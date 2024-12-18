!===========================================================================
!
!
!
!
!   /////////////////////           BEGIN            \\\\\\\\\\\\\\\\\\\\
!   \\\\\\\\\\\\\\\\\\\\\     COMMASMPI_MODULE       ////////////////////
!
!
!
!===========================================================================
!
! IN DEVELOPMENT
!

MODULE commasmpi_module

      PUBLIC

      LOGICAL,PARAMETER,PRIVATE :: debug_mpi = .false.
      logical                   :: verbose_mpi = .false.

! First grid point of tile (mpi domain)
      INTEGER,DIMENSION(:),ALLOCATABLE :: ibmpi, jbmpi, kbmpi

! Last grid point of tile (mpi domain)
      INTEGER,DIMENSION(:),ALLOCATABLE :: iempi, jempi, kempi

      INTEGER :: tileindex = 0
      INTEGER :: nxt, nyt, nzt = -1
      INTEGER :: nxbeg, nybeg, nzbeg
      INTEGER :: nxend, nyend, nzend
      INTEGER :: itile, jtile, ktile

      integer :: ncxe = -1
      integer :: ncye = -1
      integer :: ncze = -1

! First grid point of tile (full domain)
      INTEGER :: ixbeg, jybeg, kzbeg
! Last grid point of tile (full domain)
      INTEGER :: ixend, jyend, kzend

! Processor communications

      INTEGER :: process, process2
      INTEGER :: my_rank=0, number_of_processes = 1
      INTEGER :: my_comm, my_info
      INTEGER :: mpi_error_code, mpi_abort_error_code
      INTEGER :: nxpdim
      INTEGER :: nypdim
      
      INTEGER :: nproci = 1
      INTEGER :: nprocj = 1
      INTEGER :: nprock = 1
      
      INTEGER :: nxprocs = -1
      INTEGER :: nyprocs = -1
      INTEGER :: nzprocs = -1

      INTEGER :: myproci = 1
      INTEGER :: myprocj = 1
      INTEGER :: myprock = 1
      
      integer, allocatable :: procmap(:,:,:)
      
      INTEGER,DIMENSION(:),ALLOCATABLE :: iper ! offsets for periodic
      INTEGER,DIMENSION(:),ALLOCATABLE :: jper
      
      LOGICAL :: xperiodic = .false.
      LOGICAL :: yperiodic = .false.

#if defined (MPI)

! MPI variables

      REAL,DIMENSION(:,:,:),ALLOCATABLE :: awe_send, awe_recv
      REAL,DIMENSION(:,:,:),ALLOCATABLE :: asn_send, asn_recv

      REAL,DIMENSION(:,:,:),ALLOCATABLE :: kwe_send, kwe_recv
      REAL,DIMENSION(:,:,:),ALLOCATABLE :: ksn_send, ksn_recv

      REAL,DIMENSION(:,:,:),ALLOCATABLE :: swe_send, swe_recv
      REAL,DIMENSION(:,:,:),ALLOCATABLE :: ssn_send, ssn_recv

      INTEGER,DIMENSION(:),ALLOCATABLE :: w_proc, e_proc, s_proc, n_proc

      INTEGER,DIMENSION(:),ALLOCATABLE :: u_proc, d_proc  ! up,down for vertical tiling

! *********************************************************************
! *********************************************************************

CONTAINS

! *********************************************************************
! *********************************************************************

      SUBROUTINE commasmpi_startup()
          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER :: mpi_error_code, mpi_abort_error_code

          IF (debug_mpi) THEN
              WRITE (0,*) "Enter commasmpi_startup"
              WRITE (0,*) "commasmpi_startup: about to CALL MPI_Init"
          END IF !! (debug_mpi)

          CALL MPI_Init(mpi_error_code)

          IF (debug_mpi) THEN
              WRITE (0,*) "commasmpi_startup: done calling  MPI_Init,",       &
     &           " mpi_error_code=", mpi_error_code
          END IF !! (debug_mpi)

          IF (mpi_error_code /= MPI_SUCCESS) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) "commasmpi_startup: ",                          &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                  mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) "commasmpi_startup: ",                          &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (mpi_error_code /= MPI_SUCCESS)

          my_comm = MPI_COMM_WORLD
          my_info = MPI_INFO_NULL

          IF (debug_mpi) THEN
              WRITE (0,*) "commasmpi_startup: ",                              &
     &            "about to CALL MPI_Comm_rank"
          END IF !! (debug_mpi)

          CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, mpi_error_code)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "commasmpi_startup: ",                     &
     &            "done calling  MPI_Comm_rank, mpi_error_code=",             &
     &            mpi_error_code, ", my_rank=", my_rank
          END IF !! (debug_mpi)

          IF (mpi_error_code /= MPI_SUCCESS) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_startup: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_startup: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (mpi_error_code /= MPI_SUCCESS)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "commasmpi_startup: ",                     &
     &            "about to CALL MPI_Comm_size"
          END IF !! (debug_mpi)

          CALL MPI_Comm_size(MPI_COMM_WORLD, number_of_processes,             &
     &             mpi_error_code)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "commasmpi_startup: ",                     &
     &            "done calling  MPI_Comm_rank, ",                            &
     &            "mpi_error_code=",                                          &
     &            mpi_error_code,                                             &
     &            ", number_of_processes=", number_of_processes
          END IF !! (debug_mpi)

          IF (mpi_error_code /= MPI_SUCCESS) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_startup: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_startup: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (mpi_error_code /= MPI_SUCCESS)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit commasmpi_startup"
          END IF !! (debug_mpi)

      END SUBROUTINE commasmpi_startup

! *********************************************************************
! *********************************************************************

      SUBROUTINE commasmpi_bounds(nx,ny,nz)
          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER, INTENT(INOUT) :: nx,ny,nz

          INTEGER,PARAMETER :: memory_success = 0

          INTEGER :: the_process, i, j, k
          integer :: iproc, jproc, kproc
          INTEGER :: memory_status

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter commasmpi_bounds"
          END IF !! (debug_mpi)

          ALLOCATE(ibmpi(0 : number_of_processes - 1),                        &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "ibmpi of bounds 0:", number_of_processes - 1

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(jbmpi(0 : number_of_processes - 1),                        &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "jbmpi of bounds 0:", number_of_processes - 1

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(kbmpi(0 : number_of_processes - 1),                        &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "kbmpi of bounds 0:", number_of_processes - 1

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(iempi(0 : number_of_processes - 1),                        &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "iempi of bounds 0:", number_of_processes - 1

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(jempi(0 : number_of_processes - 1),                        &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "jempi of bounds 0:", number_of_processes - 1

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(kempi(0 : number_of_processes - 1),                        &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "kempi of bounds 0:", number_of_processes - 1

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(w_proc(0 : number_of_processes - 1),                       &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "w_proc of bounds 0:", number_of_processes - 1

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(e_proc(0 : number_of_processes - 1),                       &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "e_proc of bounds 0:", number_of_processes - 1

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(s_proc(0 : number_of_processes - 1),                       &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "s_proc of bounds 0:", number_of_processes - 1

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(n_proc(0 : number_of_processes - 1),                       &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "n_proc of bounds 0:", number_of_processes - 1

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)


          ALLOCATE(d_proc(0 : number_of_processes - 1),                       &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "d_proc of bounds 0:", number_of_processes - 1

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(u_proc(0 : number_of_processes - 1),                       &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "u_proc of bounds 0:", number_of_processes - 1

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_bounds: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          nxbeg = 1
          nxend = nx
          nybeg = 1
          nyend = ny
          nzbeg = 1
          nzend = nz

          IF ( debug_mpi ) THEN
            write(0,*) my_rank, ': nxprocs,nxt = ',nxprocs,nxt,nx
            write(0,*) my_rank, ': nyprocs,nyt = ',nyprocs,nyt,ny
          ENDIF

          IF ( nxprocs <= 0 ) THEN
            nx = nxt
            nxprocs = nxend/nxt
          ELSE
            nxt = nxend/nxprocs
            nx = nxt
          ENDIF
          
          IF ( nyprocs <= 0 ) THEN
            ny = nyt
            nyprocs = nyend/nyt
          ELSE
            nyt = nyend/nyprocs
            ny = nyt
          ENDIF

          IF ( nzt > nz ) nzt = nz
          
          IF ( nzprocs <= 0 ) THEN
!            nzt = nz ! for now, just setting this always -- until vertical tiling is possible
            IF ( nzt <= 0 ) THEN
              nzt = nz
            ENDIF
            nz = nzt
            nzprocs = nzend/nzt
          ELSE
            nzt = nzend/nzprocs
            nz = nzt
          ENDIF

          IF ( itile .ne. nxt ) itile = nxt
          IF ( jtile .ne. nyt ) jtile = nyt
          IF ( ktile .ne. nzt ) ktile = nzt

          IF ( .not. (nxprocs*nxt .eq. nxend .or. nxprocs*nxt .eq. nxend-1) ) THEN
             
             IF ( my_rank == 0 ) THEN
               write(0,*) 'Problem with nxt!!'
               write(0,*) 'nxend,nxt,nxprocs = ',nxend,nxt,nxprocs,nxprocs*nxt
               write(0,*) 'logicals: ',nxprocs*nxt .eq. nxend,nxprocs*nxt .eq. nxend-1
             ENDIF
             
             CALL MPI_BARRIER(my_comm, mpi_error_code)
             
             CALL COMMASMPI_ABORT()
             STOP
             
          ENDIF
          
          IF ( .not. (nyprocs*nyt .eq. nyend .or. nyprocs*nyt .eq. nyend-1) ) THEN
               
             IF ( my_rank == 0 ) THEN
               write(0,*) 'Problem with nyt!!'
               write(0,*) 'nyend,nyt,nyprocs = ',nyend,nyt,nyprocs
             ENDIF
             
             CALL MPI_BARRIER(my_comm, mpi_error_code)
             
             CALL COMMASMPI_ABORT()
             STOP
             
          ENDIF

          IF ( .not. (nzprocs*nzt .eq. nzend .or. nzprocs*nzt .eq. nzend-1) ) THEN
               
             IF ( my_rank == 0 ) THEN
               write(0,*) 'Problem with nzt!!'
               write(0,*) 'nzend,nzt,nzprocs = ',nzend,nzt,nzprocs
             ENDIF
             
             CALL MPI_BARRIER(my_comm, mpi_error_code)
             
             CALL COMMASMPI_ABORT()
             STOP
             
          ENDIF
          
          IF ( nxprocs*nyprocs*nzprocs .ne. number_of_processes ) THEN
            
             IF ( my_rank == 0 ) THEN
               write(0,*) 'Problem!  nxprocs*nyprocs .ne. number_of_processes!!'
               write(0,*) 'nxprocs,nyprocs, number_of_processes = ',nxprocs,nyprocs, number_of_processes
             ENDIF
             
             CALL MPI_BARRIER(my_comm, mpi_error_code)

             CALL COMMASMPI_ABORT()
             STOP

          ENDIF
          
!          nzt = nz
!          ktile = nzend

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "nxbeg=", nxbeg, ", nxend=", nxend
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "nybeg=", nybeg, ", nyend=", nyend
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "nzbeg=", nzbeg, ", nzend=", nzend
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "nxt=", nxt, ", nyt=", nyt, ", nzt=", nzt
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "nx=", nx, ", ny=", ny, ", nz=", nz
          END IF !! (debug_mpi)

          nproci = nxend/nxt
          nprocj = nyend/nyt
          nprock = nzend/nzt
          allocate ( procmap( nproci, nprocj, nprock) )
          allocate ( iper( nproci ) )
          allocate ( jper( nprocj ) )
          
          IF (number_of_processes == 1) THEN

              ibmpi(my_rank) = 1
              iempi(my_rank) = nxt
              jbmpi(my_rank) = 1
              jempi(my_rank) = nyt
              kbmpi(my_rank) = 1
              kempi(my_rank) = nzt
              procmap(1,1,1) = 0

          ELSE   !! (number_of_processes > 1)

              the_process = 0
              iproc = 1
              jproc = 1
              kproc = 1
              DO k = 1, nzt*nprock, nzt
                DO j = 1, nyt*nprocj, nyt
                  DO i = 1, nxt*nproci, nxt
                      ibmpi(the_process) = i
                      iempi(the_process) =                                    &
     &                    ibmpi(the_process) + nxt - 1
                      IF (iempi(the_process) > nxend) THEN
                          iempi(the_process) = nxend
                      END IF !! (iempi(the_process) > nxend)
                      jbmpi(the_process) = j
                      jempi(the_process) =                                    &
     &                    jbmpi(the_process) + nyt - 1
                      IF (jempi(the_process) > nyend) THEN
                          jempi(the_process) = nyend
                      END IF !! (jempi(the_process) > nyend)
                      kbmpi(the_process) = k
                      kempi(the_process) = nzt
                      kempi(the_process) =                                    &
     &                    kbmpi(the_process) + nzt - 1
                      IF (kempi(the_process) > nzend) THEN
                          kempi(the_process) = nzend
                      END IF !! (kempi(the_process) > nzend)
                      
                      procmap(iproc,jproc,kproc) = the_process 
                      
                      IF ( my_rank == the_process ) THEN
                        myproci = iproc
                        myprocj = jproc
                        myprock = kproc
                      ENDIF
                      
                      the_process = the_process  + 1
                      iproc = iproc + 1
                  END DO
                   iproc = 1
                   jproc = jproc + 1
                END DO
                  iproc = 1
                  jproc = 1
                  kproc = kproc + 1
              END DO

          END IF !! (number_of_processes == 1)...ELSE

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "ibmpi=", ibmpi(0:number_of_processes-1)
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "iempi=", iempi(0:number_of_processes-1)
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "jbmpi=", jbmpi(0:number_of_processes-1)
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "jempi=", jempi(0:number_of_processes-1)
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "kbmpi=", kbmpi(0:number_of_processes-1)
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "kempi=", kempi(0:number_of_processes-1)
          END IF !! (debug_mpi)

          ixbeg = ibmpi(my_rank)
          ixend = iempi(my_rank)
          jybeg = jbmpi(my_rank)
          jyend = jempi(my_rank)
          kzbeg = kbmpi(my_rank)
          kzend = kempi(my_rank)
          
          IF ( ixend .eq. nxend - 1 ) THEN
            nxt = nxt + 1
            itile = nxt
            nx = nx + 1
            ixend = nxend
!            write(0,*) 'resetting nxt for rank ',my_rank
          ENDIF
          IF ( jyend .eq. nyend - 1 ) THEN
            nyt = nyt + 1
            jtile = nyt
            ny = ny + 1
            jyend = nyend
!            write(0,*) 'resetting nyt for rank ',my_rank
          ENDIF

          IF ( kzend .eq. nzend - 1 ) THEN
            nzt = nzt + 1
            ktile = nzt
            nz = nz + 1
            kzend = nzend
!            write(0,*) 'resetting nzt for rank ',my_rank
          ENDIF

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "ixbeg=", ixbeg, ", ixend=", ixend
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "jybeg=", jybeg, ", jyend=", jyend
              WRITE (0,*) my_rank, "commasmpi_bounds: ",                     &
     &            "kzbeg=", kzbeg, ", kzend=", kzend
              WRITE (0,*) my_rank, "Exit commasmpi_bounds"
          END IF !! (debug_mpi)

      END SUBROUTINE commasmpi_bounds

! *********************************************************************
! *********************************************************************

      SUBROUTINE commasmpi_links(bcx,bcy)

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER :: process, process2
          INTEGER :: bcx,bcy
          INTEGER :: i,j

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter commasmpi_links"
          END IF !! (debug_mpi)

          DO process = 0, number_of_processes - 1

              w_proc(process) = MPI_UNDEFINED
              e_proc(process) = MPI_UNDEFINED
              s_proc(process) = MPI_UNDEFINED
              n_proc(process) = MPI_UNDEFINED
              d_proc(process) = MPI_UNDEFINED
              u_proc(process) = MPI_UNDEFINED

          END DO !! process

          DO process = 0, number_of_processes - 1

              DO process2 = 0, number_of_processes - 1

                  IF (process /= process2) THEN

                      IF ((ibmpi(process) == iempi(process2)+1)  .AND.       &
     &                    (jbmpi(process) == jbmpi(process2)  )  .AND.       &
     &                    (kbmpi(process) == kbmpi(process2)  )) THEN

                          w_proc(process) = process2

                      END IF 

                      IF ((iempi(process) == ibmpi(process2)-1)  .AND.       &
     &                    (jbmpi(process) == jbmpi(process2)  )  .AND.       &
     &                    (kbmpi(process) == kbmpi(process2)  )) THEN

                          e_proc(process) = process2

                      END IF

                      IF ((ibmpi(process) == ibmpi(process2)  )  .AND.       &
     &                    (jbmpi(process) == jempi(process2)+1)  .AND.       &
     &                    (kbmpi(process) == kbmpi(process2)  )) THEN

                          s_proc(process) = process2

                      END IF

                      IF ((ibmpi(process) == ibmpi(process2)  )  .AND.       &
     &                    (jempi(process) == jbmpi(process2)-1)  .AND.       &
     &                    (kbmpi(process) == kbmpi(process2)  )) THEN

                          n_proc(process) = process2

                      END IF

                      IF ((ibmpi(process) == ibmpi(process2)  )  .AND.       &
     &                    (jbmpi(process) == jbmpi(process2)  )  .AND.       &
     &                    (kbmpi(process) == kempi(process2)+1)) THEN

                          d_proc(process) = process2

                      END IF 

                      IF ((ibmpi(process) == ibmpi(process2)  )  .AND.       &
     &                    (jbmpi(process) == jbmpi(process2)  )  .AND.       &
     &                    (kempi(process) == kbmpi(process2)-1)) THEN

                          u_proc(process) = process2

                      END IF 


                  END IF !! (process /= process2)

              END DO !! process2
 
          END DO !! process
          
          IF ( bcy == 2 ) THEN ! North-South periodic
            
            IF ( nprocj < 2 ) THEN
              write(0,*) 'ERROR: must have nprocj >= 2 for North-South periodic!'
              CALL COMMASMPI_SHUTDOWN()
              STOP
            ENDIF
            
            yperiodic = .true.
            
            DO i = 1,nproci
              process  = procmap(i,nprocj,myprock)
              process2 = procmap(i,1,     myprock)
              
              n_proc(process)  = process2
              s_proc(process2) = process
              
            ENDDO
          ENDIF

          IF ( bcx == 2 ) THEN ! East-West periodic
            IF ( nproci < 2 ) THEN
              write(0,*) 'ERROR: must have nproci >= 2 for East-West periodic!'
              CALL COMMASMPI_SHUTDOWN()
              STOP
            ENDIF

            xperiodic = .true.

            DO j = 1,nprocj
              process  = procmap(1     ,j,myprock)
              process2 = procmap(nproci,j,myprock)
              
              w_proc(process)  = process2
              e_proc(process2) = process
            ENDDO
          ENDIF

!     WRITE (*,"("****************** my_rank:",1x,i3," ******************")") my_rank
!     WRITE (*,*) "    process  ixbeg  ixend   w_proc     e_proc"
!     WRITE (*,*) "    process  jybeg  jyend   s_proc     n_proc"
!     WRITE (*,"(1x,3(2x,i4),5x,i6,5x,i6)") process,ibmpi(process),iempi(process),w_proc(process),e_proc(process)
!     WRITE (*,"(1x,3(2x,i4),5x,i6,5x,i6)") process,jbmpi(process),jempi(process),s_proc(process),n_proc(process)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "commasmpi_links: ",                       &
     &            "w_proc=", w_proc(0:number_of_processes - 1)
              WRITE (0,*) my_rank, "commasmpi_links: ",                       &
     &            "e_proc=", e_proc(0:number_of_processes - 1)
              WRITE (0,*) my_rank, "commasmpi_links: ",                       &
     &            "s_proc=", s_proc(0:number_of_processes - 1)
              WRITE (0,*) my_rank, "commasmpi_links: ",                       &
     &            "n_proc=", n_proc(0:number_of_processes - 1)
            IF ( nprock > 1 ) THEN
              WRITE (0,*) my_rank, "commasmpi_links: ",                       &
     &            "d_proc=", d_proc(0:number_of_processes - 1)
              WRITE (0,*) my_rank, "commasmpi_links: ",                       &
     &            "u_proc=", u_proc(0:number_of_processes - 1)
            ENDIF
              WRITE (0,*) my_rank, "Exit commasmpi_links"
          END IF !! (debug_mpi)

      END SUBROUTINE commasmpi_links

! *********************************************************************
! *********************************************************************

      SUBROUTINE commasmpi_allocate(ng)
      
      implicit none
      
      integer ng

!         ALLOCATE(awe_send(-ng+1:0,-ng+1:nyt+ng,-ng+1:nzt+ng))
!         ALLOCATE(awe_recv(-ng+1:0,-ng+1:nyt+ng,-ng+1:nzt+ng))
!         ALLOCATE(asn_send(-ng+1:nxt+ng,-ng+1:0,-ng+1:nzt+ng))
!         ALLOCATE(asn_recv(-ng+1:nxt+ng,-ng+1:0,-ng+1:nzt+ng))

!         ALLOCATE(kwe_send(-ng+1:0,-ng+1:nyt+ng,-ng+1:nzt+ng))
!         ALLOCATE(kwe_recv(-ng+1:0,-ng+1:nyt+ng,-ng+1:nzt+ng))
!         ALLOCATE(ksn_send(-ng+1:nxt+ng,-ng+1:0,-ng+1:nzt+ng))
!         ALLOCATE(ksn_recv(-ng+1:nxt+ng,-ng+1:0,-ng+1:nzt+ng))

!         ALLOCATE(swe_send(-ng+1:0,-ng+1:nyt+ng,-ng+1:nzt+ng))
!         ALLOCATE(swe_recv(-ng+1:0,-ng+1:nyt+ng,-ng+1:nzt+ng))
!         ALLOCATE(ssn_send(-ng+1:nxt+ng,-ng+1:0,-ng+1:nzt+ng))
!         ALLOCATE(ssn_recv(-ng+1:nxt+ng,-ng+1:0,-ng+1:nzt+ng))

      END SUBROUTINE commasmpi_allocate

! *********************************************************************
! *********************************************************************

      SUBROUTINE commasmpi_shutdown ()

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER :: mpi_error_code, mpi_abort_error_code

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter commasmpi_shutdown"
          END IF !! (debug_mpi)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "commasmpi_shutdown: ",                    &
     &            "about to CALL MPI_Finalize"
          END IF !! (debug_mpi)

          CALL MPI_Finalize(mpi_error_code)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "commasmpi_shutdown: ",                    &
     &            "done calling  MPI_Finalize, mpi_error_code=",              &
     &            mpi_error_code
          END IF !! (debug_mpi)

          IF (mpi_error_code /= MPI_SUCCESS) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_shutdown: ",                &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_shutdown: ",                &
     &                "done calling  MPI_Abort, mpi_abort_error_code=",       &
     &                mpi_abort_error_code
              END IF !! (debug_mpi)

          END IF

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit commasmpi_shutdown"
          END IF !! (debug_mpi)

      END SUBROUTINE commasmpi_shutdown

! *********************************************************************
! *********************************************************************

      SUBROUTINE commasmpi_abort ()

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER :: mpi_error_code, mpi_abort_error_code

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter commasmpi_abort"
          END IF !! (debug_mpi)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "commasmpi_abort: ",                    &
     &            "about to CALL MPI_Abort"
          END IF !! (debug_mpi)


              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "commasmpi_abort: ",                &
     &                "done calling  MPI_Abort, mpi_abort_error_code=",       &
     &                mpi_abort_error_code
              END IF !! (debug_mpi)


          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit commasmpi_abort"
          END IF !! (debug_mpi)
          
          RETURN

      END SUBROUTINE commasmpi_abort

!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE SENDRECV_WESTWARD  <<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------
      SUBROUTINE sendrecv_westward (                                          &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                    &
     &               west_neighbor_process, east_neighbor_process,            &
     &               westward_tag, variable)

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: west_neighbor_process
          INTEGER,INTENT(IN) :: east_neighbor_process
          INTEGER,INTENT(IN) :: westward_tag

          REAL,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),              &
     &     INTENT(INOUT) :: variable

          INTEGER,PARAMETER :: memory_success = 0

          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: westward_send_buffer
          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: westward_recv_buffer

          INTEGER :: mpi_status(MPI_Status_size)
          INTEGER :: westward_send_size, westward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_western_neighbor, has_eastern_neighbor
          INTEGER :: offset

          offset = 0
          
          IF ( xperiodic .and. myproci == nproci ) THEN
              offset = 1
          ENDIF


          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_westward"
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "west_neighbor_process=", west_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "east_neighbor_process=", east_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "westward_tag=", westward_tag
          END IF !! (debug_mpi)

!........ Send process's west side to east side halo of the neighbor on
!........ the west and recv into process's east side halo the west side
!........ of the neighbor on the east.
!........ west process <-- process <-- east process 

          has_western_neighbor =                                              &
     &        (west_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (west_neighbor_process /= my_rank)
          has_eastern_neighbor =                                              &
     &        (east_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (east_neighbor_process /= my_rank)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "has_western_neighbor=", has_western_neighbor,              &
     &            "has_eastern_neighbor=", has_eastern_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(westward_send_buffer(                                      &
     &                 -ng  + 1 : 0,                                           &
     &                 -ngy + 1 : nyt + ngy,                                    &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "westward_send_buffer of bounds ",                          &
     &            -ng + 1,":",0,",",                                          &
     &            -ngy + 1,":",nyt + ngy, ",",                                  &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(westward_recv_buffer(                                      &
     &                 -ng + 1 : 0,                                           &
     &                 -ngy + 1 : nyt + ngy,                                    &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "westward_recv_buffer of bounds ",                          &
     &            -ng + 1,":",0,",",                                          &
     &            -ngy + 1,":",nyt + ngy,",",                                   &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          westward_send_buffer = 0.0
          westward_recv_buffer = 0.0
          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (has_western_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "about to prepare westward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO j = -ngy+1,nyt+ngy

                      DO ghost_zone = -ng+1, 0

                          westward_send_buffer(ghost_zone,j,k,n) =            &
     &                        variable(-ghost_zone+1,j,k,n)

                      END DO !! ghost_zone

                  END DO !! j

              END DO !! k

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "done preparing westward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_western_neighbor)

          westward_send_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nzt + ngz * 2) *                                                &
     &        (na)
          westward_recv_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nzt + ngz * 2) *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            " westward_send_size=", westward_send_size
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            " westward_recv_size=", westward_recv_size
          END IF !! (debug_mpi)

          IF (has_western_neighbor .or. has_eastern_neighbor) THEN

              IF (.NOT. has_western_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     westward_recv_buffer,                              &
     &                     westward_recv_size,                                &
     &                     MPI_REAL,                                          &
     &                     east_neighbor_process,                             &
     &                     westward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_eastern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     westward_send_buffer,                              &
     &                     westward_send_size,                                &
     &                     MPI_REAL,                                          &
     &                     west_neighbor_process,                             &
     &                     westward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_eastern_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)

                  CALL MPI_Sendrecv(                                          &
     &                     westward_send_buffer,                              &
     &                     westward_send_size,                                &
     &                     MPI_REAL,                                          &
     &                     west_neighbor_process,                             &
     &                     westward_tag,                                      &
     &                     westward_recv_buffer,                              &
     &                     westward_recv_size,                                &
     &                     MPI_REAL,                                          &
     &                     east_neighbor_process,                             &
     &                     westward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "done calling  MPI_Sendrecv, ",                     &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_eastern_neighbor)...ELSE

          END IF !! (has_western_neighbor .or. has_eastern_neighbor)


          IF (has_eastern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "about to extract westward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO j = -ngy+1,nyt+ngy

                      DO ghost_zone = -ng+1,0

                          variable(nxt-ghost_zone+1-offset,j,k,n) =                  &
     &                        westward_recv_buffer(ghost_zone,j,k,n)

                      END DO !! ghost_zone

                   END DO !! j

               END DO !! k

               END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "done extracting westward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_westward"
          END IF !! (debug_mpi)
          
          DEALLOCATE( westward_send_buffer )
          DEALLOCATE( westward_recv_buffer )

      END SUBROUTINE sendrecv_westward

!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE SENDRECV_EASTWARD  <<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------

      SUBROUTINE sendrecv_eastward (                                          &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                                   &
     &               west_neighbor_process, east_neighbor_process,            &
     &               eastward_tag, variable)

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: west_neighbor_process
          INTEGER,INTENT(IN) :: east_neighbor_process
          INTEGER,INTENT(IN) :: eastward_tag

          REAL,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),              &
     &     INTENT(INOUT) :: variable

          INTEGER,PARAMETER :: memory_success = 0

          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: eastward_send_buffer
          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: eastward_recv_buffer

          INTEGER :: mpi_status(MPI_Status_size)
          INTEGER :: eastward_send_size, eastward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_western_neighbor, has_eastern_neighbor
          INTEGER :: offset
          
          offset = 0
          
          IF ( xperiodic .and. myproci == nproci ) THEN
              offset = 1
          ENDIF
          

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_eastward"
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "west_neighbor_process=", west_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "east_neighbor_process=", east_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "eastward_tag=", eastward_tag
          END IF !! (debug_mpi)

!........ Send process's east side to west side halo of the neighbor on
!........ the east and recv into process's west side halo the east side
!........ of the neighbor on the west.
!........ west process --> process --> east process 

          has_western_neighbor =                                              &
     &        (west_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (west_neighbor_process /= my_rank)
          has_eastern_neighbor =                                              &
     &        (east_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (east_neighbor_process /= my_rank)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "has_western_neighbor=", has_western_neighbor,              &
     &            "has_eastern_neighbor=", has_eastern_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(eastward_send_buffer(                                      &
     &                 -ng + 1 : 0,                                           &
     &                 -ngy + 1 : nyt + ngy,                                    &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "eastward_send_buffer of bounds ",                          &
     &            -ng + 1,":",0,",",                                          &
     &            -ngy + 1,":",nyt + ngy,",",                                   &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(eastward_recv_buffer(                                      &
     &                 -ng + 1 : 0,                                           &
     &                 -ngy + 1 : nyt + ngy,                                    &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "eastward_recv_buffer of bounds ",                          &
     &            -ng + 1,":",0,",",                                          &
     &            -ngy + 1,":",nyt + ngy, ",",                                  &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          eastward_send_buffer = 0.0
          eastward_recv_buffer = 0.0
          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (has_eastern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "about to prepare eastward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO j = -ngy+1,nyt+ngy

                      DO ghost_zone = -ng+1,0

                          eastward_send_buffer(ghost_zone,j,k,n) =            &
                              variable(nxt+ghost_zone-offset,j,k,n)

                      END DO !! ghost_zone

                  END DO !! j

              END DO !! k

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "done preparing eastward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_western_neighbor)

          eastward_send_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nzt + ngz * 2) *                                                &
     &        (na)
          eastward_recv_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nzt + ngz * 2) *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            " eastward_send_size=", eastward_send_size
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            " eastward_recv_size=", eastward_recv_size
          END IF !! (debug_mpi)

          IF (has_western_neighbor .or. has_eastern_neighbor) THEN

              IF (.NOT. has_eastern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     eastward_recv_buffer,                              &
     &                     eastward_recv_size,                                &
     &                     MPI_REAL,                                          &
     &                     west_neighbor_process,                             &
     &                     eastward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_western_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     eastward_send_buffer,                              &
     &                     eastward_send_size,                                &
     &                     MPI_REAL,                                          &
     &                     east_neighbor_process,                             &
     &                     eastward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_eastern_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)

                  CALL MPI_Sendrecv(                                          &
     &                     eastward_send_buffer,                              &
     &                     eastward_send_size,                                &
     &                     MPI_REAL,                                          &
     &                     east_neighbor_process,                             &
     &                     eastward_tag,                                      &
     &                     eastward_recv_buffer,                              &
     &                     eastward_recv_size,                                &
     &                     MPI_REAL,                                          &
     &                     west_neighbor_process,                             &
     &                     eastward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "done calling MPI_Sendrecv, ",                      &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_eastern_neighbor)...ELSE

          END IF !! (has_western_neighbor .or. has_eastern_neighbor)


          IF (has_western_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "about to extract eastward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO j = -ngy+1,nyt+ngy

                      DO ghost_zone = -ng+1,0

                          variable(ghost_zone,j,k,n) =                        &
                              eastward_recv_buffer(ghost_zone,j,k,n)

                      END DO !! ghost_zone

                   END DO !! j

               END DO !! k

               END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "done extracting eastward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_eastward"
          END IF !! (debug_mpi)

          DEALLOCATE( eastward_send_buffer )
          DEALLOCATE( eastward_recv_buffer )

      END SUBROUTINE sendrecv_eastward

      
!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE SENDRECV_S_WRD  <<<<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------
      SUBROUTINE sendrecv_southward (                                         &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                                   &
     &               north_neighbor_process, south_neighbor_process,          &
     &               southward_tag, variable, voffset)

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: north_neighbor_process
          INTEGER,INTENT(IN) :: south_neighbor_process
          INTEGER,INTENT(IN) :: southward_tag

          REAL,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),              &
     &     INTENT(INOUT) :: variable

          INTEGER,OPTIONAL,INTENT(IN) :: voffset(na)
          
          INTEGER,PARAMETER :: memory_success = 0

          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: southward_send_buffer
          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: southward_recv_buffer

          INTEGER :: mpi_status(MPI_Status_size)
          INTEGER :: southward_send_size, southward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_northern_neighbor, has_southern_neighbor
          INTEGER :: offset(na),isoffset(na)
          
          offset(:) = 0
          isoffset(:) = 0
          
          IF ( yperiodic .and. myprocj == nprocj ) THEN
            IF ( present( voffset ) ) THEN
              offset(:) = voffset(:)
            ELSE
              offset(:) = 1
            ENDIF
          ENDIF

          IF ( yperiodic .and. myprocj == 1 ) THEN
            IF ( present( voffset ) ) THEN
              isoffset(:) = voffset(:)
            ENDIF
          ENDIF

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_southward"
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "north_neighbor_process=", north_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "south_neighbor_process=", south_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "southward_tag=", southward_tag
          END IF !! (debug_mpi)

!.......  Send process's south side to north side halo of the neighbor on the south
!.......  and recv into process's north side halo the south side of the neighbor 
!.......  on the north.
!.......  north process --> process --> south process 

          has_northern_neighbor =                                             &
     &        (north_neighbor_process /= MPI_UNDEFINED) .AND.                 &
     &        (north_neighbor_process /= my_rank)
          has_southern_neighbor =                                             &
     &        (south_neighbor_process /= MPI_UNDEFINED) .AND.                 &
     &        (south_neighbor_process /= my_rank)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "has_northern_neighbor=", has_northern_neighbor,            &
     &            "has_southern_neighbor=", has_southern_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(southward_send_buffer(                                     &
     &                 -ngx + 1 : nxt + ngx,                                    &
     &                 -ng + 1 : 0,                                           &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "southward_send_buffer of bounds ",                         &
     &            -ngx + 1,":",nxt + ngx, ",",                                  &
     &            -ng + 1,":",0,",",                                          &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(southward_recv_buffer(                                     &
     &                 -ngx + 1 : nxt + ngx,                                    &
     &                 -ng + 1 : 0,                                           &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "southward_recv_buffer of bounds ",                         &
     &            -ngx + 1,":",nxt + ngx,",",                                   &
     &            -ng + 1,":",0,",",                                          &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          southward_send_buffer = 0.0
          southward_recv_buffer = 0.0
          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (has_southern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "about to prepare southward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO i = -ngx+1,nxt+ngx

                      DO ghost_zone = -ng+1,0

                          southward_send_buffer(i,ghost_zone,k,n) =           &
                              variable(i,-ghost_zone+1+isoffset(na),k,n)

                      END DO !! ghost_zone

                  END DO !! i

              END DO !! k

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "done preparing southward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_northern_neighbor)

          southward_send_size =                                               &
     &        (nxt + ngx * 2) *                             &
     &        (ng)           *                             &
     &        (nzt + ngz * 2) *                                                &
     &        (na)
          southward_recv_size =                                               &
     &        (nxt + ngx * 2) *                             &
     &        (ng)           *                             &
     &        (nzt + ngz * 2) *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            " southward_send_size=", southward_send_size
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            " southward_recv_size=", southward_recv_size
          END IF !! (debug_mpi)

          IF (has_northern_neighbor .or. has_southern_neighbor) THEN

              IF (.NOT. has_southern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     southward_recv_buffer,                             &
     &                     southward_recv_size,                               &
     &                     MPI_REAL,                                          &
     &                     north_neighbor_process,                            &
     &                     southward_tag,                                     &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_northern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     southward_send_buffer,                             &
     &                     southward_send_size,                               &
     &                     MPI_REAL,                                          &
     &                     south_neighbor_process,                            &
     &                     southward_tag,                                     &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_southern_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)

                  CALL MPI_Sendrecv(                                          &
     &                     southward_send_buffer,                             &
     &                     southward_send_size,                               &
     &                     MPI_REAL,                                          &
     &                     south_neighbor_process,                            &
     &                     southward_tag,                                     &
     &                     southward_recv_buffer,                             &
     &                     southward_recv_size,                               &
     &                     MPI_REAL,                                          &
     &                     north_neighbor_process,                            &
     &                     southward_tag,                                     &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "done calling  MPI_Sendrecv, ",                     &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_southern_neighbor)...ELSE

          END IF !! (has_northern_neighbor .or. has_southern_neighbor)


          IF (has_northern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "about to extract southward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO i = -ngx+1,nxt+ngx

                      DO ghost_zone = -ng+1,0

                          variable(i,nyt-ghost_zone+1-offset(na),k,n) =           &
                              southward_recv_buffer(i,ghost_zone,k,n)

                      END DO !! ghost_zone

                   END DO !! i

               END DO !! k

               END DO !! n


              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "done extracting southward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_southward"
          END IF !! (debug_mpi)

          DEALLOCATE( southward_send_buffer )
          DEALLOCATE( southward_recv_buffer )

      END SUBROUTINE sendrecv_southward

!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE SENDRECV_N_WRD  <<<<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------
      SUBROUTINE sendrecv_northward (                                         &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                                   &
     &               north_neighbor_process, south_neighbor_process,          &
     &               northward_tag, variable, voffset )

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: north_neighbor_process
          INTEGER,INTENT(IN) :: south_neighbor_process
          INTEGER,INTENT(IN) :: northward_tag

          REAL,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),              &
     &     INTENT(INOUT) :: variable

          INTEGER,OPTIONAL,INTENT(IN) :: voffset(na)

          INTEGER,PARAMETER :: memory_success = 0

          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: northward_send_buffer
          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: northward_recv_buffer

          INTEGER :: mpi_status(MPI_Status_size)
          INTEGER :: northward_send_size, northward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_northern_neighbor, has_southern_neighbor

          INTEGER :: offset(na)
          
          offset(:) = 0
          
          IF ( yperiodic .and. myprocj == nprocj ) THEN
            IF ( present( voffset ) ) THEN
              offset(:) = voffset(:)
            ELSE
              offset(:) = 1
            ENDIF
          ENDIF
          
          

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_northward"
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "south_neighbor_process=", south_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "north_neighbor_process=", north_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "northward_tag=", northward_tag
          END IF !! (debug_mpi)

!.......  Send process's north side to south side halo of the neighbor on the north
!.......  and recv into process's south side halo the north side of the neighbor 
!.......  on the south.
!.......  north process <-- process <-- south process 

          has_northern_neighbor =                                             &
     &        (north_neighbor_process /= MPI_UNDEFINED) .AND.                 &
     &        (north_neighbor_process /= my_rank)
          has_southern_neighbor =                                             &
     &        (south_neighbor_process /= MPI_UNDEFINED) .AND.                 &
     &        (south_neighbor_process /= my_rank)


          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "has_northern_neighbor=", has_northern_neighbor,            &
     &            "has_southern_neighbor=", has_southern_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(northward_send_buffer(                                     &
     &                 -ngx + 1 : nxt + ngx,                                    &
     &                 -ng + 1 : 0,                                           &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "northward_send_buffer of bounds ",                         &
     &            -ngx + 1,":",nxt + ngx,",",                                   &
     &            -ng + 1,":",0,",",                                          &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(northward_recv_buffer(                                     &
     &                 -ngx + 1 : nxt + ngx,                                    &
     &                 -ng + 1 : 0,                                           &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "northward_recv_buffer of bounds ",                         &
     &            -ngx + 1,":",nxt + ngx,",",                                   &
     &            -ng + 1,":",0,",",                                          &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          northward_send_buffer = 0.0
          northward_recv_buffer = 0.0
          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (has_northern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "about to prepare northward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO i = -ngx+1,nxt+ngx

                      DO ghost_zone = -ng+1,0

                          northward_send_buffer(i,ghost_zone,k,n) =           &
                             variable(i,nyt+ghost_zone-offset(na),k,n)

                      END DO !! ghost_zone

                  END DO !! i

              END DO !! k

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "done preparing northward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_southern_neighbor)

          northward_send_size =                                               &
     &        (nxt + ngx * 2) *                             &
     &        (ng)           *                             &
     &        (nzt + ngz * 2) *                                                &
     &        (na)
          northward_recv_size =                                               &
     &        (nxt + ngx * 2) *                             &
     &        (ng)           *                             &
     &        (nzt + ngz * 2) *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            " northward_send_size=", northward_send_size
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            " northward_recv_size=", northward_recv_size
          END IF !! (debug_mpi)

          IF (has_northern_neighbor .or. has_southern_neighbor) THEN

              IF (.NOT. has_northern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     northward_recv_buffer,                             &
     &                     northward_recv_size,                               &
     &                     MPI_REAL,                                          &
     &                     south_neighbor_process,                            &
     &                     northward_tag,                                     &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_southern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     northward_send_buffer,                             &
     &                     northward_send_size,                               &
     &                     MPI_REAL,                                          &
     &                     north_neighbor_process,                            &
     &                     northward_tag,                                     &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_northern_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)

                  CALL MPI_Sendrecv(                                          &
     &                     northward_send_buffer,                             &
     &                     northward_send_size,                               &
     &                     MPI_REAL,                                          &
     &                     north_neighbor_process,                            &
     &                     northward_tag,                                     &
     &                     northward_recv_buffer,                             &
     &                     northward_recv_size,                               &
     &                     MPI_REAL,                                          &
     &                     south_neighbor_process,                            &
     &                     northward_tag,                                     &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "done calling  MPI_Sendrecv, ",                     &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_northern_neighbor)...ELSE

          END IF !! (has_southern_neighbor .or. has_northern_neighbor)


          IF (has_southern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "about to extract northward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO i = -ngx+1,nxt+ngx

                      DO ghost_zone = -ng+1,0

                          variable(i,ghost_zone,k,n) =                        &
                              northward_recv_buffer(i,ghost_zone,k,n)

                      END DO !! ghost_zone

                   END DO !! i

               END DO !! k

               END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "done extracting northward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_northward"
          END IF !! (debug_mpi)

          DEALLOCATE( northward_send_buffer )
          DEALLOCATE( northward_recv_buffer )

      END SUBROUTINE sendrecv_northward

!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE sendrecv_downward  <<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------
      SUBROUTINE sendrecv_downward (                                          &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                    &
     &               down_neighbor_process, up_neighbor_process,              &
     &               downward_tag, variable)

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: down_neighbor_process
          INTEGER,INTENT(IN) :: up_neighbor_process
          INTEGER,INTENT(IN) :: downward_tag

          REAL,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),        &
     &     INTENT(INOUT) :: variable

          INTEGER,PARAMETER :: memory_success = 0

          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: downward_send_buffer
          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: downward_recv_buffer

          INTEGER :: mpi_status(MPI_Status_size)
          INTEGER :: downward_send_size, downward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_downward_neighbor, has_upward_neighbor
          

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_downward"
              WRITE (0,*) my_rank, "sendrecv_downward: ",                     &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_downward: ",                     &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_downward: ",                     &
     &            "down_neighbor_process=", down_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_downward: ",                     &
     &            "up_neighbor_process=", up_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_downward: ",                     &
     &            "downward_tag=", downward_tag
          END IF !! (debug_mpi)

!........ Send process's down side to up side halo of the neighbor on
!........ the down and recv into process's up side halo the down side
!........ of the neighbor on the up.
!........ down process <-- process <-- up process 

          has_downward_neighbor =                                             &
     &        (down_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (down_neighbor_process /= my_rank)
          has_upward_neighbor =                                               &
     &        (up_neighbor_process /= MPI_UNDEFINED) .AND.                    &
     &        (up_neighbor_process /= my_rank)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_downward: ",                     &
     &            "has_downward_neighbor=", has_downward_neighbor,            &
     &            "has_upward_neighbor=", has_upward_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(downward_send_buffer(                                      &
     &                 -ngx + 1 : nxt + ngx,                                  &
     &                 -ngy + 1 : nyt + ngy,                                  &
     &                 -ng  + 1 : 0,                                          &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "downward_send_buffer of bounds ",                          &
     &            -ngx + 1,":",nxt + ngx,',',                                 &
     &            -ngy + 1,":",nyt + ngy, ",",                                &
     &            -ng + 1,":",0

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_downward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_downward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(downward_recv_buffer(                                      &
     &                 -ngx + 1 : nxt + ngx,                                  &
     &                 -ngy + 1 : nyt + ngy,                                  &
     &                 -ng + 1 : 0,                                           &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "downward_recv_buffer of bounds ",                          &
     &            -ngx + 1,":",nxt + ngx,',',                                 &
     &            -ngy + 1,":",nyt + ngy,",",                                 &
     &            -ng + 1,":",0

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_downward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_downward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          downward_send_buffer = 0.0
          downward_recv_buffer = 0.0
          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (has_downward_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_downward: ",                 &
     &                "about to prepare downward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO ghost_zone = -ng+1, 0

                DO j = -ngy+1,nyt+ngy

                   DO i = -ngx+1,nxt+ngx

                          downward_send_buffer(i,j,ghost_zone,n) =            &
     &                        variable(i,j,-ghost_zone+1,n)

                    END DO !! i

                  END DO !! j

              END DO !! ghost_zone

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_downward: ",                 &
     &                "done preparing downward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_downward_neighbor)

          downward_send_size =                                                &
     &        (nxt + ngx * 2) *                                               &
     &        (nyt + ngy * 2) *                                               &
     &        (ng)           *                                                &
     &        (na)
          downward_recv_size =                                                &
     &        (nxt + ngx * 2) *                                               &
     &        (nyt + ngy * 2) *                                               &
     &        (ng)           *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_downward: ",                     &
     &            " downward_send_size=", downward_send_size
              WRITE (0,*) my_rank, "sendrecv_downward: ",                     &
     &            " downward_recv_size=", downward_recv_size
          END IF !! (debug_mpi)

          IF (has_downward_neighbor .or. has_upward_neighbor) THEN

              IF (.NOT. has_downward_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_downward: ",             &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     downward_recv_buffer,                              &
     &                     downward_recv_size,                                &
     &                     MPI_REAL,                                          &
     &                     up_neighbor_process,                               &
     &                     downward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_downward: ",             &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_downward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_downward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_upward_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_downward: ",             &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     downward_send_buffer,                              &
     &                     downward_send_size,                                &
     &                     MPI_REAL,                                          &
     &                     down_neighbor_process,                             &
     &                     downward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_downward: ",             &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_downward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_downward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_upward_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_downward: ",             &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)

                  CALL MPI_Sendrecv(                                          &
     &                     downward_send_buffer,                              &
     &                     downward_send_size,                                &
     &                     MPI_REAL,                                          &
     &                     down_neighbor_process,                             &
     &                     downward_tag,                                      &
     &                     downward_recv_buffer,                              &
     &                     downward_recv_size,                                &
     &                     MPI_REAL,                                          &
     &                     up_neighbor_process,                               &
     &                     downward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_downward: ",             &
     &                    "done calling  MPI_Sendrecv, ",                     &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_downward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_downward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_upward_neighbor)...ELSE

          END IF !! (has_downward_neighbor .or. has_upward_neighbor)


          IF (has_upward_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_downward: ",                 &
     &                "about to extract downward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO ghost_zone = -ng+1,0

                  DO j = -ngy+1,nyt+ngy

                      DO i = -ngx+1,nxt+ngx

                          variable(i,j,nzt-ghost_zone+1,n) =                  &
     &                        downward_recv_buffer(i,j,ghost_zone,n)

                      END DO !! i

                   END DO !! j

               END DO !! ghost_zone

               END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_downward: ",                 &
     &                "done extracting downward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_downward"
          END IF !! (debug_mpi)
          
          DEALLOCATE( downward_send_buffer )
          DEALLOCATE( downward_recv_buffer )

      END SUBROUTINE sendrecv_downward

!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE SENDRECV_upward  <<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------

      SUBROUTINE sendrecv_upward (                                          &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                                   &
     &               down_neighbor_process, up_neighbor_process,            &
     &               upward_tag, variable)

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: down_neighbor_process
          INTEGER,INTENT(IN) :: up_neighbor_process
          INTEGER,INTENT(IN) :: upward_tag

          REAL,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),              &
     &     INTENT(INOUT) :: variable

          INTEGER,PARAMETER :: memory_success = 0

          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: upward_send_buffer
          REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: upward_recv_buffer

          INTEGER :: mpi_status(MPI_Status_size)
          INTEGER :: upward_send_size, upward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_downward_neighbor, has_upward_neighbor
          INTEGER :: offset
          
          offset = 0
          
          IF ( xperiodic .and. myproci == nproci ) THEN
              offset = 1
          ENDIF
          

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_upward"
              WRITE (0,*) my_rank, "sendrecv_upward: ",                     &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_upward: ",                     &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_upward: ",                     &
     &            "down_neighbor_process=", down_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_upward: ",                     &
     &            "up_neighbor_process=", up_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_upward: ",                     &
     &            "upward_tag=", upward_tag
          END IF !! (debug_mpi)

!........ Send process's up side to down side halo of the neighbor on
!........ the up and recv into process's down side halo the up side
!........ of the neighbor on the down.
!........ down process --> process --> up process 

          has_downward_neighbor =                                              &
     &        (down_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (down_neighbor_process /= my_rank)
          has_upward_neighbor =                                              &
     &        (up_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (up_neighbor_process /= my_rank)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_upward: ",                     &
     &            "has_downward_neighbor=", has_downward_neighbor,              &
     &            "has_upward_neighbor=", has_upward_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(upward_send_buffer(                                        &
     &                 -ngx + 1 : nxt + ngx,                                  &
     &                 -ngy + 1 : nyt + ngy,                                  &
     &                 -ng + 1 : 0,                                           &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "upward_send_buffer of bounds ",                            &
     &            -ngx + 1,":",nxt + ngx,",",                                 &
     &            -ngy + 1,":",nyt + ngy,",",                                 &
     &            -ng + 1,":",0

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_upward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_upward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(upward_recv_buffer(                                        &
     &                 -ngx + 1 : nxt + ngx,                                  &
     &                 -ngy + 1 : nyt + ngy,                                  &
     &                 -ng + 1 : 0,                                           &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "upward_recv_buffer of bounds ",                            &
     &            -ngx + 1,":",nxt + ngx,",",                                 &
     &            -ngy + 1,":",nyt + ngy, ",",                                &
     &            -ng + 1,":",0

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_upward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_upward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          upward_send_buffer = 0.0
          upward_recv_buffer = 0.0
          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (has_upward_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_upward: ",                 &
     &                "about to prepare upward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na


              DO ghost_zone = -ng+1,0

                  DO j = -ngy+1,nyt+ngy

                     DO i = -ngx+1,nxt+ngx

                          upward_send_buffer(i,j,ghost_zone,n) =            &
                              variable(i,j,nzt+ghost_zone-offset,n)

                     END DO !! i

                  END DO !! j

              END DO !! ghost_zone

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_upward: ",                 &
     &                "done preparing upward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_downward_neighbor)

          upward_send_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nxt + ngx * 2) *                                                &
     &        (na)
          upward_recv_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nxt + ngx * 2) *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_upward: ",                     &
     &            " upward_send_size=", upward_send_size
              WRITE (0,*) my_rank, "sendrecv_upward: ",                     &
     &            " upward_recv_size=", upward_recv_size
          END IF !! (debug_mpi)

          IF (has_downward_neighbor .or. has_upward_neighbor) THEN

              IF (.NOT. has_upward_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_upward: ",             &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     upward_recv_buffer,                              &
     &                     upward_recv_size,                                &
     &                     MPI_REAL,                                          &
     &                     down_neighbor_process,                             &
     &                     upward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_upward: ",             &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_upward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_upward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_downward_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_upward: ",             &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     upward_send_buffer,                              &
     &                     upward_send_size,                                &
     &                     MPI_REAL,                                          &
     &                     up_neighbor_process,                             &
     &                     upward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_upward: ",             &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_upward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_upward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_upward_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_upward: ",             &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)

                  CALL MPI_Sendrecv(                                          &
     &                     upward_send_buffer,                              &
     &                     upward_send_size,                                &
     &                     MPI_REAL,                                          &
     &                     up_neighbor_process,                             &
     &                     upward_tag,                                      &
     &                     upward_recv_buffer,                              &
     &                     upward_recv_size,                                &
     &                     MPI_REAL,                                          &
     &                     down_neighbor_process,                             &
     &                     upward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_upward: ",             &
     &                    "done calling MPI_Sendrecv, ",                      &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_upward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_upward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_upward_neighbor)...ELSE

          END IF !! (has_downward_neighbor .or. has_upward_neighbor)


          IF (has_downward_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_upward: ",                 &
     &                "about to extract upward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

                DO ghost_zone = -ng+1,0

                  DO j = -ngy+1,nyt+ngy

                      DO i = -ngx+1,nxt+ngx

                          variable(i,j,ghost_zone,n) =                        &
                              upward_recv_buffer(i,j,ghost_zone,n)

                      END DO !! i

                   END DO !! j

                 END DO !! ghost_zone

               END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_upward: ",                 &
     &                "done extracting upward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_upward"
          END IF !! (debug_mpi)

          DEALLOCATE( upward_send_buffer )
          DEALLOCATE( upward_recv_buffer )

      END SUBROUTINE sendrecv_upward


!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE SENDRECV_WESTWARD  <<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------
      SUBROUTINE sendrecv_westward_dp (                                          &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                    &
     &               west_neighbor_process, east_neighbor_process,            &
     &               westward_tag, variable)

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: west_neighbor_process
          INTEGER,INTENT(IN) :: east_neighbor_process
          INTEGER,INTENT(IN) :: westward_tag

          DOUBLE PRECISION,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),              &
     &     INTENT(INOUT) :: variable

          INTEGER,PARAMETER :: memory_success = 0

          DOUBLE PRECISION,DIMENSION(:,:,:,:),ALLOCATABLE :: westward_send_buffer
          DOUBLE PRECISION,DIMENSION(:,:,:,:),ALLOCATABLE :: westward_recv_buffer

          INTEGER :: mpi_status(MPI_Status_size)
          INTEGER :: westward_send_size, westward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_western_neighbor, has_eastern_neighbor
          INTEGER :: offset

          offset = 0
          
          IF ( xperiodic .and. myproci == nproci ) THEN
              offset = 1
          ENDIF


          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_westward"
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "west_neighbor_process=", west_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "east_neighbor_process=", east_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "westward_tag=", westward_tag
          END IF !! (debug_mpi)

!........ Send process's west side to east side halo of the neighbor on
!........ the west and recv into process's east side halo the west side
!........ of the neighbor on the east.
!........ west process <-- process <-- east process 

          has_western_neighbor =                                              &
     &        (west_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (west_neighbor_process /= my_rank)
          has_eastern_neighbor =                                              &
     &        (east_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (east_neighbor_process /= my_rank)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            "has_western_neighbor=", has_western_neighbor,              &
     &            "has_eastern_neighbor=", has_eastern_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(westward_send_buffer(                                      &
     &                 -ng  + 1 : 0,                                           &
     &                 -ngy + 1 : nyt + ngy,                                    &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "westward_send_buffer of bounds ",                          &
     &            -ng + 1,":",0,",",                                          &
     &            -ngy + 1,":",nyt + ngy, ",",                                  &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(westward_recv_buffer(                                      &
     &                 -ng + 1 : 0,                                           &
     &                 -ngy + 1 : nyt + ngy,                                    &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "westward_recv_buffer of bounds ",                          &
     &            -ng + 1,":",0,",",                                          &
     &            -ngy + 1,":",nyt + ngy,",",                                   &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          westward_send_buffer = 0.0
          westward_recv_buffer = 0.0
          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (has_western_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "about to prepare westward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO j = -ngy+1,nyt+ngy

                      DO ghost_zone = -ng+1, 0

                          westward_send_buffer(ghost_zone,j,k,n) =            &
     &                        variable(-ghost_zone+1,j,k,n)

                      END DO !! ghost_zone

                  END DO !! j

              END DO !! k

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "done preparing westward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_western_neighbor)

          westward_send_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nzt + ngz * 2) *                                                &
     &        (na)
          westward_recv_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nzt + ngz * 2) *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            " westward_send_size=", westward_send_size
              WRITE (0,*) my_rank, "sendrecv_westward: ",                     &
     &            " westward_recv_size=", westward_recv_size
          END IF !! (debug_mpi)

          IF (has_western_neighbor .or. has_eastern_neighbor) THEN

              IF (.NOT. has_western_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     westward_recv_buffer,                              &
     &                     westward_recv_size,                                &
     &                     MPI_DOUBLE_PRECISION,                              &
     &                     east_neighbor_process,                             &
     &                     westward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_eastern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     westward_send_buffer,                              &
     &                     westward_send_size,                                &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     west_neighbor_process,                             &
     &                     westward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_eastern_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)

                  CALL MPI_Sendrecv(                                          &
     &                     westward_send_buffer,                              &
     &                     westward_send_size,                                &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     west_neighbor_process,                             &
     &                     westward_tag,                                      &
     &                     westward_recv_buffer,                              &
     &                     westward_recv_size,                                &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     east_neighbor_process,                             &
     &                     westward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_westward: ",             &
     &                    "done calling  MPI_Sendrecv, ",                     &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_westward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_eastern_neighbor)...ELSE

          END IF !! (has_western_neighbor .or. has_eastern_neighbor)


          IF (has_eastern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "about to extract westward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO j = -ngy+1,nyt+ngy

                      DO ghost_zone = -ng+1,0

                          variable(nxt-ghost_zone+1-offset,j,k,n) =                  &
     &                        westward_recv_buffer(ghost_zone,j,k,n)

                      END DO !! ghost_zone

                   END DO !! j

               END DO !! k

               END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_westward: ",                 &
     &                "done extracting westward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_westward"
          END IF !! (debug_mpi)
          
          DEALLOCATE( westward_send_buffer )
          DEALLOCATE( westward_recv_buffer )

      END SUBROUTINE sendrecv_westward_dp

!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE SENDRECV_EASTWARD  <<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------

      SUBROUTINE sendrecv_eastward_dp (                                          &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                                   &
     &               west_neighbor_process, east_neighbor_process,            &
     &               eastward_tag, variable)

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: west_neighbor_process
          INTEGER,INTENT(IN) :: east_neighbor_process
          INTEGER,INTENT(IN) :: eastward_tag

          DOUBLE PRECISION,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),              &
     &     INTENT(INOUT) :: variable

          INTEGER,PARAMETER :: memory_success = 0

          DOUBLE PRECISION,DIMENSION(:,:,:,:),ALLOCATABLE :: eastward_send_buffer
          DOUBLE PRECISION,DIMENSION(:,:,:,:),ALLOCATABLE :: eastward_recv_buffer

          INTEGER :: mpi_status(MPI_Status_size)
          INTEGER :: eastward_send_size, eastward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_western_neighbor, has_eastern_neighbor
          INTEGER :: offset
          
          offset = 0
          
          IF ( xperiodic .and. myproci == nproci ) THEN
              offset = 1
          ENDIF
          

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_eastward"
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "west_neighbor_process=", west_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "east_neighbor_process=", east_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "eastward_tag=", eastward_tag
          END IF !! (debug_mpi)

!........ Send process's east side to west side halo of the neighbor on
!........ the east and recv into process's west side halo the east side
!........ of the neighbor on the west.
!........ west process --> process --> east process 

          has_western_neighbor =                                              &
     &        (west_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (west_neighbor_process /= my_rank)
          has_eastern_neighbor =                                              &
     &        (east_neighbor_process /= MPI_UNDEFINED) .AND.                  &
     &        (east_neighbor_process /= my_rank)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            "has_western_neighbor=", has_western_neighbor,              &
     &            "has_eastern_neighbor=", has_eastern_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(eastward_send_buffer(                                      &
     &                 -ng + 1 : 0,                                           &
     &                 -ngy + 1 : nyt + ngy,                                    &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "eastward_send_buffer of bounds ",                          &
     &            -ng + 1,":",0,",",                                          &
     &            -ngy + 1,":",nyt + ngy,",",                                   &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(eastward_recv_buffer(                                      &
     &                 -ng + 1 : 0,                                           &
     &                 -ngy + 1 : nyt + ngy,                                    &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "eastward_recv_buffer of bounds ",                          &
     &            -ng + 1,":",0,",",                                          &
     &            -ngy + 1,":",nyt + ngy, ",",                                  &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          eastward_send_buffer = 0.0
          eastward_recv_buffer = 0.0
          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (has_eastern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "about to prepare eastward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO j = -ngy+1,nyt+ngy

                      DO ghost_zone = -ng+1,0

                          eastward_send_buffer(ghost_zone,j,k,n) =            &
                              variable(nxt+ghost_zone-offset,j,k,n)

                      END DO !! ghost_zone

                  END DO !! j

              END DO !! k

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "done preparing eastward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_western_neighbor)

          eastward_send_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nzt + ngz * 2) *                                                &
     &        (na)
          eastward_recv_size =                                                &
     &        (ng)           *                                                &
     &        (nyt + ngy * 2) *                                                &
     &        (nzt + ngz * 2) *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            " eastward_send_size=", eastward_send_size
              WRITE (0,*) my_rank, "sendrecv_eastward: ",                     &
     &            " eastward_recv_size=", eastward_recv_size
          END IF !! (debug_mpi)

          IF (has_western_neighbor .or. has_eastern_neighbor) THEN

              IF (.NOT. has_eastern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     eastward_recv_buffer,                              &
     &                     eastward_recv_size,                                &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     west_neighbor_process,                             &
     &                     eastward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_western_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     eastward_send_buffer,                              &
     &                     eastward_send_size,                                &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     east_neighbor_process,                             &
     &                     eastward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_eastern_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)

                  CALL MPI_Sendrecv(                                          &
     &                     eastward_send_buffer,                              &
     &                     eastward_send_size,                                &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     east_neighbor_process,                             &
     &                     eastward_tag,                                      &
     &                     eastward_recv_buffer,                              &
     &                     eastward_recv_size,                                &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     west_neighbor_process,                             &
     &                     eastward_tag,                                      &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_eastward: ",             &
     &                    "done calling MPI_Sendrecv, ",                      &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_eastward: ",         &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_eastern_neighbor)...ELSE

          END IF !! (has_western_neighbor .or. has_eastern_neighbor)


          IF (has_western_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "about to extract eastward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO j = -ngy+1,nyt+ngy

                      DO ghost_zone = -ng+1,0

                          variable(ghost_zone,j,k,n) =                        &
                              eastward_recv_buffer(ghost_zone,j,k,n)

                      END DO !! ghost_zone

                   END DO !! j

               END DO !! k

               END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_eastward: ",                 &
     &                "done extracting eastward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_eastward"
          END IF !! (debug_mpi)

          DEALLOCATE( eastward_send_buffer )
          DEALLOCATE( eastward_recv_buffer )

      END SUBROUTINE sendrecv_eastward_dp

      
!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE SENDRECV_S_WRD  <<<<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------
      SUBROUTINE sendrecv_southward_dp (                                         &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                                   &
     &               north_neighbor_process, south_neighbor_process,          &
     &               southward_tag, variable, voffset)

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: north_neighbor_process
          INTEGER,INTENT(IN) :: south_neighbor_process
          INTEGER,INTENT(IN) :: southward_tag

          DOUBLE PRECISION,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),              &
     &     INTENT(INOUT) :: variable

          INTEGER,OPTIONAL,INTENT(IN) :: voffset(na)
          
          INTEGER,PARAMETER :: memory_success = 0

          DOUBLE PRECISION,DIMENSION(:,:,:,:),ALLOCATABLE :: southward_send_buffer
          DOUBLE PRECISION,DIMENSION(:,:,:,:),ALLOCATABLE :: southward_recv_buffer

          INTEGER :: mpi_status(MPI_Status_size)
          INTEGER :: southward_send_size, southward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_northern_neighbor, has_southern_neighbor
          INTEGER :: offset(na),isoffset(na)
          
          offset(:) = 0
          isoffset(:) = 0
          
          IF ( yperiodic .and. myprocj == nprocj ) THEN
            IF ( present( voffset ) ) THEN
              offset(:) = voffset(:)
            ELSE
              offset(:) = 1
            ENDIF
          ENDIF

          IF ( yperiodic .and. myprocj == 1 ) THEN
            IF ( present( voffset ) ) THEN
              isoffset(:) = voffset(:)
            ENDIF
          ENDIF

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_southward"
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "north_neighbor_process=", north_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "south_neighbor_process=", south_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "southward_tag=", southward_tag
          END IF !! (debug_mpi)

!.......  Send process's south side to north side halo of the neighbor on the south
!.......  and recv into process's north side halo the south side of the neighbor 
!.......  on the north.
!.......  north process --> process --> south process 

          has_northern_neighbor =                                             &
     &        (north_neighbor_process /= MPI_UNDEFINED) .AND.                 &
     &        (north_neighbor_process /= my_rank)
          has_southern_neighbor =                                             &
     &        (south_neighbor_process /= MPI_UNDEFINED) .AND.                 &
     &        (south_neighbor_process /= my_rank)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            "has_northern_neighbor=", has_northern_neighbor,            &
     &            "has_southern_neighbor=", has_southern_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(southward_send_buffer(                                     &
     &                 -ngx + 1 : nxt + ngx,                                    &
     &                 -ng + 1 : 0,                                           &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "southward_send_buffer of bounds ",                         &
     &            -ngx + 1,":",nxt + ngx, ",",                                  &
     &            -ng + 1,":",0,",",                                          &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(southward_recv_buffer(                                     &
     &                 -ngx + 1 : nxt + ngx,                                    &
     &                 -ng + 1 : 0,                                           &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "southward_recv_buffer of bounds ",                         &
     &            -ngx + 1,":",nxt + ngx,",",                                   &
     &            -ng + 1,":",0,",",                                          &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          southward_send_buffer = 0.0
          southward_recv_buffer = 0.0
          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (has_southern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "about to prepare southward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO i = -ngx+1,nxt+ngx

                      DO ghost_zone = -ng+1,0

                          southward_send_buffer(i,ghost_zone,k,n) =           &
                              variable(i,-ghost_zone+1+isoffset(na),k,n)

                      END DO !! ghost_zone

                  END DO !! i

              END DO !! k

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "done preparing southward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_northern_neighbor)

          southward_send_size =                                               &
     &        (nxt + ngx * 2) *                             &
     &        (ng)           *                             &
     &        (nzt + ngz * 2) *                                                &
     &        (na)
          southward_recv_size =                                               &
     &        (nxt + ngx * 2) *                             &
     &        (ng)           *                             &
     &        (nzt + ngz * 2) *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            " southward_send_size=", southward_send_size
              WRITE (0,*) my_rank, "sendrecv_southward: ",                    &
     &            " southward_recv_size=", southward_recv_size
          END IF !! (debug_mpi)

          IF (has_northern_neighbor .or. has_southern_neighbor) THEN

              IF (.NOT. has_southern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     southward_recv_buffer,                             &
     &                     southward_recv_size,                               &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     north_neighbor_process,                            &
     &                     southward_tag,                                     &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_northern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     southward_send_buffer,                             &
     &                     southward_send_size,                               &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     south_neighbor_process,                            &
     &                     southward_tag,                                     &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_southern_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)

                  CALL MPI_Sendrecv(                                          &
     &                     southward_send_buffer,                             &
     &                     southward_send_size,                               &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     south_neighbor_process,                            &
     &                     southward_tag,                                     &
     &                     southward_recv_buffer,                             &
     &                     southward_recv_size,                               &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     north_neighbor_process,                            &
     &                     southward_tag,                                     &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_southward: ",            &
     &                    "done calling  MPI_Sendrecv, ",                     &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_southward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_southern_neighbor)...ELSE

          END IF !! (has_northern_neighbor .or. has_southern_neighbor)


          IF (has_northern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "about to extract southward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO i = -ngx+1,nxt+ngx

                      DO ghost_zone = -ng+1,0

                          variable(i,nyt-ghost_zone+1-offset(na),k,n) =           &
                              southward_recv_buffer(i,ghost_zone,k,n)

                      END DO !! ghost_zone

                   END DO !! i

               END DO !! k

               END DO !! n


              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_southward: ",                &
     &                "done extracting southward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_southward"
          END IF !! (debug_mpi)

          DEALLOCATE( southward_send_buffer )
          DEALLOCATE( southward_recv_buffer )

      END SUBROUTINE sendrecv_southward_dp

!-------------------------------------------------------------------------------
!
! >>>>>>>>>>>>>>>>>>>>>>>   SUBROUTINE SENDRECV_N_WRD  <<<<<<<<<<<<<<<<<<<<<<<  !
!
!-------------------------------------------------------------------------------
      SUBROUTINE sendrecv_northward_dp (                                         &
     &               nxt, nyt, nzt, ngx, ngy, ngz, ng, na,                                   &
     &               north_neighbor_process, south_neighbor_process,          &
     &               northward_tag, variable, voffset )

          IMPLICIT NONE

          INCLUDE "mpif.h"

          INTEGER,INTENT(IN) :: nxt, nyt, nzt, ngx, ngy, ngz, ng, na
          INTEGER,INTENT(IN) :: north_neighbor_process
          INTEGER,INTENT(IN) :: south_neighbor_process
          INTEGER,INTENT(IN) :: northward_tag

          DOUBLE PRECISION,                                                               &
     &     DIMENSION(-ngx+1:nxt+ngx,-ngy+1:nyt+ngy,-ngz+1:nzt+ngz,na),              &
     &     INTENT(INOUT) :: variable

          INTEGER,OPTIONAL,INTENT(IN) :: voffset(na)

          INTEGER,PARAMETER :: memory_success = 0

          DOUBLE PRECISION,DIMENSION(:,:,:,:),ALLOCATABLE :: northward_send_buffer
          DOUBLE PRECISION,DIMENSION(:,:,:,:),ALLOCATABLE :: northward_recv_buffer

          INTEGER :: mpi_status(MPI_Status_size)
          INTEGER :: northward_send_size, northward_recv_size
          INTEGER :: n, i, j, k, ghost_zone
          INTEGER :: memory_status
          INTEGER :: mpi_error_code, mpi_abort_error_code

          LOGICAL :: has_northern_neighbor, has_southern_neighbor

          INTEGER :: offset(na)
          
          offset(:) = 0
          
          IF ( yperiodic .and. myprocj == nprocj ) THEN
            IF ( present( voffset ) ) THEN
              offset(:) = voffset(:)
            ELSE
              offset(:) = 1
            ENDIF
          ENDIF
          
          

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Enter sendrecv_northward"
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "nxt=", nxt, ",nyt=", nyt, ",nzt=", nzt, ",ng=", ng
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "na=", na
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "south_neighbor_process=", south_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "north_neighbor_process=", north_neighbor_process
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "northward_tag=", northward_tag
          END IF !! (debug_mpi)

!.......  Send process's north side to south side halo of the neighbor on the north
!.......  and recv into process's south side halo the north side of the neighbor 
!.......  on the south.
!.......  north process <-- process <-- south process 

          has_northern_neighbor =                                             &
     &        (north_neighbor_process /= MPI_UNDEFINED) .AND.                 &
     &        (north_neighbor_process /= my_rank)
          has_southern_neighbor =                                             &
     &        (south_neighbor_process /= MPI_UNDEFINED) .AND.                 &
     &        (south_neighbor_process /= my_rank)


          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            "has_northern_neighbor=", has_northern_neighbor,            &
     &            "has_southern_neighbor=", has_southern_neighbor
          END IF !! (debug_mpi)

          ALLOCATE(northward_send_buffer(                                     &
     &                 -ngx + 1 : nxt + ngx,                                    &
     &                 -ng + 1 : 0,                                           &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "northward_send_buffer of bounds ",                         &
     &            -ngx + 1,":",nxt + ngx,",",                                   &
     &            -ng + 1,":",0,",",                                          &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          ALLOCATE(northward_recv_buffer(                                     &
     &                 -ngx + 1 : nxt + ngx,                                    &
     &                 -ng + 1 : 0,                                           &
     &                 -ngz + 1 : nzt + ngz,                                    &
     &                 na),                                                   &
     &             STAT=memory_status)

          IF (memory_status /= memory_success) THEN

              WRITE (0,*) my_rank, "ERROR: can't allocate ",                  &
     &            "northward_recv_buffer of bounds ",                         &
     &            -ngx + 1,":",nxt + ngx,",",                                   &
     &            -ng + 1,":",0,",",                                          &
     &            -ngz + 1,":",nzt + ngz

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "about to CALL MPI_Abort"
              END IF !! (debug_mpi)

              mpi_error_code = -1
              CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,                  &
     &                 mpi_abort_error_code)

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "done calling MPI_Abort"
              END IF !! (debug_mpi)

          END IF !! (memory_status /= memory_success)

          northward_send_buffer = 0.0
          northward_recv_buffer = 0.0
          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (has_northern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "about to prepare northward_send_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO i = -ngx+1,nxt+ngx

                      DO ghost_zone = -ng+1,0

                          northward_send_buffer(i,ghost_zone,k,n) =           &
                             variable(i,nyt+ghost_zone-offset(na),k,n)

                      END DO !! ghost_zone

                  END DO !! i

              END DO !! k

              END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "done preparing northward_send_buffer"
              END IF !! (debug_mpi)

          END IF !! (has_southern_neighbor)

          northward_send_size =                                               &
     &        (nxt + ngx * 2) *                             &
     &        (ng)           *                             &
     &        (nzt + ngz * 2) *                                                &
     &        (na)
          northward_recv_size =                                               &
     &        (nxt + ngx * 2) *                             &
     &        (ng)           *                             &
     &        (nzt + ngx * 2) *                                                &
     &        (na)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            " northward_send_size=", northward_send_size
              WRITE (0,*) my_rank, "sendrecv_northward: ",                    &
     &            " northward_recv_size=", northward_recv_size
          END IF !! (debug_mpi)

          IF (has_northern_neighbor .or. has_southern_neighbor) THEN

              IF (.NOT. has_northern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "about to call MPI_Recv"
                  END IF !! (debug_mpi)

                  CALL MPI_Recv(                                              &
     &                     northward_recv_buffer,                             &
     &                     northward_recv_size,                               &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     south_neighbor_process,                            &
     &                     northward_tag,                                     &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "done calling MPI_Recv, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE IF (.NOT. has_southern_neighbor) THEN

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "about to call MPI_Send"
                  END IF !! (debug_mpi)

                  CALL MPI_Send(                                              &
     &                     northward_send_buffer,                             &
     &                     northward_send_size,                               &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     north_neighbor_process,                            &
     &                     northward_tag,                                     &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "done calling MPI_Send, ",                          &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              ELSE   !! (.NOT. has_northern_neighbor)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "about to call MPI_Sendrecv"
                  END IF !! (debug_mpi)

                  CALL MPI_Sendrecv(                                          &
     &                     northward_send_buffer,                             &
     &                     northward_send_size,                               &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     north_neighbor_process,                            &
     &                     northward_tag,                                     &
     &                     northward_recv_buffer,                             &
     &                     northward_recv_size,                               &
     &                     MPI_DOUBLE_PRECISION,                                          &
     &                     south_neighbor_process,                            &
     &                     northward_tag,                                     &
     &                     MPI_COMM_WORLD,                                    &
     &                     mpi_status, mpi_error_code)

                  IF (debug_mpi) THEN
                      WRITE (0,*) my_rank, "sendrecv_northward: ",            &
     &                    "done calling  MPI_Sendrecv, ",                     &
     &                    "mpi_error_code=", mpi_error_code
                  END IF !! (debug_mpi)

                  IF (mpi_error_code /= MPI_SUCCESS) THEN

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "about to CALL MPI_Abort"
                      END IF !! (debug_mpi)

                      CALL MPI_Abort(MPI_COMM_WORLD, mpi_error_code,          &
     &                         mpi_abort_error_code)

                      IF (debug_mpi) THEN
                          WRITE (0,*) my_rank, "sendrecv_northward: ",        &
     &                        "done calling MPI_Abort"
                      END IF !! (debug_mpi)

                  END IF !! (mpi_error_code /= MPI_SUCCESS)

              END IF !! (.NOT. has_northern_neighbor)...ELSE

          END IF !! (has_southern_neighbor .or. has_northern_neighbor)


          IF (has_southern_neighbor) THEN

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "about to extract northward_recv_buffer"
              END IF !! (debug_mpi)

              DO n = 1,na

              DO k = -ngz+1,nzt+ngz

                  DO i = -ngx+1,nxt+ngx

                      DO ghost_zone = -ng+1,0

                          variable(i,ghost_zone,k,n) =                        &
                              northward_recv_buffer(i,ghost_zone,k,n)

                      END DO !! ghost_zone

                   END DO !! i

               END DO !! k

               END DO !! n

              IF (debug_mpi) THEN
                  WRITE (0,*) my_rank, "sendrecv_northward: ",                &
     &                "done extracting northward_recv_buffer"
              END IF !! (debug_mpi)

          END IF

          CALL MPI_BARRIER(my_comm, mpi_error_code)

          IF (debug_mpi) THEN
              WRITE (0,*) my_rank, "Exit sendrecv_northward"
          END IF !! (debug_mpi)

          DEALLOCATE( northward_send_buffer )
          DEALLOCATE( northward_recv_buffer )

      END SUBROUTINE sendrecv_northward_dp


#else
! substitute routines for non-MPI compile

   CONTAINS
   
!-------------------------------------------------------------------------------
       SUBROUTINE commasmpi_startup()
       implicit none
              
       RETURN
      END SUBROUTINE commasmpi_startup

!-------------------------------------------------------------------------------

      SUBROUTINE commasmpi_bounds(nx,ny,nz)
       implicit none
          INTEGER, INTENT(IN) :: nx,ny,nz
       
          nxbeg = 1
          nxend = nx
          nybeg = 1
          nyend = ny
          nzbeg = 1
          nzend = nz

!          nx = nxt
!          ny = nyt

          nzt = nz
          itile = nxend
          jtile = nyend
          ktile = nzend
          
          ixbeg = 1
          jybeg = 1
          kzbeg = 1
          ixend = nxend
          jyend = nyend
          kzend = nzend
          
          nproci = 1
          nprocj = 1
          nprock = 1
          allocate ( procmap( nproci, nprocj, nprock ) )

       RETURN
       
      END SUBROUTINE commasmpi_bounds


!-------------------------------------------------------------------------------

      SUBROUTINE commasmpi_links(bcx,bcy)
         IMPLICIT NONE
         INTEGER bcx,bcy
       RETURN

      END SUBROUTINE commasmpi_links

!-------------------------------------------------------------------------------

      SUBROUTINE commasmpi_allocate(ng)

      implicit none
      
      integer ng

       RETURN

      END SUBROUTINE commasmpi_allocate

!-------------------------------------------------------------------------------

      SUBROUTINE commasmpi_shutdown ()

          IMPLICIT NONE
          
!          STOP

      END SUBROUTINE commasmpi_shutdown
!-------------------------------------------------------------------------------

      SUBROUTINE commasmpi_abort ()

          IMPLICIT NONE
          
          STOP

      END SUBROUTINE commasmpi_abort
!-------------------------------------------------------------------------------

#endif

END MODULE commasmpi_module

