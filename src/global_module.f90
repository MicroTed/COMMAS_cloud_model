!########################################################################
!########################################################################
!#########                                                      #########
!#########                  module rsa_table                    #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE rsa_table

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine read in the precalculated scattering amplitude
! from the tables, which are calcualted using T-matrix method.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 4/17/2008
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  SAVE

!-----------------------------------------------------------------------
! PARAMETER
! dsr : rain drop size,  rsa: scattering amplitude for rain drop
! dss : snow drop size,  ssa: scattering amplitude for snow aggregate
! dsh : hail drop size,  hsa: scattering amplitude for hailstone
! dsg : grpl drop size,  gsa: scattering amplitude for graupel  
!-----------------------------------------------------------------------
  INTEGER, PARAMETER :: nd = 112, ns = 8, nfw = 21
  INTEGER :: bin_opt = 0
!  INTEGER, PARAMETER :: nd_r = 100 ! 100 
!  INTEGER, PARAMETER :: nd_s = 112
!  INTEGER, PARAMETER :: nd_g = 625 ! 625 
!  INTEGER, PARAMETER :: nd_h = 875 ! 875 
  INTEGER :: nd_r = 100, nd_s = 112, nk_s = 1, nd_g = 625, nd_h = 875
  REAL, PRIVATE, ALLOCATABLE :: rsa(:,:)
  REAL, PRIVATE, ALLOCATABLE :: ssa(:,:,:)
  REAL, PRIVATE, ALLOCATABLE :: hsa(:,:,:)
  REAL, PRIVATE, ALLOCATABLE :: gsa(:,:,:)
  REAL, ALLOCATABLE :: dsr(:), dss(:), dsh(:), dsg(:)

  COMPLEX, ALLOCATABLE :: far_b(:), fbr_b(:), far_f(:), fbr_f(:)
  COMPLEX, ALLOCATABLE :: fas_b(:,:), fbs_b(:,:), fas_f(:,:), fbs_f(:,:)
  COMPLEX, ALLOCATABLE :: fah_b(:,:), fbh_b(:,:), fah_f(:,:), fbh_f(:,:)
  COMPLEX, ALLOCATABLE :: fag_b(:,:), fbg_b(:,:), fag_f(:,:), fbg_f(:,:)

  CONTAINS
  
  SUBROUTINE set_nbins(bin_opt_in)
!-----------------------------------------------------------------------
! Set number of bins for each category based on bin_opt
! bin_opt = 0: original bulk
!         = 1: Takahashi bins
!         = 2: PARSIVEL bins (not yet implemented)
!-----------------------------------------------------------------------
    INTEGER :: bin_opt_in
        
    IF(bin_opt_in == 0) THEN
      bin_opt = bin_opt_in
      nd_r = 100
      nd_s = 112
      nk_s = 1
      nd_g = 625
      nd_h = 875
    ELSEIF(bin_opt_in == 1) THEN
      bin_opt = bin_opt_in
      nd_r = 34
      nd_s = 21 ! Not used currently
      nk_s = 5  ! Not used currently
      nd_g = 45
      nd_h = 45
    ENDIF
    
  END SUBROUTINE set_nbins

  SUBROUTINE read_table (rsafndir)
!-----------------------------------------------------------------------
! Read radar scattering amplitudes tables calculated using T-matrix
! method.
!-----------------------------------------------------------------------

    INTEGER :: istatus, i, j, k
    !INTEGER, PARAMETER :: nfw = 21
    CHARACTER (LEN=256) :: rsafndir
    CHARACTER (LEN=256) :: rsafn
    CHARACTER (LEN=3), DIMENSION(nfw) :: extn = (/'000','005','010',    &
           '015','020','025','030','035','040','045','050','055','060', &
           '065','070','075','080','085','090','095','100'/)
    CHARACTER (LEN=256) :: head
    logical, save :: l_read_already = .false.

!   IF(bin_opt == 0) THEN ! Original bins for bulk schemes
!     nd_r = nd
!     nd_s = nd
!     nd_g = nd
!     nd_h = nd
!   ENDIF

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------
    
    IF ( l_read_already ) return
    
    ALLOCATE(dsr      (nd_r),stat=istatus)
    !CALL check_alloc_status(istatus, "dsr")
    ALLOCATE(dss      (nd_s),stat=istatus)
    !CALL check_alloc_status(istatus, "dss")
    ALLOCATE(dsh      (nd_h),stat=istatus)
    !CALL check_alloc_status(istatus, "dsh")
    ALLOCATE(dsg      (nd_g),stat=istatus)
    !CALL check_alloc_status(istatus, "dsg")
    ALLOCATE(rsa   (ns,nd_r),stat=istatus)
    !CALL check_alloc_status(istatus, "rsa")
    ALLOCATE(ssa(ns,nd_s,nfw),stat=istatus)
    !CALL check_alloc_status(istatus, "ssa")
    ALLOCATE(hsa(ns,nd_h,nfw),stat=istatus)
    !CALL check_alloc_status(istatus, "hsa")
    ALLOCATE(gsa(ns,nd_g,nfw),stat=istatus)
    !CALL check_alloc_status(istatus, "gsa")
    
    ALLOCATE(far_b(nd_r),stat=istatus)
    ALLOCATE(fbr_b(nd_r),stat=istatus)
    ALLOCATE(far_f(nd_r),stat=istatus)
    ALLOCATE(fbr_f(nd_r),stat=istatus)
    ALLOCATE(fas_b(nd_s,nfw),stat=istatus)
    ALLOCATE(fbs_b(nd_s,nfw),stat=istatus)
    ALLOCATE(fas_f(nd_s,nfw),stat=istatus)
    ALLOCATE(fbs_f(nd_s,nfw),stat=istatus)
    ALLOCATE(fah_b(nd_h,nfw),stat=istatus)
    ALLOCATE(fbh_b(nd_h,nfw),stat=istatus)
    ALLOCATE(fah_f(nd_h,nfw),stat=istatus)
    ALLOCATE(fbh_f(nd_h,nfw),stat=istatus)
    ALLOCATE(fag_b(nd_g,nfw),stat=istatus)
    ALLOCATE(fbg_b(nd_g,nfw),stat=istatus)
    ALLOCATE(fag_f(nd_g,nfw),stat=istatus)
    ALLOCATE(fbg_f(nd_g,nfw),stat=istatus)

!-----------------------------------------------------------------------
!  Read rain
!-----------------------------------------------------------------------
!   rsafn = TRIM(rsafndir)//'/SCTT_RAIN_fw100.dat'
    write (rsafn, '(A, A)') trim(rsafndir), '/SCTT_RAIN_fw100.dat'
    OPEN(UNIT=51,FILE=TRIM(rsafn),STATUS='old',FORM='formatted')
    READ(51,*) head

    DO j=1,nd_r
      IF(bin_opt == 0) THEN
        READ(51,'(f5.2,8e13.5)') dsr(j), (rsa(i,j),i=1,8)
      ELSEIF(bin_opt == 1) THEN
        READ(51,'(f10.6,8e13.5)') dsr(j), (rsa(i,j),i=1,8)
      ENDIF
    ENDDO
    CLOSE(51)

    far_b = CMPLX(rsa(1,:),rsa(2,:))
    fbr_b = CMPLX(rsa(3,:),rsa(4,:))
    far_f = CMPLX(rsa(5,:),rsa(6,:))
    fbr_f = CMPLX(rsa(7,:),rsa(8,:))

!-----------------------------------------------------------------------
!  Read snow
!-----------------------------------------------------------------------
    IF(bin_opt == 0) THEN ! Only read snow table for original bulk schemes for now
                          ! Eventually will add Takahashi to this
    DO k=1,nfw
!     rsafn = TRIM(rsafndir)//'/SCTT_SNOW_fw'//extn(k)//'.dat'
      write (rsafn, '(A, A, A, A)') trim(rsafndir), '/SCTT_SNOW_fw', extn(k), '.dat'
      OPEN(UNIT=51,FILE=TRIM(rsafn),STATUS='old',FORM='formatted')
      READ(51,*) head

      DO j=1,nd_s
        IF(bin_opt == 0) THEN
          READ(51,'(f5.2,8e13.5)') dss(j), (ssa(i,j,k),i=1,8)
        ELSEIF(bin_opt == 1) THEN
          READ(51,'(f10.6,8e13.5)') dss(j), (ssa(i,j,k),i=1,8)
        ENDIF
      ENDDO
      CLOSE(51)
    ENDDO

    fas_b = CMPLX(ssa(1,:,:),ssa(2,:,:))
    fbs_b = CMPLX(ssa(3,:,:),ssa(4,:,:))
    fas_f = CMPLX(ssa(5,:,:),ssa(6,:,:))
    fbs_f = CMPLX(ssa(7,:,:),ssa(8,:,:))
    
    ENDIF

!-----------------------------------------------------------------------
!  Read hail
!-----------------------------------------------------------------------
    DO k=1, nfw
!     rsafn = TRIM(rsafndir)//'/SCTT_HAIL_fw'//extn(k)//'.dat'
      write (rsafn, '(A, A, A, A)') trim(rsafndir), '/SCTT_HAIL_fw', extn(k), '.dat'
      OPEN(UNIT=51,FILE=TRIM(rsafn),STATUS='old',FORM='formatted')
      READ(51,*) head

      DO j=1,nd_h
        IF(bin_opt == 0) THEN
          READ(51,'(f5.2,8e13.5)') dsh(j), (hsa(i,j,k),i=1,8)
        ELSEIF(bin_opt == 1) THEN
          READ(51,'(f10.6,8e13.5)') dsh(j), (hsa(i,j,k),i=1,8)
        ENDIF
      ENDDO
      CLOSE(51)
    ENDDO

    fah_b = CMPLX(hsa(1,:,:),hsa(2,:,:))
    fbh_b = CMPLX(hsa(3,:,:),hsa(4,:,:))
    fah_f = CMPLX(hsa(5,:,:),hsa(6,:,:))
    fbh_f = CMPLX(hsa(7,:,:),hsa(8,:,:))

!-----------------------------------------------------------------------
!  Read graupel
!-----------------------------------------------------------------------
    DO k=1, nfw
!     rsafn = TRIM(rsafndir)//'/SCTT_GRPL_fw'//extn(k)//'.dat'
      write (rsafn, '(A, A, A, A)') trim(rsafndir), '/SCTT_GRPL_fw', extn(k), '.dat'
      OPEN(UNIT=51,FILE=TRIM(rsafn),STATUS='old',FORM='formatted')
      READ(51,*) head

      DO j=1,nd_g
        IF(bin_opt == 0) THEN
          READ(51,'(f5.2,8e13.5)') dsg(j), (gsa(i,j,k),i=1,8)
        ELSEIF(bin_opt == 1) THEN
          READ(51,'(f10.6,8e13.5)') dsg(j), (gsa(i,j,k),i=1,8)
        ENDIF
      ENDDO
      CLOSE(51)
    ENDDO

    fag_b = CMPLX(gsa(1,:,:),gsa(2,:,:))
    fbg_b = CMPLX(gsa(3,:,:),gsa(4,:,:))
    fag_f = CMPLX(gsa(5,:,:),gsa(6,:,:))
    fbg_f = CMPLX(gsa(7,:,:),gsa(8,:,:))

    l_read_already = .true.
    
    deallocate(rsa,ssa,hsa,gsa)

  END SUBROUTINE read_table

END MODULE rsa_table

