#ifndef RKIND
#define RKIND 4
#endif

#if RKIND == 4
#define FLOAT_PRECISION MPI_REAL
#else
#define FLOAT_PRECISION MPI_DOUBLE_PRECISION
#endif

#if defined (MPI) && defined (BOXMG)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mgsetupdg(nx,ny,nz,nor,
     :   dx,dy,dz,gxt,gyt,gzt,
     :   dxx,dyy,dzz,
     :   nnxs,nnys,nnzs,
     :   nbw,nbe,nbs,nbn,
     :   istretch,
     :   iunit,bcx,bcy)
      
      USE ELEC_MODULE
      USE COMMASMPI_MODULE
      USE GRID_MODULE
      USE INDEX_MODULE

      IMPLICIT NONE

! ------------------------------------------------
!     Includes
! ------------------------------------------------

      INCLUDE   'mpif.h'

#include        "BMG_constants_f77.h"
      INCLUDE   'BMG_workspace_f77.h'
      INCLUDE   'BMG_parameters_f77.h'

! ------------------------------------------------
!      Multigrid/Workspace Memory Allocation: 
! ------------------------------------------------

      integer nx,ny,nz,nor

 
      real dzz(nz),dxx(nx),dyy(ny)    

      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4)
      real  gzt(-nor+1:nz+nor,4)

      integer nnxs(2),nnys(2),nnzs(2)
      integer nbw,nbe,nbs,nbn
      real dx,dy,dz
      logical istretch
      integer iunit
      integer iestag
      logical lcalce   ! whether to calculate electric field components when iestag=1
      
      
      integer bcx, bcy
      
      integer i

      !
      ! Workspace pointers and logicals

      INTEGER  BMG_MSG_iGRIDdum(1)

      REAL*8    T1, T2, T_SOLVE_LOOP

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 
      IF ( NBMG_iPARMSdg .ne. NBMG_iPARMS .or.
     :     NBMG_rPARMSdg .ne. NBMG_rPARMS .or.
     :     NBMG_IOFLAGdg .ne. NBMG_IOFLAG .or.
     :     NBMG_pWORKdg  .ne. NBMG_pWORK  .or. 
     :     NBMG_InWORKdg .ne. NBMG_InWORK .or. 
     :     NBMG_MSG_pGRIDdg .ne. NBMG_MSG_pGRID) THEN
     
        write(0,*) 'Problem with BoxMG parameter constants! Fix in elecmodule!'
        STOP
      ENDIF
      
      BMG_iPARMSdg(id_BMG3_POINTERS) = BMG_USE_pointers

c #ifdef __apple__
c      BMG_iPARMSdg(id_BMG3_CG_COMM) = BMG_CG_ALLGATHER
c #endif
      CALL commas_SETUP_ProcGrid_file(                           
     :                 BMG_MSG_iGRIDdum, BMG_MSG_pGRIDdg,       
     :                 BMG_iPARMSdg, NBMG_MSG_iGRIDmdg, NBMG_MSG_iGRIDdg,
     :                 nbw,nbe,nbs,nbn,nnxs(1),nnys(1),nnzs(1),nx,ny,nz,1,
     :                 my_rank, number_of_processes, my_comm, mpi_error_code  
     :                 )

      allocate(BMG_MSG_iGRIDdg(1:NBMG_MSG_iGRIDdg))

      BMG_iPARMSdg(id_BMG3_POINTERS) = BMG_NO_pointers
      NBMG_MSG_iGRIDmdg = NBMG_MSG_iGRIDdg


      CALL commas_SETUP_ProcGrid_file(                           
     :                 BMG_MSG_iGRIDdg, BMG_MSG_pGRIDdg,       
     :                 BMG_iPARMSdg, NBMG_MSG_iGRIDmdg, NBMG_MSG_iGRIDdg,
     :                 nbw,nbe,nbs,nbn,nnxs(1),nnys(1),nnzs(1),nx,ny,nz,1,
     :                 my_rank, number_of_processes, my_comm, mpi_error_code  
     :                 )


      NGxdg = BMG_MSG_iGRIDdg(id_BMG_MSG_NGx)
      NGydg = BMG_MSG_iGRIDdg(id_BMG_MSG_NGy)
      NGzdg = BMG_MSG_iGRIDdg(id_BMG_MSG_NGz)

      NLxdg = BMG_MSG_iGRIDdg(id_BMG_MSG_NLx)
      NLydg = BMG_MSG_iGRIDdg(id_BMG_MSG_NLy)
      NLzdg = BMG_MSG_iGRIDdg(id_BMG_MSG_NLz)

      iGsdg = BMG_MSG_iGRIDdg(id_BMG_MSG_iGs)
      jGsdg = BMG_MSG_iGRIDdg(id_BMG_MSG_jGs)
      kGsdg = BMG_MSG_iGRIDdg(id_BMG_MSG_kGs)
      
      iGsdg0 = iGsdg - 1
      jGsdg0 = jGsdg - 1
      kGsdg0 = kGsdg - 1

      IF ( MyProcI /= BMG_MSG_iGRIDdg(id_BMG_MSG_MyProcI) .or.
     :      MyProcJ /= BMG_MSG_iGRIDdg(id_BMG_MSG_MyProcJ) .or.
     :      MyProcK /= BMG_MSG_iGRIDdg(id_BMG_MSG_MyProcK) ) THEN
           write(0,*) 'problem with myproci,j, or k! ', my_rank
           write(0,*) my_rank, MyProcI, BMG_MSG_iGRIDdg(id_BMG_MSG_MyProcI) 
           write(0,*) my_rank, MyProcJ, BMG_MSG_iGRIDdg(id_BMG_MSG_MyProcJ)
           write(0,*) my_rank, MyProcK, BMG_MSG_iGRIDdg(id_BMG_MSG_MyProcK)
       ENDIF

!      IF ( my_rank.EQ.0 ) THEN 
        write(iunit,*) 'Global grid dimensions are (nx,ny,nz) = ',   
     &       NGxdg,NGydg,NGzdg
        write(iunit,*) 'My local dims are NLxdg,NLydg,NLzdg = ',NLxdg,NLydg,NLzdg
        write(iunit,*) 'iGsdg,jGsdg,kGsdg = ',iGsdg,jGsdg,kGsdg
        write(iunit,*) 'nx,ny,nz,nnxs,nnys = ',nx,ny,nz,nnxs(1),nnys(1)
!      ENDIF

      CALL EX_SETUP_BMG_parms(                                      
     :               BMG_iPARMSdg, BMG_rPARMSdg, BMG_IOFLAGdg, 
     :               my_rank, number_of_processes, my_comm, mpi_error_code   
     :               )

      BMG_IOFLAGdg(iBMG3_BUG_RES_RELAX)    =  .FALSE.
      BMG_IOFLAGdg(iBMG3_BUG_RES_CG_SOLVE) =  .FALSE.
      BMG_IOFLAGdg(iBMG3_OUT_TIME_SETUP)   =  .FALSE.
      BMG_IOFLAGdg(iBMG3_OUT_TIME_CYCLING) =  .FALSE.
      BMG_IOFLAGdg(iBMG3_OUT_TIME_TOTAL)   =  .FALSE.

c separable fine GS speed-up doesn't seem to work with stretching for some reason,
c so make sure to use default if stretching is active.
c      IF ( istretch ) BMG_iPARMSdg(id_BMG3_SEPARABLE) = 0

! ------------------------------
!      Workspace Allocation:
! ------------------------------

      pSIdg=1
      pSRdg=1

      !
      ! Start with nothing stored in the workspace array
      !
      DO i=1, NBMG_InWORKdg
         BMG_InWORKdg(i) = .FALSE.
      ENDDO

      BMG_InWORKdg(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      BMG_InWORKdg(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      BMG_InWORKdg(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      BMG_InWORKdg(i_InWORK_RES) = .TRUE.     ! store RES in work array

      BMG_iPARMSdg(id_BMG3_POINTERS) = BMG_USE_pointers

      CALL BMG3_SymStd_SETUP_PtrWork(                     
     :          NLxdg, NLydg, NLzdg, NGxdg, NGydg, NGzdg,             
     :          iGsdg, jGsdg, kGsdg,                            
     :          BMG_iPARMSdg, BMG_rPARMSdg, BMG_IOFLAGdg,       
     :          my_rank, number_of_processes, NProcI, NProcJ, NProcK,    
     :          my_comm,                           
     :          NOGmdg, NFmdg, NSOmdg,                          
     :          NBMG_iWORKdg, NBMG_rWORKdg,                   
     :          NBMG_iWORK_PLdg, NBMG_rWORK_PLdg,             
     :          NBMG_iWORK_CSdg, NBMG_rWORK_CSdg,             
     :          BMG_pWORKdg, BMG_InWORKdg, pSRdg, pSIdg           
     :          )

      NOGdg  = BMG_iPARMSdg(id_BMG3_DIM_NOG)
      NFdg   = BMG_iPARMSdg(id_BMG3_DIM_NF)
      NCbmgdg   = BMG_iPARMSdg(id_BMG3_DIM_NC)
      NCIdg  = BMG_iPARMSdg(id_BMG3_DIM_NCI)
      NSOdg  = BMG_iPARMSdg(id_BMG3_DIM_NSO)
      NSORdg = BMG_iPARMSdg(id_BMG3_DIM_NSOR)
      NCBWdg = BMG_iPARMSdg(id_BMG3_DIM_NCBW)
      NCUdg  = BMG_iPARMSdg(id_BMG3_DIM_NCU)

      NMSGidg = BMG_iPARMSdg(id_BMG3_DIM_NMSGi)
      NMSGrdg = BMG_iPARMSdg(id_BMG3_DIM_NMSGr)

! ------------------------------
!     Report space requirements
! ------------------------------

      IF (my_rank .eq. 0) THEN
         WRITE(iunit,*)
         WRITE(iunit,*) 'Memory allocation in MB'
         WRITE(iunit,*) '======================='
         WRITE(iunit,*) 'size of real = ',RKIND
         WRITE(iunit,1210) 'SO           : ',NSOdg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'Q            : ',NFdg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'QF           : ',NFdg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'BMG_rWORK    : ',NBMG_rWORKdg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'BMG_iWORK    : ',NBMG_iWORKdg*4./(1024*1024)
         WRITE(iunit,1210) 'pMSG         : ',NBMG_pMSG*NOGdg*4./(1024*1024)
         WRITE(iunit,1210) 'pMSGSO       : ',NBMG_pMSG*NOGdg*4./(1024*1024)
         WRITE(iunit,1210) 'pLS          : ',NBMG_pLS*NOGdg*4./(1024*1024)
         WRITE(iunit,1210) 'BMG_iWORK_PL : ',NBMG_iWORK_PLdg*4./(1024*1024)
         WRITE(iunit,1210) 'BMG_rWORK_PL : ',NBMG_rWORK_PLdg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'BMG_iWORK_CS : ',NBMG_iWORK_CSdg*4./(1024*1024)
         WRITE(iunit,1210) 'BMG_rWORK_CS : ',NBMG_rWORK_CSdg*float(RKIND)/(1024*1024)
         WRITE(iunit,*)
       ENDIF

 1210 FORMAT(4X,A,1X,F18.6,1X,'MBytes')

! ------------------------------
!     Allocate the space
! ------------------------------

       allocate(SOdg(1:NSOdg))
       allocate(Qdg(1:NFdg))
       allocate(QFdg(1:NFdg))
       allocate(BMG_rWORKdg(1:NBMG_rWORKdg))
       allocate(BMG_iWORKdg(1:NBMG_iWORKdg))
       allocate(BMG_iWORK_PLdg(1:NBMG_iWORK_PLdg))
       allocate(BMG_rWORK_PLdg(1:NBMG_rWORK_PLdg))
       allocate(BMG_iWORK_CSdg(1:NBMG_iWORK_CSdg))
       allocate(BMG_rWORK_CSdg(1:NBMG_rWORK_CSdg))
       allocate(pMSGdg(1:NBMG_pMSG*NOGdg))
       allocate(pMSGSOdg(1:NBMG_pMSG*NOGdg))

         ! ------------------------------
         !     Zero the Stencil:
         ! ------------------------------
      
         DO i=1, NSOdg
            SOdg(i)=rZERO
         END DO

         DO i=1, NFdg
            Qdg(i)=rZERO
            QFdg(i)=rZERO
         END DO

         ! ------------------------------
         !     Zero the workspace:
         ! ------------------------------

         DO i=1, NBMG_rWORKdg
            BMG_rWORKdg(i)=rZERO
         END DO

         DO i=1, NBMG_iWORKdg
            BMG_iWORKdg(i)=iZERO
         END DO

         DO i=1, NBMG_rWORK_PLdg
            BMG_rWORK_PLdg(i)=rZERO
         END DO

         DO i=1, NBMG_iWORK_PLdg
            BMG_iWORK_PLdg(i)=iZERO
         END DO

         DO i=1, NBMG_rWORK_CSdg
            BMG_rWORK_CSdg(i)=rZERO
         END DO

         DO i=1, NBMG_iWORK_CSdg
            BMG_iWORK_CSdg(i)=iZERO
         END DO

      
      RETURN
      END

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mgsetuplg(nx,ny,nz,nor,
     :   dx,dy,dz,gxt,gyt,gzt,
     :   dxx,dyy,dzz,
     :   nnxs,nnys,nnzs,
     :   nbw,nbe,nbs,nbn,
     :   istretch,
     :   iunit,bcx,bcy)
      
      USE ELEC_MODULE
      USE COMMASMPI_MODULE
      USE GRID_MODULE
      USE INDEX_MODULE

      IMPLICIT NONE

! ------------------------------------------------
!     Includes
! ------------------------------------------------

      INCLUDE   'mpif.h'

#include        "BMG_constants_f77.h"
      INCLUDE   'BMG_workspace_f77.h'
      INCLUDE   'BMG_parameters_f77.h'

! ------------------------------------------------
!      Multigrid/Workspace Memory Allocation: 
! ------------------------------------------------

      integer nx,ny,nz,nor

 
      real dzz(nz),dxx(nx),dyy(ny)    

      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4)
      real  gzt(-nor+1:nz+nor,4)

      integer nnxs(2),nnys(2),nnzs(2)
      integer nbw,nbe,nbs,nbn
      real dx,dy,dz
      logical istretch
      integer iunit
      integer iestag
      
      
      integer bcx, bcy
      
      integer i

      !
      ! Workspace pointers and logicals

      INTEGER  BMG_MSG_iGRIDdum(1)

      REAL*8    T1, T2, T_SOLVE_LOOP

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 
      
      BMG_iPARMSlg(id_BMG3_POINTERS) = BMG_USE_pointers

c #ifdef __apple__
c      BMG_iPARMSlg(id_BMG3_CG_COMM) = BMG_CG_ALLGATHER
c #endif

      CALL commas_SETUP_ProcGrid_file(                           
     :                 BMG_MSG_iGRIDdum, BMG_MSG_pGRIDlg,       
     :                 BMG_iPARMSlg, NBMG_MSG_iGRIDmlg, NBMG_MSG_iGRIDlg,
     :                 nbw,nbe,nbs,nbn,nnxs(2),nnys(2),nnzs(2),nx,ny,nz,lightintx,
     :                 my_rank, number_of_processes, my_comm, mpi_error_code  
     :                 )

      allocate(BMG_MSG_iGRIDlg(1:NBMG_MSG_iGRIDlg))

      BMG_iPARMSlg(id_BMG3_POINTERS) = BMG_NO_pointers
      NBMG_MSG_iGRIDmlg = NBMG_MSG_iGRIDlg


      CALL commas_SETUP_ProcGrid_file(                           
     :                 BMG_MSG_iGRIDlg, BMG_MSG_pGRIDlg,       
     :                 BMG_iPARMSlg, NBMG_MSG_iGRIDmlg, NBMG_MSG_iGRIDlg,
     :                 nbw,nbe,nbs,nbn,nnxs(2),nnys(2),nnzs(2),nx,ny,nz,lightintx,
     :                 my_rank, number_of_processes, my_comm, mpi_error_code  
     :                 )


      NGxlg = BMG_MSG_iGRIDlg(id_BMG_MSG_NGx)
      NGylg = BMG_MSG_iGRIDlg(id_BMG_MSG_NGy)
      NGzlg = BMG_MSG_iGRIDlg(id_BMG_MSG_NGz)

      NLxlg = BMG_MSG_iGRIDlg(id_BMG_MSG_NLx)
      NLylg = BMG_MSG_iGRIDlg(id_BMG_MSG_NLy)
      NLzlg = BMG_MSG_iGRIDlg(id_BMG_MSG_NLz)

      iGslg = BMG_MSG_iGRIDlg(id_BMG_MSG_iGs)
      jGslg = BMG_MSG_iGRIDlg(id_BMG_MSG_jGs)
      kGslg = BMG_MSG_iGRIDlg(id_BMG_MSG_kGs)
      
      iGslg0 = iGslg - 1
      jGslg0 = jGslg - 1
      kGslg0 = kGslg - 1
      

      IF ( MyProcI /= BMG_MSG_iGRIDlg(id_BMG_MSG_MyProcI) .or.
     :      MyProcJ /= BMG_MSG_iGRIDlg(id_BMG_MSG_MyProcJ) .or.
     :      MyProcK /= BMG_MSG_iGRIDlg(id_BMG_MSG_MyProcK) ) THEN
           write(0,*) 'problem with myproci,j, or k! ', my_rank
           write(0,*) my_rank, MyProcI, BMG_MSG_iGRIDlg(id_BMG_MSG_MyProcI) 
           write(0,*) my_rank, MyProcJ, BMG_MSG_iGRIDlg(id_BMG_MSG_MyProcJ)
           write(0,*) my_rank, MyProcK, BMG_MSG_iGRIDlg(id_BMG_MSG_MyProcK)
       ENDIF

!      IF ( my_rank.EQ.0 ) THEN 
        write(iunit,*) 'Global lightning grid dimensions are (nx,ny,nz) = ',   
     &       NGxlg,NGylg,NGzlg
        write(iunit,*) 'lightintx = ',lightintx
        write(iunit,*) 'My local lightning dims are NLxlg,NLylg,NLzlg = ',NLxlg,NLylg,NLzlg
        write(iunit,*) 'My buffer zones are nbw,nbe,nbs,nbn = ',nbw,nbe,nbs,nbn
        write(iunit,*) 'nnxs(2),nnys(2),nnzs(2) = ',nnxs(2),nnys(2),nnzs(2)
        write(iunit,*) 'nx,ny = ',nx,ny
        write(iunit,*) 'iGslg,jGslg,kGslg = ',iGslg,jGslg,kGslg
        write(iunit,*) 'iGslg0,jGslg0,kGslg0 = ',iGslg0,jGslg0,kGslg0
!      ENDIF

      CALL EX_SETUP_BMG_parms(                                      
     :               BMG_iPARMSlg, BMG_rPARMSlg, BMG_IOFLAGlg, 
     :               my_rank, number_of_processes, my_comm, mpi_error_code   
     :               )

      BMG_IOFLAGlg(iBMG3_BUG_RES_RELAX)    =  .FALSE.
      BMG_IOFLAGlg(iBMG3_BUG_RES_CG_SOLVE) =  .FALSE.
      BMG_IOFLAGlg(iBMG3_OUT_TIME_SETUP)   =  .FALSE.
      BMG_IOFLAGlg(iBMG3_OUT_TIME_CYCLING) =  .FALSE.
      BMG_IOFLAGlg(iBMG3_OUT_TIME_TOTAL)   =  .FALSE.

c separable fine GS speed-up doesn't seem to work with stretching for some reason,
c so make sure to use default if stretching is active.
c      IF ( istretch ) BMG_iPARMSlg(id_BMG3_SEPARABLE) = 0

! ------------------------------
!      Workspace Allocation:
! ------------------------------

      pSIlg=1
      pSRlg=1

      !
      ! Start with nothing stored in the workspace array
      !
      DO i=1, NBMG_InWORKdg
         BMG_InWORKlg(i) = .FALSE.
      ENDDO

      BMG_InWORKlg(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      BMG_InWORKlg(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      BMG_InWORKlg(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      BMG_InWORKlg(i_InWORK_RES) = .TRUE.     ! store RES in work array

      BMG_iPARMSlg(id_BMG3_POINTERS) = BMG_USE_pointers

      CALL BMG3_SymStd_SETUP_PtrWork(                     
     :          NLxlg, NLylg, NLzlg, NGxlg, NGylg, NGzlg,             
     :          iGslg, jGslg, kGslg,                            
     :          BMG_iPARMSlg, BMG_rPARMSlg, BMG_IOFLAGlg,       
     :          my_rank, number_of_processes, NProcI, NProcJ, NProcK,    
     :          my_comm,                           
     :          NOGmlg, NFmlg, NSOmlg,                          
     :          NBMG_iWORKlg, NBMG_rWORKlg,                   
     :          NBMG_iWORK_PLlg, NBMG_rWORK_PLlg,             
     :          NBMG_iWORK_CSlg, NBMG_rWORK_CSlg,             
     :          BMG_pWORKlg, BMG_InWORKlg, pSRlg, pSIlg           
     :          )

      NOGlg  = BMG_iPARMSlg(id_BMG3_DIM_NOG)
      NFlg   = BMG_iPARMSlg(id_BMG3_DIM_NF)
      NCbmglg = BMG_iPARMSlg(id_BMG3_DIM_NC)
      NCIlg  = BMG_iPARMSlg(id_BMG3_DIM_NCI)
      NSOlg  = BMG_iPARMSlg(id_BMG3_DIM_NSO)
      NSORlg = BMG_iPARMSlg(id_BMG3_DIM_NSOR)
      NCBWlg = BMG_iPARMSlg(id_BMG3_DIM_NCBW)
      NCUlg  = BMG_iPARMSlg(id_BMG3_DIM_NCU)

      NMSGilg = BMG_iPARMSlg(id_BMG3_DIM_NMSGi)
      NMSGrlg = BMG_iPARMSlg(id_BMG3_DIM_NMSGr)

! ------------------------------
!     Report space requirements
! ------------------------------

      IF (my_rank .eq. 0) THEN
         WRITE(iunit,*)
         WRITE(iunit,*) 'Lightning Memory allocation in MB'
         WRITE(iunit,*) '======================='
         WRITE(iunit,*) 'size of real = ',RKIND
         WRITE(iunit,1210) 'SO           : ',NSOlg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'Q            : ',NFlg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'QF           : ',NFlg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'BMG_rWORK    : ',NBMG_rWORKlg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'BMG_iWORK    : ',NBMG_iWORKlg*4./(1024*1024)
         WRITE(iunit,1210) 'pMSG         : ',NBMG_pMSG*NOGlg*4./(1024*1024)
         WRITE(iunit,1210) 'pMSGSO       : ',NBMG_pMSG*NOGlg*4./(1024*1024)
         WRITE(iunit,1210) 'pLS          : ',NBMG_pLS*NOGlg*4./(1024*1024)
         WRITE(iunit,1210) 'BMG_iWORK_PL : ',NBMG_iWORK_PLlg*4./(1024*1024)
         WRITE(iunit,1210) 'BMG_rWORK_PL : ',NBMG_rWORK_PLlg*float(RKIND)/(1024*1024)
         WRITE(iunit,1210) 'BMG_iWORK_CS : ',NBMG_iWORK_CSlg*4./(1024*1024)
         WRITE(iunit,1210) 'BMG_rWORK_CS : ',NBMG_rWORK_CSlg*float(RKIND)/(1024*1024)
         WRITE(iunit,*)
       ENDIF

 1210 FORMAT(4X,A,1X,F18.6,1X,'MBytes')

! ------------------------------
!     Allocate the space
! ------------------------------

       allocate(SOlg(1:NSOlg))
       allocate(Qlg(1:NFlg))
       allocate(QFlg(1:NFlg))
       allocate(BMG_rWORKlg(1:NBMG_rWORKlg))
       allocate(BMG_iWORKlg(1:NBMG_iWORKlg))
       allocate(BMG_iWORK_PLlg(1:NBMG_iWORK_PLlg))
       allocate(BMG_rWORK_PLlg(1:NBMG_rWORK_PLlg))
       allocate(BMG_iWORK_CSlg(1:NBMG_iWORK_CSlg))
       allocate(BMG_rWORK_CSlg(1:NBMG_rWORK_CSlg))
       allocate(pMSGlg(1:NBMG_pMSG*NOGlg))
       allocate(pMSGSOlg(1:NBMG_pMSG*NOGlg))

         ! ------------------------------
         !     Zero the Stencil:
         ! ------------------------------
      
         DO i=1, NSOlg
            SOlg(i)=rZERO
         END DO

         DO i=1, NFlg
            Qlg(i)=rZERO
            QFlg(i)=rZERO
         END DO

         ! ------------------------------
         !     Zero the workspace:
         ! ------------------------------

         DO i=1, NBMG_rWORKlg
            BMG_rWORKlg(i)=rZERO
         END DO

         DO i=1, NBMG_iWORKlg
            BMG_iWORKlg(i)=iZERO
         END DO

         DO i=1, NBMG_rWORK_PLlg
            BMG_rWORK_PLlg(i)=rZERO
         END DO

         DO i=1, NBMG_iWORK_PLlg
            BMG_iWORK_PLlg(i)=iZERO
         END DO

         DO i=1, NBMG_rWORK_CSlg
            BMG_rWORK_CSlg(i)=rZERO
         END DO

         DO i=1, NBMG_iWORK_CSlg
            BMG_iWORK_CSlg(i)=iZERO
         END DO

      
      RETURN
      END

! ###############################################################
! ###############################################################
! ###############################################################

      SUBROUTINE commas_SETUP_ProcGrid_file(
     &                BMG_MSG_iGRID, BMG_MSG_pGRID,
     &                BMG_iPARMS, NBMG_MSG_iGRIDm, NBMG_MSG_iGRID,
     &                nbw,nbe,nbs,nbn,nnxs,nnys,nnzs,NLx,NLy,NLz,ifac,
     &                MPI_MyProc, NProc, BMG_MPI_COMM, MPI_IERROR
     &                ) 

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_SETUP_ProcGrid_file reads and communicates minimal
C     data from a file that describes the processor grid and creates
C     related information for use locally.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================
      
      USE COMMASMPI_MODULE

      IMPLICIT NONE
      
C -----------------------------
C     Includes
C
      INCLUDE 'mpif.h'

      INCLUDE   'BMG_workspace_f77.h'
      INCLUDE   'BMG_parameters_f77.h'

C -----------------------------
C     Argument Declarations
C
      INTEGER  BMG_iPARMS(NBMG_iPARMS)

      INTEGER  NBMG_MSG_iGRIDm, NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRIDm),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER   BMG_MPI_COMM, MPI_MyProc, MPI_IERROR, NProc
      
      integer  nbw,nbe,nbs,nbn,nnxs,nnys,nnzs,ifac

C -----------------------------
C     Local Declarations
C 
      INTEGER  CE, GRIDDIST, i, iGs, j, jGs, k, kGs,
     &         MyProc, ! MyProcI, MyProcJ, MyProcK,
     &         NGx, NGy, NGz, NLx, NLy, NLz, 
     &         pSI, ! NProcI, NProcJ, NProcK, 
     &         p_iGs, p_jGs, p_kGs, p_NLx, p_NLy, p_NLz, 
     &         p_PC, p_PG, p_ProcCoord, p_ProcGrid

      LOGICAL  exist

C ==========================================================================

C ------------------------------------------- 
C     MPI -> MSG indexing:
C ------------------------------------------- 

      MyProc = MPI_MyProc+1

      GRIDDIST = BMG_MSG_GRIDDIST_CUSTOM

C ------------------------------------------- 
C     Open the data file, read the header
C -------------------------------------------     

      !
      !  Master:
      !
      IF ( MyProc .EQ. 1 ) THEN
         
         !
         ! Eliminate trailing blanks
         !

      ENDIF
         
C --------------------------------------------
C     Process and distribute the grid data
C --------------------------------------------

      IF ( GRIDDIST.EQ.BMG_MSG_GRIDDIST_AUTOMATIC ) THEN

         !
         !  This section would take NGx, NGy, NGz and
         !  distribute them as evenly as possible on 
         !  the NProcI x NProcJ x NProcK processor grid.
         ! 
         
      
      ELSEIF ( GRIDDIST.EQ.BMG_MSG_GRIDDIST_CUSTOM
     &        .OR. GRIDDIST.EQ.BMG_MSG_GRIDDIST_UNIFORM ) THEN
         !
         !  Master:
         !
         IF ( MyProc.EQ.1 ) THEN
            ! 
            ! Read the number of processors in i and j
            !
!            READ(10,*) NProcI, NProcJ, NProcK
            !
         ENDIF

         !
         !  All:
         !
!         CALL MPI_Bcast( NProcI, 1, MPI_INTEGER, 
!     &                   0, BMG_MPI_COMM, MPI_IERROR )
!         
!         CALL MPI_BCast( NProcJ, 1, MPI_INTEGER, 
!     &                   0, BMG_MPI_COMM, MPI_IERROR )
!
!         CALL MPI_BCast( NProcK, 1, MPI_INTEGER, 
!     &                   0, BMG_MPI_COMM, MPI_IERROR )
         
!         IF ( NProcI*NProcJ*NProcK .NE. NProc ) THEN
!            WRITE(*,*) 'Procgrid Error : ',
!     &                 NProcI,' x ', NProcJ, 'x', NProcK, ' != ', NProc
!            STOP
!         ENDIF

         !
         ! Compute storage requirements and pointers
         !
         pSI = NBMG_MSG_iGRID_data + 1
         !
         BMG_MSG_pGRID(ip_BMG_MSG_NLx_Grid)  = pSI
         pSI = pSI + NProcI
         !
         BMG_MSG_pGRID(ip_BMG_MSG_NLy_Grid)  = pSI
         pSI = pSI + NProcJ
         !
         BMG_MSG_pGRID(ip_BMG_MSG_NLz_Grid)  = pSI
         pSI = pSI + NProcK
         !
         BMG_MSG_pGRID(ip_BMG_MSG_iGs_Grid)  = pSI
         pSI = pSI + NProcI
         !
         BMG_MSG_pGRID(ip_BMG_MSG_jGs_Grid)  = pSI
         pSI = pSI + NProcJ
         !
         BMG_MSG_pGRID(ip_BMG_MSG_kGs_Grid)  = pSI
         pSI = pSI + NProcK
         !
         BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)  = pSI
         pSI = pSI + NProcI*NProcJ*NProcK
         !
         BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord) = pSI
         pSI = pSI + 3*NProcI*NProcJ*NProcK
         !
         NBMG_MSG_iGRID = pSI-1

         IF ( BMG_iPARMS(id_BMG3_POINTERS).EQ.BMG_USE_pointers ) THEN
            
            NBMG_MSG_iGRIDm = NBMG_MSG_iGRID

!            CLOSE(10)
            RETURN

         ELSEIF ( BMG_iPARMS(id_BMG3_POINTERS).EQ.BMG_NO_pointers ) THEN
            
            IF ( NBMG_MSG_iGRID.GT.NBMG_MSG_iGRIDm ) THEN 

               IF (MyProc.EQ.1) THEN
                  WRITE(*,500) 'Space for BMG_MSG_iGRID is too small!'
                  WRITE(*,510) 'HAVE: NBMG_MSG_iGRIDm = ',
     &                         NBMG_MSG_iGRIDm
                  WRITE(*,520) 'NEED: NBMG_MSG_iGRIDm = ',
     &                         NBMG_MSG_iGRID
               ENDIF

               CALL MPI_FINALIZE(MPI_IERROR)
               STOP

            ENDIF

         ELSE

            IF (MyProc.EQ.1) THEN
               WRITE(*,500) 'Value of BMG_iPARMS(id_BMG3_pointers) is '
     &                      //'out of range ...'
               WRITE(*,510) 'HAVE: BMG_iPARMS(id_BMG3_POINTERS) = ',
     &                      BMG_iPARMS(id_BMG3_POINTERS)
               WRITE(*,505) 'NEED: BMG_iPARMS(id_BMG3_POINTERS) = '
     &                      //'0 or 1'  
            ENDIF

            CALL MPI_FINALIZE(MPI_IERROR)
            STOP

         ENDIF

         !
         !  Make local copies of pointers
         !
         p_NLx = BMG_MSG_pGRID(ip_BMG_MSG_NLx_Grid)
         p_NLy = BMG_MSG_pGRID(ip_BMG_MSG_NLy_Grid)
         p_NLz = BMG_MSG_pGRID(ip_BMG_MSG_NLz_Grid)
         p_iGs = BMG_MSG_pGRID(ip_BMG_MSG_iGs_Grid)
         p_jGs = BMG_MSG_pGRID(ip_BMG_MSG_jGs_Grid)
         p_kGs = BMG_MSG_pGRID(ip_BMG_MSG_kGs_Grid)
         !
         p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
         p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

         !
         !  Master:
         !
         IF ( MyProc.EQ.1 ) THEN
            !
            IF ( GRIDDIST.EQ.BMG_MSG_GRIDDIST_CUSTOM ) THEN
               !
               ! Read arrays of fine-grid local dimensions 
               !
               IF ( NProcI == 1 ) THEN
                 BMG_MSG_iGRID(p_NLx) = nnxs
               ELSE
               DO i=1, NProcI
                  IF ( i == 1 ) THEN
                  BMG_MSG_iGRID(p_NLx+i-1) = NLx + nbw
                  ELSEIF ( i == NProcI ) THEN
                  BMG_MSG_iGRID(p_NLx+i-1) = NLx + nbw + ifac*(nxend - nproci*nxt) ! this is for the global size, so use nbw since nbe might be zero
                                                         ! the last term, ifac*(nxend - nproci*nxt) accounts for possibility of
                                                         ! having an extra point in the east boundary tiles
                  ELSE
                  BMG_MSG_iGRID(p_NLx+i-1) = NLx
                  ENDIF
               ENDDO
               ENDIF
               !
               
               IF ( NProcJ == 1 ) THEN
                 BMG_MSG_iGRID(p_NLy) = nnys
               ELSE
               
               DO j=1, NProcJ
                 IF ( j == 1 ) THEN
                  BMG_MSG_iGRID(p_NLy+j-1) = NLy + nbs
                 
                 ELSEIF ( j == NProcJ ) THEN
                  BMG_MSG_iGRID(p_NLy+j-1) = NLy + nbs + ifac*(nyend - nprocj*nyt) ! this is for the global size, so use nbs since nbn might be zero
                                                         ! the last term, ifac*(nyend - nprocj*nyt) accounts for possibility of
                                                         ! having an extra point in the north boundary tiles
                 
                 ELSE
                  BMG_MSG_iGRID(p_NLy+j-1) = NLy
                 ENDIF
               ENDDO
               ENDIF
               !
               
               
               IF ( NProcK == 1 ) THEN
                 BMG_MSG_iGRID(p_NLz) = nnzs - 1
               ELSE
               
               DO k=1, NProcK
                 IF ( k == 1 ) THEN
                  BMG_MSG_iGRID(p_NLz+k-1) = NLz+1 ! does there need to be an extra point below ground??
                 
                 ELSEIF ( k == NProcK ) THEN
                  BMG_MSG_iGRID(p_NLz+k-1) = NLz + nbs + ifac*(nzend - nprock*nzt) ! this is for the global size, so use nbs since nbn might be zero
                                                         ! the last term, ifac*(nyend - nprocj*nyt) accounts for possibility of
                                                         ! having an extra point in the north boundary tiles
                 
                 ELSE
                  BMG_MSG_iGRID(p_NLz+k-1) = NLz
                 ENDIF
               ENDDO
               ENDIF

!               DO k=1, NProcK
!                  BMG_MSG_iGRID(p_NLz+k-1) = nnzs-1
!               ENDDO
               !
               ! Close the grid file.
               !
!               CLOSE (10)
            ELSE
               !
               ! Read fine-grid local dimensions
               !
!               READ(10,*) NLx, NLy, NLz
               !
               ! Close the grid file.
               !
!               CLOSE (10)
               !
               ! Fill arrays of fine-grid local dimensions
               !
               DO i=1, NProcI
                  BMG_MSG_iGRID(p_NLx+i-1) = NLx
               ENDDO
               !
               DO j=1, NProcJ
                  BMG_MSG_iGRID(p_NLy+j-1) = NLy
               ENDDO
               !
               DO k=1, NProcK
                  BMG_MSG_iGRID(p_NLz+k-1) = NLz
               ENDDO
               !
            ENDIF

         ENDIF
         
         !
         !  All:
         !
         CALL MPI_Bcast( BMG_MSG_iGRID(p_NLx), NProcI, MPI_INTEGER,
     &                   0, BMG_MPI_COMM, MPI_IERROR )
         
         CALL MPI_BCast( BMG_MSG_iGRID(p_NLy), NProcJ, MPI_INTEGER,
     &                   0, BMG_MPI_COMM, MPI_IERROR ) 

         CALL MPI_BCast( BMG_MSG_iGRID(p_NLz), NProcK, MPI_INTEGER,
     &                   0, BMG_MPI_COMM, MPI_IERROR ) 

         !
         ! Setup processor grid: 
         !
         !   MyProc  range is [1,NProcI*NProcJ*NProcK]
         !   MyProcI range is [1,NProcI]
         !   MyProcJ range is [1,NProcJ] 
         !   MyProcJ range is [1,NProcK] 
         !
!         MyProcK = (MyProc-1)/(NProcI*NProcJ) + 1
!         MyProcJ = (MyProc-(MyProcK-1)*NProcI*NProcJ-1)/NProcI + 1
!         MyProcI = MyProc 
!     &           - (MyProcK-1)*NProcI*NProcJ - (MyProcJ-1)*NProcI 
         
         p_PG = p_ProcGrid
         IF ( NProcK > 1 ) THEN
            write(0,*) 'STOP! Cannot handle NProcK > 1 !!'
            call commasmpi_abort()
            STOP
         ENDIF
         DO k=1, NProcK
            DO j=1, NProcJ
               DO i=1, NProcI
!                  BMG_MSG_iGRID(p_PG) = procmap(i,j) + 1
                  BMG_MSG_iGRID(p_PG) =  i + (j-1)*NProcI 
     &                                    + (k-1)*NProcJ*NProcI
!                  IF (  BMG_MSG_iGRID(p_PG) .ne. procmap(i,j)+1 ) THEN
!                    write(0,*) 'BMG_MSG_iGRID, procmap: ',BMG_MSG_iGRID(p_PG), procmap(i,j)+1
!                  ENDIF
                  
                  p_PG = p_PG + 1
               ENDDO
            ENDDO
         ENDDO

         !
         ! Store coordinates of each processor
         ! 
         p_PC = p_ProcCoord
         DO k=1, NProcK
            DO j=1, NProcJ
               DO i=1, NProcI
                  BMG_MSG_iGRID(p_PC)   = i
                  BMG_MSG_iGRID(p_PC+1) = j
                  BMG_MSG_iGRID(p_PC+2) = k
                  p_PC = p_PC + 3
               ENDDO
            ENDDO
         ENDDO

         !
         !  Global problem size
         ! 
         NGx = 0
         DO i=1, NProcI
            NGx = NGx + BMG_MSG_iGRID(p_NLx+i-1)
         ENDDO
         
         NGy = 0
         DO j=1, NProcJ
            NGy = NGy + BMG_MSG_iGRID(p_NLy+j-1)
         ENDDO
         
         NGz = 0
         DO k=1, NProcK
            NGz = NGz + BMG_MSG_iGRID(p_NLz+k-1)
         ENDDO

         !
         ! Global index (i,j) of the lower left corner
         !
         iGs = 1
         DO i=1, NprocI
            BMG_MSG_iGRID(p_iGs+i-1) = iGs
            iGs = iGs + BMG_MSG_iGRID(p_NLx+i-1)
         ENDDO

         jGs = 1
         DO j=1, NProcJ
            BMG_MSG_iGRID(p_jGs+j-1) = jGs
            jGs = jGs + BMG_MSG_iGRID(p_NLy+j-1)
         ENDDO

         kGs = 1
         DO k=1, NProcK
            BMG_MSG_iGRID(p_kGs+k-1) = kGs
            kGs = kGs + BMG_MSG_iGRID(p_NLz+k-1)
         ENDDO
         
         !
         !  Sanity check
         !
         IF ( NGx.NE.iGs-1 .OR. NGy.NE.jGs-1 .OR. NGz.NE.kGs-1 ) THEN

            IF ( MyProc.EQ.1 ) THEN
               WRITE(*,*) ' THIS IS STRANGE: NGx = ', NGx
               WRITE(*,*) '                  iGs = ', iGs
               WRITE(*,*) '                  NGy = ', NGy
               WRITE(*,*) '                  jGs = ', jGs
               WRITE(*,*) '                  NGz = ', NGz
               WRITE(*,*) '                  kGs = ', kGs
            ENDIF
               
            CALL MPI_FINALIZE(MPI_IERROR)
            STOP

         ENDIF

      ENDIF

C --------------------------------------------
C     Packup the local grid data
C --------------------------------------------

      ! 
      !  Store: NGx, NGy, NGz, NLx, NLy, NLz, iGs, jGs, kGs
      ! 
      BMG_MSG_iGRID(id_BMG_MSG_NGx) = NGx
      BMG_MSG_iGRID(id_BMG_MSG_NGy) = NGy
      BMG_MSG_iGRID(id_BMG_MSG_NGz) = NGz
      BMG_MSG_iGRID(id_BMG_MSG_NLx) = BMG_MSG_iGRID(p_NLx+MyProcI-1)
      BMG_MSG_iGRID(id_BMG_MSG_NLy) = BMG_MSG_iGRID(p_NLy+MyProcJ-1)
      BMG_MSG_iGRID(id_BMG_MSG_NLz) = BMG_MSG_iGRID(p_NLz+MyProcK-1)
      BMG_MSG_iGRID(id_BMG_MSG_iGs) = BMG_MSG_iGRID(p_iGs+MyProcI-1)
      BMG_MSG_iGRID(id_BMG_MSG_jGs) = BMG_MSG_iGRID(p_jGs+MyProcJ-1)
      BMG_MSG_iGRID(id_BMG_MSG_kGs) = BMG_MSG_iGRID(p_kGs+MyProcK-1)
      !
      !  Store: Nproc, NProcI, NProcJ, NProcK
      !
      BMG_MSG_iGRID(id_BMG_MSG_NProc)   = NProc
      BMG_MSG_iGRID(id_BMG_MSG_NProcI)  = NProcI
      BMG_MSG_iGRID(id_BMG_MSG_NProcJ)  = NProcJ
      BMG_MSG_iGRID(id_BMG_MSG_NProcK)  = NProcK
      !
      !  Store: MyProc, MyProcI, MyProcJ, MyProcK
      !
      BMG_MSG_iGRID(id_BMG_MSG_MyProc)  = MyProc
      BMG_MSG_iGRID(id_BMG_MSG_MyProcI) = MyProcI
      BMG_MSG_iGRID(id_BMG_MSG_MyProcJ) = MyProcJ
      BMG_MSG_iGRID(id_BMG_MSG_MyProcK) = MyProcK
      !
      !  Store/set communicators:
      !
      BMG_MSG_iGRID(id_BMG_MSG_COMM)    = BMG_MPI_COMM
      BMG_MSG_iGRID(id_BMG_MSG_COMM_xy) = 0
      BMG_MSG_iGRID(id_BMG_MSG_COMM_yz) = 0
      BMG_MSG_iGRID(id_BMG_MSG_COMM_xz) = 0

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SETUP_ProcGrid_file.f',
     &        //,5X,A,/)
 505  FORMAT (5X,A,/)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END

C ==========================================================================
C ==========================================================================
C ==========================================================================
C ==========================================================================

      SUBROUTINE EX_SETUP_BMG_parms( 
     &              BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &              MPI_MyProc, NProc, EX_MPI_COMM, EX_MPI_IERROR
     &              )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     EX_SETUP_BMG_parms on the master ( MPI_MyProc=0 ) first sets the
C     default cycle parameters, then reads customization parameters
C     from a data file, and then broadcasts the arrays.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      USE ELEC_MODULE
      
      IMPLICIT NONE

C -----------------------------
C     Includes
C -----------------------------

      INCLUDE 'mpif.h'

#include      "BMG_constants_f77.h"
      INCLUDE 'BMG_parameters_f77.h'

C ----------------------------
C     Argument Declarations
C ----------------------------

      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER   EX_MPI_COMM, EX_MPI_IERROR, MPI_MyProc, NProc

C ----------------------------
C     Local Declarations
C ----------------------------

      INTEGER   CE, i, IO_in, MyProc, READ_STATUS
      LOGICAL   exist

C ==========================================================================

C ------------------------------------------- 
C     MPI indexing:
C ------------------------------------------- 

      MyProc = MPI_MyProc


C -------------------------------------------
C     Set default parameter values
C -------------------------------------------
      
      IF ( MyProc.EQ.0 ) THEN
      
         CALL BMG3_SymStd_SETUP_parms( 
     &                    BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &                    ) 
         
      ENDIF

C ------------------------------------------- 
C     Read parameter data
C ------------------------------------------- 

      IF ( MyProc.EQ.0 ) THEN 



            !
            !  Open data file
            !

            !
            !  Read the parameters
            !
             BMG_iPARMS(id_BMG3_SETUP) = 0
             BMG_iPARMS(id_BMG3_STENCIL) = 1
             BMG_iPARMS(id_BMG3_SEPARABLE) = bmg_sep
!            write(0,*) 'BMG_iPARMS(id_BMG3_SEPARABLE) = ',
!     :        BMG_iPARMS(id_BMG3_SEPARABLE),id_BMG3_SEPARABLE
             BMG_iPARMS(id_BMG3_RELAX)  = 1
             BMG_iPARMS(id_BMG3_NRELAX_DOWN) = 1
             BMG_iPARMS(id_BMG3_NRELAX_UP)  = 1
             BMG_iPARMS(id_BMG3_RELAX_SYM )  = 0

             BMG_iPARMS(id_BMG3_CYCLE_CLASS) = bmg_cycletype
             BMG_iPARMS(id_BMG3_NCYCLE_TYPE) = 1

             BMG_iPARMS(id_BMG3_CG_SOLVER) = bmg_cgsolve
             BMG_iPARMS(id_BMG3_CG_MIN_DIM) = 3
         
             BMG_iPARMS(id_BMG2_CYCLE_CLASS) = 1
             BMG_iPARMS(id_BMG2_NCYCLE_TYPE) = 1

             BMG_iPARMS(id_BMG2_CG_SOLVER) = 0
             BMG_iPARMS(id_BMG2_CG_MIN_DIM) = 3

             BMG_iPARMS(id_BMG3_MAX_ITERS) = bmg_maxiter
             BMG_rPARMS(id_BMG3_STOP_TOL) = bmg_tol


            !
            !  Override default stopping criteria
            !
!            BMG_iPARMS(id_BMG3_STOP_TEST) = BMG_STOP_REL_RES_L2

            !
            !  Override debbugging flags
            !
!            BMG_IOFLAG(iBMG3_BUG_RES_RELAX)    = .TRUE.
!            BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) = .TRUE.

            !
            !  Set the Coarse Grid Comm Pattern
            !
            BMG_iPARMS(id_BMG3_CG_COMM) = BMG_CG_GATHER_SCATTER

            BMG_iPARMS(id_BMG3_SYNC_INITIAL_GUESS) 
     &           = BMG_SYNC_INITIAL_GUESS

      ENDIF

C ------------------------------------------- 
C     Broadcast read status
C ------------------------------------------- 

c      CALL MPI_Bcast( READ_STATUS, 1, MPI_INTEGER,
c     &                0, EX_MPI_COMM, EX_MPI_IERROR )


C ------------------------------------------- 
C     Broadcast parameters arrays
C ------------------------------------------- 

      CALL MPI_Bcast( BMG_iPARMS, NBMG_iPARMS, MPI_INTEGER,
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

      CALL MPI_Bcast( BMG_rPARMS, NBMG_rPARMS, FLOAT_PRECISION,
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

      CALL MPI_Bcast( BMG_IOFLAG, NBMG_IOFLAG, MPI_LOGICAL,
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

C ------------------------------------------- 
C     Override some I/O on the master
C ------------------------------------------- 

      IF ( MyProc.EQ.0 ) THEN

         !
         !  Override I/O Parameters
         !
         BMG_IOFLAG(iBMG3_OUT_WSPACE_POINT) = .FALSE.
         BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE)  = .FALSE.

         BMG_IOFLAG(iBMG3_OUT_ITERATIONS)   = .TRUE.
         BMG_IOFLAG(iBMG3_OUT_TIME_SETUP)   = .TRUE.
         BMG_IOFLAG(iBMG3_OUT_TIME_CYCLING) = .TRUE.
         BMG_IOFLAG(iBMG3_OUT_TIME_TOTAL)   = .TRUE.

         BMG_IOFLAG(iBMG3_BUG_PARAMETERS)   = .FALSE.
         !
         !  Debugging parameters?
         !
         IF ( BMG_IOFLAG(iBMG3_BUG_PARAMETERS) ) THEN
            CALL BMG3_SymStd_DUMP_parms( 
     &                       BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &                       ) 
         ENDIF

      ENDIF

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: EX_SETUP_BMG_parms.f',
     &        //,5X,A,/)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================


      RETURN
      END

C ==========================================================================
C ==========================================================================
C ==========================================================================

      SUBROUTINE PUTF( SO, QF, QFW, Q, QW,  t0, pot, nx,ny,nz,nor,
     :               nbw,nbe,nbs,nbn,
     &                 NLx, NLy, NLz, NGx, NGy, NGz, 
     &                 iGs, jGs, kGs, hx, hy, hz, isetup
     &               )
C ==========================================================================

      USE ELEC_MODULE
      USE COMMASMPI_MODULE
      
      IMPLICIT NONE

C ----------------
C     Includes
C
#include       "BMG_constants_f77.h"
      INCLUDE  'BMG_stencils_f77.h'

C ----------------
C     Parmeters
C
      REAL*RKIND   rHALF
      PARAMETER ( rHALF = 0.5D0 )

C ----------------
C     Arguments
C
      INTEGER  iGs, jGs, kGs, NLx, NLy, NLz, NGx, NGy, NGz, isetup
!      REAL*RKIND   hx, hy, hz, Q(NLx+2,NLy+2,NLz+2),
!     &         QF(NLx+2,NLy+2,NLz+2), SO(NLx+3,NLy+3,NLz+3,4)
      REAL*RKIND   hx, hy, hz, Q(0:NLx+1,0:NLy+1,0:NLz+1),
     &         QF(0:NLx+1,0:NLy+1,0:NLz+1), 
     &         QFW(-nbw:nx+1+nbe,-nbs:ny+1+nbn,0:NLz+1), 
     &         QW(-nbw:nx+1+nbe,-nbs:ny+1+nbn,0:NLz+1), 
     &         SO(0:NLx+2,0:NLy+2,0:NLz+2,4)


      integer nx,ny,nz,nor
      integer nbw,nbe,nbs,nbn
      integer, parameter    :: ng1 = 1
      
      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real pot(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

C ----------------
C     Local
C
      INTEGER  i, i1, i2, iGf, j, j1, j2, jGf, k, k1, k2, kGf, 
     &         iBEG, iEND, jBEG, jEND, kBEG, kEND, 
     &         NLx_g, NLy_g, NLz_g, NLSOx_g, NLSOy_g, NLSOz_g,
     &	       is, js, kss
     
      integer ibegw,iendw,jbegw,jendw
      REAL*RKIND   cm, fs, h2, s, xh, yh, zh, alpha
      
      integer ixe, jye


      real eperao
      parameter (eperao  = 8.8592e-12 )

C ==========================================================================

C ----------------------------------------
C     Initalize constants
C ----------------------------------------

      !
      !  Standard halo with depth = [1,1]
      !
      NLx_g = NLx + 1
      NLy_g = NLy + 1
      NLz_g = NLz + 1

      i1 = NLx + 0
      j1 = NLy + 0
      k1 = NLz + 0

      i2 = NLx - 1
      j2 = NLy - 1
      k2 = NLz - 1

      !
      !  Global indexing
      !
      iGf=iGs+i2-1+1
      jGf=jGs+j2-1+1
      kGf=kGs+k2-1+1
      
!      write(*,*) 'igf,jgf,kgf = ',igf,jgf,kgf
!      write(*,*) 'NGx,NGy,NGz = ',NGx,NGy,NGz

      !
      !  Stencil halo with depth = [1,2]
      !
      NLSOx_g = NLx + 2
      NLSOy_g = NLy + 2
      NLSOz_g = NLz + 2

      !
      !
      !  Grid Spacing
      h2=hx*hy*hz

      xh=hy*hz/hx
      yh=hx*hz/hy
      zh=hx*hy/hz

      !
      !  Loop Boundaries: iGs, jGs, kGs
      !
      IF ( iGs.EQ.1 ) THEN
         iBEG = 2
         ibegw = 1
      ELSE
         iBEG = 1
         ibegw = 1
      ENDIF

      IF ( jGs.EQ.1 ) THEN
         jBEG = 2
         jbegw = 1
      ELSE
         jBEG = 1
         jbegw = 1
      ENDIF

      IF ( kGs.EQ.1 ) THEN
         kBEG = 2
      ELSE
         kBEG = 1
      ENDIF

      !
      !  Loop Boundaries: iGf, jGf, kGf
      !
      IF ( iGf.EQ.NGx ) THEN
         iEND = i1
         iendw = nx-1
      ELSE
         iEND = NLx_g
         iendw = nx
      ENDIF

      IF ( jGf.EQ.NGy ) THEN
         jEND = j1
         jendw = ny-1
      ELSE
         jEND = NLy_g
         jendw = ny
      ENDIF

      IF ( kGf.EQ.NGz ) THEN
         kEND = k1
      ELSE
         kEND = NLz_g
      ENDIF
         
      
      IF ( isetup == 1 ) THEN
      
C ----------------------------------------
C     Zero solution
C ----------------------------------------
      
      DO k=0, NLz_g
         DO j=0, NLy_g
            DO i=0, NLx_g

               Q(i,j,k)   = rZERO
               QF(i,j,k)  = rZERO

            END DO
         END DO
      END DO

C ----------------------------------------
C     Zero Stencil
C ----------------------------------------

      DO k=0, NLSOz_g
         DO j=0, NLSOy_g
            DO i=0, NLSOx_g

               SO(i,j,k,kp)  = rZERO
               SO(i,j,k,kpw) = rZERO
               SO(i,j,k,kps) = rZERO
               SO(i,j,k,kb)  = rZERO

            END DO
         END DO
      END DO
      
      ENDIF


C ----------------------------------------
C     Compute the source (right-hand-side) multiplied by the grid point volume
C    Note that this is actually -rhs, so multiply by -1
C    For electric poisson eq., then qf = (rho/epsilon)*(grid volume)
C ----------------------------------------
      
      ixe = nx - 1
      IF ( myproci < nproci ) ixe = nx
      
      jye = ny - 1
      IF ( myprocj < nprocj ) jye = ny
      
!      write(*,*) 'putf: rank, i1 = ',my_rank,nbw,nbe,i1-1,nx+nbe-(-nbw+1)
      
      IF ( .true. ) THEN
      
      DO k=1,k1 ! k1/2-3, k1/2+3
         DO j=-nbs+1, ny+nbn
            DO i=-nbw+1, nx+nbe
               !
            IF ( k .ge. 1 .and. k .le. nz-1 .and. 
     :           j .ge. 1 .and. j .le. jye  .and.
     :           i .ge. 1 .and. i .le. ixe ) THEN
!               qfw(i,j,k)    = (t0(i,j,k)/eperao)*hx*hy*(dlz/z1d2(k,3))
               qfw(i,j,k)    = (t0(i,j,k)/eperao)*x1d2(i+nbw)*y1d2(j+nbs)*(dlz/z1d2(k,3))
               qw(i,j,k)     =  pot(i,j,k)
!               qfw(i,j,k)    = (rhofair(k)/eperao)*hx*hy*(dlz/z1d2(k,3))
            ELSE
!               qfw(i,j,k)    = (rhofair(k)/eperao)*hx*hy*(dlz/z1d2(k,3))
               qfw(i,j,k)    = (rhofair(k)/eperao)*x1d2(i+nbw)*y1d2(j+nbs)*(dlz/z1d2(k,3))

!               qw(i,j,k) = potfair(k)
            ENDIF
               !
            ENDDO
         ENDDO
      ENDDO
      
      ELSE


      DO k=1,k1 ! k1/2-3, k1/2+3
         DO j=1, j1
            DO i=1, i1
               !
               qf(i,j,k)    = 0.0 ! (rhofair(k)/eperao)*hx*hy*(dlz/z1d2(k,3))
               !
            ENDDO
         ENDDO
      ENDDO

      ENDIF
      
      
      IF ( isetup == 1 ) THEN

C ----------------------------------------
C     Stencil coefficients
C ----------------------------------------
      !
      !  South
      !
!      write(0,*) 'putf setup: rank = ', my_rank
!      write(0,*) 'ibeg,iend = ',ibeg,iend
!      write(0,*) 'jbeg,jend = ',jbeg,jend
      DO k=1, kEND
         DO j=jBEG,jEND
            DO i=1, iEND
               !
               !
               SO(i,j,k,kps) = x1d2(i)*(dlz/z1d2(k,3))/y1d2(j)
            ENDDO
         ENDDO
      ENDDO

      !
      !  West
      !
      DO k=1, kEND
         DO j=1, jEND
            DO i=iBEG, iEND
               !
               !
               SO(i,j,k,kpw) =  y1d2(j)*(dlz/z1d2(k,3))/x1d2(i)
            ENDDO
         ENDDO
      ENDDO

      !
      !  Bottom
      !
      DO k=kBEG, kEND
         DO j=1, jEND
            DO i=1, iEND
               !
               SO(i,j,k,kb) = x1d2(i)*y1d2(j)/(dlz/z1d2(k,3))
            ENDDO
         ENDDO
      ENDDO

C ----------------------------------------
C     Initial Guess
C ----------------------------------------

      DO k=1, k1
         DO j=1, j1
            DO i=1, i1
!               q(i,j,k) = 0 ! BMG_rand(r)
            ENDDO
         ENDDO
      ENDDO

C ---------------------------------------
C     Boundary Conditions
C ---------------------------------------

      !
      !  West boundary
      !
C      IF ( iGs.EQ.1 ) THEN
C
C         DO k=2, kEND
C            DO j=2, jEND
Cc               SO(2,j,k,kps)  = rHALF*SO(2,j,k,kps)
Cc               SO(2,j,k,kb)   = rHALF*SO(2,j,k,kb)
C            ENDDO
C         ENDDO
C
C
C         IF ( ibl.EQ.1) THEN
C
C            DO j=2, jEND
C               DO k=2, kEND
C!                  SO(2,j,k,kp) = SO(2,j,k,kp) + hy*hz
C                  SO(1,j,k,kpw) = xh
C               ENDDO
C            ENDDO
C
C         ENDIF
C
C      ENDIF

      !
      !  Bottom Boundary
      !
      IF ( kGs.EQ.1 ) THEN
         
         DO j=2, jEND
            DO i=2, iEND
c               SO(i,j,2,kpw)  = rHALF*SO(i,j,2,kpw)
c               SO(i,j,2,kps)  = rHALF*SO(i,j,2,kps)
            ENDDO
         ENDDO

         IF ( .true. ) THEN  ! dirichlet

           DO k=1,1
            DO j=0, nly+1
               DO i=0, nlx+1
!           DO j=0, jEND
!              DO i=0, iEND
                  Q(i,j,k-1) = 0  ! set potential at one point LOWER
!                  QF(i,j,k) = 1
               SO(i,j,k,kb) =  x1d2(i)*y1d2(j)/(dlz/z1d2(k,3))
               ENDDO
            ENDDO
           ENDDO

            DO j=1, jEND
               DO i=1, iEND
c                  SO(i,j,2,kp) = SO(i,j,2,kp) + hx*hy
               ENDDO
            ENDDO
            
         ELSE ! IF ( ibzb.EQ.0) THEN

           alpha = -1 ! gradient at edge
           DO k=1,1
!            DO j=1, nly+2
!               DO i=1, nlx+2
            DO j=1, jEND
               DO i=1, iEND
!                  Q(i,j,k) = 0
                  QF(i,j,1) = QF(i,j,1) - alpha*x1d2(i)*y1d2(j)
                 SO(i,j,k,kp) = SO(i,j,k,kp) -  x1d2(i)*y1d2(j)/(dlz/z1d2(k,3)) ! zh  ! 0.1*zh
               ENDDO
            ENDDO
           ENDDO


         ENDIF
         
      ENDIF


      !
      !  Top boundary
      !
      IF ( kGf.EQ.NGz ) THEN 
         
         DO j=1, jEND
            DO i=1, iEND
c               SO(i,j,k1,kpw) = rHALF*SO(i,j,k1,kpw)
c               SO(i,j,k1,kps) = rHALF*SO(i,j,k1,kps)
            ENDDO
         ENDDO

         IF( .true. ) THEN ! Dirichlet

!           DO k= k1-1,k1+1 ! nlz+1,nlz+2
!             write(*,*) 'k,potfair = ',k,potfair(k)
!           ENDDO

           DO k= k1+1,k1+1 ! nlz+1,nlz+2
            DO j=0, nly+1
               DO i=0, nlx+1
!           DO j=jBEG-1, jEND
!              DO i=iBEG-1, iEND
                  Q(i,j,k) = potfair(k)
!                  Q(i,j,k) = 10
!                  QF(i,j,k) = 1
                 SO(i,j,k,kb) = SO(i,j,k,kb) + x1d2(i)*y1d2(j)/(dlz/z1d2(k,3))
               ENDDO
            ENDDO
           ENDDO

            DO j=1, jEND
               DO i=1, iEND
c                  SO(i,j,k1,kp) = SO(i,j,k1,kp) + hx*hy
               ENDDO
            ENDDO

C           DO k=nlz-1,nlz+2
C            DO j=2, jEND
C               DO i=2, iEND
C                  Q(i,j,k) = 1
C               ENDDO
C            ENDDO
C           ENDDO



         ENDIF

      ENDIF




      DO k=1, k1
         DO j=1, j1
            DO i=1, i1
               cm=rONE
               IF ( ( i.EQ.2 .AND. iGs.EQ.1 ) 
     &              .OR. ( i.EQ.i1 .AND. iGf.EQ.NGx ) ) THEN
                  cm = rHALF*cm
               ENDIF
               IF ( ( j.EQ.2 .AND. jGs.EQ.1 ) 
     &              .OR. ( j.EQ.j1 .AND. jGf.EQ.NGy ) ) THEN
                  cm = rHALF*cm
               ENDIF
               IF ( ( k.EQ.2 .AND. kGs.EQ.1 ) 
     &              .OR. ( k.EQ.k1 .AND. kGf.EQ.NGz ) ) THEN
                  cm = rHALF*cm
               ENDIF
               !
               SO(i,j,k,kp) = cm*SO(i,j,k,kp)
     &                      + SO(i,j,k,kpw)+SO(i,j+1,k,kps)
     &                      + SO(i+1,j,k,kpw)+SO(i,j,k,kps)
     &                      + SO(i,j,k,kb)+SO(i,j,k+1,kb)
               !
!               qf(i,j,k) = cm*qf(i,j,k)
               !
            ENDDO
         ENDDO
      ENDDO
      
      ENDIF

C ==========================================================================

      RETURN
      END
C ==========================================================================
C ==========================================================================
C ==========================================================================

      SUBROUTINE PUTFLGT( SO, QF, QFW, Q, QW,  t0, pot, nx,ny,nz,nor,
     :               nbw,nbe,nbs,nbn,
     &                 NLx, NLy, NLz, NGx, NGy, NGz, 
     &                 iGs, jGs, kGs, hx, hy, hz, isetup
     &               )
C ==========================================================================

      USE ELEC_MODULE
      USE COMMASMPI_MODULE
      
      IMPLICIT NONE

C ----------------
C     Includes
C
#include       "BMG_constants_f77.h"
      INCLUDE  'BMG_stencils_f77.h'

C ----------------
C     Parmeters
C
      REAL*RKIND   rHALF
      PARAMETER ( rHALF = 0.5D0 )

C ----------------
C     Arguments
C
      INTEGER  iGs, jGs, kGs, NLx, NLy, NLz, NGx, NGy, NGz, isetup
!      REAL*RKIND   hx, hy, hz, Q(NLx+2,NLy+2,NLz+2),
!     &         QF(NLx+2,NLy+2,NLz+2), SO(NLx+3,NLy+3,NLz+3,4)
      REAL*RKIND   hx, hy, hz, Q(0:NLx+1,0:NLy+1,0:NLz+1),
     &         QF(0:NLx+1,0:NLy+1,0:NLz+1), 
     &         QFW(-nbw:nx+1+nbe,-nbs:ny+1+nbn,0:NLz+1), 
     &         QW(-nbw:nx+1+nbe,-nbs:ny+1+nbn,0:NLz+1), 
     &         SO(0:NLx+2,0:NLy+2,0:NLz+2,4)


      integer nx,ny,nz,nor
      integer nbw,nbe,nbs,nbn
      integer, parameter    :: ng1 = 1
      
      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)
      real pot(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

C ----------------
C     Local
C
      INTEGER  i, i1, i2, iGf, j, j1, j2, jGf, k, k1, k2, kGf, 
     &         iBEG, iEND, jBEG, jEND, kBEG, kEND, 
     &         NLx_g, NLy_g, NLz_g, NLSOx_g, NLSOy_g, NLSOz_g,
     &	       is, js, kss
     
      integer ibegw,iendw,jbegw,jendw
      REAL*RKIND   cm, fs, h2, s, xh, yh, zh, alpha
      
      integer ixe, jye


      real eperao
      parameter (eperao  = 8.8592e-12 )

C ==========================================================================

C ----------------------------------------
C     Initalize constants
C ----------------------------------------


c         write(0,*) 'sizes of qf: ',size(qf,1),size(qf,2),size(qf,3)
c         write(0,*) 'sizes of qfw: ',size(qfw,1),size(qfw,2),size(qfw,3)

      !
      !  Standard halo with depth = [1,1]
      !
      NLx_g = NLx + 1
      NLy_g = NLy + 1
      NLz_g = NLz + 1

      i1 = NLx + 0
      j1 = NLy + 0
      k1 = NLz + 0

      i2 = NLx - 1
      j2 = NLy - 1
      k2 = NLz - 1

      !
      !  Global indexing
      !
      iGf=iGs+i2-1+1
      jGf=jGs+j2-1+1
      kGf=kGs+k2-1+1
      
!      write(*,*) 'igf,jgf,kgf = ',igf,jgf,kgf
!      write(*,*) 'NGx,NGy,NGz = ',NGx,NGy,NGz

      !
      !  Stencil halo with depth = [1,2]
      !
      NLSOx_g = NLx + 2
      NLSOy_g = NLy + 2
      NLSOz_g = NLz + 2

      !
      !
      !  Grid Spacing
      h2=hx*hy*hz

      xh=hy*hz/hx
      yh=hx*hz/hy
      zh=hx*hy/hz

      !
      !  Loop Boundaries: iGs, jGs, kGs
      !
      IF ( iGs.EQ.1 ) THEN
         iBEG = 2
         ibegw = 1
      ELSE
         iBEG = 1
         ibegw = 1
      ENDIF

      IF ( jGs.EQ.1 ) THEN
         jBEG = 2
         jbegw = 1
      ELSE
         jBEG = 1
         jbegw = 1
      ENDIF

      IF ( kGs.EQ.1 ) THEN
         kBEG = 2
      ELSE
         kBEG = 1
      ENDIF

      !
      !  Loop Boundaries: iGf, jGf, kGf
      !
      IF ( iGf.EQ.NGx ) THEN
         iEND = i1
         iendw = nx-1
      ELSE
         iEND = NLx_g
         iendw = nx
      ENDIF

      IF ( jGf.EQ.NGy ) THEN
         jEND = j1
         jendw = ny-1
      ELSE
         jEND = NLy_g
         jendw = ny
      ENDIF

      IF ( kGf.EQ.NGz ) THEN
         kEND = k1
      ELSE
         kEND = NLz_g
      ENDIF
         
      
      IF ( isetup == 1 ) THEN
      
C ----------------------------------------
C     Zero solution
C ----------------------------------------
      
      DO k=0, NLz_g
         DO j=0, NLy_g
            DO i=0, NLx_g

               Q(i,j,k)   = rZERO
               QF(i,j,k)  = rZERO

            END DO
         END DO
      END DO

C ----------------------------------------
C     Zero Stencil
C ----------------------------------------

      DO k=0, NLSOz_g
         DO j=0, NLSOy_g
            DO i=0, NLSOx_g

               SO(i,j,k,kp)  = rZERO
               SO(i,j,k,kpw) = rZERO
               SO(i,j,k,kps) = rZERO
               SO(i,j,k,kb)  = rZERO

            END DO
         END DO
      END DO
      
      ENDIF


C ----------------------------------------
C     Compute the source (right-hand-side) multiplied by the grid point volume
C    Note that this is actually -rhs, so multiply by -1
C    For electric poisson eq., then qf = (rho/epsilon)*(grid volume)
C ----------------------------------------

      IF ( isetup == 0 ) THEN ! only run this if lightning grid is same as dynamics (no refinement)
      
      ixe = nx - 1
      IF ( myproci < nproci ) ixe = nx
      
      jye = ny - 1
      IF ( myprocj < nprocj ) jye = ny
      
!      write(*,*) 'putf: rank, i1 = ',my_rank,nbw,nbe,i1-1,nx+nbe-(-nbw+1)
      
      IF ( .true. ) THEN
      
      DO k=1,k1 ! k1/2-3, k1/2+3
         DO j=-nbs+1, ny+nbn
            DO i=-nbw+1, nx+nbe
               !
            IF ( k .ge. 1 .and. k .le. nz-1 .and. 
     :           j .ge. 1 .and. j .le. jye  .and.
     :           i .ge. 1 .and. i .le. ixe ) THEN
               qfw(i,j,k)    = (t0(i,j,k)/eperao)*x1d2lgt(i+nbw)*y1d2lgt(j+nbs)*(dlzlgt/z1d2lgt(k,3))
               qw(i,j,k)     =  pot(i,j,k)
!               qfw(i,j,k)    = (rhofair(k)/eperao)*hx*hy*(dlz/z1d2lgt(k,3))
            ELSE
               qfw(i,j,k)    = (rhofair(k)/eperao)*x1d2lgt(i+nbw)*y1d2lgt(j+nbs)*(dlzlgt/z1d2lgt(k,3))
!               qw(i,j,k) = potfair(k)
            ENDIF
               !
            ENDDO
         ENDDO
      ENDDO
      
      ELSE


      DO k=1,k1 ! k1/2-3, k1/2+3
         DO j=1, j1
            DO i=1, i1
               !
               qf(i,j,k)    = 0.0 ! (rhofair(k)/eperao)*hx*hy*(dlzlgt/z1d2lgt(k,3))
               !
            ENDDO
         ENDDO
      ENDDO

      ENDIF
      
      ENDIF
      
      
      IF ( isetup == 1 ) THEN

C ----------------------------------------
C     Stencil coefficients
C ----------------------------------------
      !
      !  South
      !
      DO k=1, kEND
         DO j=jBEG,jEND
            DO i=1, iEND
               !
               !
               SO(i,j,k,kps) = x1d2lgt(i)*(dlzlgt/z1d2lgt(k,3))/y1d2lgt(j)
            ENDDO
         ENDDO
      ENDDO

      !
      !  West
      !
      DO k=1, kEND
         DO j=1, jEND
            DO i=iBEG, iEND
               !
               !
               SO(i,j,k,kpw) =  y1d2lgt(j)*(dlzlgt/z1d2lgt(k,3))/x1d2lgt(i)
            ENDDO
         ENDDO
      ENDDO

      !
      !  Bottom
      !
      DO k=kBEG, kEND
         DO j=1, jEND
            DO i=1, iEND
               !
               SO(i,j,k,kb) = x1d2lgt(i)*y1d2lgt(j)/(dlzlgt/z1d2lgt(k,3))
            ENDDO
         ENDDO
      ENDDO

C ----------------------------------------
C     Initial Guess
C ----------------------------------------

      DO k=1, k1
         DO j=1, j1
            DO i=1, i1
!               q(i,j,k) = 0 ! BMG_rand(r)
            ENDDO
         ENDDO
      ENDDO

C ---------------------------------------
C     Boundary Conditions
C ---------------------------------------

      !
      !  West boundary
      !
C      IF ( iGs.EQ.1 ) THEN
C
C         DO k=2, kEND
C            DO j=2, jEND
Cc               SO(2,j,k,kps)  = rHALF*SO(2,j,k,kps)
Cc               SO(2,j,k,kb)   = rHALF*SO(2,j,k,kb)
C            ENDDO
C         ENDDO
C
C
C         IF ( ibl.EQ.1) THEN
C
C            DO j=2, jEND
C               DO k=2, kEND
C!                  SO(2,j,k,kp) = SO(2,j,k,kp) + hy*hz
C                  SO(1,j,k,kpw) = xh
C               ENDDO
C            ENDDO
C
C         ENDIF
C
C      ENDIF

      !
      !  Bottom Boundary
      !
      IF ( kGs.EQ.1 ) THEN
         
         DO j=2, jEND
            DO i=2, iEND
c               SO(i,j,2,kpw)  = rHALF*SO(i,j,2,kpw)
c               SO(i,j,2,kps)  = rHALF*SO(i,j,2,kps)
            ENDDO
         ENDDO

         IF ( .true. ) THEN  ! dirichlet

           DO k=1,1
            DO j=0, nly+1
               DO i=0, nlx+1
!           DO j=0, jEND
!              DO i=0, iEND
                  Q(i,j,k-1) = 0  ! set potential at one point LOWER
!                  QF(i,j,k) = 1
               SO(i,j,k,kb) =  x1d2lgt(i)*y1d2lgt(j)/(dlzlgt/z1d2lgt(k,3))
               ENDDO
            ENDDO
           ENDDO

            DO j=1, jEND
               DO i=1, iEND
c                  SO(i,j,2,kp) = SO(i,j,2,kp) + hx*hy
               ENDDO
            ENDDO
            
         ELSE ! IF ( ibzb.EQ.0) THEN

           alpha = -1 ! gradient at edge
           DO k=1,1
!            DO j=1, nly+2
!               DO i=1, nlx+2
            DO j=1, jEND
               DO i=1, iEND
!                  Q(i,j,k) = 0
                  QF(i,j,1) = QF(i,j,1) - alpha*x1d2lgt(i)*y1d2lgt(j)
                 SO(i,j,k,kp) = SO(i,j,k,kp) -  x1d2lgt(i)*y1d2lgt(j)/(dlzlgt/z1d2lgt(k,3)) ! zh  ! 0.1*zh
               ENDDO
            ENDDO
           ENDDO


         ENDIF
         
      ENDIF


      !
      !  Top boundary
      !
      IF ( kGf.EQ.NGz ) THEN 
         
         DO j=1, jEND
            DO i=1, iEND
c               SO(i,j,k1,kpw) = rHALF*SO(i,j,k1,kpw)
c               SO(i,j,k1,kps) = rHALF*SO(i,j,k1,kps)
            ENDDO
         ENDDO

         IF( .true. ) THEN ! Dirichlet

!           DO k= k1-1,k1+1 ! nlz+1,nlz+2
!             write(*,*) 'k,potfair = ',k,potfair(k)
!           ENDDO

           DO k= k1+1,k1+1 ! nlz+1,nlz+2
            DO j=0, nly+1
               DO i=0, nlx+1
!           DO j=jBEG-1, jEND
!              DO i=iBEG-1, iEND
                  Q(i,j,k) = potfair(k)
!                  Q(i,j,k) = 10
!                  QF(i,j,k) = 1
                 SO(i,j,k,kb) = SO(i,j,k,kb) + x1d2lgt(i)*y1d2lgt(j)/(dlzlgt/z1d2lgt(k,3))
               ENDDO
            ENDDO
           ENDDO

            DO j=1, jEND
               DO i=1, iEND
c                  SO(i,j,k1,kp) = SO(i,j,k1,kp) + hx*hy
               ENDDO
            ENDDO

C           DO k=nlz-1,nlz+2
C            DO j=2, jEND
C               DO i=2, iEND
C                  Q(i,j,k) = 1
C               ENDDO
C            ENDDO
C           ENDDO



         ENDIF

      ENDIF




      DO k=1, k1
         DO j=1, j1
            DO i=1, i1
               cm=rONE
               IF ( ( i.EQ.2 .AND. iGs.EQ.1 ) 
     &              .OR. ( i.EQ.i1 .AND. iGf.EQ.NGx ) ) THEN
                  cm = rHALF*cm
               ENDIF
               IF ( ( j.EQ.2 .AND. jGs.EQ.1 ) 
     &              .OR. ( j.EQ.j1 .AND. jGf.EQ.NGy ) ) THEN
                  cm = rHALF*cm
               ENDIF
               IF ( ( k.EQ.2 .AND. kGs.EQ.1 ) 
     &              .OR. ( k.EQ.k1 .AND. kGf.EQ.NGz ) ) THEN
                  cm = rHALF*cm
               ENDIF
               !
               SO(i,j,k,kp) = cm*SO(i,j,k,kp)
     &                      + SO(i,j,k,kpw)+SO(i,j+1,k,kps)
     &                      + SO(i+1,j,k,kpw)+SO(i,j,k,kps)
     &                      + SO(i,j,k,kb)+SO(i,j,k+1,kb)
               !
!               qf(i,j,k) = cm*qf(i,j,k)
               !
            ENDDO
         ENDDO
      ENDDO
      
      ENDIF

C ==========================================================================

      RETURN
      END

C ==========================================================================
C ==========================================================================
C ==========================================================================

      SUBROUTINE PRINTQ( SO, QF, Q,
     &                 NLx, NLy, NLz, NGx, NGy, NGz, 
     &                 iGs, jGs, kGs, hx, hy, hz, kg 
     &               )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Creates the stencil representation of the discretized PDE.
C     
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      USE ELEC_MODULE
      
      IMPLICIT NONE

C ----------------
C     Includes
C
#include       "BMG_constants_f77.h"
      INCLUDE  'BMG_stencils_f77.h'
      
C ----------------
C     Parmeters
C
      REAL*RKIND   rHALF
      PARAMETER ( rHALF = 0.5D0 )

      real eperao
      parameter (eperao  = 8.8592e-12 )

C ----------------
C     Arguments
C
      INTEGER  iGs, jGs, kGs, NLx, NLy, NLz, NGx, NGy, NGz, kg
c      REAL*RKIND   hx, hy, hz, Q(NLx+2,NLy+2,NLz+2),
c     &         QF(NLx+2,NLy+2,NLz+2), SO(NLx+3,NLy+3,NLz+3,4)
      REAL*RKIND   hx, hy, hz, Q(0:NLx+1,0:NLy+1,0:NLz+1),
     &         QF(0:NLx+1,0:NLy+1,0:NLz+1), 
     &         SO(0:NLx+2,0:NLy+2,0:NLz+2,4)

C ----------------
C     Local
C
      INTEGER  i, i1, i2, iGf, j, j1, j2, jGf, k, k1, k2, kGf, kz,
     &         iBEG, iEND, jBEG, jEND, kBEG, kEND, 
     &         NLx_g, NLy_g, NLz_g, NLSOx_g, NLSOy_g, NLSOz_g,
     &	       is, js, kss
      REAL*RKIND   cm, fs, h2, s, xh, yh, zh

       real :: maxq,minq, maxrho,minrho, rho
       integer :: ixmx, jymx, kzmx

C -----------------
C     Includes
C

C ==========================================================================

        DO k=0,nlz+1
         kz = Max(1, Min(nlz,k))
!          write(6,*) 'k,z,q = ',k,z1d2(kz,1),q(3,3,k),q(3,nly-1,k) !,q(nlx-1,nly-1,k)
          write(6,*) 'k,z,q = ',k,z1d2(kz,1),q(nlx/2,nly/2,k),qf(nlx/2,nly/2,k)*eperao !,q(nlx-1,nly-1,k)
     :        ,(q(nlx/2,nly/2,k)-q(nlx/2,nly/2,Min(k+1,nlz+1)))*z1d2(kz,3)/dlz  !   ,SO(3,3,k,kp)
c     :        ,(q(3,3,k)-q(3,3,Min(k+1,nlz+1)))*z1d2(kz,3)/dlz  !   ,SO(3,3,k,kp)
        ENDDO
        
        DO i=1,nlx+1
          write(6,*) 'i,q = ',i,q(i,nly-1,nlz/2), SO(i,nly-2,nlz/2,kpw),
     :  SO(i,nly-2,nlz/2,kp)
        ENDDO
        
        maxq = -1e20
        minq =  1e20
        maxrho = 0
        minrho = 0
        DO k = 1,nlz
         DO j = 1,nly
          DO i = 1,nlx
            IF ( q(i,j,k) .gt. maxq ) THEN
              maxq = q(i,j,k)
            ELSEIF ( q(i,j,k) .lt. minq ) THEN
              minq = q(i,j,k)
            ENDIF

            rho = qf(i,j,k)*eperao/(hx*hy*(dlzlgt/z1d2lgt(k,3)))
            IF ( rho .gt. maxrho ) THEN
              maxrho = rho
              ixmx = i
              jymx = j
              kzmx = k
            ELSEIF (rho .lt. minrho ) THEN
              minrho = rho
            ENDIF
          ENDDO
         ENDDO
        ENDDO
        
        write(6,*) 'minq,maxq,minrho,maxrho = ',minq,maxq,minrho,maxrho
        write(6,*) 'ixmx,jymx,kzmx = ',ixmx,jymx,kzmx
        DO k=1,nlz
         kz = Max(1, Min(nlz,k))
          rho = qf(ixmx,jymx,k)*eperao/(hx*hy*(dlzlgt/z1d2lgt(k,3)))
!          write(6,*) 'k,z,q = ',k,z1d2(kz,1),q(3,3,k),q(3,nly-1,k) !,q(nlx-1,nly-1,k)
          write(6,*) 'k,z,q = ',k,z1d2(kz,1),q(ixmx,jymx,k),rho !,q(nlx-1,nly-1,k)
     :        ,(q(ixmx,jymx,k)-q(ixmx,jymx,Min(k+1,nlz+1)))*z1d2(kz,3)/dlz  !   ,SO(3,3,k,kp)
c     :        ,(q(3,3,k)-q(3,3,Min(k+1,nlz+1)))*z1d2(kz,3)/dlz  !   ,SO(3,3,k,kp)
        ENDDO
        
      RETURN
      END

C ==========================================================================
C ==========================================================================
C ==========================================================================

      SUBROUTINE PUTPHI( pot2, elec, nx,ny,nz,nor,
     :               gxt,gyt,gzt,
     :               nbw,nbe,nbs,nbn,iunit,
     &                 NLx, NLy, NLz, NGx, NGy, NGz, 
     &                 iGs, jGs, kGs, hx, hy, hz
     &               )
C ==========================================================================

      USE GRID_MODULE
      USE ELEC_MODULE
      USE INDEX_MODULE

      IMPLICIT NONE

C ----------------
C     Includes
C
! #include       "BMG_constants_f77.h"
!      INCLUDE  'BMG_stencils_f77.h'


C ----------------
C     Arguments
C

      TYPE(VARIABLE)     :: elec(neelec)

      INTEGER  iGs, jGs, kGs, NLx, NLy, NLz, NGx, NGy, NGz, kg
!      REAL*RKIND   hx, hy, hz, Q(NLx+2,NLy+2,NLz+2),
!     &         QF(NLx+2,NLy+2,NLz+2), SO(NLx+3,NLy+3,NLz+3,4)
      REAL*RKIND   hx, hy, hz,
     &         pot2(-nbw:nx+1+nbe,-nbs:ny+1+nbn,0:NLz+1)

      real    :: gxt(-nor+1:nx+nor,4), gyt(-nor+1:ny+nor,4), gzt(-nor+1:nz+nor,4)

      integer nx,ny,nz,nor
      integer nbw,nbe,nbs,nbn,iunit
      integer, parameter    :: ng1 = 1
      
      real t0(-nor+ng1:nx+nor,-nor+ng1:ny+nor,-nor+ng1:nz+nor)

C ----------------
C     Local
C
      INTEGER  i, i1, i2, iGf, j, j1, j2, jGf, k, k1, k2, kGf, 
     &         iBEG, iEND, jBEG, jEND, kBEG, kEND, 
     &         NLx_g, NLy_g, NLz_g, NLSOx_g, NLSOy_g, NLSOz_g,
     &	       is, js, kss
     
      integer ibegw,iendw,jbegw,jendw
      REAL*RKIND   cm, fs, h2, s, xh, yh, zh, alpha


      real eperao
      parameter (eperao  = 8.8592e-12 )

C ==========================================================================


C ----------------------------------------
C     Initalize constants
C ----------------------------------------

      !
      !  Standard halo with depth = [1,1]
      !
      NLx_g = NLx + 1
      NLy_g = NLy + 1
      NLz_g = NLz + 1

      i1 = NLx + 0
      j1 = NLy + 0
      k1 = NLz + 0

      i2 = NLx - 1
      j2 = NLy - 1
      k2 = NLz - 1

      !
      !  Global indexing
      !
      iGf=iGs+i2-1+1
      jGf=jGs+j2-1+1
      kGf=kGs+k2-1+1
      
!      write(*,*) 'igf,jgf,kgf = ',igf,jgf,kgf
!      write(*,*) 'NGx,NGy,NGz = ',NGx,NGy,NGz

      !
      !  Stencil halo with depth = [1,2]
      !
      NLSOx_g = NLx + 2
      NLSOy_g = NLy + 2
      NLSOz_g = NLz + 2

      !
      !
      !  Grid Spacing
      h2=hx*hy*hz

      xh=hy*hz/hx
      yh=hx*hz/hy
      zh=hx*hy/hz

      !
      !  Loop Boundaries: iGs, jGs, kGs
      !
      IF ( iGs.EQ.1 ) THEN
         iBEG = 2
         ibegw = 1
      ELSE
         iBEG = 1
         ibegw = 1
      ENDIF

      IF ( jGs.EQ.1 ) THEN
         jBEG = 2
         jbegw = 1
      ELSE
         jBEG = 1
         jbegw = 1
      ENDIF

      IF ( kGs.EQ.1 ) THEN
         kBEG = 2
      ELSE
         kBEG = 1
      ENDIF

      !
      !  Loop Boundaries: iGf, jGf, kGf
      !
      IF ( iGf.EQ.NGx ) THEN
         iEND = i1
         iendw = nx-1
      ELSE
         iEND = NLx_g
         iendw = nx
      ENDIF

      IF ( jGf.EQ.NGy ) THEN
         jEND = j1
         jendw = ny-1
      ELSE
         jEND = NLy_g
         jendw = ny
      ENDIF

      IF ( kGf.EQ.NGz ) THEN
         kEND = k1
      ELSE
         kEND = NLz_g
      ENDIF
       

c
c  put pot into elec
c
      do k = 0,nz+1
      do j = 0,ny+1
      do i = 0,nx+1
         elec(ipot)%flt3d(i,j,k) = pot2(i,j,k)
      end do
      end do
      end do

!       IF ( iestag .eq. 1 .and. lcalce ) THEN
!       IF ( iestag .eq. 1 ) THEN

        call efieldstag
     >  (nx,ny,nz,nor
     >  ,elec
     :  ,gxt,gyt,gzt
     :  ,nbw,nbe+1,nbs,nbn+1,NLSOz_g
     :  ,pot2
     :  ,iunit)  
       
!       ENDIF


      
      RETURN
      END

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE commas_SymStd_SOLVE_boxmg(
     &                NLx, NLy, NLz, NGx, NGy, NGz, 
     &                iGs, jGs, kGs, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm, 
     &                SO, NSOm, SOR, NSORm, CI, NCIm,
     &                IGRD, NOGm, NOG,
     &                BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                BMG_rWORK_PL, NBMG_rWORK_PLm,
     &                BMG_iWORK_CS, NBMG_iWORK_CSm,
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                iWorkMSG, NMSGim, pMSG, pMSGSO, 
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                NProc, MSG_Buffer, NMSGrm, MPICOMM
     &                ,Qmask)

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_SOLVE_boxmg is a black box matrix solver. It takes a
C     matrix defined by the user from a given (fine) grid and constructs
C     coarser grids and their associated coefficient matrices. BMG3D
C     then performs the multigrid algorithm and returns the solution
C     vec- tor. The matrix is set by a user written subroutine. The
C     difference stencil at (i,j,k) spans three planes. The k-plane has
C     difference stencil
C
C                        xnw xn xne
C                        xw  xo xe
C                        xsw xs xse.
C
C     The (k+1)-plane has difference stencil
C
C                        xtnw xtn xtne
C                        xtw  xt  xte
C                        xtsw xts xtse.
C
C     The(k-1)-plane has difference stencil
C
C                        xbnw xbn xbne
C                        xbw  xb  xbe
C                        xbsw xbs xbse.
C
C     Here 
C           xo   = SO(i,j,k,ko)
C           xw   = - SO(i,j,k,kw)
C           xs   = - SO(i,j,k,ks)
C           xb   = - SO(i,j,k,kb)
C           xsw  = - SO(i,j,k,ksw)
C           xnw  = - SO(i,j+1,k,knw)
C           xbw  = - SO(i,j,k,kbw)
C           xbnw = - SO(i,j+1,k,kbnw)
C           xbn  = - SO(i,j+1,k,kbn)
C           xbne = - SO(i+1,j+1,k,kbne)
C           xbe  = - SO(i+1,j,k,kbe)
C           xbse = - SO(i+1,j,k,kbse)
C           xbs  = - SO(i,j,k,kbs)
C           xbsw = - SO(i,j,k,kbsw)
C
C     The remaining coefficents:
C
C     xn, xne,xe,xse, xt, xtw, xtnw, xtn, xtne, xte, xtse, xts, and xtsw 
C
C     are specified by symmetry. In the case of a seven point operator,
C     only xo, xw, xs, and xb need to be set. The difference scheme is
C     assumed to be positive definite.
C     
C     The fictitious points 
C
C         (1,j,k), j=1,jj, k=1,kk, (ii,j,k), j=1,jj, k=1,kk, 
C         (i,1,k), i=1,ii, k=1,kk, (i,jj,k), i=1,ii, k=1,kk, 
C         (i,j,1), i=1,ii, j=1,jj, (i,j,kk), i=1,ii, j=1,jj,
C
C     are assumed for ease of programming. However, it is also assumed
C     that the user has set any difference coefficient referring
C     to these points to zero. For example, if the above difference
C     template is centered at (2,2), then xb, xbw, xbnw, xbn, xbne, xbe,
C     xbse, xbs xnw, xw, xsw, xs, xse, xtnw, xtw, xtsw, xts, xtse
C     would be assumed zero.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     NXM      x-dimension of the grid, excluding fictitious points
C     NYM      y-dimension of the grid, excluding fictitious points
C     NZM      z-dimension of the grid, excluding fictitious points
C
C     TOL      Convergence tolerance.
C
C     ISTOP    Indicates the stopping criteria. At present there
C              is a relative residual test
C
C                ISTOP = BMG_STOP_REL_RES_L2 
C
C                => iteration has converged if the l2-norm of the
C                   current residual, divided by the l2-norm of the 
C                   initial residual is less than TOL.
C
C              and an absolute residual test
C  
C                ISTOP = BMG_STOP_REL_ABS_L2
C
C                => iteration has converged if the l2-norm of the 
C                  current residual is less than TOL.
C
C              NOTE: if the selected test fails and the number of
C              iterations is greater than ABS(ISTRT) then the iteration
C              has failed => set TOL=-RES_L2 and return.
C
C
C     IFD      Indicator for difference scheme. IFD=1 means a seven
C              point operator. I.e., in the template given above, xnw
C              xne, xse, xsw, xbw, xbnw, xbn, xbne, xbe, xbse, xbs, xbsw,
C              xtw, xtnw, xtn, xtne, xte, xtse, xts, xtsw.
C              are assumed zero. IFD.ne.1. means a 27 point
C              operator. (Note that 27 point operators will be
C              generated on the coarser grids even in the case of a seven
C              point operator.
C
C     IU       Number of relaxation sweeps to be performed
C              on a coarse grid before interpolation to a fine grid.
C
C     ID       Number of relaxation sweeps to be performed on a fine
C              grid before the problem is transferred to a coarse grid.
C
C     IM       Unused
C
C     ISTRT    Indicator for number of multigrid cycles to be performed
C              if ISTRT.gt.0 the algorithm will begin on the coarsest grid
C              with a full multigrid cycle (however,
C              without cubic interpolation). It will continue
C              cycling until ISTRT cycles have been performed or until
C              the error criterion is satisfied. If ISTRT.lt.0 the
C              algorithm will begin on the finest grid and will perform
C              -ISTRT cycles unless the error criterion is satisfied
C              first.
C
C     ISETUP   Indicator for whether or not to skip initial setup,
C              initial guess for Q, computation of pointers, and
C              generation of coefficients on the coarser grids.
C            
C                ISETUP=0  => perform setup, including,
C                             - compute workspace pointers
C                             - compute coarse-grid operators
C                             - compute interpolation transfer operators
C
C              If ISETUP.NE.0 then we have the following cases:
C         
C                ISETUP=1  => skip computation of workspace pointers
C                             - compute coarse-grid operators
C                             - compute interpolation operators      
C                ISETUP=2  => skip all setup 
C
C              It is important to note the interaction that ISETUP 
C              has with ISTRT.  Specifically, with ISETUP.NE.0 we have
C            
C                 ISTRT > 0 => FMG so residual is coarsened appropriately
C                 ISTRT < 0 => n-cycle so nothing extra to do here.
C
C              An example of usage would be a constant coefficient 
C              time-dependent problem where ISETUP=0 on the first time step 
C              and ISETUP=2 on subsequent time steps.
C
C     BMG_IOFLAG   Logical array of I/O and debugging flags. Admittedly
C                  this is not well documented yet, but look in the include
C                  file "BMG_parameters.h" for hints (the names are verbose.)
C
C     IRELAX   Indicator for the relaxation algorithm:
C     
C              IRELAX=1  => Colored Gauss-Seidel relaxation:
C                        -  red-black on the finest grid for 5-point stencils
C                        -  four color otherwise
C
C              IRELAX=5  => alternating red-black plane Gauss-Seidel relaxation
C                        -  ordering is xy planes, yz planes, xz planes
C                        -  for symmetric relaxation this is reversed
C                           as necessary when moving "up" or "down" levels
C
C     IRELAX_SYM  Symmetric cycles require that the relaxation ordering
C                 be reversed between the restriction and interpolation
C                 stage of the cycle.   This is important when BMG3D
C                 is being used as a preconditioner.
C              
C                    IRELAX_SYM = BMG_RELAX_NONSYM
C                    => use the stated relaxation ordering (nonsymmetric)
C
C                    IRELAX_SYM = BMG_RELAX_SYM
C                    => make the cycle symmetric.
C
C                 as defined in BMG_parameters.h
C
C     ISEP     Indicator for separable equation such that the SO array can be
C              reduced to 3 1-D arrays in the interior, saving memory bandwidth
C
C               ISEP = BMG_NONSEPARABLE (0): Always use 3d SO arrays
C               ISEP = BMG_SEPARABLE (1): only use 3d SO arrays on boundaries.
C
C
C     IVW      Indicator for type of cycle. IVW=1 means v-cycles will
C              be performed. IVW=2 means w-cycles will be performed, etc.
C
C     MCYCL    Currently disabled!
C
C     NFm      Maximum storage for a vector on all grids. This should be
C              computed by BMG3_SymStd_SETUP_space and is checked in boxmg.
C
C     NCm      Maximum storage for a vector on all coarse grids.  This should
C              be computed in BMG3_SymStd_SETUP_space and is checked in boxmg.
C
C     NSOm     See INPUT/OUTPUT
C
C     NSORm    See INPUT/OUTPUT
C
C     NCIm     See INPUT/OUTPUT
C     
C     NOGm     Maximum number of grids that can be supported by the 
C              pointer array IGRD. It should be set in calling program.
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     SO       User defined real array that contains the coefficient
C              matrix. See above description for format.  If SETUP is
C              specified then coarse-grid operators are created and returned. 
C
C     NSOm     Dimension of SO, the calling program should have computed
C              this with a call to BMG2_SymStd_SETUP_space.  Clearly,
C            
C              IFD .EQ. 1  =>    NSOm = 4*NFm + 10*NCm 
C              IFD .NE. 1  =>    NSOm = 14*NFm 
C
C              BMG3_SymStd_SOLVE_boxmg verifies that NSOm is large enough.
C
C     CI       CI is a real array which should be dimensioned
C              to (NCIm) in calling program. It is used to contain
C              the interpolation coefficients.
C
C     NCIm     Dimension of CI, set in calling program.  NCIm = 26*NCm is 
C              enough.  BOXMG checks to see if NCIm is large enough.
C
C     QF       The user defined array that contains the right hand side.
C              It must be dimensioned to at least NFMAX.
C
C     Q        the user defined array that contains the solution vector.
C              It must be dimensioned to at least NFMAX.
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C
C =======================================================================
C  --------------------
C   WORK ARRAYS:
C  --------------------
C
C     ABD      User declared two dimensional real array, which
C              is used to store the coefficient matrix for the coarsest
C              grid. It is then used by the LAPACK routines. It
C              should be dimensioned to (NCBWm,NCUm).
C
C     BBD      User declared real array of dimension NCUm for use in the
C              the LAPACK routines
C
C     NCBWm    Maximum first subscript of ABD, which needs to be > or =
C              (x-dimension+1)*(y-dimension+2) on coarsest grid. 
C
C     NCUm     Maximum second subscript of ABD, which needs to be > or =
C              (x-dimension+1)*(y-dimension+1)*(z-dimension+1) on
C              coarsest grid.
C
C     IGRD     A work array. IGRD should be dimensioned to
C              IGRD(NOGm,NBMG_pIGRD) in the calling program.
C
C     SOR      SOR is a real array, which should be dimensioned to
C              (NSORm) in the calling program. It is used to store
C              residuals, reciprocals of SO(.,1) if IRELAX=1, and lu
C              decompositions if IRELAX.gt.1.
C
C     NSORm    Dimension of SOR, set in calling program. NSORm=2*NFMAX
C              is enough for IRELAX=1 or 2. For IRELAX=3 or 4, 4*NFMAX
C              is enough. 
C
C
C =======================================================================
C  --------------------
C   REFERENCES:
C  --------------------
C
C     Dendy, J. E. Jr., "Two Multigrid Methods for Three Dimensional
C     Problems with Discountinuos and Anisotropic Coefficients", SIAM
C     Journal of Scientific and Satatistical Computing, Vol. 8, No. 2,
C     September 1987
C     
C ==========================================================================
C  --------------------
C   Error Codes:
C  --------------------
C
C     -------------------------------------------------------------------
C     CODE    DESCRIPTION                           ORIGIN
C     -------------------------------------------------------------------
C      1     initial residual has become zero     SOLVE_boxmg
C      2     NOG = 0                              SETUP_parts
C      3     NOG < 0                              SETUP_parts
C      4     NStncl out of range                  SETUP_cg_LU
C      5     Cholesky decomposition failed, look
C            in BMG_iPARMS(id_BMG3_Ext_Err_Code)
C            for the return code of the LAPACK
C            routine                              SETUP_cg_LU
C      6    NOG = 1                               SETUP_relax
C      7    NOG = 0                               SETUP_relax
C      8    NOG < 0                               SETUP_relax
C      9    fatal setup error                     SETUP_PtrGrid
C     10    IRELAX out of range                   SETUP_PtrGrid
C     11    min coarse grid dim < 3 ( NXYZc )     SETUP_space
C     12    min coarse grid dim < 3 ( NXYc )      SETUP_space
C     13    IRELAX out of range                   SETUP_space
C     14    computed number of grids too small    SETUP_space
C     15    computed number of grids < 1          SETUP_space
C     16    fatal error (check output)            SETUP_PtrWork
C     17    memory allocation mode unspecified    SETUP_PtrWork
C     18    KCF = KC                              ncycle
C     19    updown out of range                   updown
C     20    coarse grid solve failed, look in
C           BMG_iPARMS(id_BMG3_Ext_Err_Code)
C           for the return code of the LAPACK
C           routine                               SOLVE_cg
C     21    inconsistent # of stencil entries     COPY_SO_xy
C     22    inconsistent # of stencil entries     COPY_SO_yz
C     23    inconsistent # of stencil entries     COPY_SO_xz
C     24    invalid cg_operator type              SETUP_parts
C     25    invalid cg_operator construction      SETUP_parts
C
C ==========================================================================

      IMPLICIT NONE

C ------------------------------------------------
C     Includes
C
      INCLUDE 'geom_param_fort.h'
      INCLUDE 'mpi_param_fort.h'
      INCLUDE 'MSG.h'

#include    "BMG_constants_f77.h"
      INCLUDE 'BMG_stencils_f77.h'
      INCLUDE 'BMG_workspace_f77.h'
      INCLUDE 'BMG_parameters_f77.h'

C ------------------------------------------------
C     Argument Declarations
C
      INTEGER  NLx, NLy, NLz, NGx, NGy, NGz, iGs, jGs, kGs,
     &         NFm, NCm, NSOm, NSORm, NCIm, NCBWm, NCUm, NOGm, NOG,
     &         NMSGim, NMSGrm, NProc, MPICOMM, NBMG_MSG_iGRID,
     &         NBMG_iWORK_CSm, NBMG_iWORK_PLm,
     &         NBMG_rWORK_CSm, NBMG_rWORK_PLm

      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER  iWorkMSG(NMSGim),
     &         pMSG(NBMG_pMSG,NOGm),
     &         pMSGSO(NBMG_pMSG,NOGm)
      REAL*RKIND   MSG_Buffer(NMSGrm)

      INTEGER  BMG_iWORK_PL(NBMG_iWORK_PLm)
      REAL*RKIND   BMG_rWORK_PL(NBMG_rWORK_PLm)

      INTEGER  BMG_iWORK_CS(NBMG_iWORK_CSm)
      REAL*RKIND   BMG_rWORK_CS(NBMG_rWORK_CSm)

      INTEGER  IGRD(NOGm,NBMG_pIGRD)
      REAL*RKIND CI(NCIm),Q(NFm),QF(NFm),RES(NFm),SO(NSOm),SOR(NSORm)
!      REAL*RKIND, optional :: Qmask(NLx,Nly,Nlz)
      REAL*RKIND :: Qmask(NLx,Nly,Nlz)

C ------------------------------------------------
C     Local Declarations
C
      INTEGER i, iBC, id, ifd, iil, iilc, iilf, im, irelax, isep,
     &        IRELAX_SYM, ISETUP, ISTOP, istrt, istrt2, iu, ivw,
     &        jjl, jjlc, jjlf, k, kc, kfmg, kf, kkl, kklc, kklf,
     &        mcyc, mcycl, ncyc, NStncl, NMSGi, NMSGr, p_ABD, p_BBD,
     &        p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC,
     &        IIG, IIGC, IIGF, JJG, JJGC, JJGF, KKG, KKGC, KKGF,
     &        NC, NCBW, NCI, NCU, NF, NSO, NSOR, NXYZc,
     &        NBMG_iWORK_CS, NBMG_iWORK_PL, NBMG_rWORK_CS, 
     &        NBMG_rWORK_PL


      INTEGER  iGs_c, iGs_f, jGs_c, jGs_f, kGs_c, kGs_f

      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK,
     &         NProcI, NProcJ, NProcK,
     &         p_NLx, p_NLy, p_NLz,
     &         p_NLx_kg, p_NLy_kg, p_NLz_kg,
     &         p_ProcCoord, p_ProcGrid, pSI_MSG

      INTEGER  BMG_MSG_MyProc, MPI_IERR

      REAL*8  REL_RES_L2, RES_L2, RES_L2_0
      REAL*8 T, T1, T2, TOL
      REAL*8  TT1, TT2
      LOGICAL NCYCLE_FLAG, INIT_CHECK

C ==========================================================================

      IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
         RETURN
      END IF

C -------------------
C     Zero Times:
C -------------------

      t=dZERO
      t1=dZERO
      t2=dZERO

      IF ( BMG_IOFLAG(iBMG3_BUG_PARAMETERS) ) THEN
         CALL BMG3_SymStd_DUMP_parms( 
     &                    BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &                    )
      ENDIF

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------
C --------------------------------------------
C     MPI Parallel Data Decomposition:
C --------------------------------------------
      
      !
      ! Processor Grid:
      !
      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_NLx       = BMG_MSG_pGRID(ip_BMG_MSG_NLx_Grid)
      p_NLy       = BMG_MSG_pGRID(ip_BMG_MSG_NLy_Grid)
      p_NLz       = BMG_MSG_pGRID(ip_BMG_MSG_NLz_Grid)
      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)
      
C ----------------------------------------
C     MSG workspace pointers:
C ----------------------------------------

      p_NLx_kg = 1
      p_NLy_kg = p_NLx_kg + NProcI*NOGm
      p_NLz_kg = p_NLy_kg + NProcJ*NOGm

      pSI_MSG  = p_NLz_kg + NProcK*NOGm  !!NB: initalization of pSI_MSG

C ----------------------------------------
C     Multigrid Cycle parameters:
C ----------------------------------------

      IFD    = BMG_iPARMS(id_BMG3_STENCIL)
      IBC    = BMG_iPARMS(id_BMG3_BC)
      ISETUP = BMG_iPARMS(id_BMG3_SETUP)

      IRELAX = BMG_iPARMS(id_BMG3_RELAX)
      IRELAX_SYM = BMG_iPARMS(id_BMG3_RELAX_SYM )

      ISEP = BMG_iPARMS(id_BMG3_SEPARABLE)

      ID = BMG_iPARMS(id_BMG3_NRELAX_DOWN)
      IU = BMG_iPARMS(id_BMG3_NRELAX_UP)
      IM = BMG_iPARMS(id_BMG3_NRELAX_FG)

      IF (BMG_iPARMS(id_BMG3_CYCLE_CLASS).EQ.BMG_N_CYCLE) THEN
         ISTRT = - BMG_iPARMS(id_BMG3_MAX_ITERS)
      ELSEIF (BMG_iPARMS(id_BMG3_CYCLE_CLASS).EQ.BMG_FMG_CYCLE) THEN 
         ISTRT = BMG_iPARMS(id_BMG3_MAX_ITERS)
      ENDIF

      IVW   = BMG_iPARMS(id_BMG3_NCYCLE_TYPE)
      MCYCL = BMG_iPARMS(id_BMG3_FMG_NNCYCLE)

      ISTOP = BMG_iPARMS(id_BMG3_STOP_TEST)
      TOL   = BMG_rPARMS(id_BMG3_STOP_TOL)

C ================================================

      IF (BMG_iPARMS(id_BMG2_CYCLE_CLASS).EQ.BMG_N_CYCLE) THEN
         ISTRT2 = - BMG_iPARMS(id_BMG2_MAX_ITERS)
      ELSEIF (BMG_iPARMS(id_BMG2_CYCLE_CLASS).EQ.BMG_FMG_CYCLE) THEN 
         ISTRT2 = BMG_iPARMS(id_BMG2_MAX_ITERS)
      ENDIF

C ================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C ------------------------------------------------------------------- 
C     Enable MSG communications checking to make sure all initialized
C -------------------------------------------------------------------

      CALL MPI_INITIALIZED(INIT_CHECK, MPI_IERR)
      IF( .NOT. INIT_CHECK ) THEN
         PRINT *,' MPI NOT INITIALIZED!'
         RETURN
      ENDIF
      
      IF( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR. 
     &    ISETUP.EQ.BMG_SETUP_only ) THEN

         CALL MSG_set_comm_parent(MPICOMM)
         CALL MSG_enable(BMG_MSG_MyProc,NProc)
         BMG_MSG_MyProc = MSG_myproc()
         NProc  = MSG_nproc()
            
      ELSE
         
         IF( MPICOMM .EQ. MSG_COMM_PARENT) THEN
            BMG_MSG_MyProc = MSG_myproc()
            NProc  = MSG_nproc()
         ELSE
            IF( MyProc.EQ.1 ) PRINT *,' MSG NOT INITIALIZED!'
         ENDIF

      ENDIF

C ---------------------------------
C     Initialize timers
C ---------------------------------

      T  = dZERO
      T1 = dZERO
      T2 = dZERO

C ----------------------------------
C     Sanity checks!
C ----------------------------------

      !!
      !!   We need to add way more some sanity checks here  !!!
      !!   Maybe a subroutine to check all parameter values??
      !!   

      !
      ! Processor indexing
      !
      IF ( BMG_MSG_MyProc.NE.BMG_MSG_iGRID(id_BMG_MSG_MyProc) ) THEN 
         WRITE(*,*) '**** ERROR: BMG3_SymStd_SOLVE_boxmg.f ... '
         WRITE(*,*) ' MyProc =         ', MyProc
         WRITE(*,*) ' BMG_MSG_MyProc = ', BMG_MSG_MyProc
         STOP
      ENDIF


C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> BEGIN:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      IF ( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR. 
     &     ISETUP.EQ.BMG_SETUP_only .OR. 
     &     ISETUP.EQ.BMG_SETUP_ptrs ) THEN

         TT1 = MPI_Wtime()

         CALL BMG3_SymStd_SETUP_PtrGrid( 
     &             NLx, NLy, NLz, NGx, NGy, NGz,
     &             iGs, jGs, kGs, 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,    
     &             MyProc, NProc, NProcI, NProcJ, NProcK, 
     &             MPICOMM, BMG_iWORK_PL, BMG_iWORK_CS, 
     &             NOGm, NFm, NCm, NSOm, NSORm, NCIm,
     &             NBMG_iWORK_PLm, NBMG_rWORK_PLm,
     &             NBMG_iWORK_CSm, NBMG_rWORK_CSm, NMSGim, NMSGrm,
     &             NOG, NF, NC, NSO, NSOR, NCI, IGRD,
     &             NBMG_iWORK_PL, NBMG_rWORK_PL,
     &             NBMG_iWORK_CS, NBMG_rWORK_CS, NMSGi, NMSGr
     &             )
         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF


         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_SETUP_PTR_GRID) =
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_PTR_GRID) + TT2 - TT1

         ! Output      
         IF( BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE) ) THEN
            WRITE (*,260) 
     &            'Storage for a vector on all grids = ', NF
         ENDIF

         ! Output
         IF( BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE) ) THEN
            WRITE (*,270) 
     &            'Storage for a vector on all coarse grids = ', NC
         ENDIF

         !
         ! Setup of the MSG arrays (pSI_MSG initialized above).
         !
         TT1 = MPI_Wtime()


         CALL BMG3_SymStd_SETUP_MSG( 
     &             pMSG, pMSGSO, iWorkMSG, NMSGi, pSI_MSG, 
     &             IGRD, NOG, NOGm, NProc, MyProc, 
     &             iWorkMSG(p_NLx_kg), iWorkMSG(p_NLy_kg),
     &             iWorkMSG(p_NLz_kg), 
     &             BMG_MSG_iGRID(p_NLx), BMG_MSG_iGRID(p_NLy),
     &             BMG_MSG_iGRID(p_NLz),
     &             BMG_MSG_iGRID(p_ProcGrid),
     &             NProcI, NProcJ, NProcK,
     &             MPICOMM 
     &             )

         TT2 = MPI_Wtime()

         BMG_rPARMS(id_BMG3_TIME_SETUP_MSG) = 
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_MSG) + TT2 - TT1
     
      ELSE

         !
         !  Maximum dimensions are the actual dimensions.
         !

         NF   = NFm
         NC   = NCm
         NSO  = NSOm
         NSOR = NSORm
         NCI  = NCIm

         NBMG_iWORK_CS = NBMG_iWORK_CSm
         NBMG_rWORK_CS = NBMG_rWORK_CSm
 
         NMSGi = NMSGim
         NMSGr = NMSGrm

      ENDIF

      !
      ! Create local pointers for LU solve
      !
      IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN

         p_ABD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO)
         p_BBD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CU)
         
         NCBW  = BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW)
         NCU   = BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)

      ENDIF

      !
      !  No components just setup pointers, thanks!
      !
      IF ( ISETUP.EQ.BMG_SETUP_ptrs ) RETURN

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>> END:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>> BEGIN:  COMPONENT SETUP <<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------


 89   CONTINUE  ! >>>>>>>>> BEGIN: Component Setup

      IF ( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR. 
     &     ISETUP.EQ.BMG_SETUP_opers .OR. 
     &     ISETUP.EQ.BMG_SETUP_only ) THEN

         !  Start the timer
         t1 = MPI_Wtime()

         !
         ! Construct coarse-grid and interpolation operators.  In 
         ! addition, setup any necessary components for relaxation.
         !
         TT1 = MPI_Wtime()
         
         CALL BMG3_SymStd_SETUP_parts( 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm, NOG,
     &             BMG_iWORK_PL, NBMG_iWORK_PL,
     &             BMG_rWORK_PL, NBMG_rWORK_PL,
     &             BMG_iWORK_CS, NBMG_iWORK_CS,
     &             BMG_rWORK_CS, NBMG_rWORK_CS,
     &             iWorkMSG, NMSGi, pSI_MSG, pMSG, pMSGSO,
     &             MSG_Buffer, NMSGr,
     &             MyProcI, MyProcJ, MyProcK,
     &             NProcI, NProcJ, NProcK, NProc,
     &             BMG_MSG_iGRID(p_ProcGrid),
     &             BMG_MSG_iGRID(p_ProcCoord),
     &             iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &             iWorkMSG(p_NLx_kg), iWorkMSG(p_NLy_kg),
     &             iWorkMSG(p_NLz_kg),
     &             MPICOMM
     &             )
         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_SETUP_PARTS) =
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_PARTS) + TT2 - TT1

         ! Compute the setup time
         t2 = MPI_Wtime()
         t=t+t2-t1
         BMG_rPARMS(id_BMG3_TIME_SETUP_TOTAL) =
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_TOTAL) + t2 - t1

         
         ! Output the setup time
         IF( BMG_IOFLAG(iBMG3_OUT_TIME_SETUP) )  THEN
            WRITE (*,220) '(3D) SETUP TIME =', T
         ENDIF

      ENDIF

      !
      !  No solve just the setup, thanks!  
      !
      IF (ISETUP.EQ.BMG_SETUP_only) RETURN

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  COMPONENT SETUP <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN: MULTIGRID CYCLING <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

 95   CONTINUE

      ! Start the timer
      t1 = MPI_Wtime()

C ----------------------------------------------------------------
C     Initial residual:
C ----------------------------------------------------------------

      !
      ! Get pointers for the finest grid (GRID LEVEL = NOG)
      !
      CALL BMG3_SymStd_GET_pointers( 
     &          NOG, IGRD, NOGm, 
     &          p_U, p_SO, p_SOR, p_CI, 
     &          IIL, JJL, KKL, IIG, JJG, KKG
     &          )

C      FIXME: The arguments passed in should be correct, but this 
C             is a reasonable place to double check
C
C      iGs = IGRD(NOG,idL_BMG_Icoord)
C      jGs = IGRD(NOG,idL_BMG_Jcoord)
C      kGs = IGRD(NOG,idL_BMG_Kcoord)

      !
      ! Set NStncl for the finest grid
      !
      IF ( IFD.NE.1 ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF
     
      !
      ! if requested, sync the initial guess
      !
      IF ( BMG_iPARMS(id_BMG3_SYNC_INITIAL_GUESS) 
     &     .eq. BMG_SYNC_INITIAL_GUESS ) THEN
         
         CALL BMG3_SymStd_UTILS_update_ghosts(
     &             NOG, Q(p_U), IIL, JJL, KKL, iWorkMSG, pMSG,
     &             MSG_Buffer, MPICOMM
     &             )

      END IF

      CALL BMG3_SymStd_residual(
     &          NOG, NOGm, IFD,
     &          Q(p_U), QF(p_U), SO(p_SO), RES(p_U), 
     &          IIL, JJL, KKL, NStncl,
     &          iWorkMSG, NMSGi, pMSG, 
     &          MSG_Buffer, NMSGr, MPICOMM
     &          )

      CALL BMG3_SymStd_UTILS_norm_l2( 
     &           RES(p_U), IIL, JJL, KKL, RES_L2_0, MPICOMM 
     &           )

      !
      !  Output
      !
      IF ( BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
         WRITE (*,405) '*** THE INITIAL RESIDUAL LGT (l2-NORM) = ', RES_L2_0
      ENDIF

      IF ( RES_L2_0 .EQ. rZERO ) THEN
         IF ( BMG_IOFLAG(iBMG3_WARN_ZERO_RESIDUAL) ) THEN
            WRITE(*,505) '*** WARNING: bmg3d.f '
            WRITE(*,505) '    Zero initial residual !!!! '
            WRITE(*,505) '    D''ooh! Initial guess is the solution!'
         ENDIF
         RETURN
      ENDIF

C ----------------------------------------------------------------
C     Direct solve in the degenerate case of one grid!
C ----------------------------------------------------------------

      IF ( NOG.EQ.1 ) THEN

         !
         ! Solve on the coarsest grid
         !
         IF ( MyProc.EQ.1 ) THEN 
            WRITE(*,*) '*** Performing direct solve on the fine grid!'
         ENDIF
         !

         IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN
            !
            TT1 = MPI_Wtime()
            !
            CALL BMG3_SymStd_SOLVE_cg_LU( 
     &                q(p_U), qf(p_U), IIL, JJL, KKL,
     &                BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD), 
     &                NCBW, NCU, NOGm,
     &                NProcI, NProcJ, NProcK, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM, BMG_iPARMS
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ELSE
            !
            TT1 = MPI_Wtime()
            !
            
            CALL BMG3_SymStd_SOLVE_cg_boxmg( 
     &                Q(p_U), QF(p_U), NOG,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                IIL, JJL, KKL, IIG, JJG, KKG,
     &                iGs, jGs, kGs,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                NProcI, NProcJ, NProcK, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ENDIF

         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         BMG_rPARMS(id_BMG3_TIME_SOLVE_CG)  
     &       = BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) + TT2 - TT1

         !
         ! Compute the final residual
         !
         CALL BMG3_SymStd_residual( 
     &             NOG, NOG, IFD,
     &             Q(p_U), QF(p_U), SO(p_SO), RES(p_U),
     &             IIL, JJL, KKL, NStncl, 
     &             iWorkMSG, NMSGi, pMSG, 
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )
         !
         CALL BMG3_SymStd_UTILS_norm_l2( 
     &             RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM 
     &             )
         !
         IF( BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) .AND. MyProc.EQ.1 ) THEN 
            WRITE (*,405) '*** THE FINAL RESIDUAL (l2-NORM) = ', RES_L2
         ENDIF
         !
         TOL = RES_L2
         ! Jump to final I/O
         GOTO 150
         !
      ENDIF



C --------------------------------
C     Cycle Preliminaries:
C -------------------------------- 

      KF = NOG                  ! finest grid index
      KC = 1                    ! coarsest grid index
      KFMG = 1                  ! FMG current finest grid
      MCYC = ABS(ISTRT)         ! Maximum number of multigrid cycles
      NCYC = 1                  ! multigrid cycle counter
      NCYCLE_FLAG = .TRUE.      ! default to do n-cycles

C ----------------------------------------------------------------
C     Perform multigrid F-cycle:
C     (if necessary, then continue with n-cycles)
C ----------------------------------------------------------------

      IF ( ISTRT.GT.0 ) THEN

         ! Coarsen the right hand side
         DO k=KF-1, 1, -1

            CALL BMG3_SymStd_GET_pointers( k+1, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, 
     &                IILF, JJLF, KKLF, IIGF, JJGF, KKGF 
     &                )
            iGs_f = IGRD(k+1,idL_BMG_Icoord)
            jGs_f = IGRD(k+1,idL_BMG_Jcoord)
            kGs_f = IGRD(k+1,idL_BMG_Kcoord)

            CALL BMG3_SymStd_GET_pointers( k, IGRD, NOGm,
     &                p_UC, p_SOC, p_SORC, p_CIC, 
     &                IILC, JJLC, KKLC, IIGC, JJGC, KKGC  
     &                )
            iGs_c = IGRD(k,idL_BMG_Icoord)
            jGs_c = IGRD(k,idL_BMG_Jcoord)
            kGs_c = IGRD(k,idL_BMG_Kcoord)

            ! restrict the right hand side 

            TT2 = MPI_Wtime()

            CALL BMG3_SymStd_restrict( 
     &                K+1, K,  
     &                QF(p_U), QF(p_UC), CI(p_CIC),
     &                IILF, JJLF, KKLF, IILC, JJLC, KKLC,
     &                iGs_f, jGs_f, kGs_f
     &                )

            TT1 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_restrict) =
     &           BMG_rPARMS(id_BMG3_TIME_restrict) + TT2 - TT1

            ! Copy right hand side into RES(i,j,k) for interpolation
            DO i=0, IILF*JJLF*KKLF-1
               RES(p_U+i)=QF(p_U+i)
            END DO
         END DO


         ! FMG current finest grid index is set to the coarsest grid
         KFMG=KC     
         ! FMG: solve exactly on the coarsest grid
         CALL BMG3_SymStd_GET_pointers(
     &             KFMG, IGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI, 
     &             IIL, JJL, KKL, IIG, JJG, KKG 
     &             )
         iGs_c = IGRD(KFMG,idL_BMG_Icoord)
         jGs_c = IGRD(KFMG,idL_BMG_Jcoord)
         kGs_c = IGRD(KFMG,idL_BMG_Kcoord)
         
         IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU ) THEN
            !
            TT1 = MPI_Wtime()
            !
            CALL BMG3_SymStd_SOLVE_cg_LU( 
     &                q(p_U), qf(p_U), IIL, JJL, KKL, 
     &                BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD), 
     &                NCBW, NCU, NOGm,
     &                NProcI, NProcJ, NProcK, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM, BMG_iPARMS 
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ELSE
            !
            TT1 = MPI_Wtime()
            !

            CALL BMG3_SymStd_SOLVE_cg_boxmg( 
     &                Q(p_U), QF(p_U), NOG,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                IIL, JJL, KKL, IIG, JJG, KKG,
     &                iGs_c, jGs_c, kGs_c,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                NProcI, NProcJ, NProcK, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ENDIF

         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         BMG_rPARMS(id_BMG3_TIME_SOLVE_CG)  
     &       = BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) + TT2 - TT1
         
         RES_L2=dZERO
         IF ( BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
            WRITE (*,230) KFMG, RES_L2
         ENDIF

 100     CONTINUE    ! >>>>>>>> LOOP BOUNDARY: multigrid F-cycle


            KFMG=KFMG+1
            ! 
            ! FMG: Interpolate
            !
            CALL BMG3_SymStd_GET_pointers(
     &                KFMG, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, IIL, JJL, KKL,
     &                IIG, JJG, KKG 
     &                )
            iGs_f = IGRD(KFMG,idL_BMG_Icoord)
            jGs_f = IGRD(KFMG,idL_BMG_Jcoord)
            kGs_f = IGRD(KFMG,idL_BMG_Kcoord)

            CALL BMG3_SymStd_GET_pointers(
     &                KFMG-1, IGRD, NOGm,
     &                p_UC, p_SOC, p_SORC, p_CIC, IILC, JJLC, KKLC,
     &                IIGC, JJGC, KKGC  
     &                )
            iGs_c = IGRD(KFMG-1,idL_BMG_Icoord)
            jGs_c = IGRD(KFMG-1,idL_BMG_Jcoord)
            kGs_c = IGRD(KFMG-1,idL_BMG_Kcoord)

            !
            ! Set NStncl for grid KFMG
            !
            IF ( KFMG.LT.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
               NStncl=14
            ELSE
               NStncl=4
            ENDIF
            !


            TT1 = MPI_Wtime()

            IF ( KFMG < NOG ) THEN
            CALL BMG3_SymStd_interp_add( 
     &                KFMG-1, KFMG, NOG,
     &                Q(p_U) ,Q(p_UC), RES(p_U),
     &                SO(p_SO), NStncl, CI(p_CIC),
     &                IILC, JJLC, KKLC,
     &                IIL, JJL, KKL,
     &                iGs_f, jGs_f, kGs_f,
     &                iWorkMSG, NMSGi, pMSG,
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            
            ELSE ! special version for fine grid (checks mask array)
            CALL commas_SymStd_interp_add( 
     &                KFMG-1, KFMG, NOG,
     &                Q(p_U) ,Q(p_UC), RES(p_U),
     &                SO(p_SO), NStncl, CI(p_CIC),
     &                IILC, JJLC, KKLC,
     &                IIL, JJL, KKL,
     &                iGs_f, jGs_f, kGs_f,
     &                iWorkMSG, NMSGi, pMSG, Qmask,
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            ENDIF


            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_interp_add) =
     &           BMG_rPARMS(id_BMG3_TIME_interp_add) + TT2 - TT1
            !
            ! FMG; Perform an n-cycle 
            !
            IF ( BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN 
               WRITE(*,410) '*** Performing n-cycle:  Coarsest =', KC,
     &                      'Finest =', KFMG
            ENDIF


            CALL commas_SymStd_ncycle( 
     &                KC, KFMG, KF, 
     &                IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm,
     &                SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID
     &                ,Qmask)


            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            IF ( KFMG.EQ.KF ) THEN
               !
               !  Compute the final residual
               !
               CALL BMG3_SymStd_residual(
     &                   KF, KF, IFD,
     &                   Q(p_U), QF(p_U), SO(p_SO), RES(p_U),
     &                   IIL, JJL, KKL, NStncl,
     &                   iWorkMSG, NMSGi, pMSG, 
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )
               CALL BMG3_SymStd_UTILS_norm_l2( 
     &                   RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM
     &                   )

               ! Note this is still vulnerable to very small RES_L2_0
               IF ( RES_L2_0.EQ.dZERO ) THEN
                  IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR).AND.MyProc.EQ.1)
     &                 THEN
                    WRITE(*,505)
     &                    '*** FATAL ERROR: BMG3_SymStd_SOLVE_boxmg.f '
                    WRITE(*,505)'    Initial residual has become zero!'
                  END IF
                  
                  CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,1)
                  RETURN
                  
               ELSE
                  REL_RES_L2 = RES_L2/RES_L2_0
               ENDIF             

               !
               !  Output 
               !
               IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )
     &            .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
                  WRITE (*,400) '*** ITERATION =', NCYC, 
     &                          '*** RELATIVE RESIDUAL = ', REL_RES_L2
               ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 ) 
     &                 .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
                  WRITE (*,400) '*** ITERATION =', NCYC, 
     &                          '*** ABSOLUTE RESIDUAL = ', RES_L2
               ENDIF
               !
               !  Check Convergence
               !
               IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 ) 
     &            .AND. ( REL_RES_L2.LT.TOL )     ) THEN
                  ! FMG: converged
                  !      - set TOL and skip n-cycles
                  TOL = REL_RES_L2
                  NCYCLE_FLAG=.FALSE.
               ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 ) 
     &                 .AND. ( RES_L2.LT.TOL )         ) THEN
                  ! FMG: converged
                  !      - set TOL and skip n-cycles
                  TOL = RES_L2
                  NCYCLE_FLAG=.FALSE.
               ELSE IF ( NCYC.GE.MCYC ) THEN
                  ! FMG: cycle limit reached 
                  !      - set TOL and skip n-cycles
                  IF ( ISTOP.EQ.BMG_STOP_REL_RES_L2 ) THEN
                     TOL = -REL_RES_L2
                  ELSE
                     TOL = -RES_L2
                  ENDIF
                  NCYCLE_FLAG=.FALSE.
               ELSE 
                  ! FMG: didn't converge
                  !      - increase counter and continue with n-cycles
                  NCYC=NCYC+1
               ENDIF
               !
            ELSE
               ! FMG: Continue the F-cycle
               GOTO 100
            ENDIF
            
 120     CONTINUE    ! <<<<<<<< LOOP BOUNDARY: multigrid F-cycle

      ENDIF  ! ISTRT > 0



C ----------------------------------------------------------------
C     Perform multigrid n-cycle(s):
C ----------------------------------------------------------------

      IF ( NCYCLE_FLAG ) THEN
         
         CALL BMG3_SymStd_GET_pointers(
     &             KF, IGRD, NOGm, p_U, p_SO, 
     &             p_SOR, p_CI, IIL, JJL, KKL,
     &             IIG, JJG, KKG
     &             )
         

 130     CONTINUE   ! >>>>>>>> LOOP BOUNDARY: multigrid n-cycles



            !
            ! Perform a multigrid n-cycle
            !
            CALL commas_SymStd_ncycle( 
     &                KC, KF, KF, 
     &                IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm,
     &                SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr, 
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID
     &                ,Qmask)

            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            !
            !  Compute the final residual
            !
            CALL BMG3_SymStd_residual(
     &                KF, KF, IFD,
     &                Q(p_U), QF(p_U), SO(p_SO), RES(p_U), 
     &                IIL, JJL, KKL, NStncl, 
     &                iWorkMSG, NMSGi, pMSG, 
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            CALL BMG3_SymStd_UTILS_norm_l2( 
     &                RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM
     &                )

            ! Note this is still vulnerable to very small RES_L2_0
            IF ( RES_L2_0.EQ.rZERO ) THEN
               WRITE(*,505) '*** FATAL ERROR: bmg3d.f '
               WRITE(*,505) '    Initial residual has become zero!!!'
               STOP
            ELSE
               REL_RES_L2 = RES_L2/RES_L2_0
            ENDIF

            !
            !  Output 
            !
            IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )
     &           .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** ITERATION =', NCYC, 
     &                       '*** RELATIVE RESIDUAL = ', REL_RES_L2
            ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 )
     &              .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** ITERATION =', NCYC, 
     &                      '*** ABSOLUTE RESIDUAL = ', RES_L2
            ENDIF

            !
            ! Check Convergence 
            !
            IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )
     &         .AND. ( REL_RES_L2.LT.TOL )     ) THEN
               ! 
               ! n-cycles have converged in the RELATIVE RESIDUAL
               ! - set TOL and return
               !
               TOL = REL_RES_L2
            ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 ) 
     &              .AND. ( RES_L2.LT.TOL )         ) THEN
               ! 
               ! n-cycles have converged in the ABSOLUTE RESIDUAL
               ! - set TOL and return
               !
               TOL = RES_L2
            ELSE IF ( NCYC.GE.MCYC ) THEN
               !
               ! n-cycles failed to converge, cycle limit reached
               ! - set TOL and return
               !
               IF ( ISTOP.EQ.BMG_STOP_REL_RES_L2 ) THEN
                  TOL = -REL_RES_L2
               ELSE
                  TOL = -RES_L2
               ENDIF
            ELSE 
               !
               ! increase cycle count and continue n-cycling
               !
               NCYC=NCYC+1
               GOTO 130
            ENDIF

 140     CONTINUE    ! <<<<<<<< LOOP BOUNDARY: multigrid n-cycles

      ENDIF

 150  CONTINUE       ! <<<<<<<< FINAL I/O
      !
      !
      BMG_rPARMS(id_BMG3_STOP_TOL) = TOL
      !
      !  Compute the solve time and update the total time
      !
      t2 = MPI_Wtime()
      T=T+T2-T1
      BMG_rPARMS(id_BMG3_TIME_SOLVE_total) = 
     &     BMG_rPARMS(id_BMG3_TIME_SOLVE_total) + t2 - t1


      ! Output the multigrid cycling time
      IF ( BMG_IOFLAG(iBMG3_OUT_TIME_CYCLING) ) THEN
         WRITE(*,240) '(3D) MULTIGRID CYCLING TIME =', T2-T1

      ENDIF

      IF ( BMG_IOFLAG(iBMG3_OUT_TIME_TOTAL) ) THEN
         WRITE(*,240) '(3D) TOTAL TIME = ', T
      ENDIF

C -------------------------------------
C     Output:
C -------------------------------------         

C ==========================================================================

 220  FORMAT (/,1X,A,1X,F9.3,/)
 230  FORMAT ('(3D)LEVEL',I2,' RESIDUAL NORM= ',1P,E10.3)
 240  FORMAT (/,1X,A,1X,F12.3,/)
 260  FORMAT (/,/,2X,A,1X,I5)
 270  FORMAT (2X,A,1X,I5,/)

 400  FORMAT (1X,A,1X,I2,4X,A,1X,1P,E16.9)
 405  FORMAT (/,1X,A,1X,1P,E16.9,/)
 410  FORMAT (1X,A,1X,I2,4X,A,1X,I2)

 500  FORMAT (/,2X,A)
 505  FORMAT (5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
      
C ===========================================

      RETURN
      END

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE commas_SymStd_ncycle( 
     &                KC, KCF, NOG, 
     &                IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC,
     &                SO, NSO, SOR, NSOR, CI, NCI, IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr, 
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID
     &                ,Qmask)

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   BMG3_SymStd_ncycle performs a single multigrid n-cycle.  
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'mpif.h'

#include    "BMG_constants_f77.h"
      INCLUDE 'BMG_workspace_f77.h'
      INCLUDE 'BMG_stencils_f77.h'
      INCLUDE 'BMG_parameters_f77.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  KC, KCF, NOG, IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, ISEP,
     &         NF, NC, NSO, NSOR, NCI, NOGm, NMSGi, NMSGr,
     &         NBMG_MSG_iGRID, NBMG_iWORK_CS, NBMG_iWORK_PL,
     &         NBMG_rWORK_CS, NBMG_rWORK_PL

      INTEGER  BMG_iPARMS(NBMG_iPARMS),
     &         BMG_iWORK_CS(NBMG_iWORK_CS),
     &         BMG_iWORK_PL(NBMG_iWORK_PL), 
     &         iWorkMSG(NMSGi), pMSG(NBMG_pMSG,NOGm), 
     &         pMSGSO(NBMG_pMSG,NOGm), IGRD(NOGm,NBMG_pIGRD),
     &         BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS),
     &         BMG_rWORK_CS(NBMG_rWORK_CS),
     &         BMG_rWORK_PL(NBMG_rWORK_PL),
     &         CI(NCI), Q(NF), QF(NF), RES(NF),
     &         SO(NSO), SOR(NSOR), MSG_Buffer(NMSGr)
      REAL*RKIND :: Qmask(NF)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, NProc, NProcI, NProcJ, 
     &         NProcK, MPICOMM, MPI_IERROR, NStncl

      INTEGER  iGs, jGs, kGs, NLx, NLy, NLz, NGx, NGy, NGz

      INTEGER  p_ProcGrid, p_ProcCoord

      INTEGER  i, j, k, p_ABD, p_BBD, p_CI, p_SO, p_SOR, p_U,
     &         NCBW, NCU
      REAL*8   RES_L2
      REAL*8  TT1, TT2

C ==========================================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------
C --------------------------------------------
C     MPI Parallel Data Decomposition:
C --------------------------------------------
      !
      ! MPI Communicators
      !
      MPICOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM)
      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      k = KCF   ! set the current grid index

      IF (KCF.EQ.KC) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,*) '*** FATAL ERROR: BMG3_SymStd_ncycle.f ****'
            WRITE(*,*) 'The current finest grid = ', KCF
            WRITE(*,*) 'The coarest grid        = ', KC
            WRITE(*,*) 'ERROR; minimum number of grids is 2 !'
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,18)
         RETURN

      ENDIF

      ! set the n-cycle counter to zero on all grids
      DO i = KC, KCF
         IGRD(i,idL_BMG_IVW)=iZERO
      END DO


C -----------------------------------------------------------
C     Begin the n-cycle
C     (this is also the n-cycle recursive loop boundary)
C -----------------------------------------------------------

 140  CONTINUE     
         
         IF ( K.EQ.KC .AND. IGRD(k,idL_BMG_IVW).EQ.iZERO ) THEN
            !
            ! Solve on the coarsest grid
            !
            CALL BMG3_SymStd_GET_pointers( 
     &                KC, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, 
     &                NLx, NLy, NLz, NGx, NGy, NGz 
     &                )
            !
            iGs = IGRD(KC,idL_BMG_Icoord)
            jGs = IGRD(KC,idL_BMG_Jcoord)
            kGs = IGRD(KC,idL_BMG_Kcoord)
            !
            IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU ) THEN
               !
               p_ABD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO)
               p_BBD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CU)
               !
               NCBW  = BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW)
               NCU   = BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)
               !
               TT1 = MPI_Wtime()
               !
               CALL BMG3_SymStd_SOLVE_cg_LU( 
     &                   q(p_U), qf(p_U), NLx, NLy, NLz, 
     &                   BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD), 
     &                   NCBW, NCU, NOGm,
     &                   NProcI, NProcJ, NProcK, NProc, MyProc,
     &                   BMG_MSG_iGRID(p_ProcGrid),
     &                   BMG_MSG_iGRID(p_ProcCoord),
     &                   iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                   MSG_Buffer, NMSGr, MPICOMM, BMG_iPARMS
     &                   )
               !
               TT2 = MPI_Wtime()
               !
            ELSE
               !
               TT1 = MPI_Wtime()
               !
               CALL BMG3_SymStd_SOLVE_cg_boxmg( 
     &                   Q(p_U), QF(p_U), NOG,
     &                   BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                   NLx, NLy, NLz, NGx, NGy, NGz,
     &                   iGs, jGs, kGs,
     &                   BMG_iWORK_CS, NBMG_iWORK_CS,
     &                   BMG_rWORK_CS, NBMG_rWORK_CS,
     &                   NProcI, NProcJ, NProcK, NProc, MyProc,
     &                   BMG_MSG_iGRID(p_ProcGrid),
     &                   BMG_MSG_iGRID(p_ProcCoord),
     &                   iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )

               !
               TT2 = MPI_Wtime()
               !
            ENDIF

            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            BMG_rPARMS(id_BMG3_TIME_SOLVE_CG)
     &          = BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) + TT2 - TT1
            
            !
            IGRD(k,idL_BMG_IVW) = IVW
            IF( BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE)
     &          .AND. BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU
     &        ) THEN
               NStncl = 14   ! always the same on coarse grid if NOG > 1
               CALL BMG3_SymStd_residual( 
     &                   KC, NOG, IFD, 
     &                   q(p_U), qf(p_U), so(p_SO), RES(p_U),
     &                   NLx, NLy, NLz, NStncl,
     &                   iWorkMSG, NMSGi, pMSG,
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )               
               CALL BMG3_SymStd_UTILS_norm_l2( 
     &                   RES(p_U), NLx, NLy, NLz, RES_L2, MPICOMM
     &                   )
               !
               IF ( MyProc.EQ.1 ) THEN 
                  WRITE (*,100) k, RES_L2
               ENDIF
               !
            ENDIF
            !
            ! continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_IVW).LT.IVW ) THEN
            !
            ! Move down
            !            
            CALL commas_SymStd_updown( 
     &                k, KCF, NOG, BMG_DOWN, 
     &                IFD, IU, ID, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID 
     &                ,Qmask)
            !
            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            IGRD(k,idL_BMG_IVW) = IGRD(k,idL_BMG_IVW) + 1 
            k=k-1
            !
            ! Continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_IVW).EQ.IVW ) THEN
            !
            ! Move up
            !
            IGRD(k,idL_BMG_IVW) = iZERO   ! reset n-cycle counter to zero
            k=k+1
            !
            CALL commas_SymStd_updown(
     &                k, KCF, NOG, BMG_UP,
     &                IFD, IU, ID, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID 
     &                ,Qmask)
            !
            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            !
            ! Are we on the finest grid? 
            ! 
            IF (K.EQ.KCF) THEN 
               ! n-cycle is done
               RETURN
            ELSE
               ! continue the n-cycle
               GOTO 140  
            ENDIF
            !
         ENDIF

C ==========================================================================

 100     FORMAT (1X,'(3D)LEVEL ',I2,' RESIDUAL NORM = ',1P,E12.5)

C ===========================================

         RETURN
         END

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE commas_SymStd_updown ( 
     &                K, KF, NOG, UPDOWN,
     &                IFD, IU, ID, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &                Q, QF, RES, NF, NC,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID 
     &                ,Qmask)


C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_updown performs all the necessary tasks associated
C     with moving "down" to a coarser grid or moving "up" to a finer
C     grid.  It is necessary because f77 does not support recursion.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C -----------------------------

      INCLUDE 'mpif.h'

#include    "BMG_constants_f77.h" 
      INCLUDE 'BMG_stencils_f77.h'
      INCLUDE 'BMG_parameters_f77.h'
      INCLUDE 'BMG_workspace_f77.h'

C ----------------------------
C     Argument Declarations
C ----------------------------

      !
      !  BoxMG parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      !
      !  BoxMG Cycle Parameters
      !
      INTEGER  k, KC, KCF, KF, IFD, IU, ID, IVW,
     &         IRELAX, IRELAX_SYM, UPDOWN, ISEP

      !
      !  BoxMG dimensions
      !
      INTEGER  NF, NC, NSO, NSOR, NCI, NOGm, NOG,
     &         NBMG_iWORK_CS, NBMG_iWORK_PL,
     &         NBMG_rWORK_CS, NBMG_rWORK_PL

      !
      !  BoxMG internal pointers
      !
      INTEGER  IGRD(NOGm,NBMG_pIGRD)

      !
      !  BoxMG data arrays
      !
      REAL*RKIND   CI(NCI), Q(NF), QF(NF), RES(NF), SO(NSO), SOR(NSOR)
      REAL*RKIND :: Qmask(NF)

      !
      !  BoxMG plane relaxation workspace
      !
      INTEGER BMG_iWORK_PL(NBMG_iWORK_PL)
      REAL*RKIND  BMG_rWORK_PL(NBMG_rWORK_PL)

      !
      !  MSG workspace dimensions
      !
      INTEGER  NMSGi, NMSGr
      !
      !  MSG workspace and pointers
      !
      INTEGER  iWorkMSG(NMSGi), 
     &         pMSG(NBMG_pMSG,NOG),
     &         pMSGSO(NBMG_pMSG,NOG)
      REAL*RKIND   MSG_Buffer(NMSGr)
 
      !
      !  BMG_MSG grid stuff
      !
      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)
      
C ----------------------------
C     Local Declarations
C ----------------------------

      INTEGER  I, J, NLx, NLy, NLz, NLx_c, NLy_c, NLz_c, 
     &         NGx, NGy, NGz, NGx_c, NGy_c, NGz_c, 
     &         p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC
      INTEGER  MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM,
     &         iGs, jGs, kGs, NSORv, NStncl
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, MPI_IERROR,
     &         NProc, NProcI, NProcJ, NProcK, p_ProcGrid, p_ProcCoord
      REAL*8   RES_L2
      REAL*8   TT1, TT2

C ==========================================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------
C --------------------------------------------
C     MPI Parallel Data Decomposition:
C --------------------------------------------
      !
      ! MPI Communicators
      !
      MPICOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM)
      !
      ! MPI Communicatorsn for plane relaxation
      !
      XYPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xy)
      XZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xz)
      YZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_yz)
 
      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      !
      ! Processor Coordinates
      !
      iGs = IGRD(k,idL_BMG_Icoord)
      jGs = IGRD(k,idL_BMG_Jcoord)
      kGs = IGRD(k,idL_BMG_Kcoord)
      
      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      !
      ! Set NStncl for the finest grid
      !
      IF ( K.LT.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         NSORv = 2
      ELSE 
         NSORv = 2  ! the same for now but this may change
      ENDIF

      IF (UPDOWN.EQ.BMG_DOWN) THEN

         ! Get pointers for grids k and k-1
         CALL BMG3_SymStd_GET_pointers(     
     &               K, IGRD, NOGm,
     &               p_U, p_SO, p_SOR, p_CI, 
     &               NLx, NLy, NLz, NGx, NGy, NGz
     &               )
         CALL BMG3_SymStd_GET_pointers( 
     &               K-1, IGRD, NOGm,
     &               p_UC, p_SOC, p_SORC, p_CIC, 
     &               NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c 
     &               )
        
         ! Relaxation
         DO 10 i=1, ID
            !
            TT1 = MPI_Wtime()

            CALL commas_SymStd_relax( 
     &                k, KF, NOG, UPDOWN,
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,
     &                Q, QF, RES, NF, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iGs, jGs, kGs,
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
     &                ,Qmask)

            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_relax) = 
     &           BMG_rPARMS(id_BMG3_TIME_relax) + TT2 - TT1


            !
            IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) .AND. MyProc.EQ.1 ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
            !
 10      CONTINUE

         ! restrict the residual
         TT1 = MPI_Wtime()

         ! calculate the residual
         CALL BMG3_SymStd_residual( 
     &             K, NOG, IFD, 
     &             Q(p_U), QF(p_U), SO(p_SO), RES(p_U),
     &             NLx, NLy, NLz, NStncl,
     &             iWorkMSG, NMSGi, pMSG, 
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) = 
     &        BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) + TT2 - TT1
         
c$$$         CALL BMG3_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, RES(p_U), K, NOG,
c$$$     &             NLx, NLy, NLz, NGx, NGy, NGz,
c$$$     &             IGRD(K,idL_BMG_ICoord),
c$$$     &             IGRD(K,idL_BMG_JCoord),
c$$$     &             IGRD(K,idL_BMG_KCoord),
c$$$     &             'DOWN-RES-post-relax', .FALSE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProcK, NProc,
c$$$     &             MyProcI, MyProcJ, MyProcK, MPICOMM
c$$$     &             )

c         ! restrict the residual
c         TT1 = MPI_Wtime()

         CALL BMG3_SymStd_restrict( 
     &             K, K-1,  
     &             RES(p_U), QF(p_UC), CI(p_CIC),
     &             NLx, NLy, NLz, NLx_c, NLy_c, NLz_c,
     &             iGs, jGs, kGs
     &             )

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_restrict) = 
     &        BMG_rPARMS(id_BMG3_TIME_restrict) + TT2 - TT1
         

c$$$         CALL BMG3_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, QF(p_UC), K-1, NOG,
c$$$     &             NLx_c, NLy_c, NLz_c, NGx_c, NGy, NGz,
c$$$     &             IGRD(K-1,idL_BMG_ICoord),
c$$$     &             IGRD(K-1,idL_BMG_JCoord),
c$$$     &             IGRD(K-1,idL_BMG_KCoord),
c$$$     &             'DOWN-QF-post-restrict', .FALSE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProcK, NProc,
c$$$     &             MyProcI, MyProcJ, MyProcK, MPICOMM
c$$$     &             )

         ! Zero the initial guess
         CALL BMG3_SymStd_UTILS_rV_zero( 
     &               k-1, q(p_UC), NLx_c, NLy_c, NLz_c
     &               )


      ELSE IF (UPDOWN.EQ.BMG_UP) THEN

         ! Get pointers for grids k and k-1
         CALL BMG3_SymStd_GET_pointers(      
     &             k, IGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI, 
     &             NLx, NLy, NLz, NGx, NGy, NGz
     &             )
         CALL BMG3_SymStd_GET_pointers( 
     &             k-1, IGRD, NOGm,
     &             p_UC, p_SOC, p_SORC, p_CIC, 
     &             NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c 
     &             )
         

         ! Interpolate and Correct 
         TT1 = MPI_Wtime()

c        write(0,*) 'SymStd_updown: k, NOG = ',k,NOG
            IF ( k < NOG ) THEN
            CALL BMG3_SymStd_interp_add( 
     &             k-1, k, NOG,
     &             Q(p_U) ,Q(p_UC), RES(p_U),
     &             SO(p_SO), NStncl, CI(p_CIC),
     &             NLx_c, NLy_c, NLz_c,
     &             NLx, NLy, NLz,
     &             iGs, jGs, kGs,
     &             iWorkMSG, NMSGi, pMSG,
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )
            
            ELSE ! special version for fine grid (checks mask array)
c        write(0,*) 'SymStd_updown: call commas_SymStd_interp_add = ',k,NOG
         CALL commas_SymStd_interp_add( 
     &             k-1, k, NOG,
     &             Q(p_U) ,Q(p_UC), RES(p_U),
     &             SO(p_SO), NStncl, CI(p_CIC),
     &             NLx_c, NLy_c, NLz_c,
     &             NLx, NLy, NLz,
     &             iGs, jGs, kGs,
     &             iWorkMSG, NMSGi, pMSG, Qmask,
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )
         ENDIF

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_interp_add) =
     &        BMG_rPARMS(id_BMG3_TIME_interp_add) + TT2 - TT1

                  

         IF ( BMG_IOFLAG(iBMG3_BUG_RES_INTERP) ) THEN
            CALL BMG3_SymStd_residual( 
     &                K, NOG, IFD,
     &                Q(p_U), QF(p_U), SO(p_SO), RES(p_U), 
     &                NLx, NLy, NLz, NStncl,
     &                iWorkMSG, NMSGi, pMSG, 
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            CALL BMG3_SymStd_UTILS_norm_l2( 
     &                RES(p_U), NLx, NLy, NLz, RES_L2, MPICOMM
     &                )
            WRITE(*,110) K, RES_L2
         ENDIF

         ! Relaxation
         DO 20 i=1, IU
            !
            TT1 = MPI_Wtime()
            
            CALL commas_SymStd_relax(  
     &                k, KF, NOG, UPDOWN,
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,
     &                Q, QF, RES, NF, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iGs, jGs, kGs,
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
     &                ,Qmask)
            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_relax) = 
     &           BMG_rPARMS(id_BMG3_TIME_relax) + TT2 - TT1

            !
            IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) .AND. MyProc.EQ.1 ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
            !
 20      CONTINUE

      ELSE
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN         
            WRITE(*,*) '*** UPDOWN out of range: '
            WRITE(*,*) 'HAVE: UPDOWN = ', UPDOWN
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,19)
         RETURN

      ENDIF

C ==========================================================================

 100  FORMAT (1X,'(3D)LEVEL ',I2,' RESIDUAL NORM = ',1P,E12.5)
 110  FORMAT (1X,'(3D)LEVEL ',I2,
     &           ' AFTER INTERPOLATION RES_L2 = ',1P,E12.5 )

C ===========================================

      RETURN
      END

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE commas_SymStd_relax( 
     &                KG, KF, NOG, UPDOWN,
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,
     &                Q, QF, RES, NF, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iGs, jGs, kGs,
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr, 
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
     &                ,Qmask)

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_relax performs one relaxation step (either colored
C     point Gauss Seidel, or alternating plane red-black Gauss-Seidel)
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C ------------------------------------------------
C     Includes
C
#include    "BMG_constants_f77.h"
      INCLUDE 'BMG_workspace_f77.h'
      INCLUDE 'BMG_parameters_f77.h'

C ------------------------------------------------
C     Argument Declarations
C 
      INTEGER  NCI, NSO, NSOR, NSORv, NStncl, kg, kf, NBMG_iWORK_PL,
     &         NBMG_rWORK_PL, NOG, NF, NOGm, NMSGi, NMSGr

      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_iWORK_PL(NBMG_iWORK_PL),
     &         IFD, IGRD(NOGm,24), IRELAX_SYM, UPDOWN, IRELAX, ISEP
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS), BMG_rWORK_PL(NBMG_rWORK_PL),
     &         CI(NCI), q(NF), qf(NF), RES(NF),
     &         SO(NSO), SOR(NSOR), MSG_Buffer(NMSGr)

      REAL*8   RES_L2
      
      REAL*RKIND :: Qmask(NF)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  iGs, jGs, kGs,
     &         MPICOMM, iWorkMSG(NMSGi), pMSG(NBMG_pMSG,NOGm)
      INTEGER  XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, MPI_IERROR,
     &         NProc, NProcI, NProcJ, NProcK


C ----------------------------
C     Local Declarations
C
      INTEGER  NLx, NLy, NLz, NGx, NGy, NGz, 
     &         NLx_f, NLy_f, NLz_f, NGx_f, NGy_f, NGz_f, 
     &         p_CI, p_SO, p_SOR, p_U, p_ProcGrid
      
C ==========================================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      !
      ! MPI Communicators
      !
      MPICOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM)
      !
      ! MPI Communicatorsn for plane relaxation
      !
      XYPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xy)
      XZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xz)
      YZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_yz)
 
      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         NSORv = 2
      ELSE 
         NSORv = 2  ! the same for now but this may change
      ENDIF

      !
      ! Collect fine-grid dimensions 
      ! ( need to check: KF or NOG )
      !
      CALL BMG3_SymStd_GET_pointers( 
     &            KF, IGRD, NOGm,
     &            p_U, p_SO, p_SOR, p_CI, 
     &            NLx_f, NLy_f, NLz_f, NGx_f, NGy_f, NGz_f
     &            ) 

      CALL BMG3_SymStd_GET_pointers(  
     &            KG, IGRD, NOGm,
     &            p_U, p_SO, p_SOR, p_CI, 
     &            NLx, NLy, NLz, NGx, NGy, NGz
     &            )
      

      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         !
         ! Gauss-Seidel relaxation
         !
         CALL commas_SymStd_relax_GS(  
     &             KG, SO(p_SO), QF(p_U), Q(p_U), RES(p_U), SOR(p_SOR),
     &             NLx, NLy, NLz, NGx, NGy, NGz, 
     &             RES_L2, BMG_IOFLAG, NOG, NOGm,
     &             IFD, NStncl, NSORv, IRELAX, IRELAX_SYM, UPDOWN,
     &             iGs, jGs, kGs,
     &             MyProcI, MyProcJ, MyProcK, MyProc,
     &             BMG_MSG_iGRID(p_ProcGrid),
     &             NProcI, NProcJ, NProcK, NProc,
     &             iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr, MPICOMM
     &            ,ISEP,Qmask)
         !
      ELSEIF ( IRELAX.EQ.BMG_GS_RB_planes_xy_yz_xz ) THEN

         WRITE(*,*) 'Fatal Error: To busy to code planes !!!'
         CALL MPI_FINALIZE(MPI_IERROR)
         STOP

         !
      ENDIF

C ==========================================================================

      RETURN
      END


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      SUBROUTINE commas_SymStd_relax_GS( 
     &                KG, SO, QF, Q, RES, SOR,
     &                NLx, NLy, NLz, NGx, NGy, NGz, 
     &                RES_L2, BMG_IOFLAG, NOG, NOGm,
     &                IFD, NStncl, NSORv, IRELAX, IRELAX_SYM, UPDOWN,
     &                iGs, jGs, kGs,
     &                MyProcI, MyProcJ, MyProcK, MyProc,
     &                ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &                iWork, NMSGi, pMSG, MSG_Buffer, NMSGr, MPICOMM
     &               ,ISEP,Qmask)

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_relax_GS performs one sweep of Gauss Seidel (with the
C     correct ordering depending on whether we are on the way down, or
C     up in the symmetric cycling case)
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C ------------------------------------------------
C     Includes
C
      INCLUDE 'mpif.h'
      INCLUDE 'MSG.h'

#include    "BMG_constants_f77.h"
      INCLUDE 'BMG_stencils_f77.h'
      INCLUDE 'BMG_parameters_f77.h'
      INCLUDE 'BMG_workspace_f77.h'

C ------------------------------------------------
C     Argument Declarations
C
      INTEGER NGx, NGy, NGz, NLx, NLy, NLz

      INTEGER IFD, IRELAX, IRELAX_SYM, KG, NOG, NOGm,
     &        NSORv, NStncl, UPDOWN, NMSGi, NMSGr

      INTEGER iGs, jGs, kGs
      INTEGER MyProc, MyProcI, MyProcJ, MyProcK,
     &        NProc, NProcI, NProcJ, NProcK
      INTEGER ProcGrid(NProcI,NProcJ,NProcK)

      INTEGER MPICOMM, iWork(NMSGi), pMSG(NBMG_pMSG,NOGm)
      REAL*RKIND  Q(NLx,NLy,NLz), QF(NLx,NLy,NLz), 
     &        RES(NLx,NLy,NLz),
     &        SO(NLx+1,NLy+1,NLz+1,NStncl), SOR(NLx,NLy,NLz,NSORv)
     
      REAL*8 RES_L2
      
      REAL*RKIND  MSG_Buffer(NMSGr)

      REAL*RKIND :: Qmask(NLx,Nly,Nlz)

      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

C ------------------------------------------------
C     Local Declarations
C
      INTEGER i, i1, j ,j1, k, k1
      
      INTEGER ie, iw, js, jn, kbot, kt

      INTEGER pts, ibeg, iend, ptsred, ptsblk
      INTEGER ptstart, ptend, ptstride
      INTEGER ptrn, ierror
      INTEGER relaxtest
      LOGICAL comhalf
      
      REAL*RKIND SOx(NLx+1),SOy(NLy+1),SOz(NLz+1)
      
      INTEGER ISEP

C ==========================================================================

      j1 = NLy-1
      i1 = NLx-1
      k1 = NLz-1
      
c      write(0,*) 'mproc, mpi,mpj,mpk, igs,jgs,kgs = ',myproc,myproci,
c     :   myprocj,myprock,igs,jgs,kgs
c      write(0,*) 'mproc, isep, BMG_SEPARABLE = ',myproc,isep,
c     :   BMG_SEPARABLE, BMG_NONSEPARABLE
      

c      write(0,*) 'RELAX_GS: kg,nog,ifd = ',KG,NOG,ifd

      comhalf = .false.

      
      IF ( KG.LT.NOG .OR. ifd.NE.1 ) THEN 

         !
         !   27-point relaxation (8-color)
         !

C         IF ((UPDOWN.eq.BMG_UP).or.(IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) 
C     &              THEN
C            ptstart = 1
C            ptend   = 8
C            ptstride = 1
C         ELSE
C            ptstart = 8
C            ptend = 1
C            ptstride = -1
C         ENDIF         
C        
C         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<
C
C            DO k=2+mod(mod((pts-1)/4,2)+mod(kGs+1,2),2),K1,2
C               !
C               DO j=2+mod(mod((pts-1)/2,2)+mod(jGs+1,2),2),J1,2
C                  !
C                  DO i=2+mod(mod(pts-1,2)+mod(iGs+1,2),2),I1,2
                     !

         IF ( (UPDOWN.eq.BMG_UP)
     &        .OR. (IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF         
        
        
        relaxtest = 2
        
        IF ( relaxtest .eq. 0 ) THEN ! default
        
         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

c            DO k=2+mod(mod((pts-1)/4,2)+mod(kGs+1,2),2),K1,2
c               !
c               DO j=2+mod(mod((pts-1)/2,2)+mod(jGs+1,2),2),J1,2
c                  !
c                  DO i=2+mod(mod(pts-1,2)+mod(iGs+1,2),2),I1,2
             DO k = 2,k1
              
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                     !
                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !
            
            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
            !
            
            ENDIF
          ENDDO

        ELSEIF ( relaxtest .eq. 1 ) THEN ! non-redblack GS

             DO k = 2,k1
              
               DO j=2, j1
                  !
                  !
                  DO i = 2,i1
                     !
                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO

        ELSEIF ( relaxtest .eq. 2 ) THEN ! layered redblack GS

c         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

c do 'red' points on first plane
            ptsred = ptstart
            ptsblk = ptend
            DO  k=2, 2
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsred,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO

c do 'red' points on k plane and 'black' points on k-1 plane

            DO  k=3, k1

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsred,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO

            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
            !
            
            ENDIF

c black points
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-2)+ptsblk,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k-1) = ( QF(i,j,k-1) 
     &                    + SO(i,j,k-1,kpw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k-1,kpnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k-1,kps)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k-1,kpsw)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k-1,kpw)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k-1,kpnw)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k-1,kps)*Q(i,j-1,k-1)
     &                    + SO(i,j,k-1,kpsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k-1,kb)*Q(i,j,k-1-1)
     &                    + SO(i,j,k-1,kbw)*Q(i-1,j,k-1-1)
     &                    + SO(i,j+1,k-1,kbnw)*Q(i-1,j+1,k-1-1)
     &                    + SO(i,j+1,k-1,kbn)*Q(i,j+1,k-1-1)
     &                    + SO(i+1,j+1,k-1,kbne)*Q(i+1,j+1,k-1-1)
     &                    + SO(i+1,j,k-1,kbe)*Q(i+1,j,k-1-1)
     &                    + SO(i+1,j,k-1,kbse)*Q(i+1,j-1,k-1-1)
     &                    + SO(i,j,k-1,kbs)*Q(i,j-1,k-1-1)
     &                    + SO(i,j,k-1,kbsw)*Q(i-1,j-1,k-1-1)
     &                    + SO(i,j,k-1+1,kb)*Q(i,j,k-1+1)
     &                    + SO(i,j,k-1+1,kbe)*Q(i-1,j,k-1+1)
     &                    + SO(i,j+1,k-1+1,kbse)*Q(i-1,j+1,k-1+1)
     &                    + SO(i,j+1,k-1+1,kbs)*Q(i,j+1,k-1+1)
     &                    + SO(i+1,j+1,k-1+1,kbsw)*Q(i+1,j+1,k-1+1)
     &                    + SO(i+1,j,k-1+1,kbw)*Q(i+1,j,k-1+1)
     &                    + SO(i+1,j,k-1+1,kbnw)*Q(i+1,j-1,k-1+1)
     &                    + SO(i,j,k-1+1,kbn)*Q(i,j-1,k-1+1)
     &                    + SO(i,j,k-1+1,kbne)*Q(i-1,j-1,k-1+1)
     &                    )*SOR(i,j,k-1,msor)
                     !
                  ENDDO
                  !
               ENDDO

               !
            ENDDO

c do 'black' points on last plane
            DO  k=k1, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsblk,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            
       ENDIF ! relaxtest
       !
c         ENDDO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<

        IF ( relaxtest .ne. 0 ) THEN
 
            ptrn = 1
            !
            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
     
        ENDIF
            !
!         ENDDO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<


         !
      ELSE
         !  
         !  7 point relaxation (2-color)
         !         
         IF ( (UPDOWN.eq.BMG_UP)
     &        .OR. (IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF         

      IF ( ISEP .eq. BMG_SEPARABLE ) THEN

      
      DO  k=1, k1+1
        SOz(k) = SO(3,3,k,kb)
        SOx(k) = SO(3,3,k,kpw)
        SOy(k) = SO(3,3,k,kps)
c       IF ( myproc == 1 ) write(0,*) 'k, SOz,x,y = ',k,SOz(k),
c     &  SO(3,3,k,kpw),SO(3,3,k,kps)
c        IF ( myproc == 1 ) write(0,*) 'j, SOy = ',
      ENDDO
c        SOz(2) = SOz(3)
c        SOz(k1+1) = SOz(k1)
        
      
      DO j=1, j1+1
c        SOy(j) = SO(3,j,3,kps)
c        IF ( myproc == 1 ) write(0,*) 'j, SOy,x,z = ',j,SOy(j),
c     &   SO(3,j,3,kpw),SO(3,j,3,kb)
c       IF ( myproc == 1 ) write(0,*) 'k, SOz = ',j,SO(3,j,3,kb)
c        IF ( myproc == 1 ) write(0,*) 'i, SOx = ',j,SO(3,j,3,kpw)
c        write(0,*) 'j, SOy, myproc = ',j,SOy(j),myproc
      ENDDO
c        SOy(2) = SOy(3)
c        SOy(j1+1) = SOy(j1)
      
      
      DO i=1,i1+1
c        SOx(i) = SO(i,3,3,kpw)
c        IF ( myproc == 1 ) write(0,*) 'i, SOx,y,z = ',i,SOx(i),
c     &   SO(i,3,3,kps),SO(i,3,3,kb)
      ENDDO
c        SOx(2) = SOx(3)
c        SOx(i1+1) = SOx(i1)
      
      ENDIF

         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

          IF ( ISEP .eq. BMG_SEPARABLE ) THEN

          IF ( .true. ) THEN
          
          
            DO  k=2, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SOx(k)*Q(i-1,j,k)
     &                        + SOy(k)*Q(i,j+1,k)
     &                        + SOx(k)*Q(i+1,j,k)
     &                        + SOy(k)*Q(i,j-1,k)
     &                        + SOz(k)*Q(i,j,k-1)
     &                        + SOz(k+1)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)
                  END DO
                  !
               END DO
               !
            END DO
           
          ELSE ! old try

          kbot = 0
          kt = 0
          js = 0
          jn = 0
          iw = 0
          ie = 0
          
          IF (  MyProcK .eq. 1 )     kbot = 1
          IF (  MyProcK .eq. NProcK) kt = 1
          IF (  MyProcJ .eq. 1)      js = 1
          IF (  MyProcJ .eq. NProcJ) jn = 1
          IF (  MyProcI .eq. 1)      iw = 1
          IF (  MyProcI .eq. NProcI) ie = 1
          
            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  
                  IF ( ibeg .eq. 2 ) ibeg = ibeg + 2*iw
                  IF ( iend .eq. i1 ) iend = iend - 2*ie
                  
                  !
                  DO i=ibeg,iend,2
                  !
c                  IF ( .not. (( (k.eq.2      .and. MyProcK.eq.1 ) .or.
c     &                        ( k .eq. k1    .and. MyProcK .eq. NProcK) .or.
c     &                        ( j .eq. 2     .and. MyProcJ .eq. 1)  .or.
c     &                        ( j .eq. j1    .and. MyProcJ .eq. NProcJ) .or.
c     &                        ( i .eq. ibeg  .and. MyProcI .eq. 1)  .or. 
c     &                        ( i .eq. iend  .and. MyProcI .eq. NProcI) ) ) THEN
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SOx(i)*Q(i-1,j,k)
     &                        + SOy(j+1)*Q(i,j+1,k)
     &                        + SOx(i+1)*Q(i+1,j,k)
     &                        + SOy(j)*Q(i,j-1,k)
     &                        + SOz(k)*Q(i,j,k-1)
     &                        + SOz(k+1)*Q(i,j,k+1)
     &                        )/(SOx(i) + SOx(i+1) + SOy(j) + 
     &                           SOy(j+1) + SOz(k) + SOz(k+1))

                  END DO
                  !
               END DO
               !
            END DO

c edges:

           IF ( kbot .eq. 1 ) THEN
            k=2

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg 
                  !
                  DO i=ibeg,iend,2

                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)


                 END DO
              END DO
           
           ENDIF

           IF ( kt .eq. 1 ) THEN
            k=k1

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg 
                  !
                  DO i=ibeg,iend,2

                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)


                 END DO
              END DO
           
           ENDIF
           
           IF ( js .eq. 1 ) THEN
                      
            DO  k=2+kbot, k1-kt
              j = 2
              ibeg = mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend = 2*((i1-ibeg)/2)+ibeg 

               DO i=ibeg,iend,2
    
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)

                END DO
             END DO
          ENDIF

           
           IF ( jn .eq. 1 ) THEN
                      
            DO  k=2+kbot, k1-kt
              j = j1
              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg 
                  DO i=ibeg,iend,2

    
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)

                END DO
             END DO
          ENDIF

           
           IF ( ie .eq. 1 ) THEN
            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn

              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg 
              i = ibeg
               IF ( ibeg .eq. 2 ) THEN
    
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)

                ENDIF
                END DO
             END DO
          ENDIF

           IF ( iw .eq. 1 ) THEN
            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn

              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg 
              
              i = iend
              
               IF ( i .eq. i1 ) THEN
    
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)

                ENDIF
                END DO
             END DO
          ENDIF
          
         ENDIF
         

         ELSE  ! Non separable

            DO  k=2, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !
                  
                  IF ( Qmask(i,j,k) .gt. 0.9e30 ) THEN ! no channel point here
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)
                  ELSE ! there is a channel point here, so reset q(i,j,k) to channel potential (value in mask)
                       ! and solve for qf(i,j,k) instead
                      Q(i,j,k) = Qmask(i,j,k)
!                      QF(i,j,k) =  Qmask(i,j,k,1)*SO(i,j,k,kp) - 
                      QF(i,j,k) =  Qmask(i,j,k)/SOR(i,j,k,msor) - 
     &                        ( SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )
                      
                  ENDIF
                  END DO
                  !
               END DO
               !
            END DO

          ENDIF
            !

            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
            !
            
            ENDIF
            
            
            !
         END DO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<
            !
            IF ( .not. comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
            !
            
            ENDIF
            !
C         END DO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<

      ENDIF

      ptrn = 1

      CALL MSG_tbdx_close(Q, MSG_Buffer,
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &     iWork(pMSG(ipL_MSG_Proc,KG)),
     &     iWork(pMSG(ipL_MSG_Ipr,KG)),
     &     iWork(pMSG(ipL_MSG_Index,KG)),
     &     ptrn, ierror)

      IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) ) THEN
         !
         CALL BMG3_SymStd_residual( 
     &             KG, NOG, IFD,
     &             Q, QF, SO, RES, NLx, NLy, NLz, NStncl, 
     &             iWork, NMSGi, pMSG, 
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )
         CALL BMG3_SymStd_UTILS_norm_l2( 
     &             RES, NLx, NLy, NLz, RES_L2, MPICOMM
     &             )
         !
      ENDIF

C ==========================================================================

      RETURN
      END

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE commas_SymStd_interp_add(
     &                KCG, KFG, NOG,
     &                Q ,QC, RES, SO, NStncl, CI,
     &                IIC, JJC, KKC,
     &                IIF, JJF, KKF,
     &                iGs, jGs, kGs, 
     &                iWork, NMSGi, pMSG, Qmask,
     &                BUFFER, NMSGr, MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   commas_SymStd_interp_add.f interpolates Q from the coarse mesh, KCG, to 
C   the fine mesh, KFG, and adds the result to Q on fine mesh.
c
c  Special version for finest mesh: checks mask array and does not update Q 
c  in block around lightning flash 
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'mpif.h'
      INCLUDE 'MSG.h'

      INCLUDE 'BMG_stencils_f77.h'
      INCLUDE 'BMG_workspace_f77.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER IIC, IIF, JJC, JJF, KCG, KFG, KKC, KKF
      INTEGER iGs, jGs, kGs, NMSGi, NMSGr, NOG, NStncl

      INTEGER iWork(NMSGi), pMSG(NBMG_pMSG,NOG), MPICOMM
      REAL*RKIND  CI(IIC,JJC,KKC,26), Q(IIF,JJF,KKF), QC(IIC,JJC,KKC) 
      REAL*RKIND  RES(IIF,JJF,KKF), SO(IIF+1,JJF+1,KKF+1,NStncl)
      REAL*RKIND  BUFFER(NMSGr)

      REAL*RKIND :: Qmask(IIF,JJF,KKF)

C ----------------------------
C     Local Declarations
C
      INTEGER IC, I, IICF, IICF1, IIC1,
     &        JC, J, JJCF, JJCF1, JJC1,
     &        KC, K, KKCF, KKCF1, KKC1
      INTEGER ISTART, JSTART, KSTART, ICSTART, JCSTART, KCSTART
      INTEGER ierror, ptrn
      REAL*RKIND  A, AQ
      
      INTEGER imn,imx, jmn,jmx, kmn,kmx
      LOGICAL flash
      REAL    maskmn, maskmx

C ==========================================================================

C -------------------------------------------------
C     Useful index bounds:
C -------------------------------------------------

      iic1=iic-1
      jjc1=jjc-1
      kkc1=kkc-1

      iicf=(iif-2)/2+3
      jjcf=(jjf-2)/2+3
      kkcf=(kkf-2)/2+3

      IF ((mod(IIF,2).eq.0).and.(mod(iGs,2).eq.0)) THEN
         IICF = IICF - 1
      ENDIF

      IF ((mod(JJF,2).eq.0).and.(mod(jGs,2).eq.0)) THEN
         JJCF = JJCF - 1
      ENDIF
      
      IF ((mod(KKF,2).eq.0).and.(mod(kGs,2).eq.0)) THEN
         KKCF = KKCF - 1
      ENDIF

      IICF1=IICF-1
      JJCF1=JJCF-1
      KKCF1=KKCF-1

C --------------------------------------------------
C   interpolate answers from coarse to fine mesh 
C   and add to answers on fine mesh.
C --------------------------------------------------

c     odd values of iGs, jGs, or kGs indicate the
c     the case similar to the serial version


      IF (mod(iGs,2).eq.1) THEN
         ISTART = 2
         ICSTART = 2
      ELSE
         ISTART = 1
         ICSTART = 1
      ENDIF

      IF (mod(jGs,2).eq.1) THEN
         JSTART = 2
         JCSTART = 2
      ELSE
         JSTART = 1
         JCSTART = 1
      ENDIF
      
      IF (mod(kGs,2).eq.1) THEN
         KSTART = 2
         KCSTART = 2
      ELSE
         KSTART = 1
         KCSTART = 1
      ENDIF

C --------------------------------------------------
C     NB: division is only possible in the interior
C --------------------------------------------------

      imn = iif
      imx = 0
      jmn = jjf
      jmx = 0
      kmn = kkf
      kmx = 0
      flash = .false.
      
      maskmn = 1.e31
      maskmx = 0
      
      do k=2, kkf-1
         do j=2, jjf-1
            do i=2, iif-1
               RES(i,j,k)=RES(i,j,k)/so(i,j,k,kp)
               
               maskmn = Min( maskmn, Qmask(i,j,k) )
               maskmx = Max( maskmx, Qmask(i,j,k) )
               
               IF ( Qmask(i,j,k) .lt. 0.9e30 ) THEN ! mask value of 1.e30 indicates no channel
                 flash = .true.
                 imn = Min( imn, i-1 )
                 jmn = Min( jmn, j-1 )
                 kmn = Min( kmn, k-1 )

                 imx = Max( imx, i+1 )
                 jmx = Max( jmx, j+1 )
                 kmx = Max( kmx, k+1)
c           write(0,*) 'commas_SymStd_interp_add: ijk,qmask= ',i,j,k,Qmask(i,j,k)
c            write(0,*) 'commas_SymStd_interp_add: imn: ',imn,imx
               ENDIF
            enddo
         enddo
      enddo
      
      IF ( imn .le. imx ) THEN
c      write(0,*) 'flash detected in commas_SymStd_interp_add'
      flash = .true.
      imn = Max(1, imn-3)
      jmn = Max(1, jmn-3)
      kmn = Max(1, kmn-3)

      imx = Max(1, imx+3)
      jmx = Max(1, jmx+3)
      kmx = Max(1, kmx+3)
      
      ELSE
c      write(0,*) 'no flash detected in commas_SymStd_interp_add, mx,mn=',maskmn,maskmx
c      write(0,*) 'commas_SymStd_interp_add: iif,jjf,kkf = ',iif,jjf,kkf
      ENDIF
      
      
      IF ( .not. flash ) THEN
      
      k=KSTART-2
      j=JSTART
      i=ISTART

      do kc=KCSTART,kkc1

         i=ISTART
         j=JSTART
         k=k+2

         ic=ICSTART
         jc=JCSTART

         q(i,j,k)=q(i,j,k)+qc(ic,jc,kc)

         do ic=ICSTART+1,iicf1
            i=i+2
            q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)
            a = ci(ic,jc,kc,lxyr)*qc(ic,jc,kc)
     &           + ci(ic,jc,kc,lxyl)*qc(ic-1,jc,kc)
            q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)
         enddo
         
         do jc=JCSTART+1,jjcf1

            i=ISTART
            j=j+2

            ic=ICSTART

            q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)
            aq = ci(2,jc,kc,lxya)*qc(2,jc,kc)
     &         + ci(2,jc,kc,lxyb)*qc(2,jc-1,kc)
            q(i,j-1,k) = q(i,j-1,k) + aq + RES(i,j-1,k)
            
            do ic=ICSTART+1,iicf1
               
               i=i+2

               q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)
               a = ci(ic,jc,kc,lxyr)*qc(ic,jc,kc)
     &              + ci(ic,jc,kc,lxyl)*qc(ic-1,jc,kc)
               q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)
               
               aq = ci(ic,jc,kc,lxya)*qc(ic,jc,kc)
     &              + ci(ic,jc,kc,lxyb)*qc(ic,jc-1,kc)
               q(i,j-1,k) = q(i,j-1,k) + aq + RES(i,j-1,k)
               
               a = ci(ic,jc,kc,lxysw)*qc(ic-1,jc-1,kc)
     &              + ci(ic,jc,kc,lxynw)*qc(ic-1,jc,kc)
     &              + ci(ic,jc,kc,lxyne)*qc(ic,jc,kc)
     &              + ci(ic,jc,kc,lxyse)*qc(ic,jc-1,kc)
               q(i-1,j-1,k) = q(i-1,j-1,k) + a + RES(i-1,j-1,k)
               
            enddo
         enddo
      enddo


      k=KSTART-1
      do kc=KCSTART+1,kkcf1

         i=ISTART-2
         j=JSTART
         k=k+2

         jc=JCSTART

         do ic=ICSTART,iic1

            i=i+2

            q(i,j,k) = q(i,j,k) + ci(ic,jc,kc,lxza)*qc(ic,jc,kc)
     &           + ci(ic,jc,kc,lxzb)*qc(ic,jc,kc-1) + RES(i,j,k)

         enddo

         j=JSTART
         do jc=JCSTART+1,jjcf1

            i=ISTART-2
            j=j+2
            
            do ic=ICSTART,iic1

               i=i+2

               q(i,j,k) = q(i,j,k) 
     &                  + ci(ic,jc,kc,lxza)*qc(ic,jc,kc)
     &                  + ci(ic,jc,kc,lxzb)*qc(ic,jc,kc-1)
     &                  + RES(i,j,k)

               q(i,j-1,k) = q(i,j-1,k)
     &                    + ci(ic,jc,kc,lyznw)*qc(ic,jc,kc)
     &                    + ci(ic,jc,kc,lyzne)*qc(ic,jc-1,kc)
     &                    + ci(ic,jc,kc,lyzsw)*qc(ic,jc,kc-1)
     &                    + ci(ic,jc,kc,lyzse)*qc(ic,jc-1,kc-1)
     &                    + RES(i,j-1,k)

            enddo
         enddo
      enddo


      k=KSTART-1
      do kc=KCSTART+1,kkcf1
         k=k+2
         
         j=JSTART
         jc=JCSTART
         i=ISTART-1

         do ic=ICSTART+1,iicf1
            i=i+2

            q(i,j,k) = q(i,j,k)
     &               + ci(ic,jc,kc,lxznw)*qc(ic-1,jc,kc)
     &               + ci(ic,jc,kc,lxzne)*qc(ic,jc,kc)
     &               + ci(ic,jc,kc,lxzsw)*qc(ic-1,jc,kc-1)
     &               + ci(ic,jc,kc,lxzse)*qc(ic,jc,kc-1)
     &               + RES(i,j,k)

         enddo

         j=JSTART
         do jc=JCSTART+1,jjcf1
            j=j+2

            i=ISTART-1
            do ic=ICSTART+1,iicf1
               i=i+2

               q(i,j,k) = q(i,j,k) 
     &                  + ci(ic,jc,kc,lxznw)*qc(ic-1,jc,kc)
     &                  + ci(ic,jc,kc,lxzne)*qc(ic,jc,kc)
     &                  + ci(ic,jc,kc,lxzsw)*qc(ic-1,jc,kc-1)
     &                  + ci(ic,jc,kc,lxzse)*qc(ic,jc,kc-1)
     &                  + RES(i,j,k)

               q(i,j-1,k) = q(i,j-1,k) 
     &                    + ci(ic,jc,kc,ltnw)*qc(ic-1,jc,kc)
     &                    + ci(ic,jc,kc,ltne)*qc(ic,jc,kc)
     &                    + ci(ic,jc,kc,ltsw)*qc(ic-1,jc-1,kc)
     &                    + ci(ic,jc,kc,ltse)*qc(ic,jc-1,kc)
     &                    + ci(ic,jc,kc,lbnw)*qc(ic-1,jc,kc-1)
     &                    + ci(ic,jc,kc,lbne)*qc(ic,jc,kc-1)
     &                    + ci(ic,jc,kc,lbsw)*qc(ic-1,jc-1,kc-1)
     &                    + ci(ic,jc,kc,lbse)*qc(ic,jc-1,kc-1)
     &                    + RES(i,j-1,k)
            enddo

         enddo

      enddo
      
      ELSE ! there are channel points in the tile
      k=KSTART-2
      j=JSTART
      i=ISTART

      do kc=KCSTART,kkc1

         i=ISTART
         j=JSTART
         k=k+2

         ic=ICSTART
         jc=JCSTART

         IF ( .not. ( i .ge. imn .and. i .le. imx .and.
     :                j .ge. jmn .and. j .le. jmx .and.
     :                k .ge. kmn .and. k .le. kmx ) ) THEN
         q(i,j,k)=q(i,j,k)+qc(ic,jc,kc)
         ENDIF

         do ic=ICSTART+1,iicf1
            i=i+2
         IF ( .not. ( i .ge. imn .and. i .le. imx .and.
     :                j .ge. jmn .and. j .le. jmx .and.
     :                k .ge. kmn .and. k .le. kmx ) ) THEN
            q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)
            a = ci(ic,jc,kc,lxyr)*qc(ic,jc,kc)
     &           + ci(ic,jc,kc,lxyl)*qc(ic-1,jc,kc)
            q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)
         ENDIF
         enddo
         
         do jc=JCSTART+1,jjcf1

            i=ISTART
            j=j+2

            ic=ICSTART

            IF ( .not. ( i .ge. imn .and. i .le. imx .and.
     :                   j .ge. jmn .and. j .le. jmx .and.
     :                   k .ge. kmn .and. k .le. kmx ) ) THEN
            q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)
            aq = ci(2,jc,kc,lxya)*qc(2,jc,kc)
     &         + ci(2,jc,kc,lxyb)*qc(2,jc-1,kc)
            q(i,j-1,k) = q(i,j-1,k) + aq + RES(i,j-1,k)
            ENDIF
            
            do ic=ICSTART+1,iicf1
               
               i=i+2

             IF ( .not. ( i .ge. imn .and. i .le. imx .and.
     :                    j .ge. jmn .and. j .le. jmx .and.
     :                    k .ge. kmn .and. k .le. kmx ) ) THEN
               q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)
               a = ci(ic,jc,kc,lxyr)*qc(ic,jc,kc)
     &              + ci(ic,jc,kc,lxyl)*qc(ic-1,jc,kc)
               q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)
               
               aq = ci(ic,jc,kc,lxya)*qc(ic,jc,kc)
     &              + ci(ic,jc,kc,lxyb)*qc(ic,jc-1,kc)
               q(i,j-1,k) = q(i,j-1,k) + aq + RES(i,j-1,k)
               
               a = ci(ic,jc,kc,lxysw)*qc(ic-1,jc-1,kc)
     &              + ci(ic,jc,kc,lxynw)*qc(ic-1,jc,kc)
     &              + ci(ic,jc,kc,lxyne)*qc(ic,jc,kc)
     &              + ci(ic,jc,kc,lxyse)*qc(ic,jc-1,kc)
               q(i-1,j-1,k) = q(i-1,j-1,k) + a + RES(i-1,j-1,k)
               ENDIF
               
            enddo
         enddo
      enddo


      k=KSTART-1
      do kc=KCSTART+1,kkcf1

         i=ISTART-2
         j=JSTART
         k=k+2

         jc=JCSTART

         do ic=ICSTART,iic1

            i=i+2

           IF ( .not. ( i .ge. imn .and. i .le. imx .and.
     :                  j .ge. jmn .and. j .le. jmx .and.
     :                  k .ge. kmn .and. k .le. kmx ) ) THEN
            q(i,j,k) = q(i,j,k) + ci(ic,jc,kc,lxza)*qc(ic,jc,kc)
     &           + ci(ic,jc,kc,lxzb)*qc(ic,jc,kc-1) + RES(i,j,k)

            ENDIF
         enddo

         j=JSTART
         do jc=JCSTART+1,jjcf1

            i=ISTART-2
            j=j+2
            
            do ic=ICSTART,iic1

               i=i+2

            IF ( .not. ( i .ge. imn .and. i .le. imx .and.
     :                   j .ge. jmn .and. j .le. jmx .and.
     :                   k .ge. kmn .and. k .le. kmx ) ) THEN
               q(i,j,k) = q(i,j,k) 
     &                  + ci(ic,jc,kc,lxza)*qc(ic,jc,kc)
     &                  + ci(ic,jc,kc,lxzb)*qc(ic,jc,kc-1)
     &                  + RES(i,j,k)

               q(i,j-1,k) = q(i,j-1,k)
     &                    + ci(ic,jc,kc,lyznw)*qc(ic,jc,kc)
     &                    + ci(ic,jc,kc,lyzne)*qc(ic,jc-1,kc)
     &                    + ci(ic,jc,kc,lyzsw)*qc(ic,jc,kc-1)
     &                    + ci(ic,jc,kc,lyzse)*qc(ic,jc-1,kc-1)
     &                    + RES(i,j-1,k)

              ENDIF
            enddo
         enddo
      enddo


      k=KSTART-1
      do kc=KCSTART+1,kkcf1
         k=k+2
         
         j=JSTART
         jc=JCSTART
         i=ISTART-1

         do ic=ICSTART+1,iicf1
            i=i+2

         IF ( .not. ( i .ge. imn .and. i .le. imx .and.
     :                j .ge. jmn .and. j .le. jmx .and.
     :                k .ge. kmn .and. k .le. kmx ) ) THEN
            q(i,j,k) = q(i,j,k)
     &               + ci(ic,jc,kc,lxznw)*qc(ic-1,jc,kc)
     &               + ci(ic,jc,kc,lxzne)*qc(ic,jc,kc)
     &               + ci(ic,jc,kc,lxzsw)*qc(ic-1,jc,kc-1)
     &               + ci(ic,jc,kc,lxzse)*qc(ic,jc,kc-1)
     &               + RES(i,j,k)

          ENDIF
         enddo

         j=JSTART
         do jc=JCSTART+1,jjcf1
            j=j+2

            i=ISTART-1
            do ic=ICSTART+1,iicf1
               i=i+2

         IF ( .not. ( i .ge. imn .and. i .le. imx .and.
     :                j .ge. jmn .and. j .le. jmx .and.
     :                k .ge. kmn .and. k .le. kmx ) ) THEN
               q(i,j,k) = q(i,j,k) 
     &                  + ci(ic,jc,kc,lxznw)*qc(ic-1,jc,kc)
     &                  + ci(ic,jc,kc,lxzne)*qc(ic,jc,kc)
     &                  + ci(ic,jc,kc,lxzsw)*qc(ic-1,jc,kc-1)
     &                  + ci(ic,jc,kc,lxzse)*qc(ic,jc,kc-1)
     &                  + RES(i,j,k)

               q(i,j-1,k) = q(i,j-1,k) 
     &                    + ci(ic,jc,kc,ltnw)*qc(ic-1,jc,kc)
     &                    + ci(ic,jc,kc,ltne)*qc(ic,jc,kc)
     &                    + ci(ic,jc,kc,ltsw)*qc(ic-1,jc-1,kc)
     &                    + ci(ic,jc,kc,ltse)*qc(ic,jc-1,kc)
     &                    + ci(ic,jc,kc,lbnw)*qc(ic-1,jc,kc-1)
     &                    + ci(ic,jc,kc,lbne)*qc(ic,jc,kc-1)
     &                    + ci(ic,jc,kc,lbsw)*qc(ic-1,jc-1,kc-1)
     &                    + ci(ic,jc,kc,lbse)*qc(ic,jc-1,kc-1)
     &                    + RES(i,j-1,k)
            ENDIF
            
            enddo

         enddo

      enddo
      
      ENDIF





      ptrn = 1
      
      call MSG_tbdx_send(Q, buffer, 
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KFG)),
     &     iWork(pMSG(ipL_MSG_Proc,KFG)),
     &     iWork(pMSG(ipL_MSG_Ipr,KFG)),
     &     iWork(pMSG(ipL_MSG_Index,KFG)),
     &     ptrn, ierror)
      
      call MSG_tbdx_receive(Q, buffer,
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KFG)),
     &     iWork(pMSG(ipL_MSG_Proc,KFG)),
     &     iWork(pMSG(ipL_MSG_Ipr,KFG)),
     &     iWork(pMSG(ipL_MSG_Index,KFG)),
     &     ptrn, ierror)
      
      call MSG_tbdx_close(Q, buffer,
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KFG)),
     &     iWork(pMSG(ipL_MSG_Proc,KFG)),
     &     iWork(pMSG(ipL_MSG_Ipr,KFG)),
     &     iWork(pMSG(ipL_MSG_Index,KFG)),
     &     ptrn, ierror)
      
C ==========================================================================

      RETURN
      END



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE commas_SymStd_SOLVE_boxmgdg(
     &                NLx, NLy, NLz, NGx, NGy, NGz, iGs, jGs, kGs, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm, 
     &                SO, NSOm, 
     &                SOR, NSORm, 
     &                CI, NCIm,
     &                IGRD, NOGm, NOG,
     &                BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                BMG_rWORK_PL, NBMG_rWORK_PLm,
     &                BMG_iWORK_CS, NBMG_iWORK_CSm,
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                iWorkMSG, NMSGim, 
     &                pMSG, pMSGSO, 
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                NProc, MSG_Buffer, NMSGrm, 
     &                MPICOMM)
c     &                ,hx,hy,hz)

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_SOLVE_boxmg is a black box matrix solver. It takes a
C     matrix defined by the user from a given (fine) grid and constructs
C     coarser grids and their associated coefficient matrices. BMG3D
C     then performs the multigrid algorithm and returns the solution
C     vec- tor. The matrix is set by a user written subroutine. The
C     difference stencil at (i,j,k) spans three planes. The k-plane has
C     difference stencil
C
C                        xnw xn xne
C                        xw  xo xe
C                        xsw xs xse.
C
C     The (k+1)-plane has difference stencil
C
C                        xtnw xtn xtne
C                        xtw  xt  xte
C                        xtsw xts xtse.
C
C     The(k-1)-plane has difference stencil
C
C                        xbnw xbn xbne
C                        xbw  xb  xbe
C                        xbsw xbs xbse.
C
C     Here 
C           xo   = SO(i,j,k,ko)
C           xw   = - SO(i,j,k,kw)
C           xs   = - SO(i,j,k,ks)
C           xb   = - SO(i,j,k,kb)
C           xsw  = - SO(i,j,k,ksw)
C           xnw  = - SO(i,j+1,k,knw)
C           xbw  = - SO(i,j,k,kbw)
C           xbnw = - SO(i,j+1,k,kbnw)
C           xbn  = - SO(i,j+1,k,kbn)
C           xbne = - SO(i+1,j+1,k,kbne)
C           xbe  = - SO(i+1,j,k,kbe)
C           xbse = - SO(i+1,j,k,kbse)
C           xbs  = - SO(i,j,k,kbs)
C           xbsw = - SO(i,j,k,kbsw)
C
C     The remaining coefficents:
C
C     xn, xne,xe,xse, xt, xtw, xtnw, xtn, xtne, xte, xtse, xts, and xtsw 
C
C     are specified by symmetry. In the case of a seven point operator,
C     only xo, xw, xs, and xb need to be set. The difference scheme is
C     assumed to be positive definite.
C     
C     The fictitious points 
C
C         (1,j,k), j=1,jj, k=1,kk, (ii,j,k), j=1,jj, k=1,kk, 
C         (i,1,k), i=1,ii, k=1,kk, (i,jj,k), i=1,ii, k=1,kk, 
C         (i,j,1), i=1,ii, j=1,jj, (i,j,kk), i=1,ii, j=1,jj,
C
C     are assumed for ease of programming. However, it is also assumed
C     that the user has set any difference coefficient referring
C     to these points to zero. For example, if the above difference
C     template is centered at (2,2), then xb, xbw, xbnw, xbn, xbne, xbe,
C     xbse, xbs xnw, xw, xsw, xs, xse, xtnw, xtw, xtsw, xts, xtse
C     would be assumed zero.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     NXM      x-dimension of the grid, excluding fictitious points
C     NYM      y-dimension of the grid, excluding fictitious points
C     NZM      z-dimension of the grid, excluding fictitious points
C
C     TOL      Convergence tolerance.
C
C     ISTOP    Indicates the stopping criteria. At present there
C              is a relative residual test
C
C                ISTOP = BMG_STOP_REL_RES_L2 
C
C                => iteration has converged if the l2-norm of the
C                   current residual, divided by the l2-norm of the 
C                   initial residual is less than TOL.
C
C              and an absolute residual test
C  
C                ISTOP = BMG_STOP_REL_ABS_L2
C
C                => iteration has converged if the l2-norm of the 
C                  current residual is less than TOL.
C
C              NOTE: if the selected test fails and the number of
C              iterations is greater than ABS(ISTRT) then the iteration
C              has failed => set TOL=-RES_L2 and return.
C
C
C     IFD      Indicator for difference scheme. IFD=1 means a seven
C              point operator. I.e., in the template given above, xnw
C              xne, xse, xsw, xbw, xbnw, xbn, xbne, xbe, xbse, xbs, xbsw,
C              xtw, xtnw, xtn, xtne, xte, xtse, xts, xtsw.
C              are assumed zero. IFD.ne.1. means a 27 point
C              operator. (Note that 27 point operators will be
C              generated on the coarser grids even in the case of a seven
C              point operator.
C
C     IU       Number of relaxation sweeps to be performed
C              on a coarse grid before interpolation to a fine grid.
C
C     ID       Number of relaxation sweeps to be performed on a fine
C              grid before the problem is transferred to a coarse grid.
C
C     IM       Unused
C
C     ISTRT    Indicator for number of multigrid cycles to be performed
C              if ISTRT.gt.0 the algorithm will begin on the coarsest grid
C              with a full multigrid cycle (however,
C              without cubic interpolation). It will continue
C              cycling until ISTRT cycles have been performed or until
C              the error criterion is satisfied. If ISTRT.lt.0 the
C              algorithm will begin on the finest grid and will perform
C              -ISTRT cycles unless the error criterion is satisfied
C              first.
C
C     ISETUP   Indicator for whether or not to skip initial setup,
C              initial guess for Q, computation of pointers, and
C              generation of coefficients on the coarser grids.
C            
C                ISETUP=0  => perform setup, including,
C                             - compute workspace pointers
C                             - compute coarse-grid operators
C                             - compute interpolation transfer operators
C
C              If ISETUP.NE.0 then we have the following cases:
C         
C                ISETUP=1  => skip computation of workspace pointers
C                             - compute coarse-grid operators
C                             - compute interpolation operators      
C                ISETUP=2  => skip all setup 
C
C              It is important to note the interaction that ISETUP 
C              has with ISTRT.  Specifically, with ISETUP.NE.0 we have
C            
C                 ISTRT > 0 => FMG so residual is coarsened appropriately
C                 ISTRT < 0 => n-cycle so nothing extra to do here.
C
C              An example of usage would be a constant coefficient 
C              time-dependent problem where ISETUP=0 on the first time step 
C              and ISETUP=2 on subsequent time steps.
C
C     BMG_IOFLAG   Logical array of I/O and debugging flags. Admittedly
C                  this is not well documented yet, but look in the include
C                  file "BMG_parameters.h" for hints (the names are verbose.)
C
C     IRELAX   Indicator for the relaxation algorithm:
C     
C              IRELAX=1  => Colored Gauss-Seidel relaxation:
C                        -  red-black on the finest grid for 5-point stencils
C                        -  four color otherwise
C
C              IRELAX=5  => alternating red-black plane Gauss-Seidel relaxation
C                        -  ordering is xy planes, yz planes, xz planes
C                        -  for symmetric relaxation this is reversed
C                           as necessary when moving "up" or "down" levels
C
C     IRELAX_SYM  Symmetric cycles require that the relaxation ordering
C                 be reversed between the restriction and interpolation
C                 stage of the cycle.   This is important when BMG3D
C                 is being used as a preconditioner.
C              
C                    IRELAX_SYM = BMG_RELAX_NONSYM
C                    => use the stated relaxation ordering (nonsymmetric)
C
C                    IRELAX_SYM = BMG_RELAX_SYM
C                    => make the cycle symmetric.
C
C                 as defined in BMG_parameters.h
C
C     ISEP     Indicator for separable equation such that the SO array can be
C              reduced to 3 1-D arrays in the interior, saving memory bandwidth
C
C               ISEP = BMG_NONSEPARABLE (0): Always use 3d SO arrays
C               ISEP = BMG_SEPARABLE (1): only use 3d SO arrays on boundaries.
C
C
C     IVW      Indicator for type of cycle. IVW=1 means v-cycles will
C              be performed. IVW=2 means w-cycles will be performed, etc.
C
C     MCYCL    Currently disabled!
C
C     NFm      Maximum storage for a vector on all grids. This should be
C              computed by BMG3_SymStd_SETUP_space and is checked in boxmg.
C
C     NCm      Maximum storage for a vector on all coarse grids.  This should
C              be computed in BMG3_SymStd_SETUP_space and is checked in boxmg.
C
C     NSOm     See INPUT/OUTPUT
C
C     NSORm    See INPUT/OUTPUT
C
C     NCIm     See INPUT/OUTPUT
C     
C     NOGm     Maximum number of grids that can be supported by the 
C              pointer array IGRD. It should be set in calling program.
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     SO       User defined real array that contains the coefficient
C              matrix. See above description for format.  If SETUP is
C              specified then coarse-grid operators are created and returned. 
C
C     NSOm     Dimension of SO, the calling program should have computed
C              this with a call to BMG2_SymStd_SETUP_space.  Clearly,
C            
C              IFD .EQ. 1  =>    NSOm = 4*NFm + 10*NCm 
C              IFD .NE. 1  =>    NSOm = 14*NFm 
C
C              BMG3_SymStd_SOLVE_boxmg verifies that NSOm is large enough.
C
C     CI       CI is a real array which should be dimensioned
C              to (NCIm) in calling program. It is used to contain
C              the interpolation coefficients.
C
C     NCIm     Dimension of CI, set in calling program.  NCIm = 26*NCm is 
C              enough.  BOXMG checks to see if NCIm is large enough.
C
C     QF       The user defined array that contains the right hand side.
C              It must be dimensioned to at least NFMAX.
C
C     Q        the user defined array that contains the solution vector.
C              It must be dimensioned to at least NFMAX.
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C
C =======================================================================
C  --------------------
C   WORK ARRAYS:
C  --------------------
C
C     ABD      User declared two dimensional real array, which
C              is used to store the coefficient matrix for the coarsest
C              grid. It is then used by the LAPACK routines. It
C              should be dimensioned to (NCBWm,NCUm).
C
C     BBD      User declared real array of dimension NCUm for use in the
C              the LAPACK routines
C
C     NCBWm    Maximum first subscript of ABD, which needs to be > or =
C              (x-dimension+1)*(y-dimension+2) on coarsest grid. 
C
C     NCUm     Maximum second subscript of ABD, which needs to be > or =
C              (x-dimension+1)*(y-dimension+1)*(z-dimension+1) on
C              coarsest grid.
C
C     IGRD     A work array. IGRD should be dimensioned to
C              IGRD(NOGm,NBMG_pIGRD) in the calling program.
C
C     SOR      SOR is a real array, which should be dimensioned to
C              (NSORm) in the calling program. It is used to store
C              residuals, reciprocals of SO(.,1) if IRELAX=1, and lu
C              decompositions if IRELAX.gt.1.
C
C     NSORm    Dimension of SOR, set in calling program. NSORm=2*NFMAX
C              is enough for IRELAX=1 or 2. For IRELAX=3 or 4, 4*NFMAX
C              is enough. 
C
C
C =======================================================================
C  --------------------
C   REFERENCES:
C  --------------------
C
C     Dendy, J. E. Jr., "Two Multigrid Methods for Three Dimensional
C     Problems with Discountinuos and Anisotropic Coefficients", SIAM
C     Journal of Scientific and Satatistical Computing, Vol. 8, No. 2,
C     September 1987
C     
C ==========================================================================
C  --------------------
C   Error Codes:
C  --------------------
C
C     -------------------------------------------------------------------
C     CODE    DESCRIPTION                           ORIGIN
C     -------------------------------------------------------------------
C      1     initial residual has become zero     SOLVE_boxmg
C      2     NOG = 0                              SETUP_parts
C      3     NOG < 0                              SETUP_parts
C      4     NStncl out of range                  SETUP_cg_LU
C      5     Cholesky decomposition failed, look
C            in BMG_iPARMS(id_BMG3_Ext_Err_Code)
C            for the return code of the LAPACK
C            routine                              SETUP_cg_LU
C      6    NOG = 1                               SETUP_relax
C      7    NOG = 0                               SETUP_relax
C      8    NOG < 0                               SETUP_relax
C      9    fatal setup error                     SETUP_PtrGrid
C     10    IRELAX out of range                   SETUP_PtrGrid
C     11    min coarse grid dim < 3 ( NXYZc )     SETUP_space
C     12    min coarse grid dim < 3 ( NXYc )      SETUP_space
C     13    IRELAX out of range                   SETUP_space
C     14    computed number of grids too small    SETUP_space
C     15    computed number of grids < 1          SETUP_space
C     16    fatal error (check output)            SETUP_PtrWork
C     17    memory allocation mode unspecified    SETUP_PtrWork
C     18    KCF = KC                              ncycle
C     19    updown out of range                   updown
C     20    coarse grid solve failed, look in
C           BMG_iPARMS(id_BMG3_Ext_Err_Code)
C           for the return code of the LAPACK
C           routine                               SOLVE_cg
C     21    inconsistent # of stencil entries     COPY_SO_xy
C     22    inconsistent # of stencil entries     COPY_SO_yz
C     23    inconsistent # of stencil entries     COPY_SO_xz
C     24    invalid cg_operator type              SETUP_parts
C     25    invalid cg_operator construction      SETUP_parts
C
C ==========================================================================

      IMPLICIT NONE

C ------------------------------------------------
C     Includes
C
      INCLUDE 'geom_param_fort.h'
      INCLUDE 'mpi_param_fort.h'
      INCLUDE 'MSG.h'

#include    "BMG_constants_f77.h"
      INCLUDE 'BMG_stencils_f77.h'
      INCLUDE 'BMG_workspace_f77.h'
      INCLUDE 'BMG_parameters_f77.h'

C ------------------------------------------------
C     Argument Declarations
C
      INTEGER  NLx, NLy, NLz, NGx, NGy, NGz, iGs, jGs, kGs,
     &         NFm, NCm, NSOm, NSORm, NCIm, NCBWm, NCUm, NOGm, NOG,
     &         NMSGim, NMSGrm, NProc, MPICOMM, NBMG_MSG_iGRID,
     &         NBMG_iWORK_CSm, NBMG_iWORK_PLm,
     &         NBMG_rWORK_CSm, NBMG_rWORK_PLm

      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER  iWorkMSG(NMSGim),
     &         pMSG(NBMG_pMSG,NOGm),
     &         pMSGSO(NBMG_pMSG,NOGm)
      REAL*RKIND   MSG_Buffer(NMSGrm)

      INTEGER  BMG_iWORK_PL(NBMG_iWORK_PLm)
      REAL*RKIND   BMG_rWORK_PL(NBMG_rWORK_PLm)

      INTEGER  BMG_iWORK_CS(NBMG_iWORK_CSm)
      REAL*RKIND   BMG_rWORK_CS(NBMG_rWORK_CSm)

      INTEGER  IGRD(NOGm,NBMG_pIGRD)
      REAL*RKIND CI(NCIm),Q(NFm),QF(NFm),RES(NFm),SO(NSOm),SOR(NSORm)


C ------------------------------------------------
C     Local Declarations
C
      INTEGER i, iBC, id, ifd, iil, iilc, iilf, im, irelax, isep,
     &        IRELAX_SYM, ISETUP, ISTOP, istrt, istrt2, iu, ivw,
     &        jjl, jjlc, jjlf, k, kc, kfmg, kf, kkl, kklc, kklf,
     &        mcyc, mcycl, ncyc, NStncl, NMSGi, NMSGr, p_ABD, p_BBD,
     &        p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC,
     &        IIG, IIGC, IIGF, JJG, JJGC, JJGF, KKG, KKGC, KKGF,
     &        NC, NCBW, NCI, NCU, NF, NSO, NSOR, NXYZc,
     &        NBMG_iWORK_CS, NBMG_iWORK_PL, NBMG_rWORK_CS, 
     &        NBMG_rWORK_PL


      INTEGER IMdum

      INTEGER  iGs_c, iGs_f, jGs_c, jGs_f, kGs_c, kGs_f

      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK,
     &         NProcI, NProcJ, NProcK,
     &         p_NLx, p_NLy, p_NLz,
     &         p_NLx_kg, p_NLy_kg, p_NLz_kg,
     &         p_ProcCoord, p_ProcGrid, pSI_MSG

      INTEGER  BMG_MSG_MyProc, MPI_IERR

      REAL*8  REL_RES_L2, RES_L2, RES_L2_0
      REAL*8 T, T1, T2, TOL
      REAL*8  TT1, TT2
      LOGICAL NCYCLE_FLAG, INIT_CHECK, FINE_ONLY
      
      logical, parameter :: debug = .false.

C ==========================================================================

      
      IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
         RETURN
      END IF

C -------------------
C     Zero Times:
C -------------------

      t=dZERO
      t1=dZERO
      t2=dZERO

      IF ( debug )  write(0,*) 'myProc, start commas_SymStd_SOLVE_boxmgdg =',myProc
      IF ( BMG_IOFLAG(iBMG3_BUG_PARAMETERS) ) THEN
         CALL BMG3_SymStd_DUMP_parms( 
     &                    BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &                    )
      ENDIF

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------
C --------------------------------------------
C     MPI Parallel Data Decomposition:
C --------------------------------------------
      
      !
      ! Processor Grid:
      !
      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_NLx       = BMG_MSG_pGRID(ip_BMG_MSG_NLx_Grid)
      p_NLy       = BMG_MSG_pGRID(ip_BMG_MSG_NLy_Grid)
      p_NLz       = BMG_MSG_pGRID(ip_BMG_MSG_NLz_Grid)
      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)
      
C ----------------------------------------
C     MSG workspace pointers:
C ----------------------------------------

      p_NLx_kg = 1
      p_NLy_kg = p_NLx_kg + NProcI*NOGm
      p_NLz_kg = p_NLy_kg + NProcJ*NOGm

      pSI_MSG  = p_NLz_kg + NProcK*NOGm  !!NB: initalization of pSI_MSG

C ----------------------------------------
C     Multigrid Cycle parameters:
C ----------------------------------------

      IFD    = BMG_iPARMS(id_BMG3_STENCIL)
      IBC    = BMG_iPARMS(id_BMG3_BC)
      ISETUP = BMG_iPARMS(id_BMG3_SETUP)

      IRELAX = BMG_iPARMS(id_BMG3_RELAX)
      IRELAX_SYM = BMG_iPARMS(id_BMG3_RELAX_SYM )

      ISEP = BMG_iPARMS(id_BMG3_SEPARABLE)

      ID = BMG_iPARMS(id_BMG3_NRELAX_DOWN)
      IU = BMG_iPARMS(id_BMG3_NRELAX_UP)
      IM = BMG_iPARMS(id_BMG3_NRELAX_FG)
      IMdum = 0

      IF (BMG_iPARMS(id_BMG3_CYCLE_CLASS).EQ.BMG_N_CYCLE) THEN
         ISTRT = - BMG_iPARMS(id_BMG3_MAX_ITERS)
      ELSEIF (BMG_iPARMS(id_BMG3_CYCLE_CLASS).EQ.BMG_FMG_CYCLE) THEN 
         ISTRT = BMG_iPARMS(id_BMG3_MAX_ITERS)
      ENDIF

      IVW   = BMG_iPARMS(id_BMG3_NCYCLE_TYPE)
      MCYCL = BMG_iPARMS(id_BMG3_FMG_NNCYCLE)

      ISTOP = BMG_iPARMS(id_BMG3_STOP_TEST)
      TOL   = BMG_rPARMS(id_BMG3_STOP_TOL)

      IF ( BMG_IOFLAG(iBMG3_BUG_PARAMETERS) ) THEN
         write(0,*) 'myProc, nproci,nprocj,nprock =',myProc, nproci,nprocj,nprock
         
      ENDIF

C ================================================

      IF (BMG_iPARMS(id_BMG2_CYCLE_CLASS).EQ.BMG_N_CYCLE) THEN
         ISTRT2 = - BMG_iPARMS(id_BMG2_MAX_ITERS)
      ELSEIF (BMG_iPARMS(id_BMG2_CYCLE_CLASS).EQ.BMG_FMG_CYCLE) THEN 
         ISTRT2 = BMG_iPARMS(id_BMG2_MAX_ITERS)
      ENDIF

C ================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C ------------------------------------------------------------------- 
C     Enable MSG communications checking to make sure all initialized
C -------------------------------------------------------------------

      IF ( debug )  write(0,*) 'myProc, call MPI_INITIALIZED =',myProc

      CALL MPI_INITIALIZED(INIT_CHECK, MPI_IERR)
      IF( .NOT. INIT_CHECK ) THEN
         PRINT *,' MPI NOT INITIALIZED!'
         RETURN
      ENDIF
      
      IF( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR. 
     &    ISETUP.EQ.BMG_SETUP_only ) THEN

         CALL MSG_set_comm_parent(MPICOMM)
         CALL MSG_enable(BMG_MSG_MyProc,NProc)
         BMG_MSG_MyProc = MSG_myproc()
         NProc  = MSG_nproc()
            
      ELSE
         
         IF( MPICOMM .EQ. MSG_COMM_PARENT) THEN
            BMG_MSG_MyProc = MSG_myproc()
            NProc  = MSG_nproc()
         ELSE
            IF( MyProc.EQ.1 ) PRINT *,' MSG NOT INITIALIZED!'
         ENDIF

      ENDIF

C ---------------------------------
C     Initialize timers
C ---------------------------------

      T  = dZERO
      T1 = dZERO
      T2 = dZERO

C ----------------------------------
C     Sanity checks!
C ----------------------------------

      !!
      !!   We need to add way more some sanity checks here  !!!
      !!   Maybe a subroutine to check all parameter values??
      !!   

      !
      ! Processor indexing
      !
      IF ( BMG_MSG_MyProc.NE.BMG_MSG_iGRID(id_BMG_MSG_MyProc) ) THEN 
         WRITE(*,*) '**** ERROR: BMG3_SymStd_SOLVE_boxmg.f ... '
         WRITE(*,*) ' MyProc =         ', MyProc
         WRITE(*,*) ' BMG_MSG_MyProc = ', BMG_MSG_MyProc
         STOP
      ENDIF


C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> BEGIN:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      IF ( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR. 
     &     ISETUP.EQ.BMG_SETUP_only .OR. 
     &     ISETUP.EQ.BMG_SETUP_ptrs ) THEN

         TT1 = MPI_Wtime()

      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_SETUP_PtrGrid =',myProc

         CALL BMG3_SymStd_SETUP_PtrGrid( 
     &             NLx, NLy, NLz, NGx, NGy, NGz,
     &             iGs, jGs, kGs, 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,    
     &             MyProc, NProc, NProcI, NProcJ, NProcK, 
     &             MPICOMM, BMG_iWORK_PL, BMG_iWORK_CS, 
     &             NOGm, NFm, NCm, NSOm, NSORm, NCIm,
     &             NBMG_iWORK_PLm, NBMG_rWORK_PLm,
     &             NBMG_iWORK_CSm, NBMG_rWORK_CSm, NMSGim, NMSGrm,
     &             NOG, NF, NC, NSO, NSOR, NCI, IGRD,
     &             NBMG_iWORK_PL, NBMG_rWORK_PL,
     &             NBMG_iWORK_CS, NBMG_rWORK_CS, NMSGi, NMSGr
     &             )
         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF


         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_SETUP_PTR_GRID) =
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_PTR_GRID) + TT2 - TT1

         ! Output      
         IF( BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE) ) THEN
            WRITE (*,260) 
     &            'Storage for a vector on all grids = ', NF
         ENDIF

         ! Output
         IF( BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE) ) THEN
            WRITE (*,270) 
     &            'Storage for a vector on all coarse grids = ', NC
         ENDIF

         !
         ! Setup of the MSG arrays (pSI_MSG initialized above).
         !
         TT1 = MPI_Wtime()

      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_SETUP_MSG =',myProc

         CALL BMG3_SymStd_SETUP_MSG( 
     &             pMSG, pMSGSO, iWorkMSG, NMSGi, pSI_MSG, 
     &             IGRD, NOG, NOGm, NProc, MyProc, 
     &             iWorkMSG(p_NLx_kg), iWorkMSG(p_NLy_kg),
     &             iWorkMSG(p_NLz_kg), 
     &             BMG_MSG_iGRID(p_NLx), BMG_MSG_iGRID(p_NLy),
     &             BMG_MSG_iGRID(p_NLz),
     &             BMG_MSG_iGRID(p_ProcGrid),
     &             NProcI, NProcJ, NProcK,
     &             MPICOMM 
     &             )

         TT2 = MPI_Wtime()

         BMG_rPARMS(id_BMG3_TIME_SETUP_MSG) = 
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_MSG) + TT2 - TT1
     
      ELSE

         !
         !  Maximum dimensions are the actual dimensions.
         !

         NF   = NFm
         NC   = NCm
         NSO  = NSOm
         NSOR = NSORm
         NCI  = NCIm

         NBMG_iWORK_CS = NBMG_iWORK_CSm
         NBMG_rWORK_CS = NBMG_rWORK_CSm
 
         NMSGi = NMSGim
         NMSGr = NMSGrm

      ENDIF

      !
      ! Create local pointers for LU solve
      !
      IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN

         p_ABD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO)
         p_BBD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CU)
         
         NCBW  = BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW)
         NCU   = BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)

      ENDIF

      !
      !  No components just setup pointers, thanks!
      !
      IF ( ISETUP.EQ.BMG_SETUP_ptrs ) RETURN

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>> END:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>> BEGIN:  COMPONENT SETUP <<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------


 89   CONTINUE  ! >>>>>>>>> BEGIN: Component Setup

      IF ( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR. 
     &     ISETUP.EQ.BMG_SETUP_opers .OR. 
     &     ISETUP.EQ.BMG_SETUP_only ) THEN

         !  Start the timer
         t1 = MPI_Wtime()

         !
         ! Construct coarse-grid and interpolation operators.  In 
         ! addition, setup any necessary components for relaxation.
         !
         TT1 = MPI_Wtime()
         
      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_SETUP_parts =',myProc

         CALL BMG3_SymStd_SETUP_parts( 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm, NOG,
     &             BMG_iWORK_PL, NBMG_iWORK_PL,
     &             BMG_rWORK_PL, NBMG_rWORK_PL,
     &             BMG_iWORK_CS, NBMG_iWORK_CS,
     &             BMG_rWORK_CS, NBMG_rWORK_CS,
     &             iWorkMSG, NMSGi, pSI_MSG, pMSG, pMSGSO,
     &             MSG_Buffer, NMSGr,
     &             MyProcI, MyProcJ, MyProcK,
     &             NProcI, NProcJ, NProcK, NProc,
     &             BMG_MSG_iGRID(p_ProcGrid),
     &             BMG_MSG_iGRID(p_ProcCoord),
     &             iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &             iWorkMSG(p_NLx_kg), iWorkMSG(p_NLy_kg),
     &             iWorkMSG(p_NLz_kg),
     &             MPICOMM
     &             )

      IF ( debug )  write(0,*) 'myProc, done BMG3_SymStd_SETUP_parts =',myProc
      IF (debug) CALL MPI_Barrier(mpicomm, i)

         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_SETUP_PARTS) =
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_PARTS) + TT2 - TT1

         ! Compute the setup time
         t2 = MPI_Wtime()
         t=t+t2-t1
         BMG_rPARMS(id_BMG3_TIME_SETUP_TOTAL) =
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_TOTAL) + t2 - t1

         
         ! Output the setup time
         IF( BMG_IOFLAG(iBMG3_OUT_TIME_SETUP) )  THEN
            WRITE (*,220) '(3D) SETUP TIME =', T
         ENDIF

      ENDIF

      !
      !  No solve just the setup, thanks!  
      !
      IF (ISETUP.EQ.BMG_SETUP_only) RETURN

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  COMPONENT SETUP <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN: MULTIGRID CYCLING <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

 95   CONTINUE

      ! Start the timer
      t1 = MPI_Wtime()

C ----------------------------------------------------------------
C     Initial residual:
C ----------------------------------------------------------------

      !
      ! Get pointers for the finest grid (GRID LEVEL = NOG)
      !
      IF (debug) CALL MPI_Barrier(mpicomm, i)
      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_GET_pointers =',myProc

      CALL BMG3_SymStd_GET_pointers( 
     &          NOG, IGRD, NOGm, 
     &          p_U, p_SO, p_SOR, p_CI, 
     &          IIL, JJL, KKL, IIG, JJG, KKG
     &          )

C      FIXME: The arguments passed in should be correct, but this 
C             is a reasonable place to double check
C
C      iGs = IGRD(NOG,idL_BMG_Icoord)
C      jGs = IGRD(NOG,idL_BMG_Jcoord)
C      kGs = IGRD(NOG,idL_BMG_Kcoord)

      !
      ! Set NStncl for the finest grid
      !
      IF ( IFD.NE.1 ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF
     
      !
      ! if requested, sync the initial guess
      !
      IF ( BMG_iPARMS(id_BMG3_SYNC_INITIAL_GUESS) 
     &     .eq. BMG_SYNC_INITIAL_GUESS ) THEN
         
      IF (debug) CALL MPI_Barrier(mpicomm, i)
      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_UTILS_update_ghosts =',myProc

         CALL BMG3_SymStd_UTILS_update_ghosts(
     &             NOG, Q(p_U), IIL, JJL, KKL, iWorkMSG, pMSG,
     &             MSG_Buffer, MPICOMM
     &             )

      END IF

      IF (debug) CALL MPI_Barrier(mpicomm, i)
      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_residual =',myProc

      CALL BMG3_SymStd_residual(
     &          NOG, NOGm, IFD,
     &          Q(p_U), QF(p_U), SO(p_SO), RES(p_U), 
     &          IIL, JJL, KKL, NStncl,
     &          iWorkMSG, NMSGi, pMSG, 
     &          MSG_Buffer, NMSGr, MPICOMM
     &          )

      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_UTILS_norm_l2 =',myProc

      CALL BMG3_SymStd_UTILS_norm_l2( 
     &           RES(p_U), IIL, JJL, KKL, RES_L2_0, MPICOMM 
     &           )

      !
      !  Output
      !
      IF ( BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
         IF ( IM == 0 ) THEN
           WRITE (*,405) '*** THE INITIAL RESIDUAL DG (l2-NORM) = ', RES_L2_0
         ELSE
           WRITE (*,405) '*** THE INIT RES DG-FG = ', RES_L2_0
         ENDIF
      ENDIF

      IF ( RES_L2_0 .EQ. dZERO ) THEN
         IF ( BMG_IOFLAG(iBMG3_WARN_ZERO_RESIDUAL) ) THEN
            WRITE(*,505) '*** WARNING: bmg3d.f '
            WRITE(*,505) '    Zero initial residual !!!! '
            WRITE(*,505) '    D''ooh! Initial guess is the solution!'
         ENDIF
         RETURN
      ENDIF

C ----------------------------------------------------------------
C     Direct solve in the degenerate case of one grid!
C ----------------------------------------------------------------

      IF ( NOG.EQ.1 ) THEN

         !
         ! Solve on the coarsest grid
         !
         IF ( MyProc.EQ.1 ) THEN 
            WRITE(*,*) '*** Performing direct solve on the fine grid!'
         ENDIF
         !

         IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN
            !
            TT1 = MPI_Wtime()
            !
            CALL BMG3_SymStd_SOLVE_cg_LU( 
     &                q(p_U), qf(p_U), IIL, JJL, KKL,
     &                BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD), 
     &                NCBW, NCU, NOGm,
     &                NProcI, NProcJ, NProcK, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM, BMG_iPARMS
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ELSE
            !
            TT1 = MPI_Wtime()
            !
            
            CALL BMG3_SymStd_SOLVE_cg_boxmg( 
     &                Q(p_U), QF(p_U), NOG,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                IIL, JJL, KKL, IIG, JJG, KKG,
     &                iGs, jGs, kGs,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                NProcI, NProcJ, NProcK, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ENDIF

         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         BMG_rPARMS(id_BMG3_TIME_SOLVE_CG)  
     &       = BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) + TT2 - TT1

         !
         ! Compute the final residual
         !
         CALL BMG3_SymStd_residual( 
     &             NOG, NOG, IFD,
     &             Q(p_U), QF(p_U), SO(p_SO), RES(p_U),
     &             IIL, JJL, KKL, NStncl, 
     &             iWorkMSG, NMSGi, pMSG, 
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )
         !
         CALL BMG3_SymStd_UTILS_norm_l2( 
     &             RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM 
     &             )
         !
         IF( BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) .AND. MyProc.EQ.1 ) THEN 
            WRITE (*,405) '*** THE FINAL RESIDUAL (l2-NORM) = ', RES_L2
         ENDIF
         !
         TOL = RES_L2
         ! Jump to final I/O
         GOTO 150
         !
      ENDIF



C --------------------------------
C     Cycle Preliminaries:
C -------------------------------- 

      KF = NOG                  ! finest grid index
      KC = 1                    ! coarsest grid index
      KFMG = 1                  ! FMG current finest grid
      MCYC = ABS(ISTRT)         ! Maximum number of multigrid cycles
      NCYC = 1                  ! multigrid cycle counter
      NCYCLE_FLAG = .TRUE.      ! default to do n-cycles
      FINE_ONLY = .FALSE.       ! run relaxation on finest grid only
      
      IF ( IM .gt. 0 ) FINE_ONLY = .TRUE.

C ----------------------------------------------------------------
C     Perform multigrid F-cycle:
C     (if necessary, then continue with n-cycles)
C ----------------------------------------------------------------

      IF ( ISTRT.GT.0 .and. .not. FINE_ONLY ) THEN

         ! Coarsen the right hand side
         DO k=KF-1, 1, -1

            CALL BMG3_SymStd_GET_pointers( k+1, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, 
     &                IILF, JJLF, KKLF, IIGF, JJGF, KKGF 
     &                )
            iGs_f = IGRD(k+1,idL_BMG_Icoord)
            jGs_f = IGRD(k+1,idL_BMG_Jcoord)
            kGs_f = IGRD(k+1,idL_BMG_Kcoord)

            CALL BMG3_SymStd_GET_pointers( k, IGRD, NOGm,
     &                p_UC, p_SOC, p_SORC, p_CIC, 
     &                IILC, JJLC, KKLC, IIGC, JJGC, KKGC  
     &                )
            iGs_c = IGRD(k,idL_BMG_Icoord)
            jGs_c = IGRD(k,idL_BMG_Jcoord)
            kGs_c = IGRD(k,idL_BMG_Kcoord)

            ! restrict the right hand side 

            TT2 = MPI_Wtime()

            CALL BMG3_SymStd_restrict( 
     &                K+1, K,  
     &                QF(p_U), QF(p_UC), CI(p_CIC),
     &                IILF, JJLF, KKLF, IILC, JJLC, KKLC,
     &                iGs_f, jGs_f, kGs_f
     &                )

            TT1 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_restrict) =
     &           BMG_rPARMS(id_BMG3_TIME_restrict) + TT2 - TT1

            ! Copy right hand side into RES(i,j,k) for interpolation
            DO i=0, IILF*JJLF*KKLF-1
               RES(p_U+i)=QF(p_U+i)
            END DO
         END DO


         ! FMG current finest grid index is set to the coarsest grid
         KFMG=KC     
         ! FMG: solve exactly on the coarsest grid
         CALL BMG3_SymStd_GET_pointers(
     &             KFMG, IGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI, 
     &             IIL, JJL, KKL, IIG, JJG, KKG 
     &             )
         iGs_c = IGRD(KFMG,idL_BMG_Icoord)
         jGs_c = IGRD(KFMG,idL_BMG_Jcoord)
         kGs_c = IGRD(KFMG,idL_BMG_Kcoord)
         
         IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU ) THEN
            !
            TT1 = MPI_Wtime()
            !
            CALL BMG3_SymStd_SOLVE_cg_LU( 
     &                q(p_U), qf(p_U), IIL, JJL, KKL, 
     &                BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD), 
     &                NCBW, NCU, NOGm,
     &                NProcI, NProcJ, NProcK, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM, BMG_iPARMS 
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ELSE
            !
            TT1 = MPI_Wtime()
            !

            CALL BMG3_SymStd_SOLVE_cg_boxmg( 
     &                Q(p_U), QF(p_U), NOG,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                IIL, JJL, KKL, IIG, JJG, KKG,
     &                iGs_c, jGs_c, kGs_c,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                NProcI, NProcJ, NProcK, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ENDIF

         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         BMG_rPARMS(id_BMG3_TIME_SOLVE_CG)  
     &       = BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) + TT2 - TT1
         
         RES_L2=dZERO
         IF ( BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
            WRITE (*,230) KFMG, RES_L2
         ENDIF

 100     CONTINUE    ! >>>>>>>> LOOP BOUNDARY: multigrid F-cycle


            KFMG=KFMG+1
            ! 
            ! FMG: Interpolate
            !
            CALL BMG3_SymStd_GET_pointers(
     &                KFMG, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, IIL, JJL, KKL,
     &                IIG, JJG, KKG 
     &                )
            iGs_f = IGRD(KFMG,idL_BMG_Icoord)
            jGs_f = IGRD(KFMG,idL_BMG_Jcoord)
            kGs_f = IGRD(KFMG,idL_BMG_Kcoord)

            CALL BMG3_SymStd_GET_pointers(
     &                KFMG-1, IGRD, NOGm,
     &                p_UC, p_SOC, p_SORC, p_CIC, IILC, JJLC, KKLC,
     &                IIGC, JJGC, KKGC  
     &                )
            iGs_c = IGRD(KFMG-1,idL_BMG_Icoord)
            jGs_c = IGRD(KFMG-1,idL_BMG_Jcoord)
            kGs_c = IGRD(KFMG-1,idL_BMG_Kcoord)

            !
            ! Set NStncl for grid KFMG
            !
            IF ( KFMG.LT.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
               NStncl=14
            ELSE
               NStncl=4
            ENDIF
            !


            TT1 = MPI_Wtime()

            CALL BMG3_SymStd_interp_add( 
     &                KFMG-1, KFMG, NOG,
     &                Q(p_U) ,Q(p_UC), RES(p_U),
     &                SO(p_SO), NStncl, CI(p_CIC),
     &                IILC, JJLC, KKLC,
     &                IIL, JJL, KKL,
     &                iGs_f, jGs_f, kGs_f,
     &                iWorkMSG, NMSGi, pMSG,
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            

            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_interp_add) =
     &           BMG_rPARMS(id_BMG3_TIME_interp_add) + TT2 - TT1
            !
            ! FMG; Perform an n-cycle 
            !
            IF ( BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN 
               WRITE(*,410) '*** Performing n-cycle:  Coarsest =', KC,
     &                      'Finest =', KFMG
            ENDIF


            CALL commas_SymStd_ncycledg( 
     &                KC, KFMG, KF, 
     &                IFD, IU, ID, IMdum, IVW, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm,
     &                SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID
     &                ) ! ,Qmask)


            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            IF ( KFMG.EQ.KF ) THEN
               !
               !  Compute the final residual
               !
               CALL BMG3_SymStd_residual(
     &                   KF, KF, IFD,
     &                   Q(p_U), QF(p_U), SO(p_SO), RES(p_U),
     &                   IIL, JJL, KKL, NStncl,
     &                   iWorkMSG, NMSGi, pMSG, 
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )
               CALL BMG3_SymStd_UTILS_norm_l2( 
     &                   RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM
     &                   )

               ! Note this is still vulnerable to very small RES_L2_0
               IF ( RES_L2_0.EQ.dZERO ) THEN
                  IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR).AND.MyProc.EQ.1)
     &                 THEN
                    WRITE(*,505)
     &                    '*** FATAL ERROR: BMG3_SymStd_SOLVE_boxmg.f '
                    WRITE(*,505)'    Initial residual has become zero!'
                  END IF
                  
                  CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,1)
                  RETURN
                  
               ELSE
                  REL_RES_L2 = RES_L2/RES_L2_0
               ENDIF             

               !
               !  Output 
               !
               IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )
     &            .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
                  WRITE (*,400) '*** ITERATION =', NCYC, 
     &                          '*** RELATIVE RESIDUAL = ', REL_RES_L2
               ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 ) 
     &                 .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
                  WRITE (*,400) '*** ITERATION =', NCYC, 
     &                          '*** ABSOLUTE RESIDUAL = ', RES_L2
               ENDIF
               !
               !  Check Convergence
               !
               IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 ) 
     &            .AND. ( REL_RES_L2.LT.TOL )     ) THEN
                  ! FMG: converged
                  !      - set TOL and skip n-cycles
                  TOL = REL_RES_L2
                  NCYCLE_FLAG=.FALSE.
               ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 ) 
     &                 .AND. ( RES_L2.LT.TOL )         ) THEN
                  ! FMG: converged
                  !      - set TOL and skip n-cycles
                  TOL = RES_L2
                  NCYCLE_FLAG=.FALSE.
               ELSE IF ( NCYC.GE.MCYC ) THEN
                  ! FMG: cycle limit reached 
                  !      - set TOL and skip n-cycles
                  IF ( ISTOP.EQ.BMG_STOP_REL_RES_L2 ) THEN
                     TOL = -REL_RES_L2
                  ELSE
                     TOL = -RES_L2
                  ENDIF
                  NCYCLE_FLAG=.FALSE.
               ELSE 
                  ! FMG: didn't converge
                  !      - increase counter and continue with n-cycles
                  NCYC=NCYC+1
               ENDIF
               !
            ELSE
               ! FMG: Continue the F-cycle
               GOTO 100
            ENDIF
            
 120     CONTINUE    ! <<<<<<<< LOOP BOUNDARY: multigrid F-cycle

      ENDIF  ! ISTRT > 0

C ----------------------------------------------------------------
C     Perform IM relaxations on fine grid
C ----------------------------------------------------------------

      IF ( FINE_ONLY ) THEN
        NCYCLE_FLAG = .FALSE.

         CALL BMG3_SymStd_GET_pointers(
     &             KF, IGRD, NOGm, p_U, p_SO, 
     &             p_SOR, p_CI, IIL, JJL, KKL,
     &             IIG, JJG, KKG
     &             )

            CALL commas_SymStd_ncycledg( 
     &                KC, KF, KF, 
     &                IFD, IU, ID, IM, IVW, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm,
     &                SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr, 
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID
     &                ) ! ,Qmask)

            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            !
            !  Compute the final residual
            !
            CALL BMG3_SymStd_residual(
     &                KF, KF, IFD,
     &                Q(p_U), QF(p_U), SO(p_SO), RES(p_U), 
     &                IIL, JJL, KKL, NStncl, 
     &                iWorkMSG, NMSGi, pMSG, 
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            CALL BMG3_SymStd_UTILS_norm_l2( 
     &                RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM
     &                )

            ! Note this is still vulnerable to very small RES_L2_0
            IF ( RES_L2_0.EQ.dZERO ) THEN
               WRITE(*,505) '*** FATAL ERROR: bmg3d.f '
               WRITE(*,505) '    Initial residual has become zero!!!'
               STOP
            ELSE
               REL_RES_L2 = RES_L2/RES_L2_0
            ENDIF

            !
            !  Output 
            !
            IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )
     &           .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** FG_ONLY ITERATION =', NCYC, 
     &                       '*** FG_ONLY RELATIVE RESIDUAL = ', REL_RES_L2
            ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 )
     &              .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** FG_ONLY ITERATION =', NCYC, 
     &                      '*** FG_ONLY ABSOLUTE RESIDUAL = ', RES_L2
            ENDIF
      
      ENDIF

C ----------------------------------------------------------------
C     Perform multigrid n-cycle(s):
C ----------------------------------------------------------------

      IF ( NCYCLE_FLAG ) THEN
         
      IF ( debug )  write(0,*) 'myProc, call BMG3_SymStd_GET_pointers =',myProc

         CALL BMG3_SymStd_GET_pointers(
     &             KF, IGRD, NOGm, p_U, p_SO, 
     &             p_SOR, p_CI, IIL, JJL, KKL,
     &             IIG, JJG, KKG
     &             )
         

 130     CONTINUE   ! >>>>>>>> LOOP BOUNDARY: multigrid n-cycles


      IF ( debug )  write(0,*) 'myProc, call commas_SymStd_ncycledg =',myProc

            !
            ! Perform a multigrid n-cycle
            !
            CALL commas_SymStd_ncycledg( 
     &                KC, KF, KF, 
     &                IFD, IU, ID, IMdum, IVW, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm,
     &                SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr, 
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID
     &                ) ! ,Qmask)

            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            !
            !  Compute the final residual
            !
            CALL BMG3_SymStd_residual(
     &                KF, KF, IFD,
     &                Q(p_U), QF(p_U), SO(p_SO), RES(p_U), 
     &                IIL, JJL, KKL, NStncl, 
     &                iWorkMSG, NMSGi, pMSG, 
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            CALL BMG3_SymStd_UTILS_norm_l2( 
     &                RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM
     &                )

            ! Note this is still vulnerable to very small RES_L2_0
            IF ( RES_L2_0.EQ.dZERO ) THEN
               WRITE(*,505) '*** FATAL ERROR: bmg3d.f '
               WRITE(*,505) '    Initial residual has become zero!!!'
               STOP
            ELSE
               REL_RES_L2 = RES_L2/RES_L2_0
            ENDIF

            !
            !  Output 
            !
            IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )
     &           .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** ITERATION =', NCYC, 
     &                       '*** RELATIVE RESIDUAL = ', REL_RES_L2
            ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 )
     &              .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** ITERATION =', NCYC, 
     &                      '*** ABSOLUTE RESIDUAL = ', RES_L2
            ENDIF

            !
            ! Check Convergence 
            !
            IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )
     &         .AND. ( REL_RES_L2.LT.TOL )     ) THEN
               ! 
               ! n-cycles have converged in the RELATIVE RESIDUAL
               ! - set TOL and return
               !
               TOL = REL_RES_L2
            ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 ) 
     &              .AND. ( RES_L2.LT.TOL )         ) THEN
               ! 
               ! n-cycles have converged in the ABSOLUTE RESIDUAL
               ! - set TOL and return
               !
               TOL = RES_L2
            ELSE IF ( NCYC.GE.MCYC ) THEN
               !
               ! n-cycles failed to converge, cycle limit reached
               ! - set TOL and return
               !
               IF ( ISTOP.EQ.BMG_STOP_REL_RES_L2 ) THEN
                  TOL = -REL_RES_L2
               ELSE
                  TOL = -RES_L2
               ENDIF
            ELSE 
               !
               ! increase cycle count and continue n-cycling
               !
               NCYC=NCYC+1
               GOTO 130
            ENDIF

 140     CONTINUE    ! <<<<<<<< LOOP BOUNDARY: multigrid n-cycles

      ENDIF

 150  CONTINUE       ! <<<<<<<< FINAL I/O
      !
      !
      BMG_rPARMS(id_BMG3_STOP_TOL) = TOL
      !
      !  Compute the solve time and update the total time
      !
      t2 = MPI_Wtime()
      T=T+T2-T1
      BMG_rPARMS(id_BMG3_TIME_SOLVE_total) = 
     &     BMG_rPARMS(id_BMG3_TIME_SOLVE_total) + t2 - t1


      ! Output the multigrid cycling time
      IF ( BMG_IOFLAG(iBMG3_OUT_TIME_CYCLING) ) THEN
         WRITE(*,240) '(3D) MULTIGRID CYCLING TIME =', T2-T1

      ENDIF

      IF ( BMG_IOFLAG(iBMG3_OUT_TIME_TOTAL) ) THEN
         WRITE(*,240) '(3D) TOTAL TIME = ', T
      ENDIF

C -------------------------------------
C     Output:
C -------------------------------------         

C ==========================================================================

 220  FORMAT (/,1X,A,1X,F9.3,/)
 230  FORMAT ('(3D)LEVEL',I2,' RESIDUAL NORM= ',1P,E10.3)
 240  FORMAT (/,1X,A,1X,F12.3,/)
 260  FORMAT (/,/,2X,A,1X,I5)
 270  FORMAT (2X,A,1X,I5,/)

 400  FORMAT (1X,A,1X,I2,4X,A,1X,1P,E16.9)
 405  FORMAT (/,1X,A,1X,1P,E16.9,/)
 410  FORMAT (1X,A,1X,I2,4X,A,1X,I2)

 500  FORMAT (/,2X,A)
 505  FORMAT (5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
      
C ===========================================

      RETURN
      END

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE commas_SymStd_ncycledg( 
     &                KC, KCF, NOG, 
     &                IFD, IU, ID, IM, IVW, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC,
     &                SO, NSO, SOR, NSOR, CI, NCI, IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr, 
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID
     &                ) !,Qmask)

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   BMG3_SymStd_ncycle performs a single multigrid n-cycle.  
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'mpif.h'

#include    "BMG_constants_f77.h"
      INCLUDE 'BMG_workspace_f77.h'
      INCLUDE 'BMG_stencils_f77.h'
      INCLUDE 'BMG_parameters_f77.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  KC, KCF, NOG, 
     &         IFD, IU, ID, IM, IVW, IRELAX, IRELAX_SYM, ISEP,
     &         NF, NC, NSO, NSOR, NCI, NOGm, NMSGi, NMSGr,
     &         NBMG_MSG_iGRID, NBMG_iWORK_CS, NBMG_iWORK_PL,
     &         NBMG_rWORK_CS, NBMG_rWORK_PL

      INTEGER  BMG_iPARMS(NBMG_iPARMS),
     &         BMG_iWORK_CS(NBMG_iWORK_CS),
     &         BMG_iWORK_PL(NBMG_iWORK_PL), 
     &         iWorkMSG(NMSGi), pMSG(NBMG_pMSG,NOGm), 
     &         pMSGSO(NBMG_pMSG,NOGm), IGRD(NOGm,NBMG_pIGRD),
     &         BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS),
     &         BMG_rWORK_CS(NBMG_rWORK_CS),
     &         BMG_rWORK_PL(NBMG_rWORK_PL),
     &         CI(NCI), Q(NF), QF(NF), RES(NF),
     &         SO(NSO), SOR(NSOR), MSG_Buffer(NMSGr)
!      REAL*RKIND :: Qmask(NF)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, NProc, NProcI, NProcJ, 
     &         NProcK, MPICOMM, MPI_IERROR, NStncl

      INTEGER  iGs, jGs, kGs, NLx, NLy, NLz, NGx, NGy, NGz

      INTEGER  p_ProcGrid, p_ProcCoord

      INTEGER  i, j, k, p_ABD, p_BBD, p_CI, p_SO, p_SOR, p_U,
     &         NCBW, NCU
      REAL*8   RES_L2
      REAL*8  TT1, TT2

C ==========================================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------
C --------------------------------------------
C     MPI Parallel Data Decomposition:
C --------------------------------------------
      !
      ! MPI Communicators
      !
      MPICOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM)
      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      k = KCF   ! set the current grid index

      IF (KCF.EQ.KC) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,*) '*** FATAL ERROR: BMG3_SymStd_ncycle.f ****'
            WRITE(*,*) 'The current finest grid = ', KCF
            WRITE(*,*) 'The coarest grid        = ', KC
            WRITE(*,*) 'ERROR; minimum number of grids is 2 !'
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,18)
         RETURN

      ENDIF

      ! set the n-cycle counter to zero on all grids
      DO i = KC, KCF
         IGRD(i,idL_BMG_IVW)=iZERO
      END DO

      
C -----------------------------------------------------------
C     If IM > 0 then run fine grid relax only and return.
C -----------------------------------------------------------
      IF ( IM .gt. 0 ) THEN
      
            CALL commas_SymStd_updowndg(
     &                KCF, KCF, NOG, BMG_UP,
     &                IFD, IU, ID, IM, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID 
     &                ) ! ,Qmask)
            !
            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

      
      
      
         RETURN
      ENDIF


C -----------------------------------------------------------
C     Begin the n-cycle
C     (this is also the n-cycle recursive loop boundary)
C -----------------------------------------------------------

 140  CONTINUE     
         
         IF ( K.EQ.KC .AND. IGRD(k,idL_BMG_IVW).EQ.iZERO ) THEN
            !
            ! Solve on the coarsest grid
            !
            CALL BMG3_SymStd_GET_pointers( 
     &                KC, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, 
     &                NLx, NLy, NLz, NGx, NGy, NGz 
     &                )
            !
            iGs = IGRD(KC,idL_BMG_Icoord)
            jGs = IGRD(KC,idL_BMG_Jcoord)
            kGs = IGRD(KC,idL_BMG_Kcoord)
            !
            IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU ) THEN
               !
               p_ABD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO)
               p_BBD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CU)
               !
               NCBW  = BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW)
               NCU   = BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)
               !
               TT1 = MPI_Wtime()
               !
               CALL BMG3_SymStd_SOLVE_cg_LU( 
     &                   q(p_U), qf(p_U), NLx, NLy, NLz, 
     &                   BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD), 
     &                   NCBW, NCU, NOGm,
     &                   NProcI, NProcJ, NProcK, NProc, MyProc,
     &                   BMG_MSG_iGRID(p_ProcGrid),
     &                   BMG_MSG_iGRID(p_ProcCoord),
     &                   iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                   MSG_Buffer, NMSGr, MPICOMM, BMG_iPARMS
     &                   )
               !
               TT2 = MPI_Wtime()
               !
            ELSE
               !
               TT1 = MPI_Wtime()
               !
               CALL BMG3_SymStd_SOLVE_cg_boxmg( 
     &                   Q(p_U), QF(p_U), NOG,
     &                   BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                   NLx, NLy, NLz, NGx, NGy, NGz,
     &                   iGs, jGs, kGs,
     &                   BMG_iWORK_CS, NBMG_iWORK_CS,
     &                   BMG_rWORK_CS, NBMG_rWORK_CS,
     &                   NProcI, NProcJ, NProcK, NProc, MyProc,
     &                   BMG_MSG_iGRID(p_ProcGrid),
     &                   BMG_MSG_iGRID(p_ProcCoord),
     &                   iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )

               !
               TT2 = MPI_Wtime()
               !
            ENDIF

            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            BMG_rPARMS(id_BMG3_TIME_SOLVE_CG)
     &          = BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) + TT2 - TT1
            
            !
            IGRD(k,idL_BMG_IVW) = IVW
            IF( BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE)
     &          .AND. BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU
     &        ) THEN
               NStncl = 14   ! always the same on coarse grid if NOG > 1
               CALL BMG3_SymStd_residual( 
     &                   KC, NOG, IFD, 
     &                   q(p_U), qf(p_U), so(p_SO), RES(p_U),
     &                   NLx, NLy, NLz, NStncl,
     &                   iWorkMSG, NMSGi, pMSG,
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )               
               CALL BMG3_SymStd_UTILS_norm_l2( 
     &                   RES(p_U), NLx, NLy, NLz, RES_L2, MPICOMM
     &                   )
               !
               IF ( MyProc.EQ.1 ) THEN 
                  WRITE (*,100) k, RES_L2
               ENDIF
               !
            ENDIF
            !
            ! continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_IVW).LT.IVW ) THEN
            !
            ! Move down
            !            
            CALL commas_SymStd_updowndg( 
     &                k, KCF, NOG, BMG_DOWN, 
     &                IFD, IU, ID, IM, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID 
     &                ) ! ,Qmask)
            !
            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            IGRD(k,idL_BMG_IVW) = IGRD(k,idL_BMG_IVW) + 1 
            k=k-1
            !
            ! Continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_IVW).EQ.IVW ) THEN
            !
            ! Move up
            !
            IGRD(k,idL_BMG_IVW) = iZERO   ! reset n-cycle counter to zero
            k=k+1
            !
            CALL commas_SymStd_updowndg(
     &                k, KCF, NOG, BMG_UP,
     &                IFD, IU, ID, IM, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID 
     &                ) ! ,Qmask)
            !
            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            !
            ! Are we on the finest grid? 
            ! 
            IF (K.EQ.KCF) THEN 
               ! n-cycle is done
               RETURN
            ELSE
               ! continue the n-cycle
               GOTO 140  
            ENDIF
            !
         ENDIF

C ==========================================================================

 100     FORMAT (1X,'(3D)LEVEL ',I2,' RESIDUAL NORM = ',1P,E12.5)

C ===========================================

         RETURN
         END

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE commas_SymStd_updowndg ( 
     &                K, KF, NOG, UPDOWN,
     &                IFD, IU, ID, IM, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &                Q, QF, RES, NF, NC,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID 
     &                ) ! ,Qmask)


C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_updown performs all the necessary tasks associated
C     with moving "down" to a coarser grid or moving "up" to a finer
C     grid.  It is necessary because f77 does not support recursion.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C -----------------------------

      INCLUDE 'mpif.h'

#include    "BMG_constants_f77.h" 
      INCLUDE 'BMG_stencils_f77.h'
      INCLUDE 'BMG_parameters_f77.h'
      INCLUDE 'BMG_workspace_f77.h'

C ----------------------------
C     Argument Declarations
C ----------------------------

      !
      !  BoxMG parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      !
      !  BoxMG Cycle Parameters
      !
      INTEGER  k, KC, KCF, KF, IFD, IU, ID, IM, IVW,
     &         IRELAX, IRELAX_SYM, UPDOWN, ISEP

      !
      !  BoxMG dimensions
      !
      INTEGER  NF, NC, NSO, NSOR, NCI, NOGm, NOG,
     &         NBMG_iWORK_CS, NBMG_iWORK_PL,
     &         NBMG_rWORK_CS, NBMG_rWORK_PL

      !
      !  BoxMG internal pointers
      !
      INTEGER  IGRD(NOGm,NBMG_pIGRD)

      !
      !  BoxMG data arrays
      !
      REAL*RKIND   CI(NCI), Q(NF), QF(NF), RES(NF), SO(NSO), SOR(NSOR)
!      REAL*RKIND :: Qmask(NF)

      !
      !  BoxMG plane relaxation workspace
      !
      INTEGER BMG_iWORK_PL(NBMG_iWORK_PL)
      REAL*RKIND  BMG_rWORK_PL(NBMG_rWORK_PL)

      !
      !  MSG workspace dimensions
      !
      INTEGER  NMSGi, NMSGr
      !
      !  MSG workspace and pointers
      !
      INTEGER  iWorkMSG(NMSGi), 
     &         pMSG(NBMG_pMSG,NOG),
     &         pMSGSO(NBMG_pMSG,NOG)
      REAL*RKIND   MSG_Buffer(NMSGr)
 
      !
      !  BMG_MSG grid stuff
      !
      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)
      
C ----------------------------
C     Local Declarations
C ----------------------------

      INTEGER  I, J, NLx, NLy, NLz, NLx_c, NLy_c, NLz_c, 
     &         NGx, NGy, NGz, NGx_c, NGy_c, NGz_c, 
     &         p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC
      INTEGER  MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM,
     &         iGs, jGs, kGs, NSORv, NStncl
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, MPI_IERROR,
     &         NProc, NProcI, NProcJ, NProcK, p_ProcGrid, p_ProcCoord
      INTEGER  NPASS
      REAL*8   RES_L2
      REAL*8   TT1, TT2

C ==========================================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------
C --------------------------------------------
C     MPI Parallel Data Decomposition:
C --------------------------------------------
      !
      ! MPI Communicators
      !
      MPICOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM)
      !
      ! MPI Communicatorsn for plane relaxation
      !
      XYPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xy)
      XZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xz)
      YZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_yz)
 
      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      !
      ! Processor Coordinates
      !
      iGs = IGRD(k,idL_BMG_Icoord)
      jGs = IGRD(k,idL_BMG_Jcoord)
      kGs = IGRD(k,idL_BMG_Kcoord)
      
      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      !
      ! Set NStncl for the finest grid
      !
      IF ( K.LT.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         NSORv = 2
      ELSE 
         NSORv = 2  ! the same for now but this may change
      ENDIF

      IF ( IM .gt. 0 ) THEN
         DO i=1, IM
            !
            TT1 = MPI_Wtime()

            CALL commas_SymStd_relaxdg( 
     &                KF, KF, NOG, UPDOWN,
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,
     &                Q, QF, RES, NF, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iGs, jGs, kGs,
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
     &                ) ! ,Qmask)

            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_relax) = 
     &           BMG_rPARMS(id_BMG3_TIME_relax) + TT2 - TT1


            !
            IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) .AND. MyProc.EQ.1 ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
            !
          ENDDO
        
        RETURN
        
      ENDIF


      IF (UPDOWN.EQ.BMG_DOWN) THEN

         ! Get pointers for grids k and k-1
         CALL BMG3_SymStd_GET_pointers(     
     &               K, IGRD, NOGm,
     &               p_U, p_SO, p_SOR, p_CI, 
     &               NLx, NLy, NLz, NGx, NGy, NGz
     &               )
         CALL BMG3_SymStd_GET_pointers( 
     &               K-1, IGRD, NOGm,
     &               p_UC, p_SOC, p_SORC, p_CIC, 
     &               NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c 
     &               )
        
         ! Relaxation
         
         DO 10 i=1, ID
            !
            TT1 = MPI_Wtime()

            CALL commas_SymStd_relaxdg( 
     &                k, KF, NOG, UPDOWN,
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,
     &                Q, QF, RES, NF, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iGs, jGs, kGs,
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
     &                ) ! ,Qmask)

            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_relax) = 
     &           BMG_rPARMS(id_BMG3_TIME_relax) + TT2 - TT1


            !
            IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) .AND. MyProc.EQ.1 ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
            !
 10      CONTINUE

         ! restrict the residual
         TT1 = MPI_Wtime()

         ! calculate the residual
         CALL BMG3_SymStd_residual( 
     &             K, NOG, IFD, 
     &             Q(p_U), QF(p_U), SO(p_SO), RES(p_U),
     &             NLx, NLy, NLz, NStncl,
     &             iWorkMSG, NMSGi, pMSG, 
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) = 
     &        BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) + TT2 - TT1
         
c$$$         CALL BMG3_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, RES(p_U), K, NOG,
c$$$     &             NLx, NLy, NLz, NGx, NGy, NGz,
c$$$     &             IGRD(K,idL_BMG_ICoord),
c$$$     &             IGRD(K,idL_BMG_JCoord),
c$$$     &             IGRD(K,idL_BMG_KCoord),
c$$$     &             'DOWN-RES-post-relax', .FALSE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProcK, NProc,
c$$$     &             MyProcI, MyProcJ, MyProcK, MPICOMM
c$$$     &             )

c         ! restrict the residual
c         TT1 = MPI_Wtime()

         CALL BMG3_SymStd_restrict( 
     &             K, K-1,  
     &             RES(p_U), QF(p_UC), CI(p_CIC),
     &             NLx, NLy, NLz, NLx_c, NLy_c, NLz_c,
     &             iGs, jGs, kGs
     &             )

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_restrict) = 
     &        BMG_rPARMS(id_BMG3_TIME_restrict) + TT2 - TT1
         

c$$$         CALL BMG3_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, QF(p_UC), K-1, NOG,
c$$$     &             NLx_c, NLy_c, NLz_c, NGx_c, NGy, NGz,
c$$$     &             IGRD(K-1,idL_BMG_ICoord),
c$$$     &             IGRD(K-1,idL_BMG_JCoord),
c$$$     &             IGRD(K-1,idL_BMG_KCoord),
c$$$     &             'DOWN-QF-post-restrict', .FALSE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProcK, NProc,
c$$$     &             MyProcI, MyProcJ, MyProcK, MPICOMM
c$$$     &             )

         ! Zero the initial guess
         CALL BMG3_SymStd_UTILS_rV_zero( 
     &               k-1, q(p_UC), NLx_c, NLy_c, NLz_c
     &               )


      ELSE IF (UPDOWN.EQ.BMG_UP) THEN

         ! Get pointers for grids k and k-1
         CALL BMG3_SymStd_GET_pointers(      
     &             k, IGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI, 
     &             NLx, NLy, NLz, NGx, NGy, NGz
     &             )
         CALL BMG3_SymStd_GET_pointers( 
     &             k-1, IGRD, NOGm,
     &             p_UC, p_SOC, p_SORC, p_CIC, 
     &             NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c 
     &             )
         

         ! Interpolate and Correct 
         TT1 = MPI_Wtime()

            CALL BMG3_SymStd_interp_add( 
     &             k-1, k, NOG,
     &             Q(p_U) ,Q(p_UC), RES(p_U),
     &             SO(p_SO), NStncl, CI(p_CIC),
     &             NLx_c, NLy_c, NLz_c,
     &             NLx, NLy, NLz,
     &             iGs, jGs, kGs,
     &             iWorkMSG, NMSGi, pMSG,
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )
            

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_interp_add) =
     &        BMG_rPARMS(id_BMG3_TIME_interp_add) + TT2 - TT1

                  

         IF ( BMG_IOFLAG(iBMG3_BUG_RES_INTERP) ) THEN
            CALL BMG3_SymStd_residual( 
     &                K, NOG, IFD,
     &                Q(p_U), QF(p_U), SO(p_SO), RES(p_U), 
     &                NLx, NLy, NLz, NStncl,
     &                iWorkMSG, NMSGi, pMSG, 
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            CALL BMG3_SymStd_UTILS_norm_l2( 
     &                RES(p_U), NLx, NLy, NLz, RES_L2, MPICOMM
     &                )
            WRITE(*,110) K, RES_L2
         ENDIF

         ! Relaxation
         DO 20 i=1, IU
            !
            TT1 = MPI_Wtime()
            
            CALL commas_SymStd_relaxdg(  
     &                k, KF, NOG, UPDOWN,
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,
     &                Q, QF, RES, NF, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iGs, jGs, kGs,
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
     &                ) ! ,Qmask)
            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_relax) = 
     &           BMG_rPARMS(id_BMG3_TIME_relax) + TT2 - TT1

            !
            IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) .AND. MyProc.EQ.1 ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
            !
 20      CONTINUE

      ELSE
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN         
            WRITE(*,*) '*** UPDOWN out of range: '
            WRITE(*,*) 'HAVE: UPDOWN = ', UPDOWN
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,19)
         RETURN

      ENDIF

C ==========================================================================

 100  FORMAT (1X,'(3D)LEVEL ',I2,' RESIDUAL NORM = ',1P,E12.5)
 110  FORMAT (1X,'(3D)LEVEL ',I2,
     &           ' AFTER INTERPOLATION RES_L2 = ',1P,E12.5 )

C ===========================================

      RETURN
      END

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE commas_SymStd_relaxdg( 
     &                KG, KF, NOG, UPDOWN,
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,
     &                Q, QF, RES, NF, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iGs, jGs, kGs,
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr, 
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
     &                ) ! ,Qmask)

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_relax performs one relaxation step (either colored
C     point Gauss Seidel, or alternating plane red-black Gauss-Seidel)
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C ------------------------------------------------
C     Includes
C
#include    "BMG_constants_f77.h"
      INCLUDE 'BMG_workspace_f77.h'
      INCLUDE 'BMG_parameters_f77.h'

C ------------------------------------------------
C     Argument Declarations
C 
      INTEGER  NCI, NSO, NSOR, NSORv, NStncl, kg, kf, NBMG_iWORK_PL,
     &         NBMG_rWORK_PL, NOG, NF, NOGm, NMSGi, NMSGr

      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_iWORK_PL(NBMG_iWORK_PL),
     &         IFD, IGRD(NOGm,24), IRELAX_SYM, UPDOWN, IRELAX, ISEP
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS), BMG_rWORK_PL(NBMG_rWORK_PL),
     &         CI(NCI), q(NF), qf(NF), RES(NF),
     &         SO(NSO), SOR(NSOR), MSG_Buffer(NMSGr)

      REAL*8  RES_L2

!      REAL*RKIND :: Qmask(NF)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  iGs, jGs, kGs,
     &         MPICOMM, iWorkMSG(NMSGi), pMSG(NBMG_pMSG,NOGm)
      INTEGER  XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, MPI_IERROR,
     &         NProc, NProcI, NProcJ, NProcK


C ----------------------------
C     Local Declarations
C
      INTEGER  NLx, NLy, NLz, NGx, NGy, NGz, 
     &         NLx_f, NLy_f, NLz_f, NGx_f, NGy_f, NGz_f, 
     &         p_CI, p_SO, p_SOR, p_U, p_ProcGrid
      
C ==========================================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      !
      ! MPI Communicators
      !
      MPICOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM)
      !
      ! MPI Communicatorsn for plane relaxation
      !
      XYPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xy)
      XZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xz)
      YZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_yz)
 
      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         NSORv = 2
      ELSE 
         NSORv = 2  ! the same for now but this may change
      ENDIF

      !
      ! Collect fine-grid dimensions 
      ! ( need to check: KF or NOG )
      !
      CALL BMG3_SymStd_GET_pointers( 
     &            KF, IGRD, NOGm,
     &            p_U, p_SO, p_SOR, p_CI, 
     &            NLx_f, NLy_f, NLz_f, NGx_f, NGy_f, NGz_f
     &            ) 

      CALL BMG3_SymStd_GET_pointers(  
     &            KG, IGRD, NOGm,
     &            p_U, p_SO, p_SOR, p_CI, 
     &            NLx, NLy, NLz, NGx, NGy, NGz
     &            )
      

      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         !
         ! Gauss-Seidel relaxation
         !
         CALL commas_SymStd_relax_GSdg(  
     &             KG, SO(p_SO), QF(p_U), Q(p_U), RES(p_U), SOR(p_SOR),
     &             NLx, NLy, NLz, NGx, NGy, NGz, 
     &             RES_L2, BMG_IOFLAG, NOG, NOGm,
     &             IFD, NStncl, NSORv, IRELAX, IRELAX_SYM, UPDOWN,
     &             iGs, jGs, kGs,
     &             MyProcI, MyProcJ, MyProcK, MyProc,
     &             BMG_MSG_iGRID(p_ProcGrid),
     &             NProcI, NProcJ, NProcK, NProc,
     &             iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr, MPICOMM
     &            ,ISEP ) ! ,Qmask)
         !
      ELSEIF ( IRELAX.EQ.BMG_GS_RB_planes_xy_yz_xz ) THEN

         WRITE(*,*) 'Fatal Error: To busy to code planes !!!'
         CALL MPI_FINALIZE(MPI_IERROR)
         STOP

         !
      ENDIF

C ==========================================================================

      RETURN
      END


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      SUBROUTINE commas_SymStd_relax_GSdg( 
     &                KG, SO, QF, Q, RES, SOR,
     &                NLx, NLy, NLz, NGx, NGy, NGz, 
     &                RES_L2, BMG_IOFLAG, NOG, NOGm,
     &                IFD, NStncl, NSORv, IRELAX, IRELAX_SYM, UPDOWN,
     &                iGs, jGs, kGs,
     &                MyProcI, MyProcJ, MyProcK, MyProc,
     &                ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &                iWork, NMSGi, pMSG, MSG_Buffer, NMSGr, MPICOMM
     &               ,ISEP ) ! ,Qmask)

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_relax_GS performs one sweep of Gauss Seidel (with the
C     correct ordering depending on whether we are on the way down, or
C     up in the symmetric cycling case)
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      USE CPUTIME_MODULE

      IMPLICIT NONE

C ------------------------------------------------
C     Includes
C
      INCLUDE 'mpif.h'
      INCLUDE 'MSG.h'

#include    "BMG_constants_f77.h"
      INCLUDE 'BMG_stencils_f77.h'
      INCLUDE 'BMG_parameters_f77.h'
      INCLUDE 'BMG_workspace_f77.h'

C ------------------------------------------------
C     Argument Declarations
C
      INTEGER NGx, NGy, NGz, NLx, NLy, NLz

      INTEGER IFD, IRELAX, IRELAX_SYM, KG, NOG, NOGm,
     &        NSORv, NStncl, UPDOWN, NMSGi, NMSGr

      INTEGER iGs, jGs, kGs
      INTEGER MyProc, MyProcI, MyProcJ, MyProcK,
     &        NProc, NProcI, NProcJ, NProcK
      INTEGER ProcGrid(NProcI,NProcJ,NProcK)

      INTEGER MPICOMM, iWork(NMSGi), pMSG(NBMG_pMSG,NOGm)
      REAL*RKIND  Q(NLx,NLy,NLz), QF(NLx,NLy,NLz), 
     &        RES(NLx,NLy,NLz),
     &        SO(NLx+1,NLy+1,NLz+1,NStncl), SOR(NLx,NLy,NLz,NSORv)
      REAL*RKIND  MSG_Buffer(NMSGr)
      
      REAL*8 RES_L2

!      REAL*RKIND :: Qmask(NLx,Nly,Nlz)

      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

C ------------------------------------------------
C     Local Declarations
C
      INTEGER i, i1, j ,j1, k, k1
      
      INTEGER ie, iw, js, jn, kbot, kt

      INTEGER pts, ibeg, iend, ptsred, ptsblk
      INTEGER ptstart, ptend, ptstride
      INTEGER ptrn, ierror
      INTEGER relaxtest
      LOGICAL comhalf
      
      REAL*RKIND SOx(NLx+1),SOy(NLy+1),SOz(NLz+1)
      
      INTEGER ISEP

C ==========================================================================

      j1 = NLy-1
      i1 = NLx-1
      k1 = NLz-1
      
c      write(0,*) 'mproc, mpi,mpj,mpk, igs,jgs,kgs = ',myproc,myproci,
c     :   myprocj,myprock,igs,jgs,kgs
c      write(0,*) 'mproc, isep, BMG_SEPARABLE = ',myproc,isep,
c     :   BMG_SEPARABLE, BMG_NONSEPARABLE
      

c      write(0,*) 'RELAX_GS: kg,nog,ifd = ',KG,NOG,ifd

      comhalf = .true.

      CALL cld_cpu('BOXMG-GS')

      
      IF ( KG.LT.NOG .OR. ifd.NE.1 ) THEN 

         !
         !   27-point relaxation (8-color)
         !

C         IF ((UPDOWN.eq.BMG_UP).or.(IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) 
C     &              THEN
C            ptstart = 1
C            ptend   = 8
C            ptstride = 1
C         ELSE
C            ptstart = 8
C            ptend = 1
C            ptstride = -1
C         ENDIF         
C        
C         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<
C
C            DO k=2+mod(mod((pts-1)/4,2)+mod(kGs+1,2),2),K1,2
C               !
C               DO j=2+mod(mod((pts-1)/2,2)+mod(jGs+1,2),2),J1,2
C                  !
C                  DO i=2+mod(mod(pts-1,2)+mod(iGs+1,2),2),I1,2
                     !

         IF ( (UPDOWN.eq.BMG_UP)
     &        .OR. (IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF         
        
        
        relaxtest = 0
        
        IF ( relaxtest .eq. 0 ) THEN ! default
        
         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

c            DO k=2+mod(mod((pts-1)/4,2)+mod(kGs+1,2),2),K1,2
c               !
c               DO j=2+mod(mod((pts-1)/2,2)+mod(jGs+1,2),2),J1,2
c                  !
c                  DO i=2+mod(mod(pts-1,2)+mod(iGs+1,2),2),I1,2
             DO k = 2,k1
              
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                     !
                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !
            
            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
            !
            
            ENDIF
          ENDDO

        ELSEIF ( relaxtest .eq. 1 ) THEN ! non-redblack GS

             DO k = 2,k1
              
               DO j=2, j1
                  !
                  !
                  DO i = 2,i1
                     !
                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO

        ELSEIF ( relaxtest .eq. 2 ) THEN ! layered redblack GS

c         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

c do 'red' points on first plane
            ptsred = ptstart
            ptsblk = ptend
            DO  k=2, 2
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsred,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO

c do 'red' points on k plane and 'black' points on k-1 plane

            DO  k=3, k1

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsred,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO

            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
            !
            
            ENDIF

c black points
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-2)+ptsblk,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k-1) = ( QF(i,j,k-1) 
     &                    + SO(i,j,k-1,kpw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k-1,kpnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k-1,kps)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k-1,kpsw)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k-1,kpw)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k-1,kpnw)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k-1,kps)*Q(i,j-1,k-1)
     &                    + SO(i,j,k-1,kpsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k-1,kb)*Q(i,j,k-1-1)
     &                    + SO(i,j,k-1,kbw)*Q(i-1,j,k-1-1)
     &                    + SO(i,j+1,k-1,kbnw)*Q(i-1,j+1,k-1-1)
     &                    + SO(i,j+1,k-1,kbn)*Q(i,j+1,k-1-1)
     &                    + SO(i+1,j+1,k-1,kbne)*Q(i+1,j+1,k-1-1)
     &                    + SO(i+1,j,k-1,kbe)*Q(i+1,j,k-1-1)
     &                    + SO(i+1,j,k-1,kbse)*Q(i+1,j-1,k-1-1)
     &                    + SO(i,j,k-1,kbs)*Q(i,j-1,k-1-1)
     &                    + SO(i,j,k-1,kbsw)*Q(i-1,j-1,k-1-1)
     &                    + SO(i,j,k-1+1,kb)*Q(i,j,k-1+1)
     &                    + SO(i,j,k-1+1,kbe)*Q(i-1,j,k-1+1)
     &                    + SO(i,j+1,k-1+1,kbse)*Q(i-1,j+1,k-1+1)
     &                    + SO(i,j+1,k-1+1,kbs)*Q(i,j+1,k-1+1)
     &                    + SO(i+1,j+1,k-1+1,kbsw)*Q(i+1,j+1,k-1+1)
     &                    + SO(i+1,j,k-1+1,kbw)*Q(i+1,j,k-1+1)
     &                    + SO(i+1,j,k-1+1,kbnw)*Q(i+1,j-1,k-1+1)
     &                    + SO(i,j,k-1+1,kbn)*Q(i,j-1,k-1+1)
     &                    + SO(i,j,k-1+1,kbne)*Q(i-1,j-1,k-1+1)
     &                    )*SOR(i,j,k-1,msor)
                     !
                  ENDDO
                  !
               ENDDO

               !
            ENDDO

c do 'black' points on last plane
            DO  k=k1, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsblk,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            
       ENDIF ! relaxtest
       !
c         ENDDO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<

        IF ( relaxtest .ne. 0 ) THEN
 
            ptrn = 1
            !
            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
     
        ENDIF
            !
!         ENDDO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<


         !
      ELSE

      CALL cld_cpu('BOXMG-GSfine')
         !  
         !  7 point relaxation (2-color)
         !         
         IF ( (UPDOWN.eq.BMG_UP)
     &        .OR. (IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF         

      IF ( ISEP .eq. BMG_SEPARABLE ) THEN

      
      DO  k=1, k1+1
        SOz(k) = SO(3,3,k,kb)
        SOx(k) = SO(3,3,k,kpw)
        SOy(k) = SO(3,3,k,kps)
c       IF ( myproc == 1 ) write(0,*) 'k, SOz,x,y = ',k,SOz(k),
c     &  SO(3,3,k,kpw),SO(3,3,k,kps)
c        IF ( myproc == 1 ) write(0,*) 'j, SOy = ',
      ENDDO
c        SOz(2) = SOz(3)
c        SOz(k1+1) = SOz(k1)
        
      
      DO j=1, j1+1
c        SOy(j) = SO(3,j,3,kps)
c        IF ( myproc == 1 ) write(0,*) 'j, SOy,x,z = ',j,SOy(j),
c     &   SO(3,j,3,kpw),SO(3,j,3,kb)
c       IF ( myproc == 1 ) write(0,*) 'k, SOz = ',j,SO(3,j,3,kb)
c        IF ( myproc == 1 ) write(0,*) 'i, SOx = ',j,SO(3,j,3,kpw)
c        write(0,*) 'j, SOy, myproc = ',j,SOy(j),myproc
      ENDDO
c        SOy(2) = SOy(3)
c        SOy(j1+1) = SOy(j1)
      
      
      DO i=1,i1+1
c        SOx(i) = SO(i,3,3,kpw)
c        IF ( myproc == 1 ) write(0,*) 'i, SOx,y,z = ',i,SOx(i),
c     &   SO(i,3,3,kps),SO(i,3,3,kb)
      ENDDO
c        SOx(2) = SOx(3)
c        SOx(i1+1) = SOx(i1)
      
      ENDIF

         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

          IF ( ISEP .eq. BMG_SEPARABLE ) THEN

          IF ( .true. ) THEN
          
          
            DO  k=2, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SOx(k)*Q(i-1,j,k)
     &                        + SOy(k)*Q(i,j+1,k)
     &                        + SOx(k)*Q(i+1,j,k)
     &                        + SOy(k)*Q(i,j-1,k)
     &                        + SOz(k)*Q(i,j,k-1)
     &                        + SOz(k+1)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)
                  END DO
                  !
               END DO
               !
            END DO
            
          ELSE ! old try

          kbot = 0
          kt = 0
          js = 0
          jn = 0
          iw = 0
          ie = 0
          
          IF (  MyProcK .eq. 1 )     kbot = 1
          IF (  MyProcK .eq. NProcK) kt = 1
          IF (  MyProcJ .eq. 1)      js = 1
          IF (  MyProcJ .eq. NProcJ) jn = 1
          IF (  MyProcI .eq. 1)      iw = 1
          IF (  MyProcI .eq. NProcI) ie = 1
          
            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  
                  IF ( ibeg .eq. 2 ) ibeg = ibeg + 2*iw
                  IF ( iend .eq. i1 ) iend = iend - 2*ie
                  
                  !
                  DO i=ibeg,iend,2
                  !
c                  IF ( .not. (( (k.eq.2      .and. MyProcK.eq.1 ) .or.
c     &                        ( k .eq. k1    .and. MyProcK .eq. NProcK) .or.
c     &                        ( j .eq. 2     .and. MyProcJ .eq. 1)  .or.
c     &                        ( j .eq. j1    .and. MyProcJ .eq. NProcJ) .or.
c     &                        ( i .eq. ibeg  .and. MyProcI .eq. 1)  .or. 
c     &                        ( i .eq. iend  .and. MyProcI .eq. NProcI) ) ) THEN
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SOx(i)*Q(i-1,j,k)
     &                        + SOy(j+1)*Q(i,j+1,k)
     &                        + SOx(i+1)*Q(i+1,j,k)
     &                        + SOy(j)*Q(i,j-1,k)
     &                        + SOz(k)*Q(i,j,k-1)
     &                        + SOz(k+1)*Q(i,j,k+1)
     &                        )/(SOx(i) + SOx(i+1) + SOy(j) + 
     &                           SOy(j+1) + SOz(k) + SOz(k+1))

                  END DO
                  !
               END DO
               !
            END DO

c edges:

           IF ( kbot .eq. 1 ) THEN
            k=2

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg 
                  !
                  DO i=ibeg,iend,2

                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)


                 END DO
              END DO
           
           ENDIF

           IF ( kt .eq. 1 ) THEN
            k=k1

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg 
                  !
                  DO i=ibeg,iend,2

                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)


                 END DO
              END DO
           
           ENDIF
           
           IF ( js .eq. 1 ) THEN
                      
            DO  k=2+kbot, k1-kt
              j = 2
              ibeg = mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend = 2*((i1-ibeg)/2)+ibeg 

               DO i=ibeg,iend,2
    
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)

                END DO
             END DO
          ENDIF

           
           IF ( jn .eq. 1 ) THEN
                      
            DO  k=2+kbot, k1-kt
              j = j1
              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg 
                  DO i=ibeg,iend,2

    
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)

                END DO
             END DO
          ENDIF

           
           IF ( ie .eq. 1 ) THEN
            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn

              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg 
              i = ibeg
               IF ( ibeg .eq. 2 ) THEN
    
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)

                ENDIF
                END DO
             END DO
          ENDIF

           IF ( iw .eq. 1 ) THEN
            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn

              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg 
              
              i = iend
              
               IF ( i .eq. i1 ) THEN
    
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)

                ENDIF
                END DO
             END DO
          ENDIF
          
         ENDIF
         

         ELSE  ! Non separable

            DO  k=2, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !
                  
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)
                  END DO
                  !
               END DO
               !
            END DO

          ENDIF
            !

            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
            !
            
            ENDIF
            
            
            !
         END DO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<
            !
            IF ( .not. comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
            !
            
            ENDIF
            !
C         END DO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<

      CALL cld_cpu('BOXMG-GSfine')

      ENDIF

      ptrn = 1

      CALL MSG_tbdx_close(Q, MSG_Buffer,
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &     iWork(pMSG(ipL_MSG_Proc,KG)),
     &     iWork(pMSG(ipL_MSG_Ipr,KG)),
     &     iWork(pMSG(ipL_MSG_Index,KG)),
     &     ptrn, ierror)

      IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) ) THEN
         !
         CALL BMG3_SymStd_residual( 
     &             KG, NOG, IFD,
     &             Q, QF, SO, RES, NLx, NLy, NLz, NStncl, 
     &             iWork, NMSGi, pMSG, 
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )
         CALL BMG3_SymStd_UTILS_norm_l2( 
     &             RES, NLx, NLy, NLz, RES_L2, MPICOMM
     &             )
         !


      ENDIF

      CALL cld_cpu('BOXMG-GS')

C ==========================================================================

      RETURN
      END






#endif
