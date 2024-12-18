C---------------------------------------------------------------
C
C       3D ARRAY DRIVER 
C       SIXTH-ORDER LOW-PASS IMPLICIT TANGENT FILTER
C       (RAYMOND, MWR, 116, 2132-2141)
C

       SUBROUTINE RAYMOND3D_LOWPASS(XY3D,RESULT3D,EPS,NX,NY,NZ)

C
C       MODIFIED BY LJW May 2011
C---------------------------------------------------------------

        IMPLICIT NONE

C Passed variables

        INTEGER, INTENT(IN)  :: NX, NY, NZ
        REAL,    INTENT(IN)  :: XY3D(NX,NY,NZ)
        REAL,    INTENT(OUT) :: RESULT3D(NX,NY,NZ)
        REAL,    INTENT(IN)  :: EPS

C Local variables

        INTEGER :: I, J, K, NM
        REAL, ALLOCATABLE :: XY(:), RESULT(:)

        NM = MAX(NX, NY, NZ)

        ALLOCATE(XY(NM))
        ALLOCATE(RESULT(NM))

        RESULT3D(:,:,:) = XY3D(:,:,:)

C X-PASS
        IF( NX > 1 ) THEN
          DO K = 1,NZ
          DO J = 1,NY
          DO I = 1,NX
            XY(I) = RESULT3D(I,J,K)
          ENDDO
          CALL RAYMOND1D_LOWPASS(XY,RESULT,EPS,NX)
          DO I = 1,NX  
            RESULT3D(I,J,K) = RESULT(I)
          ENDDO
          ENDDO
          ENDDO
        ENDIF

C Y-PASS
        IF( NY > 1 ) THEN
          DO K = 1,NZ
          DO I = 1,NX
          DO J = 1,NY  
            XY(J) = RESULT3D(I,J,K)
          ENDDO
          CALL RAYMOND1D_LOWPASS(XY,RESULT,EPS,NY)
          DO J = 1,NY  
            RESULT3D(I,J,K) = RESULT(J)
          ENDDO
          ENDDO
          ENDDO
        ENDIF

C Z-PASS
        IF( NZ > 1 ) THEN 
          DO J = 1,NY
          DO I = 1,NX
          DO K = 1,NZ  
            XY(K) = RESULT3D(I,J,K)
          ENDDO
          CALL RAYMOND1D_LOWPASS(XY,RESULT,EPS,NZ)
          DO K = 1,NZ
            RESULT3D(I,J,K) = RESULT(K)
          ENDDO
          ENDDO
          ENDDO
        ENDIF

        DEALLOCATE(XY)
        DEALLOCATE(RESULT)

        RETURN
        END
C---------------------------------------------------------------
C
C       3D ARRAY DRIVER WHICH USES DIFFERENT COEFF FOR
C       HORIZONTAL AND VERTICAL
C       SIXTH-ORDER LOW-PASS IMPLICIT TANGENT FILTER
C       (RAYMOND, MWR, 116, 2132-2141)
C

       SUBROUTINE RAYMOND3DG_LOWPASS(XY3D,RESULT3D,EPS,NX,NY,NZ)

C
C       MODIFIED BY LJW May 2011
C---------------------------------------------------------------

        IMPLICIT NONE

C Passed variables

        INTEGER, INTENT(IN)  :: NX, NY, NZ
        REAL,    INTENT(IN)  :: XY3D(NX,NY,NZ)
        REAL,    INTENT(OUT) :: RESULT3D(NX,NY,NZ)
        REAL,    INTENT(IN)  :: EPS(2)

C Local variables

        INTEGER :: I, J, K, NM
        REAL, ALLOCATABLE :: XY(:), RESULT(:)

        NM = MAX(NX, NY, NZ)

        ALLOCATE(XY(NM))
        ALLOCATE(RESULT(NM))

        RESULT3D(:,:,:) = XY3D(:,:,:)

C X-PASS
        IF( NX > 1 .AND. EPS(1) > 0.0 ) THEN
          DO K = 1,NZ
          DO J = 1,NY
          DO I = 1,NX
            XY(I) = RESULT3D(I,J,K)
          ENDDO
          CALL RAYMOND1D_LOWPASS(XY,RESULT,EPS(1),NX)
          DO I = 1,NX  
            RESULT3D(I,J,K) = RESULT(I)
          ENDDO
          ENDDO
          ENDDO
        ENDIF

C Y-PASS
        IF( NY > 1 .AND. EPS(1) > 0.0 ) THEN
          DO K = 1,NZ
          DO I = 1,NX
          DO J = 1,NY  
            XY(J) = RESULT3D(I,J,K)
          ENDDO
          CALL RAYMOND1D_LOWPASS(XY,RESULT,EPS(1),NY)
          DO J = 1,NY  
            RESULT3D(I,J,K) = RESULT(J)
          ENDDO
          ENDDO
          ENDDO
        ENDIF

C Z-PASS
        IF( NZ > 1 .AND. EPS(2) > 0.0) THEN 
          DO J = 1,NY
          DO I = 1,NX
          DO K = 1,NZ  
            XY(K) = RESULT3D(I,J,K)
          ENDDO
          CALL RAYMOND1D_LOWPASS(XY,RESULT,EPS(2),NZ)
          DO K = 1,NZ
            RESULT3D(I,J,K) = RESULT(K)
          ENDDO
          ENDDO
          ENDDO
        ENDIF

        DEALLOCATE(XY)
        DEALLOCATE(RESULT)

        RETURN
        END
C---------------------------------------------------------------
C
C       2D ARRAY DRIVER 
C       SIXTH-ORDER LOW-PASS IMPLICIT TANGENT FILTER
C       (RAYMOND, MWR, 116, 2132-2141)
C

       SUBROUTINE RAYMOND2D_LOWPASS(XY2D,RESULT2D,EPS,NX,NY)

C
C       MODIFIED BY LJW May 2011
C---------------------------------------------------------------

        IMPLICIT NONE

C Passed variables

        INTEGER, INTENT(IN)  :: NX, NY
        REAL,    INTENT(IN)  :: XY2D(NX,NY)
        REAL,    INTENT(OUT) :: RESULT2D(NX,NY)
        REAL,    INTENT(IN)  :: EPS

C Local variables

        INTEGER :: I, J, NM
        REAL, ALLOCATABLE :: XY(:), RESULT(:)
C       REAL :: XY(513), RESULT(513)

        NM = MAX(NX, NY)

        ALLOCATE(XY(NM))
        ALLOCATE(RESULT(NM))

        RESULT2D(:,:) = XY2D(:,:)

C X-PASS
        IF( NX > 1 ) THEN
          DO J = 1,NY
          DO I = 1,NX
            XY(I) = RESULT2D(I,J)
          ENDDO
          CALL RAYMOND1D_LOWPASS(XY,RESULT,EPS,NX)
          DO I = 1,NX
            RESULT2D(I,J) = RESULT(I)
          ENDDO
          ENDDO
        ENDIF

C Y-PASS
        IF( NY > 1 ) THEN
          DO I = 1,NX
          DO J = 1,NY
            XY(J) = RESULT2D(I,J)
          ENDDO
          CALL RAYMOND1D_LOWPASS(XY,RESULT,EPS,NY)
          DO J = 1,NY
            RESULT2D(I,J) = RESULT(J)
          ENDDO
          ENDDO
        ENDIF

        DEALLOCATE(XY)
        DEALLOCATE(RESULT)

        RETURN
        END

C---------------------------------------------------------------
C

       SUBROUTINE RAYMOND1D_LOWPASS(XY,RESULT,EPS,N)
C
C       1D SIXTH-ORDER LOW-PASS IMPLICIT TANGENT FILTER
C       (RAYMOND, MWR, 116, 2132-2141)
C
C*************************************************************
C***     THIS CODE IS COPIED FROM A LISTING PROVIDED       ***
C***     BY WILLIAM H RAYMOND. SOME NOTATIONAL CHANGES     ***
C***     HAVE BEEN MADE IN THE ROUTINE LOWPAS. THE         ***
C***     ROUTINE INVLOW HAS BEEN COPIED ALMOST VERBATIM.   ***
C*************************************************************
C
C       XY     UNFILTERED VALUES ON INPUT
C       RESULT FILTERED VALUES ON OUTPUT.
C       N      NUMBER OF VALUES.
C       EPS    FILTER PARAMETER
C              (DETERMINES CUTOFF)
C
C       MODIFIED BY LJW May 2011
C---------------------------------------------------------------
C
        IMPLICIT NONE

C Passed variables
        INTEGER, INTENT(IN) :: N
        REAL, INTENT(IN)  :: XY(N)
        REAL, INTENT(OUT) :: RESULT(N)
        REAL, INTENT(IN)  :: EPS

C Local Variables
        INTEGER :: NM1, NM2, NM3, NM4, J
        REAL    :: RHS(N), XDASH(N)

        NM1 = N-1
        NM2 = N-2
        NM3 = N-3
        NM4 = N-4
C
        RHS(1) = 0
        RHS(N) = 0
        RHS(  2) = EPS*(XY(  1)-2.0*XY(  2)+XY(  3))
        RHS(NM1) = EPS*(XY(NM2)-2.0*XY(NM1)+XY(  N))
        RHS(  3) = EPS*(-1.0*(XY(  1)+XY(  5))
     +                  +4.0*(XY(  2)+XY(  4))
     +                  -6.0* XY(  3)         )
        RHS(NM2) = EPS*(-1.0*(XY(  N)+XY(NM4))
     +                  +4.0*(XY(NM1)+XY(NM3))
     +                  -6.0* XY(NM2)         )
        DO 1000 J=4,NM3
        RHS(J) = EPS*(       (XY(J-3)+XY(J+3))
     +                - 6.0*(XY(J-2)+XY(J+2))
     +                +15.0*(XY(J-1)+XY(J+1))
     +                -20.0* XY(  J)         )
 1000   CONTINUE
C
C       SOLVE THE LINEAR SYSTEM FOR XDASH
C
        CALL RAYMOND_INVLOW(RHS,N,XDASH,EPS)
C
C       ADD CORRECTION TO GET FILTERED VALUES.
C
        DO 2000 J=1,N
        RESULT(J) = XY(J) + XDASH(J)
 2000   CONTINUE
C
        RETURN
        END

C---------------------------------------------------------------
C
C
        SUBROUTINE RAYMOND_INVLOW(BB,N,XANS,EP)
C
C       GAUSSIAN ELIMINATION FOR LOW-PASS FILTER.
C
C       SIXTH-ORDER LOW-PASS IMPLICIT TANGENT FILTER.
C       (REF: WILLIAM H RAYMOND, MWR, 116, 2132-2124)
C
        IMPLICIT NONE
        INTEGER I,N

C
        REAL      A(N),B(N),C(N),D(N),E(N),
     +            DELTA(N),BETA(N),W(N),GAM(N),
     +            H(N),XANS(N),BB(N),PI(N),
     +            AP(N),F(N),Z(N)

        REAL EP
C
C---------------------------------------------------------------
C       INITIALIZE THE MATRIX
C
        DO 10 I=4,N-3
        Z(I) = 1.0-EP
        A(I) = 6.0*(1.0+EP)
        B(I) = 15.0*(1.0-EP)
        C(I) = 20.0*(1.0+EP)
        D(I) = B(I)
        E(I) = A(I)
        F(I) = Z(I)
   10   CONTINUE
C
        Z(1) = 0.
        Z(2) = 0.
        Z(3) = 0.
C
        A(1) = 0.
        A(2) = 0.
        A(3) = 1.0+EP
C
        B(1) = 0.
        B(2) = 1.0-EP
        B(3) = 4.0*(1.0-EP)
C
        C(1) = 1.0
        C(2) = 2.0*(1.0+EP)
        C(3) = 6.0*(1.0+EP)
C
        D(1) = 0.
        D(2) = 1.0-EP
        D(3) = 4.0*(1.0-EP)
C
        E(1) = 0.
        E(2) = 0.
        E(3) = 1.0+EP
C
        F(1) = 0.
        F(2) = 0.
        F(3) = 0.
C
C
        Z(N-2) = 0.
        Z(N-1) = 0.
        Z(N) = 0.
C
        A(N-2) = 1.0+EP
        A(N-1) = 0.
        A(N) = 0.
C
        B(N-2) = 4.0*(1.0-EP)
        B(N-1) = 1.0-EP
        B(N) = 0.
C
        C(N-2) = 6.0*(1.0+EP)
        C(N-1) = 2.0*(1.0+EP)
        C(N) = 1.0
C
        D(N-2) = 4.0*(1.0-EP)
        D(N-1) = 1.0-EP
        D(N) = 0.
C
        E(N-2) = 1.0+EP
        E(N-1) = 0.
        E(N) = 0.
C
        F(N-2) = 0.
        F(N-1) = 0.
        F(N) = 0.
C
C       Step One.
        BETA(1) = 0.0
        DELTA(2) = B(2)
        W(1) = C(1)
        PI(1) = 0.
        AP(1) = 0.
        AP(2) = 0.
        AP(3) = A(3)
        W(2) = C(2)-DELTA(2)*BETA(1)
        GAM(1) = 0.0
        BETA(2) = (D(2)-DELTA(2)*GAM(1))/W(2)
        GAM(2) = (E(2)-PI(1)*DELTA(2))/W(2)
        PI(2) = 0.0
        DELTA(3) = (B(3)-AP(3)*BETA(1))
        W(3) = C(3)-DELTA(3)*BETA(2)-AP(3)*GAM(1)
        BETA(3) = (D(3)-AP(3)*PI(1)-DELTA(3)*GAM(2))/W(3)
        GAM(3) = (E(3)-DELTA(3)*PI(2))/W(3)
        PI(3) = 0.0
C
C       Step Two
        DO 20 I=4,N
        AP(I) = A(I)-Z(I)*BETA(I-3)
        DELTA(I) = B(I)-AP(I)*BETA(I-2)-Z(I)*GAM(I-3)
        W(I) = C(I)-AP(I)*GAM(I-2)-DELTA(I)*BETA(I-1)
     +             -Z(I)*PI(I-3)
        BETA(I) = (D(I)-AP(I)*PI(I-2)-DELTA(I)*GAM(I-1))/W(I)
        GAM(I) = (E(I)-DELTA(I)*PI(I-1))/W(I)
        PI(I) = F(I)/W(I)
   20   CONTINUE
C
C       Step Three
        H(1) = BB(1)/W(1)
        H(2) = (BB(2)-DELTA(2)*H(1))/W(2)
        H(3) = (BB(3)-DELTA(3)*H(2)-AP(3)*H(1))/W(3)
        DO 30 I=4,N
        H(I) = (BB(I)-DELTA(I)*H(I-1)-AP(I)*H(I-2)
     +               -Z(I)*H(I-3))/W(I)
   30   CONTINUE
C
C       Step Four
        XANS(N) = H(N)
        XANS(N-1) = H(N-1)-BETA(N-1)*XANS(N)
        XANS(N-2) = H(N-2)-BETA(N-2)*XANS(N-1)-GAM(N-2)*XANS(N)
        DO 40 I=N-3,1,-1
        XANS(I) = H(I)-BETA(I)*XANS(I+1)-GAM(I)*XANS(I+2)
     +                                  -PI(I)*XANS(I+3)
   40   CONTINUE
        RETURN
        END
