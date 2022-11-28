      SUBROUTINE INVERT(N,A,RST,DET)
C
C     Subroutine to invert a matrix, and compute the determinant.
C     John Tonry, 9/2/80.
C
C     INVERT's arguments:
C
C     N - The dimension of the matrix
C     A - The NxN matrix to be inverted. Upon successful inversion, A
C         contains the inverse. A must be real*8
C     RST - A Nx1 scratch vector (the row status vector)
C     DET - The determinant of the matrix. DET is set to 0 for a
C         singular matrix, and in that case, A contains garbage.
C
      REAL*8 A(N,N), SAVE, PIVOT, ONROW, CPREV, CNOW
      INTEGER*2 RST(2,N)
C
      MRANK = 0
      ISIGN = 1
      DET = 0.
      DO 10 J = 1,N
         DO 11 I = 1,2
            RST(I,J) = 0
 11      CONTINUE
 10   CONTINUE
C
C     Loop over columns, reducing each
C
      DO 500 I = 1,N
C
C     Find the pivot element
C
         PIVOT = 0
         NROW = 0
         NCOL = 0
         DO 30 J = 1,N
            IF(RST(1,J).NE.0) GO TO 30
            DO 20 K = 1,N
               IF(RST(1,K).NE.0) GO TO 20
               IF(PIVOT.GE.DABS(A(J,K))) GO TO 20
               PIVOT = DABS(A(J,K))
               NROW = J
               NCOL = K
 20         CONTINUE
 30      CONTINUE
         PIVOT = A(NROW,NCOL)
         IF(PIVOT.EQ.0) GO TO 300
         RST(1,NCOL) = NROW
         RST(2,NCOL) = I
C
C     Swap pivot element onto the diagonal
C
         DO 40 K = 1,N
            SAVE = A(NROW,K)
            A(NROW,K) = A(NCOL,K)
            A(NCOL,K) = SAVE
 40      CONTINUE
C
C     Reduce pivot column
C
         DO 50 J = 1,N
            A(J,NCOL) = -A(J,NCOL)/PIVOT
 50      CONTINUE
         A(NCOL,NCOL) = 1/PIVOT
C
C     Reduce other columns
C
         DO 60 K = 1,N
            IF(K.EQ.NCOL) GO TO 60
C
C     Find maximum of column to check for singularity
C
            CPREV = 0
            DO 70 J = 1,N
               CPREV = DMAX1(CPREV,DABS(A(J,K)))
 70         CONTINUE
C
C     Reduce the column
C
            ONROW = A(NCOL,K)
            A(NCOL,K) = 0
            DO 80 J = 1,N
               A(J,K) = A(J,K) + ONROW*A(J,NCOL)
 80         CONTINUE
C
C     Find the new maximum of the column
C
            CNOW = 0
            DO 90 J = 1,N
               CNOW = DMAX1(CNOW,DABS(A(J,K)))
 90         CONTINUE
C
C     Quit if too many figures accuracy were lost (singular)
C
            IF(CNOW.EQ.0) GOTO 300
            DECR = CPREV / CNOW
            IF(DECR.GT.1E8) GO TO 300
C
 60      CONTINUE
C
         DET = DET + DLOG10(DABS(PIVOT))
         IF(PIVOT.LT.0) ISIGN = -ISIGN
         MRANK = MRANK + 1
 500  CONTINUE
C
C     Now untangle the mess
C
      DO 100 J = 1,N
         DO 110 K = 1,N
            IF(RST(2,K).NE.(N + 1 - J)) GO TO 110
            NCOL = RST(1,K)
            IF(NCOL.EQ.K) GO TO 100
            DO 120 L = 1,N
               SAVE = A(L,NCOL)
               A(L,NCOL) = A(L,K)
               A(L,K) = SAVE
 120        CONTINUE
            GO TO 100
 110     CONTINUE
 100  CONTINUE
      IF(ABS(DET).LT.35) DET = ISIGN*10.**DET
      RETURN
C
C     Singular exit
C
300   DET = 0
      RETURN
C
      END
