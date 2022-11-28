      SUBROUTINE MINI(N,A,FUNK,DFUNK,D2FUNK,COV,NDERIV,IERR)
*
*     This is John Tonry's black box minimization (and fitting) program.
*     Revision 2.0, 11/17/82.
*
*     MINI's arguments are as follows:
*     N - The number of parameters to be varied in searching for a
*        minimum.
*     A - A real*8 vector, it serves three purposes:
*        (1) It passes MINI the initial estimate for the parameters
*        (2) It passes arguments to the function to be minimized
*        (3) It returns the minimum parameter values
*     FUNK - The function to be minimized. FUNK should take A as an
*        argument and return the value of the function: FUNK(A).
*        FUNK must return real*8 values.
*     DFUNK - DFUNK should take A as an argument and return the first
*        derivatives of the function at A.
*        DFUNK must return real*8 values: DFUNK(A,DF)
*     D2FUNK - D2FUNK should take A as an argument and return the second
*        derivatives of the function at A: D2FUNK(A,COV)
*        If DFUNK and D2FUNK are not wanted, create dummies or else
*        disregard linker error messages.
*     COV - A NxN real*8 matrix in which the covariance matrix of the
*        fit is returned.
*     NDERIV - Variable to govern how MINI gets its derivatives:
*        NDERIV = 0 for a function of arbitrary A and numerical derivatives.
*        NDERIV = 1 for a function of arbitrary A, analytic first derivatives
*                  provided by DFUNK, and numerical second derivatives.
*        NDERIV = 2 for a function of arbitrary A, analytic first and second
*                  derivatives provided by DFUNK and by D2FUNK.
*        NDERIV = 3 for a function that is quadratic in A, and MINI
*                  will iterate once, computing the minimum exactly.
*     IERR - On input (0/1) governs whether MINI will print each iteration
*        On output, IERR indicates various error conditions:
*        IERR = 0 for no error
*        IERR = 1 for maximum iteration exceeded
*        IERR = 2 for singular matrix
*
*     Descriptions of some of the variables:
*     A - Argument for the function
*     A0 - Current guess for the minimum
*     AI - increments for a0 in computing derivatives
*     DA - Vector from A0 to new guess for the minimum
*     DF - First derivatives of the function
*     D2F - Second derivatives of the function
*     LAMBDA - Governs mix of gradient and analytic searches
*     ITER - Maximum number of iterations
*     QFRAC - Maximum fractional change for successful exit
*
*     The calling program should be as follows (eg):
*********************************************************
*     REAL*8 A(4)
*     REAL*8 COV(4,4)
*     EXTERNAL CHI, DCHI, D2CHI
*     ... (Initialize A to the guess for the minimum)
*     CALL MINI(4,A,CHI,DCHI,D2CHI,COV,NDERIV)
*     ...
*     FUNCTION CHI(A)
*     REAL*8 A(4), CHI
*     ... (define the function)
*     SUBROUTINE DCHI(A,DF)
*     REAL*8 A(4), DF(4)
*     ... (dummy or compute the derivatives)
*     SUBROUTINE D2CHI(A,COV)
*     REAL*8 A(4), COV(4,4)
*     ... (dummy or compute the derivatives)
*************************************************************
*
      REAL*8 A(N),DF(16),A0(16),DA(16),AI(16),RST(16)
      REAL*8 COV(N,N),D2F(16,16)
      REAL*8 FNOW, FTHEN, FMINUS, CURVE, FUNK

      ITEST = 0
*
*     Quit if N > 16, the space allocated
*
      IF(N.GT.16) THEN
          WRITE(6,*) ' Too many parameters, maximum = 16'
          RETURN
      END IF
*
*     Define a few parameters
*
      LAMBDA = -3
      BASE = 10
      ITERMAX = 20
      QFRAC = 1E-6
      DFRAC = .02
      PLACES = 1E-7
      VARY = 1E-5
*
*     If NDERIV = 3, compute the minimum directly and exit.
*
      IF(NDERIV.NE.3) GO TO 8
      DO 900 I = 1,N
          A(I) = 0
900   CONTINUE
      FNOW = FUNK(A)
      DO 901 I = 1,N
          A(I) = 1
          DF(I) = FUNK(A)
          DO 902 J = 1,I
              A(J) = A(J) + 1
              D2F(I,J) = FUNK(A) - DF(I) - DF(J) + FNOW
              D2F(J,I) = D2F(I,J)
              COV(I,J) = D2F(I,J)
              COV(J,I) = D2F(J,I)
              A(J) = 0
902       CONTINUE
          A(I) = 0
901   CONTINUE
      DO 903 I = 1,N
          DF(I) = DF(I) - .5 * D2F(I,I) - FNOW
903   CONTINUE
      CALL INVERT(N,COV,RST,DET)
      IF(DET.EQ.0) THEN
          IERR = 2
          RETURN
      END IF
      DO 904 I = 1,N
          A(I) = 0
          DO 905 J = 1,N
              A(I) = A(I) - COV(I,J) * DF(J)
905       CONTINUE
904   CONTINUE
      FNOW = FUNK(A)
      IF(IERR.EQ.1) CALL VPRINT(N,A,FNOW,0,LAMBDA)
      GO TO 150
8     CONTINUE
*
*     Initialize A0
*
      DO 906 I = 1,N
          A0(I) = A(I)
906   CONTINUE
      FNOW = FUNK(A)
      IF(IERR.EQ.1) CALL VPRINT(N,A0,FNOW,0,LAMBDA)
*
*     Initialize AI
*
      DO 907 I = 1,N
          AI(I) = DABS(VARY*A0(I))
          IF(AI(I).EQ.0) AI(I) = 1E-6
907   CONTINUE
*
*     Begin iteration to find minimum
*
      DO 500 ITER = 1,ITERMAX
      FTHEN = FNOW

* Compute the function derivatives.
      DO 908 J = 1,N
          A(J) = A0(J)
908   CONTINUE

      IF(NDERIV.GT.0) GOTO 10

* First case: NDERIV = 0 so entirely numerical derivatives are required
*  First the 1st derivatives
      IF(ITEST.EQ.1) THEN
         WRITE(6,*) 'Variables: ', (A0(III),III=1,N)
         WRITE(6,*) 'Offsets:   ', (AI(III),III=1,N)
      END IF
      DO 909 J = 1,N
          A(J) = A0(J) + AI(J)
          DF(J) = FUNK(A)
          A(J) = A0(J)
909   CONTINUE
* Next, get the off diagonal 2nd derivatives.
      DO 910 J = 2,N
          DO 911 K = 1,J-1
              A(K) = A0(K) + AI(K)
              A(J) = A0(J) + AI(J)
              D2F(J,K) = (FUNK(A)-DF(J)-DF(K)+FNOW)/(AI(J)*AI(K))
              D2F(K,J) = D2F(J,K)
              A(J) = A0(J)
              A(K) = A0(K)
911       CONTINUE
910   CONTINUE
      IF(ITEST.EQ.1) CALL TESTER(N,DF,D2F)
* Finally do the on diagonal 2nd derivatives, and fix the 1st ones.
      DO 912 J = 1,N
          A(J) = A0(J) - AI(J)
          FMINUS = FUNK(A)
          D2F(J,J) = (FMINUS + DF(J) - 2*FNOW)/(AI(J)*AI(J))
          DF(J) = (DF(J) - FMINUS)/(2*AI(J))
          A(J) = A0(J)
912   CONTINUE
      GOTO 50

* Second case: NDERIV = 1 so analytic first derivatives are available
10    CALL DFUNK(A,DF)
      IF(NDERIV.GT.1) GOTO 20
      DO 913 J = 1,N
          A(J) = A0(J) + AI(J)
          CALL DFUNK(A,DA)
          A(J) = A0(J)
          DO 914 I = 1,J
              D2F(I,J) = (DA(I) - DF(I)) / AI(J)
              D2F(J,I) = D2F(I,J)
914       CONTINUE
913   CONTINUE
      GOTO 50

* Third case: NDERIV = 2 so analytic derivatives are available
20    CALL D2FUNK(A,COV)
      DO 915 J = 1,N
          DO 916 I = 1,N
              D2F(I,J) = COV(I,J)
916       CONTINUE
915   CONTINUE

50    CONTINUE
      IF(ITEST.EQ.1) CALL TESTER(N,DF,D2F)
* Compute better estimates for the increments.
      DO 917 J = 1,N
          CURVE = D2F(J,J)
          IF(CURVE.EQ.0) CURVE = 1E-5
          AI(J) = DSQRT((DF(J)*DFRAC/CURVE)**2+DABS(FNOW*PLACES/CURVE))
917   CONTINUE
*
*     Begin loop to find a direction along which function decreases
*
      DO 300 II = 1,15
*
*     Get weight matrix
*
      DO 918 J = 1,N
          DO 919 K = 1,J-1
              COV(J,K) = D2F(J,K)
              COV(K,J) = D2F(J,K)
919       CONTINUE
          COV(J,J) = DABS(D2F(J,J)*(1 + BASE**LAMBDA))
918   CONTINUE

      if(itest.eq.1) then
         write(6,*) 'Hessian matrix:'
         do 2950 k = 1,n
            write(6,95) (cov(j,k),j=1,n)
 2950    continue
      end if
      CALL INVERT(N,COV,RST,DET)
      IF(DET.EQ.0) THEN
          IF(IERR.EQ.1) WRITE(6,*) 'MINI: SINGULAR MATRIX'
          DO 920 J = 1,N
              DO 921 K = 1,J-1
                  COV(J,K) = D2F(J,K)
                  COV(K,J) = D2F(J,K)
921           CONTINUE
              COV(J,J) = DABS(D2F(J,J)*(1 + BASE**LAMBDA))
920       CONTINUE
          IF(IERR.EQ.1) THEN
              DO 950 K = 1,N
                  WRITE(6,95) (COV(J,K),J=1,N)
95                FORMAT(1X,1P16E9.1)
950           CONTINUE
          END IF
          IERR = 2
          RETURN
      END IF

      if(itest.eq.1) then
         write(6,*) 'Covariance matrix:'
         do 1950 k = 1,n
            write(6,95) (cov(j,k),j=1,n)
 1950    continue
      end if

*
*     Multiply to get dA
*
      DO 922 J = 1,N
          DA(J) = 0
          DO 923 K = 1,N
              DA(J) = DA(J) - COV(J,K)*DF(K)
923       CONTINUE
922   CONTINUE
*
*     Now get new function value
*
      DO 924 J = 1,N
          A(J) = A0(J) + DA(J)
924   CONTINUE
      FNOW = FUNK(A)
      IF(ITEST.EQ.1) THEN
         WRITE(6,*) 'VARIABLES AND FUNCTION:'
         WRITE(6,*) (A(III),III=1,N),FNOW
      END IF
*
*     Test for whether the function has decreased
*     If so, adopt the new point and decrement LAMBDA
*     Else, increment LAMBDA, and get a new weight matrix
*
      IF(FNOW.LT.FTHEN) GOTO 400
C      write(6,*) lambda, fnow, fthen
C      write(6,*) (a(j),j=1,n)
      LAMBDA = LAMBDA + 1
300   CONTINUE
*
*     Normal exit, the function at A0 + DA is less than at A0
*
400   DO 925 J = 1,N
          A0(J) = A(J)
925   CONTINUE
      LAMBDA = LAMBDA - 1
*
*     Print the current status and test to see whether the function
*     has varied fractionally less than QFRAC.
*
      IF(IERR.EQ.1) CALL VPRINT(N,A0,FNOW,ITER,LAMBDA)
      IF(DABS((FTHEN - FNOW)/FNOW).LT.QFRAC) GOTO 600
500   CONTINUE
*
*     This is the final computation of the covariance matrix
*
600   CONTINUE
*
*     Quit if no minimum was found in the allowed number of iterations
*
150   CONTINUE
      IF(ITER.GE.ITERMAX) THEN
          IF(IERR.EQ.1) WRITE(6,*) '  Maximum iteration exceeded'
          IERR = 1
      ELSE
          IERR = 0
      END IF
*
*     Finally, compute the covariance matrix
*
      DO 926 J = 1,N
          DO 927 K = 1,N
              COV(J,K) = D2F(J,K) / 2
927       CONTINUE
926   CONTINUE
      CALL INVERT(N,COV,RST,DET)
      DO 928 J = 1,N
          ERR = DSQRT(DABS(COV(J,J)))
          IF(COV(J,J).LT.0) ERR = -ERR
          COV(J,J) = ERR
928   CONTINUE
      DO 929 J = 2,N
          DO 930 K = 1,J-1
              COV(J,K) = COV(J,K) / (COV(J,J)*COV(K,K))
              COV(K,J) = COV(J,K)
930       CONTINUE
929   CONTINUE

      RETURN
      END
*
*
      SUBROUTINE VPRINT(N,A,FNOW,NITER,LAMBDA)
*
*     Simple subroutine to print the current parameters
*
      REAL*8 A(N), FNOW
      WRITE(6,1000) (A(I), I = 1,N)
1000  FORMAT(' A(I) =', 1P7G11.4)
      WRITE(6,1010) FNOW,NITER,LAMBDA
1010  FORMAT(' F =',1PG15.7,'    ITER =',I4,'   LAMBDA =',I4)
      RETURN
      END

      SUBROUTINE TESTER(N,DF,D2F)
      REAL*8 D2F(16,16),DF(16)
      WRITE(6,1000) (DF(I),I=1,N)
      WRITE(6,1010) ((D2F(I,J),I=1,N),J=1,N)
1000  FORMAT(' DF: ',4(4G15.7/,'     '))
1010  FORMAT(' D2F:',64(4G15.7/,'     '))
      RETURN
      END
