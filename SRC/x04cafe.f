      SUBROUTINE X04CAF('General',' ',N,NRHS,B,LDB,'Solution(s)',IFAIL)
C     X04CAF Example Program Text
C     Mark 14 Release. NAG Copyright 1989.
C     .. Parameters ..
      INTEGER          NOUT
      PARAMETER        (NOUT=6)
      INTEGER          NMAX, LDA
      PARAMETER        (NMAX=5,LDA=NMAX)
C     .. Local Scalars ..
      INTEGER          I, IFAIL, J
C     .. Local Arrays ..
      DOUBLE PRECISION A(LDA,NMAX)
C     .. External Subroutines ..
      EXTERNAL         X04CAF
C     .. Executable Statements ..
      WRITE (NOUT,*) 'X04CAF Example Program Results'
      WRITE (NOUT,*)
C     Generate an array of data
      DO 40 J = 1, NMAX
         DO 20 I = 1, LDA
            A(I,J) = 10*I + J
   20    CONTINUE
   40 CONTINUE
C
      IFAIL = 0
C
C     Print 3 by 5 rectangular matrix
      CALL X04CAF('General',' ',3,5,A,LDA,'Example 1:',IFAIL)
C
      WRITE (NOUT,*)
C
C     Print 5 by 5 lower triangular matrix
      CALL X04CAF('Lower','Non-unit',5,5,A,LDA,'Example 2:',IFAIL)
C
      STOP
      RETURN
      END
