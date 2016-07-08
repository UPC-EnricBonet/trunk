CF07AEF Example Program Text
CMark 15 Release. NAG Copyright 1991.
C Parameters ..
      PROGRAM GLS
      INTEGER NIN, NOUT, MODE, NN
      PARAMETER (NIN=5,NOUT=6)
      INTEGER NMAX, LDA, NRHMAX, LDB
      PARAMETER (NMAX=8,LDA=NMAX,NRHMAX=NMAX,LDB=NMAX)
      CHARACTER TRANS
      PARAMETER (TRANS='N')
C     Local Scalars ..
      INTEGER I, IFAIL, INFO, J, N, NRHS
C     .. Local Arrays ..
      DOUBLE PRECISION A(LDA,NMAX), B(LDB,NRHMAX)
      INTEGER IPIV(NMAX)
      MODE =0
C     .. External Subroutines ..
      OPEN(UNIT=NIN,FILE='in.txt',ACCESS='SEQUENTIAL',STATUS='OLD')

      WRITE(*,*) 'Reading'
      READ (NIN,*)
      READ (NIN,*) N, NRHS
      WRITE (*,*) N, NRHS
      IF (N.LE.NMAX .AND. NRHS.LE.NRHMAX) THEN
C         Read A and B from data file
         DO I=1,N
            READ (NIN,100) (A(I,J),J=1,N)
            IF (MODE.EQ.0)THEN
               WRITE (*,100) (A(I,J),J=1,N)
            ENDIF
         ENDDO
         DO I=1,N
            READ (NIN,200) (B(I,J),J=1,NRHS)
            IF (MODE.EQ.0)THEN
               WRITE (*,200) (B(I,J),J=1,NRHS)
            ENDIF
         ENDDO
 100     FORMAT(3(F6.2,1x),F6.2)
 200     FORMAT((F6.2,1x),F6.2)
         NN = NRHS
         WRITE(*,*) NN
         CLOSE(NIN)
C        Factorize A
         CALL DGETRF(N,N,A,LDA,IPIV,INFO)
C
         IF (INFO.EQ.0) THEN
C        Compute solution
            CALL DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
            IFAIL = 0
            WRITE(*,*) 'Calculating'
         ELSE
            WRITE (*,*) 'The factor U is singular'
         END IF
         CALL PRINTING(A,B,N,NN,NOUT,MODE, LDA,NMAX , LDB, NRHMAX)
      END IF
C
      END

      SUBROUTINE PRINTING(A,B,N,NN,NOUT,MODE,LDA,NMAX , LDB, NRHMAX)
      
      INTEGER N, NN, NOUT, MODE, LDA,NMAX , LDB, NRHMAX 
C     Local Scalars ..
      INTEGER I, J
C     .. Local Arrays ..
      DOUBLE PRECISION A(LDA,NMAX), B(LDB,NRHMAX)
C      WRITE(*,*)' N=',N, '  NN=',NN, '  MODE=',MODE, 'NOUT= ',NOUT
C      DO I=1,N
C         WRITE (*,100) (A(I,J),J=1,N)
C      ENDDO
      DO I=1,N
         WRITE (*,200) (B(I,J),J=1,NN)
      ENDDO
      OPEN(UNIT=NOUT,FILE='out.txt',STATUS='UNKNOWN')
      WRITE (NOUT,*) 'Example Program Results'
      DO I=1,N
         WRITE (NOUT,200) (B(I,J),J=1,2)
      ENDDO
 100     FORMAT(3(F6.2,1x),F6.2)
 200     FORMAT((F6.2,1x),F6.2)
      WRITE(NOUT,*)'FIN'
      CLOSE(NOUT)
      RETURN
      END    
