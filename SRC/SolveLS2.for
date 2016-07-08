CF07AEF Example Program Text
CMark 15 Release. NAG Copyright 1991.
C Parameters ..
      PROGRAM GLS
      INTEGER NIN, NOUT
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
C     .. External Subroutines ..
      OPEN(UNIT=NIN,FILE='in.txt',ACCESS='SEQUENTIAL',STATUS='OLD')
      OPEN(UNIT=NOUT,FILE='out.txt',ACCESS='SEQUENTIAL',STATUS='OLD')
      WRITE (NOUT,*) 'F07AEF Example Program Results'
C     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N, NRHS
      IF (N.LE.NMAX .AND. NRHS.LE.NRHMAX) THEN
C         Read A and B from data file
         READ (NIN,*) ((A(I,J),J=1,N),I=1,N)
         READ (NIN,*) ((B(I,J),J=1,NRHS),I=1,N)

C        Factorize A
         CALL DGETRF(N,N,A,LDA,IPIV,INFO)
C
         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
C        Compute solution
            CALL DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
C
C
            IFAIL = 0
C            CALL X04CAF('General',' ',N,NRHS,B,LDB,'Solution(s)',IFAIL)
         ELSE
            WRITE (NOUT,*) 'The factor U is singular'
         END IF
      END IF
      CLOSE(NIN)
      CLOSE(NOUN)
C
      END
