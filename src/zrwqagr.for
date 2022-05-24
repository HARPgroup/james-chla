C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RWQAGR(IWQTAGR)
C
C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 10 april 1999
C
C**********************************************************************C
C
C Read in spatially and/or temporally varying parameters for ALGAL
C GROWTH, RESP. & PREDATION RATES, and BASE LIGHT EXTINCT. COEFF.
C (unit INWQAGR).
C
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
c
      CHARACTER TITLE(3)*79, AGRCONT*3
C
      OPEN(1,FILE=AGRFN,STATUS='UNKNOWN')
      OPEN(2,FILE='wq3d.out',STATUS='UNKNOWN',ACCESS='APPEND')
C
      IF (IWQTAGR.EQ.0) THEN
        READ(1,50) (TITLE(M),M=1,3)
        WRITE(2,999)
        WRITE(2,50) (TITLE(M),M=1,3)
      END IF
C
      WRITE(2,60)'* Algal kinetic value at', IWQTAGR,
     *  ' th day from model start'
C
      READ(1,999)
      READ(1,50) TITLE(1)
      WRITE(2,50) TITLE(1)
      DO I=1,IWQZ
c        READ(1,51) MM,WQPMC(I),WQPMD(I),WQPMG(I),WQBMRC(I),
c     *    WQBMRD(I),WQBMRG(I),WQPRRC(I),WQPRRD(I),WQPRRG(I),WQKEB(I),
c     *    wqsdcoef(i)
c        WRITE(2,51) MM,WQPMC(I),WQPMD(I),WQPMG(I),WQBMRC(I),
c     *    WQBMRD(I),WQBMRG(I),WQPRRC(I),WQPRRD(I),WQPRRG(I),WQKEB(I),
c     *    wqsdcoef(i)
        READ(1,*) MM, WQPMC(i),WQPMD(i),WQPMG(i),WQPMM(i),WQBMRC(i),
     *    WQBMRD(i),WQBMRG(i),WQBMRM(i),WQPRRC(i),WQPRRD(i),
     *    WQPRRG(i),WQPRRM(i),WQKEB(i), wqsdcoef(i)
        WRITE(2,51) MM, WQPMC(i),WQPMD(i),WQPMG(i),WQPMM(i),WQBMRC(i),
     *    WQBMRD(i),WQBMRG(i),WQBMRM(i),WQPRRC(i),WQPRRD(i),
     *    WQPRRG(i),WQPRRM(i),WQKEB(i), wqsdcoef(i)
      END DO
C
      READ(1,52) IWQTAGR, AGRCONT
      WRITE(2,52) IWQTAGR, AGRCONT
C
      IF (AGRCONT.EQ.'END') THEN
        CLOSE(1)
        IWQAGR = 0
      END IF
C
      CLOSE(2)
C
  999 FORMAT(1X)
   50 FORMAT(A79)
   51 FORMAT(I8, 14F8.3)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)
C
      RETURN
      END
