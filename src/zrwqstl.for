c
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RWQSTL(IWQTSTL)
C
C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 10 april 1999
C
C**********************************************************************C
C
C Read in spatially and/or temporally varying parameters for SETTLING
C VELOCITIES of algae, RPOM, LPOM & particulate metal (unit INWQSTL).
C
C**********************************************************************C
C
      INCLUDE 'wq.par'  
      INCLUDE 'wqcom.cmn'
C
      CHARACTER TITLE(3)*79, STLCONT*3
C
      OPEN(1,FILE=STLFN,STATUS='UNKNOWN')
      OPEN(2,FILE='wq3d.out',STATUS='UNKNOWN',ACCESS='APPEND')
C
      IF (IWQTSTL.EQ.0) THEN
        READ(1,50) (TITLE(M),M=1,3)
        WRITE(2,999)
        WRITE(2,50) (TITLE(M),M=1,3)
      END IF
      WRITE(2,60)'* Settling velocity at  ', IWQTSTL,
     *  ' th day from model start'

      READ(1,999)
      READ(1,50) TITLE(1)
      WRITE(2,50) TITLE(1)
      DO I=1,IWQZ
        READ(1,51) MM,WQWSC(I),WQWSD(I),WQWSG(I),WQWSRP(I),
     *    WQWSLP(I),WQWSS(I), WQWSM
        WRITE(2,51) MM,WQWSC(I),WQWSD(I),WQWSG(I),WQWSRP(I),WQWSLP(I),
     *    WQWSS(I), WQWSM
      END DO

      READ(1,52) IWQTSTL, STLCONT
      WRITE(2,52) IWQTSTL, STLCONT
      IF (STLCONT.EQ.'END') THEN
        CLOSE(1)
        IWQSTL = 0
      END IF
C
      CLOSE(2)
C
  999 FORMAT(1X)
   50 FORMAT(A79)
   51 FORMAT(I8, 10F8.3)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)
C
      RETURN
      END
