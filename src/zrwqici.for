C
C**********************************************************************C
C
      SUBROUTINE RWQICI(IWQTICI,LA,KC)
C
C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 10 april 1999
C     Last modified by AEE 3/16/2007
C
C**********************************************************************C
C
C Read in spatially and/or temporally varying ICs (unit INWQICI).
C
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
c
      DIMENSION XWQV(NWQVM)
      CHARACTER TITLE(3)*79, ICICONT*3
C
      OPEN(1,FILE=ICIFN,STATUS='UNKNOWN')
      OPEN(2,FILE='wq3d.out',STATUS='UNKNOWN',ACCESS='APPEND')
C
      IF (IWQTICI.EQ.0) THEN
        READ(1,50) (TITLE(M),M=1,3)
        WRITE(2,999)
        WRITE(2,50) (TITLE(M),M=1,3)
      END IF
C
      WRITE(2,60)'* Initial conditions at ', IWQTICI,
     *  ' th day from model start'
C
      READ(1,999)
      READ(1,50) TITLE(1)
      WRITE(2,50) TITLE(1)

      DO M=2,LA*KC
        READ(1,84) I,J,K,(XWQV(NW),NW=1,NWQV)
        IF (LIJW(I,J).LT.1) THEN
          PRINT*, 'i, j, k, line# = ', I,J,K,M-1
          STOP 'ERROR!! invalid (i,j) in FILE 1'
        END IF
        L=LIJW(I,J)
        DO NW=1,NWQV
          WQV(L,K,NW)=XWQV(NW)
        END DO
        WRITE(2,84) I,J,K,(WQV(L,K,NW),NW=1,NWQV)
      END DO
C
C: WQCHLx=1/WQCHLx
C
      DO L=2,LA
        DO K=1,KC
          WQCHL(L,K) = WQV(L,K,1)*WQCHLC(IWQZMAP(L,KC)) + 
     *       WQV(L,K,2)*WQCHLD(IWQZMAP(L,KC))
     *      + WQV(L,K,3)*WQCHLG(IWQZMAP(L,KC))
          IF (IWQSRP.EQ.1) THEN
            O2WQ = MAX(WQV(L,K,19), 0.0)
            WQTAMD = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ), WQV(L,K,20) )
            WQTAMP(L,K) = WQV(L,K,20) - WQTAMD
            WQPO4D(L,K) = WQV(L,K,10) / (1.0 + WQKPO4P*WQTAMP(L,K))
            WQSAD(L,K)  = WQV(L,K,17) / (1.0 + WQKSAP*WQTAMP(L,K))
           ELSE IF (IWQSRP.EQ.2) THEN
            WQPO4D(L,K) = WQV(L,K,10) / (1.0 + WQKPO4P*SEDTWQ(L,K))
            WQSAD(L,K)  = WQV(L,K,17) / (1.0 + WQKSAP*SEDTWQ(L,K))
           ELSE
            WQPO4D(L,K) = WQV(L,K,10)
            WQSAD(L,K)  = WQV(L,K,17)
          END IF
        END DO
      END DO
C
      READ(1,52) IWQTICI, ICICONT
      WRITE(2,52) IWQTICI, ICICONT
      IF (ICICONT.EQ.'END') THEN
        CLOSE(1)
        IWQICI = 0
      END IF
C
      CLOSE(2)
c
  999 FORMAT(1X)
   50 FORMAT(A79)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)
   84 FORMAT(3I5, 21E12.4)
c
      RETURN
      END
