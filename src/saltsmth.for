C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SALTSMTH
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      IF(NSBMAX.GT.10)  GO TO 1001
c     IF(NSBMAX.GT.10) THEN
c       IF(MOD(NSBMAX,10).EQ.0) GO TO 1000
c      ELSE
c       GO TO 1001
c     END IF
C
C**********************************************************************C
C
      DO K=1,KC
C
      DO L=2,LA
       TVAR3S(L)=SAL(L,K)
      END DO
C
      DO NSM=1,NSBMAX
C
      DO L=2,LA
       IF(LCT(L).GT.0.AND.LCT(L).LT.9) THEN
        I=IL(L)
        J=JL(L)
        HTN=TVAR3S(LNC(L))
        HTS=TVAR3S(LSC(L))
        HTE=TVAR3S(L+1)
        HTW=TVAR3S(L-1)
        IF (IJCT(I  ,J+1).EQ.9) HTN=TVAR3S(L)
        IF (IJCT(I  ,J-1).EQ.9) HTS=TVAR3S(L)
        IF (IJCT(I+1,J  ).EQ.9) HTE=TVAR3S(L)
        IF (IJCT(I-1,J  ).EQ.9) HTW=TVAR3S(L)
        TVAR3N(L)=(1.-WSMB)*TVAR3S(L)+0.25*WSMB*(HTN+HTS+HTE+HTW)
       END IF
      END DO
C
      DO L=2,LA
       TVAR3S(L)=TVAR3N(L)
      END DO      
C
      END DO
C
      DO L=2,LA
      SAL(L,K)=TVAR3N(L)
      SAL1(L,K)=TVAR3N(L)
      END DO
C
      END DO
C
      GO TO 2000
C
C**********************************************************************C
C
C **  IMPLEMENT SPECIAL SALINITY INITIALIZATION, VERSION 1
C
C**********************************************************************C
C
C 1000 CONTINUE
C
C      DO K=1,KC
C
C      DO L=2,LA
C       TVAR3S(L)=SAL(L,K)
C      END DO
C
C      DO NSM=1,NSBMAX
C
C      DO L=2,LA
C       IF(LCT(L).GT.0.AND.LCT(L).LT.9) THEN
C        I=IL(L)
C        J=JL(L)
C        HTN=TVAR3S(LNC(L))
C        HTS=TVAR3S(LSC(L))
C        HTE=TVAR3S(L+1)
C        HTW=TVAR3S(L-1)
C        IF (IJCT(I  ,J+1).EQ.9) HTN=TVAR3S(L)
C        IF (IJCT(I  ,J-1).EQ.9) HTS=TVAR3S(L)
C        IF (IJCT(I+1,J  ).EQ.9) HTE=TVAR3S(L)
C        IF (IJCT(I-1,J  ).EQ.9) HTW=TVAR3S(L)
C        TVAR3N(L)=(1.-WSMB)*TVAR3S(L)+0.25*WSMB*(HTN+HTS+HTE+HTW)
C       END IF
C      END DO
C
C       DO L=2,LA
C        IF(SALINIT(L,K).GT.0.0) TVAR3N(L)=SALINIT(L,K)
C       END DO
C
C      DO L=2,LA
C       TVAR3S(L)=TVAR3N(L)
C      END DO      
C
C      END DO
C
C      DO L=2,LA
C       SAL(L,K)=TVAR3N(L)
C       SAL1(L,K)=TVAR3N(L)
C      END DO
C
C      WRITE(6,6000)K,NSM
C      END DO
C
C     CALL SALPLTH(1,SAL)
C
C**********************************************************************C
C
C **  IMPLEMENT SPECIAL SALINITY INITIALIZATION, VERSION 2
C
C**********************************************************************C
C
 1001 CONTINUE
C
      DO K=1,KC
C
      DO L=2,LA
       TVAR3S(L)=SAL(L,K)
      END DO
C
      DO NSM=1,NSBMAX
C
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      TVAR3N(L)=TVAR3S(L)+(WSMB/HMP(L))
     $        *( HRU(L+1)*(TVAR3S(L+1)-TVAR3S(L  ))
     $          -HRU(L  )*(TVAR3S(L  )-TVAR3S(L-1))
     $          +HRV(LN )*(TVAR3S(LN )-TVAR3S(L  ))
     $          -HRV(L  )*(TVAR3S(L  )-TVAR3S(LS )) )
      END DO
C
      DO L=2,LA
       IF(SALINIT(L,K).GT.0.0) TVAR3N(L)=SALINIT(L,K)
      END DO
C
      DO L=2,LA
       TVAR3S(L)=TVAR3N(L)
      END DO      
C
      END DO
C
      DO L=2,LA
       SAL(L,K)=TVAR3N(L)
       SAL1(L,K)=TVAR3N(L)
      END DO
C
      WRITE(6,6001)K,NSM
      END DO
C
 2000   OPEN(1,FILE='newsalt.inp',STATUS='UNKNOWN')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='newsalt.inp',STATUS='UNKNOWN')
       IONE=1
       WRITE(1,9101)IONE
       DO L=2,LC-1
       WRITE(1,9102)L,IL(L),JL(L),(SAL(L,K),K=1,KC)
       END DO
       CLOSE(1)
C
C     CALL SALPLTH(1,SAL)
C
C**********************************************************************C
C
 6000 FORMAT(' COMPLE V1 SMOOTHING LAYER ',I5,' NSM = ',I5/)
 6001 FORMAT(' COMPLE V2 SMOOTHING LAYER ',I5,' NSM = ',I5/)
 9101 FORMAT(I5)
 9102 FORMAT(3I5,12F6.2)
 2001 CONTINUE
      RETURN
      END
