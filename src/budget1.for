C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BUDGET1
c
c    Calculate BSEDBEG for PCB 
c    modifed by J.S. 3/7/2011
C
C **  ADDED BY DON KINGERY ON 15 OCTOBER 1996
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
C **  SUBROUTINES BUDGETn CALCULATE SEDIMENT BUDGET (TOTAL SEDIMENTS)
C
C**********************************************************************C
C 
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
C
C**********************************************************************C
C
C **  INITIALIZE VOLUME, SALT MASS, SEDIMENT, AND ASSOCIATED FLUXES
C
C----------------------------------------------------------------------C
C
      VOLMBEG=0.
      SMASSBEG=0.
      BSEDBEG=0.
      SSEDBEG=0.
      SEDIN=0.
      SEDOUT=0.
      VOLMOUT=0.
      SMASSOUT=0.
      VOLMIN=0.
      SMASSIN=0.
      SEDERR=0.
      TOXOUTBC=0.
      TOXINBC=0.
      TOXBEDTOP=0 
C
C**********************************************************************C
C
      DO K=1,KB
      DO L=2,LA
        SEDBT(L,K)=0.
        SNDBT(L,K)=0.
      END DO
      END DO
C
      DO NS=1,NPCB
       DO K=1,KB
       DO L=2,LA
        BSEDBEG=BSEDBEG+SCB(L)*DXYP(L)*TOXB(L,K,NS)   !bottom PCB in ugg/m^2 
       END DO
       END DO
      END DO
C
C !volbw3 is to accumulate wc-sedbed flux
C BSEDBEG Total PCB at beging
C
      DO K=1,KB
      DO L=2,LA
       VOLBW3(L,K)=0.
!       BSEDBEG=BSEDBEG+SCB(L)*DXYP(L)*(SEDBT(L,K) )
      END DO
      END DO
C

      DO NS=1,NPCB
       DO K=1,KC
        DO L=2,LA
         SSEDBEG=SSEDBEG+SCB(L)*DXYP(L)*HWQ(L)*TOX(L,K,NS)*DZC(K)
        END DO
       END DO
      END DO
C
      SEDBEG=BSEDBEG+SSEDBEG
C
C**********************************************************************C
C

      RETURN
      END
