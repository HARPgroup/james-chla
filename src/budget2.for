C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BUDGET2
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
C **  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES
C   NCBS,NCBW,NCBE,NCBN
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)

       VOLMOUT=VOLMOUT-VHDXWQ(LN,K)*DZC(K)
C
       DO NT=1,NPCB     
       TOXOUTBC=TOXOUTBC-MIN(VHDXWQ(LN,K),0.)*TOX(LN,K,NT)*DZC(K)   !out>0
       TOXINBC =TOXINBC -MAX(VHDXWQ(LN,K),0.)*TOX(L,K,NT)*DZC(K)    !in <0
       
        SEDOUT=SEDOUT-MIN(VHDXWQ(LN,K),0.)*TOX(LN,K,NT)*DZC(K)
     $               -MAX(VHDXWQ(LN,K),0.)*TOX(L,K,NT)*DZC(K)
       END DO
C
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)

      VOLMOUT=VOLMOUT-UHDYWQ(L+1,K)*DZC(K)

      DO NT=1,NPCB
       TOXOUTBC=TOXOUTBC-MIN(UHDYWQ(L+1,K),0.)*TOX(L+1,K,NT)*DZC(K)   !out >0
       TOXINBC =TOXINBC -MAX(UHDYWQ(L+1,K),0.)*TOX(L,K,NT)*DZC(K)     !in <0
      
       SEDOUT=SEDOUT-MIN(UHDYWQ(L+1,K),0.)*TOX(L+1,K,NT)*DZC(K)
     $               -MAX(UHDYWQ(L+1,K),0.)*TOX(L,K,NT)*DZC(K)
      END DO

      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)

      VOLMOUT=VOLMOUT+UHDYWQ(L,K)*DZC(K)
C
      DO NT=1,NPCB
      TOXOUTBC=TOXOUTBC+MAX(UHDYWQ(L,K),0.)*TOX(L-1,K,NT)*DZC(K) 
      TOXINBC =TOXINBC +MIN(UHDYWQ(L,K),0.)*TOX(L,K,NT)*DZC(K)

      SEDOUT=SEDOUT+MIN(UHDYWQ(L,K),0.)*TOX(L,K,NT)*DZC(K)           ! in <0
     $               +MAX(UHDYWQ(L,K),0.)*TOX(L-1,K,NT)*DZC(K)       ! out >0
      END DO

      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)

      VOLMOUT=VOLMOUT+VHDXWQ(L,K)*DZC(K)
C
      DO NT=1,NPCB
       TOXOUTBC=TOXOUTBC+MAX(VHDXWQ(L,K),0.)*TOX(LS,K,NT)*DZC(K) !out >0
       TOXINBC =TOXINBC +MIN(VHDXWQ(L,K),0.)*TOX(L,K,NT)*DZC(K)  !in <0

       SEDOUT=SEDOUT+MIN(VHDXWQ(L,K),0.)*TOX(L,K,NT)*DZC(K)
     $               +MAX(VHDXWQ(L,K),0.)*TOX(LS,K,NT)*DZC(K)
      END DO
C
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  ACCUMULATE FLUX OF PCB in and out of top layer of sediment
C
       DO NS=1,NPCB
        DO L=2,LA 
         VOLBW3(L,2)=VOLBW3(L,2)+SCB(L)*DXYP(L)*TOXF(L,0,NS)*TOX(L,1,NS)   
     &             +TOXFB(L,2,NS)*DXYP(L)*TOXB(L,2,NS)
        END DO
       END DO 
C
C**********************************************************************C
C
      RETURN
      END
