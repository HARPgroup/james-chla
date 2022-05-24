C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BUDGET3
C
c     calculate SEDIN of budget.out, to know how much sediment into the system, Ji, 10/24/00
c
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
      INCLUDE 'wq.par'      
      INCLUDE 'efdc.cmn'
      INCLUDE 'wqcom.cmn' 
C
C**********************************************************************C
C
C     DIMENSION CONT(LCM,KCM)
C
C**********************************************************************C
C
C **  ACCUMULATE INTERNAL SOURCES AND SINKS
C
C----------------------------------------------------------------------C
C
      VOLCONT=0.
      VOLMAST=0.
C
      DO L=2,LA
      VOLMIN=VOLMIN+QSUME(L)
      VOLCONT=VOLCONT+QSUME(L)
      END DO
C
C----------------------------------------------------------------------C
C
C    ACCUMULATE INTERNAL SOURCES AND SINKS FOR PCB
C
      if(irunpcb.EQ.1) then
C     
C For total PCB
C      
       DO NS=1,NPCBpt
        L = LIJW(IPCBpt(NS),JPCBpt(NS))  
         DO k=1,kc
         DO NW=1,NPCB
	    SEDIN=SEDIN+WQWPSL(L,K,NW)           
         ENDDO
         ENDDO
       END DO
!        M=5
!       DO NS=1,NQSIJ
!        L=LQS(NS)
!        DO K=1,KC
!        SEDIN=SEDIN
!     $           +MAX(QSS(K,NS),0.)*CQS(K,NS,M)
!     $           +MIN(QSS(K,NS),0.)*TOX1(L,K,NN)
!     $           +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
!     $           +MIN(QSERT(K,NQSTMP),0.)*TOX(L,K,NN)
!        END DO
!       END DO
       
      else
      
       M=5
       DO NT=1,NTOX
       M=MSVTOX(NT)      
       DO NS=1,NQSIJ
        L=LQS(NS)
        NQSTMP=NQSERQ(NS)
        NCSTMP=NCSERQ(NS,M)       
        DO K=1,KC
        SEDIN=SEDIN
     $           +MAX(QSS(K,NS),0.)*CQS(K,NS,M)
     $           +MIN(QSS(K,NS),0.)*TOX1(L,K,NN)
     $           +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
     $           +MIN(QSERT(K,NQSTMP),0.)*TOX(L,K,NN)
        END DO
       END DO
       ENDDO

      endif
      
      RETURN
      END
