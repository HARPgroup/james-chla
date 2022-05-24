C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALPCB (ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE CALCULATES THE CONCENTRATION OF DISSOLVED AND
C **  SUSPENDED CONSTITUTENTS, INCLUDING SALINITY, TEMPERATURE, DYE AND
C **  AND SUSPENDED SEDIMENT AT TIME LEVEL (N+1). THE VALUE OF ISTL
C **  INDICATES THE NUMBER OF TIME LEVELS IN THE STEP
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      INCLUDE 'wq.par' 
      INCLUDE 'wqcom.cmn'   
!      DIMENSION EEB(LCM),CCLBTMP(LCM)
C
C**********************************************************************C
C
      DELT=DT2
      S3TL=1.0
      S2TL=0.0

      DELTD2=DELT
      
      TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCON
      
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
C     DELTD2=0.5*DELT
C
C**********************************************************************C
C
C **  3D ADVECTI0N TRANSPORT CALCULATION
C
C----------------------------------------------------------------------C

        NTOX1=NPCB
        
        DO NT=1,NTOX1
         M=MSVTOX(NT)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=TOX(L,K,NT)
          TVAR2S(L,K)=TOX(L,K,NT)
         END DO
         END DO
         
         DO K=1,KC
         DO L=2,LA
         FQC(L,K)=WQWPSL(L,K,NT)*DZIC(K)    ! get point sources ug/s WQWPSL is set in sub pcbpoint 3/1/2017 JS
         ENDDO
         ENDDO
         
         CALL CALTRPCB(M,5,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          TOX1(L,K,NT)=TVAR1S(L,K)
          TOX(L,K,NT)=TVAR2S(L,K)
         END DO
         END DO
        END DO

C

C**********************************************************************C
C
C **  VERTICAL DIFFUSION IMPLICIT HALF STEP CALCULATION
C
C----------------------------------------------------------------------C
C
!
! VERTICAL DIFFUSION CALCULATION
!
      DO L=2,LA
       HWQI(L)=1./HWQ(L)
      END DO

        IF (KC.EQ.1) GO TO 2001
!
        RCDZKK=-DELT*CDZKK(1)
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,CCUBTMP,CCMBTMP,EEB) 
        DO L=2,LA
         CCUBTMP=RCDZKK*HWQI(L)*AB(L,1)*F_AB
         CCMBTMP=1.-CCUBTMP
         EEB=1./CCMBTMP
         CU1(L,1)=CCUBTMP*EEB  
         DO NT=1,NTOX1
          TOX(L,1,NT)=TOX(L,1,NT)*EEB
         END DO
        END DO
C
        DO K=2,KS
         RCDZKMK=-DELT*CDZKMK(K)
         RCDZKK=-DELT*CDZKK(K)
         DO L=2,LA
          CCLBTMP=RCDZKMK*HWQI(L)*AB(L,K-1)*F_AB
          CCUBTMP=RCDZKK*HWQI(L)*AB(L,K)*F_AB
          CCMBTMP=1.-CCLBTMP-CCUBTMP
          EEB=1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
          CU1(L,K)=CCUBTMP*EEB
          DO NT=1,NTOX1
           TOX(L,K,NT)=(TOX(L,K,NT)-CCLBTMP*TOX(L,K-1,NT))*EEB
          END DO     
         END DO
        END DO
        K=KC
        RCDZKMK=-DELT*CDZKMK(K)
        DO L=2,LA
         CCLBTMP=RCDZKMK*HWQI(L)*AB(L,K-1)*F_AB
         CCMBTMP=1.-CCLBTMP
         EEB=1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
         DO NT=1,NTOX1
         TOX(L,K,NT)=(TOX(L,K,NT)-CCLBTMP*TOX(L,K-1,NT))*EEB
         END DO
        END DO

        DO K=KC-1,1,-1
         DO L=2,LA
          DO NT=1,NTOX1
          TOX(L,K,NT)=TOX(L,K,NT)-CU1(L,K)*TOX(L,K+1,NT)
          END DO     
         END DO
        END DO
        goto 1500
2001   continue
C=========================================
!      DO L=2,LA
!      HPI(L)=1./HP(L)
!      END DO
      
!      IF(KC.EQ.1) GO TO 1500
C
!      RCDZKK=-DELTD2*CDZKK(1)

!      DO L=2,LA
!       CCUBTMP=RCDZKK*HPI(L)*AB(L,1)
!       CCMBTMP=1.-CCUBTMP
!       EEB(L)=1./CCMBTMP
!       CU1(L,1)=CCUBTMP*EEB(L)
!      END DO

!      DO NT=1,NTOX1
!       DO L=2,LA
!        TOX(L,1,NT)=TOX(L,1,NT)*EEB(L)
!       END DO
!      END DO
!C
!      DO K=2,KS
!       RCDZKMK=-DELT*CDZKMK(K)
!       RCDZKK=-DELT*CDZKK(K)
!       DO L=2,LA
!        CCLBTMP(L)=RCDZKMK*HPI(L)*AB(L,K-1)
!        CCUBTMP=RCDZKK*HPI(L)*AB(L,K)
!        CCMBTMP=1.-CCLBTMP(L)-CCUBTMP
!        EEB(L)=1./(CCMBTMP-CCLBTMP(L)*CU1(L,K-1))
!        CU1(L,K)=CCUBTMP*EEB(L)
!       END DO
!
!       DO NT=1,NTOX1
!        DO L=2,LA
!         TOX(L,K,NT)=(TOX(L,K,NT)-CCLBTMP(L)*TOX(L,K-1,NT))*EEB(L)
!        END DO
!       END DO
!      END DO
!
!C
!      K=KC
!      RCDZKMK=-DELT*CDZKMK(K)
!      DO L=2,LA
!        CCLBTMP(L)=RCDZKMK*HPI(L)*AB(L,K-1)
!        CCMBTMP=1.-CCLBTMP(L)
!        EEB(L)=1./(CCMBTMP-CCLBTMP(L)*CU1(L,K-1))
!      END DO
!      DO NT=1,NTOX1
!        DO L=2,LA
!         TOX(L,K,NT)=(TOX(L,K,NT)-CCLBTMP(L)*TOX(L,K-1,NT))*EEB(L)
!        END DO
!      END DO
!C
!      DO K=KC-1,1,-1
!        DO NT=1,NTOX1
!         DO L=2,LA
!          TOX(L,K,NT)=TOX(L,K,NT)-CU1(L,K)*TOX(L,K+1,NT)
!         END DO
!        END DO
!      END DO
!
!C
!C**********************************************************************C
C
 1500 CONTINUE
C
C**********************************************************************C
C
C **  DATA ASSIMILATION
C
C**********************************************************************C
C
C **  VERTICAL DIFFUSION CALCULATION OPTIMIZED FOR HP 9000 S700
C
 1000 CONTINUE
C
C----------------------------------------------------------------------C
C
      RETURN
      END