C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALSFT (ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE CALSFT CALCULATES THE TRANSPORT OF SHELL FISH LARVAE 
C **  AT TIME LEVEL (N+1). 
C **  CALLED ONLY ON ODD THREE TIME LEVEL STEPS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      DIMENSION WTFKB(KCM),WTFKC(KCM)
C     DIMENSION HWQI(LCM)
CDHP  DIMENSION CUBTMP(LCM),CMBTMP(LCM),CLBTMP(LCM),EB(LCM),
CDHP $          VTMP(LCM),ABHWQI(LCM,KCM)
C
C**********************************************************************C
C
      DELT=DT2
C
C**********************************************************************C
C
C **  UPDATED TIME SERIES CONCENTRATION BOUNDARY CONDITIONS
C
C     CALL CALWQS(ISTL)
C      
C**********************************************************************C
C
C **  DETERMINE IF CURRENT TIME STEP IS DURING DAYLIGHT OR DARKNESS
C
      IF(ISSFLDN.GE.1) THEN
       ISDARK=1
       TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
      IF(NCSTEP.GT.0) TIME=(SECDLAST+TCON*TBEGIN)/86400.  !% J.S. 1/31/2014 
       ITIME=INT(TIME)
       RTIME=FLOAT(ITIME)
       TIMTMP=TIME-RTIME
       IF(TIMTMP.GE.TSRSF.AND.TIMTMP.LE.TSSSF) ISDARK=0
      END IF
C      
C**********************************************************************C
C
C **  DETERMINE IF LOCAL CONDITIONS ARE EBB OR FLOOD
C     FOR FLOOD, ISFEQ1(L,K)=1  FOR EBB,  ISFEQ1(L,K)=0 
C     UUU=0 FOR EBB AND 1 FOR FLOOD
C     VVV=1 FOR EBB AND O FOR FLOOD
C
      IF(ISSFLFE.GE.1) THEN
C
      IF(KC.EQ.1) THEN
        WTFKB(1)=1.
        WTFKC(1)=0.
      END IF
      IF(KC.EQ.2) THEN
        WTFKB(1)=1.0
        WTFKC(1)=0.0
        WTFKB(2)=0.0
        WTFKC(2)=1.0
      END IF
      IF(KC.EQ.3) THEN
       DO K=1,KC
        WTFKB(K)=FLOAT(KC-K)/FLOAT(KS)
        WTFKC(K)=1.0-WTFKB(K)
       END DO
      END IF
C
C **  SET SWITCHES TO EBB
C
      DO K=1,KC      
      DO L=2,LA
       UUU(L,K)=0.
       VVV(L,K)=1.
      END DO
      END DO
C
C **  RESET SWITCHES FOR FLOOD
C
      DO K=1,KC      
      DO L=2,LA
       LN=LNC(L)
       FANGTMP=ACCWFLD(L,1)*WTFKB(K)+ACCWFLD(L,2)*WTFKC(K)
       UTMP=0.5*STCUV(L)*(UWQ(L+1,K)+UWQ(L,K))
       VTMP=0.5*STCUV(L)*(VWQ(LN ,K)+VWQ(L,K))
       VELEKB=CUE(L)*UTMP+CVE(L)*VTMP+1.E-12
       VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
       CURANG=ATAN2(VELNKB,VELEKB)
       ANGDIF=ABS(FANGTMP-CURANG)
       IF(ANGDIF.LT.1.5708) THEN
         UUU(L,K)=1.
         VVV(L,K)=0.
       END IF
      END DO
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  SET UP ADVECTION FIELD
C
C **  SET ATTACHED TO BOTTOM AND NO ADVECTIVE TRANSPORT IN BOTTOM 
C **  LAYER DURING EBB IF APPROPRIATE
C
      IF(ISSFLFE.GE.1) THEN
C
      IF (SFNTBET.LT.1.) THEN
       K=1
       DO L=2,LA
       LN=LNC(L)
       UHDYWQ(L,K)=UUU(L,1)*UHDYWQ(L,1)+SFNTBET*VVV(L,1)*UHDYWQ(L,1)
       VHDXWQ(L,K)=UUU(L,1)*UHDYWQ(L,1)+SFNTBET*VVV(L,1)*VHDXWQ(L,1)
       UWQ(L,K)=UUU(L,1)*UWQ(L,1)+SFNTBET*VVV(L,1)*UWQ(L,1)
       VWQ(L,K)=UUU(L,1)*VWQ(L,1)+SFNTBET*VVV(L,1)*VWQ(L,1)
C          UHDYWQ(L+1,K)=SFNTBET*UHDYWQ(L+1,K)
C          VHDXWQ(LN ,K)=SFNTBET*VHDXWQ(LN ,K)
C          UWQ(L+1,K)=SFNTBET*UWQ(L+1,K)
C          VWQ(LN ,K)=SFNTBET*VWQ(LN ,K)
       END DO
      END IF
C
      END IF
C
      IF (ISTRAN(7).GE.1) CALL CALTRWQ (7,0,SFL,SFL2)
C
C**********************************************************************C
C
C **  SET UP VERTICAL MIGRATION AND SETTLING BEHAVIOR
C
C **  INITIALIZE VERTICAL VELOCTIY TO TIME DEPENDENT SETTLING VELOCITY 
C
      DO K=1,KS
      DO L=2,LA
      WWQ(L,K)=-WSFLSTT
      END DO
      END DO
C
      DO L=2,LA
      WWQ(L,KC)=0.
      WWQ(L,0)=0.
      END DO
C
      IF(ISSFLFE.GE.1.AND.ISSFLDN.GE.1) THEN
C
C **  DAYLIGHT CONDITIONS
      IF(ISDARK.EQ.0) THEN
      DO K=1,KS
      RABOVE=FLOAT(KC-K)/FLOAT(KC)
       DO L=2,LA
C **   DETERMINE DISTANCE TO SURFACE
       HABOVE=RABOVE*HWQ(L)
       IF(UUU(L,K).GT.0.) THEN 
C **    FLOOD CONDITION : SWIM UP TO MIN DIST BELOW SURFACE    
         IF(HABOVE.GT.DSFLMNT) WWQ(L,K)=WSFLSMT
        ELSE
C **    EBB CONDITION : CONTINUE TO SINK OR SWIM UP TO MAX DIST BL SURF 
         IF(HABOVE.GT.DSFLMXT) WWQ(L,K)=WSFLSMT
       END IF
       END DO
      END DO
      END IF
C
C **  DARK CONDITIONS
      IF(ISDARK.EQ.1) THEN
      DO K=1,KS
       DO L=2,LA
C **   FLOOD CONDITION : SWIM UP TO  SURFACE   
       WWQ(L,K)=VVV(L,K)*WWQ(L,K)+UUU(L,K)*WSFLSMT       
       END DO
      END DO
      END IF
C
      END IF
C
      IF(SFATBTT.GT.0.) THEN
      DO L=2,LA
      WWQ(L,0)=-WSFLSTT
      END DO
      END IF
C
C**********************************************************************C
C
C **  CALCULATE NET VERTICAL SWIMING OR SETTLING
C
      IF(WSFLSMT.EQ.0.) GO TO 100
C
C **  LIMIT VERTICAL SETTLING AND/OR SWIMMING FOR STABILITY
C
      DO K=0,KS
      DO L=2,LA
      WWW(L,K)=MIN(WWQ(L,K),0.)
      WWW(L,K)=ABS(WWW(L,K))
      WWQ(L,K)=MAX(WWQ(L,K),0.)
      END DO
      END DO
C
      TMPVAL=0.25/(DELT*FLOAT(KC))
      DO K=1,KS
      DO L=2,LA
      WMAXX=TMPVAL*HWQ(L)
      WWW(L,K)=MIN(WWW(L,K),WMAXX)
      WWQ(L,K)=MIN(WWQ(L,K),WMAXX)
      WWQ(L,K)=WWQ(L,K)-WWW(L,K)
      END DO
      END DO
C
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MAX(WWQ(L,K),0.)*SFL(L,K)
     $        +MIN(WWQ(L,K),0.)*SFL(L,K+1)
      END DO
      END DO
C
      IF(SFATBTT.GT.0.) THEN
      DO L=2,LA
      SFLSBOT(L)=SFLSBOT(L)-DELT*FWU(L,0)
      END DO
      END IF
C
      DO K=1,KC
      DO L=2,LA
      SFL(L,K)=SFL(L,K)
     $        +DELT*(FWU(L,K-1)-FWU(L,K))*DZIC(K)/HWQ(L)
      END DO
      END DO
C
      GO TO 200
C
  100 CONTINUE
C
C **  FULLY IMPLICIT SETTLING IF SWIMMING IS ZERO EVERYWHERE
C
C **  FULLY IMPLICIT SETTLING IN SURFACE LAYER
C
      TMPVAL=DELT*WSFLSTT
      DZCIT=TMPVAL/DZC(KC)
      DO L=2,LA
      TMPVAL1=DZCIT/HWQ(L)
      SFL(L,KC)=SFL(L,KC)/(1.+TMPVAL1)
      END DO
C
C **  FULLY IMPLICIT SETTLING IN REMAINING LAYERS
C
      IF (KC.GT.1) THEN
       DO K=KS,1,-1
       DZCIT=TMPVAL/DZC(K)
        DO L=2,LA
        TMPVAL1=DZCIT/HWQ(L)
        SFL(L,K)=(SFL(L,K)+TMPVAL1*SFL(L,K+1))/(1.+TMPVAL1)
        END DO
       END DO
      END IF
C
      IF(SFATBTT.GT.0.) THEN
      DO L=2,LA
      SFLSBOT(L)=SFLSBOT(L)+TMPVAL*SFL(L,1)
      END DO
      END IF
C
  200 CONTINUE

      DO L=2,LA
      FWU(L,0)=0.
      WWQ(L,0)=0.
      WWW(L,0)=0.
      END DO
C
C**********************************************************************C
C
C **  CALCULATE LINEAR DECAY
C
      IF (RKDSFLT.GE.0.) THEN
C
      CDYETMP=1./(1.+DELT*RKDSFLT)
C     CDYETMP=(1.-DELTD2*RKDYE)/(1.+DELTD2*RKDYE)
      DO K=1,KC
      DO L=2,LA
      SFL(L,K)=CDYETMP*SFL(L,K)
      END DO
      END DO
C
      END IF
C
C**********************************************************************C
C
      IF (KC.EQ.1) GO TO 2000
C
C**********************************************************************C
C
C **  VERTICAL DIFFUSION CALCULATION
C
      IF (ISHP.GE.1) GO TO 1000
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      HWQI(L)=1./HWQ(L)
      END DO
C
      RCDZKK=-DELT*CDZKK(1)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        CCUBTMP=RCDZKK*HWQI(L)*AB(L,1)
        CCMBTMP=1.-CCUBTMP
        EEB=1./CCMBTMP
        CU1(L,1)=CCUBTMP*EEB
        SFL(L,1)=SFL(L,1)*EEB
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=2,KS
        RCDZKMK=-DELT*CDZKMK(K)
        RCDZKK=-DELT*CDZKK(K)
        DO L=LF,LL
         CCLBTMP=RCDZKMK*HWQI(L)*AB(L,K-1)
         CCUBTMP=RCDZKK*HWQI(L)*AB(L,K)
         CCMBTMP=1.-CCLBTMP-CCUBTMP
         EEB=1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
         CU1(L,K)=CCUBTMP*EEB
         SFL(L,K)=(SFL(L,K)-CCLBTMP*SFL(L,K-1))*EEB
        END DO
       END DO
      END DO
C
      K=KC
      RCDZKMK=-DELT*CDZKMK(K)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        CCLBTMP=RCDZKMK*HWQI(L)*AB(L,K-1)
        CCMBTMP=1.-CCLBTMP
        EEB=1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
        SFL(L,K)=(SFL(L,K)-CCLBTMP*SFL(L,K-1))*EEB
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=KC-1,1,-1
        DO L=LF,LL
         SFL(L,K)=SFL(L,K)-CU1(L,K)*SFL(L,K+1)
        END DO
       END DO
      END DO
C
      GO TO 2000
C
C**********************************************************************C
C
C **  VERTICAL DIFFUSION CALCULATION OPTIMIZED FOR HP 9000 S700
C
 1000 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO L=1,LA
      HWQI(L)=1./HWQ(L)
      END DO
C
CDHP  DO K=1,KS
CDHP  CALL vec_$mult_vector(HWQI(1),AB(1,K),LA,ABHWQI(1,K))
CDHP  END DO
C
      K=1
      RUNIT=1.0
      RTMP=-DELT*CDZKK(1)
C     DO L=1,LA
C     CUBTMP(L)=RTMP*HWQI(L)*AB(L,1)=RTMP*ABHWQI(L,1)
C     END DO
CDHP  CALL vec_$mult_constant(ABHWQI(1,1),LA,RTMP,CUBTMP(1))
C     DO L=1,LA
C     CMBTMP(L)=1.-CUBTMP(L)
C     END DO
CDHP  CALL vec_$SUB_constant(CUBTMP(1),LA,RUNIT,CMBTMP(1))
CDHP  DO L=1,LA
CDHP  EB(L)=1./CMBTMP(L)
CDHP  END DO
C     DO L=1,LA
C     CU1(L,1)=CUBTMP(L)*EB(L)
C     END DO
CDHP  CALL vec_$mult_vector(CUBTMP(1),EB(1),LA,CU1(1,1))
C     IF (ISTRAN(1).GE.1) THEN
C       DO L=1,LA
C       SFL(L,1)=SFL(L,1)*EB(L)
C       END DO
C     END IF
CDHP  IF (ISTRAN(1).GE.1)
CDHP $  CALL vec_$mult_vector(SFL(1,1),EB(1),LA,SFL(1,1)) 
C
      DO K=2,KS
      RTMP1=-DELT*CDZKMK(K)
      RTMP2=-DELT*CDZKK(K)
C      DO L=1,LA
C      CLBTMP(L)=RTMP1*AB(L,K-1)*HWQI(L)=RTMP1*ABHWQI(L,K-1)
C      CUBTMP(L)=RTMP2*AB(L,K)*HPI(L)=RTMP2*ABHWQI(L,K)
C      END DO
CDHP  CALL vec_$mult_constant(ABHWQI(1,K-1),LA,RTMP1,CLBTMP(1))
CDHP  CALL vec_$mult_constant(ABHWQI(1,K),LA,RTMP2,CUBTMP(1))
C      DO L=1,LA
C      CMBTMP(L)=1.-CLBTMP(L)-CUBTMP(L)
C      END DO
CDHP  CALL vec_$SUB_constant(CUBTMP(1),LA,RUNIT,CMBTMP(1))
CDHP  CALL vec_$SUB_vector(CMBTMP(1),CLBTMP(1),LA,CMBTMP(1))
C      DO L=1,LA
C      CMBTMP(L)=CMBTMP(L)-CLBTMP(L)*CU1(L,K-1)
C      END DO
CDHP  CALL vec_$mult_rSUB_vector(CLBTMP(1),CU1(1,K-1),CMBTMP(1),LA,
CDHP $                           CMBTMP(1))
CDHP   DO L=1,LA
CDHP   EB(L)=1./CMBTMP(L)
CDHP   END DO
C      DO L=1,LA
C      CU1(L,K)=CUBTMP(L)*EB(L)
C      END DO
CDHP  CALL vec_$mult_vector(CUBTMP(1),EB(1),LA,CU1(1,K))  
C      DO L=1,LA
C      CLBTMP(L)=-1.*CLBTMP(L)*EB(L)
C      END DO
CDHP  CALL vec_$mult_vector(CLBTMP(1),EB(1),LA,CLBTMP(1))
CDHP  CALL vec_$neg(CLBTMP(1),LA,CLBTMP(1)) 
      IF (ISTRAN(1).GE.1) THEN
C       DO L=1,LA
C       SFL(L,K)=EB(L)*SFL(L,K)+CLBTMP(L)*SFL(L,K-1)
C       END DO
CDHP    CALL vec_$mult_vector(CLBTMP(1),SFL(1,K-1),LA,VTMP(1))
CDHP    CALL vec_$mult_add_vector(EB(1),SFL(1,K),VTMP(1),LA,SFL(1,K))  
      END IF
      END DO
C
      K=KC
      RTMP=-DELT*CDZKMK(KC)
C      DO L=1,LA
C      CLBTMP(L)=RTMP*AB(L,KC-1)*HWQI(L)=RTMP*ABHWQI(L,KS)
C      END DO
CDHP  CALL vec_$mult_constant(ABHWQI(1,KS),LA,RTMP,CLBTMP(1))
C      DO L=1,LA
C      CMBTMP(L)=1.-CLBTMP(L)
C      END DO
CDHP  CALL vec_$SUB_constant(CLBTMP(1),LA,RUNIT,CMBTMP(1))
C      DO L=1,LA
C      CMBTMP(L)=CMBTMP(L)-CLBTMP(L)*CU1(L,K-1)
C      END DO
CDHP  CALL vec_$mult_rSUB_vector(CLBTMP(1),CU1(1,KS),CMBTMP(1),LA,
CDHP $                           CMBTMP(1))
CDHP   DO L=1,LA
CDHP   EB(L)=1./CMBTMP(L)
CDHP   END DO
C      DO L=1,LA
C      CLBTMP(L)=-1.*CLBTMP(L)*EB(L)
C      END DO
CDHP  CALL vec_$mult_vector(CLBTMP(1),EB(1),LA,CLBTMP(1))
CDHP  CALL vec_$neg(CLBTMP(1),LA,CLBTMP(1))
      IF (ISTRAN(1).GE.1) THEN
C       DO L=1,LA
C       SFL(L,K)=EB(L)*SFL(L,K)+CLBTMP(L)*SFL(L,K-1)
C       END DO
CDHP    CALL vec_$mult_vector(CLBTMP(1),SFL(1,KS),LA,VTMP(1))
CDHP    CALL vec_$mult_add_vector(EB(1),SFL(1,KC),VTMP(1),LA,SFL(1,KC))    
      END IF
C
      DO K=KS,1,-1
C     IF (ISTRAN(5).GE.1) THEN
C       DO L=1,LA
C       SFL(L,K)=SFL(L,K)-CU1(L,K)*SFL(L,K+1)
C       END DO
C     END IF
CDHP  IF (ISTRAN(5).GE.1) 
CDHP $  CALL vec_$mult_rSUB_vector(CU1(1,K),SFL(1,K+1),SFL(1,K),LA,   
CDHP $                             SFL(1,K))
      END DO
C
C**********************************************************************C
C
C **  UPDATE SHELL FISH LARVAE CONCENTRATIONS
C
 2000 CONTINUE
C
      DO K=1,KC
       DO L=2,LA
        SFL2(L,K)=SFL(L,K)
       END DO
      END DO
C
C**********************************************************************C
C
      RETURN
      END
