C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALDISP3
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      DIMENSION UP(KCM),VP(KCM),CDISP(MGM,MGM),
     $          CDISPI(KCM,KCM),CCTMP(KCM,KCM),CCUTMP(KCM),CCVTMP(KCM),
     $          CSOL(MGM),INDX(MGM)
C
C**********************************************************************C
C
C **  INITIALIZE ON FIRST CALL
C
      IF (N.EQ.NDISP) THEN
C
      DO L=1,LC
      DO KK=1,KC
      DO K=1,KC
      BDISP(K,KK,L)=0.
      END DO
      END DO
      END DO
C
      DO L=1,LC
      DO K=1,KC
      BDISP(K,K,L)=1.
      FUDISP(K,L)=0.
      FVDISP(K,L)=0.
      END DO
      END DO
C
      DO L=1,LC
      DXXTCA(L)=0.
      DXYTCA(L)=0.
      DYXTCA(L)=0.
      DYYTCA(L)=0.
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  CALCULATE VERTICAL DIFFUSION MATRIX AND INVERSE AND ACCUMULATE
C **  DATA FOR DISPERSION CALCULATION
C
      DELT=FLOAT(NTSTBC)*DT
C
      DO L=2,LA
      IF (LCT(L).EQ.5.AND.SPB(L).NE.0.) THEN
      LN=LNC(L)
C
      DO K=1,KC
      UP(K)=0.5*(U(L,K)+U(L+1,K))
      VP(K)=0.5*(V(L,K)+V(LN,K))
      END DO
C
      UAVG=0.
      VAVG=0.
      DO K=1,KC
      UAVG=UAVG+DZC(K)*UP(K)
      VAVG=VAVG+DZC(K)*VP(K)
      END DO
C
      DO K=1,KC
      UP(K)=UP(K)-UAVG
      VP(K)=VP(K)-VAVG
      END DO
C
      DO KK=1,KC
      DO K=1,KC
      CDISP(K,KK)=0.
      END DO
      END DO
C
      CUTMP=-DELT*CDZKK(1)*AB(L,1)*HPI(L)
      CMTMP=1.-CUTMP
      CDISP(1,1)=CMTMP*DZC(1)
      CDISP(1,2)=CUTMP*DZC(1)
C
      DO K=2,KS
      CLTMP=-DELT*CDZKMK(K)*AB(L,K-1)*HPI(L)
      CUTMP=-DELT*CDZKK(K)*AB(L,K)*HPI(L)
      CMTMP=1.-CLTMP-CUTMP
      CDISP(K,K-1)=CLTMP*DZC(K)
      CDISP(K,K)=CMTMP*DZC(K)
      CDISP(K,K+1)=CUTMP*DZC(K)
      END DO
C
      CLTMP=-DELT*CDZKMK(KC)*AB(L,KS)*HPI(L)
      CMTMP=1.-CLTMP
      CDISP(KC,KS)=CLTMP*DZC(KC)
      CDISP(KC,KC)=CMTMP*DZC(KC)
C
      CALL LUDCMP(CDISP,KC,MGM,INDX,DDD)
C
      DO KK=1,KC
      DO K=1,KC
      CDISPI(K,KK)=0.
      END DO
      END DO
C
      DO K=1,KC
      CDISPI(K,K)=1.
      END DO
C
      DO K=1,KC
      CALL LUBKSB(CDISP,KC,MGM,INDX,CDISPI(1,K))
      END DO
C
      DO KK=1,KC
      DO K=1,KC
      CDISPI(K,KK)=CDISPI(K,KK)*DZC(K)
      END DO
      END DO
C
      DO KK=1,KC
      DO K=1,KC
      CTMP=0.
      DO KT=1,KC
      CTMP=CTMP+CDISPI(K,KT)*BDISP(KT,KK,L)
      END DO
      CCTMP(K,KK)=CTMP
      END DO
      END DO
C
      DO KK=1,KC
      DO K=1,KC
      BDISP(K,KK,L)=CCTMP(K,KK)
      END DO
      END DO
C
      DO K=1,KC
      CCUTMP(K)=FUDISP(K,L)-DT*UP(K)/HMIN
      CCVTMP(K)=FVDISP(K,L)-DT*VP(K)/HMIN
      END DO
C
      DO K=1,KC
      CCUU=0.
      CCVV=0.
      DO KK=1,KC
      CCUU=CCUU+CDISPI(K,KK)*CCUTMP(KK)
      CCVV=CCVV+CDISPI(K,KK)*CCVTMP(KK)
      END DO
      FUDISP(K,L)=CCUU
      FVDISP(K,L)=CCVV
      END DO
C
      CCUU=0.
      CCVV=0.
      CCUV=0.
      CCVU=0.
      DO K=1,KC
      CCUU=CCUU+DZC(K)*UP(K)*FUDISP(K,L)
      CCUV=CCUV+DZC(K)*UP(K)*FVDISP(K,L)
      CCVU=CCVU+DZC(K)*VP(K)*FUDISP(K,L)
      CCVV=CCVV+DZC(K)*VP(K)*FVDISP(K,L)
      END DO
      DXXTCA(L)=DXXTCA(L)+CCUU*HP(L)
      DXYTCA(L)=DXYTCA(L)+CCUV*HP(L)
      DYXTCA(L)=DYXTCA(L)+CCVU*HP(L)
      DYYTCA(L)=DYYTCA(L)+CCVV*HP(L)
C
      DO K=1,KC
      CCUU=0.
      CCVV=0.
      DO KK=1,KC
      CCUU=CCUU+DZC(KK)*UP(KK)*BDISP(KK,K,L)
      CCVV=CCVV+DZC(KK)*VP(KK)*BDISP(KK,K,L)
      END DO
      CUDISPT(K,L)=CCUU*HP(L)
      CVDISPT(K,L)=CCVV*HP(L)
      END DO
C
      END IF
      END DO
C
      IF (N.LT.NTS) RETURN
C
C**********************************************************************C
C
C **  COMPLETE CALCULATION OF DISPERSION COEFFICIENTS
C
      TPNN=TPN/FLOAT(NTSTBC)
C
      DO L=2,LA
      IF (LCT(L).EQ.5.AND.SPB(L).NE.0.) THEN
C
      DO KK=1,KC
      DO K=1,KC
      CDISP(K,KK)=-BDISP(K,KK,L)
      END DO
      END DO
C
      DO K=1,KC
      CDISP(K,K)=1.+CDISP(K,K)
      END DO
C
      CALL LUDCMP(CDISP,KC,MGM,INDX,DDD)
C
      DO K=1,KC
      CSOL(K)=FUDISP(K,L)
      END DO
      CALL LUBKSB(CDISP,KC,MGM,INDX,CSOL)
      CCUU=0.
      CCVU=0.
      DO K=1,KC
      CCUU=CCUU+CUDISPT(K,L)*CSOL(K)
      CCVU=CCVU+CVDISPT(K,L)*CSOL(K)
      END DO
      DXXTCA(L)=-(DXXTCA(L)+CCUU)*HMIN/TPNN
      DYXTCA(L)=-(DYXTCA(L)+CCVU)*HMIN/TPNN
C
      DO K=1,KC
      CSOL(K)=FVDISP(K,L)
      END DO
      CALL LUBKSB(CDISP,KC,MGM,INDX,CSOL)
      CCVV=0.
      CCUV=0.
      DO K=1,KC
      CCVV=CCVV+CVDISPT(K,L)*CSOL(K)
      CCUV=CCUV+CUDISPT(K,L)*CSOL(K)
      END DO
      DYYTCA(L)=-(DYYTCA(L)+CCVV)*HMIN/TPNN
      DXYTCA(L)=-(DXYTCA(L)+CCUV)*HMIN/TPNN
C
      END IF
      END DO
C
      DO L=2,LA
      DXXTCA(L)=DXXTCA(L)/HLPF(L)
      DXYTCA(L)=DXYTCA(L)/HLPF(L)
      DYXTCA(L)=DYXTCA(L)/HLPF(L)
      DYYTCA(L)=DYYTCA(L)/HLPF(L)
      END DO
C
C**********************************************************************C
C
C **  ADJUST DISPERSON TENSOR COMPONENTS
C
      OPEN(88,FILE='disdia.out',STATUS='UNKNOWN')
      CLOSE(88,STATUS='DELETE')
      OPEN(88,FILE='disdia.out',STATUS='UNKNOWN')
C
      DMAX=10000.
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      DXXTMP=DXXTCA(L)
      IF(DXXTMP.LT.0.OR.DXXTMP.GT.DMAX) THEN
C
      WRITE(88,8881)IL(L),JL(L),DXXTMP
      DXXTCA(L)=0.
      LN=LNC(L)
      LS=LSC(L)
      WTX=0.
      DXXWEST=0.
      DXXEAST=0.
      DXXSOUT=0.
      DXXNORT=0.
C
      IF(SUB(L).NE.0) THEN
       DXXWEST=DXXTCA(L-1)
       IF(DXXWEST.LT.0.OR.DXXWEST.GT.DMAX) THEN
        DXXWEST=0.
       ELSE
        WTX=WTX+1
       END IF
      END IF
C
      IF(SUB(L+1).NE.0) THEN
       DXXEAST=DXXTCA(L+1)
       IF(DXXEAST.LT.0.OR.DXXEAST.GT.DMAX) THEN
        DXXEAST=0.
       ELSE
        WTX=WTX+1
       END IF
      END IF
C
      IF(SVB(L).NE.0) THEN
       DXXSOUT=DXXTCA(LS)
       IF(DXXSOUT.LT.0.OR.DXXSOUT.GT.DMAX) THEN
        DXXSOUT=0.
       ELSE
        WTX=WTX+1
       END IF
      END IF
C
      IF(SVB(LN).NE.0) THEN
       DXXNORT=DXXTCA(LN)
       IF(DXXNORT.LT.0.OR.DXXNORT.GT.DMAX) THEN
        DXXNORT=0.
       ELSE
        WTX=WTX+1
       END IF
      END IF
C
      IF(WTX.NE.0) DXXTCA(L)=(DXXWEST+DXXEAST+DXXSOUT+DXXNORT)/WTX
C
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      DXYTMP=DXYTCA(L)
      IF(ABS(DXYTMP).GT.DMAX) THEN
C
      WRITE(88,8882)IL(L),JL(L),DXYTMP
      DXYTCA(L)=0.
      LN=LNC(L)
      LS=LSC(L)
      WTX=0.
      DXYWEST=0.
      DXYEAST=0.
      DXYSOUT=0.
      DXYNORT=0.
C
      IF(SUB(L).NE.0) THEN
       DXYWEST=DXYTCA(L-1)
       IF(ABS(DXYWEST).GT.DMAX) THEN
        DXYWEST=0.
       ELSE
        WTX=WTX+1
       END IF
      END IF
C
      IF(SUB(L+1).NE.0) THEN
       DXYEAST=DXYTCA(L+1)
       IF(ABS(DXYEAST).GT.DMAX) THEN
        DXYEAST=0.
       ELSE
        WTX=WTX+1
       END IF
      END IF
C
      IF(SVB(L).NE.0) THEN
       DXYSOUT=DXYTCA(LS)
       IF(ABS(DXYSOUT).GT.DMAX) THEN
        DXYSOUT=0.
       ELSE
        WTX=WTX+1
       END IF
      END IF
C
      IF(SVB(LN).NE.0) THEN
       DXYNORT=DXYTCA(LN)
       IF(ABS(DXYNORT).GT.DMAX) THEN
        DXYNORT=0.
       ELSE
        WTX=WTX+1
       END IF
      END IF
C
      DXYTCA(L)=(DXYWEST+DXYEAST+DXYSOUT+DXYNORT)/WTX
C
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      DYXTMP=DYXTCA(L)
      IF(DYXTMP.GT.DMAX) THEN
C
      WRITE(88,8883)IL(L),JL(L),DYXTMP
      DYXTCA(L)=0.
      LN=LNC(L)
      LS=LSC(L)
      WTY=0.
      DYXWEST=0.
      DYXEAST=0.
      DYXSOUT=0.
      DYXNORT=0.
C
      IF(SUB(L).NE.0) THEN
       DYXWEST=DYXTCA(L-1)
       IF(ABS(DYXWEST).GT.DMAX) THEN
        DYXWEST=0.
       ELSE
        WTY=WTY+1
       END IF
      END IF
C
      IF(SUB(L+1).NE.0) THEN
       DYXEAST=DYXTCA(L+1)
       IF(ABS(DYXEAST).GT.DMAX) THEN
        DYXEAST=0.
       ELSE
        WTY=WTY+1
       END IF
      END IF
C
      IF(SVB(L).NE.0) THEN
       DYYSOUT=DYYTCA(LS)
       IF(ABS(DYXSOUT).GT.DMAX) THEN
        DYXSOUT=0.
       ELSE
        WTY=WTY+1
       END IF
      END IF
C
      IF(SVB(LN).NE.0) THEN
       DYYNORT=DYYTCA(LN)
       IF(ABS(DYXNORT).GT.DMAX) THEN
        DYXNORT=0.
       ELSE
        WTY=WTY+1
       END IF
      END IF
C
      DYXTCA(L)=(DYXWEST+DYXEAST+DYXSOUT+DYXNORT)/WTY
C
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      DYYTMP=DYYTCA(L)
      IF(DYYTMP.LT.0.OR.DYYTMP.GT.DMAX) THEN
C
      WRITE(88,8884)IL(L),JL(L),DYYTMP
      DYYTCA(L)=0.
      LN=LNC(L)
      LS=LSC(L)
      WTY=0.
      DYYWEST=0.
      DYYEAST=0.
      DYYSOUT=0.
      DYYNORT=0.
C
      IF(SUB(L).NE.0) THEN
       DYYWEST=DYYTCA(L-1)
       IF(DYYWEST.LT.0.OR.DYYWEST.GT.DMAX) THEN
        DYYWEST=0.
       ELSE
        WTY=WTY+1
       END IF
      END IF
C
      IF(SUB(L+1).NE.0) THEN
       DYYEAST=DYYTCA(L+1)
       IF(DYYEAST.LT.0.OR.DYYEAST.GT.DMAX) THEN
        DYYEAST=0.
       ELSE
        WTY=WTY+1
       END IF
      END IF
C
      IF(SVB(L).NE.0) THEN
       DYYSOUT=DYYTCA(LS)
       IF(DYYSOUT.LT.0.OR.DYYSOUT.GT.DMAX) THEN
        DYYSOUT=0.
       ELSE
        WTY=WTY+1
       END IF
      END IF
C
      IF(SVB(LN).NE.0) THEN
       DYYNORT=DYYTCA(LN)
       IF(DYYNORT.LT.0.OR.DYYNORT.GT.DMAX) THEN
        DYYNORT=0.
       ELSE
        WTY=WTY+1
       END IF
      END IF
C
      DYYTCA(L)=(DYYWEST+DYYEAST+DYYSOUT+DYYNORT)/WTY
C
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      CLOSE (88)
C
 8881 FORMAT(1X,'I=',I5,2X,'J=',I5,2X,'DXX= ',E12.4)
 8882 FORMAT(1X,'I=',I5,2X,'J=',I5,2X,'DXY= ',E12.4)
 8883 FORMAT(1X,'I=',I5,2X,'J=',I5,2X,'DYX= ',E12.4)
 8884 FORMAT(1X,'I=',I5,2X,'J=',I5,2X,'DYY= ',E12.4)
C
C**********************************************************************C
C
C **  WRITE OUTPUT FILES
C
      OPEN(88,FILE='disten.out',STATUS='UNKNOWN')
      CLOSE(88,STATUS='DELETE')
      OPEN(88,FILE='disten.out',STATUS='UNKNOWN')
      WRITE(88,881)
C
      DO L=2,LA
      WRITE(88,2011)IL(L),JL(L),DLON(L),DLAT(L),DXXTCA(L),
     $                    DXYTCA(L),DYXTCA(L),DYYTCA(L)
      END DO
C
      CLOSE (88)
C
      OPEN(88,FILE='uvtsc.out',STATUS='UNKNOWN')
      CLOSE(88,STATUS='DELETE')
      OPEN(88,FILE='uvtsc.out',STATUS='UNKNOWN')
      WRITE(88,882)
C
      DO L=2,LA
      AMCPT=AMCP(L)*GI
      AMSPT=AMSP(L)*GI
      WRITE(88,2012)IL(L),JL(L),DLON(L),DLAT(L),AMCPT,AMSPT,
     $                    AMCUE(L),AMSUE(L),AMCVE(L),AMSVE(L)
      END DO
C
      CLOSE (88)
C
      OPEN(88,FILE='uverv.out',STATUS='UNKNOWN')
      CLOSE(88,STATUS='DELETE')
      OPEN(88,FILE='uverv.out',STATUS='UNKNOWN')
      WRITE(88,883)
C
      DO L=2,LA
      WRITE(88,2012)IL(L),JL(L),DLON(L),DLAT(L),HLPF(L),UELPF(L),
     $                    VELPF(L),SALLPF(L,1),SALLPF(L,KC)
      END DO
C
      CLOSE (88)
C
  881 FORMAT(3X,'I',3X,'J',3X,'LON',9X,'LAT',9X,'DXX',10X,'DXY',10X,
     $          'DYX',10X,'DYY')
  882 FORMAT(3X,'I',3X,'J',3X,'LON',9X,'LAT',9X,'AMCPT',8X,'AMSPT',8X,
     $          'AMCUE',8X,'AMSUE',8X,'AMCVE',8X,'AMSVE')
  883 FORMAT(3X,'I',3X,'J',3X,'LON',9X,'LAT',9X,'HLPF',9X,'UELPF',8X,
     $          'VELPF',8X,'SALLPFBOT',4X,'SALLPFSURF')
 2011 FORMAT(2I4,2X,F10.6,2X,F10.6,4(2X,E12.4))
 2012 FORMAT(2I4,2X,F10.6,2X,F10.6,6(2X,E12.4))
C
C**********************************************************************C
C
      RETURN
      END
