C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALTSXY
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
C
C**********************************************************************C
C
C ** SUBROUTINE CALTSXY UPDATES TIME VARIABLE SURFACE WIND STRESS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
      DIMENSION  PATMTT(NASERM),TATMTT(NASERM),TWETTT(NASERM),
     $           RAINTT(NASERM),EVAPTT(NASERM),SOLSWRTT(NASERM),
     $           CLOUDTT(NASERM),SVPAT(NASERM),VPAT(NASERM),
     $           RHAT(NASERM),CLEVAPT(NASERM),CCNHTTT(NASERM),
     $           WINDE(NWSERM),WINDN(NWSERM)
C
C**********************************************************************C
C
      IF(NWSER.GT.0) THEN
C
      DO NA=1,NWSER
      TIME=DT*FLOAT(N)/TCWSER(NA)+TBEGIN*(TCON/TCWSER(NA))
C
      M1=MWTLAST(NA)
  200 CONTINUE
      M2=M1+1
      IF (TIME.GT.TWSER(M2,NA)) THEN
        M1=M2
        GO TO 200
       ELSE
        MWTLAST(NA)=M1
      END IF
C
      TDIFF=TWSER(M2,NA)-TWSER(M1,NA)
      WTM1=(TWSER(M2,NA)-TIME)/TDIFF
      WTM2=(TIME-TWSER(M1,NA))/TDIFF

      DEGM1=90.-WINDD(M1,NA)
      DEGM2=90.-WINDD(M2,NA)
      WINDS1=WTM1*WINDS(M1,NA)+WTM2*WINDS(M2,NA)
      WINDS2=WTM1*WINDS(M1,NA)+WTM2*WINDS(M2,NA)
      WINDE1=WINDS(M1,NA)*COS(DEGM1/57.29578)
      WINDN1=WINDS(M1,NA)*SIN(DEGM1/57.29578)
      WINDE2=WINDS(M2,NA)*COS(DEGM2/57.29578)
      WINDN2=WINDS(M2,NA)*SIN(DEGM2/57.29578)
      WINDE(NA)=WTM1*WINDE1+WTM2*WINDE2
      WINDN(NA)=WTM1*WINDN1+WTM2*WINDN2
C
      END DO
C
      END IF
C
C**********************************************************************C
C
      IF(NWSER.GT.0) THEN
C
      DO L=2,LA
       WNDVELE(L)=0.
       WNDVELN(L)=0.
      END DO
C
      DO NA=1,NWSER
      DO L=2,LA
       WNDVELE(L)=WNDVELE(L)+WNDWHT(L,NA)*WINDE(NA)
       WNDVELN(L)=WNDVELN(L)+WNDWHT(L,NA)*WINDN(NA)
      END DO
      END DO
      if(n.eq.1)then
      open(1001,file='wind.tmp',status='unknown')
      close(1001,status='delete')
      open(1001,file='wind.tmp',status='unknown')
      write(1001,*)'       L,Time,winde,windn'
      close(1001)
      endif
      
      if(mod(n,180).eq.0)then
      open(1001,file='wind.tmp',status='unknown',access='append') 
      TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
      DO L=2,2
      write(1001,601)n,TIME,WNDVELE(L),WNDVELN(L)
      enddo !l
      close(1001)
      endif
601   format(I10,',',f12.4,',',f12.3,',',f12.3)      
C
      DO L=2,LA
       WINDST(L)=SQRT( WNDVELE(L)*WNDVELE(L)
     $                 +WNDVELN(L)*WNDVELN(L) )
C     	  IF(HP(L).LT.2.0.and.WINDST(L).GT.6.0 ) WINDST(L)=WINDST(L)*0.7
C     	  IF(HP(L).LT.1.0) WINDST(L)=WINDST(L)*HP(L)
        IF(HP(L).LT.3.0) WINDST(L)=WINDST(L)*HP(L)/3.0
        IF(HP(L).LT.1.2) WINDST(L)=WINDST(L)*0.1
      END DO
C
      DO L=2,LA
        TSEAST=1.2E-6*(0.8+0.065*WINDST(L))*WINDST(L)*WNDVELE(L)
        TSNORT=1.2E-6*(0.8+0.065*WINDST(L))*WINDST(L)*WNDVELN(L)
        DETTMP=1./( CUE(L)*CVN(L)-CUN(L)*CVE(L) )
        TSX(L)=DETTMP*(CVN(L)*TSEAST-CVE(L)*TSNORT)
        TSY(L)=DETTMP*(CUE(L)*TSNORT-CUN(L)*TSEAST)
      END DO
C
      END IF

C**********************************************************************C
C
      IF(NASER.GT.0) THEN
C
      DO NA=1,NASER
      TIME=DT*FLOAT(N)/TCASER(NA)+TBEGIN*(TCON/TCASER(NA))
C
      M1=MATLAST(NA)
  100 CONTINUE
      M2=M1+1
      IF (TIME.GT.TASER(M2,NA)) THEN
        M1=M2
        GO TO 100
       ELSE
        MATLAST(NA)=M1
       END IF
C
      TDIFF=TASER(M2,NA)-TASER(M1,NA)
      WTM1=(TASER(M2,NA)-TIME)/TDIFF
      WTM2=(TIME-TASER(M1,NA))/TDIFF
C
      PATMTT(NA)=WTM1*PATM(M1,NA)+WTM2*PATM(M2,NA)
      TATMTT(NA)=WTM1*TDRY(M1,NA)+WTM2*TDRY(M2,NA)
      TWETTT(NA)=WTM1*TWET(M1,NA)+WTM2*TWET(M2,NA)
      RAINTT(NA)=WTM1*RAIN(M1,NA)+WTM2*RAIN(M2,NA)
      EVAPTT(NA)=WTM1*EVAP(M1,NA)+WTM2*EVAP(M2,NA)
      SOLSWRTT(NA)=WTM1*SOLSWR(M1,NA)+WTM2*SOLSWR(M2,NA)
C Rico modified, by time coeffient to adjust the SOLAR SW RADIATION      
C      print*, 'CSERT(1,NCSER(2),2)=',CSERT(1,NCSER(2),2)
C      print*, 'NCSER(2)=', NCSER(2)
      SOLSWRTT(NA)=(WTM1*SOLSWR(M1,NA)+WTM2*SOLSWR(M2,NA))
 !    $ *CSERT(1,NCSER(2),2)             no adj J.S.
C Rico modified, by time coeffient to adjust the SOLAR SW RADIATION      
      CLOUDTT(NA)=WTM1*CLOUD(M1,NA)+WTM2*CLOUD(M2,NA)
      SVPAT(NA)=
     $  (10.**((0.7859+0.03477*TATMTT(NA))/(1.+0.00412*TATMTT(NA))))
C    $ *(1+1.E-6*PATMTT(NA)*(4.5+0.0006*TATMTT(NA)*TATMT(NA)))
      IF(IRELH(NA).EQ.0) THEN
        RHAT(NA)=1.
     $         -0.00066*(PATMTT(NA)/SVPAT(NA))*(TATMTT(NA)-TWETTT(NA))
       ELSE
        RHAT(NA)=TWETTT(NA)
      END IF
      VPAT(NA)=RHAT(NA)*SVPAT(NA)
      CLEVAPT(NA)=1.E-3
      CCNHTTT(NA)=1.E-3
      IF(REVC.GT.1.E-6) CLEVAPT(NA)=REVC*1.E-3
      IF(RCHC.GT.1.E-6) CCNHTTT(NA)=RCHC*1.E-3
C
      END DO
C
      END IF
C
C**********************************************************************C
C
      IF(NASER.GT.0) THEN
C
      DO L=2,LA
       PATMT(L)=0.
       TATMT(L)=0.
       RAINT(L)=0.
       EVAPT(L)=0.
       SOLSWRT(L)=0.
       CLOUDT(L)=0.
       SVPA(L)=0.
       RHA(L)=0.
       VPA(L)=0.
       CLEVAP(L)=0.
       CCNHTT(L)=0.
      END DO
C
      DO NA=1,NASER
      DO L=2,LA
       PATMT(L)=PATMT(L)+ATMWHT(L,NA)*PATMTT(NA)
       TATMT(L)=TATMT(L)+ATMWHT(L,NA)*TATMTT(NA)
       RAINT(L)=RAINT(L)+ATMWHT(L,NA)*RAINTT(NA)
       EVAPT(L)=EVAPT(L)+ATMWHT(L,NA)*EVAPTT(NA)
       SOLSWRT(L)=SOLSWRT(L)+ATMWHT(L,NA)*SOLSWRTT(NA)
       CLOUDT(L)=CLOUDT(L)+ATMWHT(L,NA)*CLOUDTT(NA)
       SVPA(L)=SVPA(L)+ATMWHT(L,NA)*SVPAT(NA)
       RHA(L)=RHA(L)+ATMWHT(L,NA)*RHAT(NA)
       VPA(L)=VPA(L)+ATMWHT(L,NA)*VPAT(NA)
       CLEVAP(L)=CLEVAP(L)+ATMWHT(L,NA)*CLEVAPT(NA)
       CCNHTT(L)=CCNHTT(L)+ATMWHT(L,NA)*CCNHTTT(NA)
      END DO
      END DO
c
c     RAINT(45)=0.0 ! Hardwired for Wister Reservoir (5,10), John & Ji, 7/30/99
C
      IF(REVC.LT.0.) THEN
        DO L=2,LA
          CLEVAP(L)=1.E-3*(0.8+0.065*WINDST(L))
        END DO
      END IF
C
      IF(RCHC.LT.0.) THEN
        DO L=2,LA
          CCNHTT(L)=1.E-3*(0.8+0.065*WINDST(L))
        END DO
      END IF
C
      END IF
C
C**********************************************************************C
C
      RETURN
      END
