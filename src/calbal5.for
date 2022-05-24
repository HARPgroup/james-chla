C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALBAL5
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINES CALBAL CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, 
C **  AND ENERGY BALANCES
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
C **  CHECK FOR END OF BALANCE PERIOD
C
      IF (NBAL.EQ.NTSMMT) THEN
cjh      WRITE(6,6666)N,NBAL
cjh 6666 FORMAT(' ACTIVE CALL TO CALBAL5, N,NBUD = ',2I5)
C
C**********************************************************************C
C
C **  CALCULATE ENDING VOLUME, SALT MASS, DYE MASS, MOMENTUM, KINETIC 
C **  ENERGY AND POTENTIAL ENERGY, AND ASSOCIATED FLUXES
C
C----------------------------------------------------------------------C
C
      VOLEND=0.
      SALEND=0.
      DYEEND=0.
      UMOEND=0.
      VMOEND=0.
      UUEEND=0.
      VVEEND=0.
      PPEEND=0.
      BBEEND=0.
C
      DO L=2,LA
      LN=LNC(L)
      VOLEND=VOLEND+SPB(L)*DXYP(L)*HP(L)
      UMOEND=UMOEND+SPB(L)*0.5*DXYP(L)*HP(L)*(DYIU(L)*HUI(L)*UHDYE(L)
     $                                 +DYIU(L+1)*HUI(L+1)*UHDYE(L+1))
      VMOEND=VMOEND+SPB(L)*0.5*DXYP(L)*HP(L)*(DXIV(L)*HVI(L)*VHDXE(L)
     $                                 +DXIV(LN)*HVI(LN)*VHDXE(LN))
      PPEEND=PPEEND+SPB(L)*0.5*DXYP(L)
     $             *(GI*P(L)*P(L)-G*BELV(L)*BELV(L))
      END DO
C
      AMOEND=SQRT(UMOEND*UMOEND+VMOEND*VMOEND)
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      SALEND=SALEND+SCB(L)*DXYP(L)*HP(L)*SAL(L,K)*DZC(K)
      DYEEND=DYEEND+SCB(L)*DXYP(L)*HP(L)*DYE(L,K)*DZC(K)
C     UUEEND=UUEEND+SPB(L)*0.25*(DXYU(L)*HU(L)*U(L,K)*U(L,K)
C    $      +DXYU(L+1)*HU(L+1)*U(L+1,K)*U(L+1,K))*DZC(K)
C     VVEEND=VVEEND+SPB(L)*0.25*(DXYV(L)*HV(L)*V(L,K)*V(L,K)
C    $      +DXYV(LN)*HV(LN)*V(LN,K)*V(LN,K))*DZC(K)
      UUEEND=UUEEND+SPB(L)*0.125*DXYP(L)*HP(L)*DZC(K)
     $      *( (U(L,K)+U(L+1,K))*(U(L,K)+U(L+1,K)) )
      VVEEND=VVEEND+SPB(L)*0.125*DXYP(L)*HP(L)*DZC(K)
     $      *( (V(L,K)+V(LN,K))*(V(L,K)+V(LN,K)) )
      BBEEND=BBEEND+SPB(L)*GP*DXYP(L)*HP(L)*DZC(K)*( BELV(L) 
     $      +0.5*HP(L)*(Z(K)+Z(K-1)) )*B(L,K)
      END DO
      END DO
C
      UUEOUT=DT*UUEOUT
      VVEOUT=DT*VVEOUT
      PPEOUT=DT*PPEOUT
      BBEOUT=DT*BBEOUT
      VOLOUT=DT*VOLOUT
      SALOUT=DT*SALOUT
      DYEOUT=DT*DYEOUT
      UMOOUT=DT*UMOOUT
      VMOOUT=DT*VMOOUT
C
      ENEBEG=UUEBEG+VVEBEG+PPEBEG+BBEBEG
      ENEEND=UUEEND+VVEEND+PPEEND+BBEEND
      ENEOUT=UUEOUT+VVEOUT+PPEOUT+BBEOUT
C
      VOLBMO=VOLBEG-VOLOUT
      SALBMO=SALBEG-SALOUT
      DYEBMO=DYEBEG-DYEOUT
      UMOBMO=UMOBEG-DYEOUT
      VMOBMO=VMOBEG-DYEOUT
      ENEBMO=ENEBEG-ENEOUT
C
      VOLERR=VOLEND-VOLBMO
      SALERR=SALEND-SALBMO
      DYEERR=DYEEND-DYEBMO
      UMOERR=UMOEND-UMOBMO
      VMOERR=VMOEND-VMOBMO
      ENEERR=ENEEND-ENEBMO
C
      RVERDE=-9999.
      RSERDE=-9999.
      RDERDE=-9999.
      RUERDE=-9999.
      RVERDE=-9999.
      REERDE=-9999.
C
      RVERDO=-9999.
      RSERDO=-9999.
      RDERDO=-9999.
      RUERDO=-9999.
      RVERDO=-9999.
      REERDO=-9999.
C
      IF(VOLEND.NE.0.) RVERDE=VOLERR/VOLEND
      IF(SALEND.NE.0.) RSERDE=SALERR/SALEND
      IF(DYEEND.NE.0.) RDERDE=DYEERR/DYEEND
      IF(UMOEND.NE.0.) RUMERDE=UMOERR/UMOEND
      IF(VMOEND.NE.0.) RVMERDE=VMOERR/VMOEND
      IF(ENEEND.NE.0.) REERDE=ENEERR/ENEEND
C
      IF(VOLOUT.NE.0.) RVERDO=VOLERR/VOLOUT
      IF(SALOUT.NE.0.) RSERDO=SALERR/SALOUT
      IF(DYEOUT.NE.0.) RDERDO=DYEERR/DYEOUT
      IF(UMOOUT.NE.0.) RUMERDO=UMOERR/UMOOUT
      IF(VMOOUT.NE.0.) RVMERDO=VMOERR/VMOOUT
      IF(ENEOUT.NE.0.) REERDO=ENEERR/ENEOUT

C**********************************************************************C
C
C **  OUTPUT BALANCE RESULTS TO FILE bal.out
C
C----------------------------------------------------------------------C
C
      IF (JSBAL.EQ.1) THEN
        OPEN(89,FILE='bal.out',STATUS='UNKNOWN')
        CLOSE(89,STATUS='DELETE')
        OPEN(89,FILE='bal.out',STATUS='UNKNOWN')
        JSBAL=0
       ELSE
        OPEN(89,FILE='bal.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
C
      WRITE(89,890)NTSMMT,N
      WRITE(89,891)
      WRITE(89,892)VOLBEG,SALBEG,DYEBEG,ENEBEG,UMOBEG,VMOBEG,AMOBEG
      WRITE(89,900)
      WRITE(89,893)
      WRITE(89,892)VOLOUT,SALOUT,DYEOUT,ENEOUT,UMOOUT,VMOOUT
      WRITE(89,900)
      WRITE(89,894)
      WRITE(89,892)VOLBMO,SALBMO,DYEBMO,ENEBMO,UMOBMO,VMOBMO
      WRITE(89,900)
      WRITE(89,895)
      WRITE(89,892)VOLEND,SALEND,DYEEND,ENEEND,UMOEND,VMOEND,AMOEND
      WRITE(89,900)
      WRITE(89,896)
      WRITE(89,892)VOLERR,SALERR,DYEERR,ENEERR,UMOERR,VMOERR
      WRITE(89,900)
      WRITE(89,897)
      WRITE(89,892)RVERDE,RSERDE,RDERDE,REERDE,RUMERDE,RVMERDE
      WRITE(89,900)
      WRITE(89,898)
      WRITE(89,892)RVERDO,RSERDO,RDERDO,REERDO,RUMERDO,RVMERDO
      WRITE(89,899)
      UUEBMO=UUEBEG-UUEOUT
      VVEBMO=VVEBEG-VVEOUT
      PPEBMO=PPEBEG-PPEOUT
      BBEBMO=BBEBEG-BBEOUT
      WRITE(89,901)UUEBEG
      WRITE(89,902)UUEOUT
      WRITE(89,903)UUEBMO
      WRITE(89,904)UUEEND
      WRITE(89,900)
      WRITE(89,905)VVEBEG
      WRITE(89,906)VVEOUT
      WRITE(89,907)VVEBMO
      WRITE(89,908)VVEEND
      WRITE(89,900)
      WRITE(89,909)PPEBEG
      WRITE(89,910)PPEOUT
      WRITE(89,911)PPEBMO
      WRITE(89,912)PPEEND
      WRITE(89,900)
      WRITE(89,913)BBEBEG
      WRITE(89,914)BBEOUT
      WRITE(89,915)BBEBMO
      WRITE(89,916)BBEEND
      WRITE(89,900)
      WRITE(89,899)
C
      CLOSE(89)
C
  890 FORMAT (' VOLUME, MASS, AND ENERGY BALANCE OVER',I5,' TIME STEPS'
     $,' ENDING AT TIME STEP',I5,//)
  891 FORMAT (' INITIAL VOLUME    INITIAL SALT    INITIAL DYE     '
     $,'INITIAL ENER    INITIAL UMO     INITIAL VMO     '
     $,'INITIAL AMO',/)
  892 FORMAT (1X,7(E14.6,2X))
  893 FORMAT (' VOLUME OUT        SALT OUT        DYE OUT         '
     $,'ENERGY OUT      UMO OUT         VMO OUT',/)
  894 FORMAT (' INITIAL-OUT VOL   INIT-OUT SALT   INIT-OUT DYE    '
     $,'INIT-OUT ENER   INIT-OUT UMO    INIT-OUT VMO',/)
  895 FORMAT (' FINAL VOLUME      FINAL SALT      FINAL DYE       '
     $,'FINAL ENERGY    FINAL UMO       FINAL VMO       '
     $,'FINAL AMO',/)
  896 FORMAT (' VOLUME ERR        SALT ERR        DYE ERR         '
     $,'ENERGY ERR      UMO ERR         VMO ERR',/)
  897 FORMAT (' R VOL/END ER      R SAL/END ER    R DYE/END ER    '
     $,'R ENE/END ER    R UMO/END ER    R VMO/END ER',/)
  898 FORMAT (' R VOL/OUT ER      R SAL/OUT ER    R DYE/OUT ER    '
     $,'R ENE/OUT ER    R UMO/OUT ER    R VMO/OUT ER',/)
  899 FORMAT (////)
  900 FORMAT (//)
  901 FORMAT(' UUEBEG =  ',E14.6)
  902 FORMAT(' UUEOUT =  ',E14.6)
  903 FORMAT(' UUEBMO =  ',E14.6)
  904 FORMAT(' UUEEND =  ',E14.6)
  905 FORMAT(' VVEBEG =  ',E14.6)
  906 FORMAT(' VVEOUT =  ',E14.6)
  907 FORMAT(' VVEBMO =  ',E14.6)
  908 FORMAT(' VVEEND =  ',E14.6)
  909 FORMAT(' PPEBEG =  ',E14.6)
  910 FORMAT(' PPEOUT =  ',E14.6)
  911 FORMAT(' PPEBMO =  ',E14.6)
  912 FORMAT(' PPEEND =  ',E14.6)
  913 FORMAT(' BBEBEG =  ',E14.6)
  914 FORMAT(' BBEOUT =  ',E14.6)
  915 FORMAT(' BBEBMO =  ',E14.6)
  916 FORMAT(' BBEEND =  ',E14.6)
C
C**********************************************************************C
C
C     RESET COUNTER
C
      NBAL=0
C
      END IF 
C
      NBAL=NBAL+1
C
C**********************************************************************C
C
      RETURN
      END
