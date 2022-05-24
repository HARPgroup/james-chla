C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE LTMT
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE LTMT EXECUTES A LONG-TERM MASS TRANSPORT
C **  TIME INTEGRATION 
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      OPEN(99,FILE='restran.inp',STATUS='UNKNOWN')
C
C**********************************************************************C
C
C **  READ FIRST MEAN MASS TRANSPORT FIELD AND INITIALIZED NEXT
C
C----------------------------------------------------------------------C
C
      CALL RESTRAN
C
      DO L=2,LA
      H1P(L)=HLPF(L)
      HP(L)=HLPF(L)
      UHDY1E(L)=0.
      UHDYE(L)=0.
      VHDX1E(L)=0.
      VHDXE(L)=0.
      END DO
C
c     IF(ISLTMT.EQ.1) THEN
        DO K=1,KS
        DO L=2,LA
        AB(L,K)=ABLPF(L,K)
        END DO
        END DO
c      ELSE
c       DO K=1,KS
c       DO L=2,LA
c       AB(L,K)=ABEFF(L,K)
c       END DO
c       END DO
c     END IF
C
      DO K=1,KC
      DO L=1,LC
      UHDY1(L,K)=UHLPF(L,K)*DYU(L)
      VHDX1(L,K)=VHLPF(L,K)*DXV(L)
      END DO
      END DO
C
      DO K=1,KC
      DO L=1,LC
      UHDY(L,K)=UHDY1(L,K)
      UHDY1E(L)=UHDY1E(L)+UHDY1(L,K)*DZC(K)
      VHDX(L,K)=VHDX1(L,K)
      VHDX1E(L)=VHDX1E(L)+VHDX1(L,K)*DZC(K)
      END DO
      END DO
C
      DO L=1,LC
      UHDYE(L)=UHDY1E(L)
      VHDXE(L)=VHDX1E(L)
      END DO
C
C----------------------------------------------------------------------C
C
      IF (ISLTMTS.EQ.0) THEN
C
      DO K=1,KC
      DO L=1,LC
      SAL(L,K)=0.
      SAL1(L,K)=0.
      END DO
      END DO
C
      ELSE
C
      DO K=1,KC
      DO L=1,LC
      SAL(L,K)=SALLPF(L,K)
      SAL1(L,K)=SALLPF(L,K)
      END DO
      END DO
C
      END IF
C
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)
      SAL1(L,K)=SALLPF(L,K)
      SAL(L,K)=SALLPF(L,K)
      END DO
      END DO
C
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)      
      SAL1(L,K)=SALLPF(L,K)
      SAL(L,K)=SALLPF(L,K)
      END DO
      END DO
C
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)
      SAL1(L,K)=SALLPF(L,K)
      SAL(L,K)=SALLPF(L,K)
      END DO
      END DO
C
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)
      SAL1(L,K)=SALLPF(L,K)
      SAL(L,K)=SALLPF(L,K)
      END DO
      END DO
C
C**********************************************************************C
C
C **  READ SECOND MEAN MASS TRANSPORT FIELD FOR UNSTEADY CASE 
C
C----------------------------------------------------------------------C
C
      IF (ISSSMMT.EQ.0) THEN
C
      CALL RESTRAN
C
      RNTCMMT=FLOAT(NTSMMT)/FLOAT(NTSPTC)
      DO L=2,LA
      HP(L)=H1P(L)+(HLPF(L)-H1P(L))/(RNTCMMT*FLOAT(NTSPTC))
      UHDYE(L)=0.
      VHDXE(L)=0.
      END DO
C
c     IF(ISLTMT.EQ.1) THEN
        DO K=1,KS
        DO L=2,LA
        AB(L,K)=SQRT(AB(L,K)*ABLPF(L,K))
        END DO
        END DO
c      ELSE
c       DO K=1,KS
c       DO L=2,LA
c       AB(L,K)=SQRT(AB(L,K)*ABEFF(L,K))
c       END DO
c       END DO
c     END IF
C
      DO K=1,KC
      DO L=1,LC
      UHDY(L,K)=UHLPF(L,K)*DYU(L)
      VHDX(L,K)=VHLPF(L,K)*DXV(L)
      SAL(L,K)=SALLPF(L,K)
      END DO
      END DO
C
      DO K=1,KC
      DO L=1,LC
      UHDYE(L)=UHDYE(L)+UHDY(L,K)*DZC(K)
      VHDXE(L)=VHDXE(L)+VHDX(L,K)*DZC(K)
      END DO
      END DO
C
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      END DO
      END DO
C
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)      
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      END DO
      END DO
C
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      END DO
      END DO
C
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      END DO
      END DO
C
      END IF
C
C**********************************************************************C
C     
C **  ADJUST THE STEADY OR INITIAL MEAN MASS TRANSPORT FIELD  
C
      CALL ADJMMT
C
C**********************************************************************C
C
C **  INITIALIZE COUNTERS
C
      NCTS=0
C
C**********************************************************************C
C**********************************************************************C
C
C **  BEGIN THE TIME INTEGRATION LOOP
C
C----------------------------------------------------------------------C
C
      DO 1000 N=1,NTS
      NCTS=NCTS+1
C
C**********************************************************************C
C
C **  ADVANCE CONCENTRATION FIELDS
C
C----------------------------------------------------------------------C
C
      IF(ISCDCA(1).NE.1) THEN
        CALL CALCONC (2)
       ELSE
        CALL CALCONC (3)
      END IF
      IF (ISHOW.GE.1) CALL SHOWVAL1
C
C----------------------------------------------------------------------C
C
C **  CHECK RANGE OF SALINITY AND DYE CONCENTRATION
C
      IF (ISMMC.EQ.1) THEN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF (SAL(L,K).GT.SALMAX) THEN
       SALMAX=SAL(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      END IF
      IF (SAL(L,K).LT.SALMIN) THEN
       SALMIN=SAL(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      END IF
      END DO
      END DO
C
      WRITE(6,6001)N
      WRITE(6,6002)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6003)SALMIN,IMIN,JMIN,KMIN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF (DYE(L,K).GT.SALMAX) THEN
       SALMAX=DYE(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      END IF
      IF (DYE(L,K).LT.SALMIN) THEN
       SALMIN=DYE(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      END IF
      END DO
      END DO
C
      WRITE(6,6004)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6005)SALMIN,IMIN,JMIN,KMIN
C
      END IF
C
 6001 FORMAT(1X,'N=',I10)
 6002 FORMAT(1X,'SALMAX=',F14.4,5X,'I,J,K=',(3I10))
 6003 FORMAT(1X,'SALMIN=',F14.4,5X,'I,J,K=',(3I10))
 6004 FORMAT(1X,'DYEMAX=',F14.4,5X,'I,J,K=',(3I10))
 6005 FORMAT(1X,'DYEMIN=',F14.4,5X,'I,J,K=',(3I10))
C
C**********************************************************************C
C
C **  ADVANCE MMT FIELDS
C
C----------------------------------------------------------------------C
C
      IF (ISSSMMT.EQ.0) THEN
C
      IF (NCTS.EQ.NTSMMT) THEN
      NCTS=0
C
      DO L=2,LA
      H1P(L)=HLPF(L)
      UHDY1E(L)=UHDYE(L)
      VHDX1E(L)=VHDXE(L)
      END DO
C
      DO K=1,KC
      DO L=2,LA
      SAL1(L,K)=SALLPF(L,K)
      UHDY1(L,K)=UHDY(L,K)
      VHDX1(L,K)=VHDX(L,K)
      END DO
      END DO
C
      CALL RESTRAN
C
      RNTCMMT=FLOAT(NTSMMT)/FLOAT(NTSPTC)
      DO L=2,LA
      HP(L)=H1P(L)+(HLPF(L)-H1P(L))/(RNTCMMT*FLOAT(NTSPTC))
      UHDYE(L)=0.
      VHDXE(L)=0.
      END DO
C
c     IF(ISLTMT.EQ.1) THEN
        DO K=1,KS
        DO L=2,LA
        AB(L,K)=SQRT(AB(L,K)*ABLPF(L,K))
        END DO
        END DO
c      ELSE
c       DO K=1,KS
c       DO L=2,LA
c       AB(L,K)=SQRT(AB(L,K)*ABEFF(L,K))
c       END DO
c       END DO
c     END IF
C
      DO K=1,KC
      DO L=1,LC
      UHDY(L,K)=UHLPF(L,K)*DYU(L)
      VHDX(L,K)=VHLPF(L,K)*DXV(L)
      END DO
      END DO
C
      DO K=1,KC
      DO L=1,LC
      UHDYE(L)=UHDYE(L)+UHDY(L,K)*DZC(K)
      VHDXE(L)=VHDXE(L)+VHDX(L,K)*DZC(K)
      END DO
      END DO
C
      CALL ADJMMT
C
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)
     $                   -SAL1(L,K))/(RNTCMMT*FLOAT(NTSPTC))
      END DO
      END DO
C
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)      
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      END DO
      END DO
C
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      END DO
      END DO
C
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))/FLOAT(NTSMMT)
      END DO
      END DO
C
      ELSE
C
      DO L=2,LA
      HP(L)=H1P(L)+(HLPF(L)-H1P(L))*FLOAT(NCTS+1)/FLOAT(NTSMMT)
      END DO
C
      DO K=1,KC
      DO LL=1,NCBS
      L=LCBS(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))*FLOAT(NCTS+1)/
     $         FLOAT(NTSMMT)
      END DO
      END DO
C
      DO K=1,KC
      DO LL=1,NCBW
      L=LCBW(LL)      
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))*FLOAT(NCTS+1)/
     $         FLOAT(NTSMMT)
      END DO
      END DO
C
      DO K=1,KC
      DO LL=1,NCBE
      L=LCBE(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))*FLOAT(NCTS+1)/
     $         FLOAT(NTSMMT)
      END DO
      END DO
C
      DO K=1,KC
      DO LL=1,NCBN
      L=LCBN(LL)
      SAL(L,K)=SAL1(L,K)+(SALLPF(L,K)-SAL1(L,K))*FLOAT(NCTS+1)/
     $         FLOAT(NTSMMT)
      END DO
      END DO
C
      END IF
C
      END IF
C
C**********************************************************************C
C
 1000 CONTINUE
C
C**********************************************************************C
C**********************************************************************C
C
C **  WRITE GRAPHICS FILES FOR RESIDUAL VARIABLES
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=1,LC
      SALLPF(L,K)=SAL(L,K)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  RESIDUAL SALINITY CONTOURING IN HORIZONTAL: SUBROUTINE RSALPLTH
C
      IF (ISRSPH(1).EQ.1.AND.ISTRAN(1).GE.1) THEN
       CALL RSALPLTH(1,SALLPF)
      END IF
C
      IF (ISRSPH(2).EQ.1.AND.ISTRAN(2).GE.1) THEN
       CALL RSALPLTH(2,TEMLPF)
      END IF
C
      IF (ISRSPH(3).EQ.1.AND.ISTRAN(3).GE.1) THEN
       CALL RSALPLTH(3,DYELPF)
      END IF
C
      IF (ISRSPH(4).EQ.1.AND.ISTRAN(4).GE.1) THEN
       CALL RSALPLTH(4,SEDTLPF)
      END IF
C
      IF (ISRSPH(4).EQ.1.AND.ISTRAN(5).GE.1) THEN
       CALL RSALPLTH(5,SNDTLPF)
      END IF
C
CFIX      IF (ISRSPH(5).EQ.1.AND.ISTRAN(6).GE.1) THEN
CFIX       CALL RSALPLTH(6,TOXLPF)
CFIX      END IF
C
      IF (ISRSPH(7).EQ.1.AND.ISTRAN(7).GE.1) THEN
       CALL RSALPLTH(7,SFLLPF)
      END IF
C
C----------------------------------------------------------------------C
C
C **  RESIDUAL VELOCITY VECTOR PLOTTING IN HORIZONTAL PLANES:
C **  SUBROUTINE RVELPLTH
C
      IF (ISRVPH.EQ.1) THEN
       CALL RVELPLTH
      END IF
C
C----------------------------------------------------------------------C
C
C **  RESIDUAL SURFACE ELEVATION PLOTTING IN HORIZONTAL PLANES:
C **  SUBROUTINE RVELPLTH
C
      IF (ISRPPH.EQ.1) THEN
       CALL RSURFPLT
      END IF
C
C----------------------------------------------------------------------C
C
C **  RESIDUAL SALINITY AND VERTICAL MASS DIFFUSIVITY CONTOURING IN 
C **  3 VERTICAL PLANES:  SUBROUTINE RSALPLTV
C
      IF (ISRSPV(1).GE.1) THEN
       CALL RSALPLTV(1)
      END IF
C
C----------------------------------------------------------------------C
C
C **  RESIDUAL NORMAL AND TANGENTIAL VELOCITY CONTOURING AND AND 
C **  TANGENTIAL VELOCITY VECTOR PLOTTING IN VERTICAL PLANES:
C **  SUBROUTINE RVELPLTV
C
      IF (ISRVPV.GE.1) THEN
       CALL RVELPLTV
      END IF
C
C**********************************************************************C
C
      CLOSE(99) 
C
C**********************************************************************C
C
      RETURN
      END
