C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE ROUT3D
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C     DIMENSION APT(KPCM,KCM), AP(KPCM), IAP(KPCM)
      DIMENSION  AKL(KPCM,LCM),  AIJ(IGM,JGM)
C
C     CHARACTER BLANK,LET(0:50),CHARY(ICM,JCM,KPCM)
C     CHARACTER *80 TITLE
      CHARACTER *12 SALFN,TEMFN,DYEFN,SEDFN,UUUFN,VVVFN,WWWFN,
     $              CMPFN,SNDFN,TOXFN
C
C     DATA BLANK/' '/
C
C     DATA LET/' ','y','Y','x','X','w','W','v','V','u','U','t','T',
C    $         's','S','r','R','q','Q','p','P','o','O','n','N',
C    $         'm','M','l','L','k','K','j','J','i','I','h','H',
C    $         'g','G','f','F','e','E','d','D','c','C','b','B',
C    $         'a','A'/    
C
C**********************************************************************C
C
C **  INITIALIZE OUTPUT FILES
C
C----------------------------------------------------------------------C
C
      IAD=I3DMAX-I3DMIN+1
      JAD=J3DMAX-J3DMIN+1
C
C     NEED TO REPLACE .asc with (KPCxJADxIAD) IN FILE NAMES BELOW
C
      NRCAL3D=NRCAL3D+1
C
      IF(NRCAL3D.EQ.1) THEN
        SALFN='rsal3d01.asc'
        TEMFN='rtem3d01.asc'
        DYEFN='rdye3d01.asc'
        SEDFN='rsed3d01.asc'
        SNDFN='rsnd3d01.asc'
        TOXFN='rtox3d01.asc'
        UUUFN='ruuu3d01.asc'
        VVVFN='rvvv3d01.asc'
        WWWFN='rwww3d01.asc'
        CMPFN='rcmp3d01.asc'
      END IF
      IF(NRCAL3D.EQ.2) THEN
        SALFN='rsal3d02.asc'
        TEMFN='rtem3d02.asc'
        DYEFN='rdye3d02.asc'
        SEDFN='rsed3d02.asc'
        SNDFN='rsnd3d02.asc'
        TOXFN='rtox3d02.asc'
        UUUFN='ruuu3d02.asc'
        VVVFN='rvvv3d02.asc'
        WWWFN='rwww3d02.asc'
        CMPFN='rcmp3d02.asc'
      END IF
      IF(NRCAL3D.EQ.3) THEN
        SALFN='rsal3d03.asc'
        TEMFN='rtem3d03.asc'
        DYEFN='rdye3d03.asc'
        SEDFN='rsed3d03.asc'
        SNDFN='rsnd3d03.asc'
        TOXFN='rtox3d03.asc'
        UUUFN='ruuu3d03.asc'
        VVVFN='rvvv3d03.asc'
        WWWFN='rwww3d03.asc'
        CMPFN='rcmp3d03.asc'
      END IF
      IF(NRCAL3D.EQ.4) THEN
        SALFN='rsal3d04.asc'
        TEMFN='rtem3d04.asc'
        DYEFN='rdye3d04.asc'
        SEDFN='rsed3d04.asc'
        SNDFN='rsnd3d04.asc'
        TOXFN='rtox3d04.asc'
        UUUFN='ruuu3d04.asc'
        VVVFN='rvvv3d04.asc'
        WWWFN='rwww3d04.asc'
        CMPFN='rcmp3d04.asc'
       END IF
       IF(NRCAL3D.EQ.5) THEN
        SALFN='rsal3d05.asc'
        TEMFN='rtem3d05.asc'
        DYEFN='rdye3d05.asc'
        SEDFN='rsed3d05.asc'
        SNDFN='rsnd3d05.asc'
        TOXFN='rtox3d05.asc'
        UUUFN='ruuu3d05.asc'
        VVVFN='rvvv3d05.asc'
        WWWFN='rwww3d05.asc'
        CMPFN='rcmp3d05.asc'
      END IF
      IF(NRCAL3D.EQ.6) THEN
        SALFN='rsal3d06.asc'
        TEMFN='rtem3d06.asc'
        DYEFN='rdye3d06.asc'
        SEDFN='rsed3d06.asc'
        SNDFN='rsnd3d06.asc'
        TOXFN='rtox3d06.asc'
        UUUFN='ruuu3d06.asc'
        VVVFN='rvvv3d06.asc'
        WWWFN='rwww3d06.asc'
        CMPFN='rcmp3d06.asc'
      END IF
      IF(NRCAL3D.EQ.7) THEN
        SALFN='rsal3d07.asc'
        TEMFN='rtem3d07.asc'
        DYEFN='rdye3d07.asc'
        SEDFN='rsed3d07.asc'
        SNDFN='rsnd3d07.asc'
        TOXFN='rtox3d07.asc'
        UUUFN='ruuu3d07.asc'
        VVVFN='rvvv3d07.asc'
        WWWFN='rwww3d07.asc'
        CMPFN='rcmp3d07.asc'
      END IF
      IF(NRCAL3D.EQ.8) THEN
        SALFN='rsal3d08.asc'
        TEMFN='rtem3d08.asc'
        DYEFN='rdye3d08.asc'
        SEDFN='rsed3d08.asc'
        SNDFN='rsnd3d08.asc'
        TOXFN='rtox3d08.asc'
        UUUFN='ruuu3d08.asc'
        VVVFN='rvvv3d08.asc'
        WWWFN='rwww3d08.asc'
        CMPFN='rcmp3d08.asc'
      END IF
      IF(NRCAL3D.EQ.9) THEN
        SALFN='rsal3d09.asc'
        TEMFN='rtem3d09.asc'
        DYEFN='rdye3d09.asc'
        SEDFN='rsed3d09.asc'
        SNDFN='rsnd3d09.asc'
        TOXFN='rtox3d09.asc'
        UUUFN='ruuu3d09.asc'
        VVVFN='rvvv3d09.asc'
        WWWFN='rwww3d09.asc'
        CMPFN='rcmp3d09.asc'
      END IF
      IF(NRCAL3D.EQ.10) THEN
        SALFN='rsal3d10.asc'
        TEMFN='rtem3d10.asc'
        DYEFN='rdye3d10.asc'
        SEDFN='rsed3d10.asc'
        SNDFN='rsnd3d10.asc'
        TOXFN='rtox3d10.asc'
        UUUFN='ruuu3d10.asc'
        VVVFN='rvvv3d10.asc'
        WWWFN='rwww3d10.asc'
        CMPFN='rcmp3d10.asc'
      END IF
      IF(NRCAL3D.EQ.11) THEN
        SALFN='rsal3d11.asc'
        TEMFN='rtem3d11.asc'
        DYEFN='rdye3d11.asc'
        SEDFN='rsed3d11.asc'
        SNDFN='rsnd3d11.asc'
        TOXFN='rtox3d11.asc'
        UUUFN='ruuu3d11.asc'
        VVVFN='rvvv3d11.asc'
        WWWFN='rwww3d11.asc'
        CMPFN='rcmp3d11.asc'
      END IF
      IF(NRCAL3D.EQ.12) THEN
        SALFN='rsal3d12.asc'
        TEMFN='rtem3d12.asc'
        DYEFN='rdye3d12.asc'
        SEDFN='rsed3d12.asc'
        SNDFN='rsnd3d12.asc'
        TOXFN='rtox3d12.asc'
        UUUFN='ruuu3d12.asc'
        VVVFN='rvvv3d12.asc'
        WWWFN='rwww3d12.asc'
        CMPFN='rcmp3d12.asc'
      END IF
      IF(NRCAL3D.EQ.13) THEN
        SALFN='rsal3d13.asc'
        TEMFN='rtem3d13.asc'
        DYEFN='rdye3d13.asc'
        SEDFN='rsed3d13.asc'
        SNDFN='rsnd3d13.asc'
        TOXFN='rtox3d13.asc'
        UUUFN='ruuu3d13.asc'
        VVVFN='rvvv3d13.asc'
        WWWFN='rwww3d13.asc'
        CMPFN='rcmp3d13.asc'
      END IF
      IF(NRCAL3D.EQ.14) THEN
        SALFN='rsal3d14.asc'
        TEMFN='rtem3d14.asc'
        DYEFN='rdye3d14.asc'
        SEDFN='rsed3d14.asc'
        SNDFN='rsnd3d14.asc'
        TOXFN='rtox3d14.asc'
        UUUFN='ruuu3d14.asc'
        VVVFN='rvvv3d14.asc'
        WWWFN='rwww3d14.asc'
        CMPFN='rcmp3d14.asc'
      END IF
      IF(NRCAL3D.EQ.15) THEN
        SALFN='rsal3d15.asc'
        TEMFN='rtem3d15.asc'
        DYEFN='rdye3d15.asc'
        SEDFN='rsed3d15.asc'
        SNDFN='rsnd3d15.asc'
        TOXFN='rtox3d15.asc'
        UUUFN='ruuu3d15.asc'
        VVVFN='rvvv3d15.asc'
        WWWFN='rwww3d15.asc'
        CMPFN='rcmp3d15.asc'
      END IF
      IF(NRCAL3D.EQ.16) THEN
        SALFN='rsal3d16.asc'
        TEMFN='rtem3d16.asc'
        DYEFN='rdye3d16.asc'
        SEDFN='rsed3d16.asc'
        SNDFN='rsnd3d16.asc'
        TOXFN='rtox3d16.asc'
        UUUFN='ruuu3d16.asc'
        VVVFN='rvvv3d16.asc'
        WWWFN='rwww3d16.asc'
        CMPFN='rcmp3d16.asc'
      END IF
      IF(NRCAL3D.EQ.17) THEN
        SALFN='rsal3d17.asc'
        TEMFN='rtem3d17.asc'
        DYEFN='rdye3d17.asc'
        SEDFN='rsed3d17.asc'
        SNDFN='rsnd3d17.asc'
        TOXFN='rtox3d17.asc'
        UUUFN='ruuu3d17.asc'
        VVVFN='rvvv3d17.asc'
        WWWFN='rwww3d17.asc'
        CMPFN='rcmp3d17.asc'
      END IF
      IF(NRCAL3D.EQ.18) THEN
        SALFN='rsal3d18.asc'
        TEMFN='rtem3d18.asc'
        DYEFN='rdye3d18.asc'
        SEDFN='rsed3d18.asc'
        SNDFN='rsnd3d18.asc'
        TOXFN='rtox3d18.asc'
        UUUFN='ruuu3d18.asc'
        VVVFN='rvvv3d18.asc'
        WWWFN='rwww3d18.asc'
        CMPFN='rcmp3d18.asc'
      END IF
      IF(NRCAL3D.EQ.19) THEN
        SALFN='rsal3d19.asc'
        TEMFN='rtem3d19.asc'
        DYEFN='rdye3d19.asc'
        SEDFN='rsed3d19.asc'
        SNDFN='rsnd3d19.asc'
        TOXFN='rtox3d19.asc'
        UUUFN='ruuu3d19.asc'
        VVVFN='rvvv3d19.asc'
        WWWFN='rwww3d19.asc'
        CMPFN='rcmp3d19.asc'
      END IF
      IF(NRCAL3D.EQ.20) THEN
        SALFN='rsal3d20.asc'
        TEMFN='rtem3d20.asc'
        DYEFN='rdye3d20.asc'
        SEDFN='rsed3d20.asc'
        SNDFN='rsnd3d10.asc'
        TOXFN='rtox3d10.asc'
        UUUFN='ruuu3d20.asc'
        VVVFN='rvvv3d20.asc'
        WWWFN='rwww3d20.asc'
        CMPFN='rcmp3d20.asc'
      END IF
      IF(NRCAL3D.EQ.21) THEN
        SALFN='rsal3d21.asc'
        TEMFN='rtem3d21.asc'
        DYEFN='rdye3d21.asc'
        SEDFN='rsed3d21.asc'
        SNDFN='rsnd3d21.asc'
        TOXFN='rtox3d21.asc'
        UUUFN='ruuu3d21.asc'
        VVVFN='rvvv3d21.asc'
        WWWFN='rwww3d21.asc'
        CMPFN='rcmp3d21.asc'
      END IF
      IF(NRCAL3D.EQ.22) THEN
        SALFN='rsal3d22.asc'
        TEMFN='rtem3d22.asc'
        DYEFN='rdye3d22.asc'
        SEDFN='rsed3d22.asc'
        SNDFN='rsnd3d22.asc'
        TOXFN='rtox3d22.asc'
        UUUFN='ruuu3d22.asc'
        VVVFN='rvvv3d22.asc'
        WWWFN='rwww3d22.asc'
        CMPFN='rcmp3d22.asc'
      END IF
      IF(NRCAL3D.EQ.23) THEN
        SALFN='rsal3d23.asc'
        TEMFN='rtem3d23.asc'
        DYEFN='rdye3d23.asc'
        SEDFN='rsed3d23.asc'
        SNDFN='rsnd3d23.asc'
        TOXFN='rtox3d23.asc'
        UUUFN='ruuu3d23.asc'
        VVVFN='rvvv3d23.asc'
        WWWFN='rwww3d23.asc'
        CMPFN='rcmp3d23.asc'
      END IF
      IF(NRCAL3D.EQ.24) THEN
        SALFN='rsal3d24.asc'
        TEMFN='rtem3d24.asc'
        DYEFN='rdye3d24.asc'
        SEDFN='rsed3d24.asc'
        SNDFN='rsnd3d24.asc'
        TOXFN='rtox3d24.asc'
        UUUFN='ruuu3d24.asc'
        VVVFN='rvvv3d24.asc'
        WWWFN='rwww3d24.asc'
        CMPFN='rcmp3d24.asc'
      END IF
C
      IF (NRCAL3D.EQ.1) THEN
        OPEN(50,FILE='rout3d.dia',STATUS='UNKNOWN')
        CLOSE(50,STATUS='DELETE')
        OPEN(50,FILE='rout3d.dia',STATUS='UNKNOWN')
        WRITE(50,520)IAD,JAD
        WRITE(50,530)NRCAL3D
        DO KP=1,KPC
        WRITE(50,502)KP,ZZP(KP)
        END DO
       ELSE
        OPEN(50,FILE='rout3d.dia',ACCESS='APPEND',STATUS='UNKNOWN')
        WRITE(50,530)NRCAL3D
       END IF             
C
      IF(IS3DSAL.GE.1) THEN
         ASALMAX=-99999999.
         ASALMIN= 99999999.
         DO K=1,KC
         DO L=2,LA
         TMPVAL=SALLPF(L,K)
         ASALMAX=MAX(ASALMAX,TMPVAL)
         ASALMIN=MIN(ASALMIN,TMPVAL)
         END DO
         END DO
         WRITE(50,521)ASALMAX,ASALMIN
         IF (JS3DSAL.EQ.0) THEN
           SAL3DMA=255.
           SAL3DMI=0.
         END IF
         IF (JS3DSAL.EQ.1) THEN
           SAL3DMA=ASALMAX
           SAL3DMI=ASALMIN
         END IF
         OPEN(51,FILE=SALFN,STATUS='UNKNOWN')
         CLOSE(51,STATUS='DELETE')
         OPEN(51,FILE=SALFN,STATUS='UNKNOWN')
      END IF
      IF(IS3DTEM.GE.1) THEN
         ATEMMAX=-99999999.
         ATEMMIN= 99999999.
         DO K=1,KC
         DO L=2,LA
         TMPVAL=TEM(L,K)
         ATEMMAX=MAX(ATEMMAX,TMPVAL)
         ATEMMIN=MIN(ATEMMIN,TMPVAL)
         END DO
         END DO
         WRITE(50,522)ATEMMAX,ATEMMIN
         IF (JS3DTEM.EQ.0) THEN
           TEM3DMA=255.
           TEM3DMI=0.
         END IF
         IF (JS3DTEM.EQ.1) THEN
           TEM3DMA=ATEMMAX
           TEM3DMI=ATEMMIN
         END IF
         OPEN(52,FILE=TEMFN,STATUS='UNKNOWN')
         CLOSE(52,STATUS='DELETE')
         OPEN(52,FILE=TEMFN,STATUS='UNKNOWN')
      END IF
      IF(IS3DDYE.GE.1) THEN
         ADYEMAX=-99999999.
         ADYEMIN= 99999999.
         DO K=1,KC
         DO L=2,LA
         TMPVAL=DYELPF(L,K)
         ADYEMAX=MAX(ADYEMAX,TMPVAL)
         ADYEMIN=MIN(ADYEMIN,TMPVAL)
         END DO
         END DO
         WRITE(50,523)ADYEMAX,ADYEMIN
         IF (JS3DDYE.EQ.0) THEN
           DYE3DMA=255.
           DYE3DMI=0.
         END IF
         IF (JS3DDYE.EQ.1) THEN
           DYE3DMA=ADYEMAX
           DYE3DMI=ADYEMIN
         END IF
         OPEN(53,FILE=DYEFN,STATUS='UNKNOWN')
         CLOSE(53,STATUS='DELETE')
         OPEN(53,FILE=DYEFN,STATUS='UNKNOWN')
      END IF
      IF(IS3DSED.GE.1) THEN
         ASEDMAX=-99999999.
         ASEDMIN= 99999999.
         DO K=1,KC
         DO L=2,LA
         TMPVAL=SEDLPF(L,K,1)
         ASEDMAX=MAX(ASEDMAX,TMPVAL)
         ASEDMIN=MIN(ASEDMIN,TMPVAL)
         END DO
         END DO
         WRITE(50,524)ASEDMAX,ASEDMIN
         IF (JS3DSED.EQ.0) THEN
           SED3DMA=255.
           SED3DMI=0.
         END IF
         IF (JS3DSED.EQ.1) THEN
           SED3DMA=ASEDMAX
           SED3DMI=ASEDMIN
         END IF
         OPEN(54,FILE=SEDFN,STATUS='UNKNOWN')
         CLOSE(54,STATUS='DELETE')
         OPEN(54,FILE=SEDFN,STATUS='UNKNOWN')
      END IF
      IF(IS3DUUU.GE.1) THEN
         AUUUMAX=-99999999.
         AUUUMIN= 99999999.
         DO K=1,KC
         DO L=2,LA
         TMPVAL=0.5*(ULPF(L,K)+ULPF(L+1,K))
         AUUUMAX=MAX(AUUUMAX,TMPVAL)
         AUUUMIN=MIN(AUUUMIN,TMPVAL)
         END DO
         END DO
         WRITE(50,525)AUUUMAX,AUUUMIN
         IF (JS3DUUU.EQ.0) THEN
           UUU3DMA=255.
           UUU3DMI=0.
         END IF
         IF (JS3DUUU.EQ.1) THEN
           UUU3DMA=AUUUMAX
           UUU3DMI=AUUUMIN
         END IF
         OPEN(55,FILE=UUUFN,STATUS='UNKNOWN')
         CLOSE(55,STATUS='DELETE')
         OPEN(55,FILE=UUUFN,STATUS='UNKNOWN')
      END IF
      IF(IS3DVVV.GE.1) THEN
         AVVVMAX=-99999999.
         AVVVMIN= 99999999.
         DO K=1,KC
         DO L=2,LA
         TMPVAL=0.5*(VLPF(L,K)+VLPF(LNC(L),K))
         AVVVMAX=MAX(AVVVMAX,TMPVAL)
         AVVVMIN=MIN(AVVVMIN,TMPVAL)
         END DO
         END DO
         WRITE(50,526)AVVVMAX,AVVVMIN
         IF (JS3DVVV.EQ.0) THEN
           VVV3DMA=255.
           VVV3DMI=0.
         END IF
         IF (JS3DVVV.EQ.1) THEN
           VVV3DMA=AVVVMAX
           VVV3DMI=AVVVMIN
         END IF
         OPEN(56,FILE=VVVFN,STATUS='UNKNOWN')
         CLOSE(56,STATUS='DELETE')
         OPEN(56,FILE=VVVFN,STATUS='UNKNOWN')
      END IF
      IF(IS3DWWW.GE.1) THEN
         AWWWMAX=-99999999.
         AWWWMIN= 99999999.
         DO K=1,KC
         DO L=2,LA
         LN=LNC(L)
         LS=LSC(L)
         TMPVAL=0.5*(WLPF(L,K)+WLPF(L,K-1))
     $       +ZZ(K)*( GI*0.0*DTI*(P(L)-P1(L))
     $       +0.5*(ULPF(L+1,K)*(HLPF(L+1)+BELV(L+1)
     $                      -HLPF(L)-BELV(L))*DXIU(L+1)
     $            +ULPF(L,K)*(HLPF(L)+BELV(L)
     $                    -HLPF(L-1)-BELV(L-1))*DXIU(L)
     $            +VLPF(LN,K)*(HLPF(LN)+BELV(LN)
     $                     -HLPF(L)-BELV(L))*DYIV(LNC(L))
     $            +VLPF(L,K)*(HLPF(L)+BELV(L)
     $                    -HLPF(LS)-BELV(LS))*DYIV(L)) )
     $    +0.5*(1.-ZZ(K))*( ULPF(L+1,K)*(BELV(L+1)-BELV(L))*DXIU(L+1)
     $                     +ULPF(L,K)*(BELV(L)-BELV(L-1))*DXIU(L)
     $                     +VLPF(LN,K)*(BELV(LN)-BELV(L))*DYIV(LN)
     $                     +VLPF(L,K)*(BELV(L)-BELV(LS))*DYIV(L) )
         AWWWMAX=MAX(AWWWMAX,TMPVAL)
         AWWWMIN=MIN(AWWWMIN,TMPVAL)
         END DO
         END DO
         WRITE(50,527)AWWWMAX,AWWWMIN
         IF (JS3DWWW.EQ.0) THEN
           WWW3DMA=255.
           WWW3DMI=0.
         END IF
         IF (JS3DWWW.EQ.1) THEN
           WWW3DMA=AWWWMAX
           WWW3DMI=AWWWMIN
         END IF
         OPEN(57,FILE=WWWFN,STATUS='UNKNOWN')
         CLOSE(57,STATUS='DELETE')
         OPEN(57,FILE=WWWFN,STATUS='UNKNOWN')
      END IF
      IF(IS3DSND.GE.1) THEN
         ASNDMAX=-99999999.
         ASNDMIN= 99999999.
         DO K=1,KC
         DO L=2,LA
         TMPVAL=SNDLPF(L,K,1)
         ASNDMAX=MAX(ASNDMAX,TMPVAL)
         ASNDMIN=MIN(ASNDMIN,TMPVAL)
         END DO
         END DO
         WRITE(50,528)ASNDMAX,ASNDMIN
         IF (JS3DSND.EQ.0) THEN
           SND3DMA=255.
           SED3DMI=0.
         END IF
         IF (JS3DSED.EQ.1) THEN
           SND3DMA=ASNDMAX
           SND3DMI=ASNDMIN
         END IF
         OPEN(58,FILE=SNDFN,STATUS='UNKNOWN')
         CLOSE(58,STATUS='DELETE')
         OPEN(58,FILE=SNDFN,STATUS='UNKNOWN')
      END IF
      IF(IS3DTOX.GE.1) THEN
         ATOXMAX=-99999999.
         ATOXMIN= 99999999.
         DO K=1,KC
         DO L=2,LA
         TMPVAL=TOXLPF(L,K,1)
         ASEDMAX=MAX(ATOXMAX,TMPVAL)
         ASEDMIN=MIN(ATOXMIN,TMPVAL)
         END DO
         END DO
         WRITE(50,529)ATOXMAX,ATOXMIN
         IF (JS3DTOX.EQ.0) THEN
           TOX3DMA=255.
           TOX3DMI=0.
         END IF
         IF (JS3DTOX.EQ.1) THEN
           TOX3DMA=ATOXMAX
           TOX3DMI=ATOXMIN
         END IF
         OPEN(59,FILE=TOXFN,STATUS='UNKNOWN')
         CLOSE(59,STATUS='DELETE')
         OPEN(59,FILE=TOXFN,STATUS='UNKNOWN')
      END IF
C
      OPEN(99,FILE=CMPFN,STATUS='UNKNOWN')
      CLOSE(99,STATUS='DELETE')
      OPEN(99,FILE=CMPFN,STATUS='UNKNOWN')
C
C**********************************************************************C
C
C **  BEGIN LOOP TO LOAD OUTPUT FILES 
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
C
      DO KP=1,KPC
      IAP(KP)=0
      AP(KP)=0.
      END DO
C
      DO K=1,KC
      DO KP=1,KPC
      APT(KP,K)=0.
      END DO
      END DO
C
      DO KP=1,KPC
      ZZPS=(ZZP(KP)-BELV(L))/HLPF(L)
      IF (ZZPS.GE.0.) THEN
        KPB(L)=KP
        GO TO 190
      END IF
      END DO
  190 CONTINUE
      DO KP=KPC,1,-1
      ZZPS=(ZZP(KP)-BELV(L))*HPI(L)
      IF (ZZPS.LE.1.) THEN
        KPS(L)=KP
        GO TO 195
      END IF
      END DO
  195 CONTINUE
C
C*DIAGNOSTIC
C     WRITE(50,500)NRCAL3D
C     WRITE(50,500)L,IL(L),JL(L),KPB(L),KPS(L)
C*DIAGNOSTIC
C
      DO KP=KPB(L),KPS(L)
      ZZPS=(ZZP(KP)-BELV(L))*HPI(L)
      IF (ZZPS.GE.0.0.AND.ZZPS.LE.1.0) THEN
        IF (ZZPS.GE.ZZ(KC)) THEN
          APT(KP,KC)= (ZZPS-ZZ(KS))/(ZZ(KC)-ZZ(KS))
          APT(KP,KS)=-(ZZPS-ZZ(KC))/(ZZ(KC)-ZZ(KS))
         ELSE
          IF (ZZPS.LE.ZZ(1)) THEN
            APT(KP,2)= (ZZPS-ZZ(1))/(ZZ(2)-ZZ(1))
            APT(KP,1)=-(ZZPS-ZZ(2))/(ZZ(2)-ZZ(1))
           ELSE
            K=1
  200       K=K+1
            IF (ZZPS.GT.ZZ(K-1).AND.ZZPS.LE.ZZ(K)) THEN
              APT(KP,K)  = (ZZPS-ZZ(K-1))/(ZZ(K)-ZZ(K-1))
              APT(KP,K-1)=-(ZZPS-ZZ(K))/(ZZ(K)-ZZ(K-1))
             ELSE
              GO TO 200
            END IF
          END IF
        END IF
      END IF
      END DO
C
C*DIAGNOSTIC
C     IF(NRCAL3D.EQ.1) THEN
C       IF(L.EQ.2) THEN
C          DO KP=1,KPC
C          WRITE(50,505)(APT(KP,K),K=1,KC)
C          END DO
C       END IF
C     END IF
C*DIAGNOSTIC
C
      DO K=1,KC
       TMP3D(K)=1.0
      END DO
      DO KP=KPB(L),KPS(L)
       AP(KP)=0.
       DO K=1,KC
        AP(KP)=AP(KP)+APT(KP,K)*TMP3D(K)
       END DO
      END DO
      DO KP=KPB(L),KPS(L)
       IAP(KP)=NINT(AP(KP))
      END DO
      WRITE(99,559)IL(L),JL(L),(IAP(K),K=1,KPC)
C
      IF(IS3DSAL.GE.1) THEN
        DO K=1,KC
        TMP3D(K)=SALLPF(L,K)
        END DO
        SCALE3D=254./(SAL3DMA-SAL3DMI)
        DO KP=KPB(L),KPS(L)
        AP(KP)=0.
        DO K=1,KC
        AP(KP)=AP(KP)+APT(KP,K)*TMP3D(K)
C*DIAGNOSTIC
C       IF(NRCAL3D.EQ.1) THEN
C         IF(L.EQ.2) THEN
C         WRITE(50,510)KP,K,AP(KP),APT(KP,K),TMP3D(K),SAL(L,K)
C         END IF
C       END IF
C*DIAGNOSTIC
        END DO
        END DO
        DO KP=KPB(L),KPS(L)
        IAP(KP)=NINT((AP(KP)-SAL3DMI)*SCALE3D)+1
        IF(IAP(KP).GT.255) IAP(KP)=255
        AP(KP)=SAL3DMA*AP(KP)
        END DO
        IF(JS3DSAL.LE.2) WRITE(51,501)(IAP(K),K=1,KPC)
        IF(JS3DSAL.EQ.3) WRITE(51,551)(AP(K),K=1,KPC)
C*DIAGNOSTIC
C       IF(NCALL3D.EQ.1) THEN
C         IF(L.EQ.2) THEN
C           DO KP=1,KPC
C           WRITE(50,506)KP,AP(KP),IAP(KP)
C           END DO
C         END IF
C       END IF
C*DIAGNOSTIC
      END IF
      IF(IS3DTEM.GE.1) THEN
        DO K=1,KC
        TMP3D(K)=TEM(L,K)
        END DO
        SCALE3D=254./(TEM3DMA-TEM3DMI)
        DO KP=KPB(L),KPS(L)
        AP(KP)=0.
        DO K=1,KC
        AP(KP)=AP(KP)+APT(KP,K)*TMP3D(K)
        END DO
        END DO
        DO KP=KPB(L),KPS(L)
        IAP(KP)=NINT((AP(KP)-TEM3DMI)*SCALE3D)+1
        IF(IAP(KP).GT.255) IAP(KP)=255
        AP(KP)=TEM3DMA*AP(KP)
        END DO
        IF(JS3DTEM.LE.2) WRITE(52,501)(IAP(K),K=1,KPC)
        IF(JS3DTEM.EQ.3) WRITE(52,551)(AP(K),K=1,KPC)
      END IF
      IF(IS3DDYE.GE.1) THEN
        DO K=1,KC
        TMP3D(K)=DYELPF(L,K)
        END DO
        SCALE3D=254./(DYE3DMA-DYE3DMI)
        DO KP=KPB(L),KPS(L)
        AP(KP)=0.
        DO K=1,KC
        AP(KP)=AP(KP)+APT(KP,K)*TMP3D(K)
        END DO
        END DO
        DO KP=KPB(L),KPS(L)
        IAP(KP)=NINT((AP(KP)-DYE3DMI)*SCALE3D)+1
        IF(IAP(KP).GT.255) IAP(KP)=255
        AP(KP)=DYE3DMA*AP(KP)
        END DO
        IF(JS3DDYE.LE.2) WRITE(53,501)(IAP(K),K=1,KPC)
        IF(JS3DDYE.EQ.3) WRITE(53,551)(AP(K),K=1,KPC)
      END IF
      IF(IS3DSED.GE.1) THEN
        DO K=1,KC
        TMP3D(K)=SEDLPF(L,K,1)
        END DO
        SCALE3D=254./(SED3DMA-SED3DMI)
        DO KP=KPB(L),KPS(L)
        AP(KP)=0.
        DO K=1,KC
        AP(KP)=AP(KP)+APT(KP,K)*TMP3D(K)
        END DO
        END DO
        DO KP=KPB(L),KPS(L)
        IAP(KP)=NINT((AP(KP)-SED3DMI)*SCALE3D)+1
        IF(IAP(KP).GT.255) IAP(KP)=255
        AP(KP)=SED3DMA*AP(KP)
        END DO
        IF(JS3DSED.LE.2) WRITE(54,501)(IAP(K),K=1,KPC)
        IF(JS3DSED.EQ.3) WRITE(54,551)(AP(K),K=1,KPC)
      END IF
      IF(IS3DUUU.GE.1) THEN
        DO K=1,KC
        TMP3D(K)=0.5*(ULPF(L,K)+ULPF(L+1,K))
        END DO
        SCALE3D=254./(UUU3DMA-UUU3DMI)
        DO KP=KPB(L),KPS(L)
        AP(KP)=0.
        DO K=1,KC
        AP(KP)=AP(KP)+APT(KP,K)*TMP3D(K)
        END DO
        END DO
        DO KP=KPB(L),KPS(L)
        IAP(KP)=NINT((AP(KP)-UUU3DMI)*SCALE3D)+1
        IF(IAP(KP).GT.255) IAP(KP)=255
        AP(KP)=UUU3DMA*AP(KP)
        END DO
        IF(JS3DUUU.LE.2) WRITE(55,501)(IAP(K),K=1,KPC)
        IF(JS3DUUU.EQ.3) WRITE(55,551)(AP(K),K=1,KPC)
      END IF
      IF(IS3DVVV.GE.1) THEN
        DO K=1,KC
        TMP3D(K)=0.5*(VLPF(L,K)+VLPF(LN,K))
        END DO
        SCALE3D=254./(VVV3DMA-VVV3DMI)
        DO KP=KPB(L),KPS(L)
        AP(KP)=0.
        DO K=1,KC
        AP(KP)=AP(KP)+APT(KP,K)*TMP3D(K)
        END DO
        END DO
        DO KP=KPB(L),KPS(L)
        IAP(KP)=NINT((AP(KP)-VVV3DMI)*SCALE3D)+1
        IF(IAP(KP).GT.255) IAP(KP)=255
        AP(KP)=VVV3DMA*AP(KP)
        END DO
        IF(JS3DVVV.LE.2) WRITE(56,501)(IAP(K),K=1,KPC)
        IF(JS3DVVV.EQ.3) WRITE(56,551)(AP(K),K=1,KPC)
      END IF
      IF(IS3DWWW.GE.1) THEN
        DO K=1,KC
        TMP3D(K)=0.5*(WLPF(L,K)+WLPF(L,K-1))
     $         +ZZ(K)*( GI*0.0*DTI*(P(L)-P1(L))
     $         +0.5*(ULPF(L+1,K)*(HLPF(L+1)+BELV(L+1)
     $                        -HLPF(L)-BELV(L))*DXIU(L+1)
     $              +ULPF(L,K)*(HLPF(L)+BELV(L)
     $                      -HLPF(L-1)-BELV(L-1))*DXIU(L)
     $              +VLPF(LN,K)*(HLPF(LN)+BELV(LN)
     $                       -HLPF(L)-BELV(L))*DYIV(LNC(L))
     $              +VLPF(L,K)*(HLPF(L)+BELV(L)
     $                      -HLPF(LS)-BELV(LS))*DYIV(L)) )
     $      +0.5*(1.-ZZ(K))*( ULPF(L+1,K)*(BELV(L+1)-BELV(L))*DXIU(L+1)
     $                       +ULPF(L,K)*(BELV(L)-BELV(L-1))*DXIU(L)
     $                       +VLPF(LN,K)*(BELV(LN)-BELV(L))*DYIV(LN)
     $                       +VLPF(L,K)*(BELV(L)-BELV(LS))*DYIV(L) )
        END DO
        SCALE3D=254./(WWW3DMA-WWW3DMI)
        DO KP=KPB(L),KPS(L)
        AP(KP)=0.
        DO K=1,KC
        AP(KP)=AP(KP)+APT(KP,K)*TMP3D(K)
        END DO
        END DO
        DO KP=KPB(L),KPS(L)
        IAP(KP)=NINT((AP(KP)-WWW3DMI)*SCALE3D)+1
        IF(IAP(KP).GT.255) IAP(KP)=255
        AP(KP)=WWW3DMA*AP(KP)
        END DO
        IF(JS3DWWW.LE.2) WRITE(57,501)(IAP(K),K=1,KPC)
        IF(JS3DWWW.EQ.3) WRITE(57,551)(AP(K),K=1,KPC)
      END IF
      IF(IS3DSND.GE.1) THEN
        DO K=1,KC
        TMP3D(K)=SNDLPF(L,K,1)
        END DO
        SCALE3D=254./(SND3DMA-SND3DMI)
        DO KP=KPB(L),KPS(L)
        AP(KP)=0.
        DO K=1,KC
        AP(KP)=AP(KP)+APT(KP,K)*TMP3D(K)
        END DO
        END DO
        DO KP=KPB(L),KPS(L)
        IAP(KP)=NINT((AP(KP)-SND3DMI)*SCALE3D)+1
        IF(IAP(KP).GT.255) IAP(KP)=255
        AP(KP)=SND3DMA*AP(KP)
        END DO
        IF(JS3DSND.LE.2) WRITE(58,501)(IAP(K),K=1,KPC)
        IF(JS3DSND.EQ.3) WRITE(58,551)(AP(K),K=1,KPC)
      END IF
      IF(IS3DTOX.GE.1) THEN
        DO K=1,KC
        TMP3D(K)=TOXLPF(L,K,1)
        END DO
        SCALE3D=254./(TOX3DMA-TOX3DMI)
        DO KP=KPB(L),KPS(L)
        AP(KP)=0.
        DO K=1,KC
        AP(KP)=AP(KP)+APT(KP,K)*TMP3D(K)
        END DO
        END DO
        DO KP=KPB(L),KPS(L)
        IAP(KP)=NINT((AP(KP)-TOX3DMI)*SCALE3D)+1
        IF(IAP(KP).GT.255) IAP(KP)=255
        AP(KP)=TOX3DMA*AP(KP)
        END DO
        IF(JS3DSND.LE.2) WRITE(59,501)(IAP(K),K=1,KPC)
        IF(JS3DSND.EQ.3) WRITE(59,551)(AP(K),K=1,KPC)
      END IF
C
      END DO
C
      IF(IS3DSAL.GE.1) CLOSE(51)
      IF(IS3DTEM.GE.1) CLOSE(52)
      IF(IS3DDYE.GE.1) CLOSE(53)
      IF(IS3DSED.GE.1) CLOSE(54)
      IF(IS3DUUU.GE.1) CLOSE(55)
      IF(IS3DVVV.GE.1) CLOSE(56)
      IF(IS3DWWW.GE.1) CLOSE(57)
      IF(IS3DSND.GE.1) CLOSE(58)
      IF(IS3DTOX.GE.1) CLOSE(59)
      CLOSE(99)
C
C**********************************************************************C
C
C **  REWRITE OUTPUT ARRAYS INTO CORRECT ORDER IF I3DRW.EQ.1
C
C----------------------------------------------------------------------C
C
      IF(I3DRW.EQ.1) THEN
C
      DO J=1,JG
      DO I=1,IG
      IAIJ(I,J)=0
      AIJ(I,J)=0
      END DO
      END DO
C
      IF(ISCLO.EQ.0.OR.NWGG.EQ.0) THEN
C  
      IF(IS3DSAL.GE.1) THEN
         OPEN(51,FILE=SALFN,STATUS='UNKNOWN')
         DO L=2,LA
          IF(JS3DSAL.LE.2) READ(51,*)(IAKL(K,L),K=1,KPC)
          IF(JS3DSAL.EQ.3) READ(51,*)(AKL(K,L),K=1,KPC)
         END DO
         CLOSE(51,STATUS='DELETE')
         OPEN(51,FILE=SALFN,STATUS='UNKNOWN')
         IF(JS3DSAL.LE.2) THEN
         DO K=1,KPC
          DO L=2,LA
          IAIJ(IL(L),JL(L))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DSAL.GT.0) THEN
            IF(IJCT(I3DMIN,J3DMIN).EQ.0) THEN
              IAIJ(I3DMIN,J3DMIN)=255
            END IF
          END IF
          DO J=J3DMAX,J3DMIN,-1         
           WRITE(51,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         END IF
         IF(JS3DSAL.GE.3.) THEN
          DO K=1,KPC
           DO L=2,LA
            AIJ(IL(L),JL(L))=AKL(K,L)
           END DO
           DO J=J3DMAX,J3DMIN,-1         
            WRITE(51,551)(AIJ(I,J),I=I3DMAX,I3DMIN,-1)
           END DO
          END DO
         END IF
         CLOSE(51)
      END IF
      IF(IS3DTEM.GE.1) THEN
         OPEN(52,FILE=TEMFN,STATUS='UNKNOWN')
         DO L=2,LA
          IF(JS3DTEM.LE.2) READ(52,*)(IAKL(K,L),K=1,KPC)
          IF(JS3DTEM.EQ.3) READ(52,*)(AKL(K,L),K=1,KPC)
         END DO
         CLOSE(52,STATUS='DELETE')
         OPEN(52,FILE=TEMFN,STATUS='UNKNOWN')
         IF(JS3DTEM.LE.2) THEN
         DO K=1,KPC
          DO L=2,LA
          IAIJ(IL(L),JL(L))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DTEM.GT.0) THEN
            IF(IJCT(I3DMIN,J3DMIN).EQ.0) THEN
              IAIJ(I3DMIN,J3DMIN)=255
            END IF
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(52,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         END IF
         IF(JS3DTEM.GE.3.) THEN
          DO K=1,KPC
           DO L=2,LA
            AIJ(IL(L),JL(L))=AKL(K,L)
           END DO
           DO J=J3DMAX,J3DMIN,-1         
            WRITE(52,551)(AIJ(I,J),I=I3DMAX,I3DMIN,-1)
           END DO
          END DO
         END IF
         CLOSE(52)
      END IF
      IF(IS3DDYE.GE.1) THEN
         OPEN(53,FILE=DYEFN,STATUS='UNKNOWN')
         DO L=2,LA
          IF(JS3DDYE.LE.2) READ(53,*)(IAKL(K,L),K=1,KPC)
          IF(JS3DDYE.EQ.3) READ(53,*)(AKL(K,L),K=1,KPC)
         END DO
         CLOSE(53,STATUS='DELETE')
         OPEN(53,FILE=DYEFN,STATUS='UNKNOWN')
         IF(JS3DDYE.LE.2) THEN
         DO K=1,KPC
          DO L=2,LA
          IAIJ(IL(L),JL(L))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DDYE.GT.0) THEN
            IF(IJCT(I3DMIN,J3DMIN).EQ.0) THEN
              IAIJ(I3DMIN,J3DMIN)=255
            END IF
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(53,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         END IF
         IF(JS3DDYE.GE.3.) THEN
          DO K=1,KPC
           DO L=2,LA
            AIJ(IL(L),JL(L))=AKL(K,L)
           END DO
           DO J=J3DMAX,J3DMIN,-1         
            WRITE(53,551)(AIJ(I,J),I=I3DMAX,I3DMIN,-1)
           END DO
          END DO
         END IF
         CLOSE(53)
      END IF
      IF(IS3DSED.GE.1) THEN
         OPEN(54,FILE=SEDFN,STATUS='UNKNOWN')
         DO L=2,LA
          IF(JS3DSED.LE.2) READ(54,*)(IAKL(K,L),K=1,KPC)
          IF(JS3DSED.EQ.3) READ(54,*)(AKL(K,L),K=1,KPC)
         END DO
         CLOSE(54,STATUS='DELETE')
         OPEN(54,FILE=SEDFN,STATUS='UNKNOWN')
         IF(JS3DSED.LE.2) THEN
         DO K=1,KPC
          DO L=2,LA
          IAIJ(IL(L),JL(L))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DSED.GT.0) THEN
            IF(IJCT(I3DMIN,J3DMIN).EQ.0) THEN
              IAIJ(I3DMIN,J3DMIN)=255
            END IF
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(54,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         END IF
         IF(JS3DSED.GE.3.) THEN
          DO K=1,KPC
           DO L=2,LA
            AIJ(IL(L),JL(L))=AKL(K,L)
           END DO
           DO J=J3DMAX,J3DMIN,-1         
            WRITE(54,551)(AIJ(I,J),I=I3DMAX,I3DMIN,-1)
           END DO
          END DO
         END IF
         CLOSE(54)
      END IF
      IF(IS3DUUU.GE.1) THEN
         OPEN(55,FILE=UUUFN,STATUS='UNKNOWN')
         DO L=2,LA
          IF(JS3DUUU.LE.2) READ(55,*)(IAKL(K,L),K=1,KPC)
          IF(JS3DUUU.EQ.3) READ(55,*)(AKL(K,L),K=1,KPC)
         END DO
         CLOSE(55,STATUS='DELETE')
         OPEN(55,FILE=UUUFN,STATUS='UNKNOWN')
         IF(JS3DUUU.LE.2) THEN
         DO K=1,KPC
          DO L=2,LA
          IAIJ(IL(L),JL(L))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DUUU.GT.0) THEN
            IF(IJCT(I3DMIN,J3DMIN).EQ.0) THEN
              IAIJ(I3DMIN,J3DMIN)=255
            END IF
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(55,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         END IF
         IF(JS3DUUU.GE.3.) THEN
          DO K=1,KPC
           DO L=2,LA
            AIJ(IL(L),JL(L))=AKL(K,L)
           END DO
           DO J=J3DMAX,J3DMIN,-1         
            WRITE(55,551)(AIJ(I,J),I=I3DMAX,I3DMIN,-1)
           END DO
          END DO
         END IF
         CLOSE(55)
      END IF
      IF(IS3DVVV.GE.1) THEN
         OPEN(56,FILE=VVVFN,STATUS='UNKNOWN')
         DO L=2,LA
          IF(JS3DVVV.LE.2) READ(56,*)(IAKL(K,L),K=1,KPC)
          IF(JS3DVVV.EQ.3) READ(56,*)(AKL(K,L),K=1,KPC)
         END DO
         CLOSE(56,STATUS='DELETE')
         OPEN(56,FILE=VVVFN,STATUS='UNKNOWN')
         IF(JS3DVVV.LE.2) THEN
         DO K=1,KPC
          DO L=2,LA
          IAIJ(IL(L),JL(L))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DVVV.GT.0) THEN
            IF(IJCT(I3DMIN,J3DMIN).EQ.0) THEN
              IAIJ(I3DMIN,J3DMIN)=255
            END IF
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(56,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         END IF
         IF(JS3DVVV.GE.3.) THEN
          DO K=1,KPC
           DO L=2,LA
            AIJ(IL(L),JL(L))=AKL(K,L)
           END DO
           DO J=J3DMAX,J3DMIN,-1         
            WRITE(56,551)(AIJ(I,J),I=I3DMAX,I3DMIN,-1)
           END DO
          END DO
         END IF
         CLOSE(56)
      END IF
      IF(IS3DWWW.GE.1) THEN
         OPEN(57,FILE=WWWFN,STATUS='UNKNOWN')
         DO L=2,LA
          IF(JS3DWWW.LE.2) READ(57,*)(IAKL(K,L),K=1,KPC)
          IF(JS3DWWW.EQ.3) READ(57,*)(AKL(K,L),K=1,KPC)
         END DO
         CLOSE(57,STATUS='DELETE')
         OPEN(57,FILE=WWWFN,STATUS='UNKNOWN')
         IF(JS3DWWW.LE.2) THEN
         DO K=1,KPC
          DO L=2,LA
          IAIJ(IL(L),JL(L))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DWWW.GT.0) THEN
            IF(IJCT(I3DMIN,J3DMIN).EQ.0) THEN
              IAIJ(I3DMIN,J3DMIN)=255
            END IF
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(57,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         END IF
         IF(JS3DWWW.GE.3.) THEN
          DO K=1,KPC
           DO L=2,LA
            AIJ(IL(L),JL(L))=AKL(K,L)
           END DO
           DO J=J3DMAX,J3DMIN,-1         
            WRITE(57,551)(AIJ(I,J),I=I3DMAX,I3DMIN,-1)
           END DO
          END DO
         END IF
         CLOSE(57)
      END IF
      IF(IS3DSND.GE.1) THEN
         OPEN(58,FILE=SNDFN,STATUS='UNKNOWN')
         DO L=2,LA
          IF(JS3DSND.LE.2) READ(58,*)(IAKL(K,L),K=1,KPC)
          IF(JS3DSND.EQ.3) READ(58,*)(AKL(K,L),K=1,KPC)
         END DO
         CLOSE(58,STATUS='DELETE')
         OPEN(58,FILE=SNDFN,STATUS='UNKNOWN')
         IF(JS3DSND.LE.2) THEN
         DO K=1,KPC
          DO L=2,LA
           IAIJ(IL(L),JL(L))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DSND.GT.0) THEN
            IF(IJCT(I3DMIN,J3DMIN).EQ.0) THEN
              IAIJ(I3DMIN,J3DMIN)=255
            END IF
          END IF
          DO J=J3DMAX,J3DMIN,-1         
           WRITE(58,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         END IF
         IF(JS3DSND.GE.3.) THEN
          DO K=1,KPC
           DO L=2,LA
            AIJ(IL(L),JL(L))=AKL(K,L)
           END DO
           DO J=J3DMAX,J3DMIN,-1         
            WRITE(58,551)(AIJ(I,J),I=I3DMAX,I3DMIN,-1)
           END DO
          END DO
         END IF
         CLOSE(58)
      END IF
      IF(IS3DTOX.GE.1) THEN
         OPEN(59,FILE=TOXFN,STATUS='UNKNOWN')
         DO L=2,LA
          IF(JS3DTOX.LE.2) READ(59,*)(IAKL(K,L),K=1,KPC)
          IF(JS3DTOX.EQ.3) READ(59,*)(AKL(K,L),K=1,KPC)
         END DO
         CLOSE(59,STATUS='DELETE')
         OPEN(59,FILE=TOXFN,STATUS='UNKNOWN')
         IF(JS3DTOX.LE.2) THEN
         DO K=1,KPC
          DO L=2,LA
           IAIJ(IL(L),JL(L))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DTOX.GT.0) THEN
            IF(IJCT(I3DMIN,J3DMIN).EQ.0) THEN
              IAIJ(I3DMIN,J3DMIN)=255
            END IF
          END IF
          DO J=J3DMAX,J3DMIN,-1         
           WRITE(54,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         END IF
         IF(JS3DTOX.GE.3.) THEN
          DO K=1,KPC
           DO L=2,LA
            AIJ(IL(L),JL(L))=AKL(K,L)
           END DO
           DO J=J3DMAX,J3DMIN,-1         
            WRITE(59,551)(AIJ(I,J),I=I3DMAX,I3DMIN,-1)
           END DO
          END DO
         END IF
         CLOSE(59)
      END IF
C
      ELSE
C  
      IF(IS3DSAL.GE.1) THEN
         OPEN(51,FILE=SALFN,STATUS='UNKNOWN')
         DO L=2,LA
         READ(51,*)(IAKL(K,L),K=1,KPC)
         END DO
         CLOSE(51,STATUS='DELETE')
         OPEN(51,FILE=SALFN,STATUS='UNKNOWN')
         DO K=1,KPC
          DO NW=1,NWGG
          L=LWGG(NW)
          IAIJ(IWGG(NW),JWGG(NW))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DSAL.GT.0) THEN
            IAIJ(I3DMIN,J3DMIN)=255
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(51,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         CLOSE(51)
      END IF
      IF(IS3DTEM.GE.1) THEN
         OPEN(52,FILE=TEMFN,STATUS='UNKNOWN')
         DO L=2,LA
         READ(52,*)(IAKL(K,L),K=1,KPC)
         END DO
         CLOSE(52,STATUS='DELETE')
         OPEN(52,FILE=TEMFN,STATUS='UNKNOWN')
         DO K=1,KPC
          DO NW=1,NWGG
          L=LWGG(NW)
          IAIJ(IWGG(NW),JWGG(NW))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DTEM.GT.0) THEN
            IAIJ(I3DMIN,J3DMIN)=255
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(52,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         CLOSE(52)
      END IF
      IF(IS3DDYE.GE.1) THEN
         OPEN(53,FILE=DYEFN,STATUS='UNKNOWN')
         DO L=2,LA
         READ(53,*)(IAKL(K,L),K=1,KPC)
         END DO
         CLOSE(53,STATUS='DELETE')
         OPEN(53,FILE=DYEFN,STATUS='UNKNOWN')
         DO K=1,KPC
          DO NW=1,NWGG
          L=LWGG(NW)
          IAIJ(IWGG(NW),JWGG(NW))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DDYE.GT.0) THEN
            IAIJ(I3DMIN,J3DMIN)=255
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(53,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         CLOSE(53)
      END IF
      IF(IS3DSED.GE.1) THEN
         OPEN(54,FILE=SEDFN,STATUS='UNKNOWN')
         DO L=2,LA
         READ(54,*)(IAKL(K,L),K=1,KPC)
         END DO
         CLOSE(54,STATUS='DELETE')
         OPEN(54,FILE=SEDFN,STATUS='UNKNOWN')
         DO K=1,KPC
          DO NW=1,NWGG
          L=LWGG(NW)
          IAIJ(IWGG(NW),JWGG(NW))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DSED.GT.0) THEN
            IAIJ(I3DMIN,J3DMIN)=255
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(54,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         CLOSE(54)
      END IF
      IF(IS3DUUU.GE.1) THEN
         OPEN(55,FILE=UUUFN,STATUS='UNKNOWN')
         DO L=2,LA
         READ(55,*)(IAKL(K,L),K=1,KPC)
         END DO
         CLOSE(55,STATUS='DELETE')
         OPEN(55,FILE=UUUFN,STATUS='UNKNOWN')
         DO K=1,KPC
          DO NW=1,NWGG
          L=LWGG(NW)
          IAIJ(IWGG(NW),JWGG(NW))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DUUU.GT.0) THEN
            IAIJ(I3DMIN,J3DMIN)=255
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(55,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         CLOSE(55)
      END IF
      IF(IS3DVVV.GE.1) THEN
         OPEN(56,FILE=VVVFN,STATUS='UNKNOWN')
         DO L=2,LA
         READ(56,*)(IAKL(K,L),K=1,KPC)
         END DO
         CLOSE(56,STATUS='DELETE')
         OPEN(56,FILE=VVVFN,STATUS='UNKNOWN')
         DO K=1,KPC
          DO NW=1,NWGG
          L=LWGG(NW)
          IAIJ(IWGG(NW),JWGG(NW))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DVVV.GT.0) THEN
            IAIJ(I3DMIN,J3DMIN)=255
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(56,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         CLOSE(56)
      END IF
      IF(IS3DWWW.GE.1) THEN
         OPEN(57,FILE=WWWFN,STATUS='UNKNOWN')
         DO L=2,LA
         READ(57,*)(IAKL(K,L),K=1,KPC)
         END DO
         CLOSE(57,STATUS='DELETE')
         OPEN(57,FILE=WWWFN,STATUS='UNKNOWN')
         DO K=1,KPC
          DO NW=1,NWGG
          L=LWGG(NW)
          IAIJ(IWGG(NW),JWGG(NW))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DWWW.GT.0) THEN
            IAIJ(I3DMIN,J3DMIN)=255
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(57,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         CLOSE(57)
      END IF
      IF(IS3DSND.GE.1) THEN
         OPEN(58,FILE=SNDFN,STATUS='UNKNOWN')
         DO L=2,LA
         READ(58,*)(IAKL(K,L),K=1,KPC)
         END DO
         CLOSE(58,STATUS='DELETE')
         OPEN(58,FILE=SNDFN,STATUS='UNKNOWN')
         DO K=1,KPC
          DO NW=1,NWGG
          L=LWGG(NW)
          IAIJ(IWGG(NW),JWGG(NW))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DSND.GT.0) THEN
            IAIJ(I3DMIN,J3DMIN)=255
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(58,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         CLOSE(58)
      END IF
      IF(IS3DTOX.GE.1) THEN
         OPEN(59,FILE=TOXFN,STATUS='UNKNOWN')
         DO L=2,LA
         READ(59,*)(IAKL(K,L),K=1,KPC)
         END DO
         CLOSE(59,STATUS='DELETE')
         OPEN(59,FILE=TOXFN,STATUS='UNKNOWN')
         DO K=1,KPC
          DO NW=1,NWGG
          L=LWGG(NW)
          IAIJ(IWGG(NW),JWGG(NW))=IAKL(K,L)
          END DO
          IF(K.EQ.1.AND.JS3DTOX.GT.0) THEN
            IAIJ(I3DMIN,J3DMIN)=255
          END IF
          DO J=J3DMAX,J3DMIN,-1         
          WRITE(59,501)(IAIJ(I,J),I=I3DMAX,I3DMIN,-1)
          END DO
          IAIJ(I3DMIN,J3DMIN)=0
         END DO
         CLOSE(59)
      END IF
C
      END IF
C
      END IF
C
C**********************************************************************C
C
  500 FORMAT(5I5)
  501 FORMAT(72I4)
  502 FORMAT(I5,F10.4)
  505 FORMAT(8F10.5)
  506 FORMAT(I5,2X,F10.5,5X,I5)
  510 FORMAT(2I5,4(2X,F10.5))
  520 FORMAT('IAD = ',I5,'  JAD = ',I5//)
  521 FORMAT('RSALMAX = ',E12.4,'  RSALMIN = ',E12.4/)
  522 FORMAT('RTEMMAX = ',E12.4,'  RTEMMIN = ',E12.4/)
  523 FORMAT('RDYEMAX = ',E12.4,'  RDYEMIN = ',E12.4/)
  524 FORMAT('RSEDMAX = ',E12.4,'  RSEDMIN = ',E12.4/)
  525 FORMAT('RUUUMAX = ',E12.4,'  RUUUMIN = ',E12.4/)
  526 FORMAT('RVVVMAX = ',E12.4,'  RVVVMIN = ',E12.4/)
  527 FORMAT('RWWWMAX = ',E12.4,'  RWWWMIN = ',E12.4/)
  528 FORMAT('RSNDMAX = ',E12.4,'  RSNDMIN = ',E12.4/)
  529 FORMAT('RTOXMAX = ',E12.4,'  RTOXMIN = ',E12.4/)
  530 FORMAT('NRCAL3D = ',I5/)
  551 FORMAT(72F7.1)
  559 FORMAT(2I4,2X,72I2)
C
C----------------------------------------------------------------------C
C
      CLOSE (50)
      RETURN
      END 
