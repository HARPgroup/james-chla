C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE OUT3D
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
      CHARACTER *11 SALFN,TEMFN,DYEFN,SEDFN,UUUFN,VVVFN,WWWFN,
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
      NCALL3D=NCALL3D+1
C
      IF(NCALL3D.EQ.1) THEN
        SALFN='sal3d01.asc'
        TEMFN='tem3d01.asc'
        DYEFN='dye3d01.asc'
        SEDFN='sed3d01.asc'
        SNDFN='snd3d01.asc'
        TOXFN='tox3d01.asc'
        UUUFN='uuu3d01.asc'
        VVVFN='vvv3d01.asc'
        WWWFN='www3d01.asc'
        CMPFN='cmp3d01.asc'
      END IF
      IF(NCALL3D.EQ.2) THEN
        SALFN='sal3d02.asc'
        TEMFN='tem3d02.asc'
        DYEFN='dye3d02.asc'
        SEDFN='sed3d02.asc'
        SNDFN='snd3d02.asc'
        TOXFN='tox3d02.asc'
        UUUFN='uuu3d02.asc'
        VVVFN='vvv3d02.asc'
        WWWFN='www3d02.asc'
        CMPFN='cmp3d02.asc'
      END IF
      IF(NCALL3D.EQ.3) THEN
        SALFN='sal3d03.asc'
        TEMFN='tem3d03.asc'
        DYEFN='dye3d03.asc'
        SEDFN='sed3d03.asc'
        SNDFN='snd3d03.asc'
        TOXFN='tox3d03.asc'
        UUUFN='uuu3d03.asc'
        VVVFN='vvv3d03.asc'
        WWWFN='www3d03.asc'
        CMPFN='cmp3d03.asc'
      END IF
      IF(NCALL3D.EQ.4) THEN
        SALFN='sal3d04.asc'
        TEMFN='tem3d04.asc'
        DYEFN='dye3d04.asc'
        SEDFN='sed3d04.asc'
        SNDFN='snd3d04.asc'
        TOXFN='tox3d04.asc'
        UUUFN='uuu3d04.asc'
        VVVFN='vvv3d04.asc'
        WWWFN='www3d04.asc'
        CMPFN='cmp3d04.asc'
       END IF
       IF(NCALL3D.EQ.5) THEN
        SALFN='sal3d05.asc'
        TEMFN='tem3d05.asc'
        DYEFN='dye3d05.asc'
        SEDFN='sed3d05.asc'
        SNDFN='snd3d05.asc'
        TOXFN='tox3d05.asc'
        UUUFN='uuu3d05.asc'
        VVVFN='vvv3d05.asc'
        WWWFN='www3d05.asc'
        CMPFN='cmp3d05.asc'
      END IF
      IF(NCALL3D.EQ.6) THEN
        SALFN='sal3d06.asc'
        TEMFN='tem3d06.asc'
        DYEFN='dye3d06.asc'
        SEDFN='sed3d06.asc'
        SNDFN='snd3d06.asc'
        TOXFN='tox3d06.asc'
        UUUFN='uuu3d06.asc'
        VVVFN='vvv3d06.asc'
        WWWFN='www3d06.asc'
        CMPFN='cmp3d06.asc'
      END IF
      IF(NCALL3D.EQ.7) THEN
        SALFN='sal3d07.asc'
        TEMFN='tem3d07.asc'
        DYEFN='dye3d07.asc'
        SEDFN='sed3d07.asc'
        SNDFN='snd3d07.asc'
        TOXFN='tox3d07.asc'
        UUUFN='uuu3d07.asc'
        VVVFN='vvv3d07.asc'
        WWWFN='www3d07.asc'
        CMPFN='cmp3d07.asc'
      END IF
      IF(NCALL3D.EQ.8) THEN
        SALFN='sal3d08.asc'
        TEMFN='tem3d08.asc'
        DYEFN='dye3d08.asc'
        SEDFN='sed3d08.asc'
        SNDFN='snd3d08.asc'
        TOXFN='tox3d08.asc'
        UUUFN='uuu3d08.asc'
        VVVFN='vvv3d08.asc'
        WWWFN='www3d08.asc'
        CMPFN='cmp3d08.asc'
      END IF
      IF(NCALL3D.EQ.9) THEN
        SALFN='sal3d09.asc'
        TEMFN='tem3d09.asc'
        DYEFN='dye3d09.asc'
        SEDFN='sed3d09.asc'
        SNDFN='snd3d09.asc'
        TOXFN='tox3d09.asc'
        UUUFN='uuu3d09.asc'
        VVVFN='vvv3d09.asc'
        WWWFN='www3d09.asc'
        CMPFN='cmp3d09.asc'
      END IF
      IF(NCALL3D.EQ.10) THEN
        SALFN='sal3d10.asc'
        TEMFN='tem3d10.asc'
        DYEFN='dye3d10.asc'
        SEDFN='sed3d10.asc'
        SNDFN='snd3d10.asc'
        TOXFN='tox3d10.asc'
        UUUFN='uuu3d10.asc'
        VVVFN='vvv3d10.asc'
        WWWFN='www3d10.asc'
        CMPFN='cmp3d10.asc'
      END IF
      IF(NCALL3D.EQ.11) THEN
        SALFN='sal3d11.asc'
        TEMFN='tem3d11.asc'
        DYEFN='dye3d11.asc'
        SEDFN='sed3d11.asc'
        SNDFN='snd3d11.asc'
        TOXFN='tox3d11.asc'
        UUUFN='uuu3d11.asc'
        VVVFN='vvv3d11.asc'
        WWWFN='www3d11.asc'
        CMPFN='cmp3d11.asc'
      END IF
      IF(NCALL3D.EQ.12) THEN
        SALFN='sal3d12.asc'
        TEMFN='tem3d12.asc'
        DYEFN='dye3d12.asc'
        SEDFN='sed3d12.asc'
        SNDFN='snd3d12.asc'
        TOXFN='tox3d12.asc'
        UUUFN='uuu3d12.asc'
        VVVFN='vvv3d12.asc'
        WWWFN='www3d12.asc'
        CMPFN='cmp3d12.asc'
      END IF
      IF(NCALL3D.EQ.13) THEN
        SALFN='sal3d13.asc'
        TEMFN='tem3d13.asc'
        DYEFN='dye3d13.asc'
        SEDFN='sed3d13.asc'
        SNDFN='snd3d13.asc'
        TOXFN='tox3d13.asc'
        UUUFN='uuu3d13.asc'
        VVVFN='vvv3d13.asc'
        WWWFN='www3d13.asc'
        CMPFN='cmp3d13.asc'
      END IF
      IF(NCALL3D.EQ.14) THEN
        SALFN='sal3d14.asc'
        TEMFN='tem3d14.asc'
        DYEFN='dye3d14.asc'
        SEDFN='sed3d14.asc'
        SNDFN='snd3d14.asc'
        TOXFN='tox3d14.asc'
        UUUFN='uuu3d14.asc'
        VVVFN='vvv3d14.asc'
        WWWFN='www3d14.asc'
        CMPFN='cmp3d14.asc'
      END IF
      IF(NCALL3D.EQ.15) THEN
        SALFN='sal3d15.asc'
        TEMFN='tem3d15.asc'
        DYEFN='dye3d15.asc'
        SEDFN='sed3d15.asc'
        SNDFN='snd3d15.asc'
        TOXFN='tox3d15.asc'
        UUUFN='uuu3d15.asc'
        VVVFN='vvv3d15.asc'
        WWWFN='www3d15.asc'
        CMPFN='cmp3d15.asc'
      END IF
      IF(NCALL3D.EQ.16) THEN
        SALFN='sal3d16.asc'
        TEMFN='tem3d16.asc'
        DYEFN='dye3d16.asc'
        SEDFN='sed3d16.asc'
        SNDFN='snd3d16.asc'
        TOXFN='tox3d16.asc'
        UUUFN='uuu3d16.asc'
        VVVFN='vvv3d16.asc'
        WWWFN='www3d16.asc'
        CMPFN='cmp3d16.asc'
      END IF
      IF(NCALL3D.EQ.17) THEN
        SALFN='sal3d17.asc'
        TEMFN='tem3d17.asc'
        DYEFN='dye3d17.asc'
        SEDFN='sed3d17.asc'
        SNDFN='snd3d17.asc'
        TOXFN='tox3d17.asc'
        UUUFN='uuu3d17.asc'
        VVVFN='vvv3d17.asc'
        WWWFN='www3d17.asc'
        CMPFN='cmp3d17.asc'
      END IF
      IF(NCALL3D.EQ.18) THEN
        SALFN='sal3d18.asc'
        TEMFN='tem3d18.asc'
        DYEFN='dye3d18.asc'
        SEDFN='sed3d18.asc'
        SNDFN='snd3d18.asc'
        TOXFN='tox3d18.asc'
        UUUFN='uuu3d18.asc'
        VVVFN='vvv3d18.asc'
        WWWFN='www3d18.asc'
        CMPFN='cmp3d18.asc'
      END IF
      IF(NCALL3D.EQ.19) THEN
        SALFN='sal3d19.asc'
        TEMFN='tem3d19.asc'
        DYEFN='dye3d19.asc'
        SEDFN='sed3d19.asc'
        SNDFN='snd3d19.asc'
        TOXFN='tox3d19.asc'
        UUUFN='uuu3d19.asc'
        VVVFN='vvv3d19.asc'
        WWWFN='www3d19.asc'
        CMPFN='cmp3d19.asc'
      END IF
      IF(NCALL3D.EQ.20) THEN
        SALFN='sal3d20.asc'
        TEMFN='tem3d20.asc'
        DYEFN='dye3d20.asc'
        SEDFN='sed3d20.asc'
        SNDFN='snd3d20.asc'
        TOXFN='tox3d20.asc'
        UUUFN='uuu3d20.asc'
        VVVFN='vvv3d20.asc'
        WWWFN='www3d20.asc'
        CMPFN='cmp3d20.asc'
      END IF
      IF(NCALL3D.EQ.21) THEN
        SALFN='sal3d21.asc'
        TEMFN='tem3d21.asc'
        DYEFN='dye3d21.asc'
        SEDFN='sed3d21.asc'
        SNDFN='snd3d21.asc'
        TOXFN='tox3d21.asc'
        UUUFN='uuu3d21.asc'
        VVVFN='vvv3d21.asc'
        WWWFN='www3d21.asc'
        CMPFN='cmp3d21.asc'
      END IF
      IF(NCALL3D.EQ.22) THEN
        SALFN='sal3d22.asc'
        TEMFN='tem3d22.asc'
        DYEFN='dye3d22.asc'
        SEDFN='sed3d22.asc'
        SNDFN='snd3d22.asc'
        TOXFN='tox3d22.asc'
        UUUFN='uuu3d22.asc'
        VVVFN='vvv3d22.asc'
        WWWFN='www3d22.asc'
        CMPFN='cmp3d22.asc'
      END IF
      IF(NCALL3D.EQ.23) THEN
        SALFN='sal3d23.asc'
        TEMFN='tem3d23.asc'
        DYEFN='dye3d23.asc'
        SEDFN='sed3d23.asc'
        SNDFN='snd3d23.asc'
        TOXFN='tox3d23.asc'
        UUUFN='uuu3d23.asc'
        VVVFN='vvv3d23.asc'
        WWWFN='www3d23.asc'
        CMPFN='cmp3d23.asc'
      END IF
      IF(NCALL3D.EQ.24) THEN
        SALFN='sal3d24.asc'
        TEMFN='tem3d24.asc'
        DYEFN='dye3d24.asc'
        SEDFN='sed3d24.asc'
        SNDFN='snd3d24.asc'
        TOXFN='tox3d24.asc'
        UUUFN='uuu3d24.asc'
        VVVFN='vvv3d24.asc'
        WWWFN='www3d24.asc'
        CMPFN='cmp3d24.asc'
      END IF
C
      IF (NCALL3D.EQ.1) THEN
        OPEN(50,FILE='out3d.dia',STATUS='UNKNOWN')
        CLOSE(50,STATUS='DELETE')
        OPEN(50,FILE='out3d.dia',STATUS='UNKNOWN')
        WRITE(50,520)IAD,JAD
        WRITE(50,530)NCALL3D
        DO KP=1,KPC
        WRITE(50,502)KP,ZZP(KP)
        END DO
       ELSE
        OPEN(50,FILE='out3d.dia',ACCESS='APPEND',STATUS='UNKNOWN')
        WRITE(50,530)NCALL3D
       END IF             
C 
      IF(IS3DSAL.GE.1) THEN
         ASALMAX=-99999999.
         ASALMIN= 99999999.
         DO K=1,KC
         DO L=2,LA
         TMPVAL=SAL(L,K)
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
         TMPVAL=DYE(L,K)
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
         TMPVAL=SEDT(L,K)
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
         TMPVAL=0.5*(U(L,K)+U(L+1,K))
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
         TMPVAL=0.5*(V(L,K)+V(LNC(L),K))
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
         TMPVAL=0.5*(W(L,K)+W(L,K-1))
     $       +GI*ZZ(K)*( DTI*(P(L)-P1(L))
     $       +0.5*(U(L+1,K)*(P(L+1)-P(L))*DXIU(L+1)
     $            +U(L,K)*(P(L)-P(L-1))*DXIU(L)
     $            +V(LN,K)*(P(LN)-P(L))*DYIV(LNC(L))
     $            +V(L,K)*(P(L)-P(LS))*DYIV(L)) )
     $       +0.5*(1.-ZZ(K))*( U(L+1,K)*(BELV(L+1)-BELV(L))*DXIU(L+1)
     $                        +U(L,K)*(BELV(L)-BELV(L-1))*DXIU(L)
     $                        +V(LN,K)*(BELV(LN)-BELV(L))*DYIV(LN)
     $                        +V(L,K)*(BELV(L)-BELV(LS))*DYIV(L) )
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
         TMPVAL=SNDT(L,K)
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
         TMPVAL=TOX(L,K,1)
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
      ZZPS=(ZZP(KP)-BELV(L))*HPI(L)
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
C     WRITE(50,500)NCALL3D
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
C     IF(NCALL3D.EQ.1) THEN
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
        TMP3D(K)=SAL(L,K)
        END DO
        SCALE3D=254./(SAL3DMA-SAL3DMI)
        DO KP=KPB(L),KPS(L)
        AP(KP)=0.
        DO K=1,KC
        AP(KP)=AP(KP)+APT(KP,K)*TMP3D(K)
C*DIAGNOSTIC
C       IF(NCALL3D.EQ.1) THEN
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
        TMP3D(K)=DYE(L,K)
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
        TMP3D(K)=SEDT(L,K)
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
        TMP3D(K)=0.5*(U(L,K)+U(L+1,K))
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
        TMP3D(K)=0.5*(V(L,K)+V(LN,K))
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
        TMP3D(K)=0.5*(W(L,K)+W(L,K-1))
     $         +GI*ZZ(K)*( DTI*(P(L)-P1(L))
     $         +0.5*(U(L+1,K)*(P(L+1)-P(L))*DXIU(L+1)
     $              +U(L,K)*(P(L)-P(L-1))*DXIU(L)
     $              +V(LN,K)*(P(LN)-P(L))*DYIV(LNC(L))
     $              +V(L,K)*(P(L)-P(LS))*DYIV(L)) )
     $         +0.5*(1.-ZZ(K))*( U(L+1,K)*(BELV(L+1)-BELV(L))*DXIU(L+1)
     $                          +U(L,K)*(BELV(L)-BELV(L-1))*DXIU(L)
     $                          +V(LN,K)*(BELV(LN)-BELV(L))*DYIV(LN)
     $                          +V(L,K)*(BELV(L)-BELV(LS))*DYIV(L) )
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
        TMP3D(K)=SNDT(L,K)
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
        TMP3D(K)=TOX(L,K,1)
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
  521 FORMAT('SALMAX = ',E12.4,'  SALMIN = ',E12.4/)
  522 FORMAT('TEMMAX = ',E12.4,'  TEMMIN = ',E12.4/)
  523 FORMAT('DYEMAX = ',E12.4,'  DYEMIN = ',E12.4/)
  524 FORMAT('SEDMAX = ',E12.4,'  SEDMIN = ',E12.4/)
  525 FORMAT('UUUMAX = ',E12.4,'  UUUMIN = ',E12.4/)
  526 FORMAT('VVVMAX = ',E12.4,'  VVVMIN = ',E12.4/)
  527 FORMAT('WWWMAX = ',E12.4,'  WWWMIN = ',E12.4/)
  528 FORMAT('SNDMAX = ',E12.4,'  SNDMIN = ',E12.4/)
  529 FORMAT('TOXMAX = ',E12.4,'  TOXMIN = ',E12.4/)
  530 FORMAT('NCALL3D = ',I5/)
  551 FORMAT(72F7.1)
  559 FORMAT(2I4,2X,72I2)
C
C----------------------------------------------------------------------C
C
      CLOSE (50)
      RETURN
      END 
