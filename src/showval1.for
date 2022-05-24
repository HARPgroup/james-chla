C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SHOWVAL1
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
      CHARACTER BLANK,ASTER,CSURF(32),CSALS(20),CSALB(20)
C
      DATA BLANK/' '/
      DATA ASTER/'*'/
C
C**********************************************************************C
C
      IF (NSHOWC.EQ.NSHOWR) THEN
        OPEN(1,FILE='show.inp',STATUS='UNKNOWN')
         DO NSKIP=1,6
         READ(1,1)
         END DO
        READ(1,*)NSHTYPE,NSHOWR,ICSHOW,JCSHOW,ISHPRT
        READ(1,*)ZSSMIN,ZSSMAX,SSALMAX
        CLOSE(1)
        IF (NSHTYPE.EQ.1) THEN
          IZSMIN=NINT(ZSSMIN)
          IZSMAX=NINT(ZSSMAX)
          ISALMAX=NINT(SSALMAX)
          NSHOWC=0
          WRITE(6,2)   
          WRITE(6,8)   
          WRITE(6,9)   
          WRITE(6,10)IZSMIN,IZSMAX,ISALMAX,ISALMAX   
          WRITE(6,2)         
         ELSE
          IF (NSHTYPE.EQ.2) THEN
            NSHOWC=0        
            WRITE(6,2)
            WRITE(6,3)
            WRITE(6,4)
            WRITE(6,5)
            WRITE(6,6)
            WRITE(6,5)
            WRITE(6,2)
           ELSE
            IF (NSHTYPE.EQ.3) THEN
            NSHOWC=0        
            WRITE(6,2)
            WRITE(6,3)
            WRITE(6,44)
            WRITE(6,5)
            WRITE(6,66)
            WRITE(6,5)
            WRITE(6,2)
            END IF
            IF (NSHTYPE.EQ.4) THEN
            NSHOWC=0        
            WRITE(6,2)
            WRITE(6,33)
            WRITE(6,44)
            WRITE(6,5)
            WRITE(6,67)
            WRITE(6,5)
            WRITE(6,2)
            END IF
            IF (NSHTYPE.EQ.5) THEN
            NSHOWC=0        
            WRITE(6,2)
            WRITE(6,34)
            WRITE(6,44)
            WRITE(6,5)
            WRITE(6,68)
            WRITE(6,5)
            WRITE(6,2)
            END IF
          END IF
        END IF
      END IF
C     
      IMODTMP=MOD(N,ISHPRT)
      IF(IMODTMP.NE.0) THEN
        NSHOWC=NSHOWC+1
        RETURN
      END IF
C
      IF (NSHTYPE.EQ.1) THEN
        NSHOWC=NSHOWC+1
        DO M=1,32
        CSURF(M)=BLANK
        END DO
        DO M=1,20
        CSALS(M)=BLANK
        CSALB(M)=BLANK
        END DO
        L=LIJ(ICSHOW,JCSHOW)
        ZSURF=(HP(L)+BELV(L))*100.
        IF(IS1DCHAN.GT.0) ZSURF=HP(L)*100.
        ZSTMP=(31.*(ZSURF-ZSSMIN)/(ZSSMAX-ZSSMIN))+1.
        IZSTMP=NINT(ZSTMP)
        IF(IZSTMP.GT.32)IZSTMP=32
        IF(IZSTMP.LT.1)IZSTMP=1
        CSURF(IZSTMP)=ASTER
        SSTMP=(19.*SAL(L,KC)/SSALMAX)+1.
        SBTMP=(19.*SAL(L,1)/SSALMAX)+1.
        ISSTMP=NINT(SSTMP)
        ISBTMP=NINT(SBTMP)
        IF(ISSTMP.GT.20)ISSTMP=20
        IF(ISSTMP.LT.1)ISSTMP=1
        IF(ISBTMP.GT.20)ISBTMP=20
        IF(ISBTMP.LT.1)ISBTMP=1
        CSALS(ISSTMP)=ASTER
        CSALB(ISBTMP)=ASTER
        WRITE(6,11)N,ICSHOW,JCSHOW,(CSURF(I),I=1,32),(CSALS(J),J=1,20),
     $                                               (CSALB(K),K=1,20)
       ELSE   
        NSHOWC=NSHOWC+1
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
        L=LIJ(ICSHOW,JCSHOW)
        LN=LNC(L)
        ZSURF=(HP(L)+BELV(L))*100.
        IF(IS1DCHAN.GT.0) ZSURF=HP(L)*100.
        UTMP=0.5*STCUV(L)*(U(L+1,KC)+U(L,KC))*100.
        VTMP=0.5*STCUV(L)*(V(LN,KC)+V(L,KC))*100.
        VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
        VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
        UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))*100.
        VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))*100.
        VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
        VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
        AVKS=AV(L,KS)*10000.*HP(L)
        AVKB=AV(L,1)*10000.*HP(L)
        ABKS=AB(L,KS)*10000.*HP(L)
        ABKB=AB(L,1)*10000.*HP(L)
        SALKC=SAL(L,KC)
        SALKB=SAL(L,1)
c       SEDKC=SEDT(L,KC)/1000.
c       SEDKB=SEDT(L,1)/1000.
c       SNDKC=SNDT(L,KC)/1000.
c       SNDKB=SNDT(L,1)/1000.
        SEDKC=SEDT(L,KC)
        SEDKB=SEDT(L,1)
        SNDKC=SNDT(L,KC)
        SNDKB=SNDT(L,1)
C       SEDKC=SEDT(L,KC)*1000.
C       SEDKB=SEDT(L,1)*1000.
C       SNDKC=SNDT(L,KC)*1000.
C       SNDKB=SNDT(L,1)*1000.
        IZSURF=NINT(ZSURF)
        IVELEKC=NINT(VELEKC)
        IVELNKC=NINT(VELNKC)
        ISALKC=NINT(SALKC)
        ISEDKC=NINT(SEDKC)
        ISNDKC=NINT(SNDKC)
        ITEMKC=NINT(TEM(L,KC))
        IAVKS=NINT(AVKS)
        IABKS=NINT(ABKS)
        IVELEKB=NINT(VELEKB)
        IVELNKB=NINT(VELNKB)
        ISALKB=NINT(SALKB)
        ISEDKB=NINT(SEDKB)
        ISNDKB=NINT(SNDKB)
        ITEMKB=NINT(TEM(L,1))
        IAVKB=NINT(AVKB)
        IABKB=NINT(ABKB)
        IF (NSHTYPE.EQ.2) THEN
          WRITE(6,7)N,ICSHOW,JCSHOW,IZSURF,
     $                            IVELEKC,IVELNKC,ISALKC,IAVKS,IABKS,
     $                            IVELEKB,IVELNKB,ISALKB,IAVKB,IABKB
         ELSE
          IF (NSHTYPE.EQ.3) THEN
          WRITE(6,77)TIME,ICSHOW,JCSHOW,IZSURF,
     $                            IVELEKC,IVELNKC,ISALKC,IAVKS,IABKS,
     $                            IVELEKB,IVELNKB,ISALKB,IAVKB,IABKB
          END IF
          IF (NSHTYPE.EQ.4) THEN
C          WRITE(6,7)N,ICSHOW,JCSHOW,IZSURF,
C     $                            IVELEKC,IVELNKC,ISEDKC,IAVKS,IABKS,
C     $                            IVELEKB,IVELNKB,ISEDKB,IAVKB,IABKB
C          WRITE(6,7)N,ICSHOW,JCSHOW,IZSURF,
C     $                            IVELEKC,IVELNKC,ISNDKC,IAVKS,IABKS,
C     $                            IVELEKB,IVELNKB,ISNDKB,IAVKB,IABKB
          WRITE(6,77)TIME,ICSHOW,JCSHOW,IZSURF,
     $                            IVELEKC,IVELNKC,ISEDKC,IAVKS,IABKS,
     $                            IVELEKB,IVELNKB,ISEDKB,IAVKB,IABKB
          WRITE(6,79)ISNDKC,ISNDKB
          END IF
          IF (NSHTYPE.EQ.5) THEN
          WRITE(6,77)TIME,ICSHOW,JCSHOW,IZSURF,
     $                            IVELEKC,IVELNKC,ITEMKC,IAVKS,IABKS,
     $                            IVELEKB,IVELNKB,ITEMKB,IAVKB,IABKB
          END IF
        END IF
      END IF
C
C**********************************************************************C
C
    1 FORMAT (80X)
    2 FORMAT('----------------------------------------------------',
     $       '-------------------------------------------')
    3 FORMAT('| TIME |  I  |  J  | ZSUR | VELE | VELN | SAL |  AV ',
     $       1X,'|  AB  | VELE | VELN | SAL |  AV  |  AB  |')
   33 FORMAT('| TIME |  I  |  J  | ZSUR | VELE | VELN | SED |  AV ',
     $       1X,'|  AB  | VELE | VELN | SED |  AV  |  AB  |')
   34 FORMAT('| TIME |  I  |  J  | ZSUR | VELE | VELN | TEM |  AV ',
     $       1X,'|  AB  | VELE | VELN | TEM |  AV  |  AB  |')
    4 FORMAT('| STEP |     |     |      | SURF | SURF | SUR | SURF',
     $       1X,'| SURF | BOTT | BOTT | BOT | BOTT | BOTT |')
    5 FORMAT('|      |     |     |      |      |      |     |     ',
     $       1X,'|      |      |      |     |      |      |')
    6 FORMAT('|      |     |     |  cm  | cm/s | cm/s | psu',
     $       1X,'|cmsq/s|cmsq/s| cm/s | cm/s | psu |cmsq/s|cmsq/s|')
    7 FORMAT('|',I6   ,'|',I4,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,
     $       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|',I5,1X,
     $       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|')
    8 FORMAT('| TIME |  I  |  J  |     FREE SURFACE ELEVATION     |',
     $       1X,' SURFACE SALINITY  |   BOTTOM SALINITY  |')
    9 FORMAT('| STEP |     |     |              (cm)              |',
     $       1X,'       (psu)       |       (psu)        |')
   10 FORMAT('|      |     |     |',I3,'                          ',
     $        I3,'| 0               ',I3,'| 0               ',I3,'|')
   11 FORMAT('|',I5,1X,'|',I4,1X,'|',I4,1X,'|',32A1,'|',20A1,'|',
     $       20A1,'|')
   44 FORMAT('|      |     |     |      | SURF | SURF | SUR | SURF',
     $       1X,'| SURF | BOTT | BOTT | BOT | BOTT | BOTT |')
   66 FORMAT('| days |     |     |  cm  | cm/s | cm/s | psu',
     $       1X,'|cmsq/s|cmsq/s| cm/s | cm/s | psu |cmsq/s|cmsq/s|')
   67 FORMAT('|      |     |     |  cm  | cm/s | cm/s |mg/l',
     $       1X,'|cmsq/s|cmsq/s| cm/s | cm/s |mg/l |cmsq/s|cmsq/s|')
   68 FORMAT('| days |     |     |  cm  | cm/s | cm/s | d:C',
     $       1X,'|cmsq/s|cmsq/s| cm/s | cm/s | d:C |cmsq/s|cmsq/s|')
   77 FORMAT('|',F6.1 ,'|',I4,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,
     $       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|',I5,1X,
     $       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|')
   79 FORMAT('|',6X ,'|',4X,1X,'|',4X,1X,'|',5X,1X,'|',5X,1X,
     $       '|',5X,1X,'|',I4,1X,'|',5X,1X,'|',5X,1X,'|',5X,1X,
     $       '|',5X,1X,'|',I4,1X,'|',5X,1X,'|',5X,1X,'|')
C
C**********************************************************************C
C
      RETURN
      END
