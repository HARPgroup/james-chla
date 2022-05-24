C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C 
      SUBROUTINE SHOWVAL2
C
C **  LAST MODIFIED BY Mike Morton ON 10 april 1999
C
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
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
        READ(1,*)NSHTYPE,NSHOWR,ICSHOW,JCSHOW, ISHPRT
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
c            WRITE(6,5)
            WRITE(6,6)
c            WRITE(6,5)
            WRITE(6,2)
           ELSE
            IF (NSHTYPE.EQ.3) THEN
            NSHOWC=0
            WRITE(6,2)
c            WRITE(6,3)
            WRITE(6,33)
            WRITE(6,44)
c            WRITE(6,5)
            WRITE(6,66)
c            WRITE(6,5)
            WRITE(6,2)
            END IF
            IF (NSHTYPE.EQ.4) THEN
            NSHOWC=0
            WRITE(6,2)
            WRITE(6,33)
            WRITE(6,44)
c            WRITE(6,5)
            WRITE(6,67)
c            WRITE(6,5)
            WRITE(6,2)
            END IF
            IF (NSHTYPE.EQ.5) THEN
            NSHOWC=0
            WRITE(6,2)
            WRITE(6,34)
            WRITE(6,44)
c            WRITE(6,5)
            WRITE(6,68)
c            WRITE(6,5)
            WRITE(6,2)
            END IF
            IF (NSHTYPE.GT.5) THEN
             NSHOWC=0
            WRITE(6,2)
            WRITE(6,300)
            WRITE(6,44)
            WRITE(6,65)
            WRITE(6,2)         
            ENDIF
          END IF
        END IF
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
       if (mod(n,ISHPRT) .eq. 0) then
        NSHOWC=NSHOWC+1
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
        IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN
        L=LIJ(ICSHOW,JCSHOW)
        LN=LNC(L)
        ZSURF=(HP(L)+BELV(L))*100.
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
        TEMKC=TEM(L,KC)
        TEMKB=TEM(L,1)
        IF(NSHTYPE.GT.5) THEN
         Alage1=WQVO(L,KC,1)/0.04
         Alage2=WQVO(L,KC,2)/0.04
         Alage3=WQVO(L,KC,3)/0.04
         Alage4=WQVO(L,KC,22)*10
         NHH4=WQVO(L,KC,14)*10
         NOO3=WQVO(L,KC,15)*10
        ENDIF
c        IZSURF=NINT(ZSURF)
c        IVELEKC=NINT(VELEKC)
c        IVELNKC=NINT(VELNKC)
c        ISALKC=NINT(SALKC)
c        ISEDKC=NINT(SEDKC)
c        ISNDKC=NINT(SNDKC)
c        ITEMKC=NINT(TEM(L,KC))
c        IAVKS=NINT(AVKS)
c        IABKS=NINT(ABKS)
c        IVELEKB=NINT(VELEKB)
c        IVELNKB=NINT(VELNKB)
c        ISALKB=NINT(SALKB)
c        ISEDKB=NINT(SEDKB)
c        ISNDKB=NINT(SNDKB)
c        ITEMKB=NINT(TEM(L,1))
c        IAVKB=NINT(AVKB)
c        IABKB=NINT(ABKB)
         IF (NSHTYPE.EQ.2) THEN
c          WRITE(6,7)N,ICSHOW,JCSHOW,IZSURF,
c     $                            IVELEKC,IVELNKC,ISALKC,IAVKS,IABKS,
c     $                            IVELEKB,IVELNKB,ISALKB,IAVKB,IABKB
c M. Morton writes real variables to PC screen:
          WRITE(6,7)N,ICSHOW,JCSHOW, ZSURF,
     +                             VELEKC, VELNKC, SALKC, AVKS, ABKS,
     +                             VELEKB, VELNKB, SALKB, AVKB, ABKB
         ELSE
         IF (NSHTYPE.EQ.3) THEN
c          WRITE(6,77)TIME,ICSHOW,JCSHOW,IZSURF,
c     $                            IVELEKC,IVELNKC,ISALKC,IAVKS,IABKS,
c     $                            IVELEKB,IVELNKB,ISALKB,IAVKB,IABKB
          WRITE(6,77) TIME, N, ICSHOW, JCSHOW, ZSURF,
     $                             VELEKC, VELNKC, SALKC,
     $                             VELEKB, VELNKB, SALKB
          END IF
          IF (NSHTYPE.EQ.4) THEN
c          WRITE(6,7)N,ICSHOW,JCSHOW,IZSURF,
c     $                            IVELEKC,IVELNKC,ISEDKC,IAVKS,IABKS,
c     $                            IVELEKB,IVELNKB,ISEDKB,IAVKB,IABKB
c          WRITE(6,7)N,ICSHOW,JCSHOW,IZSURF,
c     $                            IVELEKC,IVELNKC,ISNDKC,IAVKS,IABKS,
c     $                            IVELEKB,IVELNKB,ISNDKB,IAVKB,IABKB
          WRITE(6,77) time, N, ICSHOW, JCSHOW, ZSURF,
     $                             VELEKC, VELNKC, SEDKC,
     $                             VELEKB, VELNKB, SEDKB
c         WRITE(6,77) time, N, ICSHOW, JCSHOW, ZSURF,           ! Ji, 10/26/00
c    $                             VELEKC, VELNKC, SNDKC,
c    $                             VELEKB, VELNKB, SNDKB
          END IF
 
          IF (NSHTYPE.EQ.5) THEN
c          WRITE(6,77)TIME,ICSHOW,JCSHOW,IZSURF,
c     $                            IVELEKC,IVELNKC,ITEMKC,IAVKS,IABKS,
c     $                            IVELEKB,IVELNKB,ITEMKB,IAVKB,IABKB
          WRITE(6,77) TIME, N, ICSHOW, JCSHOW, ZSURF,
     $                             SEDKC, VELNKC, TEMKC,
     $                             SEDKB, VELNKB, TEMKB
          ENDIF
          
          IF (NSHTYPE.EQ.6) THEN
c          WRITE(6,77)TIME,ICSHOW,JCSHOW,IZSURF,
c     $                            IVELEKC,IVELNKC,ITEMKC,IAVKS,IABKS,
c     $                            IVELEKB,IVELNKB,ITEMKB,IAVKB,IABKB
          WRITE(6,77) TIME, N, ICSHOW, JCSHOW, ZSURF,
     $    Alage1,Alage2,Alage3,Alage4,NHH4,NOO3
          END IF
          IF (NSHTYPE.GT.7) THEN
c          WRITE(6,77)TIME,ICSHOW,JCSHOW,IZSURF,
c     $                            IVELEKC,IVELNKC,ITEMKC,IAVKS,IABKS,
c     $                            IVELEKB,IVELNKB,ITEMKB,IAVKB,IABKB
          L=LIJ(ICSHOW,JCSHOW)
          a_carb=(WQV(L,K1,4)+WQV(L,K1,5))
          a_agle=(WQV(L,K1,1)+WQV(L,K1,2)+WQV(L,K1,3))
          WRITE(6,77) TIME, N, ICSHOW, JCSHOW, ZSURF,
     $              TOX(L,KC,2),TOXB(L,2,2), a_carb,
     $              TOX(L,1,2), TOXB(L,1,2), a_agle
          END IF
        END IF
       end if 
      END IF 
C
C**********************************************************************C
C
    1 FORMAT (80X)
c Note: M. Morton changed format statements for cleaner appearence on PC:
    2 FORMAT(' --------------------------------------------',
     $       '-----------------------------------')
    3 FORMAT('      Time  I  J  ZSUR  VELE  VELN  SAL    AV ',
     $       '   AB  VELE  VELN  SAL    AV    AB')
    4 FORMAT('      Step              SURF  SURF SURF  SURF ',
     $       ' SURF  BOTT  BOTT BOTT  BOTT  BOTT')
    5 FORMAT(' ')
    6 FORMAT('                    cm  cm/s  cm/s  ppt cmsqs ',
     $       'cmsqs  cm/s  cm/s  ppt cmsqs cmsqs')
    7 FORMAT(' ',F10.3,1x,I2,1x,I2,1x,f5.1,1x,f5.0,1x,
     $       f5.0,1x,f4.1,1x,f5.1,1x,f5.1,1x,f5.0,1x,
     $       f5.0,1x,f4.1,1x,f5.1,1x,f5.1)

    8 FORMAT('| TIME |  I  |  J  |     FREE SURFACE ELEVATION     |',
     $       1X,' SURFACE SALINITY  |   BOTTOM SALINITY  |')
    9 FORMAT('| STEP |     |     |              (cm)              |',
     $       1X,'       (psu)       |       (psu)        |')
   10 FORMAT('|      |     |     |',I3,'                          ',
     $        I3,'| 0               ',I3,'| 0               ',I3,'|')
   11 FORMAT('|',I5,1X,'|',I4,1X,'|',I4,1X,'|',32A1,'|',20A1,'|',
     $       20A1,'|')
   33 FORMAT('     Model   Time   I   J    WSE   VELE   VELN    SED',
     $       '   VELE   VELN    SED')
   34 FORMAT('     Model   Time   I   J   ZSUR   VELE   VELN    TEM',
     +       '   VELE   VELN    TEM')
   44 FORMAT('      Time   Step                  SURF   SURF   SURF',
     $       '   BOTT   BOTT   BOTT')
   66 FORMAT('      days                    cm   cm/s   cm/s    ppt',
     $       '   cm/s   cm/s    ppt')
   65 FORMAT('      days                    cm   mg/l   degC   mg/l',
     $       '   mg/l   degC   mg/l') 
   67 FORMAT('      days                    cm   cm/s   cm/s   mg/L',
     +       '   cm/s   cm/s   mg/L')
   68 FORMAT('      days                    cm   cm/s   cm/s   degC',
     +       '   cm/s   cm/s   degC')
   77 FORMAT(' ',F9.2,1x,I6,1x,I3,1x,i3,1x,f6.1,1x,f6.1,1x,
     $       f6.1,1x,f6.1,1x,f6.1,1x,f6.1,1x,f6.1)
  300 FORMAT('    Model   Time   I   J    WSE   SED    TEMP    WQV',
     $       '   SED    TEMP    WQV')
c Note: the following format statements were John Hamrick's original code:
c   2 FORMAT('----------------------------------------------------',
c     $       '-------------------------------------------')
c    3 FORMAT('| TIME |  I  |  J  | ZSUR | VELE | VELN | SAL |  AV ',
c     $       1X,'|  AB  | VELE | VELN | SAL |  AV  |  AB  |')
c   33 FORMAT('| TIME |  I  |  J  | ZSUR | VELE | VELN | SED |  AV ',
c     $       1X,'|  AB  | VELE | VELN | SED |  AV  |  AB  |')
c   34 FORMAT('| TIME |  I  |  J  | ZSUR | VELE | VELN | TEM |  AV ',
c     $       1X,'|  AB  | VELE | VELN | TEM |  AV  |  AB  |')
c    4 FORMAT('| STEP |     |     |      | SURF | SURF | SUR | SURF',
c     $       1X,'| SURF | BOTT | BOTT | BOT | BOTT | BOTT |')
c    5 FORMAT('|      |     |     |      |      |      |     |     ',
c     $       1X,'|      |      |      |     |      |      |')
c    6 FORMAT('|      |     |     |  cm  | cm/s | cm/s | psu',
c     $       1X,'|cmsq/s|cmsq/s| cm/s | cm/s | psu |cmsq/s|cmsq/s|')
c    7 FORMAT('|',I6   ,'|',I4,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,
c     $       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|',I5,1X,
c     $       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|')
c    8 FORMAT('| TIME |  I  |  J  |     FREE SURFACE ELEVATION     |',
c     $       1X,' SURFACE SALINITY  |   BOTTOM SALINITY  |')
c    9 FORMAT('| STEP |     |     |              (cm)              |',
c     $       1X,'       (psu)       |       (psu)        |')
c   10 FORMAT('|      |     |     |',I3,'                          ',
c     $        I3,'| 0               ',I3,'| 0               ',I3,'|')
c   11 FORMAT('|',I5,1X,'|',I4,1X,'|',I4,1X,'|',32A1,'|',20A1,'|',
c     $       20A1,'|')
c   44 FORMAT('|      |     |     |      | SURF | SURF | SUR | SURF',
c     $       1X,'| SURF | BOTT | BOTT | BOT | BOTT | BOTT |')
c   66 FORMAT('| days |     |     |  cm  | cm/s | cm/s | psu',
c     $       1X,'|cmsq/s|cmsq/s| cm/s | cm/s | psu |cmsq/s|cmsq/s|')
c   67 FORMAT('|      |     |     |  cm  | cm/s | cm/s |mg/l',
c     $       1X,'|cmsq/s|cmsq/s| cm/s | cm/s |mg/l |cmsq/s|cmsq/s|')
c   68 FORMAT('| days |     |     |  cm  | cm/s | cm/s | d:C',
c     $       1X,'|cmsq/s|cmsq/s| cm/s | cm/s | d:C |cmsq/s|cmsq/s|')
c   77 FORMAT('|',F6.2 ,'|',I4,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,
c     $       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|',I5,1X,
c     $       '|',I5,1X,'|',I4,1X,'|',I5,1X,'|',I5,1X,'|')
C
C**********************************************************************C
C
      RETURN
      END
