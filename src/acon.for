C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE ACON(ITVAL)

C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
      PARAMETER (NJELM=2,NATDM=1)
C
C**********************************************************************C
C
c      DIMENSION ZAD(KCM,NATDM),  TAD(KCM,NATDM), UAGD(KCM,NATDM), 
c     $         VAGD(KCM,NATDM), WAGD(KCM,NATDM),SALAD(KCM,NATDM),
c     $        TEMAD(KCM,NATDM),DYEAD(KCM,NATDM),SFLAD(KCM,NATDM),
c     $        TOXAD(KCM,NTXM,NATDM),
c     $        SEDAD(KCM,NSCM,NATDM),SNDAD(KCM,NSNM,NATDM)
C
c      DIMENSION TOXA(NTXM),SEDA(NSCM),SNDA(NSNM)
C
c      COMMON/ADATA/ ZAD,  TAD, UAGD,VAGD, WAGD,SALAD,
c     $        TEMAD,DYEAD,SFLAD,TOXAD,SEDAD,SNDAD
C
c      COMMON/AVAL/ NAZD,NATD,ZVAL,TVAL,UAG, VAG, WAG,SALA,
c     $            TEMA,DYEA,SFLA,TOXA,SEDA,SNDA
C
C**********************************************************************C
C
       NTOX1=NTOX
      IF(NPCB.GT.0)NTOX1=NPCB
      IF(ZVAL.LT.ZAD(1,ITVAL)) THEN
        UAG=UAGD(1,ITVAL)
        VAG=VAGD(1,ITVAL)
        WAG=WAGD(1,ITVAL)
        SALA=SALAD(1,ITVAL)
        TEMA=TEMAD(1,ITVAL)
CJHFIX-REMOVE        DYEA=SALAD(1,ITVAL)
CJHFIX-ADD THE FOLLOWING 1 LINE(S)
        DYEA=DYEAD(1,ITVAL)
        SFLA=SFLAD(1,ITVAL)
        DO NT=1,NTOX1
         TOXA(NT)=TOXAD(1,NT,ITVAL)
        END DO
        DO NS=1,NSED
         SEDA(NS)=SEDAD(1,NS,ITVAL)
        END DO
        DO NX=1,NSND
         SNDA(NX)=SNDAD(1,NX,ITVAL)
        END DO
        RETURN
      END IF
C
      IF(ZVAL.GE.ZAD(NAZD,ITVAL)) THEN
        UAG=UAGD(NAZD,ITVAL)
        VAG=VAGD(NAZD,ITVAL)
        WAG=WAGD(NAZD,ITVAL)
        SALA=SALAD(NAZD,ITVAL)
        TEMA=TEMAD(NAZD,ITVAL)
        DYEA=DYEAD(NAZD,ITVAL)
        SFLA=SFLAD(NAZD,ITVAL)
        DO NT=1,NTOX1
         TOXA(NT)=TOXAD(NAZD,NT,ITVAL)
        END DO
        DO NS=1,NSED
         SEDA(NS)=SEDAD(NAZD,NS,ITVAL)
        END DO
        DO NX=1,NSND
         SNDA(NX)=SNDAD(NAZD,NX,ITVAL)
        END DO
        RETURN
      END IF
C
      NZ=1        
 1000 CONTINUE
      NZP=NZ+1
      IF(ZVAL.GE.ZAD(NZ,ITVAL).AND.ZVAL.LT.ZAD(NZP,ITVAL)) THEN
        DZI=1./(ZAD(NZP,ITVAL)-ZAD(NZ,ITVAL))
        WTNZ=DZI*(ZAD(NZP,ITVAL)-ZVAL)
        WTNZP=DZI*(ZVAL-ZAD(NZ,ITVAL))
        UAG=WTNZ*UAGD(NZ,ITVAL)+WTNZP*UAGD(NZP,ITVAL)
        VAG=WTNZ*VAGD(NZ,ITVAL)+WTNZP*VAGD(NZP,ITVAL)
        WAG=WTNZ*WAGD(NZ,ITVAL)+WTNZP*WAGD(NZP,ITVAL)
        SALA=WTNZ*SALAD(NZ,ITVAL)+WTNZP*SALAD(NZP,ITVAL)
        TEMA=WTNZ*TEMAD(NZ,ITVAL)+WTNZP*TEMAD(NZP,ITVAL)
        DYEA=WTNZ*DYEAD(NZ,ITVAL)+WTNZP*DYEAD(NZP,ITVAL)
        SFLA=WTNZ*SFLAD(NZ,ITVAL)+WTNZP*SFLAD(NZP,ITVAL)
        DO NT=1,NTOX1
         TOXA(NT)=WTNZ*TOXAD(NZ,NT,ITVAL)+WTNZP*TOXAD(NZP,NT,ITVAL)
        END DO
        DO NS=1,NSED
         SEDA(NS)=WTNZ*SEDAD(NZ,NS,ITVAL)+WTNZP*SEDAD(NZP,NS,ITVAL)
        END DO
        DO NX=1,NSND
         SNDA(NX)=WTNZ*SNDAD(NZ,NX,ITVAL)+WTNZP*SNDAD(NZP,NX,ITVAL)
        END DO
        RETURN
       ELSE
        NZ=NZ+1
        GO TO 1000
      END IF
C
C**********************************************************************C
C
      RETURN
      END
