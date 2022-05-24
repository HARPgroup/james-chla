C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALHEAT(ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
C **  SUBROUTINE CALHEAT CALCULATES SURFACE AND INTERNAL HEAT SOURCES 
C **  AND SINKS IN THE HEAT (TEM) TRANSPORT EQUATION
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
      DIMENSION TEMOLD(LCM,KCM)
C
C**********************************************************************C
C
      DELT=DT2
      S3TL=1.0
      S2TL=0.0
      IF (ISTL.EQ.2) THEN
       DELT=DT
       S3TL=0.0
       S2TL=1.0
      END IF
      
 !     return
C
C**********************************************************************C
C
      IF (ISTOPT(2).EQ.1) THEN
C
      DO K=1,KCM
      DO L=2,LA
        TEMOLD(L,K)=TEM(L,K)
      END DO
      END DO
C
      DO L=2,LA
        TVAR3S(L)=0.
      END DO
C
C **  ADSORB SW SOLR RAD TO ALL LAYERS AND BED
C
      IF (IASWRAD.EQ.0) THEN
C
       DO L=2,LA
        SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/
     $              (1.+0.00412*TEM(L,KC))))
C    $    *(1+1.E-6*PATMT(L)*(4.5+0.0006*TEM(L,KC)*TEM(L,KC)))
        HDEP=MAX(HP(L),0.)
        WNDTMP=WINDST(L)
cx        IF(HP(L).LT.HWET) WNDTMP=0.
        TVAR1S(L,KC)=HDEP*TEM(L,KC)
     $  -(DELT*DZIC(KC))*( 1.312E-14*((TEM(L,KC)+273.)**4)
     $                  *(0.39-0.05*SQRT(VPA(L)))*(1.-.8*CLOUDT(L))
     $   +5.248E-14*((TEM(L,KC)+273.)**3)*(TEM(L,KC)-TATMT(L))
     $   +1.5*0.288E-6*WNDTMP*(TEM(L,KC)-TATMT(L))
     $   +1.5*0.445E-3*WNDTMP*(SVPW-VPA(L))/PATMT(L) )
     $  +(DELT*DZIC(KC))*0.2385E-6*SOLSWRT(L)*(
     $         FSWRATF*EXP(SWRATNF*HDEP*(Z(KC)-1.))
     $        +(1.-FSWRATF)*EXP(SWRATNS*HDEP*(Z(KC)-1.))
     $        -FSWRATF*EXP(SWRATNF*HDEP*(Z(KC-1)-1.))
     $        -(1.-FSWRATF)*EXP(SWRATNS*HDEP*(Z(KC-1)-1.)) )
       END DO
C
       IF(KC.GT.1) THEN
        DO K=1,KS
         DO L=2,LA
          HDEP=MAX(HP(L),0.)
          TVAR1S(L,K)=HDEP*TEM(L,K)
     $        +(DELT*DZIC(K))*0.2385E-6*SOLSWRT(L)*(
     $         FSWRATF*EXP(SWRATNF*HDEP*(Z(K)-1.))
     $        +(1.-FSWRATF)*EXP(SWRATNS*HDEP*(Z(K)-1.))
     $        -FSWRATF*EXP(SWRATNF*HDEP*(Z(K-1)-1.))
     $        -(1.-FSWRATF)*EXP(SWRATNS*HDEP*(Z(K-1)-1.)) )
         END DO
        END DO
       END IF
C
       DO L=2,LA
        UBED=0.5*( U(L,1)+U(L+1,1) )
        VBED=0.5*( V(L,1)+V(LNC(L),1) )
        USPD=SQRT( UBED*UBED+VBED*VBED )
        TMPVAL=(HTBED1*USPD+HTBED2)*(TEM(L,1)-TEMB(L))
        TVAR1S(L,1)=TVAR1S(L,1)-DELT*DZIC(1)*TMPVAL
        TEMB(L)=TEMB(L)+DELT*TMPVAL/DABEDT
     $        +(DELT/DABEDT)*0.2385E-6*SOLSWRT(L)*(
     $        +FSWRATF*EXP(SWRATNF*HDEP*(Z(0)-1.))
     $        +(1.-FSWRATF)*EXP(SWRATNS*HDEP*(Z(0)-1.)) )
       END DO
C
       DO K=1,KC
        DO L=2,LA
         IF(HP(L).GT.0.) TEM(L,K)=(HPI(L)*TVAR1S(L,K)+TEM(L,K))*0.5
        END DO
       END DO
C
      END IF
C
C **  ADSORB SW SOLR RAD TO SURFACE LAYER
C
      IF (IASWRAD.EQ.1) THEN
C
       DO L=2,LA
        SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/
     $              (1.+0.00412*TEM(L,KC))))
C    $    *(1+1.E-6*PATMT(L)*(4.5+0.0006*TEM(L,KC)*TEM(L,KC)))
        HDEP=MAX(HP(L),0.)
        WNDTMP=WINDST(L)
cx        IF(HP(L).LT.HWET) WNDTMP=0.
        TVAR1S(L,KC)=HDEP*TEM(L,KC)
     $  -(DELT*DZIC(KC))*( 1.312E-14*((TEM(L,KC)+273.)**4)
     $         *(0.39-0.05*SQRT(VPA(L)))*(1.-.8*CLOUDT(L))
     $   +5.248E-14*((TEM(L,KC)+273.)**3)*(TEM(L,KC)-TATMT(L))
     $   +CCNHTT(L)*0.288E-3*WNDTMP*(TEM(L,KC)-TATMT(L))
     $   +CLEVAP(L)*0.445*WNDTMP*(SVPW-VPA(L))/PATMT(L) )
     $  +(DELT*DZIC(KC))*0.2385E-6*SOLSWRT(L)
       END DO
C
       DO L=2,LA
        IF(HP(L).GT.0.) TEM(L,KC)=HPI(L)*TVAR1S(L,KC)
        if(TEM(L,KC)<=0.1) TEM(L,KC)=0.1
        if(TEM(L,KC)> 31) TEM(L,KC)=31
       END DO
C
      END IF
C
C **  END IF ISOPT(2) EQ 1
C
      END IF
C
  600 FORMAT(4I5,2E12.4)
C
C**********************************************************************C
C
      IF (ISTOPT(2).EQ.2) THEN
C
C ** IMPLEMENT EXTERNALLY SPECIFIED EQUILIBRIUM TEMPERATURE FROMULATION
C
      TMPKC=FLOAT(KC)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TEM(L,KC)=TEM(L,KC)-DELT*SOLSWRT(L)*HPI(L)*(TEM(L,KC)-TATMT(L))
       END DO
      END DO
C
      END IF
C
C**********************************************************************C
C
      IF (ISTOPT(2).EQ.3) THEN
C
C ** IMPLEMENT CONSTANT COEFFICIENT EQUILIBRIUM TEMPERATURE FROMULATION
C
      DTHEQT=DELT*HEQT*FLOAT(KC)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TEM(L,KC)=TEM(L,KC)-DTHEQT*HPI(L)*(TEM(L,KC)-TEMO)
       END DO
      END DO
C
      END IF
C
C**********************************************************************C
C
      RETURN
      END
