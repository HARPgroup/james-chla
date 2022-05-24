C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALBUOY
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
C
C **  CALBUOY CALCULATES THE BUOYANCY USING MELLOR'S APPROXIMATION
C **  TO THE UNESCO EQUATION OF STATE (MELLOR, G.L., J. ATM AND OCEAN
C **  TECH, VOL 8, P 609)
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      IF(IBSC.EQ.1) GO TO 1000
C
      ISPCOR=1
C
C **  DENSITY RHOO AT P=0, S=0, AND T=TEMO
C
      RHOO=999.842594+6.793952E-2*TEMO-9.095290E-3*TEMO*TEMO
     $    +1.001685E-4*TEMO*TEMO*TEMO-1.120083E-6*TEMO*TEMO*TEMO*TEMO
     $    +6.536332E-9*TEMO*TEMO*TEMO*TEMO*TEMO
C
C **  DENSITY RHO=B(L,K) AT P=0, S=SAL(L,K), AND T=TEM(L,K)
C
      IF (ISTRAN(1).EQ.0.AND.ISTRAN(2).EQ.0) THEN
C
      DO K=1,KC
      DO L=2,LA
      B(L,K)=RHOO
      END DO
      END DO
C
      END IF
C
      IF (ISTRAN(1).GE.1.AND.ISTRAN(2).EQ.0) THEN
 !     IF (ISTRAN(1).GE.1.AND.ISTRAN(2).EQ.0) THEN
C
      DO K=1,KC
      DO L=2,LA
      SAL(L,K)=MAX(SAL(L,K),0.)
      SSTMP=SAL(L,K)
      TTMP=TEMO
      B(L,K)=RHOO+SSTMP*(0.824493-4.0899E-3*TTMP+7.6438E-5*TTMP*TTMP
     $   -8.2467E-7*TTMP*TTMP*TTMP+5.3875E-9*TTMP*TTMP*TTMP*TTMP)
     $   +SQRT(SSTMP)*SSTMP*(-5.72466E-3+1.0227E-4*TTMP
     $   -1.6546E-6*TTMP*TTMP)+4.8314E-4*SSTMP*SSTMP
      END DO
      END DO
C
      END IF
C
      IF (ISTRAN(1).EQ.0.AND.ISTRAN(2).GE.1) THEN
C
      if(mod(N,120).EQ.-2) then
       DO K=1,KC
        DO L=2,LA
         TVAR3S(L)=TEM(L,K)
        END DO
C
        DO NSM=1,1
        DO L=2,LA
        IF(LCT(L).GT.0.AND.LCT(L).LT.9) THEN
        I=IL(L)
        J=JL(L)
        HTN=TVAR3S(LNC(L))
        HTS=TVAR3S(LSC(L))
        HTE=TVAR3S(L+1)
        HTW=TVAR3S(L-1)
        IF (IJCT(I  ,J+1).EQ.9) HTN=TVAR3S(L)
        IF (IJCT(I  ,J-1).EQ.9) HTS=TVAR3S(L)
        IF (IJCT(I+1,J  ).EQ.9) HTE=TVAR3S(L)
        IF (IJCT(I-1,J  ).EQ.9) HTW=TVAR3S(L)
        TVAR3N(L)=(1.-WSMB)*TVAR3S(L)+0.25*WSMB*(HTN+HTS+HTE+HTW)
        END IF
       END DO
       END DO
C
       DO L=2,LA
        TEM(L,K)=TVAR3N(L)
       END DO
      ENDDO
      endif
C
!     if(mod(N,120).EQ.0 ) then
      DO K=1,KC
      DO L=2,LA
      TTMP=MAX(TEM(L,K),0.0);
      TTMP=MIN(TTMP,32.0);
      TTB=B(L,K)
      B(L,K)=999.842594+6.793952E-2*TTMP-9.095290E-3*TTMP*TTMP
     $    +1.001685E-4*TTMP*TTMP*TTMP-1.120083E-6*TTMP*TTMP*TTMP*TTMP
     $    +6.536332E-9*TTMP*TTMP*TTMP*TTMP*TTMP
  !    B(L,K)=(B(L,K)+TTB)/2
      END DO
      END DO
      end if
C
!      END IF
C
 !    IF (ISTRAN(1).GE.1.AND.ISTRAN(2).GE.10) THEN
      IF (ISTRAN(1).GE.1.AND.ISTRAN(2).GE.1) THEN
C
      DO K=1,KC
      DO L=2,LA
      SAL(L,K)=MAX(SAL(L,K),0.)
      SSTMP=SAL(L,K)
      TTMP=MAX(TEM(L,K),0.0);
      TTMP=MIN(TTMP,32.0);
      RHTMP=999.842594+6.793952E-2*TTMP-9.095290E-3*TTMP*TTMP
     $    +1.001685E-4*TTMP*TTMP*TTMP-1.120083E-6*TTMP*TTMP*TTMP*TTMP
     $    +6.536332E-9*TTMP*TTMP*TTMP*TTMP*TTMP
      B(L,K)=RHTMP+SSTMP*(0.824493-4.0899E-3*TTMP+7.6438E-5*TTMP*TTMP
     $   -8.2467E-7*TTMP*TTMP*TTMP+5.3875E-9*TTMP*TTMP*TTMP*TTMP)
     $   +SQRT(SSTMP)*SSTMP*(-5.72466E-3+1.0227E-4*TTMP
     $   -1.6546E-6*TTMP*TTMP)+4.8314E-4*SSTMP*SSTMP
      END DO
      END DO
C
      END IF
C
C **  APPLY MELLOR'S PRESSURE CORRECTION
C
      IF (ISPCOR.EQ.1) THEN
C
      DO K=1,KC
      DO L=2,LA
      PRES=RHOO*G*HP(L)*(1.-ZZ(K))*1.E-6
      CCON=1449.2+1.34*(SAL(L,K)-35.)+4.55*TEM(L,K)
     $    -0.045*TEM(L,K)*TEM(L,K)+0.00821*PRES+15.E-9*PRES*PRES
      TMP=PRES/(CCON*CCON)
      B(L,K)=B(L,K)+1.E+4*TMP*(1.-0.2*TMP)
      END DO
      END DO
C
      END IF
C
C **  REPLACE DENSITY B(L,K) WITH BUOYANCY B(L,K)
C
      DO K=1,KC
      DO L=2,LA
      B(L,K)=(B(L,K)/RHOO)-1.
      END DO
      END DO
C
C **  APPLY LOW SEDIMENT CONCENTRATION CORRECTION TO BUOYANCY
C
chong      IF (ISTRAN(6).GE.1) THEN
chongC
chong      DO NS=1,NSED         ! John & ji, decouple sed from density, 12/2/98
chong      DO K=1,KC
chong      DO L=2,LA
chong        B(L,K)=B(L,K)*(1.-SDEN(NS)*SED(L,K,NS))
chong     $        +(SSG(NS)-1.)*SDEN(NS)*SED(L,K,NS)
chong      END DO
chong      END DO
chong      END DO
chongC
chong      END IF
C
      IF (ISTRAN(7).GE.1) THEN
C
      DO NN=1,NSND   ! John & ji, decouple snd from density, 12/2/98
      NS=NN+NSED
      DO K=1,KC
      DO L=2,LA
        B(L,K)=B(L,K)*(1.-SDEN(NS)*SND(L,K,NN))
     $        +(SSG(NS)-1.)*SDEN(NS)*SND(L,K,NN)
      IF(ISNAN(TEM(L,K)))THEN
      TEM(L,K)=TEM1(L,K)
      ENDIF
      END DO
      END DO
      END DO
C
      END IF
C
      GO TO 2000
C
C**********************************************************************C
C
C     DENSITY AS A LINEAR FUNCTION OF SALINITY ONLY.  FOR DIAGNOSTIC
C     PURPOSES ONLY
C
 1000 CONTINUE
C
      DO K=1,KC
      DO L=2,LA
      B(L,K)=0.00075*SAL(L,K)
      END DO
      END DO
C
C**********************************************************************C
C
 2000 CONTINUE
C
      RETURN
      END
