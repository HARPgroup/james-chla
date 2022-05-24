C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BEDINIT 
C
C**********************************************************************C
C
C **  SUBROUTINE BEDINIT INITIALIZES SEDIMENT AND TOXIC VARIABLES
C **  IT SEDIMENT BED FOR HOT AND COLD START CONDITIONS
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C**********************************************************************C
C
C **  HOT START INITIALIZATION
C
C**********************************************************************C
C
      IF(ISRESTI.NE.0) THEN   ! #1
C
C----------------------------------------------------------------------C
C
C **  SET POROSITY
C
      DO K=1,KB
      DO L=2,LA
       PORBED(L,K)=VDRBED(L,K)/(1.+VDRBED(L,K))
       PORBED1(L,K)=VDRBED1(L,K)/(1.+VDRBED1(L,K))
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  SET BULK DENSITY
C
      DO K=1,KB
      DO L=2,LA
       SEDBT(L,K)=0.
       SNDBT(L,K)=0.
      END DO
      END DO
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NS=1,NSND
       DO K=1,KB
       DO L=2,LA
        SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      DO K=1,KB
      DO L=2,LA
        BDENBED(L,K)=0.
        IF(HBED(L,K).GT.0.) BDENBED(L,K)=1000.*PORBED(L,K)
     &            +0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED(L,K)
      END DO
      END DO
C
      DO K=1,KB
      DO L=2,LA
       SEDBT(L,K)=0.
       SNDBT(L,K)=0.
      END DO
      END DO
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        SEDBT(L,K)=SEDBT(L,K)+SEDB1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NS=1,NSND
       DO K=1,KB
       DO L=2,LA
        SNDBT(L,K)=SNDBT(L,K)+SNDB1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      DO K=1,KB
      DO L=2,LA
        BDENBED1(L,K)=0.
        IF(HBED1(L,K).GT.0.) BDENBED1(L,K)=1000.*PORBED1(L,K)
     &              +0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED1(L,K)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  SET TOP BED LAYER
C
      DO L=2,LA
       KBT(L)=1
      END DO
C
      DO K=1,KB
      DO L=2,LA
       IF(HBED(L,K).GT.0.) KBT(L)=K
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES
C
      IF(ISTRAN(6).GE.1) THEN
      IF(IWRSP(1).EQ.0) THEN  !  efdc.inp, C40, IWRSP(1)=0
        DO K=1,KB
        DO L=2,LA
c         TAURS(L,K)=TAUR(1)    !@, Ji, 2/7/00
          TAURS(L,K)=xLTAUR(L,K)       ! surface layer uses spatially varying xltaur
          if(kb.ge.2.and.k.lt.kb) taurs(L,k)=taurc(k) ! deeper bed layers use uniform taurc on each layer
          TAURB(L,K)=1.E6                  ! bulk erosion shear stress, to make it never happens
c         WRSPS(L,K)=WRSPO(1)   !@
          WRSPS(L,K)=XLWRSPO(L)
          WRSPB(L,K)=0.0        ! => no mass erosion
        END DO
        END DO
      END IF
      IF(IWRSP(1).GE.1) THEN
        DO K=1,KB
        DO L=2,LA
          TAURS(L,K)=CSEDTAUS(BDENBED(L,K),IWRSP(1))
          TAURB(L,K)=CSEDTAUB(BDENBED(L,K),IWRSP(1))
          WRSPS(L,K)=CSEDRESS(BDENBED(L,K),IWRSP(1))
	    aatt=0.0
          WRSPB(L,K)=CSEDRESS(aatt,IWRSP(1))
        END DO
        END DO
      END IF
      END IF
C
C----------------------------------------------------------------------C
C
C **  SET SEDIMENT VOLUME FRACTIONS
C
      DO K=1,KB
      DO L=2,LA
       BEDLINIT(L,K)=0.
       BEDDINIT(L,K)=0.
      END DO
      END DO
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=SDEN(NS)*SEDB(L,K,NS)
        VFRBED1(L,K,NS)=SDEN(NS)*SEDB1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
      NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=SDEN(NS)*SNDB(L,K,NX)
        VFRBED1(L,K,NS)=SDEN(NS)*SNDB1(L,K,NX)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
        BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
      NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
        BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
        VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
      NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
        VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
       END DO
       END DO
      END DO
      END IF
C
C----------------------------------------------------------------------C
C
      RETURN
C
      END IF           ! #1
C
C**********************************************************************C
C**********************************************************************C
C
C **  COLD START INITIALIZATION: IBMECH=0
C
C**********************************************************************C
C
      IF(IBMECH.EQ.0) THEN    ! #2
C
C----------------------------------------------------------------------C
C
C **  SET POROSITY AND VOID RATIO
C
      DO K=1,KB
      DO L=2,LA
       PORBED(L,K)=BEDPORC
       PORBED1(L,K)=BEDPORC
       VDRBED(L,K)=PORBED(L,K)/(1.-PORBED(L,K))
       VDRBED1(L,K)=PORBED1(L,K)/(1.-PORBED1(L,K))
       HBED(L,K)=0.
       HBED1(L,K)=0.
       KBT(L)=1
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  UNIFORM SEDIMENT MASS PER UNIT AREA ALL CELLS, ALL BED LAYERS
C **  CALCULATE LAYER THICKNESS AND BULK DENSITY
C
      IF(ISEDINT.LE.1) THEN
C
      DO K=1,KB
      DO L=2,LA
       SEDBT(L,K)=0.
       SNDBT(L,K)=0.
      END DO
      END DO
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        HBED(L,K)=HBED(L,K)+SDEN(NS)*SEDB(L,K,NS)
        SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
       NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        HBED(L,K)=HBED(L,K)+SDEN(NS)*SNDB(L,K,NX)
        SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NX)
       END DO
       END DO
      END DO
      END IF
C
      DO K=1,KB
      DO L=2,LA
       HBED(L,K)=(1.+VDRBED(L,K))*HBED(L,K)
       IF(HBED(L,K).GT.0.) KBT(L)=K
      END DO
      END DO
C
      DO K=1,KB
      DO L=2,LA
        BDENBED(L,K)=0.
        IF(HBED(L,K).GT.0.) BDENBED(L,K)=1000.*PORBED(L,K)
     &              +0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED(L,K)
      END DO
      END DO
C
      DO K=1,KB
      DO L=2,LA
       HBED1(L,K)=HBED(L,K)
       BDENBED1(L,K)=BDENBED(L,K)
      END DO
      END DO
C
      END IF
C
C----------------------------------------------------------------------C
C
C **  NONUNIFORM SEDIMENT MASS PER UNIT AREA ALL CELLS, ALL BED LAYERS
C **  AND INITIAL CONDITIONS ARE IN MASS PER UNIT AREA
C **  CALCULATE LAYER THICKNESS AND BULK DENSITY
C
      IF(ISEDINT.GE.2) THEN
      IF(ISEDBINT.EQ.0) THEN
C
      DO K=1,KB
      DO L=2,LA
       SEDBT(L,K)=0.
       SNDBT(L,K)=0.
      END DO
      END DO
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        HBED(L,K)=HBED(L,K)+SDEN(NS)*SEDB(L,K,NS)
        SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
c       write(182,1821) K,L,HBED(L,K),SEDBT(L,K),sedb(L,K,NS)    ! Ji, 10/22/00
1821  format(2i6,999e12.4)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
       NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        HBED(L,K)=HBED(L,K)+SDEN(NS)*SNDB(L,K,NX)
        SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NX)
       END DO
       END DO
      END DO
      END IF
C
      DO K=1,KB
      DO L=2,LA
       HBED(L,K)=(1.+VDRBED(L,K))*HBED(L,K)
       IF(HBED(L,K).GT.0.) KBT(L)=K
      END DO
      END DO
C
      DO K=1,KB
      DO L=2,LA
        BDENBED(L,K)=0.
        IF(HBED(L,K).GT.0.) BDENBED(L,K)=1000.*PORBED(L,K)
     &              +0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED(L,K)
      END DO
      END DO
C
      DO K=1,KB
      DO L=2,LA
       HBED1(L,K)=HBED(L,K)
       BDENBED1(L,K)=BDENBED(L,K)
      END DO
      END DO
C
      END IF
      END IF
C
C----------------------------------------------------------------------C
C
C **  NONUNIFORM SEDIMENT MASS PER UNIT AREA ALL CELLS, ALL BED LAYERS
C **  AND INITIAL CONDITIONS ARE IN MASS FRACTION
C **  CALCULATE LAYER THICKNESS AND BULK DENSITY
C **  THIS OPTION REQUIRES INITIAL LAYER THICKNESSES
C
      IF(ISEDINT.GE.2) THEN
      IF(ISEDBINT.EQ.1) THEN
C
      IF(IBEDLAYU.EQ.1) THEN
        DO K=1,KB
        DO L=2,LA
         BEDLINIT(L,K)=0.1*BEDLINIT(L,K)
        END DO
        END DO
      END IF
C
      DO K=1,KB
      DO L=2,LA
       HBED(L,K)=BEDLINIT(L,K)
       HBED1(L,K)=BEDLINIT(L,K)
       IF(HBED(L,K).GT.0.) KBT(L)=K
      END DO
      END DO
C
      DO K=1,KB
      DO L=2,LA
       BDENBED(L,K)=0.
      END DO
      END DO
C
      DO K=1,KB
      DO L=2,LA
       BDENBED(L,K)=0.
      END DO
      END DO
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        BDENBED(L,K)=BDENBED(L,K)+1000.*SSG(NS)*SEDB(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
       NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        BDENBED(L,K)=BDENBED(L,K)+1000.*SSG(NS)*SNDB(L,K,NX)
       END DO
       END DO
      END DO
      END IF
C
      DO K=1,KB
      DO L=2,LA
        BDENBED(L,K)=1000.*PORBED(L,K)+(1.-PORBED(L,K))*BDENBED(L,K)
        BDENBED1(L,K)=BDENBED(L,K)
      END DO
      END DO
C
      DO K=1,KB
      DO L=2,LA
       SEDBT(L,K)=1000.*HBED(L,K)*(BDENBED(L,K)-1000.*PORBED(L,K))
      END DO
      END DO
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        SEDB(L,K,NS)=SEDB(L,K,NS)*SEDBT(L,K)
        SEDB1(L,K,NS)=SEDB(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
       DO K=1,KB
       DO L=2,LA
        SNDB(L,K,NX)=SNDB(L,K,NX)*SEDBT(L,K)
        SNDB1(L,K,NX)=SNDB(L,K,NX)
       END DO
       END DO
      END DO
      END IF
C
      END IF
      END IF
C
C----------------------------------------------------------------------C
C
C **  SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES
C
      IF(ISTRAN(6).GE.1) THEN
      IF(IWRSP(1).EQ.0) THEN
        DO K=1,KB
        DO L=2,LA
c         TAURS(L,K)=TAUR(1)    !@, Ji, 2/7/00
          TAURS(L,K)=xLTAUR(L,K)       ! surface layer uses spatially varying xltaur
      if(kb.ge.2.and.k.lt.kb) taurs(L,k)=taurc(k) ! deeper bed layers use uniform taurc on each layer
          TAURB(L,K)=1.E6
c         WRSPS(L,K)=WRSPO(1)   !@
          WRSPS(L,K)=XLWRSPO(L)
          WRSPB(L,K)=0.0
        END DO
        END DO
      END IF
      IF(IWRSP(1).GE.1) THEN
        DO K=1,KB
        DO L=2,LA
          TAURS(L,K)=CSEDTAUS(BDENBED(L,K),IWRSP(1))
          TAURB(L,K)=CSEDTAUB(BDENBED(L,K),IWRSP(1))
          WRSPS(L,K)=CSEDRESS(BDENBED(L,K),IWRSP(1)) 
		aatt=0.0   
          WRSPB(L,K)=CSEDRESS(aatt,IWRSP(1))
!          WRSPB(L,K)=CSEDRESS(BDENBED(L,K),IWRSP(1))
        END DO
        END DO
      END IF
      END IF
C
C----------------------------------------------------------------------C
C
C **  SET SEDIMENT VOLUME FRACTIONS
C
      DO K=1,KB
      DO L=2,LA
       BEDLINIT(L,K)=0.
       BEDDINIT(L,K)=0.
      END DO
      END DO
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=SDEN(NS)*SEDB(L,K,NS)
        VFRBED1(L,K,NS)=SDEN(NS)*SEDB1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
      NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=SDEN(NS)*SNDB(L,K,NX)
        VFRBED1(L,K,NS)=SDEN(NS)*SNDB1(L,K,NX)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
        BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
      NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
        BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
        VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
      NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
        VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
       END DO
       END DO
      END DO
      END IF
C
C----------------------------------------------------------------------C
C
      END IF          ! #2
C
C**********************************************************************C
C**********************************************************************C
C
C **  COLD START INITIALIZATION: IBMECH.GE.1
C
C**********************************************************************C
C
      IF(IBMECH.GE.1) THEN   ! #3
C
C----------------------------------------------------------------------C
C
C **  CONVERT AND INITIALIZE BED LAYER THICKNESS AND DEFINE
C **  INITIAL TOP LAYER
C
      IF(IBEDLAYU.EQ.1) THEN
        DO K=1,KB
        DO L=2,LA
         BEDLINIT(L,K)=0.1*BEDLINIT(L,K)
        END DO
        END DO
      END IF
C
      DO L=2,LA
       KBT(L)=1
      END DO
C
      DO K=1,KB
      DO L=2,LA
       HBED(L,K)=BEDLINIT(L,K)
       HBED1(L,K)=BEDLINIT(L,K)
       IF(HBED(L,K).GT.0.) KBT(L)=K
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  CONVERT AND INITIALIZE BED BULK DENSITY
C
      IF(IBEDBDNU.EQ.1) THEN
        DO K=1,KB
        DO L=2,LA
         BEDBINIT(L,K)=1000.*BEDBINIT(L,K)
        END DO
        END DO
      END IF
C
      DO K=1,KB
      DO L=2,LA
       BDENBED(L,K)=BEDBINIT(L,K)
       BDENBED1(L,K)=BEDBINIT(L,K)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  CONVERT AND DRY DENSITY OF BED
C **  IBEDDDNU=0,1 ACTUAL DRY DENSITY, =2  POROSITY, =3 VOID RATIO
C
      IF(IBEDDDNU.EQ.1) THEN
        DO K=1,KB
        DO L=2,LA
          BEDDINIT(L,K)=1000.*BEDDINIT(L,K)
        END DO
        END DO
      END IF
C
C----------------------------------------------------------------------C
C
C **  CALCULATE POROSITY AND VOID RATIO
C
      IF(IBEDDDNU.LE.1) THEN
        DO K=1,KB
        DO L=2,LA
         PORBED(L,K)=0.001*(BEDBINIT(L,K)-BEDDINIT(L,K))
         VDRBED(L,K)=PORBED(L,K)/(1.-PORBED(L,K))
        END DO
        END DO
      END IF
C
      IF(IBEDDDNU.EQ.2) THEN
        DO K=1,KB
        DO L=2,LA
         PORBED(L,K)=BEDDINIT(L,K)
         VDRBED(L,K)=PORBED(L,K)/(1.-PORBED(L,K))
        END DO
        END DO
      END IF
C
      IF(IBEDDDNU.EQ.3) THEN
        DO K=1,KB
        DO L=2,LA
         VDRBED(L,K)=BEDDINIT(L,K)
         PORBED(L,K)=VDRBED(L,K)/(1.+VDRBED(L,K))
        END DO
        END DO
      END IF
C
      DO K=1,KB
      DO L=2,LA
        VDRBED1(L,K)=VDRBED(L,K)
        PORBED1(L,K)=PORBED(L,K)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  INITIALIZE BED SEDIMENT FOR MASS FRACTION INPUT BY CACLUALTING
C **  AND STORING TOTAL MASS OF SED/AREA IN BEDDINIT(L,K)
C
      DO K=1,KB
      DO L=2,LA
        BEDDINIT(L,K)=HBED(L,K)*(BDENBED(L,K)-1000.*PORBED(L,K))
      END DO
      END DO
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       IF(ISEDBU(NS).EQ.1) THEN
         DO K=1,KB
         DO L=2,LA
          SEDB(L,K,NS)=SEDBINIT(L,K,NS)*BEDDINIT(L,K)
          SEDB1(L,K,NS)=SEDB(L,K,NS)
         END DO
         END DO
       END IF
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NS=1,NSND
       IF(ISNDBU(NS).EQ.1) THEN
         DO K=1,KB
         DO L=2,LA
          SNDB(L,K,NS)=SNDBINIT(L,K,NS)*BEDDINIT(L,K)
          SNDB1(L,K,NS)=SNDB(L,K,NS)
         END DO
         END DO
       END IF
      END DO
      END IF
C
C----------------------------------------------------------------------C
C
C **  SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES
C
      IF(ISTRAN(6).GE.1) THEN
      IF(IWRSP(1).EQ.0) THEN
        DO K=1,KB
        DO L=2,LA
          TAURS(L,K)=TAUR(1)
          TAURB(L,K)=1.E6
          WRSPS(L,K)=WRSPO(1)
          WRSPB(L,K)=0.0
        END DO
        END DO
      END IF
      IF(IWRSP(1).GE.1) THEN
        DO K=1,KB
        DO L=2,LA
          TAURS(L,K)=CSEDTAUS(BDENBED(L,K),IWRSP(1))
          TAURB(L,K)=CSEDTAUB(BDENBED(L,K),IWRSP(1))
          WRSPS(L,K)=CSEDRESS(BDENBED(L,K),IWRSP(1))
	    aatt=0.0
          WRSPB(L,K)=CSEDRESS(aatt,IWRSP(1))
        END DO
        END DO
      END IF
      END IF
C
C----------------------------------------------------------------------C
C
C **  SET SEDIMENT VOLUME FRACTIONS
C
      DO K=1,KB
      DO L=2,LA
       BEDLINIT(L,K)=0.
       BEDDINIT(L,K)=0.
      END DO
      END DO
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=SDEN(NS)*SEDB(L,K,NS)
        VFRBED1(L,K,NS)=SDEN(NS)*SEDB1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
      NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=SDEN(NS)*SNDB(L,K,NX)
        VFRBED1(L,K,NS)=SDEN(NS)*SNDB1(L,K,NX)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
        BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
      NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
        BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
        VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
      NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
        VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
       END DO
       END DO
      END DO
      END IF
C
C----------------------------------------------------------------------C
C
      END IF      ! #3
C
C**********************************************************************C
C**********************************************************************C
C
      OPEN(1,file='bedinit.dia')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,file='bedinit.dia')
      WRITE(1,102)
      DO L=2,LA
       WRITE(1,101)IL(L),JL(L),SEDB(L,1,1),SNDB(L,1,1),HBED(L,1),
     $             PORBED(L,1),VDRBED(L,1),BDENBED(L,1),TAURS(L,1)
      END DO
      CLOSE(1)
C
  101 FORMAT(2I5,8E12.4)
  102 FORMAT('   IL   JL    SEDB         SNDB        HBED',
     $       '      PORBED      VDRBED       BDENBED    TAURS')
c
C**********************************************************************C
C
      RETURN
      END
