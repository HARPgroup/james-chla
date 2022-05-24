C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      REAL FUNCTION CSEDRESS(DENBULK,IOPT)
C
C **  CALCULATES SURFACE EROSION RATE OF COHESIVE 
C **  SEDIMENT AS A FUNCTION OF BED BULK DENSITY
C
C **  IOPT=1  BASED ON 
C **
C **  HWANG, K. N., AND A. J. MEHTA, 1989: FINE SEDIMENT DRODIBILITY
C **  IN LAKE OKEECHOBEE FLORIDA. COASTAL AND OCEANOGRAPHIC ENGINEERING
C **  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661
C
C **  IOPT=2  BASED ON 
C **
C **  HWANG, K. N., AND A. J. MEHTA, 1989: FINE SEDIMENT DRODIBILITY
C **  IN LAKE OKEECHOBEE FLORIDA. COASTAL AND OCEANOGRAPHIC ENGINEERING
C **  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661
C
      IF(IOPT.EQ.1) THEN
       DENBULK=DENBULK*0.001
       IF(DENBULK.LE.1.065) THEN
         CSEDRESS=0.62
        ELSE
         TMP=0.198/(DENBULK-1.0023)
         TMP=EXP(TMP)
         CSEDRESS=6.4E-4*(10.**TMP)
       END IF
      END IF
C
      IF(IOPT.EQ.2) THEN
       DENBULK=0.001*DENBULK
       IF(DENBULK.LE.1.065) THEN
         CSEDRESS=0.62
        ELSE
         TMP=0.198/(DENBULK-1.0023)
         TMP=EXP(TMP)
         CSEDRESS=6.4E-4*(10.**TMP)
       END IF
      END IF
C
      RETURN
      END
