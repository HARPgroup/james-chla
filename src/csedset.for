C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      REAL FUNCTION CSEDSET(SED,SHEAR,IOPT)
C
C **  CALCULATES CONCENTRATION DEPENDENT SETTLING VELOCITY OF COHESIVE 
C **  SEDIMENT
C
C **  IOPT=1  BASED ON 
C **
C **  HWANG, K. N., AND A. J. MEHTA, 1989: FINE SEDIMENT DRODIBILITY
C **  IN LAKE OKEECHOBEE FLORIDA. COASTAL AND OCEANOGRAPHIC ENGINEERING
C **  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661
C
      IF(IOPT.EQ.1) THEN
       TMP=SED/2000.
       TMP=LOG10(TMP)
       TMP=-16.*TMP*TMP/9.
       TMP=10.**TMP
       CSEDSET=8.E-4*TMP
      END IF
C
C **  IOPT=2  BASED ON 
C **
C **  SHRESTA, P. L., AND G. T. ORLOB, 1996: MULTIPHASE DISTRIBUTION   
C **  OF COHESIVE SEDIMENTS AND HEAVY METALS IN ESTUARINE SYSTEMS,   
C **  J. ENVIRONMENTAL ENGINEERING, 122, 730-740. 
C
C     CONVERT SED CONCENTRATION FORM MG/L TO GM/L BY X 1.E-3
C     SHEAR IN 1/SEC
C     CONVERT SETTLING FROM M/HOUR TO M/SEC BY / 3600
C
      IF(IOPT.EQ.2) THEN
       SED=1.E-3*SED
       RNG=0.11075+0.0386*SHEAR
       BG=EXP(-4.20706+0.1465*SHEAR)
       WTMP=BG*( SED**RNG )
       CSEDSET=WTMP/3600.
      END IF
C
C **  IOPT=3  BASED ON 
C **
C **  ZIEGLER, C. K., AND B. S. NESBIT, 1995: LONG-TERM SIMULAITON   
C **  OF FINE GRAIN SEDIMENT TRANSPORT IN LARGE RESERVOIR,   
C **  J. HYDRAULIC ENGINEERING, 121, 773-781. 
C
C     CONVERT SED CONCENTRATION FORM MG/L (GM/M**3) TO GM/CM**3 BY
C      X 1.E-6
C     SHEAR IS STRESS IN DY/CM**2. CONVERT KINEMATIC STRESS (M/S)**2
C      TO N/M**2 BY X 1.E+3, CONVERT N/M**2 TO DY/CM**2 BY X 10
C     CONVERT SETTLING FROM CM/SEC TO M/SEC BY X 0.01
C
      IF(IOPT.EQ.3) THEN
       SED=1.E-6*SED
       GG=1.E4*SHEAR
       CG=GG*SED
       CG=MAX(CG,7.51E-6)
       BD2=-0.4-0.25*LOG10(CG-7.5E-6)
       CON=9.6E-4*( (1.E-8)**BD2 )
       VAL=CG**(-0.85-BD2)
       CSEDSET=0.01*CON*VAL
      END IF
C
C **  IOPT=4  BASED ON 
C **
      IF(IOPT.EQ.4) THEN
       GG=1.E4*SHEAR
       TMP=GG*SED
       CSEDSET=8.E-5
       IF(TMP.LT.40.0) CSEDSET=1.51E-5*(TMP**0.45)
       IF(TMP.GT.400.0) CSEDSET=0.893E-6*(TMP**0.75)
      END IF
C
C **  IOPT=5  BASED ON 
C **
      IF(IOPT.EQ.5) THEN
       GG=1.E4*SHEAR
       TMP=GG*SED
       CSEDSET=3.82E-5*(TMP**0.12)
      END IF
C
C **  IOPT=6  BASED ON 
C **
      IF(IOPT.EQ.6) THEN
       GG=1.E4*SHEAR
       TMP=GG*SED
       IF(TMP.LT.100.0) CSEDSET=2.*1.16E-5*(TMP**0.5)
       IF(TMP.GE.100.0) CSEDSET=2.*1.84E-5*(TMP**0.4)
      END IF
C
C **  IOPT=6  BASED ON 
C **
      IF(IOPT.EQ.7) THEN
       TMP=SHEAR*SED
       CSEDSET=0.0052*(TMP**0.470138)
      END IF
C
      RETURN
      END