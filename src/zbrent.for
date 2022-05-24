C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      FUNCTION ZBRENT(ISMERR)
C
C**********************************************************************C
C
C Using Brent's method, find the root of a FUNC SEDFLUX known to lie
C   betweenRMIN & RMAX within an accuracy of TOL (p. 253 in Numerical
c   Recipe).
C
C**********************************************************************C
C
      EXTERNAL SEDFLUX
C
      PARAMETER (IZMAX=100,EPS=3.0E-8,TOL=1.0E-5,
     $             RMIN=1.0E-4,RMAX=100.0)
C
      ISMERR = 0
      A = RMIN
      B = RMAX
      FA = SEDFLUX(A)
      FB = SEDFLUX(B)
      IF (FA*FB.GT.0.0) THEN
        ISMERR = 1
        RETURN
      END IF
C
      FC = FB
      DO II=1,IZMAX
        IF (FB*FC.GT.0.0) THEN
          C = A
          FC = FA
          D = B-A
          E = D
        END IF
        IF (ABS(FC).LT.ABS(FB)) THEN
          A = B
          B = C
          C = A
          FA = FB
          FB = FC
          FC = FA
        END IF
        TOL1 = 2.0*EPS*ABS(B) + 0.5*TOL
        XM = 0.5 * (C-B)
        IF (ABS(XM).LE.TOL1 .OR. FB.EQ.0.0) THEN
          ZBRENT = B
          RETURN
        END IF
        IF (ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S = FB / FA
          IF (A.EQ.C) THEN
            P = 2.0 * XM * S
            Q = 1.0 - S
           ELSE
            Q = FA / FC
            R = FB / FC
            P = S * (2.0*XM*Q*(Q-R) - (B-A)*(R-1.0))
            Q = (Q-1.0) * (R-1.0) * (S-1.0)
          END IF
          IF (P.GT.0.0) Q = -Q
          P = ABS(P)
          IF (2.0*P .LT. MIN(3.0*XM*Q-ABS(TOL1*Q), ABS(E*Q))) THEN
            E = D
            D = P / Q
           ELSE
            D = XM
            E = D
          END IF
         ELSE
          D = XM
          E = D
        END IF
C
        A = B
        FA = FB
        IF (ABS(D).GT.TOL1) THEN
          B = B + D
         ELSE
          B = B + SIGN(TOL1,XM)
        END IF
        FB = SEDFLUX(B)
      END DO
C
      ISMERR = 2
      ZBRENT = B
C
      RETURN
      END
