C 
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WQzero3(LA,KC)
C
C**********************************************************************C
C M. Morton  29 Apr 1999
C Initializes the limitation and D.O. component analysis arrays:
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
C
      DO LL=2,LA
        DO K=1,KC
c
c zero the Diurnal DO variables:
c
          xLimIc(LL,K) = 0.0
          xLimId(LL,K) = 0.0
          xLimIg(LL,K) = 0.0
          xLimIm(LL,K) = 0.0
          xLimNc(LL,K) = 0.0
          xLimNd(LL,K) = 0.0
          xLimNg(LL,K) = 0.0
          xLimNm(LL,K) = 0.0
          xLimPc(LL,K) = 0.0
          xLimPd(LL,K) = 0.0
          xLimPg(LL,K) = 0.0
          xLimPm(LL,K) = 0.0
          xLimTc(LL,K) = 0.0
          xLimTd(LL,K) = 0.0
          xLimTg(LL,K) = 0.0
          xLimTm(LL,K) = 0.0
          xDOsat(LL,K) = 0.0
          xDOdef(LL,K) = 0.0
          xDOpsl(LL,K) = 0.0
          xDOsod(LL,K) = 0.0
          xDOkar(LL,K) = 0.0
          xDOdoc(LL,K) = 0.0
          xDOnit(LL,K) = 0.0
          xDOcod(LL,K) = 0.0
          xDOppB(LL,K) = 0.0
          xDOrrB(LL,K) = 0.0
          xDOppM(LL,K) = 0.0
          xDOrrM(LL,K) = 0.0
          xDOtrn(LL,K) = 0.0
          xDOall(LL,K) = 0.0
          xDOdz(LL,K)  = 0.0
          xDOacoef(LL,K)=0.0
        END DO
      END DO
      timesum3 = 0.0
      NLIM = 0
c
      RETURN
      END
