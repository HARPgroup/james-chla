C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WQzero2(LA,KC)
C
C**********************************************************************C
C M. Morton  12 Apr 1999
C Initializes the diurnal DO summation arrays:
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
          SODsum(LL,K) = 0.0
          RKAsum(LL,K) = 0.0
          SWQsum(LL,K) = 0.0
          TEMsum(LL,K) = 0.0
          DZsum(LL,K)  = 0.0
          DOOsum(LL,K) = 0.0
          DOSsum(LL,K) = 0.0
          CYAsum(LL,K) = 0.0
          DIAsum(LL,K) = 0.0
          GRNsum(LL,K) = 0.0
          XMACsum(LL,K) = 0.0
          DO I=1,4
            RESPsum(LL,K,I) = 0.0
            PRODsum(LL,K,I) = 0.0
          END DO
        END DO
      END DO
      timesum2 = 0.0
      ndocnt = 0
c
      RETURN
      END
