C 
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WQzero4(LA)
C
C**********************************************************************C
C M. Morton  02 Jun 1999
C Initializes the benthic flux arrays to 0.0
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
C
      DO LL=2,LA
c
c zero the benthic flux arrays:
c
        BFO2sum(LL)  = 0.0
        BFNH4sum(LL) = 0.0
        BFNO3sum(LL) = 0.0
        BFPO4sum(LL) = 0.0
        BFSADsum(LL) = 0.0
        BFCODsum(LL) = 0.0
      END DO
      timebf = 0.0
      nbfcnt = 0

c
      RETURN
      END
