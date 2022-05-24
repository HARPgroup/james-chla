c
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RSMICI(ISMTICI,LA)
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 10 april 1999
C
C**********************************************************************C
C
C Read in spatially and/or temporally varying ICs (unit INSMICI).
C
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
c
      DIMENSION XSMPON(NSMGM),XSMPOP(NSMGM),XSMPOC(NSMGM)
      CHARACTER TITLE(3)*79, ICICONT*3
C
      OPEN(1,FILE='wqsdici.inp',STATUS='OLD')
C
      OPEN(2,FILE='wq3d.out',STATUS='UNKNOWN',ACCESS='APPEND')
C
      IF (ISMTICI.EQ.0) THEN
        READ(1,50) (TITLE(M),M=1,3)
        WRITE(2,999)
        WRITE(2,50) (TITLE(M),M=1,3)
      END IF
C
      WRITE(2,60)'* Initial conditions at ', ISMTICI,
     *  ' th day from model start'
C
      READ(1,999)
      READ(1,50) TITLE(1)
      WRITE(2,50) TITLE(1)
C
      DO M=2,LA
cQUESTION        READ(INSMRST,90) I,J,(XSMPON(NW),NW=1,NSMG),
        READ(1,90) I,J,(XSMPON(NW),NW=1,NSMG),
     *    (XSMPOP(NW),NW=1,NSMG),(XSMPOC(NW),NW=1,NSMG),XSM1NH4,
     *    XSM2NH4,XSM2NO3,XSM2PO4,XSM2H2S,XSMPSI,XSM2SI,XSMBST,XSMT
        WRITE(2,90) I,J,(XSMPON(NW),NW=1,NSMG),
     *    (XSMPOP(NW),NW=1,NSMG),(XSMPOC(NW),NW=1,NSMG),XSM1NH4,
     *    XSM2NH4,XSM2NO3,XSM2PO4,XSM2H2S,XSMPSI,XSM2SI,XSMBST,XSMT
        IF (LIJW(I,J).LT.1 ) THEN
          PRINT*, 'i, j, line# = ', I,J,M-1
          STOP 'ERROR!! invalid (i,j) in FILE wqsdici.inp'
        END IF
        L=LIJW(I,J)
        DO MM=1,NSMG
          SMPON(L,MM)=XSMPON(MM)
          SMPOP(L,MM)=XSMPOP(MM)
          SMPOC(L,MM)=XSMPOC(MM)
        END DO
        SM1NH4(L)=XSM1NH4
        SM2NH4(L)=XSM2NH4
        SM2NO3(L)=XSM2NO3
        SM2PO4(L)=XSM2PO4
        SM2H2S(L)=XSM2H2S
        SMPSI(L) =XSMPSI
        SM2SI(L) =XSM2SI
        SMBST(L) =XSMBST
        SMT(L)   =XSMT
      END DO
C
      READ(1,52) ISMTICI, ICICONT
      WRITE(2,52) ISMTICI, ICICONT
C
      IF (ICICONT.EQ.'END') THEN
        CLOSE(1)
        ISMICI = 0
      END IF
C
      CLOSE(2)
C
  999 FORMAT(1X)
   50 FORMAT(A79)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)
   84 FORMAT(3I5, 20F8.4, F8.2)
   90 FORMAT(2I5, 18E12.4)
C
      RETURN
      END
