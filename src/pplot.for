C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C  
      SUBROUTINE PPLOT (IPT)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C  
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C  
      DIMENSION BNDU(51),BNDL(51)
C
      CHARACTER BLANK,ASTER,LET1(51),LET2(51),CHARY(ICM,JCM)
C
      DATA BLANK/' '/
      DATA ASTER/'*'/
C
      DATA LET1/'A',' ','B',' ','C',' ','D',' ','E',' ','F',' ',
     $         'G',' ','H',' ','I',' ','J',' ','K',' ','L',' ',
     $         'M',' ','N',' ','O',' ','P',' ','Q',' ','R',' ',
     $         'S',' ','T',' ','U',' ','V',' ','W',' ','X',' ',
     $         'Y',' ','Z'/                                    
C
      DATA LET2/'A','a','B','b','C','c','D','d','E','e','F','f',
     $         'G','g','H','h','I','i','J','j','K','k','L','l',
     $         'M','m','N','n','O','o','P','p','Q','q','R','r',
     $         'S','s','T','t','U','u','V','v','W','w','X','x',
     $         'Y','y','Z'/                                    
C
      PMAX=-99999.
      PMIN= 99999.
C
      DO L=2,LA
      IF (PAM(L) .GT. PMAX) PMAX=PAM(L)
      IF (PAM(L) .LT. PMIN) PMIN=PAM(L)
      END DO
C
      RNBAN=FLOAT(NBAN)
      PINV=(PMAX-PMIN)/RNBAN
C
C **  SET BOUND ARRAYS FOR REFERENCE TABLE
C
      BNDU(1)=PMAX
      BNDL(1)=PMAX-PINV
C
      DO M=2,NBAN
      MM=M-1
      BNDU(M)=BNDL(MM)
      BNDL(M)=BNDU(M)-PINV
      END DO
C
      IF (IPT.EQ.1) THEN
       DO M=1,NBAN
       WRITE (7,10) BNDU(M),LET1(M),BNDL(M)
       END DO
      ELSE
       DO M=1,NBAN
       WRITE (7,10) BNDU(M),LET2(M),BNDL(M)
       END DO
      END IF
C
   10 FORMAT (5X,E12.4,5X,A1,5X,E12.4)
C
C     WRITE (7,11)
   11 FORMAT (////)
      WRITE(7,12)
   12 FORMAT(1H1)
C
C **  LOAD CHARACTER ARRAY
C
      DO J=1,JC
      DO I=1,IC
      IF (IJCT(I,J) .NE. 9) CHARY(I,J)=BLANK
      IF (IJCT(I,J) .EQ. 9) CHARY(I,J)=ASTER
      END DO
      END DO
C
      BNDU(1)=BNDU(1)+1.
      BNDL(NBAN)=BNDL(NBAN)-1.
C
      DO L=2,LA
      I=IL(L)
      J=JL(L)
       IF (IPT.EQ.1) THEN
        DO M=1,NBAN
        IF(PAM(L).LT.BNDU(M).AND.PAM(L).GE.BNDL(M)) CHARY(I,J)=LET1(M)
        END DO
       ELSE
        DO M=1,NBAN
        IF(PAM(L).LT.BNDU(M).AND.PAM(L).GE.BNDL(M)) CHARY(I,J)=LET2(M)
        END DO
       END IF
      END DO
C
      DO JJ=1,JC,120
      JS=JJ
      JE=JJ+119
      IF(JE.GT.JC) JE=JC
      WRITE(7,22)JS,JE
      DO I=1,IC
      WRITE (7,20) I,(CHARY(I,J),J=JS,JE)
      END DO 
      END DO
C
C  20 FORMAT (10X,80A1)
   20 FORMAT (1X,I3,2X,120A1)
   22 FORMAT(1H1,'VALUES FOR J=',I5,2X,'TO J=',I5,//)
C
      RETURN
      END
