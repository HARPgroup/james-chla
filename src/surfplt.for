C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SURFPLT
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE SURFPLT WRITES FILES TO CONTOUR FREE SURFACE 
C **  ELEVATION
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      CHARACTER*80 TITLE
C
C**********************************************************************C
C
      IF (JSPPH.NE.1) GO TO 300
C
      OPEN(10,FILE='surfcon.out',STATUS='UNKNOWN')
      CLOSE(10,STATUS='DELETE')
      OPEN(10,FILE='surfcon.out',STATUS='UNKNOWN')
      TITLE='INSTANTANEOUS SURFACE ELEVATION CONTOURS'
C
      LINES=LA-1
      LEVELS=1
      DBS=0.
C
      WRITE (10,99) TITLE
      WRITE (10,100)LINES,LEVELS
      WRITE (10,250)DBS
      CLOSE(10)
      JSPPH=0
C
  300 CONTINUE
C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON    
C
      OPEN(10,FILE='surfcon.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (10,100)N,TIME
C
      IF(IS1DCHAN.EQ.0) THEN
        DO L=2,LA
         SURFEL=GI*P(L)
         WRITE(10,200)IL(L),JL(L),DLON(L),DLAT(L),SURFEL,BELV(L),HP(L)
        END DO
       ELSE
        DO L=2,LA
         SURFEL=GI*P(L)
         WRITE(10,200)IL(L),JL(L),DLON(L),DLAT(L),SURFEL,BELV(L),HP(L),
     $                FADXP(L),FADYP(L)
        END DO
      END IF
C
      CLOSE(10)
C
C**********************************************************************C
C
   99 FORMAT(A80)
  100 FORMAT(I10,F12.4)
  101 FORMAT(2I10)
  200 FORMAT(2I5,1X,8E13.4)
  250 FORMAT(12E12.4)
cmrm  200 FORMAT(2I5,1X,1p,8E13.5) 
cmrm  250 FORMAT(1p,12E11.3)
C
C**********************************************************************C
C
      RETURN
      END
