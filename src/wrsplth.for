C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WRSPLTH
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE RSALPLTH WRITES FILES FOR RESIDUAL SCALAR FIELD 
C **  CONTOURING IN HORIZONTAL PLANES       
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      REAL DBS(10)
      CHARACTER*80 TITLE
C
C**********************************************************************C
C
      IF (JSWRPH.NE.1) GO TO 300
C
C----------------------------------------------------------------------C
C
      LINES=LA-1
      LEVELS=2
      LEVELSS=3
      LEVEL1=1
      DBS(1)=0.
      DBS(2)=99.
      DBS(3)=-99.
C
        TITLE=' WAVE REYNOLDS STRESS RXX '
        LUN=21
        OPEN(LUN,FILE='wvrxxh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='wvrxxh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELS
        WRITE (LUN,250)(DBS(L),L=1,LEVELS)
        CLOSE(LUN)
C
        TITLE=' WAVE REYNOLDS STRESS RYY '
        LUN=22
        OPEN(LUN,FILE='wvryyh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='wvryyh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELS
        WRITE (LUN,250)(DBS(L),L=1,LEVELS)
        CLOSE(LUN)
C
        TITLE=' WAVE REYNOLDS STRESS RXY '
        LUN=23
        OPEN(LUN,FILE='wvrxyh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='wvrxyh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELS
        WRITE (LUN,250)(DBS(L),L=1,LEVELS)
        CLOSE(LUN)
C
        TITLE=' WAVE DISSIPATION, (m/s)**3 '
        LUN=53
        OPEN(LUN,FILE='wvdisp.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='wvdisp',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVEL1
        WRITE (LUN,250)(DBS(L),L=1,LEVEL1)
        CLOSE(LUN)
C
        JSWRPH=0
C
C----------------------------------------------------------------------C
C
  300 CONTINUE
C
      OPEN(11,FILE='wvrxxh.col',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(12,FILE='wvryyh.col',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(13,FILE='wvrxyh.col',ACCESS='APPEND',STATUS='UNKNOWN')
      CLOSE(11,STATUS='DELETE')
      CLOSE(12,STATUS='DELETE')
      CLOSE(13,STATUS='DELETE')
      OPEN(11,FILE='wvrxxh.col',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(12,FILE='wvryyh.col',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(13,FILE='wvrxyh.col',ACCESS='APPEND',STATUS='UNKNOWN')
C
      OPEN(31,FILE='wvfxh.col',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(32,FILE='wvfyh.col',ACCESS='APPEND',STATUS='UNKNOWN')
      CLOSE(31,STATUS='DELETE')
      CLOSE(32,STATUS='DELETE')
      OPEN(31,FILE='wvfxh.col',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(32,FILE='wvfyh.col',ACCESS='APPEND',STATUS='UNKNOWN')
C
      OPEN(41,FILE='wvubh.col',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(42,FILE='wvvbh.col',ACCESS='APPEND',STATUS='UNKNOWN')
      CLOSE(41,STATUS='DELETE')
      CLOSE(42,STATUS='DELETE')
      OPEN(41,FILE='wvubh.col',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(42,FILE='wvvbh.col',ACCESS='APPEND',STATUS='UNKNOWN')
C
      OPEN(51,FILE='wvkx.col',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(52,FILE='wvky.col',ACCESS='APPEND',STATUS='UNKNOWN')
      CLOSE(51,STATUS='DELETE')
      CLOSE(52,STATUS='DELETE')
      OPEN(51,FILE='wvkx.col',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(52,FILE='wvky.col',ACCESS='APPEND',STATUS='UNKNOWN')
C
      OPEN(54,FILE='wvdisp.col',ACCESS='APPEND',STATUS='UNKNOWN')
      CLOSE(54,STATUS='DELETE')
      OPEN(54,FILE='wvdisp.col',ACCESS='APPEND',STATUS='UNKNOWN')
C
      OPEN(21,FILE='wvrxxh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(22,FILE='wvryyh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(23,FILE='wvrxyh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(53,FILE='wvdisp.out',STATUS='UNKNOWN')
C
      WRITE (21,100)N
      WRITE (22,100)N
      WRITE (23,100)N
C
      DO L=2,LA
      WVRSTOP=WVHUU(L,KC)+WVPP(L,KC)
      WVRSBOT=WVHUU(L,1)+WVPP(L,1)
      WRITE(11,200)IL(L),JL(L),WVRSTOP
      WRITE(21,200)IL(L),JL(L),DLON(L),DLAT(L),WVRSTOP,WVRSBOT
      WVRSTOP=WVHVV(L,KC)+WVPP(L,KC)
      WVRSBOT=WVHVV(L,1)+WVPP(L,1)
      WRITE(12,200)IL(L),JL(L),WVRSTOP
      WRITE(22,200)IL(L),JL(L),DLON(L),DLAT(L),WVRSTOP,WVRSBOT
      WRITE(13,200)IL(L),JL(L),WVHUV(L,KC)
      WRITE(23,200)IL(L),JL(L),DLON(L),DLAT(L),WVHUV(L,KC),
     $                                         WVHUV(L,1)
      WRITE(53,200)IL(L),JL(L),DLON(L),DLAT(L),WVDISP(L,KC)
      FXWTMP=0.5*(FXWAVE(L,KC)+FXWAVE(L+1   ,KC))
      FYWTMP=0.5*(FYWAVE(L,KC)+FYWAVE(LNC(L),KC))
      WRITE(31,200)IL(L),JL(L),FXWTMP
      WRITE(32,200)IL(L),JL(L),FYWTMP
C     FXWTMP=0.5*(WVUBU(L)+WVUBU(L+1   ))
C     FYWTMP=0.5*(WVVBV(L)+WVVBV(LSC(L)))
C     WRITE(41,200)IL(L),JL(L),FXWTMP
C     WRITE(42,200)IL(L),JL(L),FYWTMP
C     WRITE(51,200)IL(L),JL(L),WVACOS(L)
C     WRITE(52,200)IL(L),JL(L),WVASIN(L)
      WRITE(54,200)IL(L),JL(L),WVDISP(L,KC)
      END DO
C
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)
      CLOSE(31)
      CLOSE(32)
      CLOSE(41)
      CLOSE(42)
      CLOSE(51)
      CLOSE(52)
      CLOSE(53)
      CLOSE(54)
C
C**********************************************************************C
C
   99 FORMAT(A80)
  100 FORMAT(I10)
  101 FORMAT(2I10)
  200 FORMAT(2I5,1X,6E14.6)
  250 FORMAT(12E12.4)
cmrm  200 FORMAT(2I5,1X,1p,6E13.5) 
cmrm  250 FORMAT(1p,12E11.3)
C
C**********************************************************************C
C
      RETURN
      END
