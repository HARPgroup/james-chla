C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RWQATM(LA,KC,DT,N,TBEGIN,TCON)
C
C**********************************************************************C
C
C mrm **  added by Mike Morton  8 June 1998
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 10 april 1999
C
C**********************************************************************C
C ** Computes wet atmospheric deposition using constant concentrations
c ** for the 21 state variables multiplied by the rainfall flow rate
c ** entering each grid cell.  Computed loads are in g/day.
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn' 
C
C**********************************************************************C
c  cv2 = conversion to get units of g/day
c  wqatm(nw) has units of mg/L
c  raint(L) has units of m/sec
c  dxyp(L) has units of m2
c  WQATML(L,kc,nw) has units of g/day
      cv2=86400.0
!      DO NW=1,NWQV
!        DO L=2,LA
!         WQATML(L,kc,nw)=wqatm(nw)*raintwq(L)*dxypwq(L)*cv2
!        END DO
!      END DO
C
        L0 = LIJW(ICPSL(IWQPS), JCPSL(IWQPS)) ! Last point source
      DO NW=1,NWQV
        DO L=2,LA
         WQATML(L,KC,NW)=WQWPSL(L0,KC,NW)*dxypwq(L)*KC  ! Input used equal weight for each layer  Otherwise, don't time 8.0
         WQWPSL(L0,KC,NW)=0
        END DO
      END DO
      

      
      IF(ITNWQ.EQ.0) THEN
C
       OPEN(1,FILE='wqatm.dia',STATUS='UNKNOWN')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='wqatm.dia',STATUS='UNKNOWN')
C
       TIME = (DT*FLOAT(N) + TBEGIN*TCON)/86400.0
       WRITE(1,*) N,TIME
C
       DO L=2,LA
         WRITE(1,110) ILW(L),JLW(L),(WQATML(L,kc,nw),NW=1,NWQV)
       END DO
C
      CLOSE(1)
C
      END IF
C
  110 FORMAT(1X,2I5,2X,1p,7E11.3,/,15X,7E11.3,/,15X,7E11.3)
  112 FORMAT('# Wet atmospheric deposition diagnostic file',/,
     $   ' N, TIME = ', I10, F12.5/)
C
C**********************************************************************C
C
      RETURN
      END
