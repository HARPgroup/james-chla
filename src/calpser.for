C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALPSER (ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
C ** SUBROUTINE CALPSER UPDATES TIME VARIABLE SURFACE ELEVATION 
C ** BOUNDARY CONDITIONS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      PSERT(0)=0.
C
      DO NS=1,NPSER
C
      TIME=DT*FLOAT(N)/TCPSER(NS)+TBEGIN*(TCON/TCPSER(NS))
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCPSER(NS)+TBEGIN*(TCON/TCPSER(NS))  !% J.S. 1/31/2014 
!      write(*,*)'P time = ',time
      M1=MPTLAST(NS)
  100 CONTINUE
      M2=M1+1
      IF (TIME.GT.TPSER(M2,NS)) THEN
       M1=M2
       GO TO 100
      ELSE
       MPTLAST(NS)=M1
      END IF      
C
      TDIFF=TPSER(M2,NS)-TPSER(M1,NS)
      WTM1=(TPSER(M2,NS)-TIME)/TDIFF
      WTM2=(TIME-TPSER(M1,NS))/TDIFF
      PSERT(NS)=WTM1*PSER(M1,NS)+WTM2*PSER(M2,NS)
C     WRITE(6,6000)N,PSERT(NS)
C
      END DO
C
 6000 FORMAT('N, PSERT = ',I6,4X,F12.4)
C
C**********************************************************************C
C
      RETURN
      END
