C**********************************************************************C
C 
      SUBROUTINE RWQSUN(DT,N,TBEGIN,TCON,NTSPTC,NCSTEP,SECDLAST)
C
C**********************************************************************C
C
C **  NEW VERSION BY J. M. HAMRICK  7 APRIL 1997
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 10 april 1999
C
C**********************************************************************C
C
C **  READS AND INTERPOLATES DAILY AVERAGE SOLAR RADIATION AND
C **  DAYLIGHT FRACTION
C
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
      INTEGER, SAVE :: NSUNDAY
      REAL*8 SECDLAST 
C
C**********************************************************************C
C
      IF(ITNWQ.GT.0) GO TO 1000
C
C**********************************************************************C
C
C **  READ IN DAILY AVERAGE SOLAR SW RAD SERIES FROM FILE 'sunday.inp'
C
C----------------------------------------------------------------------C
C
      OPEN(1,FILE='wqsunser.inp',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES WIHT 'C'
C
      CALL SKIPCOMM(1,'C')
C
      M=0
      ISPAR=1
      READ(1,*,IOSTAT=ISO)NSUNDAY,TCSUNDAY,
     $                   TASUNDAY,RMULADJ,ADDADJ
      IF(ISO.GT.0) GO TO 900
      DO M=1,NSUNDAY
        READ(1,*,IOSTAT=ISO)TSSRD(M),SOLSRD(M) 
        IF(ISO.GT.0) GO TO 900
        TSSRD(M)=TCSUNDAY*( TSSRD(M)+TASUNDAY )
        SOLSRD(M)=RMULADJ*(SOLSRD(M)+ADDADJ)
      END DO
C
      CLOSE(1)
	open(221,file='sor.dia')
C
      GO TO 901
C
  900 CONTINUE
      WRITE(6,601)M
      STOP
C
  901 CONTINUE
C
    1 FORMAT(120X)
  601 FORMAT(' READ ERROR FILE sunser.inp ')
C
C**********************************************************************C
C
 1000 CONTINUE
C
C**********************************************************************C
C
C **  DAILY AVERAGE SOLAR SW RADIATION INTERPOLTATION FOR WATER QUALITY
C
       TIME=DT*FLOAT(N)
     &           +TBEGIN*(TCON)
     
      IF(NCSTEP.GT.0) TIME=SECDLAST+TBEGIN*TCON  ! Note the time is in second 6/16/2014  
!
       M1=ISPAR
!
  100  CONTINUE
       M2=M1+1
       IF (TIME.GT.TSSRD(M2).and.M2.LT.NSUNDAY) THEN
         M1=M2
         GO TO 100
        ELSE
         ISPAR=M1
       END IF
!
       TDIFF=TSSRD(M2)-TSSRD(M1)
       WTM1=(TSSRD(M2)-TIME)/TDIFF
       WTM2=(TIME-TSSRD(M1))/TDIFF
       SUNDT=WTM1*SOLSRD(M1)+WTM2*SOLSRD(M2)
!       write(*,*)'Time Sun ',TIME,SUNDT,TIME/3600./24.
!
!**********************************************************************C
!

      RETURN
      END
