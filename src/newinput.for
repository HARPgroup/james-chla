c       
	 SUBROUTINE newinput
	   !(iyear,tplus)
c----------------------------------------------------
c  Bo Hong, 2009-Jun.
c  for running the long-term case
c  the input file should be prepared year by year.
c----------------------------------------------------
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
	character*2 mypser, myqser, mywin, mysser  
C	character*4 iyear
C	real tplus

       print*,'iyear= ',iyear

	   mypser='ye'
	   myqser='no'
	   mywin='ye'
         mysser='ye'
         
         if(ifed_inc.LT.0) then
 	   mypser='no'
	   myqser='no'
         mysser='no'        
         endif

    1 FORMAT (120X)

       JSTMSR=1

      DO NS=1,NWSER
        MWTLAST(NS)=1
      END DO
C
      DO NS=1,NPSER
       MPTLAST(NS)=1
      END DO
C
      DO NS=1,NQSER
       MQTLAST(NS)=1
      END DO
C
      NTMP=4+NSED+NSND+NTOX
      DO NC=1,NTMP
      DO NN=1,NCSER(NC)
      MCTLAST(NN,NC)=1
      END DO
      END DO
c-------------- pser ------------------             
c 
                if(mypser.eq.'ye') then

      IF (NPSER.GE.1) THEN	  
        OPEN(1,FILE='pser_'//iyear//'.inp',STATUS='UNKNOWN')
c       OPEN(1,FILE='pser.inp',STATUS='UNKNOWN')

C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,6
      READ(1,1)
      END DO
C
        DO NS=1,NPSER
        READ (1,*,IOSTAT=ISO) MPSER(NS),TCPSER(NS),TAPSER(NS),
     $                   RMULADJ,ADDADJ
        IF(ISO.GT.0) GO TO 850
         DO M=1,MPSER(NS)
         READ(1,*,IOSTAT=ISO)TPSER(M,NS),PSERTMP
        IF(ISO.GT.0) GO TO 850
c         TPSER(M,NS)=TPSER(M,NS)+TAPSER(NS)
         TPSER(M,NS)=TPSER(M,NS)+TAPSER(NS)+tplus*(TCON/TCPSER(NS))
         PSER(M,NS)=G*(PSERTMP+ADDADJ)*RMULADJ
         END DO
        END DO
        CLOSE(1)
      END IF

chong----- output for test ------------------
!         open(200, file='testpser'//iyear//'.dat',status='unknown')
!        DO NS=1,NPSER
!          write (200,*) MPSER(NS),TCPSER(NS),TAPSER(NS),
!     $                   RMULADJ,ADDADJ
!          DO M=1,MPSER(NS)
!            write(200,*)TPSER(M,NS),PSER(M,NS)
!          END DO
!        END DO
!	   close(200)
c-------------------------------------------------	   	         
	           endif           ! end of mypser
C	         pause 'pser.inp'
chong------------qser--------------------------
C
               if(myqser.eq.'ye') then

      IF(NQSER.GE.1) THEN
 !     OPEN(1,FILE='qser_'//iyear//'.inp',STATUS='UNKNOWN')
       OPEN(1,FILE='qser.inp',STATUS='UNKNOWN')

C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
        DO IS=1,15
        READ(1,1)
        END DO
C
        DO NS=1,NQSER
C        NTMP=0
        READ(1,*,IOSTAT=ISO)ISTYP, MQSER(NS),TCQSER(NS),TAQSER(NS),
     $                   RMULADJ,ADDADJ,ICHGQS
        IF(ISO.GT.0) GO TO 860
        IF(ISTYP.EQ.1) THEN
!        Read(1,*,IOSTAT=ISO)
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GO TO 860
           DO K=1,KC
            WKQ(K)=0
           ENDDO
           WKQ(KC)=1.0

           DO M=1,MQSER(NS)
C           NTMP=NTMP+1
           READ(1,*,IOSTAT=ISO)TQSER(M,NS),QSERTMP
!           write(*,*)NS,TQSER(M,NS)
C	   WRITE(2,2222)NS,NTMP,TQSER(M,NS),QSERTMP
           IF(ISO.GT.0) GO TO 860
           TQSER(M,NS)=TQSER(M,NS)+TAQSER(NS)+tplus  !hong
           QSERTMP=(RMULADJ*(QSERTMP+ADDADJ))
           IF (ICHGQS.EQ.1) QSERTMP=MAX(QSERTMP,0.0)
           IF (ICHGQS.EQ.-1) QSERTMP=MIN(QSERTMP,0.0)
            DO K=1,KC
            QSER(M,K,NS)=QSERTMP*WKQ(K)
            END DO
           END DO
         ELSE
          DO M=1,MQSER(NS)
          READ(1,*,IOSTAT=ISO)TQSER(M,NS),(QSER(M,K,NS), K=1,KC)
          IF(ISO.GT.0) GO TO 860
          TQSER(M,NS)=TQSER(M,NS)+TAQSER(NS)+tplus  !hong
           DO K=1,KC
           QSER(M,K,NS)=RMULADJ*(QSER(M,K,NS)+ADDADJ)
           IF (ICHGQS.EQ.1) QSER(M,K,NS)=MAX(QSER(M,K,NS),0.0)
           IF (ICHGQS.EQ.-1) QSER(M,K,NS)=MIN(QSER(M,K,NS),0.0)
           END DO
          END DO
        END IF
        END DO
C
      CLOSE(1)
C      CLOSE(2)
      END IF
c--output---
!       open(200, file='testqser'//iyear//'.dat',status='unknown')
!
!       DO NS=1,NQSER
!           WRITE(200,*)ISTYP, MQSER(NS),TCQSER(NS),TAQSER(NS),
!     $                   RMULADJ,ADDADJ,ICHGQS
!           DO M=1,MQSER(NS)
!              WRITE(200,*)TQSER(M,NS),(QSER(M,K,NS), K=1,KC)
!           END DO
!        END DO
!
!	   close(200)
	                endif    ! myqser
C	PAUSE 'QSER HONG'
C
chong---------------------wser----------------              
             if (mywin.eq.'ye')  then

      IF(NWSER.GT.0) THEN
	IF(NWSER.GT.100) then
      OPEN(1,FILE='wser1_'//iyear//'.inp',STATUS='UNKNOWN')
	ELSE
      OPEN(1,FILE='wser_'//iyear//'.inp',STATUS='UNKNOWN')
	ENDIF
c      OPEN(1,FILE='wser.inp',STATUS='UNKNOWN')

C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,17
      READ(1,1)
      END DO
C
      DO NS=1,NWSER
C
      READ (1,*,IOSTAT=ISO) MWSER(NS),TCWSER(NS),TAWSER(NS),WINDSCT,
     $                      ISWDINT(NS),WCUT
      IF(ISO.GT.0) GO TO 940
c      write(*,*)n,MWSER(NS),TCWSER(NS),TAWSER(NS),WINDSCT
C
      DO M=1,MWSER(NS)
      READ(1,*,IOSTAT=ISO)TWSER(M,NS),WINDS(M,NS),WINDD(M,NS)
      if(WINDS(M,NS).GT.WCUT) WINDS(M,NS)=WCUT
      IF(ISO.GT.0) GO TO 940
      END DO
C
      DO M=1,MWSER(NS)
       TWSER(M,NS)=TWSER(M,NS)+TAWSER(NS)+tplus !hong
      END DO
C
      IF(ISWDINT(NS).LE.1) THEN
      DO M=1,MWSER(NS)
       WINDS(M,NS)=WINDSCT*WINDS(M,NS)
      END DO
      END IF
C
      IF(ISWDINT(NS).EQ.1) THEN
      DO M=1,MWSER(NS)
        IF(WINDD(M,NS).LE.180.) THEN
          WINDD(M,NS)=WINDD(M,NS)+180.
          IF(WINDD(M,NS).EQ.360.) WINDD(M,NS)=0.
         ELSE
          WINDD(M,NS)=WINDD(M,NS)-180.
          IF(WINDD(M,NS).EQ.360.) WINDD(M,NS)=0.
        END IF
      END DO
      END IF
C
      IF(ISWDINT(NS).EQ.2) THEN
      DO M=1,MWSER(N)
       WINDS(M,NS)=WINDSCT*WINDS(M,NS)
       WINDD(M,NS)=WINDSCT*WINDD(M,NS)
      END DO
      END IF
C
      END DO
C
      CLOSE(1)
      END IF
c--- output for check ---------
!          open(200, FILE='testwser'//iyear//'.dat',status='unknown') 
!c	1   form='unformatted')
!      DO NS=1,NWSER
!         WRITE (200,*) MWSER(NS),TCWSER(NS),TAWSER(NS),WINDSCT,
!     $                      ISWDINT(NS),WCUT
!         DO M=1,MWSER(NS)
!            WRITE(200,*)TWSER(M,NS),WINDS(M,NS),WINDD(M,NS)
!         END DO
!      END DO
!          close(200)
c-------------------------wndmap-------------
C
      IF(NWSER.GT.1) THEN
	 IF(NWSER.LT.100) then
        OPEN(1,FILE='wndmap_'//iyear//'.inp',STATUS='UNKNOWN')
c        OPEN(1,FILE='wndmap.inp',STATUS='UNKNOWN')
        DO IS=1,4
          READ(1,1)
        END DO
        DO L=2,LA
          READ(1,*)LD,ID,JD,(WNDWHT(LD,NS),NS=1,NWSER)
        END DO
        CLOSE(1)
	 ELSE
        OPEN(1,FILE='wndmap1_'//iyear//'.inp',STATUS='UNKNOWN')
c        OPEN(1,FILE='wndmap.inp',STATUS='UNKNOWN')
        DO IS=1,4
          READ(1,1)
        END DO
        DO L=2,LA
          READ(1,*)LD,ID,JD,IDT
          DO NS=1,NWSER
            READ(1,*)IDWNDW(LD,NS),WNDWHT(LD,NS)     ! Read weight for 10 nearby station
          END DO
        END DO
        CLOSE(1)
	 ENDIF
      END IF

	      endif ! mywin

c	pause 'wser '
chong-----------------------------------------
               
C---------- sser.inp -------------------

      if (mysser.eq.'ye')  then
      IF(NCSER(1).GE.1) THEN
        OPEN(1,FILE='sser_'//iyear//'.inp',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,15
      READ(1,1)
      END DO
C
        NC=1
        DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),
     $                   TACSER(NS,NC),RMULADJ,ADDADJ
        IF(ISO.GT.0) GO TO 870
        IF(ISTYP.EQ.1) THEN   !hong (ISTYP=0 in sser.inp)
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GO TO 870
           DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),CSERTMP
           IF(ISO.GT.0) GO TO 870
           TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)+tplus   !Rico
            DO K=1,KC
            CSER(M,K,NS,NC)=(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
            END DO
           END DO
         ELSE
          DO M=1,MCSER(NS,NC)
          READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),(CSER(M,K,NS,NC), K=1,KC)
 ! 1999-2000
 !         WKQ(1)=1.2
 !         WKQ(2)=1.175
 !         WKQ(3)=1.15  
 !         WKQ(4)=1.125 
 !         WKQ(5)=1.1   
 !         WKQ(6)=1.09  
 !         WKQ(7)=1.08  
 !         WKQ(8)=1.07  
 ! 2006-2012
          WKQ(1)=1.
          WKQ(2)=1.
          WKQ(3)=1. 
          WKQ(4)=1. 
          WKQ(5)=1. 
          WKQ(6)=1. 
          WKQ(7)=1. 
          WKQ(8)=1.
           
          IF(ISO.GT.0) GO TO 870
          TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)+tplus   !Rico
           DO K=1,KC
           CSER(M,K,NS,NC)=RMULADJ*(CSER(M,K,NS,NC)+ADDADJ)*WKQ(K)
           END DO
          END DO
        END IF
        END DO
        CLOSE(1)
      END IF


C
      IF(NCSER(2).GE.1) THEN
!        OPEN(1,FILE='tser.inp',STATUS='UNKNOWN')
        OPEN(1,FILE='tser_'//iyear//'.inp',STATUS='UNKNOWN')
C
C **  SKIP OVER TITLE AND AND HEADER LINES
C
      DO IS=1,15
      READ(1,1)
      END DO
C
        NC=2
        DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),
     $                   TACSER(NS,NC),RMULADJ,ADDADJ
        IF(ISO.GT.0) GO TO 880
        IF(ISTYP.EQ.1) THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF(ISO.GT.0) GO TO 880
           DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),CSERTMP
           IF(ISO.GT.0) GO TO 880
           TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)
            DO K=1,KC
            CSER(M,K,NS,NC)=(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
            END DO
           END DO
         ELSE
          DO M=1,MCSER(NS,NC)
          READ(1,*,IOSTAT=ISO)TCSER(M,NS,NC),(CSER(M,K,NS,NC), K=1,KC)
          IF(ISO.GT.0) GO TO 880
          TCSER(M,NS,NC)=TCSER(M,NS,NC)+TACSER(NS,NC)+tplus  
           DO K=1,KC
           CSER(M,K,NS,NC)=RMULADJ*(CSER(M,K,NS,NC)+ADDADJ)
           END DO
          END DO
        END IF
        END DO
        CLOSE(1)
      END IF
C

      endif  ! endif sser


      GO TO 3000

  850 WRITE(6,851)
      WRITE(8,851)
  851 FORMAT(1X,'READ ERROR FOR FILE pser.inp ')
      STOP

  860 WRITE(6,861)
      WRITE(8,861)
  861 FORMAT(1X,'READ ERROR FOR FILE qser.inp ')
      STOP

  870 WRITE(6,871)
  871 FORMAT(1X,'READ ERROR FOR FILE sser.inp ')
      STOP

  880 WRITE(6,881)
  881 FORMAT(1X,'READ ERROR FOR FILE tser.inp ')
      STOP

  940 WRITE(6,941)
      WRITE(8,941)
  941 FORMAT(1X,'READ ERROR FOR FILE aser.inp ')
      STOP

 3000 CONTINUE

      return
	END
