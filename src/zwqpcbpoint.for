        SUBROUTINE pcbpoint
C
C **  LAST MODIFIED BY J.Shen 2016
C
C**********************************************************************C
C
C **  SUBROUTINE Read point source for pcb as mass
C     Using water quality arrary
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'wq.par'      
      INCLUDE 'efdc.cmn'
      INCLUDE 'wqcom.cmn' 
      INTEGER, SAVE :: IN_PCB = 0 
      DIMENSION RLDTMP(21) !,IACP(NACP),IACP1(NACP1)
      
      IF(IN_PCB.EQ.0) THEN   ! Read input files
      
      IN_PCB=1
!  
! Initialize combined PSL array to zero:
!
      DO NW=1,21
       DO K=1,KC
        DO L=2,LA
         WQWPSL(L,K,NW) = 0.0
        END DO
       END DO
      END DO  
! 
       OPEN(10,FILE='pcb1.dia',STATUS='UNKNOWN')
       CLOSE(10,STATUS='DELETE')
      
      OPEN(3,FILE='pcbser.inp',STATUS='UNKNOWN') 
      read(3,*) NPCBpt
      if(NPCBpt.LE.0) then
      write(*,*)'No point source inputs !'
      close(3)
      return
      else
 ! Read in point source location of PCB. It can be differ from hydrodynamis model location
 ! It is more flaxible. Total point source set 200 (change (eddc.cmn))
      DO i=1,NPCBpt  
      read(3,*) tmp1,IPCBpt(i),JPCBpt(i)
      ENDDO
      endif
         
       DO NW=1,NPCB
        RLDTMP(NW)=0.0 
       END DO
C
C **  SKIP OVER TITLE AND HEADER LINES
C
       DO IS=1,13
         READ(3,1)
       END DO
   1   FORMAT(120X)
C
       DO NS=1,NPCBpt
        MWQPTLT(NS)=1
        READ(3,*,IOSTAT=ISO)MWQPSR(NS),TCWQPSR(NS),
     $                  TAWQPSR(NS),RMULADJ,ADDADJ
        IF(ISO.GT.0) GO TO 900

        READ(3,*,IOSTAT=ISO) (WKWQ(NS,K),K=1,KC) ! weiting for each layer
        IF(ISO.GT.0) GO TO 900
c 
        if(MWQPSR(NS).gt.NDWQPSR) then   ! check array dimensions
         write(*,*) " Error:  MWQPSR(NS)>NDWQPSR ", MWQPSR(NS),NDWQPSR
         stop
        endif 
        
        RMULADJ=1.0*RMULADJ   !  g per day
        ADDADJ=ADDADJ
        
         DO M=1,MWQPSR(NS)
          READ(3,*,IOSTAT=ISO)TWQPSER(M,NS),(RLDTMP(NW),NW=1,NPCB)
          IF(ISO.GT.0) then
          write(*,*)'PCB point source input error ',M, TWQPSER
          GO TO 900
          ENDIF
          !      Temporally use C,N, P as reduction test
          RLDTMP(1)=RLDTMP(1)*(1-RC_R)
          RLDTMP(2)=RLDTMP(2)*(1-RN_R)
          RLDTMP(3)=RLDTMP(3)*(1-RP_R)
          TWQPSER(M,NS)=TWQPSER(M,NS)+TAWQPSR(NS)
          WQPSSER(M,1,NS)=RMULADJ*RLDTMP(1) 
          WQPSSER(M,2,NS)=RMULADJ*RLDTMP(2) 
          WQPSSER(M,3,NS)=RMULADJ*RLDTMP(3) 
          DO kk=4,21
            WQPSSER(M,kk,NS)= 0
          ENDDO        
         END DO    
      ENDDO
 900  continue          
      close(3)
      
      ELSE       ! get input 
 
C
C**********************************************************************C
C
C **  LOADING SERIES INTERPOLTATION
C
    !  write(*,*)'Total nonpint source ',NPSTMSR
       DO NS=1,NPCBpt

        TIME=DT*FLOAT(N-1)/TCWQPSR(NS)
     &           +TBEGIN*(TCON/TCWQPSR(NS))
     
       IF(NCSTEP.GT.0) TIME=(SECDLAST-DT)/TCON+TBEGIN  !% J.S. 6/16/2014  
C
       M1=MWQPTLT(NS)
  100  CONTINUE
       M2=M1+1
       if(M2.GT.MWQPSR(NS))then
       write(*,*)'Point source problem !', M,NS,MWQPSR(NS)
       pause
       endif
       IF (TIME.GT.TWQPSER(M2,NS)) THEN
         M1=M2
         GO TO 100
        ELSE
         MWQPTLT(NS)=M1
       ENDIF
C
       TDIFF=TWQPSER(M2,NS)-TWQPSER(M1,NS)
       WTM1=(TWQPSER(M2,NS)-TIME)/TDIFF
       WTM2=(TIME-TWQPSER(M1,NS))/TDIFF
       DO NW=1,NPCB
         WQPSSRT(NW,NS)=WTM1*WQPSSER(M1,NW,NS)
     &                   +WTM2*WQPSSER(M2,NW,NS)
       ENDDO
C
      ENDDO
C
C**********************************************************************C
!
! Initialize combined PSL array to zero:
!
       DO NW=1,NPCB
        DO K=1,KC
         DO L=2,LA
         WQWPSL(L,K,NW) = 0.0
         END DO
        END DO
       END DO
!
!  WQWPSLC:  is for constant point source input. Not supported (AEE, 3/20/2007)
!  WQWPSL:   final point loads go to model
!  WQPSSRT : input point time series
!

!       WRITE(1,112)N,TIME

      IF(N.LE.10) THEN
       OPEN(10,FILE='pcb.dia',STATUS='UNKNOWN')
       CLOSE(10,STATUS='DELETE')
      ENDIF 

       IF(N.LE.10) WRITE(900,112)N,TIME
!       OPEN(10,FILE='pcb.dia',ACCESS='APPEND')
!       WRITE(900,112)N,TIME

      DO NS=1,NPCBpt
        L = LIJW(IPCBpt(NS),JPCBpt(NS))
        ITMP = MVPSL(NS)
        DO k=1,kc
         DO NW=1,NPCB
	    WQWPSL(L,K,NW) = WQWPSL(L,K,NW)
     +         +  WQPSSRT(NW,NS)*WkWQ(NS,K)
         END DO
         IF(N.LE.10) THEN
          WRITE(900,110)ILW(L),JLW(L),K,(WQWPSL(L,K,NW),NW=1,NPCB)
	   ENDIF   
        ENDDO
      END DO
 !     ENDIF
705   format(i4, f12.4,7i4,99f12.4)


!      IF(N.LE.10) CLOSE(10)

!      END IF

  110 FORMAT(1X,3I4,2X,7E12.4,/,15X,7E12.4,/,15X,7E12.4)
  111 FORMAT(1X,I4,2X,7E12.4,/,7X,7E12.4,/,7X,7E12.4)
  112 FORMAT(' N, TIME = ', I10, F12.5/)
       
       ENDIF
       
       return
       end