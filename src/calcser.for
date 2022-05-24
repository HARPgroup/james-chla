C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALCSER (ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
C ** SUBROUTINE CALPSER UPDATES TIME VARIABLE SALINITY, TEMPERATURE
C ** DYE, SEDIMENT, AND SHELL FISH LARVAE
C ** BOUNDARY CONDITIONS AND INFLOW CONCENTRATIONS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
C **  INITIALIZE NULL SERIES CONCENTRATIONS
C
      NTOX1=NTOX
      IF(NTOX.EQ.0) NTOX1=NPCB
      NTT=4+NTOX1+NSED+NSND
      DO NT=1,NTT
       CQWRSERT(0,NT)=0.
       DO K=1,KC
        CSERT(K,0,NT)=0.
       END DO
      END DO
C
C**********************************************************************C
C
C **  CONCENTRATION SERIES INTERPOLTATION, SAL,TEM,DYE,SFL
C
      DO NC=1,4
      IF (ISTRAN(NC).EQ.0) GO TO 200
C
       DO NS=1,NCSER(NC)
C
       IF (ISTL.EQ.2) THEN
         TIME=DT*(FLOAT(N)-0.5)/TCCSER(NS,NC)
     &            +TBEGIN*(TCON/TCCSER(NS,NC))
        ELSE
         TIME=DT*FLOAT(N-1)/TCCSER(NS,NC)
     &           +TBEGIN*(TCON/TCCSER(NS,NC))
       END IF
       
! 
! %change timestep
!      
        IF(NCSTEP.GT.0) THEN   ! J.S. 12/24/2010
        IF (ISTL.EQ.2) THEN
         TIME=(SECDLAST-0.5*DT)/TCCSER(NS,NC)
     &     +TBEGIN*(TCON/TCCSER(NS,NC))
        ELSE
         TIME=(SECDLAST-DT)/TCCSER(NS,NC)+TBEGIN*(TCON/TCCSER(NS,NC))
        END IF        
       ENDIF
              
       M1=MCTLAST(NS,NC)
  100  CONTINUE
       M2=M1+1
c     write(612,619) N,ns,nc,dt,tccser(NS,NC),Tcon,tbegin,time
c     write(611,619) M2,NS,NC,time,tcser(m2,ns,nc)
c619	format(3i6,9f10.2)
       IF (TIME.GT.TCSER(M2,NS,NC)) THEN
         M1=M2
         GO TO 100
        ELSE
         MCTLAST(NS,NC)=M1
       END IF
C
       TDIFF=TCSER(M2,NS,NC)-TCSER(M1,NS,NC)
       WTM1=(TCSER(M2,NS,NC)-TIME)/TDIFF
       WTM2=(TIME-TCSER(M1,NS,NC))/TDIFF
        DO K=1,KC
        CSERT(K,NS,NC)=WTM1*CSER(M1,K,NS,NC)+WTM2*CSER(M2,K,NS,NC)
        END DO
C
C      WRITE(6,6000)N,CSERT(1,NS,NC),CSERT(KC,NS,NC)

C    
 
       IF(ISCDA(1).GT.0) THEN                 !% this part need rework
       
        if(NS.GE.NCSERA(1,1).and.NC.EQ.1) then
         if(abs(TCSER(M2,NS,NC)-TIME).LT.STLAGCDA)then   ! STLAGCDA in days   
         CSERTT(NS)=TCSER(M2,NS,NC)    
!        DO K=1,KC
!          CSERT(K,NS,NC)=CSER(M2,K,NS,NC)
!        ENDDO
         elseif(abs(TIME-TCSER(M1,NS,NC)).LT.STLAGCDA)then
c         CSERTT(NS)=TCSER(M1,NS,NC) 
!         DO K=1,KC
!          CSERT(K,NS,NC)=CSER(M1,K,NS,NC)
!         ENDDO           
         else
          DO K=1,KC
c          CSERT(K,NS,NC)=0.0   
          ENDDO 
         endif
        endif
c       
       ENDIF
      
      END DO      
      
  200 CONTINUE
      END DO
C
C **  CONCENTRATION SERIES INTERPOLTATION FOR  TOX
C
      IF (ISTRAN(5).GE.1) THEN
C
      DO NT=1,NTOX1
      NC=MSVTOX(NT)
       DO NS=1,NCSER(NC)
C
       IF (ISTL.EQ.2) THEN
         TIME=DT*(FLOAT(N)-0.5)/TCCSER(NS,NC)
     &            +TBEGIN*(TCON/TCCSER(NS,NC))
        ELSE
         TIME=DT*FLOAT(N-1)/TCCSER(NS,NC)
     &           +TBEGIN*(TCON/TCCSER(NS,NC))
       END IF
 ! 
! %change timestep
!      
        IF(NCSTEP.GT.0) THEN   ! J.S. 12/24/2010
        IF (ISTL.EQ.2) THEN
         TIME=(SECDLAST-0.5*DT)/TCCSER(NS,NC)
     &     +TBEGIN*(TCON/TCCSER(NS,NC))
        ELSE
         TIME=(SECDLAST-DT)/TCCSER(NS,NC)+TBEGIN*(TCON/TCCSER(NS,NC))
        END IF        
       ENDIF
       
       M1=MCTLAST(NS,NC)
  101  CONTINUE
       M2=M1+1
       IF (TIME.GT.TCSER(M2,NS,NC)) THEN
         M1=M2
         GO TO 101
        ELSE
         MCTLAST(NS,NC)=M1
       END IF
C
       TDIFF=TCSER(M2,NS,NC)-TCSER(M1,NS,NC)
       WTM1=(TCSER(M2,NS,NC)-TIME)/TDIFF
       WTM2=(TIME-TCSER(M1,NS,NC))/TDIFF
        DO K=1,KC
        CSERT(K,NS,NC)=WTM1*CSER(M1,K,NS,NC)+WTM2*CSER(M2,K,NS,NC)
        END DO
C
C      WRITE(6,6000)N,CSERT(1,NS,NC),CSERT(KC,NS,NC)
       END DO
C
      END DO
      END IF
C
C **  CONCENTRATION SERIES INTERPOLTATION FOR  SED
C
      IF (ISTRAN(6).GE.1) THEN
C
      DO NT=1,NSED
      NC=MSVSED(NT)
       DO NS=1,NCSER(NC)
C
       IF (ISTL.EQ.2) THEN
         TIME=DT*(FLOAT(N)-0.5)/TCCSER(NS,NC)
     &            +TBEGIN*(TCON/TCCSER(NS,NC))
        ELSE
         TIME=DT*FLOAT(N-1)/TCCSER(NS,NC)
     &           +TBEGIN*(TCON/TCCSER(NS,NC))
       END IF
! 
! %change timestep
!      
        IF(NCSTEP.GT.0) THEN   ! J.S. 12/24/2010
        IF (ISTL.EQ.2) THEN
         TIME=(SECDLAST-0.5*DT)/TCCSER(NS,NC)
     &     +TBEGIN*(TCON/TCCSER(NS,NC))
        ELSE
         TIME=(SECDLAST-DT)/TCCSER(NS,NC)+TBEGIN*(TCON/TCCSER(NS,NC))
        END IF        
       ENDIF
       
c      write(908,9009) ns,nc,n,TIME,DT,TCCSER(NS,NC),TBEGIN,TCON
       M1=MCTLAST(NS,NC)
  102  CONTINUE
       M2=M1+1
c	write(909,9009) m2,ns,nc, time,TCSER(M2,NS,NC)
9009	format(3i6,9e12.4)
       IF (TIME.GT.TCSER(M2,NS,NC)) THEN
         M1=M2
         GO TO 102
        ELSE
         MCTLAST(NS,NC)=M1
       END IF
C
       TDIFF=TCSER(M2,NS,NC)-TCSER(M1,NS,NC)
       WTM1=(TCSER(M2,NS,NC)-TIME)/TDIFF
       WTM2=(TIME-TCSER(M1,NS,NC))/TDIFF
        DO K=1,KC
        CSERT(K,NS,NC)=WTM1*CSER(M1,K,NS,NC)+WTM2*CSER(M2,K,NS,NC)
        END DO
C
C      WRITE(6,6000)N,CSERT(1,NS,NC),CSERT(KC,NS,NC)
       END DO
C
      END DO
      END IF
C
C **  CONCENTRATION SERIES INTERPOLTATION FOR  SND
C
      IF (ISTRAN(7).GE.1) THEN
C
      DO NT=1,NSND
      NC=MSVSND(NT)
       DO NS=1,NCSER(NC)
C
       IF (ISTL.EQ.2) THEN
         TIME=DT*(FLOAT(N)-0.5)/TCCSER(NS,NC)
     &            +TBEGIN*(TCON/TCCSER(NS,NC))
        ELSE
         TIME=DT*FLOAT(N-1)/TCCSER(NS,NC)
     &           +TBEGIN*(TCON/TCCSER(NS,NC))
       END IF
! 
! %change timestep
!      
        IF(NCSTEP.GT.0) THEN   ! J.S. 12/24/2010
        IF (ISTL.EQ.2) THEN
         TIME=(SECDLAST-0.5*DT)/TCCSER(NS,NC)
     &     +TBEGIN*(TCON/TCCSER(NS,NC))
        ELSE
         TIME=(SECDLAST-DT)/TCCSER(NS,NC)+TBEGIN*(TCON/TCCSER(NS,NC))
        END IF        
       ENDIF
       
       M1=MCTLAST(NS,NC)
  103  CONTINUE
       M2=M1+1
       IF (TIME.GT.TCSER(M2,NS,NC)) THEN
         M1=M2
         GO TO 103
        ELSE
         MCTLAST(NS,NC)=M1
       END IF
C
       TDIFF=TCSER(M2,NS,NC)-TCSER(M1,NS,NC)
       WTM1=(TCSER(M2,NS,NC)-TIME)/TDIFF
       WTM2=(TIME-TCSER(M1,NS,NC))/TDIFF
        DO K=1,KC
        CSERT(K,NS,NC)=WTM1*CSER(M1,K,NS,NC)+WTM2*CSER(M2,K,NS,NC)
        END DO
C
C      WRITE(6,6000)N,CSERT(1,NS,NC),CSERT(KC,NS,NC)
       END DO
C
      END DO
      END IF
C
C **  WRITE DIAGNOSTIC FILE FOR CSER INTERPOLTATION
C
      IF(ISDIQ.GE.1.AND.N.EQ.1) THEN
      OPEN(1,FILE='cdiag.out',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='cdiag.out',STATUS='UNKNOWN')
      DO NC=1,NTT
       WRITE(1,1001)NC
       DO NS=1,NCSER(NC)
         WRITE(1,1002)NS,(CSERT(K,NS,NC),K=1,KC)
       END DO
      END DO
      CLOSE(1)
      END IF
 1001 FORMAT(/' TRANSPORT VARIABLE ID =',I5/)
 1002 FORMAT(I5,2X,12E12.4)
C
C **  SHELL FISH LARVAE BEHAVIOR TIME SERIES INTERPOLTATION 
C
      IF (ISTRAN(4).EQ.0) GO TO 400
C
       IF (ISTL.EQ.2) THEN
         TIME=DT*(FLOAT(N)-0.5)/TCSFSER
     &            +TBEGIN*(TCON/TCSFSER)
        ELSE
         TIME=DT*FLOAT(N-1)/TCSFSER
     &           +TBEGIN*(TCON/TCSFSER)
       END IF
! 
! %change timestep
!      
        IF(NCSTEP.GT.0) THEN   ! J.S. 12/24/2010
        IF (ISTL.EQ.2) THEN
         TIME=(SECDLAST-0.5*DT)/TCSFSER
     &     +TBEGIN*(TCON/TCSFSER)
        ELSE
         TIME=(SECDLAST-DT)/TCSFSER+TBEGIN*(TCON/TCSFSER)
        END IF        
       ENDIF
       
       M1=MSFTLST
  300  CONTINUE
       M2=M1+1
       IF (TIME.GT.TSFSER(M2)) THEN
         M1=M2
         GO TO 300
        ELSE
         MSFTLST=M1
       END IF      
C
       TDIFF=TSFSER(M2)-TSFSER(M1)
       WTM1=(TSFSER(M2)-TIME)/TDIFF
       WTM2=(TIME-TSFSER(M1))/TDIFF
       RKDSFLT=WTM1*RKDSFL(M1)+WTM2*RKDSFL(M2)
       WSFLSTT=WTM1*WSFLST(M1)+WTM2*WSFLST(M2)
       WSFLSMT=WTM1*WSFLSM(M1)+WTM2*WSFLSM(M2)
       DSFLMNT=WTM1*DSFLMN(M1)+WTM2*DSFLMN(M2) 
       DSFLMXT=WTM1*DSFLMX(M1)+WTM2*DSFLMX(M2)
       SFNTBET=WTM1*SFNTBE(M1)+WTM2*SFNTBE(M2)
       SFATBTT=WTM1*SFATBT(M1)+WTM2*SFATBT(M2)
C
  400 CONTINUE
C
 6000 FORMAT('N, CSERT(1),CSERT(KC) = ',I6,4X,2F12.2)
C
C**********************************************************************C
C
      RETURN
      END
