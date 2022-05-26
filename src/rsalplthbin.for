C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RSALPLTH(ICON,CONC)
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
      DIMENSION CONC(LCM,KCM)
      REAL DBS(10)
      CHARACTER*80 TITLE
C
C**********************************************************************C
C
      IF (JSRSPH(ICON).NE.1) GO TO 300
C
C----------------------------------------------------------------------C
C
      LINES=LA-1
      LEVELS=2
      LEVELSS=3
      DBS(1)=0.
      DBS(2)=99.
      DBS(3)=-99.
C
      IF (ISTRAN(1).GE.1)CLOSE(531)
      IF (ISTRAN(2).GE.1)CLOSE(532)
      IF (ISTRAN(3).GE.1)CLOSE(533) 
      IF (ISTRAN(5).GE.1)CLOSE(535)
      IF (ISTRAN(4).GE.1)CLOSE(534)        
      IF (ISTRAN(6).GE.1)CLOSE(536)   
              
      IF (ISTRAN(1).GE.1) THEN
        TITLE='RESIDUAL HORIZONTAL SALINITY CONTOURS'
        LUN=531
        OPEN(LUN,FILE='rsalcnh'//iyear//'.bin',STATUS='UNKNOWN',
     &   form='unformatted')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='rsalcnh'//iyear//'.bin',form='unformatted')
        WRITE (LUN) LC-1,KC+2
      END IF
C
      IF (ISTRAN(2).GE.1) THEN
        TITLE='RESIDUAL HORIZONTAL TEMPERATURE CONTOURS'
        LUN=532
        OPEN(LUN,FILE='rtemcnh'//iyear//'.bin',STATUS='UNKNOWN',
     &   form='unformatted')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='rtemcnh'//iyear//'.bin',form='unformatted')
        WRITE (LUN) LC-1,KC+2
      END IF
C
      IF (ISTRAN(3).GE.1) THEN
        TITLE='RESIDUAL HORIZONTAL DYE CONC CONTOURS'
        LUN=533
        OPEN(LUN,FILE='rdyecnh'//iyear//'.bin',STATUS='UNKNOWN',
     &  form='unformatted')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='rdyecnh'//iyear//'.bin',form='unformatted')
        WRITE (LUN) LC-1,KC+2
      END IF
C
      IF (ISTRAN(6).GE.1) THEN
        TITLE='RESIDUAL HORIZ COHESIVE SEDIMENT CONC CONTOURS'
        LUN=536
        OPEN(LUN,FILE='rsedcnh'//iyear//'.bin',STATUS='UNKNOWN',
     &  form='unformatted')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='rsedcnh'//iyear//'.bin',form='unformatted')
        WRITE (LUN) LC-1,KC+2
      END IF
C
!      IF (ISTRAN(7).GE.1) THEN
!        TITLE='RESIDUAL HORIZ NONCOH SEDIMENT CONC CONTOURS'
!        LUN=15
!        OPEN(LUN,FILE='rsndcnh.out',STATUS='UNKNOWN')
!        CLOSE(LUN,STATUS='DELETE')
!        OPEN(LUN,FILE='rsndcnh.out',STATUS='UNKNOWN')
!        WRITE (LUN,99) TITLE
!        WRITE (LUN,101)LINES,LEVELSS
!        WRITE (LUN,250)(DBS(L),L=1,LEVELSS)
!        CLOSE(LUN)
!      END IF
C
      IF (ISTRAN(5).GE.1) THEN
        TITLE='RESIDUAL HORIZ TOXIC CONTAM CONC CONTOURS'
        LUN=535
        OPEN(LUN,FILE='rtoxcnh'//iyear//'.bin',STATUS='UNKNOWN',
     & form='unformatted')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='rtoxcnh'//iyear//'.bin',form='unformatted')
        WRITE (LUN) LC-1,KC+2,NTOX

 !       TITLE='RESIDUAL HORIZ TOXIC PART FRAC CONTOURS'
 !       LUNF=26
 !       OPEN(LUNF,FILE='rtxpcnh.out',STATUS='UNKNOWN')
 !       CLOSE(LUNF,STATUS='DELETE')
 !       OPEN(LUNF,FILE='rtxpcnh.out',STATUS='UNKNOWN')
 !       WRITE (LUNF,99) TITLE
 !       WRITE (LUNF,101)LINES,LEVELSS
 !       WRITE (LUNF,250)(DBS(L),L=1,LEVELSS)
 !       CLOSE(LUNF)
      END IF
C
      IF (ISTRAN(4).GE.1) THEN
        TITLE='RESIDUAL HORIZONTAL SFL CONC CONTOURS'
        LUN=534
        OPEN(LUN,FILE='rsflcnh'//iyear//'.bin',STATUS='UNKNOWN',
     &  form='unformatted')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='rsflcnh'//iyear//'.bin',form='unformatted')
        WRITE (LUN) LC-1,KC+2
      END IF
C
      DO ITMP=1,7
       JSRSPH(ITMP)=0
      END DO 
C
C----------------------------------------------------------------------C
C
  300 CONTINUE
C

      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON    
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
C
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014
C
      IF(ICON.EQ.1) LUN=531  ! salt
      IF(ICON.EQ.2) LUN=532  ! temp
      IF(ICON.EQ.3) LUN=533  ! dye
      IF(ICON.EQ.4) LUN=534  ! Shellfish      
      IF(ICON.EQ.6) LUN=536  ! sed
      IF(ICON.EQ.7) LUN=537  
      IF(ICON.EQ.5) LUN=535  ! tox     
      WRITE (LUN)TIME
                   
      IF(ICON.LE.3) THEN
       DO L=2,LA      
         write(LUN)IL(L),JL(L),(CONC(L,K), k=1,kc)
       ENDDO
      END IF
     
      IF(ICON.EQ.6) THEN
       DO L=2,LA
C      SEDBT=1000*SSG*SEDB(L)
chong       WRITE(LUN,200)IL(L),JL(L),CONC(L,KC),CONC(L,1),
chong     $               SEDBTLPF(L)
       WRITE(LUN)IL(L),JL(L),(CONC(L,K),k=1,kc)
       END DO
      END IF
      
      IF(ICON.EQ.7) THEN
       DO L=2,LA
C      SEDBT=1000*SSG*SNDB(L,KBT(L))
       WRITE(LUN)IL(L),JL(L),CONC(L,KC),CONC(L,1),
     $               SNDBTLPF(L,1)
       END DO
      END IF
      
 !     IF(ICON.EQ.5) THEN
 !      DO J=1,NTOX
 !      DO L=2,LA
 !      TOXBT=1000.*TOXBLPF(L,1)
 !      WRITE(LUN,200)IL(L),JL(L),CONC(L,KC),CONC(L,1),
 !    $               TOXBT
 !      END DO

 !      DO L=2,LA
C      SEDBT=1000*SSG*SEDB(L,KBT(L),1)
!       WRITE(LUNF,200)IL(L),JL(L),TXPFLPF(L,KC,1,1),
!     $               TXPFLPF(L,1,1,1),TOXBLPF(L,1)
!       END DO
!       CLOSE(LUNF)
!      END IF
      
      IF(ICON.EQ.4) THEN
       DO L=2,LA
       WRITE(LUN)IL(L),JL(L),(CONC(L,K),k=1,kc),
     $                 SFLSBOT(L)
       END DO
      END IF
      

C
C

C
C**********************************************************************C
C
   99 FORMAT(A80)
  100 FORMAT(I10,F12.4)
  101 FORMAT(2I10)
  200 FORMAT(2I5,1X,6E14.6)
  299 FORMAT(2I5,1X,100E14.6)
  250 FORMAT(20E12.4)
  400 FORMAT(1X,6E14.6)
  420 FORMAT(1X,13E11.3)
cmrm  200 FORMAT(2I5,1X,1p,6E13.5) 
cmrm  250 FORMAT(1p,12E11.3)
C
C**********************************************************************C
C
      RETURN
      END
