C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RVELPLTH
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE RVELPLTH WRITES HORIZONTAL EULERIAN RESIDUAL, VECTOR 
C **  POTENTIAL AND MEAN MASS TRANSPORT VELOCITY VECTOR FILES 
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      REAL DBS(10)
      CHARACTER*80 TITLE1,TITLE2,TITLE3
	real RVELE(KC),RVELN(KC)
C
C**********************************************************************C
C
      IF (JSRVPH.NE.1) GO TO 300
C
C----------------------------------------------------------------------C
C
C **  WRITE HEADINGS
C
      CLOSE(541)
      TITLE1='HORIZ EULERIAN MEAN TRANSPORT VELOCITY'
      TITLE2='HORIZ VECTOR POTENTIAL TRANSPORT VELOCITY'
      TITLE3='HORIZ MEAN MASS TRANSPORT VELOCITY'
      IF(ISRVPH.EQ.1) LINES=LA-1
      IF(ISRVPH.EQ.2) LINES=NRC
      IF(ISRVPH.EQ.3) LINES=NBC
      LEVELS=2
      DBS(1)=0.
      DBS(2)=99.
C
      OPEN(541,FILE='rvelvch'//iyear//'.bin',STATUS='UNKNOWN',
     &  form='binary')

!      OPEN(542,FILE='pvelvch'//iyear//'.bin',STATUS='UNKNOWN')
!      OPEN(543,FILE='mvelvch'//iyear//'.bin',STATUS='UNKNOWN')
      CLOSE(541,STATUS='DELETE')
!      CLOSE(542,STATUS='DELETE')
!      CLOSE(543,STATUS='DELETE')
      OPEN(541,FILE='rvelvch'//iyear//'.bin',form='binary')
!      OPEN(542,FILE='pvelvch'//iyear//'.bin',form='binary')
!      OPEN(543,FILE='mvelvch'//iyear//'.bin',form='binary')
      write(541)LA-1,KC
C
      JSRVPH=0
C
C----------------------------------------------------------------------C
C
  300 CONTINUE
C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON 
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014          
      
      WRITE (541)TIME
c      WRITE (12,100)N,TIME
c      WRITE (13,100)N,TIME
C
      IF(ISRVPH.EQ.1) THEN
       DO L=2,LA
       LN=LNC(L)
        DO k=1,kc
        UTMP=50.*STCUV(L)*(UHLPF(L+1,K)+UHLPF(L,K))/HMP(L)
        VTMP=50.*STCUV(L)*(VHLPF(LN,K)+VHLPF(L,K))/HMP(L)
        RVELE(K)=UTMP   
        RVELN(K)=VTMP    !hb velocity without projection
        ENDDO
       WRITE(541)IL(L),JL(L),(RVELE(K),k=1,kc)
       WRITE(541)(RVELN(K),k=1,kc)   
       END DO
      END IF


      IF(ISRVPH.EQ.2) THEN
       DO L=2,LA
       LN=LNC(L)

       DO k=1,kc
       UTMP=50.*STCUV(L)*(UHLPF(L+1,K)+UHLPF(L,K))/HMP(L)
       VTMP=50.*STCUV(L)*(VHLPF(LN,K)+VHLPF(L,K))/HMP(L)
       RVELE(K)=CUE(L)*UTMP+CVE(L)*VTMP
       RVELN(K)=CUN(L)*UTMP+CVN(L)*VTMP
       ENDDO
       WRITE(541)IL(L),JL(L),(RVELE(K),k=1,kc)
       WRITE(541)(RVELN(K),k=1,kc) 
       ENDDO     
      ENDIF
       
C
C
C**********************************************************************C
C
      RETURN
      END
