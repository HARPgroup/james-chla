C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RVELPLTH_ac
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
      OPEN(11,FILE='rvelvch'//iyear//'.out',STATUS='UNKNOWN')
      OPEN(12,FILE='pvelvch'//iyear//'.out',STATUS='UNKNOWN')
      OPEN(13,FILE='mvelvch'//iyear//'.out',STATUS='UNKNOWN')
      CLOSE(11,STATUS='DELETE')
      CLOSE(12,STATUS='DELETE')
      CLOSE(13,STATUS='DELETE')
      OPEN(11,FILE='rvelvch'//iyear//'.out',STATUS='UNKNOWN')
      OPEN(12,FILE='pvelvch'//iyear//'.out',STATUS='UNKNOWN')
      OPEN(13,FILE='mvelvch'//iyear//'.out',STATUS='UNKNOWN')
      WRITE (11,99) TITLE1
      WRITE (11,101)LINES,LEVELS
      WRITE (11,250)(DBS(L),L=1,LEVELS)
      WRITE (12,99) TITLE2
      WRITE (12,101)LINES,LEVELS
      WRITE (12,250)(DBS(L),L=1,LEVELS)
      WRITE (13,99) TITLE3
      WRITE (13,101)LINES,LEVELS
      WRITE (13,250)(DBS(L),L=1,LEVELS)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
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
C
      OPEN(11,FILE='rvelvch'//iyear//'.out',ACCESS='APPEND',
	1STATUS='UNKNOWN')
c      OPEN(12,FILE='pvelvch'//iyear//'.out',ACCESS='APPEND', !hong comment file 12 and 13
c	1STATUS='UNKNOWN')
c      OPEN(13,FILE='mvelvch'//iyear//'.out',ACCESS='APPEND',
c	1STATUS='UNKNOWN')
      WRITE (11,100)N,TIME
c      WRITE (12,100)N,TIME
c      WRITE (13,100)N,TIME
C
      IF(ISRVPH.EQ.1) THEN
       DO L=2,LA
       LN=LNC(L)
chong
c       UTMP=50.*STCUV(L)*(UHLPF(L+1,KC)+UHLPF(L,KC))/HMP(L)
c       VTMP=50.*STCUV(L)*(VHLPF(LN,KC)+VHLPF(L,KC))/HMP(L)
c       RVELEKC=CUE(L)*UTMP+CVE(L)*VTMP
c       RVELNKC=CUN(L)*UTMP+CVN(L)*VTMP
c       UTMP=50.*STCUV(L)*(UHLPF(L+1,1)+UHLPF(L,1))/HMP(L)
c       VTMP=50.*STCUV(L)*(VHLPF(LN,1)+VHLPF(L,1))/HMP(L)
c       RVELEKB=CUE(L)*UTMP+CVE(L)*VTMP
c       RVELNKB=CUN(L)*UTMP+CVN(L)*VTMP
c       WRITE(11,200)IL(L),JL(L),DLON(L),DLAT(L),RVELEKC,RVELNKC,
c     $              RVELEKB,RVELNKB
         do k=1,kc
       UTMP=50.*STCUV(L)*(UHLPF(L+1,K)+UHLPF(L,K))/HMP(L)
       VTMP=50.*STCUV(L)*(VHLPF(LN,K)+VHLPF(L,K))/HMP(L)
c       RVELE(K)=CUE(L)*UTMP+CVE(L)*VTMP
c       RVELN(K)=CUN(L)*UTMP+CVN(L)*VTMP
       RVELE(K)=UTMP   
       RVELN(K)=VTMP    !hb velocity without projection
         enddo
       WRITE(11,299)IL(L),JL(L),(RVELE(K),k=1,kc),(RVELN(K),k=1,kc)
     

c       UTMP=50.*STCUV(L)*(UVPT(L+1,KC)+UVPT(L,KC))/HMP(L)
c       VTMP=50.*STCUV(L)*(VVPT(LN,KC)+VVPT(L,KC))/HMP(L)
c       PVELEKC=CUE(L)*UTMP+CVE(L)*VTMP
c       PVELNKC=CUN(L)*UTMP+CVN(L)*VTMP
c       UTMP=50.*STCUV(L)*(UVPT(L+1,1)+UVPT(L,1))/HMP(L)
c       VTMP=50.*STCUV(L)*(VVPT(LN,1)+VVPT(L,1))/HMP(L)
c       PVELEKB=CUE(L)*UTMP+CVE(L)*VTMP
c       PVELNKB=CUN(L)*UTMP+CVN(L)*VTMP
c       WRITE(12,200)IL(L),JL(L),DLON(L),DLAT(L),PVELEKC,PVELNKC,
c     $              PVELEKB,PVELNKB
c       PVELEKC=PVELEKC+RVELEKC
c       PVELNKC=PVELNKC+RVELNKC
c       PVELEKB=PVELEKB+RVELEKB
c       PVELNKB=PVELNKB+RVELNKB
c       WRITE(13,200)IL(L),JL(L),DLON(L),DLAT(L),PVELEKC,PVELNKC,
c     $              PVELEKB,PVELNKB
       END DO
      END IF
C
      IF(ISRVPH.EQ.2) THEN
       DO LR=1,NRC
       L=LRC(LR)
       LN=LNC(L)
       UTMP=50.*STCUV(L)*(UHLPF(L+1,KC)+UHLPF(L,KC))/HMP(L)
       VTMP=50.*STCUV(L)*(VHLPF(LN,KC)+VHLPF(L,KC))/HMP(L)
       RVELEKC=CUE(L)*UTMP+CVE(L)*VTMP
       RVELNKC=CUN(L)*UTMP+CVN(L)*VTMP
       UTMP=50.*STCUV(L)*(UHLPF(L+1,1)+UHLPF(L,1))/HMP(L)
       VTMP=50.*STCUV(L)*(VHLPF(LN,1)+VHLPF(L,1))/HMP(L)
       RVELEKB=CUE(L)*UTMP+CVE(L)*VTMP
       RVELNKB=CUN(L)*UTMP+CVN(L)*VTMP
       WRITE(11,200)IL(L),JL(L),DLON(L),DLAT(L),RVELEKC,RVELNKC,
     $              RVELEKB,RVELNKB
       UTMP=50.*STCUV(L)*(UVPT(L+1,KC)+UVPT(L,KC))/HMP(L)
       VTMP=50.*STCUV(L)*(VVPT(LN,KC)+VVPT(L,KC))/HMP(L)
       PVELEKC=CUE(L)*UTMP+CVE(L)*VTMP
       PVELNKC=CUN(L)*UTMP+CVN(L)*VTMP
       UTMP=50.*STCUV(L)*(UVPT(L+1,1)+UVPT(L,1))/HMP(L)
       VTMP=50.*STCUV(L)*(VVPT(LN,1)+VVPT(L,1))/HMP(L)
       PVELEKB=CUE(L)*UTMP+CVE(L)*VTMP
       PVELNKB=CUN(L)*UTMP+CVN(L)*VTMP
       WRITE(12,200)IL(L),JL(L),DLON(L),DLAT(L),PVELEKC,PVELNKC,
     $              PVELEKB,PVELNKB
       PVELEKC=PVELEKC+RVELEKC
       PVELNKC=PVELNKC+RVELNKC
       PVELEKB=PVELEKB+RVELEKB
       PVELNKB=PVELNKB+RVELNKB
       WRITE(13,200)IL(L),JL(L),DLON(L),DLAT(L),PVELEKC,PVELNKC,
     $              PVELEKB,PVELNKB
       END DO
      END IF
C
      IF(ISRVPH.EQ.3) THEN
       DO LB=1,NBC
       L=LBC(LB)
       LN=LNC(L)
       UTMP=50.*STCUV(L)*(UHLPF(L+1,KC)+UHLPF(L,KC))/HMP(L)
       VTMP=50.*STCUV(L)*(VHLPF(LN,KC)+VHLPF(L,KC))/HMP(L)
       RVELEKC=CUE(L)*UTMP+CVE(L)*VTMP
       RVELNKC=CUN(L)*UTMP+CVN(L)*VTMP
       UTMP=50.*STCUV(L)*(UHLPF(L+1,1)+UHLPF(L,1))/HMP(L)
       VTMP=50.*STCUV(L)*(VHLPF(LN,1)+VHLPF(L,1))/HMP(L)
       RVELEKB=CUE(L)*UTMP+CVE(L)*VTMP
       RVELNKB=CUN(L)*UTMP+CVN(L)*VTMP
       WRITE(11,200)IL(L),JL(L),DLON(L),DLAT(L),RVELEKC,RVELNKC,
     $              RVELEKB,RVELNKB
       UTMP=50.*STCUV(L)*(UVPT(L+1,KC)+UVPT(L,KC))/HMP(L)
       VTMP=50.*STCUV(L)*(VVPT(LN,KC)+VVPT(L,KC))/HMP(L)
       PVELEKC=CUE(L)*UTMP+CVE(L)*VTMP
       PVELNKC=CUN(L)*UTMP+CVN(L)*VTMP
       UTMP=50.*STCUV(L)*(UVPT(L+1,1)+UVPT(L,1))/HMP(L)
       VTMP=50.*STCUV(L)*(VVPT(LN,1)+VVPT(L,1))/HMP(L)
       PVELEKB=CUE(L)*UTMP+CVE(L)*VTMP
       PVELNKB=CUN(L)*UTMP+CVN(L)*VTMP
       WRITE(12,200)IL(L),JL(L),DLON(L),DLAT(L),PVELEKC,PVELNKC,
     $              PVELEKB,PVELNKB
       PVELEKC=PVELEKC+RVELEKC
       PVELNKC=PVELNKC+RVELNKC
       PVELEKB=PVELEKB+RVELEKB
       PVELNKB=PVELNKB+RVELNKB
       WRITE(13,200)IL(L),JL(L),DLON(L),DLAT(L),PVELEKC,PVELNKC,
     $              PVELEKB,PVELNKB
       END DO
      END IF
C
      CLOSE(11)
chong      CLOSE(12)
chong      CLOSE(13)
C
C**********************************************************************C
C
   99 FORMAT(A80)
  100 FORMAT(I10,F12.4)
  101 FORMAT(2I10)
  200 FORMAT(2I5,1X,6E14.6)
  299 FORMAT(2I5,1X,20E14.6,20E14.6)

  250 FORMAT(20E12.4)
cmrm  200 FORMAT(2I5,1X,1p,6E13.5) 
cmrm  250 FORMAT(1p,12E11.3)
C
C**********************************************************************C
C
      RETURN
      END
