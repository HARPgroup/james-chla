C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WSMRST(LA,DT,N,TCON,TBEGIN,NCSTEP,SECDLAST)
C
C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 10 april 1999
C
C**********************************************************************C
C
C Write spatial distributions at the end of simulation to unit ISMORST.
C
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
      logical fexist
      character*4 dayt
      REAL*8 SECDLAST  
C
c Write ASCII restart file:
c
       TIME=TBEGIN+DT*FLOAT(N)/TCON
       IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN 
       ITIME=int(TIME)
!       IF (mod(ITIME,365).EQ.0) THEN
       OPEN(1,FILE='wqsdrst.out',STATUS='UNKNOWN')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='wqsdrst.out',STATUS='UNKNOWN')
C

       WRITE(1,101) N,TIME
       WRITE(1,888)
C
       DO L=2,LA
        WRITE(1,90) L,(SMPON(L,NW),NW=1,NSMG),
     *    (SMPOP(L,NW),NW=1,NSMG),(SMPOC(L,NW),NW=1,NSMG),SM1NH4(L),
     *    SM2NH4(L),SM2NO3(L),SM2PO4(L),SM2H2S(L),SMPSI(L),SM2SI(L),
     *    SMBST(L),SMT(L)
      END DO
C
      CLOSE(1)
!      ENDIF
      

!      IF(mod(ITIME,IWQRST).EQ.0) THEN
       write(dayt,'(i4.4)')ITIME
      OPEN(1,FILE='wqsdrst'//dayt//'.out',STATUS='UNKNOWN')       
      WRITE(1,*) N,TIME
      WRITE(1,888)
      DO L=2,LA
        WRITE(1,90) L,(SMPON(L,NW),NW=1,NSMG),
     *    (SMPOP(L,NW),NW=1,NSMG),(SMPOC(L,NW),NW=1,NSMG),SM1NH4(L),
     *    SM2NH4(L),SM2NO3(L),SM2PO4(L),SM2H2S(L),SMPSI(L),SM2SI(L),
     *    SMBST(L),SMT(L)
      END DO
      close(1)
!      ENDIF
c
c Also write Binary restart file:
c
c      INQUIRE(FILE='WQSDrst.bin', EXIST=fexist)
c      if (fexist) then
c        OPEN(UNIT=1, FILE='WQSDrst.bin', ACCESS='transparent',
c     +    FORM='unformatted', STATUS='unknown')
c        CLOSE(UNIT=1, STATUS='DELETE')
c      end if
c      open(UNIT=1, FILE='WQSDrst.bin', ACCESS='transparent',
c     +   FORM='unformatted', STATUS='unknown')
c      write(1) N, TIME
c      DO L=2,LA
c        write(1) L
c        WRITE(1) (SMPON(L,NW),NW=1,NSMG),
c     *    (SMPOP(L,NW),NW=1,NSMG),(SMPOC(L,NW),NW=1,NSMG),SM1NH4(L),
c     *    SM2NH4(L),SM2NO3(L),SM2PO4(L),SM2H2S(L),SMPSI(L),SM2SI(L),
c     *    SMBST(L),SMT(L)
c      END DO
c      CLOSE(1)
C
c   90 FORMAT(I5, 18E12.4)
   90 FORMAT(I5, 1p, 18E12.4)
  101 FORMAT('CC  SM RESTART FILE TIME STEP, TIME = ',I10,F13.5)
  888 FORMAT('    L',
     *'       GPON1       GPON2       GPON3       GPOP1       GPOP2',
     *'       GPOP3       GPOC1       GPOC2       GPOC3       G1NH4',
     *'       G2NH4       G2NO3       G2PO4       G2H2S        GPSI',
     *'        G2SI        GBST          GT')
C
      RETURN
      END
