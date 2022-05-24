C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WWQRST(LA,KC,DT,N,TCON,TBEGIN,NCSTEP,SECDLAST)
C
C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 10 april 1999
C
C**********************************************************************C
C
C Write spatial distributions at the end of simulation to unit IWQORST.
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
!       IF (mod(ITIME,385).EQ.0) THEN 
       OPEN(1,FILE='wqwcrst.out',STATUS='UNKNOWN')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='wqwcrst.out',STATUS='UNKNOWN')
C

!      WRITE(1,101) N,TIME
       WRITE(1,*) N,TIME
       WRITE(1,102)
C
C J.S.
      NWQV0=NWQV
      IF(IDNOTRVA.GT.0) NWQV0=NWQV0+1
      DO L=2,LA
        DO K=1,KC
          WRITE(1,90) L,K,(WQV(L,K,NW),NW=1,NWQV0)
        END DO
      END DO
C
      CLOSE(1)
!      ENDIF    

!      IF(mod(ITIME,IWQRST).EQ.0) THEN
       write(dayt,'(i4.4)')ITIME
      OPEN(1,FILE='wqwcrst'//dayt//'.out',STATUS='UNKNOWN')
      WRITE(1,*) N,TIME
      WRITE(1,102)
      NWQV0=NWQV
      IF(IDNOTRVA.GT.0) NWQV0=NWQV0+1
      DO L=2,LA
        DO K=1,KC
          WRITE(1,90) L,K,(WQV(L,K,NW),NW=1,NWQV0)
        END DO
      END DO
      CLOSE(1)
!      ENDIF
c
c Also write Binary restart file:
c
c      INQUIRE(FILE='WQWCrst.bin', EXIST=fexist)
c      if (fexist) then
c        OPEN(UNIT=1, FILE='WQWCrst.bin', ACCESS='transparent',
c     +    FORM='unformatted', STATUS='unknown')
c        CLOSE(UNIT=1, STATUS='DELETE')
c      end if
c      open(UNIT=1, FILE='WQWCrst.bin', ACCESS='transparent',
c     +   FORM='unformatted', STATUS='unknown')
c      write(1) N, TIME
c      NWQV0=NWQV
c      IF(IDNOTRVA.GT.0) NWQV0=NWQV0+1
c      DO L=2,LA
c        DO K=1,KC
c          write(1) L, K
c          do nw=1,nwqv0
c            write(1) wqv(L,K,NW)
c          end do
c        END DO
c      END DO
c      CLOSE(1)
C
   90 FORMAT(2I5, 1p, 22E12.4)
  101 FORMAT('CC  WQ RESTART FILE TIME STEP, TIME = ',I10,F12.5)
  102 FORMAT('c   L    K  bc          bd          bg          ',
     $       'rpoc        lpoc        doc         ',
     $       'rpop        lpop        dop         pto4        ',
     $       'rpon        lpon        don         amn         ',
     $       'nit         su          sa          cod         ',
     $       'do          tam         fcb        Malg')
C
      RETURN
      END
