C 
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RWQRST(LA,KC)
C
C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 10 april 1999
C
C**********************************************************************C
C
C Read ICs from restart file from INWQRST.
C
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
      logical fexist
      integer nn
      real xtime
C
c Check first to see if binary restart file exists.  If not, use
c the ASCII file instead.
c
      LK=(LA-1)*KC
C
      INQUIRE(FILE='WQWCrst.bin', EXIST=fexist)
      if (.not. fexist) then
        OPEN(1,FILE='wqwcrst.inp',STATUS='UNKNOWN')
        OPEN(2,FILE='wqwcrst_n.inp',STATUS='UNKNOWN')        
        READ(1,999)
        READ(1,999)
C J.S.
        NWQV0=NWQV
        write(*,*)'Read restart wqwcrst.inp ' 
        IF(IDNOTRVA.GT.0) NWQV0=NWQV0+1
        DO L1=2,LA
         DO K1=1,KC 
          READ(1,* ) L,K,(WQV(L,K,NW),NW=1,NWQV0)
         END DO
!         write(*,*)'L= ', L1
          DO K1=1,KC
          write(2,121 ) L,K,(WQV(L,K,NW),NW=1,NWQV0)
          END DO              
	  ENDDO
	   IF(IDNOTRVA.GT.0)THEN
	   DO L=2,LA
         IF(SMAC(L).GT.0.1.and.WQV(L,1,22).LT.0.0001)THEN
         WQV(L,1,22)=10
         ENDIF
         ENDDO
         ENDIF
C
121      format(2I5,24E12.4)
        CLOSE(1)
         CLOSE(2)       
      else
c        open(UNIT=1, FILE='WQWCrst.bin', ACCESS='transparent',
        open(UNIT=1, FILE='WQWCrst.bin', ACCESS='sequential',
     +     FORM='unformatted', STATUS='unknown')
        read(1) nn, xtime
        xtime=xtime
        write(0,911) nn, xtime
911     format(' Reading binary WQWCRST.BIN file ...    NN, TIME = ',
     +     I7, F11.5)
        NWQV0=NWQV
        IF(IDNOTRVA.GT.0) NWQV0=NWQV0+1
        DO M=1,LK
          READ(1) L, K
          do nw=1,nwqv0
            read(1) WQV(L, K, NW)
          end do
        END DO
        CLOSE(1)
      end if
C
c Initialize macroalgae so biomass only exists in bottom layer:
      IF(IDNOTRVA.GT.0) THEN
C      DO L=2,LA
C       WQV(L,1,22)=SMAC(LL)*WQV(L,K,22)
C      END DO
        IF(KC.GT.1) THEN
          DO K=2,KC
           DO L=2,LA
            WQV(L,k,22)=0.
           END DO
          END DO
        END IF
      END IF
C
   90 FORMAT(2I5, 21E12.4)
  999 FORMAT(1X)
C
      RETURN
      END
