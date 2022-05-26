C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WSMTSbin(N,LA,iyear1)
C
C**********************************************************************C
C **  LAST MODIFIED BY M.R. MORTON ON 02 JUNE 1999
C**********************************************************************C
C Write sediment time-series output to binary file.
C Averages benthic flux rates over ISMTSDT time steps (e.g., daily avg).
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
      integer flength
      INTEGER, SAVE :: inisd=1,iyearss=0,inibins=1,inisd1=1
      character*4 iyears
c
      IF(iyearss.ne.iyear1) then
        write(iyears,'(I4)')iyear1
        write(*,*)'WQ SAVE YEAR ',iyear
        iyearss=iyear1
        inisd=1
      ENDIF

	flength=1 ! Ji, 9/17/99
c
      if(inisd.eq.1) then
       open(1,file='WQSD3D'//iyears//'.bin',form='unformatted') 
       CLOSE(1,STATUS='DELETE')
       open(1,file='WQSD3D'//iyears//'.bin',form='unformatted') 
       CLOSE(1)
       inisd=0
      endif
      if(inisd1.eq.1)then
       open(1,file='WQSDTS.out') 
       CLOSE(1,STATUS='DELETE')
       open(1,file='WQSDTS.out') 
       write(1,*)'TIM,BFO2,BFNH4,BFNO3,BFPO4,BFSAD,BFCOD'
       close(1)
       inisd1=0
      endif
      
      if (isSDBIN .gt. 0) then
        IF ( MOD(ITNWQ,IAVGBIN) .EQ. 0 ) THEN
          nrec4 = nrec4+1
          TIMTMP = timebf / float(nbfcnt)
          write(*,*)'Out flux time= ',TIMTMP,ITNWQ
          OPEN(UNIT=2, FILE='WQSD3D'//iyears//'.bin',ACCESS='append',
     +       form='unformatted',STATUS='unknown')
          OPEN(UNIT=3, FILE='WQSDTS.out',ACCESS='append',
     +        STATUS='unknown')
c         INQUIRE(UNIT=2, FLEN=flength)
          write(2) TIMTMP
          DO LL=2,LA
            BFO2sum(LL)  = BFO2sum(LL)  / float(nbfcnt)
            BFNH4sum(LL) = BFNH4sum(LL) / float(nbfcnt)
            BFNO3sum(LL) = BFNO3sum(LL) / float(nbfcnt)
            BFPO4sum(LL) = BFPO4sum(LL) / float(nbfcnt)
            BFSADsum(LL) = BFSADsum(LL) / float(nbfcnt)
            BFCODsum(LL) = BFCODsum(LL) / float(nbfcnt)
!            BFO2sum(LL)  = BFO2sum(LL)  
!            BFNH4sum(LL) = BFNH4sum(LL) 
!            BFNO3sum(LL) = BFNO3sum(LL) 
!            BFPO4sum(LL) = BFPO4sum(LL) 
!            BFSADsum(LL) = BFSADsum(LL) 
!            BFCODsum(LL) = BFCODsum(LL) 
          END DO
          WRITE(2) BFO2sum
		  write(2) BFNH4sum
		  write(2) BFNO3sum
          write(2) BFPO4sum
		  write(2) BFSADsum
		  write(2) BFCODsum
          CLOSE(2)
          DO M=1,IWQTS                              
            LL=LWQTS(M)
            write(3,3000)TIMTMP,BFO2sum(LL),BFNH4sum(LL),BFNO3sum(LL),
     +       BFPO4sum(LL),BFSADsum(LL), BFCODsum(LL)
          ENDDO
          close(3)
          call WQzero4(LA)
 
        END IF
      end if
3000  Format(7f10.3)
C
      RETURN
      END
