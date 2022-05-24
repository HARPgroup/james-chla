C
C**********************************************************************C
C 
      SUBROUTINE WQ3DINP(LA,KC,IC,JC,IWQDT,DT,NTSPTC,IWQS)
C
C**********************************************************************C
C
C  READ WATER QUALITY SUBMODEL INPUT FILES
C
C  ORGINALLY CODED BY K.-Y. PARK
C  OPTIMIZED AND MODIFIED BY J. M. HAMRICK
C
C  LAST MODIFIED BY J. M. HAMRICK  7 APRIL 1997
C
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
      INTEGER, SAVE :: IOUT=1  
C
      CHARACTER*11  HHMMSS
c
      DATA IWQTICI,IWQTAGR,IWQTSTL,IWQTSUN,IWQTBEN,IWQTPSL,IWQTNPL/7*0/
      DATA ISMTICI/0/
      iwqtici=iwqtici
      iwqtagr=iwqtagr
      iwqtstl=iwqtstl
      iwqtsun=iwqtsun
      iwqtben=iwqtben
      iwqtpsl=iwqtpsl
      iwqtnpl=iwqtnpl
      ismtici=ismtici

      OPEN(1,FILE='wq3d.out',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
c
C **  HARDWIRE BY PASS OF RATE COEFFICIENT MAPS
C
c mrm      ISWQCMAP=0
c mrm      ISWQSMAP=0
C
      NWQKCNT=0
      NWQKDPT=1
C
!      UHEQ(1)=0.0
!      UHEQ(LC)=0.0
!      DO ND=1,NDMWQ
!       LF=2+(ND-1)*LDMWQ
!       LL=LF+LDM-1
!       DO L=LF,LL
!        UHEQ(L)=1.0
!       END DO
!      END DO
Ch-
      TINDAY = 0.0
      tinday=tinday
      ITNWQ = 0
c
      RKCWQ = 1.0/REAL(KC)
      DO K=1,KC
        WQHT(K)=REAL(KC-K)*RKCWQ
      END DO
!
      WQTSNAME(1)  = 'chl'  !Not used, Ji, 9/17/99
      WQTSNAME(2)  = 'tpc'
      WQTSNAME(3)  = 'doc'
      WQTSNAME(4)  = 'tpp'
      WQTSNAME(5)  = 'dop'
      WQTSNAME(6)  = 'p4t'
      WQTSNAME(7)  = 'p4d'
      WQTSNAME(8)  = 'apc'
      WQTSNAME(9)  = 'tnn'
      WQTSNAME(10) = 'don'
      WQTSNAME(11) = 'nh4'
      WQTSNAME(12) = 'no3'
      WQTSNAME(13) = 'tsi'
      WQTSNAME(14) = 'suu'
      WQTSNAME(15) = 'saa'
      WQTSNAME(16) = 'sad'
      WQTSNAME(17) = 'cod'
      WQTSNAME(18) = 'doo'
      WQTSNAME(19) = 'tam'
      WQTSNAME(20) = 'tmp'
      WQTSNAME(21) = 'fcb'
C

!      endif
c 1969 FORMAT('c   I    J    K    TIME       chl        tpc',
c     $       '        doc        tpp        dop        p4t',
c     $       '        p4d        apc        tnn        don',
c     $       '        nh4        no3        tsi        suu',
c     $       '        saa        sad        cod        doo',
c     $       '        tam        tmp        fcb        Malg')
c1969 FORMAT('c   I    J    K    TIME       nh4        no3',   ! for Wister, Ji, 7/27/99
c    $       '        ogn        tpp        p4t        chl',
c    $       '        bod        doo                      ',
c    $       '                                            ',
c    $       '                                            ',
c    $       '                                             ')
 1969 FORMAT('c   I    J    K    TIME       DO         CHL',   ! for TestWQ, Ji, 9/19/02
     $       '        TP         SRP        TKN        NH4',
     $       '        NOX        SI         Macroalg   TOC',
     $       '        BOD                                 ',
     $       '                                            ',
     $       '                                             ')
C

      DO M=0,NWQPSM
       WQPSQ(M)=0.0
       WQPSQC(M)=0.0
       DO J=1,NWQVM
        WQWPSLC(M,J)=0.0
       END DO
      END DO
C
      DO K=1,KC
       IWQPSC(1,K)=0
       WQDSQ(1,K)=0.0
       IWQPSC(LA+1,K)=0
       WQDSQ(LA+1,K)=0.0
      END DO
C

       DO K=1,KC
        DO L=2,LA
         IWQPSC(L,K)=0
         IWQPSV(L,K)=0
         WQDSQ(L,K)=0.0
        END DO
       END DO

C
      DO J=1,NWQVM
       DO K=1,KC
        DO L=1,LA+1
         WQWDSL(L,K,J)=0.0
         WQWPSL(L,K,J)=0.0
        END DO
       END DO
      END DO
C
	CALL WQMAPPING
	
      IF(IWQS.GT.0) THEN
         Write(*,*)'Runing simplified model'
	   CALL RWQINPS(IWQDT,DT,IC,JC,KC,LA,NTSPTC,IWQS)
	ELSE
	   Write(*,*)'Runing full version of the model'
	   CALL RWQINP(IWQDT,DT,IC,JC,KC,LA,NTSPTC)
	ENDIF
c
      write(*,*)'Write binay time series ',isTBIN
      
c      IF(isTBIN.EQ.0) THEN
c       OPEN(1,FILE='wqwcts.out',STATUS='UNKNOWN')
c       CLOSE(1,STATUS='DELETE')
c       OPEN(1,FILE='wqwcts.out',STATUS='UNKNOWN')   
c       WRITE(1,1969)
c       CLOSE(1)
c       if(isWQMIN==1) then
c        OPEN(1,FILE='wqwtsmin.out',STATUS='UNKNOWN')
c        CLOSE(1,STATUS='DELETE')
c        OPEN(1,FILE='wqwtsmin.out',STATUS='UNKNOWN')   
c        WRITE(1,1969)
c        CLOSE(1)
      
c        OPEN(1,FILE='wqwtsmax.out',STATUS='UNKNOWN')
c        CLOSE(1,STATUS='DELETE')
c        OPEN(1,FILE='wqwtsmax.out',STATUS='UNKNOWN')   
c        WRITE(1,1969)
c        CLOSE(1)     
c       endif
     
c      ELSE
C    
        OPEN(1,FILE='wqwcts.bin',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(551,FILE='wqwcts.bin',STATUS='UNKNOWN',form='binary')   
C 
       if(isWQMIN==1) then      
        OPEN(1,FILE='wqwtsmin.bin',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(552,FILE='wqwtsmin.bin',STATUS='UNKNOWN',form='binary')   
      
        OPEN(1,FILE='wqwtsmax.bin',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(553,FILE='wqwtsmax.bin',STATUS='UNKNOWN',form='binary')      

        OPEN(1,FILE='wqwtsavg.bin',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(554,FILE='wqwtsavg.bin',STATUS='UNKNOWN',form='binary')              
      
        OPEN(1,FILE='Algaeavg.bin',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(555,FILE='Algaeavg.bin',STATUS='UNKNOWN',form='binary')                     
        write(555)LCMWQ
        
       endif
c      ENDIF
C
C **  INITIALIZE DIURNAL DO ANALYSIS
C
      IF(NDDOAVG.GE.1) THEN
        OPEN(1,FILE='diurndo.out')
        CLOSE(1,STATUS='DELETE')
        DO K=1,KC
         DO L=2,LA
          DDOMAX(L,K)=-1.E6
          DDOMIN(L,K)=1.E6
         END DO
        END DO
      END IF
!
! **  INITIALIZE LIGHT EXTINCTION ANALYSIS
!
      IF(NDLTAVG.GE.1) THEN
        OPEN(1,FILE='light.out')
        CLOSE(1,STATUS='DELETE')
        NDLTCNT=0
        DO K=1,KC
         DO L=2,LA
          RLIGHTT(L,K)=0.
          RLIGHTC(L,K)=0.
         END DO
        END DO
      END IF
!
! ** Initialize water quality averaging summation arrays:
!
      call wqzero(LA,KC)

      IF (IWQICI.GE.2) CALL RWQRST(LA,KC)

      IF (IWQBEN.EQ.1) THEN
        CALL SMINIT(LA)
	  CALL SMRIN1(IWQDT,DT,IC,JC,KC,LA,ITIMES,NTSPTC)
        IF (ISMICI.EQ.2) CALL RSMRST(LA)
      END IF

      RETURN
      END
