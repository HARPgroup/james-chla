C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C 
      SUBROUTINE WSMTS(TINDAY,DT,N,TCON,TBEGIN,NTSPTC,NCSTEP,SECDLAST)
C
C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK ON  7 APRIL 1997
C
C**********************************************************************C
C
C Write time-series output
C
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
      DIMENSION TMPOUT(30)
      REAL*8 SECDLAST  
      TIMTMP=(DT*FLOAT(N-1)+TCON*TBEGIN)/86400  
      IF(NCSTEP.GT.0) TIMTMP=SECDLAST/TCON+TBEGIN 
          
      tinday=tinday
C
      IF(isTBIN.EQ.0) THEN
      OPEN(1,FILE='wqsdts1.out',STATUS='UNKNOWN',ACCESS='APPEND')
      OPEN(2,FILE='wqsdts2.out',STATUS='UNKNOWN',ACCESS='APPEND')
C
!      TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
!      TIMTMP=(DT*FLOAT(N-1)+TCON*TBEGIN)/86400
C
      DO M=1,ISMTS
        LL=LSMTS(M)
        TSMPON = SMPON(LL,1)+SMPON(LL,2)+SMPON(LL,3)
        TSMPOP = SMPOP(LL,1)+SMPOP(LL,2)+SMPOP(LL,3)
        TSMPOC = SMPOC(LL,1)+SMPOC(LL,2)+SMPOC(LL,3)
        WRITE(1,71)ILW(LL),JLW(LL),TIMTMP,SM1NH4(LL),SM2NH4(LL),
     *   WQBFNH4(LL),SM1NO3(LL),SM2NO3(LL),WQBFNO3(LL),SM1PO4(LL),
     *   SM2PO4(LL),WQBFPO4D(LL),SM1H2S(LL),SM2H2S(LL),WQBFO2(LL),
     *   WQBFCOD(LL),SM1SI(LL),SM2SI(LL),WQBFSAD(LL),SMT(LL),SMBST(LL),
     *   TSMPON,TSMPOP,TSMPOC
        WRITE(2,71)ILW(LL),JLW(LL),TIMTMP,SMCSOD(LL),SMNSOD(LL),
     *   SMD1PO4(LL),SMD1SI(LL),SMSS(LL),SMJNIT(LL),SMJDEN(LL),
     *   SMJAQH2S(LL),SMJGCH4(LL),SMDGFN(LL),SMDGFP(LL),SMDGFC(LL),
     *   SMDFN(LL,1),SMDFN(LL,2),SMDFN(LL,3),SMDFP(LL,1),SMDFP(LL,2),
     *   SMDFP(LL,3),SMDFC(LL,1),SMDFC(LL,2),SMDFC(LL,3)       ! Why silica fluxes are not calculated and not given? 10/9/02
      END DO
C
cmrm   71 FORMAT(2I5, F11.5, 23E12.4)
   71 FORMAT(2I5, F11.5, 1p, 23E11.3)
C
      CLOSE(1)
      CLOSE(2)
      
      ELSE
      write(550)TIMTMP
      DO M=1,ISMTS
      LL=LSMTS(M)  
      TSMPON = SMPON(LL,1)+SMPON(LL,2)+SMPON(LL,3)
      TSMPOP = SMPOP(LL,1)+SMPOP(LL,2)+SMPOP(LL,3)
      TSMPOC = SMPOC(LL,1)+SMPOC(LL,2)+SMPOC(LL,3)
      TMPOUT(1)=SMPON(LL,1)
      TMPOUT(2)=SMPON(LL,2)
      TMPOUT(3)=SMPON(LL,3) 
      TMPOUT(4)=SMPOP(LL,1)
      TMPOUT(5)=SMPOP(LL,2)
      TMPOUT(6)=SMPOP(LL,3)
      TMPOUT(7)=SMPOC(LL,1)
      TMPOUT(8)=SMPOC(LL,2)
      TMPOUT(9)=SMPOC(LL,3)
      TMPOUT(10)=SM1NH4(LL)
      TMPOUT(11)=SM2NH4(LL)
      TMPOUT(12)=WQBFNH4(LL)
      TMPOUT(13)=SM1NO3(LL)
      TMPOUT(14)=SM2NO3(LL)
      TMPOUT(15)=WQBFNO3(LL)
      TMPOUT(16)=SM1PO4(LL)
      TMPOUT(17)=SM2PO4(LL)
      TMPOUT(18)=WQBFPO4D(LL)
      TMPOUT(19)=SM1H2S(LL)
      TMPOUT(20)=SM2H2S(LL)
      TMPOUT(21)=WQBFO2(LL)
      TMPOUT(22)=WQBFCOD(LL)
      TMPOUT(23)=SM1SI(LL)
      TMPOUT(24)=SM2SI(LL)
      TMPOUT(25)=WQBFSAD(LL)
      TMPOUT(26)=SMT(LL)
      TMPOUT(27)=SMBST(LL)
      TMPOUT(28)=TSMPON
      TMPOUT(29)=TSMPOP
      TMPOUT(30)=TSMPOC  
      write(550)TMPOUT     
      ENDDO
      ENDIF 
C
      RETURN
      END
