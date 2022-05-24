!**********************************************************************C
! 
      SUBROUTINE WWQTS(TINDAY,KC,DT,N,TCON,TBEGIN,NTSPTC,
     & NCSTEP,SECDLAST)
!
!**********************************************************************C
!
! **  LAST MODIFIED BY JOHN HAMRICK ON  7 APRIL 1997
! **  LAST MODIFIED BY AEE 3/5/2007
!
!**********************************************************************C
!
! Write time-series output: WQCHLx=1/WQCHLx
!
!
!  B  B  B  R  L  D  R  L  D  P  R  L  D  N  N  S  S  C  D  T  F  M
!  c  d  g  P  P  O  P  P  O  O  P  P  O  H  O  U  A  O  O  A  C  A
!           O  O  C  O  O  P  4  O  O  N  4  3        D     M     C
!           C  C     P  P     t  N  N                             A
!  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
!**********************************************************************C
!
      INCLUDE 'wq.par' 
      INCLUDE 'wqcom.cmn'
      DIMENSION TMPOUT(30),CHtemp(LCMWQ)
      REAL*8 SECDLAST         
      tinday=tinday
      
 !      ITRYCHLA=1  ! Use cary C:Chla ratio (Cerco, 2004) HAB 4/9/2015
      IF(ITRYCHLA.EQ.3) THEN
	KECST=30           ! also active respiration with growth
	ELSEIF(ITRYCHLA.EQ.2)THEN
	KECST=40
	ENDIF
! 
      TIMTMP=(DT*FLOAT(N-1)+TCON*TBEGIN)/86400
      IF(NCSTEP.GT.0) TIMTMP=SECDLAST/TCON+TBEGIN 
        write(*,*)'Time to write time series ',  TIMTMP      
      IF(isTBIN.EQ.0) THEN
       OPEN(1,FILE='wqwcts.out',STATUS='UNKNOWN',ACCESS='APPEND')
      ELSE
      WRITE(551)TIMTMP
       IF (isWQMIN==1) THEN
        if(isWQAVG.EQ.0) then
        WRITE(552)TIMTMP
        WRITE(553)TIMTMP
        endif
        WRITE(554)TIMTMP
       ENDIF
      ENDIF

      DO L=2,LCMWQ
       CHtemp(L) =WQ3DA(L,KC,1)+
     &            WQ3DA(L,KC,2) +  WQ3DA(L,KC,3) 
      ENDDO
      write(555)TIMTMP
      write(555)CHtemp
!      write(*,*)'Chl = ',TIMTMP,CHtemp(100),WQ3DA(100,KC,2), 
!     *    WQCHLG(IWQZMAP(100,KC))
!
      DO M=1,IWQTS                             ! #200
        DO K=1,KC                              ! #100
        LL=LWQTS(M)
        
        IF(ITRYCHLA.GE.2) THEN             !JS  2018
        WQCHLD(IWQZMAP(LL,KC))=1/(KECST+90*exp(-1.19*WQketot(LL,KC)))
     &  *1000      
        WQCHLC(IWQZMAP(LL,KC))=1/(KECST+150*exp(-1.18*WQketot(LL,KC)))
     &  *1000  
        WQCHLG(IWQZMAP(LL,KC))=WQCHLC(IWQZMAP(LL,KC)) 
        ENDIF
             
         CHLWQ = WQVO(LL,K,1)*WQCHLC(IWQZMAP(LL,KC)) +
     *       WQVO(LL,K,2)*WQCHLD(IWQZMAP(LL,KC))             ! total algae as chlorophyll, ug/l
     *      + WQVO(LL,K,3)*WQCHLG(IWQZMAP(LL,KC))
        
          TBWQ = WQVO(LL,K,1)+WQVO(LL,K,2)+WQVO(LL,K,3)                 ! total algae as C, mg/l
          TOCWQ = WQVO(LL,K,4)+WQVO(LL,K,5)+WQVO(LL,K,6) + TBWQ         ! TOC
          IF (IWQSRP.EQ.1) THEN
            O2WQ = MAX(WQVO(LL,K,19), 0.0)
            TAMDWQ = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ), WQVO(LL,K,20) )
            TAMPWQ = WQVO(LL,K,20) - TAMDWQ
            PO4DWQ = WQVO(LL,K,10) / (1.0 + WQKPO4P*TAMPWQ)
            SADWQ = WQVO(LL,K,17) / (1.0 + WQKSAP*TAMPWQ)
           ELSE IF (IWQSRP.EQ.2) THEN                                   ! SEDT
            PO4DWQ = WQVO(LL,K,10) / (1.0 + WQKPO4P*SEDTWQ(LL,K))       !=> Dissolved PO4t, Eq.(3-8c)
            SADWQ = WQVO(LL,K,17) / (1.0 + WQKSAP*SEDTWQ(LL,K))         !=> Dissolved SA  , Eq.(3-15c)
           ELSE                                                         ! SEDT
            PO4DWQ = WQVO(LL,K,10)
            SADWQ  = WQVO(LL,K,17)
          END IF
c
          XPO4DWQ = MAX(PO4DWQ,0.0)
          APCWQ = 1.0 / (WQCP1PRM1(IMWQZT(LL))
     &      + WQCP2PRM*EXP(-WQCP3PRM*XPO4DWQ))
          TPWQ = WQVO(LL,K,7)+WQVO(LL,K,8)+WQVO(LL,K,9)+WQVO(LL,K,10)    ! included algae's
     *      + APCWQ*TBWQ
          TNWQ = WQVO(LL,K,11)+WQVO(LL,K,12)+WQVO(LL,K,13)+WQVO(LL,K,14) ! included algae's
     *    +WQVO(LL,K,15) + 
     *   WQANCC1(IMWQZT(LL))*WQVO(LL,K,1)
     *   +WQANCD1(IMWQZT(LL))*WQVO(LL,K,2)
     *   +WQANCG1(IMWQZT(LL))*WQVO(LL,K,3)
         TKN=TNWQ-WQVO(LL,K,15)                                         ! TKN=TN-NOX
!
! BOD5 based on HydroQual, LISS report (1991), P4-39, Eq. (4-24)
!
         BOD5 = 2.67 * ( WQVO(LL,K,5)*(1.0-exp(-5.0*WQKLC))
     +           + WQVO(LL,K,6)*(1.0-EXP(-5.0*WQKDC))
     +           + WQVO(LL,K,18)*(1.0-EXP(-5.0*WQKCD))
     +           + WQVO(LL,K,1)*(1.0-EXP(-5.0*WQBMRc(1)))
     +           + WQVO(LL,K,2)*(1.0-EXP(-5.0*WQBMRD(1)))
     +           + WQVO(LL,K,3)*(1.0-EXP(-5.0*WQBMRG(1))) )
     +           + 4.57 * WQVO(LL,K,14)*(1.0-EXP(-5.0*WQNITM))
!
          DO KK=1,19
          TMPOUT(KK)=WQVO(LL,K,KK)
          ENDDO
          TMPOUT(20)=CHLWQ
          TMPOUT(21)=TOCWQ
          TMPOUT(22)=TNWQ
          TMPOUT(23)=TPWQ          
          TMPOUT(24)=PO4DWQ
          TMPOUT(25)=TKN
          TMPOUT(26)=SEDTWQ(LL,K) !!BOD5
          TMPOUT(27)=WQketot(LL,K)   
          TMPOUT(28)=TEMWQ(LL,K)
          TMPOUT(29)=wqvo(LL,K,22)  
          TMPOUT(30)=SALWQ(LL,K)          
          write(551)TMPOUT
         
         END DO                                      
      END DO     
!	                                    
! reset the concentration for output in other subroutine    
C
   71 FORMAT(3I5,F11.5, 1p, 23E11.3)
C
C  Output min and max
C
      if(isWQMIN==0) return
           
       DO M=1,IWQTS                           ! Require active save binary file (DUMP3)  
       DO K=1,KC                              
        LL=LWQTS(M)      
         BOD5 = 2.67 * ( WQVmax(LL,K,5)*(1.0-exp(-5.0*WQKLC))
     +           + WQVmax(LL,K,6)*(1.0-EXP(-5.0*WQKDC))
     +           + WQVmax(LL,K,18)*(1.0-EXP(-5.0*WQKCD))
     +           + WQVmax(LL,K,1)*(1.0-EXP(-5.0*WQBMRc(1)))
     +           + WQVmax(LL,K,2)*(1.0-EXP(-5.0*WQBMRD(1)))
     +           + WQVmax(LL,K,3)*(1.0-EXP(-5.0*WQBMRG(1))) )
     +           + 4.57 * WQVmax(LL,K,14)*(1.0-EXP(-5.0*WQNITM))
     
       
          APCWQ = 1.0 / (WQCP1PRM1(IMWQZT(LL))
     &      + WQCP2PRM*EXP(-WQCP3PRM*XPO4DWQ))
     
          TPWQ1 = APCWQ*(WQVmax(LL,K,1)+WQVmax(LL,K,2)+WQVmax(LL,K,3))    ! included algae's
   
          TNWQ1 = WQANCC1(IMWQZT(LL))*WQVmax(LL,K,1)
     *            + WQANCD1(IMWQZT(LL))*WQVmax(LL,K,2)
     *            + WQANCG1(IMWQZT(LL))*WQVmax(LL,K,3)     
     
          TPWQ2 = APCWQ*(WQVmin(LL,K,1)+WQVmin(LL,K,2)+WQVmin(LL,K,3))    ! included algae's
   
          TNWQ2 = WQANCC1(IMWQZT(LL))*WQVmin(LL,K,1)
     *            + WQANCD1(IMWQZT(LL))*WQVmin(LL,K,2)
     *            + WQANCG1(IMWQZT(LL))*WQVmin(LL,K,3)  
  
          AAAA=WQ3DA(LL,K,1)/WQCHLC(IWQZMAP(LL,KC))
     &          +WQ3DA(LL,K,2)/WQCHLD(IWQZMAP(LL,KC))
     &          +WQ3DA(LL,K,3)/WQCHLG(IWQZMAP(LL,KC))
     
        TPWQ3 = APCWQ*AAAA    ! included algae's
   
        TNWQ3 = WQANCC1(IMWQZT(LL))*WQ3DA(LL,K,1)/WQCHLC(IWQZMAP(LL,KC))
     *        + WQANCD1(IMWQZT(LL))*WQ3DA(LL,K,2)/WQCHLD(IWQZMAP(LL,KC))
     *        + WQANCG1(IMWQZT(LL))*WQ3DA(LL,K,3)/WQCHLG(IWQZMAP(LL,KC))
        
     
          BOD51 = 2.67 * ( WQVmin(LL,K,5)*(1.0-exp(-5.0*WQKLC))
     +           + WQVmin(LL,K,6)*(1.0-EXP(-5.0*WQKDC))
     +           + WQVmin(LL,K,18)*(1.0-EXP(-5.0*WQKCD))
     +           + WQVmin(LL,K,1)*(1.0-EXP(-5.0*WQBMRc(1)))
     +           + WQVmin(LL,K,2)*(1.0-EXP(-5.0*WQBMRD(1)))
     +           + WQVmin(LL,K,3)*(1.0-EXP(-5.0*WQBMRG(1))) )
     +           + 4.57 * WQVmin(LL,K,14)*(1.0-EXP(-5.0*WQNITM))
     
       
        IF(ITRYCHLA.GE.2) THEN             !JS  2018
        WQCHLD(IWQZMAP(LL,KC))=1/(KECST+90*exp(-1.19*WQketot(LL,KC)))
     &  *1000      
        WQCHLC(IWQZMAP(LL,KC))=1/(KECST+150*exp(-1.18*WQketot(LL,KC)))
     &  *1000  
        WQCHLG(IWQZMAP(LL,KC))=WQCHLC(IWQZMAP(LL,KC)) 
        ENDIF
        
c          CHLWQ1 = WQVmin(LL,K,1)*WQCHLC(IWQZMAP(LL,KC)) + 
c     *     WQVmin(LL,K,2)*WQCHLD(IWQZMAP(LL,KC))             ! total algae as chlorophyll, ug/l
c     *      + WQVmin(LL,K,3)*WQCHLG(IWQZMAP(LL,KC))  
c     
c           CHLWQ2 = WQ3DA(LL,K,1)*WQCHLC(IWQZMAP(LL,KC))
c     &           + WQ3DA(LL,K,2)*WQCHLD(IWQZMAP(LL,KC)) 
c     &           + WQ3DA(LL,K,3)*WQCHLG(IWQZMAP(LL,KC))       
c           CHLWQ = WQVmax(LL,K,1)*WQCHLC(IWQZMAP(LL,KC)) +
c     *        WQVmax(LL,K,2)*WQCHLD(IWQZMAP(LL,KC))             ! total algae as chlorophyll, ug/l
c     *      + WQVmax(LL,K,3)*WQCHLG(IWQZMAP(LL,KC))  
    
           CHLWQ1 =WQVmin(LL,K,1)+ 
     *             WQVmin(LL,K,2)+ WQVmin(LL,K,3)            ! total algae as chlorophyll, ug/l   
           CHLWQ2 =WQ3DA(LL,K,1) +
     &             WQ3DA(LL,K,2) +  WQ3DA(LL,K,3)          
           CHLWQ = WQVmax(LL,K,1)+
     *             WQVmax(LL,K,2)+ WQVmax(LL,K,3)           ! total algae as chlorophyll, ug/l
      
          BOD52 = 2.67 * ( WQ3DA(LL,K,5)*(1.0-exp(-5.0*WQKLC))
     +           + WQ3DA(LL,K,6)*(1.0-EXP(-5.0*WQKDC))
     +           + WQ3DA(LL,K,18)*(1.0-EXP(-5.0*WQKCD))
     +           + WQ3DA(LL,K,1)*(1.0-EXP(-5.0*WQBMRc(1)))
     +           + WQ3DA(LL,K,2)*(1.0-EXP(-5.0*WQBMRD(1)))
     +           + WQ3DA(LL,K,3)*(1.0-EXP(-5.0*WQBMRG(1))) )
     +           + 4.57 * WQ3DA(LL,K,14)*(1.0-EXP(-5.0*WQNITM))

          TMPOUT(1)=WQVmax(LL,K,1) !*WQCHLC(IWQZMAP(LL,KC))
          TMPOUT(2)=WQVmax(LL,K,2) !*WQCHLD(IWQZMAP(LL,KC)) 
          TMPOUT(3)=WQVmax(LL,K,3) !*WQCHLG(IWQZMAP(LL,KC))  
          DO KK=4,19  
          TMPOUT(KK)=WQVmax(LL,K,KK)
          ENDDO        
          TMPOUT(20)=CHLWQ
          TMPOUT(21)=WQTOTmax(LL,K,1)
          TMPOUT(22)=WQTOTmax(LL,K,2)
          TMPOUT(23)=WQTOTmax(LL,K,3)          
          TMPOUT(24)=WQTOTmax(LL,K,4)
          TMPOUT(25)=WQTOTmax(LL,K,5)
          TMPOUT(26)=SEDTWQ(LL,K) !BOD5
          TMPOUT(27)=WQketot(LL,K)   
          TMPOUT(28)=TEMWQ(LL,K)
          TMPOUT(29)=WQVmax(LL,K,22)  
          TMPOUT(30)=SALWQ(LL,K)           
          WRITE(553)TMPOUT
          
          TMPOUT(1)=WQVmin(LL,K,1) !*WQCHLC(IWQZMAP(LL,KC))
          TMPOUT(2)=WQVmin(LL,K,2) !*WQCHLD(IWQZMAP(LL,KC)) 
          TMPOUT(3)=WQVmin(LL,K,3) !*WQCHLG(IWQZMAP(LL,KC))  
          DO KK=4,19
          TMPOUT(KK)=WQVmin(LL,K,KK)
          ENDDO
          TMPOUT(20)=CHLWQ1
          TMPOUT(21)=WQTOTmin(LL,K,1)
          TMPOUT(22)=WQTOTmin(LL,K,2)
          TMPOUT(23)=WQTOTmin(LL,K,3)          
          TMPOUT(24)=WQTOTmin(LL,K,4)
          TMPOUT(25)=WQTOTmin(LL,K,5)
          TMPOUT(26)=SEDTWQ(LL,K) !BOD51
          TMPOUT(27)=WQVmin(LL,K,22)   
          TMPOUT(28)=TEMWQ(LL,K)  
          TMPOUT(29)=WQVmin(LL,K,22)  
          TMPOUT(30)=SALWQ(LL,K)          
          WRITE(552)TMPOUT
c        ENDIF
          TMPOUT(1)=WQ3DA(LL,K,1)!*WQCHLC(IWQZMAP(LL,KC))
          TMPOUT(2)=WQ3DA(LL,K,2)!*WQCHLD(IWQZMAP(LL,KC)) 
          TMPOUT(3)=WQ3DA(LL,K,3)!*WQCHLG(IWQZMAP(LL,KC))  
          DO KK=4,19
          TMPOUT(KK)=WQ3DA(LL,K,KK)
          ENDDO
          TMPOUT(20)=CHLWQ2 
          TMPOUT(21)=WQTOT(LL,K,1)
          TMPOUT(22)=WQTOT(LL,K,2)
          TMPOUT(23)=WQTOT(LL,K,3)          
          TMPOUT(24)=WQTOT(LL,K,4)
          TMPOUT(25)=WQTOT(LL,K,5)
          TMPOUT(26)=SEDTWQ(LL,K) !BOD52
          TMPOUT(27)=WQketot(LL,K)     
          TMPOUT(28)=TEMWQ(LL,K) 
          TMPOUT(29)=WQ3DA(LL,K,22)  
          TMPOUT(30)=SALWQ(LL,K)           
          WRITE(554)TMPOUT
                   
c        ENDIF               
      END DO                                      
      END DO     
!	                                    
! reset the concentration for output in other subroutine    
!
C
      RETURN
      END
