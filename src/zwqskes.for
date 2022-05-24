C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WQSKES(LA,KC,DT,N,TCON,TBEGIN,TIDALP,NTSPTC,IWQS,
     &  SECDLAST)
C
C**********************************************************************C
C
C  This is simplified version of subroutin to solve all kinetic Eqs. 
C  from K=KC (surface layer) to K=1 (bottom).
C: after computing new values, store WQVO+WQV into WQVO(L,K,NWQV).  
C
C
C  ORGINALLY CODED BY K.-Y. PARK
C  OPTIMIZED AND MODIFIED BY J.M. HAMRICK
C  LAST MODIFIED BY J.M. HAMRICK  7 APRIL 1997
C
C  MODIFIED ON 7/25/2008 BY AEE
C
C**********************************************************************C
C
      INCLUDE 'wq.par'  
      INCLUDE 'wqcom.cmn'
      REAL*8 SECDLAST 
      DIMENSION WQDOS(LCMWQ), WQI0BOT(LCMWQ),WQPM1(LCMWQ)

      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014  
      CNS1=2.718                                                         ! (3-1e)

      NS=1

	NDMWQ=1
	LDM=LA-1
C
C
C Initial solar radiation at top of surface layer  
C
      DO L=2,LA
        WQI0BOT(L)=WQI0 * pSHADE(L)                   ! Light (pSHADE is forced to  =1)
        xI00(L)= xI00(L)+WQI0 * pSHADE(L) *DTWQ       ! J.S. 8/16/08
        xHXY(L)=xHXY(L)+HPWQ(L)        
      END DO
C
C Start of DO ND=1,NDMWQ loop
C
                             
      LF=2      
      LL=LA        

C Start of DO K=KC,1,-1  loop

      DO K=KC,1,-1       ! END #200                         
C
C DZWQ=1/h, VOLWQ=1/VOL
C
       DO L=LF,LL
        TWQ(L)=TEMWQ(L,K)
        SWQ(L)=MAX(SALWQ(L,K), 0.0)
        DZWQ(L) = 1.0 / (DZCWQ(K)*HPWQ(L))
        VOLWQ(L) = DZWQ(L) / DXYPWQ(L)
        IMWQZT(L)=IWQZMAP(L,K)
       END DO
C
       DO L=LF,LL
        WQBGSET(L,1) = WQWSG(IMWQZT(L))*DZWQ(L)
        WQLPSET(L,1) = WQWSLP(IMWQZT(L))*DZWQ(L)
       END DO
C
       IF (K.NE.KC) THEN
C
        DO L=LF,LL
         IMWQZT1(L)=IWQZMAP(L,K+1)
        END DO
C
        DO L=LF,LL
         WQBGSET(L,2) = WQWSG(IMWQZT1(L))*DZWQ(L)
         WQLPSET(L,2) = WQWSLP(IMWQZT1(L))*DZWQ(L)
        END DO
C
       END IF
C
C find an index for look-up table for temperature dependency
C
       DO L=LF,LL
        IWQT(L) = 10.0*TWQ(L) + 51    !must be consistent with the look-up table:-5C-> 50C in WQINPUTS
        IF (IWQT(L).LT.1 .OR. IWQT(L).GT.NWQTD) THEN
          timtmp = (dt*float(n) + tcon*tbegin)/86400.0
          OPEN(1,FILE='error.log',ACCESS='APPEND',STATUS='UNKNOWN')
          write(1,911) timtmp, L, ILW(L), JLW(L), K, twq(L)
911       format(/,'ERROR!! invalid water temperature, sub WQSKE',/,
     +     'TIME, L, I, J, K, TWQ(L) = ', f10.5, 4i4, f10.4)
          IWQT(L)=NWQTD
          close(1)
          WRITE(6,600)ILW(L),JLW(L),K,TWQ(L)
        !  STOP 'ERROR!! invalid water temperature'
        END IF
       END DO
C
  600 FORMAT(' I,J,K,TEM = ',3I5,E13.4)
C
C Note: MRM 04/29/99  Added arrays to keep track of
C       nitrogen, phosphorus, light, and temperature limits
C       for algae growth for cyanobacteria, diatoms, greens,         ! Ji, useful!!, 9/18/02
C       and macroalgae.  These are the arrays:
C        xLimNx(L,K) = nitrogen    limitation for algae group x
C        xLimPx(L,K) = phosphorus  limitation for algae group x
C        xLimIx(L,K) = light       limitation for algae group x
C        xLimTx(L,K) = temperature limitation for algae group x
C
C1-3 algal growth: nutrient
C
       DO L=LF,LL                              ! END #100
        RNH4WQ = MAX (WQVO(L,K,14), 0.0)
        RNO3WQ = MAX (WQVO(L,K,15), 0.0)
        PO4DWQ = MAX (WQPO4D(L,K), 0.0)
        RNH4NO3 = RNH4WQ + RNO3WQ
        WQGNG = max(RNH4NO3 / (WQKHNG+RNH4NO3+ 1.E-18),AN_Lim)   
        WQGPG = max(PO4DWQ / (WQKHPG+PO4DWQ+ 1.E-18),P_Lim)
        xlimNg(L,K) = xlimNg(L,K) + WQGNG
        xlimPg(L,K) = xlimPg(L,K) + WQGPG
        WQF1NG = MIN(WQGNG, WQGPG)                       ! f(N) Eq. (3-1c)
        
C
        IF(IDNOTRVA.GT.0 .AND. K.EQ.1) THEN              ! macroalgae
          WQGNM = RNH4NO3 / (WQKHNM+RNH4NO3 + 1.E-18)    ! N
          WQGPM = PO4DWQ / (WQKHPM+PO4DWQ + 1.E-18)      ! P
          WQF1NM = MIN(WQGNM, WQGPM)                     ! N limit fun.
 !     WQF1NM = MAX(WQF1NM,WQFNCm)           ! 2010 modify nutrient limiting 
          xlimNm(L,K) = xlimNm(L,K) + WQGNM
          xlimPm(L,K) = xlimPm(L,K) + WQGPM
        END IF
C
!

!        IF(IDNOTRVA.GT.0) PO4DWQ = MAX (WQPO4D(L,K), 0.0)
 
C
C algal growth: light, WQHT(K)=REAL(KC-K)/REAL(KC)
C In C&C, F2Ic=F2Ic/FCYAN, factor to allow cyanobacteria mat formation
C
C MRM 05/12/1999 use Riley (1956) equation to compute light extinction
C     as a function of CHL conc. if WQKECHL is less than zero:
C
        WQCHL(L,K) =  WQV(L,K,3)*WQCHLG(I)
!        xmrm = WQKEMA
        IF (WQKECHL(IMWQZT(L)) .LT. 0.0) THEN
          xmrm = 0.054*WQCHL(L,K)**0.6667 + 0.0088*WQCHL(L,K)
        END IF
        WQKESS = WQKEB(IMWQZT(L))+WQKETSS(IMWQZT(L))*SEDTWQ(L,K) + xmrm
        WQKESS1 = WQKESS
        IF (K.NE.KC) THEN
          xmrm = WQKECHL(IMWQZT(L))*WQCHL(L,KC)
          IF (WQKECHL(IMWQZT(L)).LT. 0.0) THEN
            xmrm = 0.054*WQCHL(L,KC)**0.6667 + 0.0088*WQCHL(L,KC)
          END IF
         WQKESS1=WQKEB(IMWQZT(L))+WQKETSS(IMWQZT(L))*SEDTWQ(L,KC) + xmrm
        END IF

         WQAVGIO = WQCIA*WQI0 + WQCIB*WQI1 + WQCIC*WQI2       ! Average light
C
        IF (IWQSUN .EQ. 2) THEN                               !use hourly solar radiation from aser.inp
          WQAVGIO = WQCIA*WQI1 + WQCIB*WQI2 + WQCIC*WQI3
        END IF

        xKe(L)=xKe(L)+WQKESS1/Real(KC)   ! Record Ke  J.S. 8/16/08

        WQAVGIO = WQAVGIO * pSHADE(L)

        WQISG = MAX( WQAVGIO*EXP(-WQKESS1*WQDOPG), WQISMIN )
        WQTT1 = (CNS1 * WQFD * DZWQ(L)) / WQKESS                         ! (3-1e)
        WQFDI0 = - WQI0BOT(L) / (WQFD + 1.E-18)
C
        WQFDG = WQFDI0 / (WQISG + 1.E-18)
        WQHTT = WQHT(K) * HPWQ(L)
C
        WQTTB = EXP( -WQKESS * (WQHTT+1.0/DZWQ(L)) )                     ! (3-1f)
        WQTTT = EXP( -WQKESS * WQHTT )                                   ! (3-1g)
        WQF2IG = WQTT1 * (EXP(WQFDG*WQTTB) - EXP(WQFDG*WQTTT))

!
! J.S. Test new light calculation (12/18/04)
!
        WQIZ=WQI0*EXP(-WQKESS1 * (WQHTT+0.5/DZWQ(L)) )                      ! light at Z
!	  WQI_0=  WQI0*WQI0+ PARADJ*2.065*350                     ! Is 
	  
	  if(WQISMIN.GE.100.0) then
	  WQF2IG =WQIZ/sqrt( WQI0*WQI0+WQISMIN*WQISMIN)           ! 2010 J.S. Using WQISMIN read from Card 10 (set to 200)
        endif
        
        xlimIg(L,K) = xlimIg(L,K) + WQF2IG
C
C update solar radiation at bottom of this layer
C
        WQI0BOT(L)=WQI0BOT(L)*exp(-WQKESS*(1.0/DZWQ(L)))
C
        IF (IDNOTRVA.GT.0 .AND. K.EQ.1) THEN                    ! Macalgae
          
          WQFDI0 = - WQI0BOT(L) / (WQFD + 1.E-18)
          WQISM = MAX( WQAVGIO*EXP(-WQKESS1*WQDOPM), WQISMIN )
          WQFDM = WQFDI0 / (WQISM + 1.E-18)
          WQF2IM = WQTT1 * (EXP(WQFDM*WQTTB) - EXP(WQFDM*WQTTT))
          
          if(WQISMIN.GE.100.0) then                  !use new light function; con. nomorlized by 30
           if(WQVO(L,1,IDNOTRVA).lt.30)then
            WQ_t=1
           else
            WQ_t=WQVO(L,1,IDNOTRVA)/30.0
           endif
           
           if(WQKEMA(IMWQZT(L)).GT.0)then
             WQKESS1=WQKETSS(IMWQZT(L))*SEDTWQ(L,KC) 
     &       +abs(WQKEMA(IMWQZT(L)))*WQCHL(L,KC)
     &       +WQKEB(IMWQZT(L))*WQVO(L,1,IDNOTRVA)
           else
            WQ_t=WQVO(L,1,IDNOTRVA)
            WQKESS1=WQKEB(IMWQZT(L))+WQKETSS(IMWQZT(L))*SEDTWQ(L,KC)
     &      +abs(WQKEMA(IMWQZT(L)))*WQ_t
           endif
          
           WQIZ=WQI0*EXP(-WQKESS1 * (WQHTT) ) 
           WQF2IG =WQIZ/sqrt( WQI0*WQI0+WQISMIN*WQISMIN)      
           WQF2IM =WQF2IG *WQCIM     
                                                     ! WQF2IM:light limit,test new light for micalgae
          endif
          
       WQPM(L)= WQPMM(IMWQZT(L))*WQF1NM*WQF2IM*WQTDGM(IWQT(L),IMWQZT(L))       ! WQF1NM :nutrient limit;WQTDGM: temp. limit;
       WQPM1(L)= WQPMM(IMWQZT(L))*max(WQF1NM,WQFNCm)                ! Modify nutrient limist
     &        *WQF2IM*WQTDGM(IWQT(L),IMWQZT(L))
          xlimIm(L,K) = xlimIm(L,K) + WQF2IM
          xlimTm(L,K) = xlimTm(L,K) + WQTDGM(IWQT(L),IMWQZT(L))                    ! (3-1k)
        END IF
        xlimTg(L,K) = xlimTg(L,K) + WQTDGG(IWQT(L),IMWQZT(L))
C 

C
C: WQSTOX=WQSTOX**2
C
        WQPG(L)=WQPMG(IMWQZT(L))*WQF1NG*WQF2IG*WQTDGG(IWQT(L),IMWQZT(L))        !         green
        IF(iTMgop.EQ.1)
     &  WQPG(L) = WQPMG(IMWQZT(L))*WQF1NG*WQF2IG*
     &    max(WQTDGG(IWQT(L),IMWQZT(L)),WQTDGM(IWQT(L),IMWQZT(L))) 

c
c Using hourly solar radiation, shut down photosynthesis
c at night, i.e., when solar radiation is less than 0.001.
c
        if (iwqsun .eq. 2) then
          if (wqi0 .le. 0.100) then     
            wqpc(L) = 0.0
            wqpd(L) = 0.0
            wqpg(L) = 0.0
            wqpm(L) = 0.0
          end if
        end if
C
C algal basal metabolism & predation
C
        WQBMG(L) = WQBMRG(IMWQZT(L)) * WQTDRG(IWQT(L))       ! Eq. (3-1m)         
        WQPRG(L) = WQPRRG(IMWQZT(L)) * WQTDRG(IWQT(L))       ! Eq. (3-1n)
C
C  Macalgae
C
      IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQBMM(L) = WQBMRM(IMWQZT(L)) * WQTDRM(IWQT(L))
        WQPRM(L) = WQPRRM(IMWQZT(L)) * WQTDRM(IWQT(L))
      END IF

C
C4-6 organic carbon: WQAANOX=WQAANOX*WQKHORDO
C
        WQOBTOT = WQVO(L,K,3)        ! Total algae
        WQKLPC(L) = (WQKLC + WQKLCALG*WQOBTOT) * WQTDHDR(IWQT(L))      ! WQKLCALG is set to zero for simple model
        xmrm = 0.0
        IF (IDNOTRVA.GT.0 .AND. K.EQ.1) THEN
          xmrm = WQKDCALM * WQVO(L,K,IDNOTRVA)
        END IF
        WQKDOC = (WQKDC + WQKDCALG*WQOBTOT + xmrm) * WQTDMNL(IWQT(L))  ! DOC decay
        O2WQ = MAX(WQVO(L,K,19), 0.0)
        WQTT1 = WQKDOC / (WQKHORDO + O2WQ+ 1.E-18)
        WQKHR(L) = WQTT1 * O2WQ
        WQDENIT(L) = WQTT1 * WQAANOX * RNO3WQ/(WQKHDNN+RNO3WQ+ 1.E-18)
C
C 7-10 phosphorus
C
C    wqKHPC = KHPg = phosphorus half-saturation for algae greens algae (mg/L)

        WQAPC(L) = 1.0 / (WQCP1PRM1(IMWQZT(L))
     &     + WQCP2PRM*EXP(-WQCP3PRM*PO4DWQ)) ! Eq.(3-8e), wq3dwc.inp, C21

        WQKHP = WQKHPG                                               ! Eq.(3-8i)
        WQTT1 = WQKHP / (WQKHP+PO4DWQ+ 1.E-18) * WQOBTOT             !    (3-8f)
        WQKLPP(L) = (WQKLP + WQKLPALG*WQTT1) * WQTDHDR(IWQT(L))      !    (3-8g)
        WQKDOP(L) = (WQKDP + WQKDPALG*WQTT1) * WQTDMNL(IWQT(L))      !    (3-8h)
C
C10 PO4t
C
        IF (IWQSRP.EQ.1) THEN
          WQTTM = WQKPO4P*WQTAMP(L,K)
          WQH10(L) = - WQWSSET(L,1) * WQTTM / (1.0+WQTTM)
          IF (K.NE.KC) THEN
            WQTTM = WQKPO4P*WQTAMP(L,K+1)
            WQT10(L) = WQWSSET(L,2) * WQTTM / (1.0+WQTTM)
          END IF
         ELSE IF (IWQSRP.EQ.2) THEN                   ! SEDT
C   KPO4p = partition coefficient for sorbed/dissolved PO4
          WQTTS = WQKPO4P*SEDTWQ(L,K)
          WQH10(L) = - WQSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)    ! Eq. (3-8), settling term
          IF (K.NE.KC) THEN
            WQTTS = WQKPO4P*SEDTWQ(L,K)
            WQT10(L) = WQSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)
          END IF                                      ! SEDT
         ELSE
          WQH10(L) = 0.0
          WQT10(L) = 0.0
        END IF
C
C 11-15 nitrogen
C
        WQKHN = WQKHNG                                                  ! (3-13e)
        WQTT1 = WQKHN / (WQKHN+RNH4NO3+ 1.E-18) * WQOBTOT
        WQKLPN(L) = (WQKLN + WQKLNALG*WQTT1) * WQTDHDR(IWQT(L))         ! (3-13c)
        WQKDON(L) = (WQKDN + WQKDNALG*WQTT1) * WQTDMNL(IWQT(L))         ! (3-13d)
C
C14 NH4: WQFTNIT=WQNITM*WQFTNIT
C
        IF (RNH4NO3.LE.1.0E-12) THEN             ! =NH4+NOx
          WQPNG(L)=0.0
          WQPNM(L)=0.0
        ELSE
          WQTTG = RNH4WQ/(WQKHNG+RNO3WQ+ 1.E-18)
          WQTTM = RNH4WQ/(WQKHNM+RNO3WQ+ 1.E-18)
          WQPNG(L) = (RNO3WQ/(WQKHNG+RNH4WQ+ 1.E-18)                    ! green
     $              + WQKHNG/(RNH4NO3+ 1.E-18)) * WQTTG
          WQPNM(L) = (RNO3WQ/(WQKHNM+RNH4WQ+ 1.E-18)                    ! macroalgae
     $              + WQKHNM/(RNH4NO3+ 1.E-18)) * WQTTM

        END IF

          WQNIT(L) = O2WQ * RNH4WQ * WQTDNIT(IWQT(L)) /                 ! (3-13g)
     *      ( (WQKHNDO+O2WQ) * (WQKHNN+RNH4WQ) + 1.E-18)

C
C The following arrays were added to keep track of the various components
C of dissolved oxygen.  The instantaneous values for each component are
C summed in the arrays and then dumped to the WQDOCOMP.bin file at the
C same time interval as for the WQWCAVG.bin files (i.e., IWQTSDT intervals,
C usually daily averages).  The array descriptions are:
C
C  xDOsat(L,K) = D.O. saturation for cell L, layer K (mg/L)
C  xDOdef(L,K) = D.O. deficit for cell L, layer K (mg/L)
C  xDOpsl(L,K) = D.O. component for external loads (mg/L/day)
C  xDOsod(L,K) = D.O. component for sediment oxygen demand
C  xDOkar(L,K) = D.O. component for reaeration
C  xDOdoc(L,K) = D.O. component for Diss. Org. Carbon decay
C  xDOnit(L,K) = D.O. component for ammonia nitrification
C  xDOcod(L,K) = D.O. component for Chem. Oxy. Demand oxidation
C  xDOppB(L,K) = D.O. component for photosynthesis of total chlorophyll
C  xDOrrB(L,K) = D.O. component for respiration of total chlorophyll
C  xDOppM(L,K) = D.O. component for photosynthesis of macroalgae
C  xDOrrM(L,K) = D.O. component for respiration of macroalgae
C  xDOall(L,K) = sum of the above 10 D.O. components
C  xDOdz (L,K) = layer thickness (meters)
C  NLIM = counter for number of items summed in each array slot
C
C18-19 COD, O2: WQO18(L)=DTWQO2*WQO18, WQKRDOS(L)=-WQP19(L)*WQDOS
C
        WQO18(L)=- DTWQO2*WQKCOD(IWQT(L))*O2WQ/(WQKHCOD+O2WQ+ 1.E-18)    ! (3-16)
C
C Tt The following modification to the D.O. saturation calculation made
C Tt by J.M. Hamrick / M.R. Morton on 03/08/97.  See Chapra (1997) pg. 361-364.
C
        TVAL1=1./(TWQ(L)+273.15)
        TVAL2=TVAL1*TVAL1
        TVAL3=TVAL1*TVAL2
        TVAL4=TVAL2*TVAL2
        RLNSAT1=-139.3441+(1.575701E+5*TVAL1)-(6.642308E+7*TVAL2)
     $                   +(1.2438E+10*TVAL3)-(8.621949E+11*TVAL4)
        RLNSAT2=RLNSAT1-SWQ(L)*( 1.7674E-2-(1.0754E+1*TVAL1)
     $                           +(2.1407E+3*TVAL2) )
        WQDOS(L) = EXP(RLNSAT2)
        xDOsat(L,K) = xDOsat(L,K) + WQDOS(L)*DTWQ*DZCWQ(K)*HPWQ(L)

        IF (K.EQ.KC) THEN
c
c Do not allow wind speeds above 11 m/sec in the following equation:

        windrea = WINDSTWQ(L)
        if (WINDSTWQ(L) .GT. 11.0) windrea = 11.0
         WQWREA = 0.728*SQRT(windrea) + (0.0372*windrea - 0.317)*windrea
c
c MRM 04/29/1999  User specifies constant reaeration WQKRO:
c
          IF (IWQKA .EQ. 0) THEN           
            WQVREA = WQKRO
            WQWREA = 0.0
          END IF
c
c MRM 04/12/1999  Constant reaeration due to water velocity,
c                 wind velocity computed above:
c
          IF (IWQKA .EQ. 1) THEN
            WQVREA = WQKRO
          END IF
c
c MRM 03/06/1999  O'Connor-Dobbins (1958) equation for reaeration is:
c    WQKRO = 3.933 typically
c
          IF (IWQKA .EQ. 2) THEN
            umrm = UWQS(L)                               !max(U(L,K), U(L+1,K))
            vmrm = VWQS(L)                               !max(V(L,K), V(LNC(L),K))
            xmrm = SQRT(umrm*umrm + vmrm*vmrm)           !SQRT(U(L,K)*U(L,K) + V(L,K)*V(L,K))
            WQVREA = WQKRO * xmrm**0.5 / HPWQ(L)**1.5
          END IF
c
c MRM 04/12/1999  Owens and Gibbs (1964) reaeration equation:
c    WQKRO = 5.32 typically
c
          IF (IWQKA .EQ. 3) THEN
            umrm = UWQS(L)                               ! max(U(L,K), U(L+1,K))
            vmrm = VWQS(L)                               ! max(V(L,K), V(LNC(L),K))
            xmrm = SQRT(umrm*umrm + vmrm*vmrm)
            WQVREA = WQKRO * xmrm**0.67 / HPWQ(L)**1.85
          END IF
c
c now combine reaeration due to water velocity and wind stress:
c
          WQVREA = WQVREA * REAC(IWQZMAP(L,K))
          WQWREA = WQWREA * REAC(IWQZMAP(L,K))
          WQP19(L) = - (WQVREA + WQWREA) * DZWQ(L)* WQTDKR(IWQT(L))
          WQKRDOS(L) = - WQP19(L)*WQDOS(L)
          LL0=LIJW(46,103)
          IF(L.EQ.LL0.AND.MOD(N,500).EQ.0)write(*,'(8F8.3)')WQVREA,
     &  WQWREA,DZWQ(L),WQTDKR(IWQT(L)),WQKRDOS(L),WQP19(L),WQDOS(L),
     &  WQVO(L,K,19)
        ELSE
          WQP19(L) = 0.0
        END IF
C
       END DO                ! #100
C
C trapezoidal solution of kinetic eqs: After computing new values, store
C WQVO+WQV into WQVO(L,K,NWQV)

C J.S
C
       DO L=LF,LL
C
C3 Bg: WQT3D=WQBGSET(L,2)
C
        agsett=0.05
        if(IWQNC.eq.0)then                !2015 HAB add capacity and aggregation of settling
        a1 = (WQVO(L,K,14)+WQVO(L,K,15))
        a2 = min(a1/(WQPG(L)*WQV(L,K,3)*WQANCG+1.0e-8),1.0)
        else
        a2=1
        endif
        WQA3G = (WQPG(L)*a2 - WQBMG(L) - WQPRG(L) - WQBGSET(L,1)
     &       -agsett*WQVO(L,K,3))*DTWQO2
        WQKK(L) = 1.0 / (1.0 - WQA3G)
        WQR3G = (WQWDSL(L,K,3) + WQWPSL(L,K,3)) * VOLWQ(L)
        WQRR(L) = WQVO(L,K,3) + DTWQ*WQR3G + WQA3G*WQVO(L,K,3)
       END DO
C
       IF (K.NE.KC) THEN
         DO L=LF,LL
          WQRR(L) = WQRR(L) + DTWQO2*WQBGSET(L,2)*WQVO(L,K+1,3)
         END DO
       END IF
C
       DO L=LF,LL
       WQV(L,K,3)=SCBWQ(L)*( WQRR(L)*WQKK(L) )+(1.-SCBWQ(L))*WQVO(L,K,3)
       WQ_SET_G(L,K)=WQV(L,K,3)*WQBGSET(L,1)                    ! J.S. save settling algae
       WQVO(L,K,3) = WQVO(L,K,3)+WQV(L,K,3)
       END DO
!
C
C Modified by J.S.(5/5/98) for Macalgae
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        DO L=LF,LL
        if(IWQNC.eq.0)then
        a1 = (WQVO(L,K,14)+WQVO(L,K,15))-WQPG(L)*WQV(L,K,3)*WQANCG   ! estimate available N; alg uptake first, then macalage 
        a2 = min(a1/(WQPM1(L)*WQV(L,K,IDNOTRVA)*WQANCM+1.0e-8),1.0)
        else
        a2=1.0
        endif
        WQA1C = (WQPM1(L)*a2 - WQBMM(L) - WQPRM(L)-WQWSM*DZWQ(L))*DTWQO2               !Macalgae growth Growth is modified 1/25/2010
        WQVA1C = 1.0 / (1.0 - WQA1C)
        WQV(L,K,IDNOTRVA)=(WQVO(L,K,IDNOTRVA)+WQA1C*WQVO(L,K,IDNOTRVA))
     *                   *WQVA1C*SMAC(L)
        WQV(L,K,IDNOTRVA) = max(WQV(L,K,IDNOTRVA),WQMCMIN)*SMAC(L)
        WQVO(L,K,IDNOTRVA) = WQVO(L,K,IDNOTRVA)+WQV(L,K,IDNOTRVA)
        END DO
       END IF

!
! Record algae con. as Chl-a at n+1 and Production, J.S. 11/25/2014
!
       DO L=LF,LL
        xAlg(L)=xAlg(L)+0.5*( WQVO(L,K,3))/(Real(KC)) 
        xPro(L)=xPro(L)+ 
     &    WQV(L,K,3)*( (1+ DTWQO2*WQPG(L))/(1- DTWQO2*WQPG(L))-1 )
     &  *(DZCWQ(K)*HPWQ(L))
 
       ENDDO
!

C
C5 LPOC: WQT5=WQLPSET(L,2)
C 
      if(IWQS.EQ.1) then   ! J.S.  Simplified model I
C
       DO L=LF,LL
        WQC5 = - (WQKLPC(L)+WQLPSET(L,1))
        WQKK(L) = 1.0 / (1.0 - DTWQO2*WQC5)
        WQA5 = WQFCLP * WQPRG(L)*WQVO(L,K,3)

        IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQA5 =WQA5 + WQFCLPM * WQPRM(L)*WQVO(L,K,IDNOTRVA)  ! add macalgal source 
        END IF

        WQR5 = (WQWDSL(L,K,5) + WQWPSL(L,K,5)) * VOLWQ(L)
        WQRR(L) = WQVO(L,K,5) + DTWQ*WQR5 + DTWQO2*( WQA5
     *    + WQC5*WQVO(L,K,5) )
       END DO
C
       IF (K.NE.KC) THEN
         DO L=LF,LL
          WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,5)
         END DO
       END IF
C
       DO L=LF,LL
       WQV(L,K,5)=SCBWQ(L)*( WQRR(L)*WQKK(L) )+(1.-SCBWQ(L))*WQVO(L,K,5)
       WQ_SET_C(L,K)=WQV(L,K,5)*WQLPSET(L,2)                              ! J.S. save LPC
       WQVO(L,K,5) = WQVO(L,K,5)+WQV(L,K,5)
       ENDDO
C
C6 DOC: CFCDxWQ=1-WQFCDx, WQB6=WQKRPC(L),WQC6=WQKLPC(L)
C
       DO L=LF,LL
        WQD6 = - (WQKHR(L)+WQDENIT(L))
        WQKK(L) = 1.0 / (1.0 - DTWQO2*WQD6)
        O2WQ = MAX(WQVO(L,K,19), 0.0)
        WQA6G=(WQFCDG+CFCDGWQ*WQKHRG/(WQKHRG+O2WQ+ 1.E-18))*WQBMG(L)
        WQA6 = ( WQA6G + WQFCDP*WQPRG(L) )*WQVO(L,K,3)
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
       WQA6M=(WQFCDM+(1-WQFCDM)*WQKHRM/(WQKHRM+O2WQ+ 1.E-18))*WQBMM(L)
       WQA6 =WQA6+ (WQA6M+ WQFCDPM*WQPRM(L))*WQVO(L,K,IDNOTRVA)
       END IF
C
        WQR6 = (WQWDSL(L,K,6) + WQWPSL(L,K,6)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
c        if (k.eq.kc) then
c          WQR6 = WQR6 + WQATML(L,kc,6) * VOLWQ(L)
c        end if
        WQRR(L) = WQVO(L,K,6) + DTWQ*WQR6 + DTWQO2*( WQA6 + 
     *  WQKLPC(L)*WQVO(L,K,5) + WQD6*WQVO(L,K,6) )
       ENDDO

       DO L=LF,LL
        WQV(L,K,6)=SCBWQ(L)*( WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,6)
        WQVO(L,K,6) = WQVO(L,K,6)+WQV(L,K,6)
       END DO
C
      else
!
! Simplified model II. J.S.
!
C----------------------------------------------------------------------------
C6 DOC: CFCDxWQ=1-WQFCDx, WQB6=WQKRPC(L),WQC6=WQKLPC(L)
C
       DO L=LF,LL
        WQD6 = - (WQKHR(L)+WQDENIT(L)+WQLPSET(L,1))  ! add seltting WQLPSET
        WQKK(L) = 1.0 / (1.0 - DTWQO2*WQD6)
        O2WQ = MAX(WQVO(L,K,19), 0.0)
        WQA6G=(WQFCDG+CFCDGWQ*WQKHRG/(WQKHRG+O2WQ+ 1.E-18))*WQBMG(L)
        WQA6 = ( WQA6G + WQPRG(L) )*WQVO(L,K,3)     ! due predation. no partiction for simple model II
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
       WQA6M=(WQFCDM+(1-WQFCDM)*WQKHRM/(WQKHRM+O2WQ+ 1.E-18))*WQBMM(L)
       WQA6 =WQA6+ (WQA6M+ WQPRM(L))*WQVO(L,K,IDNOTRVA)
       END IF
C
        WQR6 = (WQWDSL(L,K,6) + WQWPSL(L,K,6)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
c        if (k.eq.kc) then
c          WQR6 = WQR6 + WQATML(L,kc,6) * VOLWQ(L)
c        end if
        WQRR(L) = WQVO(L,K,6) + DTWQ*WQR6 + DTWQO2*( WQA6 + 
     *  WQD6*WQVO(L,K,6) )
       ENDDO
C       
C Add Seltting
       IF (K.NE.KC) THEN
         DO L=LF,LL
          WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,6)
         END DO
       END IF
       
       DO L=LF,LL
        WQV(L,K,6)=SCBWQ(L)*( WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,6)
        WQVO(L,K,6) = WQVO(L,K,6)+WQV(L,K,6)
       END DO
C            
      endif
!
      if(IWQS.EQ.1) then   ! J.S.  Simplified model I
C
C8 LPOP: WQT8=WQT5=WQLPSET(L,2)
C 
      DO L=LF,LL
        WQF8 = - (WQKLPP(L)+WQLPSET(L,1))
        WQKK(L) = 1.0 / (1.0 - DTWQO2*WQF8)
        WQA8G = (WQFPLG*WQBMG(L) + WQFPLP*WQPRG(L)) * WQVO(L,K,3)
        WQA8 = (WQA8G) * WQAPC(L)
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
       WQA8 = WQA8 +     (WQFPLM*WQBMM(L) + WQFPLPM*WQPRM(L))
     *         * WQVO(L,K,IDNOTRVA)* WQAPC(L)*WQAPCM
       END IF
C J.S.
        WQR8 = (WQWDSL(L,K,8) + WQWPSL(L,K,8)) * VOLWQ(L)
        WQRR(L) = WQVO(L,K,8) + DTWQ*WQR8 + DTWQO2*( WQA8
     *     + WQF8*WQVO(L,K,8) )
       END DO
C
       IF (K.NE.KC) THEN
         DO L=LF,LL
          WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,8)
         END DO
       END IF
C
       DO L=LF,LL
        WQV(L,K,8)=SCBWQ(L)*( WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,8)
        WQVO(L,K,8) = WQVO(L,K,8)+WQV(L,K,8)
       END DO
C
C9 DOP: WQE9=WQKRPP(L),WQF9=WQKLPP(L),WQG9=-WQKDOP(L)
C
       DO L=LF,LL
        WQR9=-WQKDOP(L)
        WQKK(L) = 1.0 / (1.0 - DTWQO2*WQR9)
        WQA9G = (WQFPDG*WQBMG(L) + WQFPDP*WQPRG(L)) * WQVO(L,K,3)
        WQA9 = (WQA9G) * WQAPC(L)
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
       WQA9 = WQA9 + (WQFPDM*WQBMM(L) + WQFPDPM*WQPRM(L))
     *       * WQVO(L,K,IDNOTRVA) * WQAPC(L)*WQAPCM
       END IF
C 
        WQR9 = (WQWDSL(L,K,9) + WQWPSL(L,K,9)) * VOLWQ(L)
        WQRR(L) = WQVO(L,K,9) + DTWQ*WQR9 + DTWQO2*(WQA9 
     *   + WQKLPP(L)*WQVO(L,K,8) - WQKDOP(L)*WQVO(L,K,9))
        WQV(L,K,9)=SCBWQ(L)*( WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,9)
        WQVO(L,K,9) = WQVO(L,K,9)+WQV(L,K,9)
C
C10 PO4t: WQG10=WQKDOP(L),WQA10=WQKK(L)
C
        WQA10G = (WQFPIG*WQBMG(L)+WQFPIP*WQPRG(L)-WQPG(L))*WQVO(L,K,3)
        WQKK(L) = (WQA10G) * WQAPC(L)
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQKK(L) =WQKK(L)+(WQFPIM*WQBMM(L)+WQFPIP*WQPRM(L)-WQPM(L))
     *           *WQVO(L,K,IDNOTRVA) * WQAPC(L)*WQAPCM
       END IF
C 
        WQRR(L) = (WQWDSL(L,K,10)+WQWPSL(L,K,10)) * VOLWQ(L)

       END DO
C
       IF (K.EQ.1) THEN
         DO L=LF,LL
          WQRR(L) = WQRR(L) + WQBFPO4D(L)*DZWQ(L)
         END DO
       END IF
C
       DO L=LF,LL
        WQRR(L) = WQVO(L,K,10) + DTWQ*WQRR(L) + DTWQO2*( WQKK(L)
     *    + WQKDOP(L)*WQVO(L,K,9) + WQH10(L)*WQVO(L,K,10) )
       END DO
C
       IF (K.NE.KC) THEN
         DO L=LF,LL
          WQRR(L) = WQRR(L) + DTWQO2*WQT10(L)*WQVO(L,K+1,10)
         END DO
       END IF
C
       DO L=LF,LL
       WQKK(L) = 1.0 / (1.0 - DTWQO2*WQH10(L))
       WQV(L,K,10)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,10)
       WQVO(L,K,10) = WQVO(L,K,10)+WQV(L,K,10)
       END DO
C
       else   ! simplified model 2 TOP and PO4  J.S.
C-----------------------------------------------------------------
C9 DOP: WQE9=WQKRPP(L),WQF9=WQKLPP(L),WQG9=-WQKDOP(L)
C
       DO L=LF,LL
        WQR9=-WQKDOP(L)-WQLPSET(L,1)      ! 
        WQKK(L) = 1.0 / (1.0 - DTWQO2*WQR9)
        WQA9G = ( (WQFPDG+WQFPLG)*WQBMG(L) + 
     *     (WQFPDP+WQFPLP)*WQPRG(L)) * WQVO(L,K,3)
        WQA9 = (WQA9G) * WQAPC(L)
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
       WQA9 = WQA9 + ( (WQFPDM+ WQFPLM )*WQBMM(L) + 
     *  (WQFPDPM+WQFPLPM)*WQPRM(L))
     *       * WQVO(L,K,IDNOTRVA) * WQAPC(L)*WQAPCM
       END IF

        WQR9 = (WQWDSL(L,K,9) + WQWPSL(L,K,9)) * VOLWQ(L)
        WQRR(L) = WQVO(L,K,9) + DTWQ*WQR9 + DTWQO2*(WQA9 
     *    - WQKDOP(L)*WQVO(L,K,9)-WQLPSET(L,1)*WQVO(L,K,9))
     
      ENDDO
      
       IF (K.NE.KC) THEN
         DO L=LF,LL
          WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,9)
         END DO
       END IF
       
       DO L=LF,LL    
        WQV(L,K,9)=SCBWQ(L)*( WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,9)
        WQVO(L,K,9) = WQVO(L,K,9)+WQV(L,K,9)
       ENDDO

C
C10 PO4t: WQG10=WQKDOP(L),WQA10=WQKK(L)
C 
       DO L=LF,LL 
       
        WQA10G = (WQFPIG*WQBMG(L)+WQFPIP*WQPRG(L)-WQPG(L))*WQVO(L,K,3)
        WQKK(L) = (WQA10G) * WQAPC(L)
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQKK(L) =WQKK(L)+(WQFPIM*WQBMM(L)+WQFPIP*WQPRM(L)-WQPM(L))
     *           *WQVO(L,K,IDNOTRVA) * WQAPC(L)*WQAPCM
       END IF
C 
        WQRR(L) = (WQWDSL(L,K,10)+WQWPSL(L,K,10)) * VOLWQ(L)

       END DO
C
       IF (K.NE.KC) THEN
         DO L=LF,LL
          WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,9)
         END DO
       END IF
C
       IF (K.EQ.1) THEN
         DO L=LF,LL
          WQRR(L) = WQRR(L) + WQBFPO4D(L)*DZWQ(L)
         END DO
       END IF
C
       DO L=LF,LL
        WQRR(L) = WQVO(L,K,10) + DTWQ*WQRR(L) + DTWQO2*( WQKK(L)
     *    + WQKDOP(L)*WQVO(L,K,9) + WQH10(L)*WQVO(L,K,10) )
       END DO
C
       IF (K.NE.KC) THEN
         DO L=LF,LL
          WQRR(L) = WQRR(L) + DTWQO2*WQT10(L)*WQVO(L,K+1,10)
         END DO
       END IF
C
       DO L=LF,LL
       WQKK(L) = 1.0 / (1.0 - DTWQO2*WQH10(L))
       WQV(L,K,10)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,10)
       WQVO(L,K,10) = WQVO(L,K,10)+WQV(L,K,10)
       END DO
C      
       
       endif
!
C12 LPON: WQT12=WQT5=WQLPSET(L,2)
C
! J.S
      if(IWQS.EQ.1) then  ! Simplified version I
!
       DO L=LF,LL
         WQJ12 = - (WQKLPN(L)+WQLPSET(L,1))
         WQKK(L) = 1.0 / (1.0 - DTWQO2*WQJ12)
         WQA12G = (WQFNLG*WQBMG(L)+WQFNLP*WQPRG(L))*WQANCG*WQVO(L,K,3)
         WQA12 = WQA12G
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
       WQA12 =WQA12 +(WQFNLM*WQBMM(L)+WQFNLPM*WQPRM(L))
     *      *WQANCM*WQVO(L,K,IDNOTRVA)
       END IF
C J.S.
         WQR12 = (WQWDSL(L,K,12)+WQWPSL(L,K,12)) * VOLWQ(L)
         WQRR(L) = WQVO(L,K,12) + DTWQ*WQR12 + DTWQO2*( WQA12
     *     + WQJ12*WQVO(L,K,12) )
       END DO
C
       IF (K.NE.KC) THEN
         DO L=LF,LL
           WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,12)
         END DO
       END IF
C
       DO L=LF,LL
       WQV(L,K,12)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,12)
       WQVO(L,K,12) = WQVO(L,K,12)+WQV(L,K,12)
	 ENDDO
C
C13 DON: WQI13=WQKRPN(L),WQJ13=WQKLPN(L),WQK13=-WQKDON(L)
C
       DO L=LF,LL
        WQKK(L) = 1.0 / (1.0 + DTWQO2*WQKDON(L))
        WQA13G = (WQFNDG*WQBMG(L)+WQFNDP*WQPRG(L)) *WQANCG*WQVO(L,K,3)
        WQA13 = WQA13G
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQA13 =WQA13 + (WQFNDM*WQBMM(L)+WQFNDPM*WQPRM(L))
     *         *WQANCM*WQVO(L,K,IDNOTRVA)
       END IF
C J.S.
        WQR13 = (WQWDSL(L,K,13) + WQWPSL(L,K,13)) * VOLWQ(L)
        WQRR(L) = WQVO(L,K,13) + DTWQ*WQR13 + DTWQO2*( WQA13
     *     + WQKLPN(L)*WQVO(L,K,12)
     *     - WQKDON(L)*WQVO(L,K,13) )
       WQV(L,K,13)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,13)
       WQVO(L,K,13) = WQVO(L,K,13)+WQV(L,K,13)
C
C14 NH4: WQL14=-WQNIT(L),WQK14=WQKDON(L),WQR14=WQRR(L)
C

        WQRR(L) = (WQWDSL(L,K,14)+WQWPSL(L,K,14)) * VOLWQ(L)

       END DO
C
C J.S. 5/20/03  Test nutrient input due to sediment erosion
c
       IF (K.EQ.1) THEN
         DO L=LF,LL
           WQRR(L) = WQRR(L) + WQBFNH4(L)*DZWQ(L)   
!     &    +SMTHKN(2)**(TWQ(L)-20.) *SEDFF
         END DO
       END IF
C
       DO L=LF,LL
        WQKK(L) = 1.0 / (1.0 + DTWQO2*WQNIT(L))
        WQA14C = 0
        WQA14D = 0
        WQA14G = WQFNIG*WQBMG(L) + WQFNIP*WQPRG(L) - WQPNG(L)*WQPG(L)
        WQA14 = WQA14C*WQANCC*WQVO(L,K,1)
     *     + WQA14D*WQANCD*WQVO(L,K,2) + WQA14G*WQANCG*WQVO(L,K,3)
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
              WQA14 = WQA14 + (WQFNIM*WQBMM(L)+WQFNIPM*WQPRM(L)
     *      - WQPNM(L)*WQPM(L))*WQANCM*WQVO(L,K,IDNOTRVA)
       END IF
C J.S.
        WQRR(L) = WQVO(L,K,14) + DTWQ*WQRR(L) + DTWQO2*( WQA14
     *     + WQKDON(L)*WQVO(L,K,13) - WQNIT(L)*WQVO(L,K,14) )
       WQV(L,K,14)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,14)
       WQVO(L,K,14) = WQVO(L,K,14)+WQV(L,K,14)
C
C15 NO3: WQKK(L)=1,WQR15=WQRR(L),WQD15=-WQANDC*WQDENIT(L),WQL15=WQNIT(L)
C
         WQRR(L) = (WQWDSL(L,K,15)+WQWPSL(L,K,15)) * VOLWQ(L)
       END DO
C
       IF (K.EQ.1) THEN
         DO L=LF,LL
           WQRR(L) = WQRR(L) + WQBFNO3(L)*DZWQ(L)
!     &    +SMTHKN(2)**(TWQ(L)-20.) *SEDFF                 ! J.S
         END DO
       END IF
C
       DO L=LF,LL
         WQA15C = 0.0
         WQA15D = 0.0
         WQA15G = (WQPNG(L)-1.0)*WQPG(L) * WQANCG * WQVO(L,K,3)
         WQA15 = WQA15C+WQA15D+WQA15G
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQA15 =WQA15 + (WQPNM(L)-1.0)*WQPM(L)*WQANCM
     *         *WQVO(L,K,IDNOTRVA)
       END IF
C J.S.
        WQV(L,K,15)=SCBWQ(L)*( WQVO(L,K,15) + DTWQ*WQRR(L)
     *     + DTWQO2*( WQA15
     *     - WQANDC*WQDENIT(L)*WQVO(L,K,6) + WQNIT(L)*WQVO(L,K,14) ))
     *     +(1.-SCBWQ(L))*WQVO(L,K,15)
         WQVO(L,K,15) = WQVO(L,K,15)+WQV(L,K,15)
       END DO
!
! Simplified version II
!------------------------------------------------------
      else 
C
C13 DON: WQI13=WQKRPN(L),WQJ13=WQKLPN(L),WQK13=-WQKDON(L)
C
       DO L=LF,LL
        WQKK(L) = 1.0 / (1.0 + DTWQO2*(WQKDON(L)+WQLPSET(L,1)))
        WQA13G = ((WQFNDG+WQFNLG) *WQBMG(L)+
     *        (WQFNDP+WQFNLP )*WQPRG(L)) *WQANCG*WQVO(L,K,3)
        WQA13 = WQA13G
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQA13 =WQA13 + ( (WQFNDM+WQFNLM )*WQBMM(L)+
     *  (WQFNDPM+ WQFNLPM)*WQPRM(L))
     *         *WQANCM*WQVO(L,K,IDNOTRVA)
       END IF
C J.S.
        WQR13 = (WQWDSL(L,K,13) + WQWPSL(L,K,13)) * VOLWQ(L)
        WQRR(L) = WQVO(L,K,13) + DTWQ*WQR13 + DTWQO2*( WQA13
     *     - (WQKDON(L)+WQLPSET(L,1))*WQVO(L,K,13) )
     
       ENDDO
        IF (K.NE.KC) THEN         ! add settling J.S
         DO L=LF,LL
           WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,13)
         END DO
       END IF      
       DO L=LF,LL       
       WQV(L,K,13)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,13)
       WQVO(L,K,13) = WQVO(L,K,13)+WQV(L,K,13)
       ENDDO
C
C14 NH4: WQL14=-WQNIT(L),WQK14=WQKDON(L),WQR14=WQRR(L)
C

       DO L=LF,LL  
        WQRR(L) = (WQWDSL(L,K,14)+WQWPSL(L,K,14)) * VOLWQ(L)
       END DO
       
       IF (K.EQ.1) THEN
         DO L=LF,LL
           WQRR(L) = WQRR(L) + WQBFNH4(L)*DZWQ(L)   
         END DO
       END IF
C
       DO L=LF,LL
        WQKK(L) = 1.0 / (1.0 + DTWQO2*WQNIT(L))
        WQA14C = 0
        WQA14D = 0
        WQA14G = WQFNIG*WQBMG(L) + WQFNIP*WQPRG(L) - WQPNG(L)*WQPG(L)
        WQA14 = WQA14C*WQANCC*WQVO(L,K,1)
     *     + WQA14D*WQANCD*WQVO(L,K,2) + WQA14G*WQANCG*WQVO(L,K,3)
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
              WQA14 = WQA14 + (WQFNIM*WQBMM(L)+WQFNIPM*WQPRM(L)
     *      - WQPNM(L)*WQPM(L))*WQANCM*WQVO(L,K,IDNOTRVA)
       END IF
C J.S.
        WQRR(L) = WQVO(L,K,14) + DTWQ*WQRR(L) + DTWQO2*( WQA14
     *     + WQKDON(L)*WQVO(L,K,13) - WQNIT(L)*WQVO(L,K,14) )
       WQV(L,K,14)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,14)
       WQVO(L,K,14) = WQVO(L,K,14)+WQV(L,K,14)
C
C15 NO3: WQKK(L)=1,WQR15=WQRR(L),WQD15=-WQANDC*WQDENIT(L),WQL15=WQNIT(L)
C
         WQRR(L) = (WQWDSL(L,K,15)+WQWPSL(L,K,15)) * VOLWQ(L)
       END DO
C
       IF (K.EQ.1) THEN
         DO L=LF,LL
           WQRR(L) = WQRR(L) + WQBFNO3(L)*DZWQ(L)
         END DO
       END IF
C
       DO L=LF,LL
         WQA15C = 0.0
         WQA15D = 0.0
         WQA15G = (WQPNG(L)-1.0)*WQPG(L) * WQANCG * WQVO(L,K,3)
         WQA15 = WQA15C+WQA15D+WQA15G
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQA15 =WQA15 + (WQPNM(L)-1.0)*WQPM(L)*WQANCM
     *         *WQVO(L,K,IDNOTRVA)
       END IF
C J.S.
        WQV(L,K,15)=SCBWQ(L)*( WQVO(L,K,15) + DTWQ*WQRR(L)
     *     + DTWQO2*( WQA15
     *     - WQANDC*WQDENIT(L)*WQVO(L,K,6) + WQNIT(L)*WQVO(L,K,14) ))
     *     +(1.-SCBWQ(L))*WQVO(L,K,15)
         WQVO(L,K,15) = WQVO(L,K,15)+WQV(L,K,15)
       END DO      
      
      endif
!
C
C18 COD: WQO18(L)=DTWQO2*WQO18
C
       DO L=LF,LL
         WQKK(L) = 1.0 / (1.0 - WQO18(L))
         WQRR(L) = (WQWDSL(L,K,18)+WQWPSL(L,K,18)) * VOLWQ(L)
       END DO
C
       IF (K.EQ.1) THEN
         DO L=LF,LL
           WQRR(L) = WQRR(L) + WQBFCOD(L)*DZWQ(L)
         END DO
       END IF
C
       DO L=LF,LL
       WQRR(L) = WQVO(L,K,18) + DTWQ*WQRR(L) + WQO18(L)*WQVO(L,K,18)
       WQV(L,K,18)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,18)
        WQVO(L,K,18) = WQVO(L,K,18)+WQV(L,K,18)
C
C19 O2: WQD19=-WQAOCR*WQKHR(L),WQL19=-WQAONT*WQNIT(L),WQO19=WQO18
C: WQKRDOS(L)=-WQP19(L)*WQDOS
C
        WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP19(L))
        WQRR(L) = (WQWDSL(L,K,19)+WQWPSL(L,K,19)) * VOLWQ(L)
        xDOpsl(L,K) = xDOpsl(L,K) + WQRR(L)*DTWQ*DZCWQ(K)*HPWQ(L)
        xDOall(L,K) = xDOall(L,K) + WQRR(L)*DTWQ*DZCWQ(K)*HPWQ(L)
       END DO
C
       IF (K.EQ.KC) THEN
         DO L=LF,LL
           WQRR(L) = WQRR(L) + WQKRDOS(L)
             LL0=LIJW(46,103)
  !        IF(L.EQ.LL0)write(*,*)WQKRDOS(L),WQVO(L,K,19)      
         END DO
       END IF
C
       IF (K.EQ.1) THEN
         DO L=LF,LL
           WQRR(L) = WQRR(L) + WQBFO2(L)*DZWQ(L)
           xDOsod(L,K) = xDOsod(L,K) + WQBFO2(L)*DTWQ
           xDOall(L,K) = xDOall(L,K) + WQBFO2(L)*DTWQ
         END DO
       END IF
C
       DO L=LF,LL
         O2WQ = MAX(WQVO(L,K,19), 0.0)
         WQTTC = 0.0
         WQTTD = 0.0
         WQTTG = (1.3 - 0.3*WQPNG(L)) * WQPG(L)*DO_G
         xDOppB(L,K) = xDOppB(L,K) +  
     *     ( WQTTG*WQVO(L,K,3) ) * WQAOCR*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         xDOall(L,K) = xDOall(L,K) +  
     *     ( WQTTG*WQVO(L,K,3) ) * WQAOCR*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)

         xmrm = CFCDGWQ*O2WQ*WQBMG(L)/(WQKHRG+O2WQ+ 1.E-18)
         WQA19G = WQTTG - xmrm
         xDOrrB(L,K) = xDOrrB(L,K) - xmrm*WQVO(L,K,3) * WQAOCR*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         xDOall(L,K) = xDOall(L,K) - xmrm*WQVO(L,K,3) * WQAOCR*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         WQA19 = (  WQA19G*WQVO(L,K,3) ) * WQAOCR
C
C add macalgal source     J.S.
C
C Modified by MRM 05/23/99 to allow different AOCR constants to be applied
C   to photosynthesis and respiration terms for macroalgae:
C     WQAOCRpm = AOCR applied to macroalgae photosynthesis term
C     WQAOCRrm = AOCR applied to macroalgae respiration term
C
         IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
           WQTTM = (1.3 - 0.3*WQPNM(L)) * WQPM1(L)                  !Modify Macalgae growth and DO production by reduce nutrient limits
           xmrm = (1.0-WQFCDM)*o2wq*WQBMM(L)/(WQKHRM+O2WQ+ 1.E-18)
           WQA19A = WQTTM * WQVO(L,K,IDNOTRVA) * WQAOCRpm -
     *       xmrm *  WQVO(L,K,IDNOTRVA) * WQAOCRrm
           WQA19 = WQA19 + WQA19A
           xDOppM(L,K) = xDOppM(L,K) +
     *         WQTTM*WQVO(L,K,IDNOTRVA)*WQAOCRpm*DTWQO2*DZCWQ(K)*HPWQ(L)
           xDOall(L,K) = xDOall(L,K) +
     *         WQTTM*WQVO(L,K,IDNOTRVA)*WQAOCRpm*DTWQO2*DZCWQ(K)*HPWQ(L)
           xDOrrM(L,K) = xDOrrM(L,K) -
     *         xmrm*WQVO(L,K,IDNOTRVA)*WQAOCRrm*DTWQO2*DZCWQ(K)*HPWQ(L)
           xDOall(L,K) = xDOall(L,K) -
     *         xmrm*WQVO(L,K,IDNOTRVA)*WQAOCRrm*DTWQO2*DZCWQ(K)*HPWQ(L)
         END IF
C J.S.
         WQRR(L) = WQVO(L,K,19) + DTWQ*WQRR(L) + DTWQO2*( WQA19
     *     - WQAOCR*WQKHR(L)*WQVO(L,K,6)*DO_R
     *     - WQAONT*WQNIT(L)*WQVO(L,K,14)
     *     + WQO18(L)*WQVO(L,K,18) + WQP19(L)*WQVO(L,K,19) )
       WQV(L,K,19)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,19)
C
C MRM do not allow D.O. to go negative:
C
         WQV(L,K,19) = max (WQV(L,K,19), 0.0)
         WQVO(L,K,19) = WQVO(L,K,19)+WQV(L,K,19)
C
C compute and save D.O. deficit:
C
         xmrm = WQDOS(L) - WQV(L,K,19)
         xDOdef(L,K) = xDOdef(L,K) + xmrm*DTWQ*DZCWQ(K)*HPWQ(L)
         IF (K.EQ.KC) THEN
          xDOkar(L,K) = xDOkar(L,K) + WQKRDOS(L)*DTWQ*DZCWQ(K)*HPWQ(L)+
     *        WQP19(L)*WQVO(L,K,19)*DTWQO2*DZCWQ(K)*HPWQ(L)
          xDOall(L,K) = xDOall(L,K) + WQKRDOS(L)*DTWQ*DZCWQ(K)*HPWQ(L)+
     *    WQP19(L)*WQVO(L,K,19)*DTWQO2*DZCWQ(K)*HPWQ(L)
         END IF
         xDOdoc(L,K)=xDOdoc(L,K) - WQAOCR*WQKHR(L)*WQVO(L,K,6)*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         xDOall(L,K)=xDOall(L,K) - WQAOCR*WQKHR(L)*WQVO(L,K,6)*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         xDOnit(L,K)=xDOnit(L,K) - WQAONT*WQNIT(L)*WQVO(L,K,14)*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         xDOall(L,K)=xDOall(L,K) - WQAONT*WQNIT(L)*WQVO(L,K,14)*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         xDOcod(L,K)=xDOcod(L,K) - WQO18(L)*WQVO(L,K,18)*DTWQO2
c        xDOcod(L,K)=xDOcod(L,K) + WQO18(L)*WQVO(L,K,18)             ! codemike.011701, why?
     *     *DZCWQ(K)*HPWQ(L)
         xDOall(L,K)=xDOall(L,K) - WQO18(L)*WQVO(L,K,18)*DTWQO2
c        xDOall(L,K)=xDOall(L,K) + WQO18(L)*WQVO(L,K,18)             ! codemike.011701, why?
     *     *DZCWQ(K)*HPWQ(L)
         xDOdz(L,K) = xDOdz(L,K) + DZCWQ(K)*HPWQ(L)
       END DO

C

C
C end of DO K=KC,1,-1  loop           
C end of DO ND=1,NDMWQ loop
C
      END DO                          ! #200

C
C increment counter for limitation and xDOxxx DO component arrays:
C
      TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/TIDALP
      timesum3 = timesum3 + TIMTMP
      NLIM = NLIM + 1
C
C compute WQCHL,WQTAMP,WQPO4D,WQSAD at a new time step: WQCHLx=1/WQCHLx
C
      DO K=1,KC
        DO L=2,LA
          WQCHL(L,K) =  WQV(L,K,3)*WQCHLG(I)
        END DO
      END DO
      
      IF (IWQSRP.EQ.2) THEN          ! SEDT
        DO K=1,KC
          DO L=2,LA
            WQPO4D(L,K) = WQV(L,K,10) / (1.0 + WQKPO4P*SEDTWQ(L,K))
            WQSAD(L,K)  = WQV(L,K,17) / (1.0 + WQKSAP*SEDTWQ(L,K))
          END DO
        END DO
       ELSE                                ! SEDT
        DO K=1,KC
          DO L=2,LA
            WQPO4D(L,K) = WQV(L,K,10)
            WQSAD(L,K)  = WQV(L,K,17)
          END DO
        END DO
      END IF
C
C  Coupling to sediment model
C: evaluate dep. flux using new values cause implicit scheme is used in  
C  SPM
C
 
      Fflux=1.0                       
c
      IF (IWQBEN.EQ.1) THEN           ! #400, Ji, 9/18/02

       LF=2
       LL=LA
       DO L=LF,LL
         IMWQZ = IWQZMAP(L,1)
         WQDFBC(L) = 0.0
         WQDFBD(L) = 0.0 
         WQDFBG(L) = SCBWQ(L)*WQWSG(IMWQZ)*WQV(L,1,3)*Fflux
     +             +WQWSM*DZWQ(L)*WQV(L,1,IDNOTRVA)*Fflux
         WQDFRC(L) = 0.0
         if(IWQS.EQ.1) then
         WQDFLC(L) = SCBWQ(L)*WQWSLP(IMWQZ)*WQV(L,1,5)*(1.0-FSTOC)   
         else
         WQDFLC(L) = SCBWQ(L)*WQWSLP(IMWQZ)*WQV(L,1,6)*(1.0-FSTOC)      
         endif
         WQDFRP(L) = 0.0
         if(IWQS.EQ.1) then     
         WQDFLP(L) = SCBWQ(L)*WQWSLP(IMWQZ)*WQV(L,1,8)*(1.0-FSTOP) 
         else
         WQDFLP(L) = SCBWQ(L)*WQWSLP(IMWQZ)*WQV(L,1,9)*(1.0-FSTOP)         
         endif
         WQDFRN(L) = 0.0
         if(IWQS.EQ.1) then 
         WQDFLN(L) = SCBWQ(L)*WQWSLP(IMWQZ)*WQV(L,1,12)*(1.0-FSTON)
         else
          WQDFLN(L) = SCBWQ(L)*WQWSLP(IMWQZ)*WQV(L,1,13)*(1.0-FSTON)        
         endif
         IF (IWQSI.EQ.1) WQDFSI(L) = SCBWQ(L)*WQWSD(IMWQZ)*WQV(L,1,16)
       END DO

         IF (IWQSRP.EQ.2) THEN             ! SEDT
          DO ND=1,NDMWQ
          LF=2+(ND-1)*LDMWQ
          LL=LF+LDM-1
          DO L=LF,LL
           WQDFLP(L) = SCBWQ(L)*( WQDFLP(L)+WQSEDO(NS)*( WQV(L,1,10)
     *                         -WQPO4D(L,1) ) )
           IF (IWQSI.EQ.1) WQDFSI(L) = SCBWQ(L)*( WQDFSI(L)
     *        + WQSEDO(NS)*( WQV(L,1,17)-WQSAD(L,1) ) )
          END DO
          END DO
         END IF                                  ! SEDT
C
      END IF                                       ! #400
C
C DIURNAL DO ANALYSIS
C
      IF(NDDOAVG.GE.1) THEN                     ! =0, wq3dwc.inp, C4
        OPEN(1,FILE='diurndo.out',ACCESS='APPEND')
        NDDOCNT=NDDOCNT+1
        NSTPTMP=NDDOAVG*NTSPTC/2
        RMULTMP=1./FLOAT(NSTPTMP)
C
        DO K=1,KC
         DO L=2,LA
          DDOMAX(L,K)=MAX(DDOMAX(L,K),WQV(L,K,19))
          DDOMIN(L,K)=MIN(DDOMIN(L,K),WQV(L,K,19))
         END DO
        END DO
C
        IF(NDDOCNT.EQ.NSTPTMP) THEN
          NDDOCNT=0
          TIME=DT*FLOAT(N)+TCON*TBEGIN
          TIME=TIME/TCON
          WRITE(1,1111)N,TIME
          DO L=2,LA
           WRITE(1,1112)ILW(L),JLW(L),(DDOMIN(L,K),K=1,KC),
     $                               (DDOMAX(L,K),K=1,KC)
          END DO
          DO K=1,KC
           DO L=2,LA
            DDOMAX(L,K)=-1.E6
            DDOMIN(L,K)=1.E6
           END DO
          END DO
        END IF
C
        CLOSE(1)
      END IF
C
C LIGHT EXTINCTION ANALYSIS
C
      IF(NDLTAVG.GE.1) THEN                 
        OPEN(1,FILE='light.out',ACCESS='APPEND')
        NDLTCNT=NDLTCNT+1
        NSTPTMP=NDLTAVG*NTSPTC/2
        RMULTMP=1./FLOAT(NSTPTMP)
C
        DO K=1,KC
         DO L=2,LA
          RLIGHT1=WQKEB(IMWQZT(L))+WQKETSS(IMWQZT(L))*SEDTWQ(L,K)
C
C MRM 05/12/1999 use Riley (1956) equation to compute light extinction
C     as a function of CHL conc. if WQKECHL is less than zero:
C
C          RLIGHT2=WQKECHL*WQCHL(L,K)
          xmrm = WQKECHL(IMWQZT(L))*WQCHL(L,K)
          IF (WQKECHL(IMWQZT(L)).LT. 0.0) THEN
            xmrm = 0.054*WQCHL(L,K)**0.6667 + 0.0088*WQCHL(L,K)
          END IF
          RLIGHT2 = xmrm
          RLIGHTT(L,K)=RLIGHTT(L,K)+RLIGHT1
          RLIGHTC(L,K)=RLIGHTC(L,K)+RLIGHT1+RLIGHT2
         END DO
        END DO
C
        IF(NDLTCNT.EQ.NSTPTMP) THEN
          NDLTCNT=0
          TIME=DT*FLOAT(N)+TCON*TBEGIN
          TIME=TIME/TCON
           DO K=1,KC
           DO L=2,LA
            RLIGHTT(L,K)=RMULTMP*RLIGHTT(L,K)
            RLIGHTC(L,K)=RMULTMP*RLIGHTC(L,K)
           END DO
          END DO
          WRITE(1,1111)N,TIME
          DO L=2,LA
           WRITE(1,1113)ILW(L),JLW(L),(RLIGHTT(L,K),K=1,KC),
     $                               (RLIGHTC(L,K),K=1,KC)
          END DO
          DO K=1,KC
           DO L=2,LA
            RLIGHTT(L,K)=0.
            RLIGHTC(L,K)=0.
           END DO
          END DO
        END IF
C
        CLOSE(1)
      END IF
C
C
 1111 FORMAT(I12,F10.4)
 1112 FORMAT(2I5,12F7.2)
 1113 FORMAT(2I5,20E12.4)
c
      RETURN
      END
