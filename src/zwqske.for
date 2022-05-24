C 
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WQSKE(LA,KC,DT,N,TCON,TBEGIN,TIDALP,NTSPTC,SECDLAST)
C
C**********************************************************************C
C
C  This is main subroutin to solve all kinetic Eqs. 
C  from K=KC (surface layer) to K=1 (bottom).
C: after computing new values, store WQVO+WQV into WQVO(L,K,NWQV).  
C
C
C  ORGINALLY CODED BY K.-Y. PARK
C  OPTIMIZED AND MODIFIED BY J.M. HAMRICK
C  LAST MODIFIED BY J.M. HAMRICK  7 APRIL 1997
C
C  MODIFIED ON 4/1/2007 BY AEE
C
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
      REAL*8 SECDLAST 
      DIMENSION WQDOS(LCMWQ), WQI0BOT(LCMWQ),
     &  DNH4c(LCMWQ),DNH4g(LCMWQ),DNH4d(LCMWQ),
     &  DNO3c(LCMWQ),DNO3g(LCMWQ),DNO3d(LCMWQ)

      TIME=DT*FLOAT(N)+TCON*TBEGIN
      IF(NCSTEP.GT.0) TIMTMP=SECDLAST/TCON+TBEGIN 
      CNS1=2.718                                                         ! (3-1e)

      NS=1

	NDMWQ=1
	LDM=LA-1
	DO L=2,LA
	DNH4c(L)=1
	DNH4d(L)=1
	DNH4g(L)=1
	DNO3c(L)=1
	DNO3d(L)=1
	DNO3g(L)=1
	ENDDO
	
!	ITRYCHLA=1   ! JS 2018  use vary C/Chl ratio
      IF(ITRYCHLA.EQ.3) THEN
	KECST=30           ! also active respiration with growth
	ELSEIF(ITRYCHLA.EQ.2)THEN
	KECST=40
	ENDIF
C
C Initial solar radiation at top of surface layer  (J.S. 5/4/99)
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
      DO L=2,LA
        WQI0BOT(L)=WQI0 * pSHADE(L)
        xI00(L)= xI00(L)+WQI0 * pSHADE(L) *DTWQ          ! J.S. 8/16/08
        xHXY(L)=xHXY(L)+HPWQ(L)
      END DO
C
C Start of DO ND=1,NDMWQ loop
C
      DO ND=1,NDMWQ                               
      LF=2+(ND-1)*LDMWQ          
      LL=LF+LDM-1        

C Start of DO K=KC,1,-1  loop

      DO K=1,KC
       DO L=2,LA
       IF(ITRYCHLA.GE.2) THEN                             ! JS 01/2018 using verying C:Chl a ratio
       WQCHLD(IWQZMAP(L,KC))=1/(KECST+90*exp(-1.19*WQketot(L,KC)))*1000      
       WQCHLC(IWQZMAP(L,KC))=1/(KECST+150*exp(-1.18*WQketot(L,KC)))*1000  
       WQCHLG(IWQZMAP(L,KC))=WQCHLC(IWQZMAP(L,KC)) 
       ENDIF
          WQCHL(L,K) = WQV(L,K,1)*WQCHLC(IWQZMAP(L,KC)) + 
     *    WQV(L,K,2)*WQCHLD(IWQZMAP(L,KC))
     *      + WQV(L,K,3)*WQCHLG(IWQZMAP(L,KC))
       END DO
      END DO

      DO K=KC,1,-1                                
C
C DZWQ=1/h, VOLWQ=1/VOL
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
       DO L=LF,LL
        TWQ(L)=TEMWQ(L,K)
        SWQ(L)=MAX(SALWQ(L,K), 0.0) 
        DZWQ(L) = 1.0 / (DZCWQ(K)*HPWQ(L))
        VOLWQ(L) = DZWQ(L) / DXYPWQ(L)
        IMWQZT(L)=IWQZMAP(L,K)
       END DO

       IF(WQPMC(4).GT.0.4) THEN           ! If Laf. R. HAB is simulated 2018
       DO L=LF,LL                         ! Add HAB seeds    JS 2018
       IF(TWQ(L).GT.18) THEN
       if(ILW(L)>314.and.JLW(L)>74) then
        if(ILW(L)<323.and.JLW(L)<82) then
          WQV(L,K,1)=max(WQV(L,K,1),WQMCMIN)
          WQVO(L,K,1)=max(WQVO(L,K,1),WQMCMIN)
        endif
       endif
        if(ILW(L)>288.and.JLW(L)>47) then
        if(ILW(L)<291.and.JLW(L)<81) then
          WQV(L,K,1)=min(WQV(L,K,1),1.0)
        endif
        endif
       ENDIF
 !      If(IMWQZT(L).EQ.5) then
 !      if((WQV(L,K,3))>50*0.06)then               !JS 2018
 !      WQV(L,K,1)=max(WQV(L,K,1),5.0*0.06)
 !      WQVO(L,K,1)=max(WQVO(L,K,1),5.0*0.06)
 !      endif
 !      endif                     ! Add freshwater HAB seed when Chl-a > 45 ug/L  JS 2018
       ENDDO
       ENDIF
C 
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
       IPPTC=0
       ISWIMM=0        
       IF(abs(WQWSC(3)) >2.0.or. abs(WQWSC(4))>2.0) ISWIMM=1   ! HAB  add swimming 4/22/2015
       DO L=LF,LL
       IF(ISWIMM.EQ.1) THEN
        IF(IMWQZT(L).EQ.4.or.IMWQZT(L).EQ.3) then              ! J.S. 2015, HAB add to using very settling vel. due to swimming
        WQAVGIO = WQCIA*WQI0 + WQCIB*WQI1 + WQCIC*WQI2
 !       HTEM=HPWQ(L)/KC
 !       HTEM=HTEM/(DTWQO2*2)
         TSWIM3=abs(WQWSC(3)) !0.000050*86400.0 !0.000015*86400.0
         TSWIM4=abs(WQWSC(4)) 
 !        TSWIM=min(TSWIM,HTEM)
         if(IPPTC.EQ.1) then
         write(*,*)'Set settling vel= ',TSWIM,WQAVGIO
         IPPTC=0
         endif
        if(WQAVGIO > 30) then 
        WQWSC(4)=-TSWIM4
        WQWSC(3)=-TSWIM3       
        else
        WQWSC(4)=TSWIM4
        WQWSC(3)=TSWIM3
        endif
        ENDIF
!       WQWSC(4)=0.04
        ENDIF
        
        ACCSET=0
        WQBCSET(L,1) = WQWSC(IMWQZT(L))*DZWQ(L)
     %          + ACCSET*(WQVO(L,K,2)/1.2)          ! HAB modify acceleartion settling. ACCSET=0.0 turn off   J.S. 2018
        WQBDSET(L,1) = WQWSD(IMWQZT(L))*DZWQ(L)
  !   $          +ACCSET*(WQVO(L,K,2)/1.2)  
        WQBGSET(L,1) = WQWSG(IMWQZT(L))*DZWQ(L)
  !   &          +ACCSET*(WQVO(L,K,3)/1.2)           
        WQRPSET(L,1) = WQWSRP(IMWQZT(L))*DZWQ(L)
        WQLPSET(L,1) = WQWSLP(IMWQZT(L))*DZWQ(L)
       END DO
C
       IF (IWQSRP.EQ.1) THEN 
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
        DO L=LF,LL
         WQWSSET(L,1) = WQWSS(IMWQZT(L))*DZWQ(L)
        END DO
       END IF
C
       IF (K.NE.KC) THEN
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
        DO L=LF,LL
         IMWQZT1(L)=IWQZMAP(L,K+1)
        END DO
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
        DO L=LF,LL
         WQBCSET(L,2) = WQWSC(IMWQZT1(L))*DZWQ(L)
         WQBDSET(L,2) = WQWSD(IMWQZT1(L))*DZWQ(L)+
     &        ACCSET*(WQVO(L,K+1,2)/1.2) 
         WQBGSET(L,2) = WQWSG(IMWQZT1(L))*DZWQ(L)
     &        +ACCSET*(WQVO(L,K+1,3)/1.2) 
         WQRPSET(L,2) = WQWSRP(IMWQZT1(L))*DZWQ(L)
         WQLPSET(L,2) = WQWSLP(IMWQZT1(L))*DZWQ(L)
        END DO


        IF (IWQSRP.EQ.1) THEN
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
         DO L=LF,LL
         WQWSSET(L,2) = WQWSS(IMWQZT1(L))*DZWQ(L)
         END DO
        END IF
C
       END IF
C
C find an index for look-up table for temperature dependency
C
       DO L=LF,LL
c       IWQT(L) = 10.0*TWQ(L) +151
        IWQT(L) = 10.0*TWQ(L) + 51               !must be consistent with the look-up table:-5C-> 50C in rwqc1 & wqwin, Ji, 8/29/02
        IF (IWQT(L).LT.1 .OR. IWQT(L).GT.NWQTD) THEN
c MRM +++++++++ added by M. Morton 08/05/98
          timtmp = (dt*float(n) + tcon*tbegin)/86400.0
          OPEN(1,FILE='error.log',ACCESS='APPEND',STATUS='UNKNOWN')
          write(1,911) timtmp, L, ILW(L), JLW(L), K, twq(L)
911       format(/,'ERROR!! invalid water temperature, sub WQSKE',/,
     +     'TIME, L, I, J, K, TWQ(L) = ', f10.5, 4i4, f10.4)
          close(1)
c MRM +++++++++ added by M. Morton 07/24/98
C         PRINT*, 'L, K, TEM(L,K) = ', L,K,TWQ(L)
          WRITE(6,600)ILW(L),JLW(L),K,TWQ(L)
          STOP 'ERROR!! invalid water temperature'
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

       DO L=LF,LL                              
        RNH4WQ = MAX (WQVO(L,K,14), 0.0)
        RNO3WQ = MAX (WQVO(L,K,15), 0.0)
        PO4DWQ = MAX (WQPO4D(L,K),  0.0)
        RDONWQ = MAX (WQVO(L,K,13)/2.0, 0.0) 
        IF(IMWQZT(L).EQ.3.or.IMWQZT(L).EQ.4.or.IMWQZT(L).EQ.7) THEN          
        RNH4NO3c = RNH4WQ + RNO3WQ + RDONWQ                        ! J.S. 2017 HAB. add DON limits
        WQGNC = max(RNH4NO3c/ (WQKHNC+RNH4NO3c+ 1.E-18),AN_Lim)  
        RNH4NO3 = RNH4WQ + RNO3WQ                                  ! no change or algal D and G
        WQGND = max(RNH4NO3 / (WQKHND+RNH4NO3+ 1.E-18),AN_Lim)        
        WQGNG = max(RNH4NO3 / (WQKHNG+RNH4NO3+ 1.E-18),AN_Lim)        
        ELSE
        RNH4NO3 = RNH4WQ + RNO3WQ  
        WQGNC = max(RNH4NO3 / (WQKHNC+RNH4NO3+ 1.E-18),AN_Lim)     ! EFDC WQ manual, P10, Eq. (3-1c)          
        WQGND = max(RNH4NO3 / (WQKHND+RNH4NO3+ 1.E-18),AN_Lim)
        WQGNG = max(RNH4NO3 / (WQKHNG+RNH4NO3+ 1.E-18),AN_Lim)         
        ENDIF
        WQGPC = max(PO4DWQ / (WQKHPC+PO4DWQ+ 1.E-18),P_Lim)
        WQGPD = max(PO4DWQ / (WQKHPD+PO4DWQ+ 1.E-18),P_Lim)
        WQGPG = max(PO4DWQ / (WQKHPG+PO4DWQ+ 1.E-18),P_Lim)
        xlimNc(L,K) = xlimNc(L,K) + WQGNC
        xlimNd(L,K) = xlimNd(L,K) + WQGND
        xlimNg(L,K) = xlimNg(L,K) + WQGNG
        xlimPc(L,K) = xlimPc(L,K) + WQGPC
        xlimPd(L,K) = xlimPd(L,K) + WQGPD
        xlimPg(L,K) = xlimPg(L,K) + WQGPG
C
C Modified by J.S.(5/5/98) for Macalgae
C
        IF(IDNOTRVA.GT.0 .AND. K.EQ.1) THEN              ! macroalgae
          WQGNM = RNH4NO3 / (WQKHNM+RNH4NO3 + 1.E-18)
          WQGPM = PO4DWQ / (WQKHPM+PO4DWQ + 1.E-18)
          WQF1NM = MIN(WQGNM, WQGPM)
          xlimNm(L,K) = xlimNm(L,K) + WQGNM
          xlimPm(L,K) = xlimPm(L,K) + WQGPM
        END IF
C J.S.
        WQF1NC = MIN(WQGNC, WQGPC)                   ! WQ manual, Eq. (3-1c)
C  IWQSI   = switch to activate silica state variables (0=off; 1=activated), wq3dwc.inp
        IF (IWQSI.EQ.1) THEN
          SADWQ = MAX (WQSAD(L,K), 0.0)              ! dissolved SA
          WQGSD = SADWQ / (WQKHS+SADWQ+ 1.E-18)
          WQF1ND = MIN(WQGND, WQGPD, WQGSD)          ! Eq. (3-1d)
         ELSE
          WQF1ND = MIN(WQGND, WQGPD)
        END IF
        WQF1NG = MIN(WQGNG, WQGPG)                   ! Eq. (3-1c)

       IF(IDNOTRVA.GT.0) THEN
       PO4DWQ = MAX (WQPO4D(L,K), 0.0)
       END IF
C
C algal growth: light, WQHT(K)=REAL(KC-K)/REAL(KC)
C In C&C, F2Ic=F2Ic/FCYAN, factor to allow cyanobacteria mat formation
C
C MRM 05/12/1999 use Riley (1956) equation to compute light extinction
C     as a function of CHL conc. if WQKECHL is less than zero:
C
        xmrm = WQKECHL(IMWQZT(L))*WQCHL(L,K)
        IF (WQKECHL(IMWQZT(L)) .LT. 0.0) THEN
          xmrm = 0.054*WQCHL(L,K)**0.6667 + 0.0088*WQCHL(L,K)
        END IF
        
C        WQKESS = WQKEB(IMWQZT(L))+WQKETSS*SEDT(L,K)+WQKECHL*WQCHL(L,K)
          WQKESS = WQKEB(IMWQZT(L))+WQKETSS(IMWQZT(L))*SEDTWQ(L,K) 
     &     +WQKESAT(IMWQZT(L))*SALWQ(L,K)
     &     +WQKECHL(IMWQZT(L))*WQCHL(L,K) ! Use surface Chl a   
c	write(903,9112) L, K, IMWQZT(L), WQKESS, WQKEB(IMWQZT(L)),
c     *  WQKETSS,SEDT(L,K),xmrm 
        WQKESS1 = WQKESS

        IF (K.NE.KC) THEN
C         WQKESS1=WQKEB(IMWQZT(L))+WQKETSS*SEDT(L,KC)+WQKECHL*WQCHL(L,KC)
          xmrm = WQKECHL(IMWQZT(L))*WQCHL(L,KC)
          IF (WQKECHL(IMWQZT(L)).LT. 0.0) THEN
            xmrm = 0.054*WQCHL(L,KC)**0.6667 + 0.0088*WQCHL(L,KC)
          END IF
9112  format(3i6,999e12.4)
!        WQKESS1=WQKEB(IMWQZT(L))+WQKETSS(IMWQZT(L))*SEDTWQ(L,KC) + xmrm
c      write(902,9112) L, K,IMWQZT(L),WQKEB(IMWQZT(L)),WQKETSS,SEDT(L,KC)
c     *  ,xmrm, WQKESS1
        END IF
C 
C compute Secchi depth for use as output variable:
C
C        WQketot(L,K) = wqsdcoef(imwqzt(L)) / wqkess
        WQketot(L,K) = WQketot(L,K)+wqkess

       IF(ITRYCHLA.GE.2) THEN                             ! JS 01/2018 using verying C:Chl a ratio
       WQCHLD(IWQZMAP(L,KC))=1/(KECST+90*exp(-1.19*WQketot(L,KC)))*1000      
       WQCHLC(IWQZMAP(L,KC))=1/(KECST+150*exp(-1.18*WQketot(L,KC)))*1000  
       WQCHLG(IWQZMAP(L,KC))=WQCHLC(IWQZMAP(L,KC)) 
       ENDIF
        
        WQAVGIO = WQCIA*WQI0 + WQCIB*WQI1 + WQCIC*WQI2    ! J.S. 2014 Use daily averaged data

C
! J.S. 2014 
!        IF (IWQSUN .EQ. 2) THEN                !use hourly solar radiation from aser.inp
!          WQAVGIO = WQCIA*WQI1 + WQCIB*WQI2 + WQCIC*WQI3
!        END IF
!        WQAVGIO = WQAVGIO * pSHADE(L)        
        WQAVGIO= WQI0* pSHADE(L)                                      ! J.S. 2014 Use hourly data         

C        IF (IWQSUN.EQ.2) WQAVGIO=WQIS0
C
c      write(901,9111) L, K, WQAVGIO,WQKESS1,WQDOPC,WQISMIN
9111  format(2i6,999e12.4)

        if(WQISMIN.LE.5) then
        WQISMIN1=40      ! J.S. 2014 force it to 40 
        WQISC = MAX( WQAVGIO*EXP(-WQKESS1*WQDOPC), WQISMIN1 )
        WQISD = MAX( WQAVGIO*EXP(-WQKESS1*WQDOPD), WQISMIN1 )
        WQISG = MAX( WQAVGIO*EXP(-WQKESS1*WQDOPG), WQISMIN1 )
        WQTT1 = (CNS1 * WQFD * DZWQ(L)) / WQKESS                         ! (3-1e)
C        WQFDI0 = - WQI0 / (WQFD+ 1.E-18)
        WQFDI0 = - WQI0BOT(L) / (WQFD + 1.E-18)
C
        WQFDC = WQFDI0 / (WQISC + 1.E-18)                                 ! (3-1f)
        WQFDD = WQFDI0 / (WQISD + 1.E-18)
        WQFDG = WQFDI0 / (WQISG + 1.E-18)
        WQHTT = WQHT(K) * HPWQ(L)
C
        WQTTB = EXP( -WQKESS * (WQHTT+1.0/DZWQ(L)) )                     ! (3-1f)
        WQTTT = EXP( -WQKESS * WQHTT )                                   ! (3-1g)
        WQF2IC = WQTT1 * (EXP(WQFDC*WQTTB) - EXP(WQFDC*WQTTT))           ! (3-1e)
        WQF2ID = WQTT1 * (EXP(WQFDD*WQTTB) - EXP(WQFDD*WQTTT))
        WQF2IG = WQTT1 * (EXP(WQFDG*WQTTB) - EXP(WQFDG*WQTTT))

        elseif(WQISMIN.LE.10) then
! Try hourly data J.S. 7-9-2014
        WQHTT = WQHT(K) * HPWQ(L)    ! Current depth
        WQAVGIO= WQI0* pSHADE(L) 
        WQISC = MAX( WQAVGIO*EXP(-WQKESS1*WQDOPC), WQISMIN )
        WQISD = MAX( WQAVGIO*EXP(-WQKESS1*WQDOPD), WQISMIN )
        WQISG = MAX( WQAVGIO*EXP(-WQKESS1*WQDOPG), WQISMIN )
                
        WQISC1 = WQAVGIO*EXP(-WQKESS1*(WQHTT+0.5/DZWQ(L)))/ WQISC     ! 
        WQISD1 = WQAVGIO*EXP(-WQKESS1*(WQHTT+0.5/DZWQ(L)))/ WQISD
        WQISG1 = WQAVGIO*EXP(-WQKESS1*(WQHTT+0.5/DZWQ(L)))/ WQISG
        WQF2IC = WQISC1 * EXP(-WQISC1+1)     
        WQF2ID = WQISD1 * EXP(-WQISD1+1)  
        WQF2IG = WQISG1 * EXP(-WQISG1+1) 
        
        else
        ! J.S. Test new light calculation 2014  HAB
!
        WQHTT = WQHT(K) * HPWQ(L)    ! Current depth
        WQIZ=WQI0*EXP(-WQKESS1 * (WQHTT+0.5/DZWQ(L)) )  ! light at Z
!	  WQI_0=  WQI0*WQI0+ PARADJ*2.065*350             
c        WQF2IC = WQIZ/sqrt( WQI0*WQI0+WQISMIN*WQISMIN)  ! Using WQISMIN read from Card 10 
c        WQF2ID = WQIZ/sqrt( WQI0*WQI0+WQISMIN*WQISMIN)  
c        WQF2IG = WQIZ/sqrt( WQI0*WQI0+WQISMIN*WQISMIN)       

!        WQF2IC = WQIZ/sqrt( WQIZ*WQIZ+WQISMIN*WQISMIN)  ! Using WQISMIN read from Card 10 
!        WQF2ID = WQIZ/sqrt( WQIZ*WQIZ+WQISMIN*WQISMIN)  
!        WQF2IG = WQIZ/sqrt( WQIZ*WQIZ+WQISMIN*WQISMIN) 
!        if(IMWQZT(L).EQ.4.or.IMWQZT(L).eq.3)then       
        WQF2IC = WQIZ/sqrt( WQIZ*WQIZ+WQISMINC*WQISMINC)   ! Using WQISMIN read from Card 10 4/22/15, 2018
!        else
!        WQF2IC = WQIZ/sqrt( WQIZ*WQIZ+WQISMIN*2*WQISMIN)  ! Using WQISMIN read from Card 10 4/22/15
!        endif
        WQF2ID = WQIZ/sqrt( WQIZ*WQIZ+WQISMIN*WQISMIN)  
        WQF2IG = WQIZ/sqrt( WQIZ*WQIZ+WQISMING*WQISMING) 
        
!	  WQF2IC = WQIZ/sqrt( WQIZ*WQIZ+(WQPMC(IMWQZT(L))/0.01)**2)  ! Using WQISMIN read from Card 10 
!        WQF2ID = WQIZ/sqrt( WQIZ*WQIZ+(WQPMD(IMWQZT(L))/0.01)**2)  
!        WQF2IG = WQIZ/sqrt( WQIZ*WQIZ+(WQPMG(IMWQZT(L))/0.01)**2)     
!----------------------------------------------------
        endif
        
        xlimIc(L,K) = xlimIc(L,K) + WQF2IC
        xlimId(L,K) = xlimId(L,K) + WQF2ID
        xlimIg(L,K) = xlimIg(L,K) + WQF2IG
!
        xKe(L)=xKe(L)+WQKESS1/Real(KC)   ! Record Ke  J.S. 8/16/08
!        IF(K.EQ.KC) xKe(L)=WQKESS1       ! Record Ke  at surface J.S. 8/15/2014       
C
C update solar radiation at bottom of this layer
C
        WQI0BOT(L)=WQI0BOT(L)*exp(-WQKESS*(1.0/DZWQ(L)))
C
C Modified by J.S.(5/5/98) for Macalgae
C
       IF (IDNOTRVA.GT.0 .AND. K.EQ.1) THEN
        WQHTT = WQHT(K) * HPWQ(L)    ! Current depth
        WQIZ=WQI0*EXP(-WQKESS1 * (WQHTT+0.5/DZWQ(L)) )  
        WQF2IM =WQIZ/sqrt(WQIZ*WQIZ+40)
      
      ! WQFDI0 = - WQI0BOT(L) / (WQFD + 1.E-18)
      ! WQISM = MAX( WQAVGIO*EXP(-WQKESS1*WQDOPM), WQISMINC )
      ! WQFDM = WQFDI0 / (WQISM + 1.E-18)
      ! WQF2IM = WQTT1 * (EXP(WQFDM*WQTTB) - EXP(WQFDM*WQTTT))
       
       WQPM(L)= WQPMM(IMWQZT(L))*WQF1NM*WQF2IM*WQTDGM(IWQT(L),IMWQZT(L))
       xlimIm(L,K) = xlimIm(L,K) + WQF2IM
       xlimTm(L,K) = xlimTm(L,K) + WQTDGM(IWQT(L),IMWQZT(L))                    ! (3-1k)
       END IF
        xlimTc(L,K) = xlimTc(L,K) + WQTDGC(IWQT(L),IMWQZT(L))
        xlimTd(L,K) = xlimTd(L,K) + WQTDGD(IWQT(L),IMWQZT(L))
        xlimTg(L,K) = xlimTg(L,K) + WQTDGG(IWQT(L),IMWQZT(L))
C
C J.S. Salinity toxisity 2014
C
C: WQSTOX=WQSTOX**2
C
C  IWQSTOX = cyanobacteria salinity toxicity switch (0=no toxicity; 1=toxicity)
        TOX_f=1                       ! HAB change salinity related motality JS. 2015/2018
        TOX_m=1
        IF (IWQSTOX.EQ.1) THEN
         TOX_f1= SWQ(L)**2  / (1.0 + SWQ(L)**2 +1.E-12)
         TOX_f2= WQSTOX**2 / (WQSTOX**2 + SWQ(L)+1.E-12)       
         if (IMWQZT(L).EQ.2) then        ! J.S. 2015 HAB apply salinty toxcity to middl section only
         TOX_m=1.0      
         TOX_f=TOX_f2
         endif
         if (IMWQZT(L).EQ.7) then        ! J.S. 2015 HAB apply salinty toxcity to middl section only
         TOX_m=1.0      
         TOX_f=TOX_f1
         endif 
         
         ass0=(SWQ(L)-34.0)*(SWQ(L)-28.0)    ! J.S. add salinity effect for HAB
         if( SWQ(L).LE.34.0) then
         ass0=exp(-0.0024*ass0)
         else
          ass0=exp(-0.0222*ass0)  
         endif      
!       IF (IWQSTOX.EQ.1) THEN
!       WQF4SC = WQSTOX / (WQSTOX + SWQ(L)*SWQ(L)+1.E-12)
!       WQPC(L)=WQPMC(IMWQZT(L))*WQF1NC*WQF2IC*WQTDGC(IWQT(L),IMWQZT(L))
!     & *WQF4SC ! (3-1b)
!       ELSE
!       WQPC(L)= WQPMC(IMWQZT(L))*WQF1NC*WQF2IC*WQTDGC(IWQT(L),IMWQZT(L))      ! (3-1a), cyanobacteria
!       END IF
       IF(ITRYCHLA.GE.2) THEN                ! JS 01/2018
       WQCHLD(IWQZMAP(L,KC))=(KECST+ 90*exp(-1.19*WQketot(L,KC)))      
       WQCHLC(IWQZMAP(L,KC))=(KECST+150*exp(-1.18*WQketot(L,KC)))   ! summer
       WQCHLG(IWQZMAP(L,KC))=WQCHLC(IWQZMAP(L,KC)) 
       PmaxD=WQPMD(IMWQZT(L))/WQCHLD(IWQZMAP(L,KC))             ! Change growth rate based on C:Chla ratio
       PmaxC=WQPMC(IMWQZT(L))/WQCHLC(IWQZMAP(L,KC))
       PmaxG=WQPMG(IMWQZT(L))/WQCHLG(IWQZMAP(L,KC))  
       WQPD(L)= PmaxD*WQF1ND*WQF2ID*WQTDGD(IWQT(L),IMWQZT(L))
     #          *TOX_f        !         diatom
       WQPG(L)= PmaxG*WQF1NG*WQF2IG*WQTDGG(IWQT(L),IMWQZT(L))        !         green
     &          *TOX_f !TOX_m
       WQPC(L)= PmaxC*WQF1NC*WQF2IC*WQTDGC(IWQT(L),IMWQZT(L))
     +          *TOX_f*ass0           
       ELSE                                   ! JS 01/2018
       WQPD(L)= WQPMD(IMWQZT(L))*WQF1ND*WQF2ID*WQTDGD(IWQT(L),IMWQZT(L))
     #          *TOX_f        !         diatom
       WQPG(L)= WQPMG(IMWQZT(L))*WQF1NG*WQF2IG*WQTDGG(IWQT(L),IMWQZT(L))        !         green
     &          *TOX_f !TOX_m
       WQPC(L)= WQPMC(IMWQZT(L))*WQF1NC*WQF2IC*WQTDGC(IWQT(L),IMWQZT(L))
     +          *TOX_f*ass0
       ENDIF
c
!       if (IMWQZT(L).EQ.7) then        ! JS. 2015 HAB apply salinty toxcity to middl section only
!       TOX_f2= SWQ(L)/(SWQ(L)+1) !WQSTOX*2 / (WQSTOX*2 + SWQ(L)+1.E-12)     
!       WQPC(L)= WQPMC(IMWQZT(L))*WQF1NC*WQF2IC*WQTDGC(IWQT(L),IMWQZT(L))
!     $       *(1-TOX_f2)
!       endif
      ENDIF
C
c--Check capacity for growth JS 2015 HAB
C
        ICAPACITY=0
        IF(ICAPACITY.EQ.1) THEN
        DIN_mx= (MAX (WQVO(L,K,14), 0.0)+MAX (WQVO(L,K,15), 0.0))   ! DIN
        DIN_c=WQANCC1(IMWQZT(L))*WQPC(L)*WQVO(L,K,1)*DTWQO2 !DTWQ   
        DIN_g=WQANCG1(IMWQZT(L))*WQPG(L)*WQVO(L,K,2)*DTWQO2 !DTWQ   
        DIN_d=WQANCD1(IMWQZT(L))*WQPD(L)*WQVO(L,K,3)*DTWQO2 !DTWQ  
        DIN_t=DIN_c+DIN_g+DIN_d
        DNH4c(L)=1
        DNH4d(L)=1
        DNH4g(L)=1
        IF(DIN_t.GT.DIN_mx) THEN  ! Limiting capacity
        afac=DIN_mx/DIN_t
        a_t=WQVO(L,K,1)+WQVO(L,K,2)+WQVO(L,K,3)+1.0e-8
        DNH4c(L)=afac*WQVO(L,K,1)/a_t
        DNH4d(L)=afac*WQVO(L,K,2)/a_t
        DNH4g(L)=afac*WQVO(L,K,3)/a_t
        WQPD(L)=DNH4d(L)*WQPD(L)     
        WQPG(L)=DNH4g(L)*WQPG(L)
        WQPC(L)=DNH4c(L)*WQPC(L)
        ENDIF
 !       DIN_mx= MAX (WQVO(L,K,15), 0.0)
 !      DNO3c(L)=1
 !      DNO3d(L)=1
 !      DNO3g(L)=1      
 !      IF(DIN_t.GT.DIN_mx) THEN
 !      afac=DIN_mx/DIN_t
 !      a_t=WQVO(L,K,1)+WQVO(L,K,2)+WQVO(L,K,3)+1.0e-8
 !      DNO3c(L)=afac*WQVO(L,K,1)/a_t
 !      DNO3d(L)=afac*WQVO(L,K,2)/a_t
 !      DNO3g(L)=afac*WQVO(L,K,3)/a_t       
 !      ENDIF
        ENDIF
     &                    ! 
     
        WQA14 = WQA14C*WQANCC1(IMWQZT(L))*WQVO(L,K,1)
     *     + WQA14D*WQANCD1(IMWQZT(L))*WQVO(L,K,2) + 
     &       WQA14G*WQANCG1(IMWQZT(L))*WQVO(L,K,3)

c MRM: When using hourly solar radiation, shut down photosynthesis
c      at night, i.e., when solar radiation is less than 0.001 (05/11/99).
C   IWQSUN = solar radiation switch
C            0=use constant solar radiation (I0) and FD from card C10
C            1=use daily average solar rad. and FD from file SUNDAY.INP
C            2=use hourly solar rad. from ASER.INP file
        if (iwqsun .eq. 2) then
c         if (wqi0 .le. 0.001) then
          if (wqi0 .le. 0.100) then     ! Ji, 9/19/02
            wqpc(L) = 0.0
            wqpd(L) = 0.0
            wqpg(L) = 0.0
            wqpm(L) = 0.0
          end if
        end if
C
C algal basal metabolism & predation
C    
        F_CC=1.0                                             !J.S. 2014 motality
!        IIday=mod(int(TIME),400)
!        IF(IIday.GT.240)then
!        IF(IMWQZT(L).GE.2.and.IMWQZT(L).LE.4)F_CC=20
!        endif
        
        WQBMC(L) = WQBMRC(IMWQZT(L)) * WQTDRC(IWQT(L))       ! Eq. (3-1m), cyanobacteria
        WQPRC(L) = WQPRRC(IMWQZT(L)) * WQTDRC(IWQT(L))       ! Eq. (3-1n), cyanobacteria
c
        WQBMD(L) = WQBMRD(IMWQZT(L)) * WQTDRD(IWQT(L)) *
     &  WQTDGp(iwqt(L),IMWQZT(L))    ! diatom
        WQPRD(L) = WQPRRD(IMWQZT(L)) * WQTDRD(IWQT(L)) * 
     &  WQTDGp(iwqt(L),IMWQZT(L))*F_CC
        WQBMG(L) = WQBMRG(IMWQZT(L)) * WQTDRG(IWQT(L))                    ! green algae
        WQPRG(L) = WQPRRG(IMWQZT(L)) * WQTDRG(IWQT(L))
        IF(ITRYCHLA.EQ.3) THEN      ! JS 2018 add respiration with growth Cerco 2010
        WQBMC(L)=WQBMC(L)+0.12*WQPC(L)
        WQBMD(L)=WQBMD(L)+0.12*WQPD(L)
        WQBMG(L)=WQBMG(L)+0.12*WQPG(L)
        ENDIF
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
        WQOBTOT = WQVO(L,K,1)+WQVO(L,K,2)+WQVO(L,K,3)        ! Total algae
        WQKRPC(L) = (WQKRC + WQKRCALG*WQOBTOT) * WQTDHDR(IWQT(L))
        WQKLPC(L) = (WQKLC + WQKLCALG*WQOBTOT) * WQTDHDR(IWQT(L))
        xmrm = 0.0
        IF (IDNOTRVA.GT.0 .AND. K.EQ.1) THEN
          xmrm = WQKDCALM * WQVo(L,K,IDNOTRVA)
        END IF
        WQKDOC = (WQKDC + WQKDCALG*WQOBTOT + xmrm) * WQTDMNL(IWQT(L))
        O2WQ = MAX(WQVO(L,K,19), 0.0)
        WQTT1 = WQKDOC / (WQKHORDO + O2WQ+ 1.E-18)
        WQKHR(L) = WQTT1 * O2WQ
        WQDENIT(L) = WQTT1 * WQAANOX * RNO3WQ/(WQKHDNN+RNO3WQ+ 1.E-18)
C
C 7-10 phosphorus
C
C    wqKHPC = KHPc = phosphorus half-saturation for cyanobacteria (mg/L), wqwin.inp, C7
C             KHPd = phosphorus half-saturation for algae diatoms (mg/L)
C             KHPg = phosphorus half-saturation for algae greens algae (mg/L)

        WQAPC(L) = 1.0 / (WQCP1PRM1(IMWQZT(L)) +
     &          WQCP2PRM*EXP(-WQCP3PRM*PO4DWQ)) ! Eq.(3-8e), wq3dwc.inp, C21
c     if(L.eq.911.and.k.eq.kc.and.MOD(N,ITIMES).eq.0) then  ! LZ40, surface, check APC ratio, Ji, 10/6/02
c     write(151,1501) N,time,wqapc(L),WQCP1PRM,WQCP2PRM,WQCP3PRM,PO4DWQ
c1501  format(i8,99f12.4)
c     endif
        WQKHP = (WQKHPC+WQKHPD+WQKHPG) / 3.0                         ! Eq.(3-8i)
        WQTT1 = WQKHP / (WQKHP+PO4DWQ+ 1.E-18) * WQOBTOT             !    (3-8f)
        WQKRPP(L) = (WQKRP + WQKRPALG*WQTT1) * WQTDHDR(IWQT(L))      ! Eq.(3-8f)
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
        WQKHN = (WQKHNC+WQKHND+WQKHNG) / 3.0                            ! (3-13e)
        WQTT1 = WQKHN / (WQKHN+RNH4NO3+ 1.E-18) * WQOBTOT
        WQKRPN(L) = (WQKRN + WQKRNALG*WQTT1) * WQTDHDR(IWQT(L))         ! (3-13b)
        WQKLPN(L) = (WQKLN + WQKLNALG*WQTT1) * WQTDHDR(IWQT(L))         ! (3-13c)
        WQKDON(L) = (WQKDN + WQKDNALG*WQTT1) * WQTDMNL(IWQT(L))         ! (3-13d)
C
C14 NH4: WQFTNIT=WQNITM*WQFTNIT
C
        IF (RNH4NO3.EQ.0.0) THEN             ! =NH4+NOx
          WQPNC(L)=0.0
          WQPND(L)=0.0
          WQPNG(L)=0.0
          WQPNM(L)=0.0
        ELSE
          WQTTC = RNH4WQ/(WQKHNC+RNO3WQ+ 1.E-18)
          WQTTD = RNH4WQ/(WQKHND+RNO3WQ+ 1.E-18)
          WQTTG = RNH4WQ/(WQKHNG+RNO3WQ+ 1.E-18)
          WQTTM = RNH4WQ/(WQKHNM+RNO3WQ+ 1.E-18)      ! DIN PREFER HAB
          WQPNC(L) = (RNO3WQ/(WQKHNC+RNH4WQ+ 1.E-18)                    ! (3-13a), cyanobacteria
     $              + WQKHNC/(RNH4NO3+ 1.E-18)) * WQTTC
          WQPND(L) = (RNO3WQ/(WQKHND+RNH4WQ+ 1.E-18)                    !          diatom
     $              + WQKHND/(RNH4NO3+ 1.E-18)) * WQTTD
          WQPNG(L) = (RNO3WQ/(WQKHNG+RNH4WQ+ 1.E-18)                    !          green
     $              + WQKHNG/(RNH4NO3+ 1.E-18)) * WQTTG
          WQPNM(L) = (RNO3WQ/(WQKHNM+RNH4WQ+ 1.E-18)                    !          macroalgae
     $              + WQKHNM/(RNH4NO3+ 1.E-18)) * WQTTM
          WQPDONC(L)=0                    ! WQPDONC(L):  perference for DIN
          WQPNO3C(L)=1-WQPNC(L)           ! WQPNO3C(L):  Add Perferene NO3 (Orginai 1-WQPNC(L)  
        END IF
        IF(WQPMC(4).GT.0.4) THEN          ! IF HAB is simulalated in region 4, Elz. River. i.e. growth rate>0.1
        IF(IMWQZT(L).EQ.4.or.IMWQZT(L).EQ.3) THEN
          IF (RNH4NO3c.EQ.0.0) THEN  
           WQPDONC(L)=0
          ELSE   
          DIN_TT=RNH4WQ+RNO3WQ
          WQPDIN = DIN_TT*RDONWQ/(WQKHNC+RDONWQ)/(WQKHNC+DIN_TT)+
     &        DIN_TT*WQKHNC/(DIN_TT+ RDONWQ)/( WQKHNC+ RDONW)
          WQPDONC(L) = 1-WQPDIN                 ! Perference for DON   HAB
          WQPNC(L) =(1- WQPDONC(L))* WQPNC(L)   ! Perference for NH4   HAB
          WQPNO3C(L)=1-(WQPNC(L)+WQPDONC(L))    ! Perference for NO3   HAB
          ENDIF
        ENDIF
        ENDIF
c          WQNIT(L) = O2WQ * WQTDNIT(IWQT(L)) /                         ! Mike Morton found this bug!, 9/25/00
          WQNIT(L) = O2WQ * RNH4WQ * WQTDNIT(IWQT(L)) /                 ! (3-13g)
     *      ( (WQKHNDO+O2WQ) * (WQKHNN+RNH4WQ) + 1.E-18)
C   
C16-17 silica
C
C  IWQSI   = switch to activate silica state variables (0=off; 1=activated), wq3dwc.inp
        IF (IWQSI.EQ.1) THEN
          IF (IWQSRP.EQ.1) THEN
            WQTTM = WQKSAP*WQTAMP(L,K)
            WQN17(L) = - WQWSSET(L,1) * WQTTM / (1.0+WQTTM)
            IF (K.NE.KC) THEN
              WQTTM = WQKSAP*WQTAMP(L,K+1)
              WQT17(L) = WQWSSET(L,2) * WQTTM / (1.0+WQTTM)
            END IF
           ELSE IF (IWQSRP.EQ.2) THEN                       ! SEDT
            WQTTS = WQKSAP*SEDTWQ(L,K)
            WQN17(L) = - WQSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)      ! (3-15), settling term, but what about d/dz? Ji, 9/18/02
            IF (K.NE.KC) THEN
              WQTTS = WQKSAP*SEDTWQ(L,K+1)
              WQT17(L) = WQSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)      ! (3-15), settling term
            END IF
           ELSE                                             ! SEDT
            WQN17(L) = 0.0
            WQT17(L) = 0.0
            END IF
        END IF
C
c Good approach, I can do the same for other viriables, Ji, 9/18/02
c
C 04/29/99 MRM:
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
c MRM: 06/20/98
c In the following equation, salinity must be in mg/L, hence, swq(L)
c is multiplied by 1000.
c        tval1 = twq(L)
c        tval2 = tval1*tval1
c        wqdos = 14.5532 -  0.38217*tval1 + 5.4258e-3*tval2 -
c     *     (swq(L)*1000.0/1.80655) * (1.665e-4 - 5.866e-6*tval1 +
c     *     9.796e-8*tval2)
CTt          WQDOS = 1.45532E1 + (5.4258E-3*TWQ(L)- 3.8217E-1)*TWQ(L) -
CTt  *  (SWQ(L)/1.80655)
CTt  $           * ( 1.665E-4 + (9.796E-8*TWQ(L) - 5.866E-6)*TWQ(L) )
CTt  uncomment and change line below to activate wind effect on reareation
c MRM: 06/20/98
c Do not allow wind speeds above 11 m/sec in the following equation:
        windrea = windstwq(L)
        if (windstwq(L) .gt. 11.0) windrea = 11.0
        WQWREA = 0.728*SQRT(windrea) + (0.0372*windrea - 0.317)*windrea
c        WQWREA = 0.728*SQRT(WINDST(L))
c     $           +(0.0372*WINDST(L) - 0.317)*WINDST(L)
CTt       WQWREA = 0.0

C J.S. 2/21/99 Correct surface velocity
c          WQP19(L) = - (WQKRO*SQRT(U(L,KC)**2+V(L,KC)**2) + WQWREA)
c          WQP19(L) = - (WQKRO*SQRT(0.5) + WQWREA) * DZWQ(L)
c     *      * WQTDKR(IWQT(L))
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
c           xmrm = SQRT(U(L,K)*U(L,K) + V(L,K)*V(L,K))
            umrm = UWQS(L) !max(U(L,K), U(L+1,K))
            vmrm = VWQS(L)  !max(V(L,K), V(LNC(L),K))
            xmrm = SQRT(umrm*umrm + vmrm*vmrm)
            WQVREA = WQKRO * xmrm**0.5 / HPWQ(L)**1.5
          END IF
c
c MRM 04/12/1999  Owens and Gibbs (1964) reaeration equation:
c    WQKRO = 5.32 typically
c
          IF (IWQKA .EQ. 3) THEN
c           xmrm = SQRT(U(L,K)*U(L,K) + V(L,K)*V(L,K))
            umrm = UWQS(L) ! max(U(L,K), U(L+1,K))
            vmrm = VWQS(L) ! max(V(L,K), V(LNC(L),K))
            xmrm = SQRT(umrm*umrm + vmrm*vmrm)
            WQVREA = WQKRO * xmrm**0.67 / HPWQ(L)**1.85
          END IF
c
c -- from codemike.011701, other options  ----
c
c Modified Owens and Gibbs reaeration equation:
c Note: Normalized to a depth of 1.0 ft, i.e., this equation gives the
c       same reaeration as Owens & Gibbs at 1.0 ft depth; at higher
c       depths it gives larger reaeration than Owens & Gibbs.
c WQKRO = 5.32 typically
c
c             IF (IWQKA(IZ) .EQ. 4) THEN
c                xmrm = xmrm * SQRT(HP(L))/0.3048
c                xmrm = SQRT(U(L,K)*U(L,K) + V(L,K)*V(L,K))
c               umrm = max(U(L,K), U(L+1,K))
c               vmrm = max(V(L,K), V(LNC(L),K))
c               xmrm = SQRT(umrm*umrm + vmrm*vmrm)
c               ymrm = HP(L)*3.0*(1.0 - HP(L)/(HP(L)+0.1524))
c               WQVREA = WQKRO(IZ) * xmrm**0.67 / ymrm**1.85
c             END IF
c
c Melching & Flores (1999), ASCE J. Environ. Enrg. 125(5):407-414
c Pool and Riffle type streams
c
c             IF (IWQKA(IZ) .EQ. 5) THEN
c               umrm = max(U(L,K), U(L+1,K))
c               vmrm = max(V(L,K), V(LNC(L),K))
c               xvel = SQRT(umrm*umrm + vmrm*vmrm)
c               UTMP = MAX( UHDYE(L), UHDYE(L+1) )
c               VTMP = MAX( VHDXE(L), VHDXE(LNC(L)) )
c               xflo = SQRT(utmp*utmp + vtmp*vtmp)
c               xvs = xvel*slope(iz)
c               if (xflo .lt. 0.556) then
c                 WQVREA = 517.0 * xvs**0.524 * xflo**(-0.242)
c               else
c                 WQVREA = 596.0 * xvs**0.528 * xflo**(-0.136)
c               end if
c             END IF
c
c Melching & Flores (1999), ASCE J. Environ. Enrg. 125(5):407-414
c Channel control type streams
c
c             IF (IWQKA(IZ) .EQ. 6) THEN
c               umrm = max(U(L,K), U(L+1,K))
c               vmrm = max(V(L,K), V(LNC(L),K))
c               xvel = SQRT(umrm*umrm + vmrm*vmrm)
c               UTMP = MAX( UHDYE(L), UHDYE(L+1) )
c               VTMP = MAX( VHDXE(L), VHDXE(LNC(L)) )
c               xflo = SQRT(utmp*utmp + vtmp*vtmp)
c               xvs = xvel*slope(iz)
c               if (xflo .lt. 0.556) then
c                 WQVREA = 88.0 * xvs**0.313 * HP(L)**(-0.353)
c               else
c                 xwidth = min (dxp(L), dyp(L))
c                 WQVREA = 142.0 * xvs**0.333 * HP(L)**(-0.66) *
c    +                     xwidth**(-0.243)
c               end if
c             END IF
c
c
c
c now combine reaeration due to water velocity and wind stress:
c
!          WQWREA = 0.0
          WQVREA = WQVREA * REAC(IWQZMAP(L,K))
          WQWREA = WQWREA * REAC(IWQZMAP(L,K))
          WQP19(L) = - (WQVREA + WQWREA) * DZWQ(L)* WQTDKR(IWQT(L))
          WQKRDOS(L) = - WQP19(L)*WQDOS(L)
          xDOacoef(L,K) = xDOacoef(L,K) +  
     &          (WQVREA + WQWREA) *  WQTDKR(IWQT(L))
        ELSE
          WQP19(L) = 0.0
        END IF
C
C20 TAM: WQTDTAM=WQKHBMF*BFTAM*EXP()
c
c  => TAM is NOT used when IWQSRP=2, which is the case of using TSS to carry the particulates, Ji, 9/18/02
C
        IF (IWQSRP.EQ.1) THEN
          WQR20(L) = (WQWDSL(L,K,20)+WQWPSL(L,K,20))*VOLWQ(L)
     *      + (WQVO(L,K,20) - WQTAMP(L,K)) * WQWSSET(L,1)
c mrm  add wet atmospheric deposition:
          if (k.eq.kc) then
            WQR20(L) = WQR20(L) + WQATML(L,kc,20) * VOLWQ(L)
          end if
          IF (K.EQ.1) WQR20(L) = WQR20(L)
     *      + WQTDTAM(IWQT(L))*DZWQ(L)/(WQKHBMF+O2WQ+ 1.E-18)
          IF (K.NE.KC) WQR20(L) = WQR20(L)
     *      + (WQVO(L,K+1,20) - WQTAMP(L,K+1)) * WQWSSET(L,2)
        END IF
C
       END DO                ! #100
C
C trapezoidal solution of kinetic eqs: After computing new values, store
C WQVO+WQV into WQVO(L,K,NWQV)
C
C Modified by J.S.(5/5/98) for Macalgae
C
      IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
       DO L=LF,LL
        WQA1C = (WQPM(L) - WQBMM(L) - WQPRM(L)-WQWSM*DZWQ(L))*DTWQO2
        WQVA1C = 1.0 / (1.0 - WQA1C)
        WQV(L,K,IDNOTRVA)=(WQVO(L,K,IDNOTRVA)+WQA1C*WQVO(L,K,IDNOTRVA))
     *                   *WQVA1C*SMAC(L)
        WQV(L,K,IDNOTRVA) = max(WQV(L,K,IDNOTRVA),WQMCMIN)*SMAC(L)
        WQVO(L,K,IDNOTRVA) = WQVO(L,K,IDNOTRVA)+WQV(L,K,IDNOTRVA)
       END DO
      END IF
      if(mod(N,480).eq.0)then
      write(971,*)WQPMM(6),WQPM(6),WQI0BOT(6),WQBMM(6), WQPRM(6),WQMCMIN
      endif
C J.S
C
c FEC (=21) is NOT calculated here, Ji, 9/18/02
c
c-----------------------------------
c--- loadings & integration --------
c
C1 Bc: WQT1C=WQBCSET(L,2)
C
       DO L=LF,LL                              ! settling to lower layer
        WQA1C = (WQPC(L) - WQBMC(L) - WQPRC(L) - WQBCSET(L,1))*DTWQO2   ! (3-1), the first 4 terms
        WQKK(L) = 1.0 / (1.0 - WQA1C)
        WQR1C = (WQWDSL(L,K,1) + WQWPSL(L,K,1)) * VOLWQ(L)
              ! Dry deposition
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQR1C = WQR1C + WQATML(L,kc,1) * VOLWQ(L)
        end if
        WQRR(L) = WQVO(L,K,1) + DTWQ*WQR1C + WQA1C*WQVO(L,K,1)
       END DO
C
       IF (K.NE.KC) THEN
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)        
         DO L=LF,LL           ! settling from upper layer
          WQRR(L) = WQRR(L) + DTWQO2*WQBCSET(L,2)*WQVO(L,K+1,1)         ! (3-1)
         END DO
       END IF
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,WQA2D) 

       DO L=LF,LL
CTt update modified next line added line below commented out
       WQV(L,K,1)=SCBWQ(L)*( WQRR(L)*WQKK(L) )+(1.-SCBWQ(L))*WQVO(L,K,1)   ! SCBWQ=0=boundary cells, =1=interia cells, for concentrations, Ji, 9/18/02
       WQV(L,K,1)=max(WQV(L,K,1),WQMCMIN)  ! !J.S. mini. algae 2018
ctt        WQV(L,K,1) = WQRR(L)*WQKK(L)                                 ! => Bounadry cells use old WQV0, which should be specified by BC (?)
        WQVO(L,K,1) = WQVO(L,K,1)+WQV(L,K,1)
C
C2 Bd: WQT2D=WQBDSET(L,2)
C
        WQA2D = (WQPD(L) - WQBMD(L) - WQPRD(L) - WQBDSET(L,1))*DTWQO2
        WQKK(L) = 1.0 / (1.0 - WQA2D)
        WQR2D = (WQWDSL(L,K,2) + WQWPSL(L,K,2)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQR2D = WQR2D + WQATML(L,kc,2) * VOLWQ(L)
        end if
        WQRR(L) = WQVO(L,K,2) + DTWQ*WQR2D + WQA2D*WQVO(L,K,2)
                          
       END DO
C
       IF (K.NE.KC) THEN
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
         DO L=LF,LL
          WQRR(L) = WQRR(L) + DTWQO2*WQBDSET(L,2)*WQVO(L,K+1,2)
         END DO
       END IF
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,WQA3G,WQR3G)

       DO L=LF,LL
CTt update modified next line added line below commented out
      WQV(L,K,2)=SCBWQ(L)*( WQRR(L)*WQKK(L) )+(1.-SCBWQ(L))*WQVO(L,K,2)
ctt        WQV(L,K,2) = WQRR(L)*WQKK(L)
        WQV(L,K,2)=max(WQV(L,K,2),WQMCMIN)             !J.S. mini. algae 2018
        WQVO(L,K,2) = WQVO(L,K,2)+WQV(L,K,2)
C
C3 Bg: WQT3D=WQBGSET(L,2)
C
        WQA3G = (WQPG(L) - WQBMG(L) - WQPRG(L) - WQBGSET(L,1))*DTWQO2
        WQKK(L) = 1.0 / (1.0 - WQA3G)
        WQR3G = (WQWDSL(L,K,3) + WQWPSL(L,K,3)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQR3G = WQR3G + WQATML(L,kc,3) * VOLWQ(L)
        end if
        WQRR(L) = WQVO(L,K,3) + DTWQ*WQR3G + WQA3G*WQVO(L,K,3)
       END DO
C
       IF (K.NE.KC) THEN
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
         DO L=LF,LL
          WQRR(L) = WQRR(L) + DTWQO2*WQBGSET(L,2)*WQVO(L,K+1,3)
         END DO
       END IF
C  
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)     
       DO L=LF,LL
       WQV(L,K,3)=SCBWQ(L)*( WQRR(L)*WQKK(L) )+(1.-SCBWQ(L))*WQVO(L,K,3)
        WQV(L,K,3)=max(WQV(L,K,3),WQMCMIN)          !JS 2018 min. algae
       WQVO(L,K,3) = WQVO(L,K,3)+WQV(L,K,3)
       ENDDO
!
! Record algae con. as Chl-a at n+1 and Production, J.S. 8/16/08
!
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
       DO L=LF,LL
        xAlg(L)=xAlg(L)+0.5*( WQVO(L,K,1)*WQCHLC(IWQZMAP(L,KC))          !J.S. 11/27/2014
     &   + WQCHLD(IWQZMAP(L,KC))*WQVO(L,K,2)+
     &  WQCHLG(IWQZMAP(L,KC))*WQVO(L,K,3))/(Real(KC)) 
 !       xPro(L)=xPro(L)+DTWQO2*(
 !    &  WQVO(L,K,1)*WQPC(L)+WQVO(L,K,2)*WQPG(L)+WQVO(L,K,3)*WQPG(L) )  
 
 !       xPro(L)=xPro(L)+DTWQO2*(
 !    &  WQVO(L,K,1)*WQPC(L)+WQVO(L,K,2)*WQPG(L)+WQVO(L,K,3)*WQPG(L) )
 !    *  *(DZCWQ(K)*HPWQ(L))  
     
 !    & ( WQVO(L,K,1)*( (1+ DTWQO2*WQPC(L))/(1- DTWQO2*WQPC(L))-1 ) +
 !    &       WQVO(L,K,2)*( (1+ DTWQO2*WQPD(L))/(1- DTWQO2*WQPD(L))-1 ) +  
 !    &       *( (1+ DTWQO2*)/(1- DTWQO2*WQPG(L))-1 ) )
 
 ! save production
         xPro(L)=xPro(L)+
     &  ( WQV(L,K,1)*MAX(exp(2*DTWQO2*WQPC(L))-1,0.0 ) + 
     &    WQV(L,K,2)*MAX(exp(2*DTWQO2*WQPD(L))-1,0.0 ) +  
     &    WQV(L,K,3)*MAX(exp(2*DTWQO2*WQPG(L))-1,0.0 ) )
     &  *(DZCWQ(K)*HPWQ(L))
 
       ENDDO
!

C
C4 RPOC: WQA4=WQRR(L),WQT4=WQRPSET(L,2)
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,WQB4,WQA4,WQR4)
       DO L=LF,LL
        WQB4 = - (WQKRPC(L)+WQRPSET(L,1))
        WQKK(L) = 1.0 / (1.0 - DTWQO2*WQB4)
        WQA4 = WQFCRP * (WQPRC(L)*WQVO(L,K,1)
     *     + WQPRD(L)*WQVO(L,K,2) + WQPRG(L)*WQVO(L,K,3))
C
C add macalgal source     J.S.
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQA4 = WQA4+WQFCRPM*WQPRM(L)*WQVO(L,K,IDNOTRVA)
      END IF
C J.S.
        WQR4 = (WQWDSL(L,K,4) + WQWPSL(L,K,4)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQR4 = WQR4 + WQATML(L,kc,4) * VOLWQ(L)
        end if
        WQRR(L) = WQVO(L,K,4) + DTWQ*WQR4 + DTWQO2*( WQA4
     *    + WQB4*WQVO(L,K,4) )
       END DO
C
       IF (K.NE.KC) THEN
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
         DO L=LF,LL
          WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQVO(L,K+1,4)
         END DO
       END IF
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,WQC5,WQA5,WQR5)
       DO L=LF,LL
CTt update modified next line added line below commented out
      WQV(L,K,4)=SCBWQ(L)*( WQRR(L)*WQKK(L) )+(1.-SCBWQ(L))*WQVO(L,K,4)
ctt        WQV(L,K,4) = WQRR(L)*WQKK(L)
        WQVO(L,K,4) = WQVO(L,K,4)+WQV(L,K,4)
C
C5 LPOC: WQT5=WQLPSET(L,2)
C
        IF(K>1) THEN                       ! J.S. 2018 mofilfy resuspension
        WQC5 = - (WQKLPC(L)+WQLPSET(L,1))
        ELSE
        WQC5 = - (WQKLPC(L)+WQLPSET(L,1)*0.3)        
        ENDIF
        WQKK(L) = 1.0 / (1.0 - DTWQO2*WQC5)
        WQA5 = WQFCLP * (WQPRC(L)*WQVO(L,K,1)
     *    + WQPRD(L)*WQVO(L,K,2) + WQPRG(L)*WQVO(L,K,3))
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQA5 =WQA5 + WQFCLPM * WQPRM(L)*WQVO(L,K,IDNOTRVA)
      END IF
C J.S.
        WQR5 = (WQWDSL(L,K,5) + WQWPSL(L,K,5)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQR5 = WQR5 + WQATML(L,kc,5) * VOLWQ(L)
        end if
        WQRR(L) = WQVO(L,K,5) + DTWQ*WQR5 + DTWQO2*( WQA5
     *    + WQC5*WQVO(L,K,5) )
       END DO
C
       IF (K.NE.KC) THEN
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)       
         DO L=LF,LL
          WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,5)
         END DO
       END IF
C
       DO L=LF,LL
CTt update modified next line added line below commented out
       WQV(L,K,5)=SCBWQ(L)*( WQRR(L)*WQKK(L) )+(1.-SCBWQ(L))*WQVO(L,K,5)
ctt        WQV(L,K,5) = WQRR(L)*WQKK(L)
        WQVO(L,K,5) = WQVO(L,K,5)+WQV(L,K,5)
C
C6 DOC: CFCDxWQ=1-WQFCDx, WQB6=WQKRPC(L),WQC6=WQKLPC(L)
C
        WQD6 = - (WQKHR(L)+WQDENIT(L))
        WQKK(L) = 1.0 / (1.0 - DTWQO2*WQD6)
        O2WQ = MAX(WQVO(L,K,19), 0.0)
        WQA6C=(WQFCDC+CFCDCWQ*WQKHRC/(WQKHRC+O2WQ+ 1.E-18))*WQBMC(L)
        WQA6D=(WQFCDD+CFCDDWQ*WQKHRD/(WQKHRD+O2WQ+ 1.E-18))*WQBMD(L)
        WQA6G=(WQFCDG+CFCDGWQ*WQKHRG/(WQKHRG+O2WQ+ 1.E-18))*WQBMG(L)
        WQA6 = ( WQA6C + WQFCDP*WQPRC(L) )*WQVO(L,K,1)
     *         + ( WQA6D + WQFCDP*WQPRD(L) )*WQVO(L,K,2)
     *         + ( WQA6G + WQFCDP*WQPRG(L) )*WQVO(L,K,3)
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
       WQA6M=(WQFCDM+(1-WQFCDM)*WQKHRM/(WQKHRM+O2WQ+ 1.E-18))*WQBMM(L)
       WQA6 =WQA6+ (WQA6M+ WQFCDPM*WQPRM(L))*WQVO(L,K,IDNOTRVA)
      END IF
C J.S.
C
        WQR6 = (WQWDSL(L,K,6) + WQWPSL(L,K,6)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQR6 = WQR6 + WQATML(L,kc,6) * VOLWQ(L)
        end if
        WQRR(L) = WQVO(L,K,6) + DTWQ*WQR6 + DTWQO2*( WQA6 + WQKRPC(L)*
     *    WQVO(L,K,4) + WQKLPC(L)*WQVO(L,K,5) + WQD6*WQVO(L,K,6) )
CTt update modified next line added line below commented out
        WQV(L,K,6)=SCBWQ(L)*( WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,6)
ctt        WQV(L,K,6) = WQRR(L)*WQKK(L)
        WQVO(L,K,6) = WQVO(L,K,6)+WQV(L,K,6)
C
C7 RPOP: WQT7=WQT4=WQRPSET(L,2)
C
        WQE7 = - (WQKRPP(L)+WQRPSET(L,1))
        WQKK(L) = 1.0 / (1.0 - DTWQO2*WQE7)
        WQA7C = (WQFPRC*WQBMC(L) + WQFPRP*WQPRC(L)) * WQVO(L,K,1)
        WQA7D = (WQFPRD*WQBMD(L) + WQFPRP*WQPRD(L)) * WQVO(L,K,2)
        WQA7G = (WQFPRG*WQBMG(L) + WQFPRP*WQPRG(L)) * WQVO(L,K,3)
        WQA7 = (WQA7C+WQA7D+WQA7G) * WQAPC(L)
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQA7 = WQA7 + (WQFPRM*WQBMM(L) + WQFPRPM*WQPRM(L))
     *         * WQVO(L,K,IDNOTRVA)* WQAPC(L)*WQAPCM

      END IF
C J.S.
        WQR7 = (WQWDSL(L,K,7) + WQWPSL(L,K,7)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQR7 = WQR7 + WQATML(L,kc,7) * VOLWQ(L)
        end if
        WQRR(L) = WQVO(L,K,7) + DTWQ*WQR7 + DTWQO2*( WQA7
     *     + WQE7*WQVO(L,K,7) )
       END DO
C
       IF (K.NE.KC) THEN
         DO L=LF,LL
          WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQVO(L,K+1,7)
         END DO
       END IF
C
       DO L=LF,LL
CTt update modified next line added line below commented out
        WQV(L,K,7)=SCBWQ(L)*( WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,7)
ctt        WQV(L,K,7) = WQRR(L)*WQKK(L)
        WQVO(L,K,7) = WQVO(L,K,7)+WQV(L,K,7)
C
C8 LPOP: WQT8=WQT5=WQLPSET(L,2)
C
        WQF8 = - (WQKLPP(L)+WQLPSET(L,1))
        WQKK(L) = 1.0 / (1.0 - DTWQO2*WQF8)
        WQA8C = (WQFPLC*WQBMC(L) + WQFPLP*WQPRC(L)) * WQVO(L,K,1)
        WQA8D = (WQFPLD*WQBMD(L) + WQFPLP*WQPRD(L)) * WQVO(L,K,2)
        WQA8G = (WQFPLG*WQBMG(L) + WQFPLP*WQPRG(L)) * WQVO(L,K,3)
        WQA8 = (WQA8C+WQA8D+WQA8G) * WQAPC(L)
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
       WQA8 = WQA8 +     (WQFPLM*WQBMM(L) + WQFPLPM*WQPRM(L))
     *         * WQVO(L,K,IDNOTRVA)* WQAPC(L)*WQAPCM
       END IF
C J.S.
        WQR8 = (WQWDSL(L,K,8) + WQWPSL(L,K,8)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQR8 = WQR8 + WQATML(L,kc,8) * VOLWQ(L)
        end if
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
CTt update modified next line added line below commented out
        WQV(L,K,8)=SCBWQ(L)*( WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,8)
ctt        WQV(L,K,8) = WQRR(L)*WQKK(L)
        WQVO(L,K,8) = WQVO(L,K,8)+WQV(L,K,8)
C
C9 DOP: WQE9=WQKRPP(L),WQF9=WQKLPP(L),WQG9=-WQKDOP(L)
C
        WQKK(L) = 1.0 / (1.0 + DTWQO2*WQKDOP(L))
        WQA9C = (WQFPDC*WQBMC(L) + WQFPDP*WQPRC(L)) * WQVO(L,K,1)
        WQA9D = (WQFPDD*WQBMD(L) + WQFPDP*WQPRD(L)) * WQVO(L,K,2)
        WQA9G = (WQFPDG*WQBMG(L) + WQFPDP*WQPRG(L)) * WQVO(L,K,3)
        WQA9 = (WQA9C+WQA9D+WQA9G) * WQAPC(L)
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
       WQA9 = WQA9 + (WQFPDM*WQBMM(L) + WQFPDPM*WQPRM(L))
     *       * WQVO(L,K,IDNOTRVA) * WQAPC(L)*WQAPCM
       END IF
C J.S.
        WQR9 = (WQWDSL(L,K,9) + WQWPSL(L,K,9)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQR9 = WQR9 + WQATML(L,kc,9) * VOLWQ(L)
        end if
        WQRR(L) = WQVO(L,K,9) + DTWQ*WQR9 + DTWQO2*(WQA9 + WQKRPP(L)*
     *   WQVO(L,K,7) + WQKLPP(L)*WQVO(L,K,8) - WQKDOP(L)*WQVO(L,K,9))
CTt update modified next line added line below commented out
        WQV(L,K,9)=SCBWQ(L)*( WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,9)
ctt        WQV(L,K,9) = WQRR(L)*WQKK(L)
        WQVO(L,K,9) = WQVO(L,K,9)+WQV(L,K,9)
C
C10 PO4t: WQG10=WQKDOP(L),WQA10=WQKK(L)
C
        WQA10C = (WQFPIC*WQBMC(L)+WQFPIP*WQPRC(L)-WQPC(L))*WQVO(L,K,1)
        WQA10D = (WQFPID*WQBMD(L)+WQFPIP*WQPRD(L)-WQPD(L))*WQVO(L,K,2)
        WQA10G = (WQFPIG*WQBMG(L)+WQFPIP*WQPRG(L)-WQPG(L))*WQVO(L,K,3)
        WQKK(L) = (WQA10C+WQA10D+WQA10G) * WQAPC(L)
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQKK(L) =WQKK(L)+(WQFPIM*WQBMM(L)+WQFPIP*WQPRM(L)-WQPM(L))
     *           *WQVO(L,K,IDNOTRVA) * WQAPC(L)*WQAPCM
       END IF
C J.S.
        WQRR(L) = (WQWDSL(L,K,10)+WQWPSL(L,K,10)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQRR(L) = WQRR(L) + WQATML(L,kc,10) * VOLWQ(L)
        end if
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
CTt update modified next line added line below commented out
      WQV(L,K,10)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,10)
ctt         WQV(L,K,10) = WQRR(L)*WQKK(L)
         WQVO(L,K,10) = WQVO(L,K,10)+WQV(L,K,10)
C
C11 RPON: WQT11=WQT4=WQRPSET(L,2)
C
      WQI11 = - (WQKRPN(L)+WQRPSET(L,1))
      WQKK(L) = 1.0 / (1.0 - DTWQO2*WQI11)
      WQA11C = (WQFNRC*WQBMC(L)+WQFNRP*WQPRC(L))*
     &   WQANCC1(IMWQZT(L))*WQVO(L,K,1)
         WQA11D = (WQFNRD*WQBMD(L)+WQFNRP*WQPRD(L))*
     &   WQANCD1(IMWQZT(L))*WQVO(L,K,2)
         WQA11G = (WQFNRG*WQBMG(L)+WQFNRP*WQPRG(L))*
     &   WQANCG1(IMWQZT(L))*WQVO(L,K,3)
         WQA11 = WQA11C+WQA11D+WQA11G
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
      WQA11 =WQA11 +     (WQFNRM*WQBMM(L)+WQFNRPM*WQPRM(L))
     *        *WQANCM*WQVO(L,K,IDNOTRVA)
      END IF
C J.S.
         WQR11 = (WQWDSL(L,K,11)+WQWPSL(L,K,11)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQR11 = WQR11 + WQATML(L,kc,11) * VOLWQ(L)
        end if
         WQRR(L) = WQVO(L,K,11) + DTWQ*WQR11 + DTWQO2*( WQA11
     *     + WQI11*WQVO(L,K,11) )
       END DO
C
       IF (K.NE.KC) THEN
         DO L=LF,LL
           WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQVO(L,K+1,11)
         END DO
       END IF
C
       DO L=LF,LL
CTt update modified next line added line below commented out
      WQV(L,K,11)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,11)
ctt         WQV(L,K,11) = WQRR(L)*WQKK(L)
         WQVO(L,K,11) = WQVO(L,K,11)+WQV(L,K,11)
C
C12 LPON: WQT12=WQT5=WQLPSET(L,2)
C
         WQJ12 = - (WQKLPN(L)+WQLPSET(L,1))
         WQKK(L) = 1.0 / (1.0 - DTWQO2*WQJ12)
         WQA12C = (WQFNLC*WQBMC(L)+WQFNLP*WQPRC(L))*
     &            WQANCC1(IMWQZT(L))*WQVO(L,K,1)
         WQA12D = (WQFNLD*WQBMD(L)+WQFNLP*WQPRD(L))*
     &            WQANCD1(IMWQZT(L))*WQVO(L,K,2)
         WQA12G = (WQFNLG*WQBMG(L)+WQFNLP*WQPRG(L))*
     &            WQANCG1(IMWQZT(L))*WQVO(L,K,3)
         WQA12 = WQA12C+WQA12D+WQA12G
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
       WQA12 =WQA12 +(WQFNLM*WQBMM(L)+WQFNLPM*WQPRM(L))
     *      *WQANCM*WQVO(L,K,IDNOTRVA)
       END IF
C J.S.
         WQR12 = (WQWDSL(L,K,12)+WQWPSL(L,K,12)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQR12 = WQR12 + WQATML(L,kc,12) * VOLWQ(L)
        end if
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
CTt update modified next line added line below commented out
       WQV(L,K,12)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,12)
ctt        WQV(L,K,12) = WQRR(L)*WQKK(L)
        WQVO(L,K,12) = WQVO(L,K,12)+WQV(L,K,12)
C
C13 DON: WQI13=WQKRPN(L),WQJ13=WQKLPN(L),WQK13=-WQKDON(L)
C
        WQKK(L) = 1.0 / (1.0 + DTWQO2*WQKDON(L))
        WQA13C = (WQFNDC*WQBMC(L)+WQFNDP*WQPRC(L)) 
     &   *WQANCC1(IMWQZT(L))*WQVO(L,K,1)
        WQA13D = (WQFNDD*WQBMD(L)+WQFNDP*WQPRD(L))
     &   *WQANCD1(IMWQZT(L))*WQVO(L,K,2)
        WQA13G = (WQFNDG*WQBMG(L)+WQFNDP*WQPRG(L))
     &   *WQANCG1(IMWQZT(L))*WQVO(L,K,3)
        WQA13 = WQA13C+WQA13D+WQA13G
        IF(IMWQZT(L).EQ.3.or.IMWQZT(L).EQ.4) THEN
         WQA13 =WQA13 - 
     &   (WQPDONC(L)*WQPC(L))*WQANCC1(IMWQZT(L))*WQVO(L,K,1)    ! HAB uptake DON by remove DON
        ENDIF
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQA13 =WQA13 + (WQFNDM*WQBMM(L)+WQFNDPM*WQPRM(L))
     *         *WQANCM*WQVO(L,K,IDNOTRVA)
       END IF
C J.S.
        WQR13 = (WQWDSL(L,K,13) + WQWPSL(L,K,13)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQR13 = WQR13 + WQATML(L,kc,13) * VOLWQ(L)
        end if
        WQRR(L) = WQVO(L,K,13) + DTWQ*WQR13 + DTWQO2*( WQA13
     *     + WQKRPN(L)*WQVO(L,K,11) + WQKLPN(L)*WQVO(L,K,12)
     *     - WQKDON(L)*WQVO(L,K,13) )
CTt update modified next line added line below commented out
       WQV(L,K,13)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,13)
ctt        WQV(L,K,13) = WQRR(L)*WQKK(L)
        WQVO(L,K,13) = WQVO(L,K,13)+WQV(L,K,13)
C
C14 NH4: WQL14=-WQNIT(L),WQK14=WQKDON(L),WQR14=WQRR(L)
C

        WQRR(L) = (WQWDSL(L,K,14)+WQWPSL(L,K,14)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQRR(L) = WQRR(L) + WQATML(L,kc,14) * VOLWQ(L)
        end if
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
        WQA14C = WQFNIC*WQBMC(L) + WQFNIP*WQPRC(L) - WQPNC(L)*WQPC(L)
        WQA14D = WQFNID*WQBMD(L) + WQFNIP*WQPRD(L) - WQPND(L)*WQPD(L)
        WQA14G = WQFNIG*WQBMG(L) + WQFNIP*WQPRG(L) - WQPNG(L)*WQPG(L)
        WQA14 = WQA14C*WQANCC1(IMWQZT(L))*WQVO(L,K,1)
     *     + WQA14D*WQANCD1(IMWQZT(L))*WQVO(L,K,2) + 
     &       WQA14G*WQANCG1(IMWQZT(L))*WQVO(L,K,3)
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

CTt update modified next line added line below commented out
       WQV(L,K,14)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,14)
ctt        WQV(L,K,14) = WQRR(L)*WQKK(L)
        WQVO(L,K,14) = WQVO(L,K,14)+WQV(L,K,14)
C
C15 NO3: WQKK(L)=1,WQR15=WQRR(L),WQD15=-WQANDC*WQDENIT(L),WQL15=WQNIT(L)
C
         WQRR(L) = (WQWDSL(L,K,15)+WQWPSL(L,K,15)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQRR(L) = WQRR(L) + WQATML(L,kc,15) * VOLWQ(L)
        end if
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
!      WQA15C = (WQPNC(L)-1.0)*WQPC(L) * WQANCC1(IMWQZT(L)) * WQVO(L,K,1)
      WQA15C = -(WQPNO3C(L)*WQPC(L)) * WQANCC1(IMWQZT(L)) * WQVO(L,K,1)           ! HAB remove DON
      WQA15D = (WQPND(L)-1.0)*WQPD(L) * WQANCD1(IMWQZT(L)) * WQVO(L,K,2)
      WQA15G = (WQPNG(L)-1.0)*WQPG(L) * WQANCG1(IMWQZT(L)) * WQVO(L,K,3)
      WQA15 = WQA15C+WQA15D+WQA15G
C
C add macalgal source     J.S.
C
       IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
        WQA15 =WQA15 + (WQPNM(L)-1.0)*WQPM(L)*WQANCM
     *         *WQVO(L,K,IDNOTRVA)
       END IF
C J.S.
CTt update modified next line added line below commented out
        WQV(L,K,15)=SCBWQ(L)*( WQVO(L,K,15) + DTWQ*WQRR(L)
     *     + DTWQO2*( WQA15
     *     - WQANDC*WQDENIT(L)*WQVO(L,K,6) + WQNIT(L)*WQVO(L,K,14) ) )
     *     +(1.-SCBWQ(L))*WQVO(L,K,15)
ctt         WQV(L,K,15) = WQVO(L,K,15) + DTWQ*WQRR(L) + DTWQO2*( WQA15
ctt     *     - WQANDC*WQDENIT(L)*WQVO(L,K,6) + WQNIT(L)*WQVO(L,K,14) )
         WQVO(L,K,15) = WQVO(L,K,15)+WQV(L,K,15)
       END DO
C
C16 SU: WQT16=WQBDSET(L,2)
C
C  IWQSI   = switch to activate silica state variables (0=off; 1=activated), wq3dwc.inp
       IF (IWQSI.EQ.1) THEN
C
         DO L=LF,LL
           WQM16 = - (WQKSUA(IWQT(L)) + WQBDSET(L,1))
           WQKK(L) = 1.0 / (1.0 - DTWQO2*WQM16)
           WQA16D = (WQFSPD*WQBMD(L) + WQFSPP*WQPRD(L)) * WQASCD
     *       * WQVO(L,K,2)
           WQR16 = (WQWDSL(L,K,16)+WQWPSL(L,K,16)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQR16 = WQR16 + WQATML(L,kc,16) * VOLWQ(L)
        end if
           WQRR(L) = WQVO(L,K,16) + DTWQ*WQR16 + DTWQO2*( WQA16D
     *       + WQM16*WQVO(L,K,16) )
         END DO
C
         IF (K.NE.KC) THEN
           DO L=LF,LL
             WQRR(L) = WQRR(L) + DTWQO2*WQBDSET(L,2)*WQVO(L,K+1,16)
           END DO
         END IF
C
         DO L=LF,LL
CTt update modified next line added line below commented out
       WQV(L,K,16)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,16)
ctt           WQV(L,K,16) = WQRR(L)*WQKK(L)
           WQVO(L,K,16) = WQVO(L,K,16)+WQV(L,K,16)
C
C17 SA: WQM17=WQKSUA(IWQT(L)),WQA17D=WQKK(L)
C
           WQKK(L) = (WQFSID*WQBMD(L) + WQFSIP*WQPRD(L) - WQPD(L))
     *       * WQASCD * WQVO(L,K,2)
           WQRR(L) = (WQWDSL(L,K,17)+WQWPSL(L,K,17)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
           if (k.eq.kc) then
             WQRR(L) = WQRR(L) + WQATML(L,kc,17) * VOLWQ(L)
           end if
         END DO
C
         IF (K.EQ.1) THEN
           DO L=LF,LL
             WQRR(L) = WQRR(L) + WQBFSAD(L)*DZWQ(L)
           END DO
         END IF
C
         DO L=LF,LL
           WQRR(L) = WQVO(L,K,17) + DTWQ*WQRR(L) + DTWQO2*( WQKK(L)
     *       + WQKSUA(IWQT(L))*WQVO(L,K,16) + WQN17(L)*WQVO(L,K,17) )
         END DO
C
         IF (K.NE.KC) THEN
           DO L=LF,LL
             WQRR(L) = WQRR(L) + DTWQO2*WQT17(L)*WQVO(L,K+1,17)
           END DO
         END IF
C
         DO L=LF,LL
           WQKK(L) = 1.0 / (1.0 - DTWQO2*WQN17(L))
CTt update modified next line added line below commented out
       WQV(L,K,17)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,17)
ctt           WQV(L,K,17) = WQRR(L)*WQKK(L)
           WQVO(L,K,17) = WQVO(L,K,17)+WQV(L,K,17)
         END DO
C
       END IF
C
C18 COD: WQO18(L)=DTWQO2*WQO18
C
       DO L=LF,LL
         WQKK(L) = 1.0 / (1.0 - WQO18(L))
         WQRR(L) = (WQWDSL(L,K,18)+WQWPSL(L,K,18)) * VOLWQ(L)
c mrm  add wet atmospheric deposition:
         if (k.eq.kc) then
           WQRR(L) = WQRR(L) + WQATML(L,kc,18) * VOLWQ(L)
         end if
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
CTt update modified next line added line below commented out
       WQV(L,K,18)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,18)
ctt        WQV(L,K,18) = WQRR(L)*WQKK(L)
        WQVO(L,K,18) = WQVO(L,K,18)+WQV(L,K,18)
C
C19 O2: WQD19=-WQAOCR*WQKHR(L),WQL19=-WQAONT*WQNIT(L),WQO19=WQO18
C: WQKRDOS(L)=-WQP19(L)*WQDOS
C
        WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP19(L))
        WQRR(L) = (WQWDSL(L,K,19)+WQWPSL(L,K,19)) * VOLWQ(L)
        xDOpsl(L,K) = xDOpsl(L,K) + WQRR(L)*DTWQ*DZCWQ(K)*HPWQ(L)
        xDOall(L,K) = xDOall(L,K) + WQRR(L)*DTWQ*DZCWQ(K)*HPWQ(L)
c mrm  add wet atmospheric deposition:
        if (k.eq.kc) then
          WQRR(L) = WQRR(L) + WQATML(L,kc,19) * VOLWQ(L)
        end if
       END DO
C
       IF (K.EQ.KC) THEN
         DO L=LF,LL
           WQRR(L) = WQRR(L) + WQKRDOS(L)
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
         WQTTC = (1.3 - 0.3*WQPNC(L)) * WQPC(L)*DO_G
         WQTTD = (1.3 - 0.3*WQPND(L)) * WQPD(L)*DO_G
         WQTTG = (1.3 - 0.3*WQPNG(L)) * WQPG(L)*DO_G
         xDOppB(L,K) = xDOppB(L,K) + ( WQTTC*WQVO(L,K,1)
     *     + WQTTD*WQVO(L,K,2) + WQTTG*WQVO(L,K,3) ) * WQAOCR*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         xDOall(L,K) = xDOall(L,K) + ( WQTTC*WQVO(L,K,1)
     *     + WQTTD*WQVO(L,K,2) + WQTTG*WQVO(L,K,3) ) * WQAOCR*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)

         xmrm = CFCDCWQ*O2WQ*WQBMC(L)/(WQKHRC+O2WQ+ 1.E-18)
         WQA19C = WQTTC - xmrm
         xDOrrB(L,K) = xDOrrB(L,K) - xmrm*WQVO(L,K,1) * WQAOCR*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         xDOall(L,K) = xDOall(L,K) - xmrm*WQVO(L,K,1) * WQAOCR*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         xmrm = CFCDDWQ*O2WQ*WQBMD(L)/(WQKHRD+O2WQ+ 1.E-18)
         WQA19D = WQTTD - xmrm
         xDOrrB(L,K) = xDOrrB(L,K) - xmrm*WQVO(L,K,2) * WQAOCR*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         xDOall(L,K) = xDOall(L,K) - xmrm*WQVO(L,K,2) * WQAOCR*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         xmrm = CFCDGWQ*O2WQ*WQBMG(L)/(WQKHRG+O2WQ+ 1.E-18)
         WQA19G = WQTTG - xmrm
         xDOrrB(L,K) = xDOrrB(L,K) - xmrm*WQVO(L,K,3) * WQAOCR*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         xDOall(L,K) = xDOall(L,K) - xmrm*WQVO(L,K,3) * WQAOCR*DTWQO2
     *     *DZCWQ(K)*HPWQ(L)
         WQA19 = ( WQA19C*WQVO(L,K,1) + WQA19D*WQVO(L,K,2)
     *     + WQA19G*WQVO(L,K,3) ) * WQAOCR
C
C add macalgal source     J.S.
C
C Modified by MRM 05/23/99 to allow different AOCR constants to be applied
C   to photosynthesis and respiration terms for macroalgae:
C     WQAOCRpm = AOCR applied to macroalgae photosynthesis term
C     WQAOCRrm = AOCR applied to macroalgae respiration term
C
         IF(IDNOTRVA.GT.0.AND.K.EQ.1)THEN
           WQTTM = (1.3 - 0.3*WQPNM(L)) * WQPM(L)
           xmrm = (1.0-WQFCDM)*o2wq*WQBMM(L)/(WQKHRM+O2WQ+ 1.E-18)
c MRM ++++ start new code
c           WQA19A =  WQTTM - xmrm
c           WQA19 = WQA19 + WQA19A * WQVO(L,K,IDNOTRVA) * WQAOCR
           WQA19A = WQTTM * WQVO(L,K,IDNOTRVA) * WQAOCRpm -
     *       xmrm *  WQVO(L,K,IDNOTRVA) * WQAOCRrm
           WQA19 = WQA19 - WQA19A
c           WQA19 = WQA19 + WQA19A     ! Sign difference!, codemike.011701 has this! figure out why! 9/18/02
c MRM ++++ end new code
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
     *     - WQAOCR*WQKHR(L)*WQVO(L,K,6) - WQAONT*WQNIT(L)*WQVO(L,K,14)
     *     + WQO18(L)*WQVO(L,K,18) + WQP19(L)*WQVO(L,K,19) )
C
CTt update modified next line added line below commented out
       WQV(L,K,19)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,19)
ctt         WQV(L,K,19) = WQRR(L)*WQKK(L)
C
C MRM do not allow D.O. to go negative:
C
         wqv(L,K,19) = max (wqv(L,K,19), 0.0)
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
  !        xDOacoef(L,K) = xDOacoef(L,K) + WQP19(L)
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
C TAM: WQQ20=-WQWSSET(L,1),WQT20=WQWSSET(L,2)
C
       IF (IWQSRP.EQ.1) THEN       ! TAM is calculated only if IWQSRP=1!, Ji, 9/18/02
C
         DO L=LF,LL
           WQT20 = - DTWQO2*WQWSSET(L,1)
           WQKK(L) = 1.0 / (1.0 - WQT20)
           WQRR(L) = WQVO(L,K,20) + DTWQ*WQR20(L) + WQT20*WQVO(L,K,20)
         END DO
C
         IF (K.NE.KC) THEN
           DO L=LF,LL
             WQRR(L) = WQRR(L) + DTWQO2*WQWSSET(L,2)*WQVO(L,K+1,20)
           END DO
         END IF
C
         DO L=LF,LL
CTt update modified next line added line below commented out
       WQV(L,K,20)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,20)
ctt           WQV(L,K,20) = WQRR(L)*WQKK(L)
           WQVO(L,K,20) = WQVO(L,K,20)+WQV(L,K,20)
         END DO
C
       END IF
C
C FCB: WQTDFCB(IWQT(L))=WQKFCB*WQTFCB**(T-20),S21=-KFCB*TFCB^(T-20)
C WQTD1FCB=1+DTWQO2*WQS21,WQTD2FCB=1/(1-DTWQO2*S21)
C
       IF (IWQFCB.EQ.1) THEN
C
         DO L=LF,LL
           WQKK(L) = WQTD2FCB(IWQT(L))
           WQR21= (WQWDSL(L,K,NWQV)+WQWPSL(L,K,NWQV))*VOLWQ(L)          ! (3-19)
c mrm  add wet atmospheric deposition:
           if (k.eq.kc) then
             WQR21 = WQR21 + WQATML(L,kc,21) * VOLWQ(L)
           end if
           WQRR(L) = WQVO(L,K,NWQV)*WQTD1FCB(IWQT(L)) + DTWQ*WQR21
CTt update modified next line added line below commented out
       WQV(L,K,21)=SCBWQ(L)*(WQRR(L)*WQKK(L))+(1.-SCBWQ(L))*WQVO(L,K,21)
ctt           WQV(L,K,21) = WQRR(L)*WQKK(L)
           WQVO(L,K,21) = WQVO(L,K,21)+WQV(L,K,21)
         END DO
C
       END IF
C
C end of DO K=KC,1,-1  loop
C end of DO ND=1,NDMWQ loop
C
      END DO                          ! #200
      END DO                          ! #300
C
C increment counter for limitation and xDOxxx DO component arrays:
C
      TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/TIDALP
      timesum3 = timesum3 + TIMTMP
      NLIM = NLIM + 1
C
C compute WQCHL,WQTAMP,WQPO4D,WQSAD at a new time step: WQCHLx=1/WQCHLx
C  
! J.S. 2014    
!         IF(TWQ(L)>15) THEN
!         WQCHLAD=1000.0/(exp(-1.19*xKe(L))*90+30)
!         ELSE
!         WQCHLAD=1000.0/(exp(-1.18*xKe(L))*150+30 )        
!         ENDIF
!      DO K=1,KC
!        DO L=2,LA
!          WQCHL(L,K) = WQV(L,K,1)*WQCHLAD + 
!     *    WQV(L,K,2)*WQCHLAD
!     *      + WQV(L,K,3)*WQCHLAD
!        END DO
!      END DO        
      DO K=1,KC
        DO L=2,LA
          WQCHL(L,K) = WQV(L,K,1)*WQCHLC(IWQZMAP(L,KC)) + 
     *    WQV(L,K,2)*WQCHLD(IWQZMAP(L,KC))
     *      + WQV(L,K,3)*WQCHLG(IWQZMAP(L,KC))
        END DO
      END DO
      IF (IWQSRP.EQ.1) THEN
        DO K=1,KC
          DO L=2,LA
            O2WQ = MAX(WQV(L,K,19), 0.0)
            WQTAMD = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ), WQV(L,K,20) )
            WQTAMP(L,K) = WQV(L,K,20) - WQTAMD
            WQPO4D(L,K) = WQV(L,K,10) / (1.0 + WQKPO4P*WQTAMP(L,K))
            WQSAD(L,K)  = WQV(L,K,17) / (1.0 + WQKSAP*WQTAMP(L,K))
          END DO
        END DO
       ELSE IF (IWQSRP.EQ.2) THEN          ! SEDT
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
 
      Fflux=1.0                       ! J.S. 2007 This is a buge in the previous SLR model. Set to 1. for correct settling
c
      IF (IWQBEN.EQ.1) THEN           ! #400, Ji, 9/18/02

       DO ND=1,NDMWQ
       LF=2+(ND-1)*LDMWQ
       LL=LF+LDM-1
       DO L=LF,LL
         IMWQZ = IWQZMAP(L,1)
         SET_TMP=WQWSC(IMWQZ)         ! J.S. 2015 HAB
         IF(IMWQZ.EQ.4)SET_TMP=0.0005 ! forced to use low sellting for HAB
         WQDFBC(L) = SCBWQ(L)*SET_TMP*WQV(L,1,1)*Fflux
         WQDFBD(L) = SCBWQ(L)*WQWSD(IMWQZ)*WQV(L,1,2)*Fflux
         WQDFBG(L) = SCBWQ(L)*WQWSG(IMWQZ)*WQV(L,1,3)*Fflux
     +             +WQWSM*DZWQ(L)*WQV(L,1,IDNOTRVA)*Fflux
         WQDFRC(L) = SCBWQ(L)*WQWSRP(IMWQZ)*WQV(L,1,4)*Fflux
         WQDFLC(L) = SCBWQ(L)*WQWSLP(IMWQZ)*WQV(L,1,5)*Fflux
         WQDFRP(L) = SCBWQ(L)*WQWSRP(IMWQZ)*WQV(L,1,7)*Fflux
         WQDFLP(L) = SCBWQ(L)*WQWSLP(IMWQZ)*WQV(L,1,8)*Fflux
         WQDFRN(L) = SCBWQ(L)*WQWSRP(IMWQZ)*WQV(L,1,11)*Fflux
         WQDFLN(L) = SCBWQ(L)*WQWSLP(IMWQZ)*WQV(L,1,12)*Fflux
         IF (IWQSI.EQ.1) WQDFSI(L) = SCBWQ(L)*WQWSD(IMWQZ)*WQV(L,1,16)
       END DO
       END DO

        IF (IWQSRP.EQ.1) THEN
          DO ND=1,NDMWQ
          LF=2+(ND-1)*LDMWQ
          LL=LF+LDM-1
          DO L=LF,LL
            IMWQZ = IWQZMAP(L,1)
            WQDFLP(L) = SCBWQ(L)*( WQDFLP(L)
     *        + WQWSS(IMWQZ)*( WQV(L,1,10)-WQPO4D(L,1) ) )
            IF (IWQSI.EQ.1) WQDFSI(L) = SCBWQ(L)*( WQDFSI(L)
     *        + WQWSS(IMWQZ)*( WQV(L,1,17)-WQSAD(L,1) ) )
          END DO
          END DO
         ELSE IF (IWQSRP.EQ.2) THEN             ! SEDT
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
         ! RLIGHT1=WQKEB(IMWQZT(L))+WQKETSS(IMWQZT(L))*SEDTWQ(L,K)
          
           RLIGHT1= WQKEB(IMWQZT(L))+WQKETSS(IMWQZT(L))*SEDTWQ(L,K) 
     &     +WQKESAT(IMWQZT(L))*SALWQ(L,K)
     &     +WQKECHL(IMWQZT(L))*WQCHL(L,K)
C
C MRM 05/12/1999 use Riley (1956) equation to compute light extinction
C     as a function of CHL conc. if WQKECHL is less than zero:
C
C          RLIGHT2=WQKECHL*WQCHL(L,K)
          xmrm = WQKECHL(IMWQZT(L))*WQCHL(L,K)
          IF (WQKECHL(IMWQZT(L)) .LT. 0.0) THEN
            xmrm = 0.054*WQCHL(L,K)**0.6667 + 0.0088*WQCHL(L,K)
          END IF
      !     RLIGHT2 = xmrm
      !     RLIGHTT(L,K)=RLIGHTT(L,K)+RLIGHT1
      !     RLIGHTC(L,K)=RLIGHTC(L,K)+RLIGHT1+RLIGHT2
         
          RLIGHTC(L,K)=RLIGHT1
         END DO
        END DO
C
        IF(NDLTCNT.EQ.NSTPTMP) THEN  ! J.S. 2015
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
!           WRITE(1,1113)ILW(L),JLW(L),(RLIGHTT(L,K),K=1,KC),
!     $                               (RLIGHTC(L,K),K=1,KC)
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
 1113 FORMAT(2I5,12E12.4)
c
      RETURN
      END
