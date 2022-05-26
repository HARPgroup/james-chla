C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SSEDTOX_PCB(ISTL,CORDT)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
C **  SUBROUTINE SSEDTOX CALCULATES SETTLING AND WATER COLUMN-BED
C **  EXCHANGE OF SEDIMENT AND SORBED TOXIC CONTAMINANTS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'wq.par'      
      INCLUDE 'efdc.cmn'
      INCLUDE 'wqcom.cmn'
      LOGICAL, SAVE :: init=.true.
      DIMENSION TOXPFTWW(LCM,10)
      DIMENSION CTMPDRY(LCM),SEDFPS(LCM),RESBS1_S(LCM)
      DIMENSION WSETA(LCM,0:KSM,NSTM)
      DIMENSION  TOXS(LCM,KCM,NTXM)
      DIMENSION TOXBS(LCM,KBM,NTXM)
      DIMENSION  CDECAYW(LCM,KCM),CDECAYB(LCM,KBM)
      DIMENSION  ALOW(LCM,KBM+1),BMNN(LCM,KBM+1),CUPP(LCM,KBM+1),
     $         RRHS(LCM,KBM+1),TOXTMP(LCM,KBM+1),GAMTMP(LCM,KBM+1)
      DIMENSION GASC(LCM,3),TOX_TMP(KC)

C
C**********************************************************************C
C
      II_binary=2          ! =1 save3D binary file, =2 out binary 2D for homolog
      NTSORB =3            ! Only 3 sorbed material
!      RGASA(1)=0.273
!      RGASA(2)=0.198  
!      RGASA(3)=0.13	            
      
	IF(KBP.NE.2) THEN
	KBP = 2              ! Fix 2-layer
      write(*,*)'Check sediment layer in PCB model'
      ENDIF
      
      DELT=DT2
      S3TL=1.0
      S2TL=0.0
      ISUD=1
      IF (ISTL.NE.3) THEN
        DELT=DT
        S3TL=0.0
        S2TL=1.0
        ISUD=0
      END IF
      DELTI=1./DELT
C
      SEDMDGM=SQRT(SEDMDMX*SEDMDMN)

      BEDEX=1.

      NVAL=MOD(N,2)

C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON
      FOURDPI=4./PI
C
C**********************************************************************C
C
      IF(init)THEN
      init=.false.

      if(II_binary.eq.0) then
      
        OPEN(1,FILE='PCB.dia',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        
        open(2,file='pcbcon1.out',STATUS='UNKNOWN')
        CLOSE(2,STATUS='DELETE')
        open(2,file='pcbcon1.out',STATUS='UNKNOWN')  
        write(2,*)'I  J PCBb PCBs KWb  KWs GAs'
        write(2,*)LA-1
        close(2)

        open(2,file='pcbcon2.out',STATUS='UNKNOWN')
        CLOSE(2,STATUS='DELETE')
        open(2,file='pcbcon2.out',STATUS='UNKNOWN')  
        write(2,*)'I  J PCBb PCBs KWb  KWs GAs'
        write(2,*)LA-1
        close(2)
        
        open(2,file='pcbcon3.out',STATUS='UNKNOWN')
        CLOSE(2,STATUS='DELETE')
        open(2,file='pcbcon3.out',STATUS='UNKNOWN')  
        write(2,*)'I  J PCBb PCBs KWb  KWs GAs'
        write(2,*)LA-1
        close(2)
     
        open(2,file='pcbtcon.out',STATUS='UNKNOWN')
        CLOSE(2,STATUS='DELETE')
        open(2,file='pcbtcon.out',STATUS='UNKNOWN')  
        write(2,*)'L, TOX TOXB, Kbw'
        write(2,*)LA-1
        close(2)
 
! Save as binary
        elseif(II_binary.eq.2) then
        
        open(2,file='pcbcon1.bin',STATUS='UNKNOWN')
        CLOSE(2,STATUS='DELETE')
        open(801,file='pcbcon1.bin',STATUS='UNKNOWN',
     &    form='unformatted')    
        write(801)LA-1

        open(2,file='pcbcon2.bin',STATUS='UNKNOWN')
        CLOSE(2,STATUS='DELETE')
        open(802,file='pcbcon2.bin',STATUS='UNKNOWN',
     &    form='unformatted')     
        write(802)LA-1
        
        open(2,file='pcbcon3.bin',STATUS='UNKNOWN')
        CLOSE(2,STATUS='DELETE')
        open(803,file='pcbcon3.bin',STATUS='UNKNOWN',
     &    form='unformatted')      
        write(803)LA-1
     
!        open(2,file='pcbtcon.out',STATUS='UNKNOWN')
!        CLOSE(2,STATUS='DELETE')
!        open(2,file='pcbtcon.out',STATUS='UNKNOWN')  
!        write(2,*)'L, TOX TOXB, Kbw'
!        write(2,*)LA-1
!        close(2)       

        else
        
        OPEN(1,FILE='PCB.dia',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        
        open(2,file='PCB_water.bin',STATUS='UNKNOWN')
        CLOSE(2,STATUS='DELETE')
        open(801,file='PCB_water.bin',STATUS='UNKNOWN',
     &    form='unformatted')  
 !       write(2,*)'I  J PCBb PCBs KWb  KWs GAs'
        write(801)LA-1
 !       close(801)

        open(2,file='PCB_sed.bin',STATUS='UNKNOWN')
        CLOSE(2,STATUS='DELETE')
        open(802,file='PCB_sed.bin',STATUS='UNKNOWN',
     &    form='unformatted')  
 !       write(802,*)'I  J PCBb PCBs KWb  KWs GAs'
        write(802)LA-1
 !       close(2)
        
        open(2,file='PCB_parm1.bin',STATUS='UNKNOWN')
        CLOSE(2,STATUS='DELETE')
        open(803,file='PCB_parm.bin',STATUS='UNKNOWN',
     &    form='unformatted')  
 !       write(2,*)'I  J PCBb PCBs KWb  KWs GAs'
        write(803)LA-1
 !       close(2)
     
        open(2,file='PCB_parm2.bin',STATUS='UNKNOWN')
        CLOSE(2,STATUS='DELETE')
        open(804,file='PCB_parm2.bin',STATUS='UNKNOWN',
     &    form='unformatted')  
 !       write(2,*)'L, TOX TOXB, Kbw'
        write(804)LA-1
 !       close(2)
        
       endif
!                  
      DO L=2,LA
!	 DO K=1,KBP
!	  LL=IWQZMAP(L,1)
        SEDPOC(L,2)=SMPOC(L,1)+SMPOC(L,2)+SMPOC(L,3)  ! SMPOC are read from input or restart file
        SEDPOC(L,1)=2500000*(1.-BPRO2)*0.01                            
!       ENDDO
        RESBS1_S(L)=RESBS1
      ENDDO
      
	DO L=2,LA
        HBEDP(L,2)=HTOP         ! SMHSED(ISMZMAP(L))*2   !BBH2  Top layer use water wquality layer thickness
        HBEDP(L,1)=HBOT         
        PORBEDP(L,2)=BPRO1              !prosity
        PORBEDP(L,1)=BPRO2              !prosity
        KBTP(L) =KBP    
	ENDDO
C
C**********************************************************************C
C
C **  IF N=1 AND ISTRAN(5)=1 CHECK INITIAL TOXIC CONCENTRATIONS IN
C **  BED AND REINITILIZE IF NECESSARY
C
      IF(ISRESTI.EQ.0.OR.ISCI(5).EQ.0) THEN    ! not restart 
C
C **  CALCULATE TOTAL PARTICULATE FRACTION OF EACH TOXIC IN THE BED
C
C     TOXPFB(L,K,NS,NT)
C     NS =1, NTSORB, #of sorbed substance i.e. 3 carbons
C     NT = # of homologes of PCB 
c     TOXPARB(NS,NT) partition coefficient
c     NS= # of sorption substance

      DO NT=1,NPCB
       DO NS=1,NTSORB
        DO K=1,KBP
        DO L=2,LA
         TOXPFB(L,K,NS,NT)=0.
        END DO
        END DO
       END DO
      END DO
C
      DO NT=1,NPCB
        DO K=1,KBP
        DO L=2,LA
        TOXPFB(L,K,1,NT)= SEDPOC(L,K)*TOXPARB(1,NT)         ! Only one POC
        END DO
        END DO
      END DO
C
      DO NT=1,NPCB
       DO K=1,KBP
       DO L=2,LA
        TOXPFTB(L,K,NT)=0.
       END DO
       END DO
       DO K=1,KBP
       DO L=2,LA
         TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)+TOXPFB(L,K,1,NT)   ! Only one POC
       END DO
       END DO
      END DO
C
      DO NT=1,NPCB
       DO K=1,KBP
       DO L=2,LA
        IF(TOXPFTB(L,K,NT).GT.0.0) THEN                      !Total particulate fraction
          TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)
     &                   /(PORBEDP(L,K)*HBEDP(L,K)+TOXPFTB(L,K,NT))
         ELSE
          TOXPFTB(L,K,NT)=1.
        END IF
       END DO
       END DO
      END DO
C
C **  CONVERT MASS TOX/MASS SED INITIAL CONDITION TO TOTAL TOXIC
C **  CONCENTRATION IN BED 0.001 CONVERTS TOXINTB UNITS OF MG/KG
C **  TO TOXB UNITS OF OF MG/M**2 (or ug/kg-> ug/m^2)
C
C     FOR PCB input ug/kg or ug/M^3 convert to ug/m^2
C
      if(RESPCB.eq.0) then   ! not restart 
      
      DO NT=1,NPCB
      IF(ITXBDUT(NT).EQ.0) THEN
        DO K=1,KBP
        DO L=2,LA
         TOXB(L,K,NT)=HBEDP(L,K)*TOXB(L,K,NT)                  ! input ug/m^3 
         TOXB1(L,K,NT)=TOXB(L,K,NT)
        END DO
        END DO
      END IF
      IF(ITXBDUT(NT).EQ.1) THEN
!        DO K=1,KBP
        K=KBP
        DO L=2,LA
!         TOXB(L,K,NT)=0.001*TOXB(L,K,NT)*(SEDBT(L,K))    
!     $               /TOXPFTB(L,K,NT)
!         TOXB1(L,K,NT)=TOXB(L,K,NT)
        TOXB(L,K,NT)=HBEDP(L,K)*TOXB(L,K,NT)*2500.0*(1-PORBEDP(L,K))   ! input ug/kg  ->ug/m^2
     $               /TOXPFTB(L,K,NT)                                  ! use fixed partition coeff. for sed 
        TOXB(L,K-1,NT)=TOXB(L,K,NT)                               
        TOXB1(L,K,NT)=TOXB(L,K,NT)        ! force 2 layers equal
        TOXB1(L,K-1,NT)=TOXB(L,K,NT)
!        END DO
        END DO
      END IF
      END DO
      
      endif
C
C ** DIAGNOSTICS OF INITIALIZATION
C
      OPEN(2,file='toxbed.dia')
      CLOSE(2,STATUS='DELETE')
      OPEN(2,file='toxbed.dia')
	WRITE(2,*)'IL, JL, TOXPFTB,TOXB,TOXB/Hd,TOX'
      DO L=2,LA
       TMP1=TOXB(L,2,1)/HBEDP(L,2)/(2500.0*(1-PORBEDP(L,2))) 
       WRITE(2,2222)IL(L),JL(L),TOXPFTB(L,1,1),TOXB(L,2,1),
     $              TMP1,TOX(L,1,1)
      END DO
      CLOSE(2)
C
      END IF
      
      CALL BUDGET1   ! Ini. budget
      
      END IF  ! Complete initical setup if it is code start
c
 2222 FORMAT(2I5,7E13.4)
C
C**********************************************************************C
C
C **  SAVE OLD VALUES (Srart here for each call)
C

      DO NT=1,NPCB
       DO K=1,KC
        DO L=2,LA
         TOX1(L,K,NT)=TOX(L,K,NT)
        END DO
       END DO
      END DO
      DO NT=1,NPCB
       DO K=1,KBP
        DO L=2,LA
         TOXB1(L,K,NT)=TOXB(L,K,NT)
        END DO
       END DO
      END DO

C
C
C
C**********************************************************************
c
C      Update POC in the sediemtn. POC in th uplayer is from wq model.
c      Decay and burial considered in lower layer. Use WQ model parameter Card#9 and #19 
c      SMKPOC(3) = carbon decay m/day (group 3) 
c      SMW2(IZ)  = POC bural rate (m/day)
c
c      Top layer carbon input from water quality model. Group 3-carbon into one group
c
      DO L=2,LA 
      SEDPOC(L,KBP)=(SMPOC(L,1)+SMPOC(L,2)+SMPOC(L,3))            ! use input
 !     SEDPOC(L,KBP)=(SMPOC(L,1)+SMPOC(L,2)+SMPOC(L,3))/           ! POC at top sediment layer
 !    *  (1+RESPOC*DELT )                                   ! resuspension
 !     IZ = ISMZMAP(L)
 !	WWW1=SMW2(IZ)/86400.0/HBEDP(L,2)  ! m/d ->1/s
 !	DECY1=(SMKPOC(3)+SMKPOC(2))*0.1/86400.0/HBEDP(L,2)         !need temp correction  
 !     SEDPOC(L,KBP-1)=(SEDPOC(L,KBP-1)+WWW1*DELT*SEDPOC(L,KBP) )/   ! use group 3 decay
 !    &  (1+ DELT*(WWW1+DECY1))
 !      SEDPOC(L,KBP)=max(SEDPOC(L,KBP),5.0)
 !      SEDPOC(L,KBP-1)=max( SEDPOC(L,KBP-1),5.0)
      ENDDO

C
C
C**********************************************************************C
C
C **  COHESIVE SEDIMENT, KC=1 (SINGLE LAYER IN VERTICAL)
C
         DSEDGMM=1./(1.E6*SSG(1))   ! =sden, efdc.inp, C39

C
C----------------------------------------------------------------------C
C
C **  set settling velocities
c     WQWSLP and WQWSG are input from WQ model
C
      DO K=1,KC
       DO L=2,LA
	  LL=IWQZMAP(L,K)
        WSETA(L,K,1)= WQWSLP(LL)/86400.0 *F_SET    ! POC  (m/second)
        WSETA(L,K,2)= WQWSG(LL)/86400.0 *F_SET     ! PAC  
       END DO
      ENDDO
C
C----------------------------------------------------------------------C
C
c **  horizontal loop 
c     Set vertical carbon flux due to settling
c
      DO K=1,KC
      DO NS=1,NTSORB
      DO L=2,LA
        SEDFP(L,K,NS)=0.
      END DO
      ENDDO
	ENDDO

	DO K=0,KC-1
	DO L=2,LA
        SEDFP(L,K,1) = -WSETA(L,K+1,1)*(WQV(L,K+1,4)+WQV(L,K+1,5))      !POC
	  SEDFP(L,K,2) = -WSETA(L,K+1,2)*                                 !PAC 
     *                (WQV(L,K+1,1)+WQV(L,K+1,2)+WQV(L,K+1,3)) 
      END DO
      ENDDO
C 
C     Considering resuspension from the bottom layer
C 
!      tau_c=1.0e-5
      IF(RESPOC.GT.0.0) THEN
       DO L=2,LA
        if(T_cr.GT.1.0e-8)then
         tau_1=sqrt(UWQ(L,1)*UWQ(L,1)+VWQ(L,1)*VWQ(L,1))
         HVDZBR=HP(L)*DZC(1)*0.5/0.002
         cd_1=.16/( LOG(HVDZBR)**2)
         tau_1=cd_1*tau_1*tau_1
           if(N. GE.1440*4.and.N.LE.1440*7)then 
             if(mod(N,60).EQ.0) then
              write(193,'(F12.2,I5,E12.4)')Time,L,tau_1
             endif  
           endif
         if( (tau_1-T_cr*TAU_MD(L)).GT.0) then
               garm_1=min((tau_1-T_cr*TAU_MD(L))**Tauexp/
     &         (T_cr*TAU_MD(L)),100.0)  
         else
               garm_1=0.0
         endif
        else
         garm_1=1.0
        endif 
       DO K=1,NTSORB
!	  SEDFP(L,0,K)=SEDFP(L,0,K)+garm_1*SEDPOC(L,KBP)*RESPOC        ! constant resuspension of POC
	 ENDDO
	  if(T_cr.GT.1.0e-8) then
	   TT00=WQV(L,1,1)+WQV(L,1,2)+WQV(L,1,3)
	   TT1=(WSETA(L,1,1)*WQV(L,1,5)+WSETA(L,1,2)*TT00)/
     &   (WQV(L,1,5)+TT00+1.0e-8)
	   TT0=garm_1*RESPOC/HBEDP(L,KBP)*TAU_MD1(L)
	   TT0=min(TT0,TT1) 
	   SEDFPS(L)=TT0*SEDPOC(L,KBP)
!	   TT0=(SEDFP(L,0,1)+SEDFP(L,0,2)-SEDFPS(L))/HBEDP(L,KBP)*DELT
!     *      *SEDPOC(L,KBP)              ! This is max can do now
	 else
	   SEDFPS(L)=garm_1*SEDPOC(L,KBP)*RESPOC/HBEDP(L,KBP)*TAU_MD(L)  ! Constant erosion 
!	   SEDFPS(L)=TT0     ! if no erosion, no depostion to balance the mass
!	   WSETA(L,1,2)=0
	 endif	  
        if(mod(N,60).EQ.0.and.N<=1440*4) then
          if(LIJ(ICSHOW,JCSHOW).eq.L ) then
          write(991,9991)Time,tau_1,garm_1,TT0,TT1,
     &    SEDFPS(L),SEDPOC(L,KBP),T_cr*TAU_MD(L),TAU_MD1(L)
          endif  
        endif 
 9991   format(F12.2,8E12.4)  
	 RESBS1_S(L)= garm_1*RESBS1/HBEDP(L,KBP)*TAU_MD1(L) !SEDFPS(L) ! use same rate as POC. POC/SED =0.0x fix ratio                                ! seperate settling and resuspension
 !       if(LIJ(49,95 ).eq.L) then
 !        write(991,*)Time,SEDFPS(L),1/TAU_MD(L),SEDPOC(L,KBP),RESPOC
 !       endif
      ENDDO
      ENDIF
        
      If(N>2880) close(51)
        
      IF(idiaout.eq.3) THEN 
       if(mod(N,ltimstep).eq.0) then
        write(188,*)'WSETA 1 WSETA 2 POC  Algae'
        write(188,*) WSETA(lcadia,1,1), WSETA(lcadia,2,2),
     &   WQV(lcadia,1,4)+WQV(lcadia,1,5),
     &   WQV(lcadia,1,2)+WQV(lcadia,1,3)
        write(188,*) 'SEDEP 0-1 RESPOC SEDPOC '
        write(188,*) SEDFP(lcadia,0,1),
     &   SEDFP(lcadia,1,1),RESPOC,SEDPOC(lcadia,2)
       endif
      ENDIF      
C
C**********************************************************************C
C
C **  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS
C **  IN WATER COLUMN
C
C **  TOXPFW(L,K,NS,NT) = PARTICULATE FRACTION IN WATER COLUMN
C **  TOXPARW(NS,NT) =    PARTITION COEFFICIENT IN WATER COLUMN
C **  TOXPFTW(L,K,NT) =   TOTAL PARTICULATE FRACTION IN WATER COLUMN
C                         USED AS TEMPORARY VARIBLE IN THIS AND
C                         FOLLOWING CODE BLOCK
C
C
      DO NT=1,NPCB
       DO NS=1,NTSORB
        DO K=1,KC
         DO L=2,LA
          TOXPFW(L,K,NS,NT)=0.
         END DO
        END DO
       END DO
      END DO
C
      L00=2132
 1880 format(8F10.4)
      DO NT=1,NPCB
!       IF(ITXPARW(NS,NT).EQ.0) THEN
        DO K=1,KC
         DO L=2,LA

          a_carb=(WQV(L,K,4)+WQV(L,K,5))
          TOXPFW(L,K,1,NT)=a_carb*TOXPARW(1,NT)        ! POC*Fp

          a_agle=(WQV(L,K,1)+WQV(L,K,2)+WQV(L,K,3))
          TOXPFW(L,K,2,NT)=a_agle * TOXPARW(2,NT)      ! PAC*Fp                

          TOXPFW(L,K,3,NT)=WQV(L,K,6)*TOXPARW(3,NT)                     ! DOC*Fp
         IF(idiaout.eq.3) THEN 
          if(mod(N,ltimstep).eq.0.and.L.eq.lcadia.and.K.eq.1) then
          write(188,*)'TEOXPFW 1-2  POC Alage  DOC'
          write(188,1880)TOXPFW(L,1,1,NT),TOXPFW(L,1,2,NT)
     &    ,TOXPFW(L,1,3,NT),a_carb,a_agle,WQV(L,K,6)
          endif
         ENDIF
         END DO
        END DO
!       END IF
      END DO
C
      DO NT=1,NPCB
       DO K=1,KC
        DO L=2,LA
         TOXPFTW(L,K,NT)=0.
        END DO
       END DO
       DO K=1,KC
        DO L=2,LA
         DO NS=1,NTSORB  !=3
          TOXPFTW(L,K,NT)=TOXPFTW(L,K,NT)+TOXPFW(L,K,NS,NT)
         END DO
        END DO
       END DO
      END DO
 
!       if(mod(N,480).eq.0)then
!       open(2,file='pcbcon1.out',ACCESS='APPEND')
!       write(2,*)TIME
!       NT=1
!       DO L=2,LA
!       write(2,201)IL(L),JL(L)
!     &  ,TOXPFW(L,KC,1,NT),TOXPFW(L,KC,2,NT),TOXPFTW(L,KC,NT)
!       ENDDO
!       close(2) 
!       endif     
      DO NT=1,NPCB
      DO L=2,LA
       TOXPFTWW(L,NT)=TOXPFW(L,1,1,NT)+TOXPFW(L,1,2,NT)
      ENDDO
      ENDDO
C
C**********************************************************************C
C
C **  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS
C **  IN SEDIMENT BED
C
C **  TOXPFB(L,NS,NT) = PARTICULATE FRACTION IN SEDIMENT BED
C **  TOXPARB(NS,NT) = PARTITION COEFFICIENT IN SEDIMENT BED
C **  TOXPFTB(L,NT) = TOTAL PARTICULATE FRACTION IN SEDIMENT BED
C                       USED AS TEMPORARY VARIBLE IN THIS AND
C                       FOLLOWING CODE BLOCK
C     SEDPOC(L,K)  = POC in sedimement layer K
C
      DO NT=1,NPCB
	 DO NS=1,NTSORB
        DO K=1,KBP
        DO L=2,LA
         TOXPFB(L,K,NS,NT)=0.
        END DO
        END DO
       END DO
      END DO
C
      DO NT=1,NPCB
        DO K=1,KBP
        DO L=2,LA
         TOXPFB(L,K,1,NT)=SEDPOC(L,K)*TOXPARB(1,NT)          ! POC only
         TOXPFTB(L,K,NT)=TOXPFB(L,K,1,NT)
        END DO
        END DO
      END DO
C
      IF(idiaout.eq.3) THEN 
       if(mod(N,ltimstep).eq.0) then
        write(188,*)'TOXPFB_2 SEDPOC', TOXPFB(lcadia ,2,1,1),
     &     SEDPOC(lcadia,2)
        write(188,*)'TOXPFB_1 SEDPOC', TOXPFB(lcadia ,1,1,1),
     &     SEDPOC(lcadia,1)
       endif
      ENDIF
         
      DO NT=1,NPCB
       DO NS=1,1
        DO K=1,KBP
        DO L=2,LA
           HBEDTMP=HBEDP(L,K)+1.E-9
           TOXPFB(L,K,NS,NT)=TOXPFB(L,K,NS,NT)/
     $                   (PORBEDP(L,K)*HBEDTMP+TOXPFTB(L,K,NT))
        END DO
        END DO
       END DO
      END DO
C
C**********************************************************************C
C
C **  CALCULATE PARTICULATE TOXIC CONTAMINANT SETTLING
C **  AND BED EXCHANGE FLUXES
C
C **  TOXF(L,KC,NT) = TOXIC CONTAMINANT SETTLING AND BED EXCHANGE
C                        FLUX.  USED AS TEMPORARY VARIABLE IN THIS
C                        AND FOLLOWING CODE BLOCK
C

      DO NT=1,NPCB
       DO K=0,KC
        DO L=2,LA
         TOXF(L,K,NT)=0.
        END DO
       END DO
      END DO
C
      IF(KC.GE.2) THEN
      DO NT=1,NPCB
       DO NS=1,NTSORB
          DO K=1,KS
           DO L=2,LA
            TOXF(L,K,NT)=TOXF(L,K,NT)+SEDFP(L,K,NS)*TOXPARW(NS,NT)    ! settling from lop layer
           END DO
          END DO
       END DO
      END DO
      END IF
C
      IF(KC.GE.2) THEN
      DO NT=1,NPCB
       DO K=1,KS
        DO L=2,LA
         TOXF(L,K,NT)=TOXF(L,K,NT)/(1.+TOXPFTW(L,K+1,NT))             ! Get fraction Pfd_i
        END DO
       END DO
      END DO
      END IF
C
      DO NT=1,NPCB
       DO NS=1,3
          DO L=2,LA
           TOXF(L,0,NT)=TOXF(L,0,NT) +
     $                   MIN(SEDFP(L,0,NS),0.)*TOXPARW(NS,NT)
!     $    /(1.+TOXPFTW(L,1,NT))                              ! Add settling (<0)
          END DO
       END DO
	ENDDO
C
      IF(idiaout.eq.3) THEN 
       if(mod(N,ltimstep).eq.0) then
        write(188,*)'TOXF(0,1) SEDFP1-2  TOXPARW 1-2',
     &   TOXF(lcadia,0,1),SEDFP(lcadia,0,1),
     &   SEDFP(lcadia,0,2),TOXPARW(1,1)  
        write(188,*)'TOXF(1,1)', TOXF(lcadia,1,1),SEDFP(lcadia,1,2)
       endif
      ENDIF
      
      DO NT=1,NPCB
       DO L=2,LA
        TOXF(L,0,NT)=TOXF(L,0,NT)/(1.+TOXPFTW(L,1,NT))
       END DO
      END DO
C
      DO NT=1,NPCB
       DO L=2,LA
        TOXFB(L,KBP,NT)=0.
       END DO
      END DO
C
      DO NT=1,NPCB
       DO NS=1,1
        DO L=2,LA
!         TOXFB(L,KBP,NT)=TOXFB(L,KBP,NT)                                  !fixed layers, otherwise use KBT(L)
!     &                 +MAX(SEDFP(L,0,NS),0.)*TOXPARB(NS,NT)             ! net lost due to Resuspension > 0 

         TOXFB(L,KBP,NT)=TOXFB(L,KBP,NT)                                  !fixed layers, otherwise use KBT(L)
     &                 +MAX(SEDFPS(L),0.)*TOXPARB(NS,NT)/
     &    (PORBEDP(L,KBP)*HBEDP(L,KBP)+TOXPFTB(L,KBP,NT))
     
	    WWW1=SETT_S1  !WWW1=SMW2(ISMZMAP(L))/86400.0 !
          TOXFB(L,KBP-1,NT)=0.  !-SEDPOC(L,KBP-1)*WWW1*TOXPARB(NS,NT)           ! Burial
     
        if(mod(N,60).EQ.0.and.N<=1440*4) then
                  if(LIJ(ICSHOW,JCSHOW).eq.L ) then
          !   if(LIJ(310,71).eq.L ) then
            write(992,*)Time,NT,TOXFB(L,KBP,NT),SEDFPS(L),TOXPARB(NS,NT)
             endif  
        endif           
          
        END DO                                                          
       END DO
      END DO
C
C ** DIAGNOSTICS OF FLUX
C
       IF(N.EQ.1)THEN
       OPEN(2,file='toxflx.dia')
       CLOSE(2,STATUS='DELETE')
       OPEN(2,file='toxflx.dia')
       DO L=2,LA
        WRITE(2,2222)IL(L),JL(L),HBEDP(L,KBP),TOXPFTB(L,KBP,1),
     $         TOXFB(L,KBTP(L),1),TOXF(L,0,1)
       END DO
       CLOSE(2)
       END IF
c
C
C**********************************************************************C
C
C **  CALCULATE TOTAL PARTICULATE FRACTIONS IN WATER COLUMN AND BED
C
C
!      DO NT=1,NPCB
!       DO K=1,KC
!        DO L=2,LA
!          TOXPFTW(L,K,NT)=TOXPFTW(L,K,NT)/(1.+TOXPFTW(L,K,NT))
!        END DO
!       END DO
!      END DO
C
      DO NT=1,NPCB
      DO L=2,LA
       TOXPFTWW(L,NT)=TOXPFTWW(L,NT)/(1.+TOXPFTWW(L,NT))
      ENDDO
      ENDDO
!      if(mod(N,20).eq.0) write(*,*)TOXPFTWW(100,2),TOXPFTW(100,1,2)
c      
      DO NT=1,NPCB
       DO K=1,KBP
       DO L=2,LA
          HBEDTMP=HBEDP(L,K)+1.E-12
          TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)
     $                   /(PORBEDP(L,K)*HBEDTMP+TOXPFTB(L,K,NT))
       END DO
       END DO
      END DO
C
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT, KC=1 (SINGLE LAYER IN VERTICAL)
C
      
      IF(KC.EQ.1) THEN
        DO NT=1,NPCB
C
C----------------------------------------------------------------------C
C
         DO L=2,LA
          WVEL=DELTI*HP(L)*DZC(1)
          AA11=WVEL-TOXF(L,0,NT)
          AA12=-TOXFB(L,KBTP(L),NT)
          AA21=TOXF(L,0,NT)
          AA22=DELTI+TOXFB(L,KBTP(L),NT)
          BB11=WVEL*TOX(L,1,NT)
          BB22=DELTI*TOXB1(L,KBTP(L),NT)
          DETI=1./(AA11*AA22-AA12*AA21)
          TOX(L,1,NT)=DETI*(BB11*AA22-BB22*AA12)
          TOXB1(L,KBTP(L),NT)=S3TL*TOXB(L,KBTP(L),NT)
     $                      +S2TL*TOXB1(L,KBTP(L),NT)
          TOXB(L,KBTP(L),NT)=DETI*(AA11*BB22-AA21*BB11)
         END DO
C
  676 FORMAT('N,I,J,TW,TB,TB1,TF,TFB=',3I5,5E13.4)
C
C----------------------------------------------------------------------C
C
        END DO
      END IF
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT, KC=2 (TWO LAYERS IN VERTICAL)
C
      IF(KC.EQ.2) THEN
        DO NT=1,NPCB
C
C----------------------------------------------------------------------C
C
         K=2
         DO L=2,LA
          WVEL=DELTI*HP(L)*DZC(K)
          CLEFT=WVEL-TOXF(L,K-1,NT)
          CRIGHT=WVEL*TOX(L,K,NT)
          TOX(L,K,NT)=CRIGHT/CLEFT
         END DO
C
         DO L=2,LA
	    KBTP(L) =KBM                     ! Fix at top layer
          WVEL=DELTI*HP(L)*DZC(1)
          AA11=WVEL-TOXF(L,0,NT)
          AA12=-TOXFB(L,KBTP(L),NT)
          AA21=TOXF(L,0,NT)
          AA22=DELTI+TOXFB(L,KBTP(L),NT)
          BB11=WVEL*TOX(L,1,NT)-TOXF(L,1,NT)*TOX(L,2,NT)
          BB22=DELTI*TOXB1(L,KBTP(L),NT)
          DETI=1./(AA11*AA22-AA12*AA21)
          TOX(L,1,NT)=DETI*(BB11*AA22-BB22*AA12)
          TOXB1(L,KBTP(L),NT)=S3TL*TOXB(L,KBTP(L),NT)
     $                      +S2TL*TOXB1(L,KBTP(L),NT)
          TOXB(L,KBTP(L),NT)=DETI*(AA11*BB22-AA21*BB11)
         END DO
C
C----------------------------------------------------------------------C
C
        END DO
      END IF
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT, KC=3 (THREE LAYERS IN VERTICAL)
C
      IF(KC.EQ.3) THEN
        DO NT=1,NPCB
C
C----------------------------------------------------------------------C
C
         K=3
         DO L=2,LA
          WVEL=DELTI*HP(L)*DZC(K)
          CLEFT=WVEL-TOXF(L,K-1,NT)
          CRIGHT=WVEL*TOX(L,K,NT)
          TOX(L,K,NT)=CRIGHT/CLEFT
         END DO
C
         K=2
         DO L=2,LA
          WVEL=DELTI*HP(L)*DZC(K)
          CLEFT=WVEL-TOXF(L,K-1,NT)
          CRIGHT=WVEL*TOX(L,K,NT)-TOXF(L,K,NT)*TOX(L,K+1,NT)
          TOX(L,K,NT)=CRIGHT/CLEFT
         END DO
C
         DO L=2,LA
	    KBTP(L) =KBM                     ! Fix at top layer
          WVEL=DELTI*HP(L)*DZC(1)
          AA11=WVEL-TOXF(L,0,NT)
          AA12=-TOXFB(L,KBTP(L),NT)
          AA21=TOXF(L,0,NT)
          AA22=DELTI+TOXFB(L,KBTP(L),NT)
          BB11=WVEL*TOX(L,1,NT)-TOXF(L,1,NT)*TOX(L,2,NT)
          BB22=DELTI*TOXB1(L,KBTP(L),NT)
          DETI=1./(AA11*AA22-AA12*AA21)
          TOX(L,1,NT)=DETI*(BB11*AA22-BB22*AA12)
          TOXB1(L,KBTP(L),NT)=S3TL*TOXB(L,KBTP(L),NT)
     $                      +S2TL*TOXB1(L,KBTP(L),NT)
          TOXB(L,KBTP(L),NT)=DETI*(AA11*BB22-AA21*BB11)
         END DO
C
C----------------------------------------------------------------------C
C
        END DO
      END IF
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT, KC.GE.3 (THREE OR MORE LAYERS IN VERTICAL)
C  

      IF(KC.GT.3) THEN
        DO NT=1,NPCB
C
C----------------------------------------------------------------------C
C
       
         K=KC
         DO L=2,LA
          WVEL=DELTI*HP(L)*DZC(K)
          CLEFT=WVEL-TOXF(L,K-1,NT)
          CRIGHT=WVEL*TOX(L,K,NT)
          TOX(L,K,NT)=CRIGHT/CLEFT
         END DO
C
          
         DO K=KS,2,-1
          DO L=2,LA
           WVEL=DELTI*HP(L)*DZC(K)
           CLEFT=WVEL-TOXF(L,K-1,NT)
           CRIGHT=WVEL*TOX(L,K,NT)-TOXF(L,K,NT)*TOX(L,K+1,NT)
           TOX(L,K,NT)=CRIGHT/CLEFT
          END DO
         END DO
C

         DO L=2,LA
          if(HBEDP(L,KBP).GT.0.01) THEN
!	    KBTP(L)=2                        ! Fixed top layer for 2 bottom layer model
          WVEL=DELTI*HP(L)*DZC(1)
          AA11=WVEL-TOXF(L,0,NT)
          AA12=-TOXFB(L,KBP,NT)
          AA21=TOXF(L,0,NT)
          AA22=DELTI+TOXFB(L,KBP,NT) ! -TOXFB(L,KBP-1,NT)           ! add settling to the next lower layer
          BB11=WVEL*TOX(L,1,NT)-TOXF(L,1,NT)*TOX(L,2,NT)
          BB22=DELTI*TOXB1(L,KBTP(L),NT)
          DETI=1./(AA11*AA22-AA12*AA21+1.0e-14)
          TOX(L,1,NT)=DETI*(BB11*AA22-BB22*AA12)
          TOXB1(L,KBTP(L),NT)=S3TL*TOXB(L,KBTP(L),NT)
     $                      +S2TL*TOXB1(L,KBTP(L),NT)
          TOXBTMP=DETI*(AA11*BB22-AA21*BB11)
          IF(TOXBTMP.LT.0.0) THEN
            TOXBTMP=TOXB1(L,KBTP(L),NT)                               ! no resuspension/settling
            TOX(L,KBTP(L),NT)=TOXS(L,1,NT)-TOXF(L,1,NT)*TOX(L,2,NT)/WVEL   ! Only consider settling
          END IF

          TOXB(L,KBTP(L),NT)=TOXBTMP
!
! Second step to update bottom toxics
!
          AA11=1+DELT*TOXPFTB(L,KBTP(L),NT)*SETT_S1/HBEDP(L,KBP)
          AA12=-DELT*RESBS2* TOXPFTB(L,KBTP(L)-1,NT)/HBEDP(L,KBP-1) 
          AA21=-DELT*TOXPFTB(L,KBTP(L),NT)*SETT_S1/HBEDP(L,KBP) 
          AA22= 1+DELT*(RESBS2+SETT_S2)/HBEDP(L,KBP-1)*
     &    TOXPFTB(L,KBTP(L)-1,NT)
          BB11=TOXB(L,KBTP(L),NT)
          BB22=TOXB(L,KBTP(L)-1,NT)          
          DETI=1./(AA11*AA22-AA12*AA21+1.0e-14)
          TOXB(L,KBTP(L),NT)=DETI*(BB11*AA22-BB22*AA12)
          TOXB(L,KBTP(L)-1,NT)=DETI*(AA11*BB22-AA21*BB11) 
          SEDERR=SEDERR+SETT_S2*TOXB(L,KBTP(L)-1,NT)*DXYP(L)
     &        *TOXPFTB(L,KBTP(L)-1,NT)/HBEDP(L,KBP-1) 
         ELSE
          TOXB(L,KBTP(L),NT)=0.0
         ENDIF       
         END DO
C
C----------------------------------------------------------------------C
C
        END DO
      END IF

        
C
C**********************************************************************C
C
C **  UPDATE SEDIMENT BED PHYSICAL PROPERTIES
C
C
C
c------------------------------------------------------------------
C  Constant layer thickness model: 1 One-layer model
C
!      DO L=2,LA
!	 A_TT=HBEDP(L,KBTP(L))+DELT*SETT_S1+RESBS1
!	 PORBEDP(L,K)=1.0-HBEDP(L,KBTP(L))*(1-PORBED1(L,K))/A_TT+
!    & DELT*(SEDFP(L,K,1)+SEDFP(L,K,2))/A_TT * DSEDGMM *
!     & (CSRITO+1)/CSRITO
!      ENDDO

c------------------------------------------------------------------
C   Constant sediment layer thinckness 2-layer model   (Top layer k=2, bottom layer k=1)
C   model (1-pro). i,e, Fis depth and vary prosity
C  
      dcorr=0.002
      if(iBSEDOP.NE.0)then
      
      if(iBSEDOP.EQ.1) then
      
       DO L=2,LA
	  S_TT=(RPOC2IS+1)/RPOC2IS
        AA11=(1+DELT*(RESBS1_S(L)+SETT_S1)/HBEDP(L,KBP) )
	  AA12=-DELT*RESBS2/HBEDP(L,KBP-1)
        AA21=-DELT*SETT_S1/HBEDP(L,KBP)
        AA22=1+DELT*(RESBS2+SETT_S2)/HBEDP(L,KBP-1)
        BB11=-DELT*(SEDFP(L,0,1)+SEDFP(L,0,2))/HBEDP(L,2)*DSEDGMM*S_TT         ! POC->SED  
     &   +(1-PORBEDP(L,2))
	  BB22=(1-PORBEDP(L,1))
        DETI=AA11*AA22-AA21*AA12
        ax1=(BB11*AA22-BB22*AA21)/DETI
	  ax2=(AA11*BB22-AA21*BB11)/DETI
        PORBEDP(L,2)=(1-ax1)
        PORBEDP(L,1)=(1-ax2)
       ENDDO
      endif
      
       if(iBSEDOP.EQ.2.or.iBSEDOP.EQ.3) then   ! Constant prosity, but change depth
     
       DO L=2,LA
c        SETT_SS1=SETT_S1/HBEDP(L,KBP)
c        RESBSS2=RESBS2/HBEDP(L,KBP-1)
c        SETT_SS2=SETT_S2/HBEDP(L,KBP-1)
        S_TT=(RPOC2IS+1)/RPOC2IS
        S_TT=-(SEDFP(L,0,1)+SEDFP(L,0,2))*DSEDGMM*S_TT
        HBEDP(L,KBP)=HBEDP(L,KBP)+DELT*(S_TT/(1-PORBEDP(L,2)) 
     &   -(RESBS1_S(L)+SETT_S1-RESBS2*(1-PORBEDP(L,1))
     &  /(1-PORBEDP(L,2))))
        HBEDP(L,KBP)=max(0.01,HBEDP(L,KBP)) 
        HBEDP(L,KBP-1)=HBEDP(L,KBP-1)+DELT*( 
     &  SETT_S1 *(1-PORBEDP(L,2))/(1-PORBEDP(L,1))-(RESBS2+SETT_S2) )
        HBEDP(L,KBP-1)=max(0.01,HBEDP(L,KBP-1)) 
      
       IF(iBSEDOP.EQ.3) then
        if(HBEDP(L,KBP).LE.(HTOP-dcorr)) then
         if(HBEDP(L,KBP-1).GT.dcorr) then
          DH0=HBEDP(L,KBP)
          DHH=HTOP-HBEDP(L,KBP)
          DHH1=min(DHH,HBEDP(L,KBP-1)-dcorr)
          HBEDP(L,KBP)=HBEDP(L,KBP)+DHH1
          HBEDP(L,KBP-1)=max(HBEDP(L,KBP-1)-DHH1,dcorr)
          DO NT=1,NPCB
          TOXB(L,KBP,NT)=(TOXB(L,KBP,NT)*DH0+TOXB(L,KBP-1,NT)*DHH1)
     &          /(DH0+DHH1)  
    
          ENDDO
         else
          if(HBEDP(L,KBP).LT.0.01)HBEDP(L,KBP)=min(HBEDP(L,KBP),dcorr)
         endif
        endif
       ENDIF
       ENDDO
           
      ENDIF
     
      ENDIF
C**********************************************************************C
C
C **  UPDATE TOTAL PARTICULATE FRACTION OF EACH TOXIC IN THE BED
C
      DO NT=1,NPCB
       DO NS=1,NSED+NSND
        DO K=1,KBP
        DO L=2,LA
         TOXPFB(L,K,NS,NT)=0.
        END DO
        END DO
       END DO
      END DO
C
      DO NT=1,NPCB
        DO K=1,KBP
        DO L=2,LA
         TOXPFB(L,K,1,NT)=SEDPOC(L,K)*TOXPARB(1,NT)
	   TOXPFTB(L,K,NT)=TOXPFB(L,K,1,NT)
        END DO
        END DO
      END DO
C
      DO NT=1,NPCB
       DO K=1,KBP
       DO L=2,LA
        IF(HBEDP(L,K).GT.0.0) THEN
          TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)
     &                   /(PORBEDP(L,K)*HBEDP(L,K)+TOXPFTB(L,K,NT))
         ELSE
          TOXPFTB(L,K,NT)=0.99
        END IF
       END DO
       END DO
      END DO
C
C**********************************************************************C
C
C **  DIFFUSE TOXICS IN BED AND INTO BOTTOM WATER COLUMN LAYER
C
! Try new 2-layer model
!
       
       LL=LIJ(26,89)
       
       DO NT=1,NPCB
       DO L=2,LA
       BT_1=DELTI*DZC(1)*HP(L)+DIFTOX(NT)*(1.-TOXPFTWW(L,NT)) 
       GM_1=-DIFTOX(NT)/(PORBEDP(L,2)*HBEDP(L,2))*(1.-TOXPFTB(L,2,NT))
       AF_2=-DIFTOX(NT)*(1.-TOXPFTWW(L,NT)) 
       BT_2=(DELTI-GM_1)
       GM_2=-DIFTOXB(NT)/(PORBEDP(L,1)*HBEDP(L,1))*(1.-TOXPFTB(L,1,NT))
       AF_3=-DIFTOXB(NT)/(PORBEDP(L,2)*HBEDP(L,2))*(1.-TOXPFTB(L,2,NT)) 
       BT_3=(DELTI-GM_2)
       BB_1=DELTI*DZC(1)*HP(L)*TOX(L,1,NT)
       BB_2=DELTI*TOXB(L,2,NT)
       Bb_3=DELTI*TOXB(L,1,NT) 
       alm_1=GM_1/BT_1
       d_1=BB_1/BT_1
       alm_2=GM_2/(BT_2-AF_2*alm_1)
       d_2=(BB_2-AF_2*d_1)/(BT_2-AF_2*alm_1)
       TOXB(L,1,NT)=(BB_3-AF_3*d_2)/(BT_3-AF_3*alm_2)
       TOXB(L,2,NT)=d_2-alm_2*TOXB(L,1,NT)
        TOX(L,1,NT)=d_1-alm_1*TOXB(L,2,NT)
       
       SMASSIN=SMASSIN+(AF_2*TOX(L,1,NT)-TOXB(L,2,NT)*GM_1)*DXYP(L) ! >0 difuse to water column
              
       
       IF(idiaout.eq.4) THEN 
       IF(L.EQ.LL.and.mod(N,240).eq.0.and.NT.eq.2) then
       write(501,501)TIME,TOX(LL,1,2), 
     &  TOX(LL,1,2)*(1.-TOXPFTWW(LL,2))
       write(501,501)TIME,TOXB(LL,2,2),
     &  TOXB(LL,2,2)*(1.-TOXPFTB(LL,2,2))/(PORBEDP(L,2)*HBEDP(L,2)),
     &  (1.-TOXPFTB(LL,2,2)),PORBEDP(L,2),SEDPOC(LL,2)
       write(501,501)TIME,TOXB(LL,1,2),
     &  TOXB(LL,2,2)*(1.-TOXPFTB(LL,1,2))/(PORBEDP(L,1)*HBEDP(L,1)),
     & (1.-TOXPFTB(LL,1,2)),PORBEDP(LL,1),SEDPOC(LL,1)
       ENDIF  
       ENDIF
                
       ENDDO
       ENDDO
  501 format(F7.2,8F12.4)
!-------------------------------------------------------------------------

        goto 330
        
        DO NT=1,NPCB
C
              ! fixed 2 layer model
        DO L=2,LA
         KBTP1=KBTP(L)+1
         ALOW(L,1)=0.
         CUPP(L,KBTP1)=0.
         DO K=1,KBTP(L)-1
 !        CUPP(L,K)=-DIFTOX(NT)*(PORBEDP(L,K)+PORBEDP(L,K+1))/
 !    $                     (HBEDP(L,K)+HBEDP(L,K+1))
 !        CUPP(L,K)=-DIFTOXB(NT)*(PORBEDP(L,K)+PORBEDP(L,K+1))/
 !    $                     (HBEDP(L,K)+HBEDP(L,K+1)) 
     
          CUPP(L,K)=-DIFTOX(NT)*(0.82+0.82)/
     $                     (HBEDP(L,K)+HBEDP(L,K+1)) 
         END DO
!         CUPP(L,KBTP(L))=-DIFTOX(NT)*PORBEDP(L,KBTP(L))
          CUPP(L,KBTP(L))=-DIFTOX(NT)*0.82
     &	   /HBEDP(L,KBTP(L))
         DO K=2,KBTP1
          ALOW(L,K)=CUPP(L,K-1)
         END DO
         DO K=1,KBTP(L)
          BMNN(L,K)=DELTI*HBEDP(L,K)/(1.-TOXPFTB(L,K,NT))
         END DO
         BMNN(L,KBTP1)=DELTI*DZC(1)*HP(L)/(1.-TOXPFTWW(L,NT))         !(1.-TOXPFTW(L,1,NT))        !(1.-TOXPFTW(L,1,NT)) is dissolved
         DO K=1,KBTP1
          BMNN(L,K)=BMNN(L,K)-ALOW(L,K)-CUPP(L,K)
         END DO
         DO K=1,KBTP(L)
          RRHS(L,K)=DELTI*TOXB(L,K,NT)
         END DO
         RRHS(L,KBTP1)=DELTI*DZC(1)*HP(L)*TOX(L,1,NT)
        END DO
C
        DO L=2,LA
         KBTP1=KBTP(L)+1
         BETTMP=BMNN(L,1)
         TOXTMP(L,1)=RRHS(L,1)/(BETTMP)  ! J.S.
         DO KK=2,KBTP1
          GAMTMP(L,KK)=CUPP(L,KK-1)/(BETTMP)
          BETTMP=BMNN(L,KK)-ALOW(L,KK)*GAMTMP(L,KK)
          TOXTMP(L,KK)=(RRHS(L,KK)-ALOW(L,KK)*TOXTMP(L,KK-1))
     &		/(BETTMP)
         END DO
         DO KK=KBTP(L),1,-1
          TOXTMP(L,KK)=TOXTMP(L,KK)-GAMTMP(L,KK+1)*TOXTMP(L,KK+1)
         END DO
        END DO
C
        DO L=2,LA
         KBTP1=KBTP(L)+1
         DO K=1,KBTP(L)
          TOXB(L,K,NT)=HBEDP(L,K)*TOXTMP(L,K)/(1.-TOXPFTB(L,K,NT))
         END DO
         TOX(L,1,NT)=TOXTMP(L,KBTP1)/(1.-TOXPFTWW(L,NT))     !(1.-TOXPFTW(L,1,NT))
        END DO
C
        END DO
500     continue
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT REACTIONS
C 
        DO NT=1,NPCB
C
C----------------------------------------------------------------------C
C
C **     NOTES:
C
C        BULK DECAY COEFFICIENT
C
C          RKTOXWT=RKTOXW(NT)    !*(TEM(L,K)-TKTOXW(NT))**ALPTOX
C          RKTOXBT=RKTOXB(NT)    !*(TEM(L,1)-TRTOXB(NT))**ALPTOX
C
C        VOLITIZATION
C
C          VOLTOX(NT)
C          RMOLTX(NT)=MOLECULAR WEIGHT
C
C        PHOTOLOSIS
C
C          RKTOXP(NT)=BASE RATE
C          SKTOXP(NT)=SOLAR RADIATION AT BASE RATE
C
         DO K=1,KC
         DO L=2,LA
          CDECAYW(L,K)=1./(1.+DELT*RKTOXW(NT))
         END DO
         END DO
C
         DO K=1,KC
         DO L=2,LA
          TOX(L,K,NT)=CDECAYW(L,K)*TOX(L,K,NT)
         END DO
         END DO
C
         DO K=1,KBP
         DO L=2,LA
          CDECAYB(L,K)=1./(1.+DELT*RKTOXB(NT))
         END DO
         END DO
C
         DO K=1,KBP
         DO L=2,LA
          TOXB(L,K,NT)=CDECAYB(L,K)*TOXB(L,K,NT)
         END DO
         END DO
C
C----------------------------------------------------------------------C
C
        END DO
C**********************************************************************C
C
C **  TOXIC CONTAMINANT REACTIONS SURFACE
C     AirDiff = diffusion
C     TOXS0   = area dissolved PCB
C
                
330	continue

      IWQKA=2
!	OXA0=0.5*1.0e-3   ! input 0.5 ng/m^3
	DK_h=35
!	che=3.92e-4 
!	chelog=log(che)/2.302585 !3.92e-4   ! (atm m^3 mol, for Kh=40 need change for each one)

      DO L=2,LA
 
         umrm = UWQS(L)
         vmrm = VWQS(L)
         wind_t=max(sqrt(umrm*umrm+vmrm*vmrm),0.2)
	   CALL doboundary(Os,Asss,TEM(L,KC),SAL(L,KC),WINDST(L),
     &       umrm,vmrm,
	1   HP(L),DZC(KC),0,IWQKA)
	   AirD=AirDiff*Asss/86400.0	   
!	    if(L.EQ.100) write(*,*)AirDiff*86400.,Gastox,
!     * TOXPFTW(L,KC,1)*TOX(L,KC,1)
       DO NT=1,NPCB 

        chelog=HENC(NT) ! log(HENC(NT))/2.302585
        che=10**(chelog+2.3+26.39-7868/(273.15+TEM(L,KC)))   ! Nelson's eq.F-14a       
  
        ckd=1-(TOXPFW(L,KC,1,NT)+TOXPFW(L,KC,2,NT))/(1+TOXPFTW(L,KC,NT))       ! dissolved fraction
        ckg=168*wind_t*(18.0/RMOLTX(NT))**0.25                !Eq. F-12
        ckg=ckg/86400.0
        ckl=AirD*(32.0/RMOLTX(NT))**0.25                      !Eq. F-10
 	  Gastox=RGASA(NT)*TOXA0*8.206e-5*(273+TEM(L,KC)) /che *1.0e-6 !Gas*RT/Kh   change pg/m^3 to ug/M^3
 	  AirD=ckl*che/(che+8.206e-5*(273+TEM(L,KC-1))*ckl/ckg) 
 	  TAMP=TOX(L,KC,NT) 
! Implicity
        CDECAYW(L,KC)=1./(1.+DELT*AirD*ckd )
        TOX(L,KC,NT)=CDECAYW(L,KC)*(TOX(L,KC,NT)+DELT*AirD*Gastox
     &   +DELT*RGASA(NT)*Atmdep/(HP(L)*DZC(KC)))
       
        VOLBW3(L,1)=VOLBW3(L,1)+SCB(L)*DXYP(L)*AirD                 ! exchange from air
     &     *(Gastox-TOX(L,KC,NT)*ckd)*HP(L)*DZC(KC)
 !    &      +Atmdep*DXYP(L)*RGASA(NT)
        GASC(L,NT)=Gastox
        ADEP=ADEP+ Atmdep*DXYP(L)*RGASA(NT)*DELT  
  !      CDECAYW(L,KC)=1./(1.+0.5*DELT*AirD*ckd )
  !      TOX(L,KC,NT)=CDECAYW(L,KC)*(TOX(L,KC,NT)+DELT*AirD*
  !   $   (Gastox-0.5*TOX(L,KC,NT)*ckd +DELT*Atmdep/(HP(L)*DZC(KC)) ))
  !      VOLBW3(L,1)=VOLBW3(L,1)+SCB(L)*DXYP(L)*AirD                 ! exchange from air
  !   &        *(Gastox-0.5*(TOX(L,KC,NT)+TAMP)*ckd)
  !   $        *HP(L)*DZC(KC)+Atmdep*DXYP(L)

C  Explicity     
C         TOX(L,KC,NT)=TOX(L,KC,NT)+
C     $    DELT*AirD*(Gastox-TAMP*ckd)   
C         VOLBW3(L,1)=VOLBW3(L,1)+SCB(L)*DXYP(L)*
C     $    AirD*(Gastox-TAMP*ckd)
       END DO
      ENDDO 
C
C
C--------------------
!      IPCON=2
      if(mod(N,ISVIN).eq.0.and.IPCON.gt.0)then
      
      if(II_binary.eq.0) then
      
       open(2,file='pcbtcon.out',ACCESS='APPEND')
       write(2,*)TIME
       DO L=2,LA
        DO K=1,KC
        TOX_TMP(K)=0
        ENDDO
        TOXB_TMP=0.
        TOXA_TMP=0.
        DO NT=1,NPCB
        DO K=1,KC
        TOX_TMP(K)=TOX_TMP(K)+TOX(L,K,NT)
        TOXA_TMP=TOXA_TMP+TOX(L,K,NT)
        ENDDO
        TOXB_TMP=TOXB_TMP+TOXB(L,KBP,NT)*TOXPFTB(L,KBP,NT)/HBEDP(L,KBP)/
     &    (2500.0*(1-PORBEDP(L,KBP)))   ! model unit ug/m^2 convert tot ug/kg  (ng/g)       
        ENDDO
        write(2,202)L,(TOX_TMP(K),K=1,KC),TOXB_TMP,TOXA_TMP/KC

       ENDDO
       close(2)
 202   format(I7,54F9.2) 
C 

          
       open(2,file='pcbcon1.out',ACCESS='APPEND')
       write(2,*)TIME
C
       IF(IPCON.EQ.1) THEN
       NT=1
       DO L=2,LA
       ckdb=1-(TOXPFW(L,1,1,NT)+TOXPFW(L,1,2,NT))/(1+TOXPFTW(L,1,NT))        
       ckds=1-(TOXPFW(L,KC,1,NT)+TOXPFW(L,KC,2,NT))/(1+TOXPFTW(L,KC,NT))
       write(2,201)IL(L),JL(L),TOX(L,K,NT),TOX(L,KC,NT),
     &    ckdb,ckds,GASC(L,NT)
     & ,TOXPFW(L,KC,1,NT),TOXPFW(L,KC,2,NT),TOXPFTW(L,KC,NT)
       ENDDO
       ELSE
       NT=1
       DO L=2,LA
       write(2,202)L,(TOX(L,K,NT),K=1,KC),
     $  TOXB(L,KBP,NT)*TOXPFTB(L,KBP,NT)/HBEDP(L,KBP)/
     &    (2500.0*(1-PORBEDP(L,KBP))) ,
     &  TOXB(L,KBP-1,NT)*TOXPFTB(L,KBP-1,NT)/HBEDP(L,KBP)/
     &    (2500.0*(1-PORBEDP(L,KBP))),
     & (TOXPFTW(L,K,NT),K=1,KC),(TOXPFTB(L,K,NT),K=1,2)
       ENDDO
       ENDIF
       close(2)
        
       open(2,file='pcbcon2.out',ACCESS='APPEND')
       write(2,*)TIME
       NT=2
       IF(IPCON.EQ.1)THEN
       DO L=2,LA
       ckdb=1-(TOXPFW(L,1,1,NT)+TOXPFW(L,1,2,NT))/(1+TOXPFTW(L,1,NT))        
       ckds=1-(TOXPFW(L,KC,1,NT)+TOXPFW(L,KC,2,NT))/(1+TOXPFTW(L,KC,NT))
       write(2,201)IL(L),JL(L),TOX(L,1,NT),TOX(L,KC,NT),
     &     ckdb,ckds,GASC(L,NT) 
       ENDDO
       ELSE
       NT=2
       DO L=2,LA
       write(2,202)L,(TOX(L,K,NT),K=1,KC),
     $  TOXB(L,KBP,NT)*TOXPFTB(L,KBP,NT)/HBEDP(L,KBP)/
     &    (2500.0*(1-PORBEDP(L,KBP))) ,
     &  TOXB(L,KBP-1,NT)*TOXPFTB(L,KBP-1,NT)/HBEDP(L,KBP-1)/
     &   ( 2500.0*(1-PORBEDP(L,KBP))),
     & (TOXPFTW(L,K,NT),K=1,KC),(TOXPFTB(L,K,NT),K=1,2)
       ENDDO
       ENDIF
       close(2) 
       
       open(2,file='pcbcon3.out',ACCESS='APPEND')
       write(2,*)TIME      
       IF(IPCON.EQ.1)THEN      
       NT=3
       DO L=2,LA
       ckdb=1-(TOXPFW(L,1,1,NT)+TOXPFW(L,1,2,NT))/(1+TOXPFTW(L,1,NT))        
       ckds=1-(TOXPFW(L,KC,1,NT)+TOXPFW(L,KC,2,NT))/(1+TOXPFTW(L,KC,NT)) 
       write(2,201)IL(L),JL(L),TOX(L,1,NT),TOX(L,KC,NT),
     &    ckdb,ckds,GASC(L,NT) 
       ENDDO
       ELSE
       NT=3
       DO L=2,LA
       write(2,202)L,(TOX(L,K,NT),K=1,KC),
     $  TOXB(L,KBP,NT)*TOXPFTB(L,KBP,NT)/HBEDP(L,KBP)/
     &    (2500.0*(1-PORBEDP(L,KBP))) ,
     &  TOXB(L,KBP-1,NT)*TOXPFTB(L,KBP-1,NT)/HBEDP(L,KBP-1)/
     &   ( 2500.0*(1-PORBEDP(L,KBP))),
     & (TOXPFTW(L,K,NT),K=1,KC),(TOXPFTB(L,K,NT),K=1,2)   
       ENDDO
       ENDIF
       close(2) 
            
201   format(2I6,4F9.3,E12.4,3F8.3)

      elseif(II_binary.eq.2) then
         
       write(801)TIME
C
 
       NT=1
       DO L=2,LA
       write(801)L,(TOX(L,K,NT),K=1,KC),      
     $  TOXB(L,KBP,NT)*TOXPFTB(L,KBP,NT)/HBEDP(L,KBP)/
     &    (2500.0*(1-PORBEDP(L,KBP))),      
     &  TOXB(L,KBP-1,NT)*TOXPFTB(L,KBP-1,NT)/HBEDP(L,KBP)/
     &    (2500.0*(1-PORBEDP(L,KBP))),
     & (TOXPFTW(L,K,NT),K=1,KC),(TOXPFTB(L,K,NT),K=1,2)
       ENDDO   
       IF(IPCON.GT.1)THEN
       write(802)TIME
       NT=2
       DO L=2,LA
       write(802)L,(TOX(L,K,NT),K=1,KC),
     $  TOXB(L,KBP,NT)*TOXPFTB(L,KBP,NT)/HBEDP(L,KBP)/
     &    (2500.0*(1-PORBEDP(L,KBP))) ,
     &  TOXB(L,KBP-1,NT)*TOXPFTB(L,KBP-1,NT)/HBEDP(L,KBP-1)/
     &   ( 2500.0*(1-PORBEDP(L,KBP))),
     & (TOXPFTW(L,K,NT),K=1,KC),(TOXPFTB(L,K,NT),K=1,2)
       ENDDO
       ENDIF
       IF(IPCON.GT.2)THEN              
       write(803)TIME      
       NT=3
       DO L=2,LA
       write(803)L,(TOX(L,K,NT),K=1,KC),
     $  TOXB(L,KBP,NT)*TOXPFTB(L,KBP,NT)/HBEDP(L,KBP)/
     &    (2500.0*(1-PORBEDP(L,KBP))) ,
     &  TOXB(L,KBP-1,NT)*TOXPFTB(L,KBP-1,NT)/HBEDP(L,KBP-1)/
     &   ( 2500.0*(1-PORBEDP(L,KBP))),
     & (TOXPFTW(L,K,NT),K=1,KC),(TOXPFTB(L,K,NT),K=1,2)   
       ENDDO
       ENDIF
       
       
       else
          
       write(801)TIME,TOX
       write(802)TIME,TOXB
       write(803)TIME,TOXPFW,HBEDP
       write(804)TIME,TOXPFTB,TOXPFTW,PORBEDP
     
      endif
      endif
C
C **  CALCULATE TOTAL PARTICULATE FRACTIONS IN WATER COLUMN FOR output
C
C
      DO NT=1,NPCB
       DO K=1,KC
        DO L=2,LA
          TOXPFTW(L,K,NT)=TOXPFTW(L,K,NT)/(1.+TOXPFTW(L,K,NT))
        END DO
       END DO
      END DO
C**********************************************************************C
C
910      CLOSE(1)
C
C**********************************************************************C
C
      RETURN
      END

