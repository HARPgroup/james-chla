C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE BUDGET5
C
C **  ADDED BY DON KINGERY ON 15 OCTOBER 1996
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
C **  SUBROUTINES BUDGETn CALCULATE SEDIMENT BUDGET (TOTAL SEDIMENTS)
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      
C**********************************************************************C
C
C **  CALCULATE ENDING SUSPENDED AND BOTTOM SEDIMENT IN THE MODEL DOMAIN
C
C----------------------------------------------------------------------C
      TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
C
      SDFLUX=0.
      AIRFUX=0.
      SSEDEND=0.
      BSEDEND=0.
      VOLMEND=0.
      SMASSEND=0.
      TOXBEDTOP=0.
      SEDIN1=0.
c
       DO L=2,LA
        VOLMEND=VOLMEND+SPB(L)*DXYP(L)*HWQ(L)
       END DO
C
C  Initialize bottom and suspended sediment mass  ---  DLK 9/26
C
      DO L=2,LA
       SDFLUX=SDFLUX+SCB(L)*VOLBW3(L,2)
       AIRFUX=AIRFUX+SCB(L)*VOLBW3(L,1)
      END DO
C
      DO K=1,KB
      DO L=2,LA
        SEDBT(L,K)=0.
        SNDBT(L,K)=0.
      END DO
      END DO
C
      DO NS=1,NPCB
       DO K=1,KB
       DO L=2,LA
        BSEDEND=BSEDEND+SCB(L)*DXYP(L)*TOXB(L,K,NS)                    !bottom PCB in ug/m^2
       END DO
       END DO 
       DO L=2,LA
        TOXBEDTOP=TOXBEDTOP+SCB(L)*DXYP(L)*TOXB(L,KB,NS)  
       ENDDO
      END DO
      
C
!      DO K=1,KB
!      DO L=2,LA
!       BSEDEND=BSEDEND+SCB(L)*DXYP(L)*(SEDBT(L,K) )           !bottom PCB in ug
!      END DO
!      END DO
C

      DO NS=1,NPCB
       DO K=1,KC
        DO L=2,LA
         SSEDEND=SSEDEND+SCB(L)*DXYP(L)*HWQ(L)*TOX(L,K,NS)*DZC(K)  !bottom PCB in ug
        END DO
       END DO
      END DO

      DT2=2*DT    ! All budget was called at even number.
      
      SDFLUX=DT2*SDFLUX    ! flux from sediment
      AIRFUX=DT2*AIRFUX    ! flux from air
      SEDOUT1=DT2*SEDOUT   ! flux from boundary     >0 out
      SEDIN1=DT2*SEDIN     ! load from point source >0 in
      TOXOUTBC=-DT2*TOXOUTBC  ! note out >0 , revers sign here
      TOXINBC =-DT2*TOXINBC   ! note in <0  , revers sign here

!      VOLMOUT=DT2*VOLMOUT
!      VOLMIN=DT2*VOLMIN
!      SMASSIN=  Diffusion between w&s
C
      S_NETW=SEDIN1-SEDOUT1 + SDFLUX+ AIRFUX +SMASSIN*DT2 +ADEP      ! SEDOUT1> 0 out, SMASSIN >0 flux to water)      
      S_NETW1=SSEDBEG+S_NETW
      S_NETS=-SDFLUX-SEDERR*DT2-SMASSIN*DT2 
      S_NETS1=BSEDBEG+S_NETS
      ERRO_W=S_NETW1-SSEDEND
      ERRO_E=S_NETS1-BSEDEND
      ERRO_WR=ERRO_W/SSEDBEG *100
      ERRO_ER=ERRO_E/BSEDBEG *100
      VOLMERR=VOLMEND-VOLMBMO
  
C
C**********************************************************************C
C
C **  OUTPUT BALANCE RESULTS TO FILE budget.out
C
C----------------------------------------------------------------------C
C

        OPEN(89,FILE='budget.out',ACCESS='APPEND',STATUS='UNKNOWN')
        write(89,*)'Time =              ',TIME
        write(89,*)'Beg in water & sed. ',SSEDBEG,BSEDBEG
        write(89,*)'End in water & sed. ',SSEDEND,BSEDEND
        write(89,*)'Net from boundary   ',-SEDOUT1
        write(89,*)'Input from boundary ',TOXINBC
        write(89,*)'Transport out at BC ',TOXOUTBC                
        write(89,*)'From Point source   ',SEDIN1 
        write(89,*)'From air exchange   ',AIRFUX
        write(89,*)'Atm. depostion      ',ADEP
        write(89,*)'From sediment       ',SDFLUX+SMASSIN*DT2
        write(89,*)'Bural               ',SEDERR*DT2
        write(89,*)'Diff from sedimen   ',SMASSIN*DT2
        write(89,*)'Settle to sediment  ',SDFLUX        !<0 settle to sed
        write(89,*)'Total to  w & sed.  ',S_NETW,S_NETS
        write(89,*)'Total tox on top bed',TOXBEDTOP,TOXBEDTOP
        write(89,*)'Err in water & sed. ',ERRO_W,ERRO_E
        WRITE(89,*)'Re Err in w & S %   ',ERRO_WR,ERRO_SR
        CLOSE(89)
  189   format(A20,2E15.8)
        OPEN(89,FILE='budgetime.out',ACCESS='APPEND',STATUS='UNKNOWN')       
        write(89,'(F10.2,22E12.4)') TIME,SSEDBEG,BSEDBEG,
     &   SSEDEND,BSEDEND,-SEDOUT1,AIRFUX,SDFLUX+SMASSIN*DT2,
     $   SEDERR*DT2,SMASSIN*DT2,SDFLUX,S_NETW,S_NETS,
     $   ERRO_W,ERRO_E,ERRO_WR,ERRO_SR,TOXINBC,TOXOUTBC,ADEP,TOXBEDTOP,
     &   SEDIN1 
        CLOSE(89)
C
      
      SSEDBEG=SSEDEND
      BSEDBEG=BSEDEND
      SDFLUX=0.
      AIRFUX=0.
      SEDOUT=0.
      SEDIN=0.
      SMASSIN=0.
      SEDERR=0.
      ADEP=0
      TOXOUTBC=0.
      TOXINBC=0. 
      DO L=2,LA
       VOLBW3(L,2)=0
       VOLBW3(L,1)=0
      END DO    
      
      RETURN
      END
