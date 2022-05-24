C Preprocess of ADV
C
      SUBROUTINE PREADV
      INCLUDE 'efdc.par'
      INCLUDE 'wq.par'
      INCLUDE 'efdc.cmn'
C
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
        IF(NCSTEP.GT.0) TIME=(SECDLAST+TCON*TBEGIN)/86400. !% J.S. 1/31/2014
C
C **  ACCUMULTATE COURANT NUMBERS (CFL)
C
      DELT=DT2
      CCFLUM=10000
      CCFLVM=10000
      
      DO K=1,KC
      DO L=2,LA
        CKAX(L,K)=1.-DELT*ABS(DXIU(L)*UWQ(L,K))
        CKAY(L,K)=1.-DELT*ABS(DYIV(L)*VWQ(L,K))
        CKAW(L,K)=1.-DELT*ABS(HPI(L)*DZIG(K)*WWQ(L,K))
        if(CKAX(L,K).LT.0.01) then
        CCFLUM=min(CCFLUM,CKAX(L,K))
        II_LL=IL(L)
        JJ_LL=JL(L)
 !       write(*,*)'CFL condition U',IL(L),JL(L),UWQ(L,K),CKAX(L,K)
 !       CKAX(L,K)=0.2
        endif  
        if(CKAY(L,K).LT.0.01) then
        CCFLVM=min(CCFLVM,CKAY(L,K))
        II_LLv=IL(L)
        JJ_LLv=JL(L)
 !       write(*,*)'CFL condition V',IL(L),JL(L),K,VWQ(L,K),CKAY(L,K)
 !       CKAY(L,K)=0.2
        endif       
        BETAX1(L,K)= 0.25*CKAX(L,K)+0.5-1./(12.*CKAX(L,K))      !A1
        BETAX2(L,K)=-0.25*CKAX(L,K)+0.5+1./(12.*CKAX(L,K))      !B1
        BETAY1(L,K)= 0.25*CKAY(L,K)+0.5-1./(12.*CKAY(L,K))      
        BETAY2(L,K)=-0.25*CKAY(L,K)+0.5+1./(12.*CKAY(L,K))  
        BETAZ1(L,K)= 0.25*CKAW(L,K)+0.5-1./(12.*CKAW(L,K))      
        BETAZ2(L,K)=-0.25*CKAW(L,K)+0.5+1./(12.*CKAW(L,K))  

      END DO
      END DO
      if(CCFLUM.LT.0.or.CCFLVM.LT.0)then
      if(II_LL.NE.332) then
      write(500,*)II_LL,JJ_LL,TIME,CCFLUM
!      write(500,*)II_LLv,JJ_LLv,TIME,CCFLVM  
      endif   
      endif
C
      RETURN
      END