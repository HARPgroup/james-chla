        SUBROUTINE ZDYE_inv
       !(ISTL,N,DT,TCON,TBEGIN,TIDALP,
       ! &            NTSPTC,IWQDT,LA,KC,IC,JC,IWQS,IOP_SAVE,iyear1)
C------------------------------------------------------------------------     
C Comput dye or sedimetn using saved dynamic fields
C  J. Shen 3/16/2014
C  NSTATEV = number of variable used for dye
C  Use WQ  input as boundary condition 
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn' 
      INCLUDE 'wq.par' 
      INCLUDE 'wqcom.cmn'  
      INTEGER, SAVE :: inidy=1, NSTATEV=2

      DELT=DT2
      DELTA=DT2
      DELTD2=DT    
      
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON  
               
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
      DO L=2,LA
      HWQI(L)=1./HWQ(L)
      END DO

! **  UPDATE TIME SERIES BOUNDARY CONDITIONS
!
      CALL RWQCSR(KC,DT,N,TCON,TBEGIN,NTRNVA,NTSPTC,NCSTEP,SECDLAST) 
      
      ITNWQ = ITNWQ + IWQDT
      
!      CALL WQMAPPING

!      CALL CALWQCN(2,WQV,LCMWQ,KCWM,NWQVM,IWQFCB,NWQV,IWTRC,F_AB)  
      
      DO NW=1,NTRNVA
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
       
         DO K=1,KC
         DO L=2,LA
          CWQ(L,K)=WQV(L,K,NW)
          CWQ2(L,K)=WQV(L,K,NW)
         END DO
         END DO
!  
!     Transport
!       
	   CALL CALTRWQG (8,NW,CWQ,CWQ2)
	   
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
          DO K=1,KC
	    DO L=2,LA
	     WQV(L,K,NW)=CWQ(L,K)         
          END DO
          END DO         


!
! VERTICAL DIFFUSION CALCULATION
!
        IF (KC.EQ.1) GO TO 2001
!
        RCDZKK=-DELT*CDZKK(1)
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,CCUBTMP,CCMBTMP,EEB) 
        DO L=2,LA
         CCUBTMP=RCDZKK*HWQI(L)*AB(L,1)*F_AB
         CCMBTMP=1.-CCUBTMP
         EEB=1./CCMBTMP
         CU1(L,1)=CCUBTMP*EEB
         WQV(L,1, NW)=WQV(L,1, NW)*EEB
        END DO
C
        DO K=2,KS
         RCDZKMK=-DELT*CDZKMK(K)
         RCDZKK=-DELT*CDZKK(K)
         DO L=2,LA
          CCLBTMP=RCDZKMK*HWQI(L)*AB(L,K-1)*F_AB
          CCUBTMP=RCDZKK*HWQI(L)*AB(L,K)*F_AB
          CCMBTMP=1.-CCLBTMP-CCUBTMP
          EEB=1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
          CU1(L,K)=CCUBTMP*EEB
          WQV(L,K, NW)=(WQV(L,K, NW)-CCLBTMP*WQV(L,K-1, NW))*EEB
         END DO
        END DO
        K=KC
        RCDZKMK=-DELT*CDZKMK(K)
        DO L=2,LA
         CCLBTMP=RCDZKMK*HWQI(L)*AB(L,K-1)*F_AB
         CCMBTMP=1.-CCLBTMP
         EEB=1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
         WQV(L,K, NW)=(WQV(L,K, NW)-CCLBTMP*WQV(L,K-1, NW))*EEB
        END DO

        DO K=KC-1,1,-1
         DO L=2,LA
          WQV(L,K, NW)=WQV(L,K, NW)-CU1(L,K)*WQV(L,K+1, NW)
         END DO
        END DO

2001   continue
 !      if(mod(N,100).EQ.0) write(*,*)'Dye = ',WQV(306,1,1)
 
      ENDDO
      
      fct=400
      if(TIME.LT.fct)fct=TIME
      DO NW=1,NTRNVA
      
       DO K=1,KC
        DO L=2,LA
         WQV(L,K,NW)=WQV(L,K,NW)+abs(DELT)/86400*WQ3DA(L,K,NW)
         if(WQV(L,K,NW).GT.fct) THEN
 !         write(500,'(I8,I5,F10.2,E12.4)')L,K,TIME,WQV(L,K,1)
          WQV(L,K,NW)=fct
 !         write(500,'(I8,I5,F10.2,E12.4)')L,K,TIME,WQV(L,K,1)     
         endif
         if(WQV(L,K,NW).LT.0.0)  WQV(L,K,NW)=0
        END DO
       END DO
       
      ENDDO
 !      L1=LIJ(160,94)
 !      L2=LIJ(160,95)
 !      DO K=1,KC
 !      WQV(L1,K,1)=0
 !      WQV(L2,K,1)=0
 !      ENDDO
 !      L1=LIJ(161,94)
 !      L2=LIJ(161,95)
 !      DO K=1,KC
 !      WQV(L1,K,1)=0
 !      WQV(L2,K,1)=0
 !      ENDDO
       
 !        DO L=2,LA
 !          WQV(L,KC,1)=0
 !        ENDDO
               
      RETURN
      END