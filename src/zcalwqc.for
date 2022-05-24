C 
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALWQCN(ISTL,WQV,LCMWQ,KCWM,NWQVM,IWQFCB,NWQV,
     &           IWTRC,F_AB)  
C
C     Trasport subrouting called by water quality model
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE CALWQC CALCULATES THE CONCENTRATION OF DISSOLVED AND 
C **  SUSPENDED WATER QUALITY CONSTITUTENTS AT TIME LEVEL (N+1). 
C **  CALLED ONLY ON ODD THREE TIME LEVEL STEPS
C
!  B  B  B  R  L  D  R  L  D  P  R  L  D  N  N  S  S  C  D  T  F  M
!  c  d  g  P  P  O  P  P  O  O  P  P  O  H  O  U  A  O  O  A  C  A
!           O  O  C  O  O  P  4  O  O  N  4  3        D     M     C
!           C  C     P  P     t  N  N                             A
!  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
C**********************************************************************C
C
      INCLUDE 'efdc.par'  
      INCLUDE 'efdc.cmn'
 !     INTEGER FUNCTION OMP_GET_NUM_THREADS         
	DIMENSION WQV(LCMWQ,KCWM,0:NWQVM)
      DIMENSION ICON_T(21),IWTRC(21)
      DATA ICON_T/1,2,3,0,5,6,0,8,9,10,0,12,13,14,15,0,17,18,19,0,0/      
!      DATA ICON_T/1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/  
!      INCLUDE 'wq.par'
!      INCLUDE 'wqcom.cmn'
C
C**********************************************************************C
C
   !   CALL OMP_SET_NUM_THREADs(4) 
      DELT=DT2
      IF (ISADAC(8).EQ.2)THEN
      HSIMT=1
      ELSE
      HSIMT=0      
      ENDIF
      
C
C**********************************************************************C
C
C **  UPDATED TIME SERIES CONCENTRATION BOUNDARY CONDITIONS
C
C     CALL CALWQS(ISTL)
C      
C**********************************************************************C
C
C **  3D ADVECTI0N TRANSPORT CALCULATION 
C

      IF(HSIMT.EQ.1) CALL PREADV
      
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
      DO L=2,LA
      HWQI(L)=1./HWQ(L)
      END DO

      NVQVT=NWQV
      IF (IWQFCB.EQ.0) NVQVT=NWQV-1
C

      DO NW=1,NVQVT
       if(ICON_T(NW).EQ.0) goto 2001
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
!	 
!         IF(HSIMT.EQ.0)THEN
         CALL CALTRWQ (8,NW,CWQ,CWQ2)
!         ELSE
!	   CALL CALTRWQH (8,NW,CWQ,CWQ2)
!	   ENDIF
	   
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
          DO K=1,KC
	    DO L=2,LA
	     WQV(L,K,NW)=CWQ(L,K)         
          END DO
          END DO       
          
          
C        IF(IWQNC.GE.2) THEN
!          A_tmp=1.0e5
!          DO K=1,KC
!	    DO L=2,LA
!	    if(A_tmp>CWQ(L,K))A_tmp=CWQ(L,K )
!	    ENDDO
!	    ENDDO
!	    if( A_tmp.GT.0.0) THEN    
!          DO K=1,KC
!	    DO L=2,LA
!	     WQV(L,K,NW)=CWQ(L,K)         
!          END DO
!          END DO
!          Else
!          if(abs(A_tmp).LT.0.1) THEN 
!           DO K=1,KC
!	    DO L=2,LA
!	     WQV(L,K,NW)=CWQ(L,K)         
!          END DO
!          END DO
!          endif         
!          ENDIF

C         ELSE
C!$OMP          PARALLEL DO DEFAULT(SHARED)
C!$OMP&         PRIVATE(L) 
C          DO K=1,KC
C	    DO L=2,LA
C	     IF(CWQ(L,K).GE.0) WQV(L,K,NW)=CWQ(L,K)         
C          END DO
C          END DO
C        ENDIF         
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

      ENDDO
 
C**********************************************************************C
C
      RETURN
      END
