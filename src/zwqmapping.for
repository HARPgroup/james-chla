
C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WQMAPPING 
C
C**********************************************************************C
C
C  Mapping 2D arry to 1D arry. This subroutine needs to be revides to 
C  line CH3D or other model.
C  
C  LIJW(I,J)  : 1D L value corresponding to (i,j). Active water cell start at 2
C  ILW(L) and JLW(L): returns the i and j values corresponding to L
C  DXDYPW(L)  : surface area of cell L
C  SEDTWQ(L,K): sediment concentration
C  TEMWQ(L,K) : Temperature
C  HPWQ(L)    : Total depth at L
C  SALWQ(L,K) : salinity
C  DZCWQ(K)   : layer thickness
C  UWQS(L)   : Surface U velocity
C  VWQS(L)   : Surface V velocity
C  WINDSTWQ(L): wind speed
C  SCBWQ(L)   : = 1 open dary 
C               = 0 other cell
C  WQSEDO(I)  : constant sediment settling velocity
C  RAINTWQ(L) : rainfull (for wet depsition
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      include 'wq.par'
      INCLUDE 'efdc.cmn'
      INCLUDE 'wqcom.cmn'
      INTEGER, SAVE :: inimap = 1 
           
!      CALL OMP_SET_NUM_THREADs(iCPU)
 !     CALL OMP_SET_NUM_THREADs(4)

      IF(inimap.EQ.1) THEN
      inimap=0
      DO I=1,IC
	DO J=1,JC
	LIJW(I,J)=LIJ(I,J)
	ENDDO
	ENDDO

!$OMP          PARALLEL DO DEFAULT(SHARED)
	DO I=1,LA
	ILW(I)=IL(I)
	ENDDO

!$OMP          PARALLEL DO DEFAULT(SHARED)
	DO J=1,LA
	JLW(J)=JL(J)
	ENDDO
 
 !$OMP          PARALLEL DO DEFAULT(SHARED)       
	DO L=2,LA
      DXYPWQ(L)=DXYP(L)   ! for use in Atmospheric Deposition
     	SCBWQ(L)=SCB(L)        ! Open BC = 0
      ENDDO
      
      DO K=1,KC
      DZCWQ(K)=DZC(K)
	ENDDO  
	    
      ENDIF
 ! ----------------------------end initical--------------------------
      
      if(ifed_inc.GE.0)then
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)       
	 DO K=1,KC
	 DO L=2,LA
	  SEDTWQ(L,K)=SEDT(L,K)
       ENDDO
	 ENDDO
       if(ISTRAN(6).eq.0) then
 	  DO K=1,KC
	  DO L=2,LA
	   SEDTWQ(L,K)=15
        ENDDO
	  ENDDO     
       endif
      endif
      
	IF(IFOCETMP.EQ.1) THEN
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)    
	DO K=1,KC
	DO L=2,LA
	TEMWQ(L,K)=(TPMAG)*abs(cos(3.1415/365*(DT*N/86400-TPHPSE)))+TMPADD
	TEM(L,K)=TEMWQ(L,K)
	TEM1(L,K)=TEMWQ(L,K)
	TEMWQ(L,K)=min(TEM(L,K),40.0)
	TEMWQ(L,K)=max(TEM(L,K),-4.0)
      ENDDO
	ENDDO
!      write(*,*)'Call simple temp calculation '	,TEMWQ(10,KC)
	ELSE
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
	 DO K=1,KC
	 DO L=2,LA	
       TEMWQ(L,K)=min(TEM(L,K),40.0)
	 TEMWQ(L,K)=max(TEM(L,K),-4.0)
       ENDDO
	 ENDDO
C
	ENDIF
	
      if(ifed_inc.GE.0)then
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
	DO K=1,KC
	DO L=2,LA
	SALWQ(L,K)=SAL(L,K)
      ENDDO
	ENDDO
      DO L=1,LA
      HPWQ(L)=HWQ(L)
! 	HPWQ(L)=HP(L)  
!	UWQS(L)=max(U(L,KC), U(L+1,KC))
!	VWQS(L)=max(V(L,KC), V(LNC(L),KC))  ! surface velocity 	   
	UWQS(L)=max(UWQ(L,KC), UWQ(L+1,KC))
	VWQS(L)=max(VWQ(L,KC), VWQ(LNC(L),KC))  ! surface velocity 		
      ENDDO
      endif
      
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)     
      DO L=1,LA
!	DXYPWQ(L)=DXYP(L)      ! surface area
	WINDSTWQ(L)=WINDST(L)  ! Wind
!	SCBWQ(L)=SCB(L)        ! Open BC = 0
	ENDDO
	
      DO K=1,KC
      DZCWQ(K)=DZC(K)
	ENDDO

      WSEWQ=WSEDO(1)

C**********************************************************************C
C
C **  SET WEIGHTS FOR SALINITY AND TEMPERATURE BOUNDARY INTERPOLATION
C
C----------------------------------------------------------------------C
C
      IF (KC.GT.1) THEN
      DO K=1,KC
      WTCII(K,1)=FLOAT(K-KC)/FLOAT(1-KC)
      WTCII(K,2)=FLOAT(K-1)/FLOAT(KC-1)
      END DO
      ELSE
      WTCII(1,1)=0.5 
      WTCII(1,2)=0.5
      END IF

!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
	DO L=2,LA
	RAINTWQ(L)=raint(L)   ! Input rain fall 
	ENDDO
C
      DO I=1,NSTM
      WQSEDO(I)=WSEDO(I)
	ENDDO

	RETURN
      END



