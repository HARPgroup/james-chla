C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WWQNC2(LA,KC)   ! for checking negatives, Jeff Ji, 7/29/99
C
C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK ON  7 APRIL 1997
C
C**********************************************************************C
C
C Write information of negative WQ state variables (unit IWQONC).
C
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
C
      CHARACTER*5 WQVN(22)
C
       OPEN(1,FILE=NCOFN,STATUS='UNKNOWN',ACCESS='APPEND')
C      OPEN(1,FILE=NCOFN,STATUS='UNKNOWN')      
C
      DATA WQVN/
     * 'Bc   ','Bd   ','Bg   ','RPOC ','LPOC ','DOC  ','RPOP ','LPOP ',
     * 'DOP  ','PO4t ','RPON ','LPON ','DON  ','NH4  ','NO3  ','SU   ',
     * 'SA   ','COD  ','O2   ','TAM  ','FCB  ','MALG '/
C
      KK=0
      KK1=0
      IF (IWQNC.GE.1) THEN
       DO L=2,LA
        DO K=1,KC
         DO NW=1,NWQV
            IF (WQV(L,K,NW).LT.0.0) then
		     if(mod(kk,100).eq.0)WRITE(1,90) 'TRN0 ',WQVN(NW),
     *        ITNWQ,L,ILW(L),JLW(L),K,WQV(L,K,NW)
             wqv(L,k,nw)=0.0  ! force it to zero 
             KK=KK+1   
            ENDIF	
	      IF (IWQNC.GT.1) THEN	            
	       IF(abs(WQV(L,K,NW)).GT.IWQNC) then
              if(mod(KK1,100).eq.0)WRITE(1,90) 'TRN+ ',
     &         WQVN(NW),ITNWQ,L,ILW(L),JLW(L),K,WQV(L,K,NW)
               WQV(L,K,NW)=IWQNC !real(IWQNC)/10.0
               KK1=KK1+1
             ENDIF
	      ENDIF
         END DO
        END DO
       END DO
      ENDIF
C
      CLOSE(1)
C
   90 FORMAT(A5,A5, I8, 4I5, E11.3, "   WWQNC2")
C
      RETURN
      END
