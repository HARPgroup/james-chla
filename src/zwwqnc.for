C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WWQNC(LA,KC)
C
C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 10 april 1999
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
!      OPEN(1,FILE=NCOFN,STATUS='UNKNOWN',ACCESS='APPEND')
      OPEN(1,FILE=NCOFN,STATUS='UNKNOWN')
C
      DATA WQVN/
     * 'Bc   ','Bd   ','Bg   ','RPOC ','LPOC ','DOC  ','RPOP ','LPOP ',
     * 'DOP  ','PO4t ','RPON ','LPON ','DON  ','NH4  ','NO3  ','SU   ',
     * 'SA   ','COD  ','O2   ','TAM  ','FCB  ','MALG '/
C
       NW1=0
       LL=0
      IF (IWQNC.GE.1) THEN  
       DO L=2,LA
        DO K=1,KC
         DO NW=1,NWQV
           IF (WQV(L,K,NW).LT.0.0) then
!		      WRITE(1,90) WQVN(NW),
!     *        ITNWQ,L,ILW(L),JLW(L),K,WQV(L,K,NW)
	        WQV(L,K,NW)=0     
	        NW1=NW1+1   
	        LL=L  
           ENDIF
		 END DO
        END DO
       END DO  
      ENDIF
C
      if(NW1.GT.0)WRITE(1,*) LL,NW1
      CLOSE(1)
C
   90 FORMAT(A5, I8, 4I5, E11.3)
C
      RETURN
      END
