C                                                                                                    
C**********************************************************************C                            
C**********************************************************************C                            
C**********************************************************************C                            
C                                                                                                   
      SUBROUTINE RSMRST(LA)                                                                             
C                                                                                                   
C**********************************************************************C                            
C                                                                                                   
C **  LAST MODIFIED BY JOHN HAMRICK ON  7 APRIL 1997                                                
C                                                                                                   
C**********************************************************************C                            
C                                                                                                   
C Read ICs from restart file from INSMRST.                                                          
C                                                                                                   
C**********************************************************************C                            
C                                                                                                   
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'                                                                          
      logical fexist                                                                                
      integer nn                                                                                    
      real xtime                                                                                    
C                                                                                                   
c Check first to see if binary restart file exists.  If not, use                                    
c the ASCII file instead.                                                                           
c                                                                                                   
      INQUIRE(FILE='WQSDrst.bin', EXIST=fexist)                                                     
      if (.not. fexist) then                                                                        
        OPEN(1,FILE='wqsdrst.inp',STATUS='UNKNOWN')                                                 
        READ(1,999)                                                                                 
        READ(1,999)                                                                                 
C                                                                                                   
        DO M=2,LA                                                                                   
          READ(1,*) L,(SMPON(L,NW),NW=1,NSMG),                                                      
     *      (SMPOP(L,NW),NW=1,NSMG),(SMPOC(L,NW),NW=1,NSMG),SM1NH4(L),                              
     *      SM2NH4(L),SM2NO3(L),SM2PO4(L),SM2H2S(L),SMPSI(L),SM2SI(L),                              
     *      SMBST(L),SMT(L)                                                             
        END DO                                                                                      
c Bug: SOD (=wqbfo2) is discontinous when hot start, jeff ji, 11/1/99                               
C                                                                                                   
        CLOSE(1)                                                                                    
      else                                                                                          
        open(UNIT=1, FILE='WQSDrst.bin', ACCESS='append',                                           
     +     FORM='unformatted', STATUS='unknown')                                                    
        read(1) nn, xtime                                                                           
        xtime=xtime                                                                                 
        write(0,911) nn, xtime                                                                      
911     format(' Reading binary WQSDRST.BIN file ...    NN, TIME = ',                               
     +      I7, F11.5)                                                                              
        DO M=2,LA                                                                                   
          READ(1) L                                                                                 
          READ(1) (SMPON(L,NW),NW=1,NSMG),                                                          
     *      (SMPOP(L,NW),NW=1,NSMG),(SMPOC(L,NW),NW=1,NSMG),SM1NH4(L),                              
     *      SM2NH4(L),SM2NO3(L),SM2PO4(L),SM2H2S(L),SMPSI(L),SM2SI(L),                              
     *      SMBST(L),SMT(L)                                                                         
        END DO                                                                                      
        CLOSE(1)                                                                                    
      end if                                                                                        
C                                                                                                   
   90 FORMAT(I5, 18E12.4)                                                                           
  999 FORMAT(1X)                                                                                    
C                                                                                                   
      RETURN                                                                                        
      END                                                                                           
