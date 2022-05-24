C                                                                                                   
C**********************************************************************C                            
C**********************************************************************C                            
C**********************************************************************C                            
C                                                                                                   
      SUBROUTINE WQzero(LA,KC)                                                                             
C                                                                                                   
C**********************************************************************C                            
C M. Morton  28 Jun 1998                                                                            
C Initializes the water quality averaging summation arrays:                                         
C**********************************************************************C                            
C                                                                                                   
      INCLUDE 'wq.par' 
      INCLUDE 'wqcom.cmn'                                                                          
                                                                              
      DO LL=2,LA                                                                                    
        DO K=1,KC                                                                                   
c            LL=LWQTS(M)                                                                            
          DO NW=1,NWQV                                                                                                                                             
            WQVmin(LL,K,NW) = 9.99e+21                                                              
            WQVmax(LL,K,NW) = 0.0                                                                   
          END DO                                                                                    
          TOCWQsum(LL,K) = 0.0                                                                      
          PO4DWQsum(LL,K) = 0.0                                                                     
          SADWQsum(LL,K) = 0.0                                                                      
          TPWQsum(LL,K) = 0.0                                                                       
          TNWQsum(LL,K) = 0.0                                                                       
          TSIWQsum(LL,K) = 0.0                                                                      
          BOD5sum(LL,K) = 0.0                                                                       
          CHLmsum(LL,K) = 0.0                                                                       
          POCsum(LL,K) = 0.0                                                                        
          POPsum(LL,K) = 0.0                                                                        
          PONsum(LL,K) = 0.0                                                                        
          TOCWQmin(LL,K) = 9.99e+21                                                                 
          TOCWQmax(LL,K) = 0.0                                                                      
          TNWQmin(LL,K) = 9.99e+21                                                                  
          TNWQmax(LL,K) = 0.0                                                                       
          TPWQmin(LL,K) = 9.99e+21                                                                  
          TPWQmax(LL,K) = 0.0                                                                       
          SADWQmin(LL,K) = 9.99e+21                                                                 
          SADWQmax(LL,K) = 0.0                                                                      
          POCmin(LL,K) = 9.99e+21                                                                   
          POCmax(LL,K) = 0.0                                                                        
          POPmin(LL,K) = 9.99e+21                                                                   
          POPmax(LL,K) = 0.0                                                                        
          PONmin(LL,K) = 9.99e+21                                                                   
          PONmax(LL,K) = 0.0                                                                        
          CHLmmin(LL,K) = 9.99e+21                                                                  
          CHLmmax(LL,K) = 0.0                                                                       
          SALsum(LL,K) = 0.0                                                                        
          SALmn(LL,K) = 9.99e+21                                                                    
          SALmx(LL,K) = 0.0                                                                         
          TSSsum(LL,K) = 0.0                                                                        
          TSSmn(LL,K) = 9.99e+21                                                                    
          TSSmx(LL,K) = 0.0                                                                         
          WQketsum(LL,K) = 0.0                                                                      
          WQketmn(LL,K) = 9.99e+21                                                                  
          WQketmx(LL,K) = 0.0                                                                       
        END DO                                                                                      
      END DO                                                                                        
      timesum = 0.0                                                                                 
      nwqcnt = 0                                                                                    
c                                                                                                   
      RETURN                                                                                        
      END                                                                                           
