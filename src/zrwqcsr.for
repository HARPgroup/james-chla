C                                                                                                   
C**********************************************************************C                            
C**********************************************************************C                            
C**********************************************************************C                            
C                                                                                                   
      SUBROUTINE RWQCSR(KC,DT,N,TCON,TBEGIN,NWQV_1,NTSPTC,
     &  NCSTEP,SECDLAST)                                               
C                                                                                                   
C**********************************************************************C                            
C                                                                                                   
C **  LAST MODIFIED BY JOHN HAMRICK ON  7 APRIL 1997        
C **  Last modified by AEE 3/21/2007                                        
C                                                                                                   
C**********************************************************************C                            
C                                                                                                   
C     READ AND UPDATE WATER COLUMN WATER QUALITY VARIABLE TIME SERIES                               
C                                                                                                   
C**********************************************************************C                            
C                                                                                                   
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn' 
      REAL*8 SECDLAST 
	dimension WKQ(KCWM)                                                    
      CHARACTER*11 FNWQSR(21)                                          
C                                                                       
C**********************************************************************C
C                                                                       
      IF(ITNWQ.GT.0) GO TO 1000                                         
C                                                                       
C**********************************************************************C
C                                                                       
       FNWQSR( 1)='cwqsr01.inp'                                         
       FNWQSR( 2)='cwqsr02.inp'                                         
       FNWQSR( 3)='cwqsr03.inp'                                         
       FNWQSR( 4)='cwqsr04.inp'                                         
       FNWQSR( 5)='cwqsr05.inp'                                         
       FNWQSR( 6)='cwqsr06.inp'                                         
       FNWQSR( 7)='cwqsr07.inp'                                         
       FNWQSR( 8)='cwqsr08.inp'                                         
       FNWQSR( 9)='cwqsr09.inp'                                         
       FNWQSR(10)='cwqsr10.inp'                                         
       FNWQSR(11)='cwqsr11.inp'                                         
       FNWQSR(12)='cwqsr12.inp'                                         
       FNWQSR(13)='cwqsr13.inp'                                         
       FNWQSR(14)='cwqsr14.inp'                                         
       FNWQSR(15)='cwqsr15.inp'                                         
       FNWQSR(16)='cwqsr16.inp'                                         
       FNWQSR(17)='cwqsr17.inp'                                         
       FNWQSR(18)='cwqsr18.inp'                                         
       FNWQSR(19)='cwqsr19.inp'                                         
       FNWQSR(20)='cwqsr20.inp'                                         
       FNWQSR(21)='cwqsr21.inp'                                         
C                                                                       
C**********************************************************************C
C                                                                       
C **  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE WQ CONCENTRATION       
C **  TIME SERIES FROM THE FILES cwqsrNN.inp                            
C                                                                       
C----------------------------------------------------------------------C
C                                                                       
      DO NW=1,NWQV_1                                                    
C                                                                       
      IF(NWQCSR(NW).GE.1) THEN                                          
        OPEN(1,FILE=FNWQSR(NW),STATUS='UNKNOWN')                        
C                                                                       
C **  SKIP OVER TITLE AND AND HEADER LINES                              
C                                                                       
      DO IS=1,15                                                        
      READ(1,1)                                                         
      END DO                                                            
C                                                                       
        DO NS=1,NWQCSR(NW)                                              
        MWQCTLT(NS,NW)=1                                                
        READ(1,*,IOSTAT=ISO)ISTYP,MWQCSR(NS,NW),TCWQCSR(NS,NW),         
     $                   TAWQCSR(NS,NW),RMULADJ,ADDADJ                  
        IF(ISO.GT.0) GO TO 900                                          
                                                        
        if(MWQCSR(NS,NW).gt.NDWQCSR) then   ! check array dimensions    
        write(6,*)" Error:  MWQCSR(NS,NW)>NDWQCSR", NS,NW,MWQCSR(NS,NW),
     +  NDWQCSR                                                         
        stop                                                            
        endif                                                           
                                                           
        IF(ISTYP.EQ.1) THEN                                             
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)                          
          IF(ISO.GT.0) GO TO 900                                        
           DO M=1,MWQCSR(NS,NW)                                         
           READ(1,*,IOSTAT=ISO)TWQCSER(M,NS,NW),CSERTMP                 
           IF(ISO.GT.0) GO TO 900                                       
           TWQCSER(M,NS,NW)=TWQCSER(M,NS,NW)+TAWQCSR(NS,NW)             
            DO K=1,KC                                                   
            WQCSER(M,K,NS,NW)=(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)         
            END DO                                                      
           END DO                                                       
         ELSE                                                           
          DO M=1,MWQCSR(NS,NW)                                          
          READ(1,*,IOSTAT=ISO)TWQCSER(M,NS,NW),                         
     $                       (WQCSER(M,K,NS,NW), K=1,KC)                
          IF(ISO.GT.0) GO TO 900                                        
          TWQCSER(M,NS,NW)=TWQCSER(M,NS,NW)+TAWQCSR(NS,NW)              
           DO K=1,KC                                                    
           WQCSER(M,K,NS,NW)=RMULADJ*(WQCSER(M,K,NS,NW)+ADDADJ)         
           END DO                                                       
          END DO                                                        
        END IF                                                          
        END DO                                                          
        CLOSE(1)                                                        
      END IF                                                            
C                                                                       
      END DO                                                            
                                                                        
      GO TO 901                                                         
C                                                                       
  900 CONTINUE                                                          
      WRITE(6,601)NW,NS,M                                               
      STOP                                                              
C                                                                       
  901 CONTINUE                                                          
C                                                                       
    1 FORMAT(120X)                                                      
  601 FORMAT(' READ ERROR WQ TIME SERIES, NWQ,NSER,MDATA = ',3I5)       
  602 FORMAT(' READ OF FILES cwqsrNN.inp SUCCESSFUL'/)                  
C                                                                       
C**********************************************************************C
C                                                                       
 1000 CONTINUE                                                          
C                                                                       
C**********************************************************************C
C                                                                       
C **  INITIALIZE NULL SERIES CONCENTRATIONS                             
C                                                                       
      DO NW=1,NWQV_1                                                      
       DO K=1,KC                                                        
        CSERTWQ(K,0,NW)=0.                                              
       END DO                                                           
      END DO                                                            
C                                                                       
C**********************************************************************C
C                                                                       
C **  CONCENTRATION SERIES INTERPOLTATION FOR WATER QUALITY VARIABLES   
C                                                                       
      DO NW=1,NWQV_1                                                      
C                                                                       
       DO NS=1,NWQCSR(NW)                                               
       TIME=DT*FLOAT(N-1)/TCWQCSR(NS,NW)                                
     &           +TBEGIN*(TCON/TCWQCSR(NS,NW))                          
       IF(NCSTEP.GT.0) TIME=(SECDLAST-DT)/TCON+TBEGIN   
C                                                                       
       M1=MWQCTLT(NS,NW)                                                
  100  CONTINUE                                                         
       M2=M1+1                                                          
       IF (TIME.GT.TWQCSER(M2,NS,NW)) THEN                              
         M1=M2                                                          
         GO TO 100                                                      
        ELSE                                                            
         MWQCTLT(NS,NW)=M1                                              
       END IF                                                           
C                                                                       
       TDIFF=TWQCSER(M2,NS,NW)-TWQCSER(M1,NS,NW)                        
       WTM1=(TWQCSER(M2,NS,NW)-TIME)/TDIFF                              
       WTM2=(TIME-TWQCSER(M1,NS,NW))/TDIFF                              
        DO K=1,KC                                                       
         CSERTWQ(K,NS,NW)=WTM1*WQCSER(M1,K,NS,NW)                       
     &                   +WTM2*WQCSER(M2,K,NS,NW)                       
        END DO                                                          
                                          
       END DO                                                           
C                                                                       
      END DO                                                            
C                                                                       
C**********************************************************************C
C                                                                       
C **  ON FIRST CALL, INITIALIZE OUTFLOW CONCENTRATIONS                  
C                                                                       
C----------------------------------------------------------------------C
C                                                                       
      IF (ITNWQ.EQ.0) THEN                                              
C                                                                       
      DO NW=1,NWQV_1                                                      
      DO K=1,KC                                                         
       DO LL=1,NWQOBS                                                   
        NSID=IWQOBS(LL,NW)                                              
        L=LIJW( IWQCBS(LL),JWQCBS(LL) )                                  
!        CWQLOS(LL,K,NW)=WTCI(K,1)*WQOBCS(LL,1,NW)                       
!     $    +WTCI(K,2)*WQOBCS(LL,2,NW)+CSERTWQ(K,NSID,NW) 
        CWQLOS(LL,K,NW)=CSERTWQ(K,NSID,NW)                              
        NWQLOS(LL,K,NW)=0                                               
       END DO                                                           
       DO LL=1,NWQOBW                                                   
        NSID=IWQOBW(LL,NW)                                              
        L=LIJW( IWQCBW(LL),JWQCBW(LL) )                                  
!        CWQLOW(LL,K,NW)=WTCI(K,1)*WQOBCW(LL,1,NW)                       
!     $    +WTCI(K,2)*WQOBCW(LL,2,NW)+CSERTWQ(K,NSID,NW)  
        CWQLOW(LL,K,NW)=CSERTWQ(K,NSID,NW)                              
        NWQLOW(LL,K,NW)=0                                               
       END DO                                                           
       DO LL=1,NWQOBE                                                   
        NSID=IWQOBE(LL,NW)                                              
        L=LIJW( IWQCBE(LL),JWQCBE(LL) )                                  
        CWQLOE(LL,K,NW)=CSERTWQ(K,NSID,NW)
!        CWQLOE(LL,K,NW)=WTCI(K,1)*WQOBCE(LL,1,NW)                       
!     $    +WTCI(K,2)*WQOBCE(LL,2,NW)+CSERTWQ(K,NSID,NW)                       
        NWQLOE(LL,K,NW)=0                                               
       END DO                                                           
       DO LL=1,NWQOBN                                                   
        NSID=IWQOBN(LL,NW)                                              
        L=LIJW( IWQCBN(LL),JWQCBN(LL) )                                  
!        CWQLON(LL,K,NW)=WTCI(K,1)*WQOBCN(LL,1,NW)                       
!     $    +WTCI(K,2)*WQOBCN(LL,2,NW)+CSERTWQ(K,NSID,NW)
        CWQLON(LL,K,NW)=CSERTWQ(K,NSID,NW)                              
        NWQLON(LL,K,NW)=0                                               
       END DO                                                           
      END DO                                                            
      END DO                                                            
C                                                                       
      END IF                                                            
C                                                                       
C**********************************************************************C
C                                                                       
      RETURN                                                            
      END                                                               
