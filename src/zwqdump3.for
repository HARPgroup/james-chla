      SUBROUTINE WQDUMP3(N,KC,LA,DT,NTSPTC,TCON,TBEGIN,iyear1,
     & NCSTEP,SECDLAST,IWQPT,ITTTT)
!  
!  save EFDC WQ output 
!
! These are the index numbers of the parameters in WQVo array:
! WQVo(LL,K, 1) = cyanobacteria as C  WQVo(LL,K,12) = labile PON
! WQVo(LL,K, 2) = diatoms as C        WQVo(LL,K,13) = diss. org. nitrogen
! WQVo(LL,K, 3) = green algae as C    WQVo(LL,K,14) = ammonia nitrogen
! WQVo(LL,K, 4) = refractory POC      WQVo(LL,K,15) = nitrate nitrogen
! WQVo(LL,K, 5) = labile POC          WQVo(LL,K,16) = part. biogenic silica
! WQVo(LL,K, 6) = diss. org. carbon   WQVo(LL,K,17) = available silica
! WQVo(LL,K, 7) = refractory POP      WQVo(LL,K,18) = COD
! WQVo(LL,K, 8) = labile POP          WQVo(LL,K,19) = dissolved oxygen
! WQVo(LL,K, 9) = diss. org. phos.    WQVo(LL,K,20) = total active metal
! WQVo(LL,K,10) = tot. inorg. phos.   WQVo(LL,K,21) = fecal coliform bacteria
! WQVo(LL,K,11) = refractory PON      WQVo(LL,K,22) = macroalgae
!
! WQTOTA(L ,K,1)= Total Algae
! WQTOTA(L ,K,2)= Total N 
! WQTOTA(L ,K,3)= Total P
! WQTOTA(L ,K,4)= PO4
! WQTOTA(L ,K,5)= Si
!-----------------------------------------------------------------------------------------
!
      INCLUDE 'wq.par' 
      INCLUDE 'wqcom.cmn'
      INTEGER, SAVE :: inibin=1,iavg=0,iyears=0
      character*4 iyear
       REAL*8 SECDLAST  
 !     write(*,*)'wq write year= ',iyear1,iyears
 !       ITRYCHLA=1  ! Use Cercl vary C:Chla ration HAB 4/9/2015
      IF(iyears.ne.iyear1) then
        write(iyear,'(I4)')iyear1
        write(*,*)'WQ SAVE YEAR ',iyear
        iyears=iyear1
        inibin=1
      ENDIF
!     write(*,*)'Year == ',iyears,iyear,iyear1

      if(inibin.EQ.1) THEN
!	 IAVGBIN=NTSPTC/24 * 24
!       IF (IWQICI.LE.2) THEN
!	 write(*,*)'Continue to write to the exising wq3d.bin is turned off '
!       IWWQS =0
!       ELSEIF(IWQICI.EQ.3) THEN
!       IWWQS=1
!       write(*,*)'Continue to write to the exising wq3d.bin '
!       ENDIF
       
!	 if(IWWQS.LE.2.and.isBIN.GE.1) then

	 if(isBIN.GE.1) then
  
        open(1,file='wq3d'//iyear//'.bin',form='unformatted') 
        CLOSE(1,STATUS='DELETE')
        open(1,file='wq3d'//iyear//'.bin',form='unformatted')
	  write(1)NWQVM,KCWM,LCMWQ,ICMWQ,JCMWQ
	  write(1)LIJW,ILW,JLW
	  close(1)
	  open(1,file='wqsed3d'//iyear//'.bin',form='unformatted')	
	  CLOSE(1,STATUS='DELETE')
 	  
        IF(isWQMIN.GE.1) then
         open(1,file='wq3dmax'//iyear//'.bin',form='unformatted') 
         CLOSE(1,STATUS='DELETE')
         open(1,file='wq3dmax'//iyear//'.bin',form='unformatted')
	   write(1)NWQVM,KCWM,LCMWQ,ICMWQ,JCMWQ
	   write(1)LIJW,ILW,JLW
	   close(1)
         open(1,file='wq3dmin'//iyear//'.bin',form='unformatted') 
         CLOSE(1,STATUS='DELETE')
         open(1,file='wq3dmin'//iyear//'.bin',form='unformatted')
	   write(1)NWQVM,KCWM,LCMWQ,ICMWQ,JCMWQ
	   write(1)LIJW,ILW,JLW
	   close(1)  	 
        ENDIF
        
       endif
!       if(IWWQS.LE.2) then
        OPEN(1,FILE='production'//iyear//'.out',STATUS='UNKNOWN') 
        CLOSE(1,STATUS='DELETE') 
        OPEN(1,FILE='production'//iyear//'.out',STATUS='UNKNOWN')        
         write(1,*)'I  J  Time  H  I0  Alg   Ke  Pro Pro_e' 
        close(1)
!       endif
!
! Arry dimension 
!
! LIJW(ICMWQ,JCMWQ ),ILW(LCMWQ),JLW(LCMWQ)
! Time,WQ3DA(LCMWQ,KCWM,NWQVM),WQTOT(LCMWQ,KCWM,5)
!

       DO NW=1,NWQVM
        DO K=1,KCWM
         DO L=1,LCMWQ
          WQ3DA(L,K,NW)=0.0
         ENDDO
        ENDDO
       ENDDO
   
       DO NW=1,NWQVM
        DO K=1,KCWM
         DO L=1,LCMWQ
          WQVmin(L,K,NW)=1.0E6
          WQVmax(L,K,NW)=0.0      
         ENDDO
        ENDDO
       ENDDO
        DO I=1,5
         DO K=1,KC
          DO L=1,LA
          WQTOTmax(L,K,I)=0.0
          WQTOTmin(L,K,I)=1.0E6
          WQTOT(L,K,I)=0.0          
          ENDDO
         ENDDO
        ENDDO 
!
       DO K=1,KCWM
       DO L=2,LA       
       WQketot(L,K)=0 
       ENDDO
       ENDDO
       
       DO L=2,LA
       xAlg(L)=0.0
       xPro(L)=0.0
        xKe(L)=0.0
       xI00(L)  = 0.0
       xHXY(L)  = 0.0        
       ENDDO  
!      
      inibin=0
!
	 IF(N.LE.2) RETURN
	ENDIF      ! end initical for each yar

!     Algae: convert to chl a

C       DO K=1,KC
C        DO L=2,LA
C         WQ3DA(L ,K,1)=WQ3DA(L ,K,1)+WQVO(L,K,1)*WQCHLC(IWQZMAP(L,KC))   
C         WQ3DA(L ,K,2)=WQ3DA(L ,K,2)+WQVO(L,K,2)*WQCHLD(IWQZMAP(L,KC))
C         WQ3DA(L ,K,3)=WQ3DA(L ,K,3)+WQVO(L,K,3)*WQCHLG(IWQZMAP(L,KC))
C        ENDDO
C       ENDDO
       
       IF(ITRYCHLA.EQ.3) THEN
      	KECST=30           ! also active respiration with growth
	 ELSEIF(ITRYCHLA.EQ.2)THEN
	  KECST=40
	 ENDIF
       
       
        IF(ITRYCHLA.GE.2) THEN             !JS  2018
        DO LL=2,LA
        WQCHLD(IWQZMAP(LL,KC))=1/(KECST+90*exp(-1.19*WQketot(LL,KC)))
     &  *1000      
        WQCHLC(IWQZMAP(LL,KC))=1/(KECST+150*exp(-1.18*WQketot(LL,KC)))
     &  *1000  
        WQCHLG(IWQZMAP(LL,KC))=WQCHLC(IWQZMAP(LL,KC)) 
        ENDDO
        ENDIF
       
 !      IF(ITRYCHLA.GE.2) THEN
 !       DO L=2,LA
 !        CCHa_1=CK_b+90*exp(-1.19*WQketot(L,KC))/1000      
 !        CCHa_2=CK_b+150*exp(-1.18*WQketot(L,KC))/1000 
 !       DO K=1,KC          
 !        WQ3DA(L ,K,1)=WQ3DA(L ,K,1)+WQVO(L,K,1)*CCHa_1   
 !        WQ3DA(L ,K,2)=WQ3DA(L ,K,2)+WQVO(L,K,2)*CCHa_2
 !        WQ3DA(L ,K,3)=WQ3DA(L ,K,3)+WQVO(L,K,3)*CCHa_1
 !       ENDDO
 !      ENDDO
       
c      ELSE
      
!        DO K=1,KC
!         DO L=2,LA
!         WQ3DA(L ,K,1)=WQ3DA(L ,K,1)+WQVO(L,K,1)*WQCHLC(IWQZMAP(L,KC))   
!         WQ3DA(L ,K,2)=WQ3DA(L ,K,2)+WQVO(L,K,2)*WQCHLD(IWQZMAP(L,KC))
!         WQ3DA(L ,K,3)=WQ3DA(L ,K,3)+WQVO(L,K,3)*WQCHLG(IWQZMAP(L,KC))
!         ENDDO
!        ENDDO
     
c       ENDIF
      
c       DO K=1,KC
c        DO L=2,LA
c         WQ3DA(L ,K,1)=WQ3DA(L ,K,1)+WQVO(L,K,1)   
c         WQ3DA(L ,K,2)=WQ3DA(L ,K,2)+WQVO(L,K,2)
c         WQ3DA(L ,K,3)=WQ3DA(L ,K,3)+WQVO(L,K,3)
c        ENDDO
c       ENDDO
       
!   other variables

       DO NW=1,NWQV
        DO K=1,KC
         DO L=2,LA
          WQ3DA(L ,K,NW)=WQ3DA(L ,K,NW)+WQVO(L,K,NW)
         ENDDO
        ENDDO
       ENDDO
      
! Total C

       DO K=1,KC
        DO L=2,LA
	   DO I=1,6
         WQTOT(L ,K,1)=WQTOT(L ,K,1)+WQVO(L,K,i)
         ENDDO
        ENDDO
       ENDDO
!
!  Total N
!
       DO K=1,KC
        DO L=2,LA
	   DO I=11,15
         WQTOT(L ,K,2)=WQTOT(L ,K,2)+WQVO(L,K,I)
         ENDDO
        ENDDO
       ENDDO
       DO K=1,KC
        DO L=2,LA
         WQTOT(L,K,2)=WQTOT(L,K,2)+WQVO(L,K,1)*WQANCC1(IMWQZT(L))+
     &   WQVO(L,K,2)*WQANCD1(IMWQZT(L))+WQVO(L,K,3)*WQANCG1(IMWQZT(L))
        ENDDO
       ENDDO
!
!  Total P
!
       DO K=1,KC
        DO L=2,LA
	   DO I=7,10
         WQTOT(L ,K,3)=WQTOT(L ,K,3)+WQVO(L,K,I)
         ENDDO
        ENDDO
       ENDDO

        DO K=1,KC
        DO L=2,LA
         APCWQ = 1.0 / (WQCP1PRM1(IMWQZT(L))
     &      + WQCP2PRM*EXP(-WQCP3PRM*XPO4DWQ))
         TPWQ1 = APCWQ*(WQVO(L,K,1)+WQVO(L,K,2)+WQVO(L,K,3))
         WQTOT(L ,K,3)=WQTOT(L ,K,3)+TPWQ1
        ENDDO
       ENDDO      
       
       DO K=1,KC
        DO L=2,LA
         WQTOT(L ,K,4)=WQTOT(L,K,4)+
     &                  WQVO(L,K,10) / (1.0 + WQKPO4P*SEDTWQ(L,K))   ! dissoved PO4
         WQTOT(L ,K,5)=WQTOT(L,K,5)+ 
     &  	              WQVO(L,K,17) / (1.0 + WQKSAP*SEDTWQ(L,K))  ! dissoved DSi
        ENDDO
       ENDDO
!
! Find min and max
!
!      write(*,*)'Check min/max IWQPT,IAVGBIN,iavg= ',
!     &       IWQPT,IAVGBIN,iavg
      if(isWQMIN==1) then
      
!      write(*,*)'Check min/max IWQPT,IAVGBIN,iavg= ',
!     &       IWQPT,IAVGBIN,iavg
       
       DO NW=1,NWQVM
       DO L=2,LA
       DO K=1,KC
         WQVmin(L,K,NW)=min(WQVmin(L,K,NW),WQVO(L,K,NW))
         WQVmax(L,K,NW)=max(WQVmax(L,K,NW),WQVO(L,K,NW))     
       ENDDO
       ENDDO
       ENDDO
! TC
       DO L=2,LA
        DO K=1,KC	
         Ctmp=0.0   
         DO I=1,6
          Ctmp=Ctmp+WQVO(L,K,I)
         ENDDO
          WQTOTmax(L,K,1)=max(WQTOTmax(L,K,1),Ctmp )
          WQTOTmin(L,K,1)=min(WQTOTmin(L,K,1),Ctmp)        
        ENDDO
       ENDDO
!  Total N
        DO L=2,LA
         DO K=1,KC	
           Ctmp=0.0   
 	    DO I=11,15 
 	     Ctmp=Ctmp+WQVO(L,K,I)   
          ENDDO 
          Ctmp=Ctmp+WQVO(L,K,1)*WQANCC1(IMWQZT(L))
     &        +WQVO(L,K,2)*WQANCD1(IMWQZT(L))
     &        +WQVO(L,K,3)*WQANCG1(IMWQZT(L))   	                     
          WQTOTmax(L,K,2)=max(WQTOTmax(L,K,2), Ctmp)   
          WQTOTmin(L,K,2)=min(WQTOTmin(L,K,2), Ctmp)               
        ENDDO
       ENDDO
!TP       
        DO K=1,KC
         DO L=2,LA
          Ctmp=0.0
          DO I=7,10
          Ctmp=Ctmp+WQVO(L,K,I) 
          ENDDO
          APCWQ = 1.0 / (WQCP1PRM1(IMWQZT(L))
     &      + WQCP2PRM*EXP(-WQCP3PRM*XPO4DWQ))      
          Ctmp=Ctmp+APCWQ*(WQVO(L,K,1)+WQVO(L,K,2)+WQVO(L,K,3))       
          WQTOTmin(L,K,3)=min(WQTOTmin(L,K,3), Ctmp)      ! TP  
          WQTOTmax(L,K,3)=max(WQTOTmax(L,K,3), Ctmp)                      
          Ctmp= WQVO(L,K,10)/(1.0 + WQKPO4P*SEDTWQ(L,K))  ! PO4
          WQTOTmin(L,K,4)=min(WQTOTmin(L,K,4), Ctmp)         
          WQTOTmax(L,K,4)=max(WQTOTmax(L,K,4), Ctmp)
          Ctmp= WQVO(L,K,17) / (1.0 + WQKSAP*SEDTWQ(L,K))  ! DSi                  
 !         WQTOTmin(L,K,5)=min(WQTOTmin(L,K,5), Ctmp)         
 !         WQTOTmax(L,K,5)=max(WQTOTmax(L,K,5), Ctmp)
         ENDDO
        ENDDO
               
       endif 
!
      iavg=iavg+1
      IF(MOD(IWQPT,IAVGBIN).EQ.0) THEN
!      write(*,*)'** Check save max N,IWQPT,IAVGBIN,iavg= '
!      write(*,*)'       ', N,IWQPT,IAVGBIN,iavg
C
C ---- Average results----------------
C
        DO NW=1,NWQV
         DO K=1,KC
          DO L=2,LA
           WQ3DA(L,K,NW)=WQ3DA(L,K,NW)/iavg
          ENDDO
         ENDDO
        ENDDO

        DO NW=1,5
         DO K=1,KC
          DO L=2,LA
           WQTOT(L,K,NW)=WQTOT(L,K,NW)/iavg      
          ENDDO
         ENDDO
        ENDDO
 
       
       DO K=1,KCWM
       DO L=2,LA       
       WQketot(L,K)=WQketot(L,K)/iavg/ITTTT  
       WQTOT(L,K,5)=WQketot(L,K) 
       ENDDO
       ENDDO
        
       write(iyear,'(I4)')iyears
    
       DO L=2,LA     
       xAlg(L)=xAlg(L)/iavg
       xKe(L)=xKe(L)/iavg
       xI00(L)=xI00(L)/10.45/PARADJ  ! check
 !      if(xI00(L).GT.60.0) xI00(L)=60.0
       xHXY(L)=xHXY(L)/iavg      
       ENDDO  
 101  format(2I6,3F9.2,F12.5,F10.3,8F12.3) 
 
      OPEN(1,FILE='production'//iyear//'.out',
     &  STATUS='UNKNOWN')
       TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/86400
       IF(NCSTEP.GT.0) TIMTMP=SECDLAST/TCON+TBEGIN 
      DO M=1,IWQTS                              
          LL=LWQTS(M)
          att=3.8*(xAlg(LL)*xI00(LL)/(xKe(LL)+1.0E-8))*0.001
          write(1,101)ILW(LL),JLW(LL),TIMTMP,xHXY(LL),xI00(LL),
     &    xAlg(LL),xKe(LL),xPro(LL),att,real(iavg)    
      ENDDO
      CLOSE(1)
!     
       DO K=1,KC
        DO L=2,LA
         WQ3DA(L ,K,1)=WQ3DA(L ,K,1)*WQCHLC(IWQZMAP(L,KC))   
         WQ3DA(L ,K,2)=WQ3DA(L ,K,2)*WQCHLD(IWQZMAP(L,KC))
         WQ3DA(L ,K,3)=WQ3DA(L ,K,3)*WQCHLG(IWQZMAP(L,KC))
        ENDDO
       ENDDO
 
      DO K=1,KC
       DO L=2,LA
        WQVmin(L,K,1)=WQVmin(L,K,1)*WQCHLC(IWQZMAP(L,KC))   ! convert to chl a
        WQVmin(L,K,2)=WQVmin(L,K,2)*WQCHLD(IWQZMAP(L,KC)) 
        WQVmin(L,K,3)=WQVmin(L,K,3)*WQCHLG(IWQZMAP(L,KC)) 
        WQVmax(L,K,1)=WQVmax(L,K,1)*WQCHLC(IWQZMAP(L,KC))  
        WQVmax(L,K,2)=WQVmax(L,K,2)*WQCHLD(IWQZMAP(L,KC))
        WQVmax(L,K,3)=WQVmax(L,K,3)*WQCHLG(IWQZMAP(L,KC))     
       ENDDO
      ENDDO   
                 
      if(isBIN.GE.1) then    
                       
	 open(1,file='wq3d'//iyear//'.bin',form='unformatted')	
        TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/86400
        IF(NCSTEP.GT.0) TIMTMP=SECDLAST/TCON+TBEGIN 
        WRITE(1) TIMTMP,WQ3DA,WQTOT  
 	  close(1)  
!        
! Save sediment flux here

	  open(1,file='wqsed3d'//iyear//'.bin',form='unformatted')	
        WRITE(1) TIMTMP,WQBFO2,WQBFNH4,WQBFNO3,WQBFPO4D  
 	  close(1)  
                           
       if(isWQMIN.GE.1) then
   
 	  open(1,file='wq3dmax'//iyear//'.bin',form='unformatted',
     &  position='append')	
        TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/86400
        IF(NCSTEP.GT.0) TIMTMP=SECDLAST/TCON+TBEGIN 
        WRITE(1) TIMTMP,WQVmax,WQTOTmax  
 	  close(1)
      
        open(1,file='wq3dmin'//iyear//'.bin',form='unformatted',
     &  position='append')	
        WRITE(1) TIMTMP,WQVmin,WQTOTmin  
	  close(1)
	  	  	  
       endif
      endif	  
   
      IWQPT=0
      CALL WWQTS(TINDAY,KC,DT,N,TCON,TBEGIN,NTSPTC,NCSTEP,SECDLAST)

		           	  
        DO NW=1,NWQV
         DO K=1,KC
          DO L=1,LA
          WQVmin(L,K,NW)=1.0E6
          WQVmax(L,K,NW)=0.0         
          ENDDO
         ENDDO
        ENDDO	 

        DO I=1,5
         DO K=1,KC
          DO L=1,LA
          WQTOTmax(L,K,I)=0.0
          WQTOTmin(L,K,I)=1.0E6
          ENDDO
         ENDDO
        ENDDO                                    
!
       DO NW=1,NWQV
        DO K=1,KC
         DO L=1,LA
          WQ3DA(L,K,NW)=0.0
         ENDDO
        ENDDO
       ENDDO

       DO NW=1,5
        DO K=1,KCWM
         DO L=1,LCMWQ
          WQTOT(L,K,NW)=0.0
         ENDDO
        ENDDO
       ENDDO
      
      !
       DO L=2,LA
       xAlg(L)=0.0
       xPro(L)=0.0
        xKe(L)=0.0
       xI00(L)  = 0.0
       xHXY(L)  = 0.0    
       ENDDO  

       DO K=1,KCWM
       DO L=2,LA       
       WQketot(L,K)=0 
       ENDDO
       ENDDO
           
       iavg=0

	ENDIF !IF(MOD(IWQPT,IAVGBIN).EQ.0)
C
 2000 CONTINUE
C
      RETURN
      END
