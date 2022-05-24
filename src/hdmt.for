C 
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE HDMT
C
C **  SUBROUTINE HDMT EXECUTES THE FULL HYDRODYNAMIC AND MASS TRANSPORT
C **  TIME INTERGATION
C **  MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C     Last modified by Jian Shen on July 6, 2016
C
C     if(irunpcb.EQ.1.and.ISTRAN(5).GE.1) THEN  ! Run PCB
C     if(ifed_inc.LT.0.and.ISTRAN(8).ge.1)
C     IF(NTRNVA.GT.0) run dye using saved field
C     if(idiaout.LT.0)  out out resuts
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn' 
      INCLUDE 'wq.par' 
      INCLUDE 'wqcom.cmn'     
      DIMENSION TOXAV(LCM),TOXAVS(LCM)  
      DIMENSION TOX_TMP(9),PTMDL(20),PTMDLS(20),MAPTM(LCM),NPTMDL(20)
      character*2 IDYEE
      
      INTEGER NTMDLAVG
      
      isb_grid=ISUBMODE  ! RUN Subgrid model C93
      write(*,*)'Submodel on/off =',isb_grid
      write(*,*)'Run residence time ',NTRNVA 
 !     CALL OMP_SET_NUM_THREADs(4)
            
      ITMDL=0
      IF(ITMDL.GT.0) THEN
       NTMDLAVG =0   
       DO L=1,LCM
        TOXAV(L)=0
        TOXAVS(L)=0
       ENDDO
       DO L=1,20
        PTMDL(L)=0.
        PTMDLS(L)=0.
        NPTMDL(L)=0
       ENDDO
       NGTMDL=0
       OPEN(1,file='TMDLgroup.inp')
       read(1,*)
       DO L=2,LA
        read(1,*)L0,MAPTM(L0)
        NPTMDL(MAPTM(L0))=NPTMDL(MAPTM(L0))+1
        IF(MAPTM(L0).gt.NGTMDL) NGTMDL=MAPTM(L0)
       ENDDO
 109   CLOSE(1)
      ENDIF
      
      IF(ISCRAY.EQ.0) THEN
        TTMP=SECNDS(0.0)
       ELSE
        TTMP=SECOND( )
        CALL TIMEF(WTTMP)
      END IF
!
      ITRWQOP=1
C---------------------------------------------------------------------c
C Dye study use saved resutls and wq transport
C---------------------------------------------------------------------c
      IF(NTRNVA.GT.0) THEN
      
      ITNWQ=0
      
      DO L=1,20
      NWQCSR(L)=1
      ENDDO
      
      NWQOBS=NCBS
      NWQOBW=NCBW
      NWQOBE=NCBE
      NWQOBN=NCBN
 !     write(*,*)'CK BC =',NWQOBE
      DO L=1,NCBS
      IWQCBS(L)=ICBS(L)
      JWQCBS(L)=JCBS(L)
      ENDDO
      
      DO L=1,NCBW
      IWQCBW(L)=ICBW(L)
      JWQCBW(L)=JCBW(L)
      ENDDO
     
      DO L=1,NCBE
      IWQCBE(L)=ICBE(L)
      JWQCBE(L)=JCBE(L)
      ENDDO  
               
      DO L=1,NCBN
      IWQCBN(L)=ICBN(L)
      JWQCBN(L)=JCBN(L)
      ENDDO   
 !-------------------------------------------------------------------
 !    Run dye
 !------------------------------------------------------------------- 
      IF(IRELSOP.EQ.2) THEN   ! Run dye
        
        OPEN(501,FILE='rdy1cnh.bin',STATUS='UNKNOWN',
     &   form='binary')
        CLOSE(501,STATUS='DELETE')
        OPEN(501,FILE='rdy1cnh.bin',STATUS='UNKNOWN',
     &   form='binary')
        WRITE (501) LC-1,KC+2
  
        OPEN(502,FILE='rdy2cnh.bin',STATUS='UNKNOWN',
     &   form='binary')
        CLOSE(502,STATUS='DELETE')
        OPEN(502,FILE='rdy2cnh.bin',STATUS='UNKNOWN',
     &   form='binary')
        WRITE (502) LC-1,KC+2
        
        OPEN(503,FILE='dyets1.bin',STATUS='UNKNOWN',
     &   form='binary')
         CLOSE(503,STATUS='DELETE')    
         OPEN(503,FILE='dyets1.bin',STATUS='UNKNOWN',
     &   form='binary')
          WRITE (503) KC+2
 
         OPEN(504,FILE='agets1.bin',STATUS='UNKNOWN',
     &   form='binary')
         CLOSE(504,STATUS='DELETE')    
         OPEN(504,FILE='agets1.bin',STATUS='UNKNOWN',
     &   form='binary')
         WRITE (504) KC+2
         
       IF(NTRNVA.GT.2) THEN     

        OPEN(601,FILE='rdy3cnh.bin',STATUS='UNKNOWN',
     &   form='binary')
        CLOSE(601,STATUS='DELETE')
        OPEN(601,FILE='rdy3cnh.bin',STATUS='UNKNOWN',
     &   form='binary')
        WRITE (601) LC-1,KC+2
  
        OPEN(602,FILE='rdy4cnh.bin',STATUS='UNKNOWN',
     &   form='binary')
        CLOSE(602,STATUS='DELETE')
        OPEN(602,FILE='rdy4cnh.bin',STATUS='UNKNOWN',
     &   form='binary')
        WRITE (602) LC-1,KC+2
        
        OPEN(603,FILE='dyets2.bin',STATUS='UNKNOWN',
     &   form='binary')
         CLOSE(603,STATUS='DELETE')    
         OPEN(603,FILE='dyets2.bin',STATUS='UNKNOWN',
     &   form='binary')
          WRITE (603) KC+2
 
         OPEN(604,FILE='agets2.bin',STATUS='UNKNOWN',
     &   form='binary')
         CLOSE(604,STATUS='DELETE')    
         OPEN(604,FILE='agets2.bin',STATUS='UNKNOWN',
     &   form='binary')
         WRITE (604) KC+2
         
         ENDIF 
         
         ELSE        !IRELSOP=3
         
         OPEN(501,FILE='rdy1cnh.bin',STATUS='UNKNOWN',
     &   form='binary')
         CLOSE(501,STATUS='DELETE')
         OPEN(501,FILE='rdy1cnh.bin',STATUS='UNKNOWN',
     &   form='binary')
         WRITE (501) LC-1,KC+2       
         
         ENDIF
                               
        DO NW=1,NTRNVA
        DO K=1,KC
         DO L=1,LA+1
          WQV(L,K,NW)=0
          WQVO(L,K,NW)=0
         END DO
        END DO 
        ENDDO
        
        NAVGCT=0
        write(*,*)'DT DT2= ',DT, DT2
C-----------------------------------------------------
C   Inverse computing residence time
C----------------------------------------------------       
      IF(IRELSOP.EQ.5) THEN
        write(*,*)'Computing residence time using inverse model'
   
      IF(NTRNVA.GT.2) THEN   ! Read locat residence time set up
        write(*,*)'Computing residence time using inverse model' 
       
        DO NW=1,NTRNVA 
         write(IDYEE,'(I2)')NW
         IDDD=700+NW
         OPEN(IDDD,FILE='residenc'//IDYEE//'.bin',STATUS='UNKNOWN',
     &   form='binary')
         CLOSE(IDDD,STATUS='DELETE')
         OPEN(IDDD,FILE='residenc'//IDYEE//'.bin',STATUS='UNKNOWN',
     &   form='binary')
         WRITE (IDDD) LC-1,KC+2       
        ENDDO    
             
        WQ3DA=0
        open(1,file='localRse.txt')
        
        DO i=1,99999
         read(1,*,end=919)La,Lb,Lc,Ld
         DO K=1,KC
          L0=LIJ(Lb,Lc)
          WQ3DA(L0,K,Ld+1)=1
         ENDDO
        ENDDO
 919    close(1)
        DO L=1,LA
        DO J=1,KC
        WQ3DA(L,J,1)=1
        ENDDO
        ENDDO
       ENDIF
        
  !      DO N=1,NTS
         DO N=1,(TES-TBS)*NTSPTC 
  !      TIME=DT*FLOAT(N)+TCON*TBEGIN
         TIME=TCON*TES-DT*FLOAT(N)
         TIME=TIME/TCON    
 !        write(*,*) TIME,TCON,TES
!        IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014             
        
        if(mod(N,2).EQ.0) then
         if(isb_grid.eq.1) then
	    !CALL flow_in_sub
	    CALL flow_in_inv_sub
	   else
 	   IF(IVRDOP.EQ.0)THEN
 	    CALL flow_in
 	   ELSE
 	    CALL flow_in_inv  !use direct access
 	   ENDIF
 	   endif
	   CALL WQMAPPING
         CALL ZDYE_inv               
        
         ISHPRT=NTSMMT
         if(mod(N,ISHPRT).EQ.0) then
          Open(1,file='out.txt')
	    DO L=2,LA
	    if(NTRNVA.LE.2) then
	     write(1,118)L
     &    ,WQV(L,KC,1)
     &    ,WQV(L, 1,1)
     &    ,U(L,KC),HWQ(L),TIME 
          else
   	     write(1,'(I6,4E12.3,F8.2)')L
     &    ,WQV(L,KC,1)
     &    ,WQV(L,KC,2),WQV(L,KC,3),WQV(L,KC,4),TIME        
          endif
	    ENDDO 
	   endif 

         IF(mod(N,NTSMMT).EQ.0)then   
         DO NW=1,NTRNVA
          IDDD=700+NW  
          write(IDDD)TIME   
          DO L=2,LA      
           write(IDDD)IL(L),JL(L),(WQVO(L,K,NW)/NAVGCT, K=1,KC) 
 !          write(*,*)'IL JL= ', IL(L),JL(L),LA  
          ENDDO 
         ENDDO      

         DO NW=1,NTRNVA
          DO K=1,KC
          DO L=2,LA 
           WQVO(L,K,NW)=0
          ENDDO
          ENDDO
         ENDDO
          NAVGCT=0    
C         
         ELSE
C
         DO NW=1,NTRNVA
          DO K=1,KC
          DO L=2,LA 
           WQVO(L,K,NW)=WQVO(L,K,NW)+WQV(L,K,NW)
          ENDDO
          ENDDO
        ENDDO
           NAVGCT=NAVGCT+1              
C
         END IF      

       endif  ! mode(N,2)       	           
      ENDDO  ! Ned N loop
      
      ENDIF
C            
C------End inverse computing residence time-------------
C-----------------------------------------------------
C   Inverse computing residence at selected locations 
C----------------------------------------------------       
      IF(IRELSOP.EQ.6) THEN

        write(*,*)'Computing residence time at selected location' 
         OPEN(500,FILE='residtime.bin',STATUS='UNKNOWN',
     &   form='binary')
         CLOSE(500,STATUS='DELETE')
         OPEN(500,FILE='residtime.bin',STATUS='UNKNOWN',
     &   form='binary')            
         WQ3DA=0
        
         IF(MLTMSR>20)MLTMSR=20
         NTRNVA=MLTMSR
         DO MS=1,MLTMSR 
          I0=ILTMSR(MS)
          J0=JLTMSR(MS)
          LL=LIJ(I0,J0)
          DO K=1,KC
          WQ3DA(LIJ(I0,J0),K,MS)=1
          WQ3DA(LIJ(I0+1,J0),K,MS)=1
          WQ3DA(LIJ(I0-1,J0),K,MS)=1
          WQ3DA(LIJ(I0,J0+1),K,MS)=1
          WQ3DA(LIJ(I0,J0-1),K,MS)=1
          ENDDO        
         ENDDO 
   
  !      DO N=1,NTS
         DO N=1,(TES-TBS)*NTSPTC 
  !      TIME=DT*FLOAT(N)+TCON*TBEGIN
         TIME=ICON*TES-DT*FLOAT(N)
         TIME=TIME/TCON    

!        IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014             
        
        if(mod(N,2).EQ.0) then
         if(isb_grid.eq.1) then
	    !CALL flow_in_sub
	    CALL flow_in_inv_sub
	   else
 	   IF(IVRDOP.EQ.0)THEN
 	    CALL flow_in
 	   ELSE
 	    CALL flow_in_inv  !use direct access
 	   ENDIF
 	   endif
	   CALL WQMAPPING
         CALL ZDYE_inv               
        
         ISHPRT=NTSMMT
         if(mod(N,ISHPRT).EQ.0) then
          Open(1,file='out.txt')
	    DO L=2,LA
	    if(NTRNVA.LE.2) then
	     write(1,118)L
     &    ,WQV(L,KC,1)
     &    ,WQV(L, 1,1)
     &    ,U(L,KC),HWQ(L),TIME 
          else
   	     write(1,'(I5,4E12.3,F7.2)')L
     &    ,WQV(L,KC,1)
     &    ,WQV(L,KC,2),WQV(L,KC,3),WQV(L,KC,4),TIME        
          endif
	    ENDDO 
	   endif 

         IF(mod(N,NTSMMT).EQ.0)then   

          write(500)TIME  
          DO NW=1,NTRNVA
           I0=ILTMSR(NW)
           J0=JLTMSR(NW)
           LL=LIJ(I0,J0)  
           write(500)I0,J0,(WQVO(LL,K,NW)/NAVGCT, k=1,kc)    
          ENDDO    

          DO NW=1,NTRNVA
           DO K=1,KC
            DO L=2,LA 
            WQVO(L,K,NW)=0
            ENDDO
           ENDDO
          ENDDO
          NAVGCT=0    
C         
         ELSE
C
          DO NW=1,NTRNVA
           DO K=1,KC
            DO L=2,LA 
             WQVO(L,K,NW)=WQVO(L,K,NW)+WQV(L,K,NW)
            ENDDO
           ENDDO
          ENDDO
          NAVGCT=NAVGCT+1              
C
         END IF      
c
       endif  ! mode(N,2)       	           
      ENDDO  ! Ned N loop
      
      ENDIF
C===================================================================
C Compute age
C===================================================================

      IF(IRELSOP.EQ.1) THEN
      
       DO NW=1,NTRNVA 
       if(NW<10) then
         write(IDYEE,'(I1)')NW
       else
         write(IDYEE,'(I2)')NW       
       endif
         IDDD=500+NW
         OPEN(IDDD,FILE='Age'//IDYEE//'.bin',STATUS='UNKNOWN',
     &   form='binary')
         CLOSE(IDDD,STATUS='DELETE')
         OPEN(IDDD,FILE='Age'//IDYEE//'.bin',STATUS='UNKNOWN',
     &   form='binary')
         WRITE (IDDD) LC-1,KC+2       
       ENDDO 
        WQ3DA=0
        DO L=1,LA
        DO J=1,KC
        WQ3DA(L,J,1)=1  ! default total age
        ENDDO
        ENDDO  
             
!      IF(NTRNVA.GT.5) THEN   ! Read locat residence time set up
!        write(*,*)'Computing local age'      
!        open(1,file='localRse.txt')       
!        DO i=1,99999
!         read(1,*,end=918)La,Lb,Lc,Ld
!         DO K=1,KC
!          L0=LIJ(Lb,Lc)
!          WQ3DA(L0,K,Ld*2+1)=1
!         ENDDO
!        ENDDO
! 918    close(1)
!       ENDIF
          
      DO N=1,NTS

       TIME=DT*FLOAT(N)+TCON*TBEGIN
       TIME=TIME/TCON    
       IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 

 !       if(NCSTEP.GT.0)then
 !        SECDLAST=SECDLAST+DT        ! % time in second
 !        DAYLAST=SECDLAST/86400.0    ! % time in days  J.S. 12/24/2010
 !        TIME=secdlast/tcon+TBEGIN   ! % Time in day
 !       endif 

        iyr_tm=365
        if(mod(iyear1,4).eq.0)iyr_tm=366
        if (N.eq.int((tplus+iyr_tm-TBEGIN)*NTSPTC) )  then  ! at day 365
         iyear1=iyear1+1
         write(iyear,'(I4)')iyear1
         write(6,*)'YEAR ',iyear
	   tplus=tplus+365
	   if(mod(iyear1-1,4).eq.0)tplus=tplus+1
!	   call newinput  !(iyear, tplus)
        endif
 	     
 !	  CALL CALTSXY
        IF(IRELSOP.EQ.3) CALL CALQVS (3)     

       if(mod(N,2).EQ.0) then

	   if(isb_grid.eq.1) then
	    CALL flow_in_sub
	   else
	    CALL flow_in
	   endif
	   CALL WQMAPPING
         CALL ZDYE
C
        ISHPRT=NTSMMT
        if(mod(N,ISHPRT).EQ.0) then
         Open(1,file='out.txt',err=1009)
!         DO L=2,LA
!  	   write(1,118)L,WQV(L,KC,1),WQV(L,1,1),
!     &    WQV(L,KC,2)/(1.0e-6+WQV(L,KC,1))
!     &   ,WQV(L, 1,2)/(1.0e-6+WQV(L,1,1 )),U(L,KC),TIME 
!	   ENDDO
	   DO L=2,LA
	    IF(NTRNVA.LT.5) THEN 
	     	   write(1,118)L
     &   ,WQV(L, 1,2)/(1.0e-6+WQV(L,1,1))
     &   ,WQV(L, KC,2)/(1.0e-6+WQV(L,KC,1))
     &   ,WQV(L, 1,1)
     &   ,WQV(L, 1,2),U(L,KC),TIME 
          ELSE
 	     	   write(1,118)L
     &   ,WQV(L, 1,2)/(1.0e-6+WQV(L,1,1))
     &   ,WQV(L, 1,6)/(1.0e-6+WQV(L,1,5))
     &   ,WQV(L, 1,8)/(1.0e-6+WQV(L,1,7))
     &   ,WQV(L, 1,10)/(1.0e-6+WQV(L,1,9)),U(L,KC),TIME        
          ENDIF
 ! 	   write(1,118)L
 !    &   ,WQV(L,KC/2,4)/(1.0e-6+WQV(L,KC/2,3))
 !    &   ,WQV(L, 1,4)/(1.0e-6+WQV(L,1,3 ))
 !    &   ,WQV(L,KC,2)/(1.0e-6+WQV(L,KC,1))
 !    &   ,WQV(L, 1,2)/(1.0e-6+WQV(L,1,1 )),U(L,KC),TIME 
	   ENDDO
	   
1009    close(1)          
       endif
118   format(I5,2E14.5, 4F12.3)	      
      
      IF(mod(N,NTSMMT).EQ.0)then
      
       DO NW=1,NTRNVA
        DO K=1,KC
         DO L=2,LA 
          WQVO(L,K,NW)=WQVO(L,K,NW)/NAVGCT
         ENDDO
        ENDDO
       ENDDO     
        
       DO NW=1,NTRNVA/2 
        IDDD=500+NW*2-1  
        IDDD1=500+NW*2    
        write(IDDD)TIME 
       
       DO L=2,LA      
 !        write(IDDD)IL(L),JL(L), 
 !    &        ( WQV(L,K,NW*2)/(WQV(L,K,NW*2-1)+0.0001),k=1,kc)  
         write(IDDD)IL(L),JL(L), 
     &        ( WQVO(L,K,NW*2)/(WQVO(L,K,NW*2-1)+0.0001),k=1,kc)      
       ENDDO
      ENDDO
            
       DO NW=1,NTRNVA
        DO K=1,KC
        DO L=2,LA 
         WQVO(L,K,NW)=0
        ENDDO
        ENDDO
       ENDDO
       NAVGCT=0
       
      ELSE

       DO NW=1,NTRNVA
       DO K=1,KC
       DO L=2,LA 
        WQVO(L,K,NW)=WQVO(L,K,NW)+WQV(L,K,NW)
       ENDDO
       ENDDO
       ENDDO
       NAVGCT=NAVGCT+1
       
      END IF      
 
      endif  !END mod(N,2)
              
      ENDDO  ! N loop
C
      DO i=501,504
      close(i)
      enddo
      IF(NTRNVA.GT.2) THEN     
      DO i=601,604
      close(i)
      enddo
      ENDIF
      
      ENDIF  ! IRELSOP NE. 5 
      
      RETURN
      ENDIF   ! End dye study 
C
C---------------------------------------------------------------------C
C Read saved hydronamics and ouput timeseries
C---------------------------------------------------------------------C
      if(idiaout.LT.0) then
       NCSTEP=0
       NCTMSR=1  
       JSTMSR=1 
       write(*,*)'Read saved resutls and output timeseries for check!'
       DO N=1,NTS
 
        iyr_tm=365
        if (N.eq.int((tplus+iyr_tm-TBEGIN)*NTSPTC) )  then  ! at day 365
         iyear1=iyear1+1
         write(iyear,'(I4)')iyear1
         write(6,*)'YEAR ',iyear
	   tplus=tplus+365
	   if(mod(iyear1-1,4).eq.0)tplus=tplus+1
	  NCSTEP=0
        NCTMSR=1  
        JSTMSR=1 
	  endif
	  if(isb_grid.eq.1) then
	  CALL flow_in_sub
	  else
	  CALL flow_in
	  endif
C
C **  WRITE TO TIME SERIES FILES
C
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON    
        IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
C
        IF (ISTMSR.GE.1) THEN
         IF (N.GE.NBTMSR) THEN
          IF (NCTMSR.EQ.NWTMSR) THEN
             CALL TMSR
             ! CALL TMSR_asc
            ICALLTP=1
            NCTMSR=1  
           ELSE 
            NCTMSR=NCTMSR+1
          END IF
         END IF
        END IF
       ENDDO
      
      IF(ISCRAY.EQ.0) THEN
        THDMT=THDMT+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        THDMT=T2TMP-TTMP
        WTHDMT=(WT2TMP-WTTMP)*0.001
      END IF

      WRITE(6,6995)THDMT/60.0,TCONG
      return
      
      endif 
C -------------!end(idiaout.LT.0)
C
C----------------------------------------------------------------------------------c
C Run water quality model from saved dynamic field
C
C----------------------------------------------------------------------------------c
      ITRWQOP=0               
C             
      if(ifed_inc.LT.0.and.ISTRAN(8).ge.1) then  ! J.S. Set switch for transport WQ
      open(610,file='cflmax.out')
      
       NCTMSR=1  
       JSTMSR=1 
       ISHPRT=1440
 !      IWQDT0=6  !Read in from 'wq3dwc.inp'
       DO N=1,NTS
       
        TIME=DT*FLOAT(N)/86400+TBEGIN 
        SECDLAST=SECDLAST+DT         ! % time in second
        DAYLAST=SECDLAST/86400.0    ! % time in days  J.S. 12/24/2010       
        if(NCSTEP.GT.0)then
         TIME=SECDLAST/TCON+TBEGIN    ! % Time in day
        endif 
        
  !     call SHOWVAL2
                                 ! Temporaly use NQCTLT as controal
        if(NQCTLT>0) THEN        ! This is for Linganore and ConowingoLake only   JS 10/13/2018
         CALL CALQVS (3)      ! Get out flow and compute outflow loading
        endif
                    
        iyr_tm=365
         if(mod(iyear1,4).eq.0)iyr_tm=366
 !       if (N.eq.int((tplus+iyr_tm-TBEGIN)*NTSPTC) )  then  ! at day 365
         AT=DAYLAST-(tplus+real(iyr_tm)-TBEGIN)
!         write(*,*)'AT ==',TIME, AT,DAYLAST,SECDLAST,tplus,
!     &    (tplus+real(iyr_tm)-TBEGIN)
         if(abs(AT).LT.1.0e-6) then
         iyear1=iyear1+1
         write(iyear,'(I4)')iyear1
         write(6,*)'YEAR ',iyear
	   tplus=tplus+iyr_tm
	!   if(mod(iyear1-1,4).eq.0)tplus=tplus+1
	   call newinput  !(iyear, tplus)
        endif
 	     
	  CALL CALTSXY  ! get wind and rain
 !      CALL CALQVS (ISTL)     
        if(mod(N,2).EQ.0) then
	  if(isb_grid.eq.1) then
	  CALL flow_in_sub
	  else
	  CALL flow_in
	  endif
         CALL WQ3D(3,N,DT,TCON,TBEGIN,TIDALP,NTSPTC,
     $   IWQDT,LA,KC,IC,JC,IWQS,IOP_SAVE,iyear1,NCSTEP,SECDLAST)
    
         if(IOP_SAVE.GT.0) Call zwqfield_bin  ! save WQ resutls for PCB C# 92
 
	   endif
	                
      if(mod(N,ISHPRT).EQ.0) then
       OPEN(1,FILE='show.inp',STATUS='UNKNOWN')
         DO NSKIP=1,6
         READ(1,*)
         END DO
       READ(1,*)NSHTYPE,NSHOWR,ICSHOW,JCSHOW,ISHPRT
       READ(1,*)ZSSMIN,ZSSMAX,SSALMAX
       READ(1,*)IVV1,IVV2,ILY
       CLOSE(1) 
            
      Open(1,file='out.txt',err=1008)
	DO L=2,LA
!	write(1,11)L,DLON(L),DLAT(L),HP(L)+BELV(L),SAL(L,1),SAL(L,KC)
	Ag_1=WQVO(L,KC,1)/0.04
	Ag_2=WQVO(L,KC,2)/0.04
	Ag_3=WQVO(L,KC,3)/0.04
	Ag_4=WQVO(L,1,22)/0.1
	write(1,119)L,Ag_1,Ag_2,Ag_3,Ag_4,WQVO(L,ILY,IVV1),WQVO(L,ILY,IVV2)
     $,TIME 
	enddo
1008  close(1)

	endif
119   format(I5,8E12.4)	

C
C **  WRITE TO TIME SERIES FILES
C
       if(idiaout.eq.2) then
        TIME=DT*FLOAT(N)+TCON*TBEGIN
        TIME=TIME/TCON    
        IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
C
        IF (ISTMSR.GE.1) THEN
         IF (N.GE.NBTMSR) THEN
          IF (NCTMSR.EQ.NWTMSR) THEN
            CALL TMSR
            ICALLTP=1
            NCTMSR=1  
           ELSE 
            NCTMSR=NCTMSR+1
          END IF
         END IF
        END IF
       endif
C test change timestep
c
c  % Change timestep here J.S. 12/24/2010
c     
c   NCDAY(NLAST)=last timestep, it is set in efdc.for
C
      IF(NCSTEP.GT.0) THEN
       IF(N.EQ.NCDAY(NLAST))THEN
        NTSPTC=NTSPTCC(NLAST)           ! change to new timestep
        DT=TIDALP*FLOAT(NFLTMT)/FLOAT(NTSPTCC(NLAST))
        DTI=1./DT
        DT2=2.*DT
        NSTDAY=N                        ! record the last N before changing timestep
        DTD = DT/86400.0                ! wq step in day   
        IWQDT=IWQDT0*NTSPTC/NTSPTCC(1)  ! rescale IWQDT timestep
        if(mod(IWQDT,2).NE.0)IWQDT=IWQDT+1 
        NWQKDPT=IWQDT /2 
!        write(*,*) 'IWQDT NWQKDPT ',IWQDT, NWQKDPT
! Boundary
        DO LL=1,NWQOBW
        NTSCRW(LL)=NTSCRW(LL)*NTSPTC/NTSPTCC(NLAST-1)  
        ENDDO
        DO LL=1,NWQOBS
        NTSCRS(LL)=NTSCRS(LL)*NTSPTC/NTSPTCC(NLAST-1)  
        ENDDO   
        DO LL=1,NWQOBE
        NTSCRE(LL)=NTSCRE(LL)*NTSPTC/NTSPTCC(NLAST-1)  
        ENDDO
        DO LL=1,NWQOBN
        NTSCRN(LL)=NTSCRN(LL)*NTSPTC/NTSPTCC(NLAST-1) 
        ENDDO
        nrefx=NTSPTC*nref
        ISDRY=ISDRYC(NLAST)  !Reset wet-dry
        NLAST=NLAST+1
        Iaverage=NTSPTC*aver3d/24.0+0.5
        write(841,*) "2  ",i,Iaverage
        write(*,*)'Change timestep ',NSTDAY,Iaverage,nrefx,NLAST-1
        write(*,*)IAVGBIN ,NTSPTC,DTWQ,SMHSED(1)
!
! Average results CALMMT   ! Fix 2 chose i.e hourly or daily
! Note to ensure to get correct average results, only change timestep by day
!
        if(NTSMMT.LT.NTSPTCC(NLAST-1))THEN
         NTSMMT=NTSPTC/24
        else
         NTSMMT=NTSPTC       
        endif

! Reset TMSR
        Itimes=NTSPTC/24  ! % force output in every hour *savets/24.0+0.5
        NWTMSR=ITIMES
        ICALLTP=1
        NCTMSR=1
! Reset DUMP   ! not used when using efdcwin.inp
!       NSDUMP 
        ICALLTP=1
        NCDUMP=1 
! Reset DUMP2 and DUMP3
!   DUMP2 & DUMP3 use new control (MOD((N-NSTDAY),ITIMES)),
!   NSTDAY is the last N before changing timestep, ITIMES is updated   
c
! WQ   
        IF(ISTRAN(8).GE.1) THEN
         DTWQ = DTD*REAL(IWQDT)  !( IWQDT1) ! this will   
         if(mod(NTSPTC,IWQDT).NE.0) then
         write(*,*)'Check timestep  mod(NTSPTC,IWQDT).NE.0'
         stop
         endif
         DTWQO2 = DTWQ*0.5
         IWQTSDT=Itimes  
         WQTSDT=24.0*ITIMES/NTSPTC
         ISMTSDT=NTSPTC             ! sediment force to output daily
         IAVGBIN=NTSPTC/24*IAVGBIN0 ! output interval

         write(*,*)'WQ Change timestep IAVGBIN,ISMTSDT IWQDT NWQKDPT ='
         write(*,*) IAVGBIN,ISMTSDT,IWQDT, NWQKDPT
         ITNWQ=NTSPTC
         
          IF(IWQBEN.GT.0)then 
          ISMTDMBS1=ISMTDMBS
          ISMTCMBS1=ISMTCMBS
          IWQTSDT = NINT(WQTSDT*3600.0/DT) 
       !   ISMTSDT = NINT(savets*3600.0/DT)
          ISMTDMBS = NINT(SMTDMBS/DTWQ)
          ISMTCMBS = NINT(SMTCMBS/DTWQ)
          SM1OKMDP = 1.0/SMKMDP
          SMBST1 = 1.0 / (1.0 + SMKBST*DTWQ)
          R1=ISMTDMBS/ISMTDMBS1
          R2=ISMTCMBS/ISMTCMBS1
          DO LL=2,LA
           IF(ISMHYPD(L).GT.0)ISMHYPD(L)= ISMHYPD(L)*R1
          ENDDO
          write(*,*)'ISMH ',ISMTDMBS1,ISMTDMBS,ISMTCMBS1,ISMTCMBS,
     &      ISMHYPD(20)
        
          DO I=1,ISMZ
           SMDTOH(I) = DTWQ/SMHSED(I)
           SMHODT(I) = SMHSED(I)/DTWQ
           SM1DIFT(I) = SMDIFT * SMDTOH(I)/(SMHSED(I)+ 1.E-18)
           SM2DIFT(I) = 1.0 / (1.0 + SM1DIFT(I))
           SMW2DTOH(I) = 1.0 + SMW2(I)*SMDTOH(I)
           SMW2PHODT(I) = SMW2(I) + SMHODT(I)
          ENDDO
         ENDIF
        ENDIF   ! end if WQ parameter
c        
       ENDIF
      ENDIF
C---------------------------------end change timestep---------------------------- 
      ENDDO  ! End major loop

      IF(ISCRAY.EQ.0) THEN
        THDMT=THDMT+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        THDMT=T2TMP-TTMP
        WTHDMT=(WT2TMP-WTTMP)*0.001
      END IF

      WRITE(6,6995)THDMT/60.0,TCONG

      close(550)
      RETURN
	endif
C--------------------------------------------------------------------------------c
C                                                                                c 
C    Run PCB model from saved dynamic and water quality data                     c
C                                                                                c
C--------------------------------------------------------------------------------c
      if(irunpcb.EQ.1.and.ISTRAN(5).GE.1) THEN  ! Run PCB
      ! NTOX=NPCB
      write(*,*)'Run toxic model using saved transport fields !'
       ISTRAN(1)=0
       ISTRAN(2)=0
       ISTRAN(3)=0
       ISTRAN(4)=0    
      
      I_BUDGET=1 
      
      OPEN(89,FILE='budget.out',STATUS='UNKNOWN')
      CLOSE(89,STATUS='DELETE')
      OPEN(89,FILE='budget.out',STATUS='UNKNOWN')
      WRITE(89,*)'Budget'
      CLOSE(89)
      OPEN(89,FILE='budgetime.out',STATUS='UNKNOWN')
      CLOSE(89,STATUS='DELETE')
      OPEN(89,FILE='budgetime.out',STATUS='UNKNOWN') 
       write(89,*)'Total PCB during out period'
       write(89,*)'> 0 flux into system, <0 out from system'
       write(89,*)'Unit: ug/m^3/#day Current set 30 day'
       write(89,*)'Daily loading x 0.000001/30'       
       write(89,*)'End in w & sed. 2 3 SSEDEND,BSEDEND'
       write(89,*)'From boundary     4  -SEDOUT1'
       write(89,*)'From Point source 5  SEDIN1' 
       write(89,*)'From air          6  AIRFUX'
       write(89,*)'From sediment     7  SDFLUX+SMASSIN*DT2'
       write(89,*)'Bural             8  SEDERR*DT2'
       write(89,*)'Diff from sedimen 9  SMASSIN*DT2'
       write(89,*)'Settle to sed.    10 SDFLUX'
       write(89,*)'Total to w & sed. 11  S_NETW S_NETS'
       write(89,*)'Err in w & sed.   12 13 ERRO_W,ERRO_E'
       WRITE(89,*)'Re Err in w & S % 14 15   ERRO_WR,ERRO_SR'
         
      write(89,*)'1 TIME'
      write(89,*)'2 SSEDBEG'
      write(89,*)'3 BSEDBEG'
      write(89,*)'4 SSEDEND'
      write(89,*)'5 BSEDEND'
      write(89,*)'6 -SEDOUT1'
      write(89,*)'7 AIRFUX'
      write(89,*)'8 SDFLUX+SMASSIN*DT2 '
      write(89,*)'9 SEDERR*DT2'
      write(89,*)'10 SMASSIN*DT2'
      write(89,*)'11 SDFLUX'
      write(89,*)'12 S_NETW'
      write(89,*)'13 S_NETS'
      write(89,*)'14 ERRO_W'
      write(89,*)'15 ERRO_E'
      write(89,*)'16 ERRO_WR'
      write(89,*)'17 ERRO_SR'
      write(89,*)'18 TOXINBC'
      write(89,*)'19 TOXINBC'
      write(89,*)'20 ADEP'  
      write(89,*)'21 TOXBEDTOP ug'
      write(89,*)'22 NP/Poin souces ug'
      write(89,'(I10,22I12)')(I,I=1,22)          
      close(89)

        
      NCTMSR=1  
      JSTMSR=1            
         CALL WQ3D(3,N,DT,TCON,TBEGIN,TIDALP,NTSPTC,
     $   IWQDT,LA,KC,IC,JC,IWQS,IOP_SAVE,iyear1,NCSTEP,SECDLAST)
     
!      Open(1,file='outini.txt',err=977)
!	 write(1,*)'L,X,Y,Z,BB,SS'
!	DO L=2,LA
!	write(1,11)L,DLON(L),DLAT(L),HP(L)+BELV(L),SAL(L,1),SAL(L,KC)
!	write(1,11)L,DLON(L),DLAT(L),TOXB(L,2,1),TOX(L,1,1),TOX(L,KC,1)
!	enddo
!	close(1)

 !     CALL BUDGET1   ! Ini. budget
      
      ACV=86400.0
      ACV1=1/2500.0
!
! Restart
      if(RESPCB.ge.1) then
       if(RESPCB.eq.1) then
        open(1,file='toxrestart.inp')
       else
        open(1,file='toxrestart.out')
       endif
      DO NT=1,NPCB
      DO L=2,LA
      read(1,*)LL,(TOX(L,K,NT),K=1,KC),
     &  (TOXB(L,K,NT),K=1,2),(PORBEDP(L,K),K=1,2)
      ENDDO
      ENDDO
      DO K=1,KC
      DO NT=1,NPCB
      DO L=2,LA
      TOX1(L,K,NT)=TOX(L,K,NT)
      ENDDO
      ENDDO
      ENDDO
      DO K=1,2
      DO NT=1,NPCB
      DO L=2,LA 
      TOXB1(L,K,NT)=TOXB(L,K,NT)
      ENDDO
      ENDDO
      ENDDO   
      CLOSE(1)
      endif

      call SHOWVAL2
      
      DO N=1,NTS
      
        call SHOWVAL2
        
        TIME=DT*FLOAT(N)/86400+TBEGIN 
        SECDLAST=SECDLAST+DT         ! % time in second
        DAYLAST=SECDLAST/86400.0    ! % time in days  J.S. 12/24/2010       
        if(NCSTEP.GT.0)then
         TIME=SECDLAST/TCON+TBEGIN    ! % Time in day
        endif 

        iyr_tm=365
        AT=DAYLAST-(tplus+real(iyr_tm)-TBEGIN)
        if(abs(AT).LT.1.0e-6) then
         iyear1=iyear1+1
         write(iyear,'(I4)')iyear1
         write(6,*)'YEAR ',iyear
	   tplus=tplus+365
	   if(mod(iyear1-1,4).eq.0)tplus=tplus+1
	   call newinput  !(iyear, tplus)
        endif
 	     
	  CALL CALTSXY    ! get wind and rain
 !      CALL CALQVS (ISTL)  
    
      IF(mod(N,2).EQ.0) THEN
	   if(isb_grid.eq.1) then
	    CALL flow_in_sub
	   else
	    CALL flow_in
	   endif  
	   
	 CALL zwqfield_bin        ! Change zwqfield to zwqfield_bin using sequency read   
	 CALL WQMAPPING

	 CALL CALCSER (ISTL)      ! Boundary. use toxser.inp   J.S. 9/20/2016 
       CALL CALQVS (ISTL)   
       CALL pcbpoint            ! PCB point source input. Location read from "PCBser.inp"  J.S. 7/6/2016
 !
       CALL CALPCB(3)                ! transport
       CALL BUDGET2     
       CALL BUDGET3
       
	 CALL SSEDTOX_PCB(3,CORDT)     ! PCB dynamics

       IF(MOD(N,1440*30).EQ.0)then
 !       write(*,*)' SEDIN ',SEDIN
        CALL BUDGET5   ! BUDGET5 sum loads in this period.
 !       L = LIJW(IPCBpt(99),JPCBpt(99))
 !       write(*,*)' WQWPSL ',WQWPSL(L,1,1) ,SEDIN 
       ENDIF

       if(idiaout.gt.0) then
        if(mod(N,480).EQ.0) then
         Open(1,file='outpcb.txt',err=977)
         Open(2,file='outsd.txt',err=977)
         Open(3,file='outcon.txt',err=977)
	   write(1,*)'L,X,Y,SEDPOC,TOXb,TOXs,TOXFc, TOXFa,SEDFP'
	   write(2,*)'L,X,Y,TOXB2,TOXB1,PRO2,PRO1,SEDFP1,SEDFP2'
	   write(3,*)'L,X,Y,BPCB1,BPCB2,HB,PCBb,PCBs,PFb,PFs'
	   DO L=2,LA
	    IF(ITXBDUT(1).EQ.1) THEN
           T_TMP=TOXB(L,KBP,2)/HBEDP(L,KBP)/
     &     (2500.0*(1-PORBEDP(L,KBP)))
           T_TMP1=TOXB(L,KBP-1,2)/HBEDP(L,KBP-1)
     &      /(2500.0*(1-PORBEDP(L,KBP-1)))
          ENDIF
!	write(1,11)L,DLON(L),DLAT(L),HP(L)+BELV(L),SAL(L,1),SAL(L,KC)
	   write(1,111)L,DLON(L),DLAT(L),SEDPOC(L,2),TOX(L,1,1),TOX(L,KC,1),
     &	 TOXF(L,0,1)*ACV,TOXF(L,0,2)*ACV,(SEDFP(L,0,1)+SEDFP(L,0,2))*ACV
         write(2,112)L,DLON(L),DLAT(L),T_TMP,T_TMP1,
     &   PORBEDP(L,2),PORBEDP(L,1),SEDFP(L,0,2)*ACV,SEDFP(L,0,2)*ACV

	  FA_1=TOXPFTB(L,2,1)/HBEDP(L,2)/2500.0/(1-PORBEDP(L,2))
!	 FA_2=(1-TOXPFTB(L,2,1))/HBEDP(L,2)
	   FA_2=1.0/HBEDP(L,2)/2500.0/(1-PORBEDP(L,2))
	   write(3,113)L,DLON(L),DLAT(L),TOXB(L,2,1)*FA_1,            ! PPCB: ug/kg  DPCB:ug/m^3
     &   TOXB(L,2,1)*FA_2, HBEDP(L,2), TOX(L,1,1),
     &   TOX(L,KC,1), TOXPFTW(L,1,1),TOXPFTW(L,KC,1)                            ! PCB : ug/m^3
         ENDDO
	   close(1)
	   close(2)
	   close(3)
977	   continue
	  endif     
	 endif
c	
	ENDif  ! end (mod(N,2).EQ.0)
c
111   format(I8,',',F12.3,',',F12.3,',',5(E12.4,','),E12.4)
112   format(I8,',',F12.3,',',F12.3,',',5(E12.4,','),E12.4)
113   format(I8,',',F12.3,',',F12.3,',',6(E12.4,','),E12.4)
C**********************************************************************C
C
C **  WRITE TO TIME SERIES FILES
C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON    
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
C
 !     if(idiaout.gt.0)then
      IF (ISTMSR.GE.1) THEN
        IF (N.GE.NBTMSR) THEN
          IF (NCTMSR.EQ.NWTMSR) THEN
            CALL TMSR
      !      CALL TMSR_asc
            ICALLTP=1
            NCTMSR=1  
           ELSE
            NCTMSR=NCTMSR+1
          END IF
        END IF
      END IF
 !     endif
C
 !     N00=365*2*720
      N00=10
! Get Mean value
      IF(N.GT.N00) THEN
       NTMDLAVG=NTMDLAVG+1
       DO L=2,LA
        DO K=1,KC
        TOX_TMP(K)=0
        ENDDO
        TOXB_TMP=0.
        TOXA_TMP=0.
        DO NT=1,NPCB
        DO K=1,KC
        TOX_TMP(K)=TOX_TMP(K)+TOX(L,K,NT)
        TOXA_TMP=TOXA_TMP+TOX(L,K,NT)
        ENDDO
        TOXB_TMP=TOXB_TMP+TOXB(L,KBP,NT)*TOXPFTB(L,KBP,NT)/HBEDP(L,KBP)/
     &    (2500.0*(1-PORBEDP(L,KBP)))   ! model unit ug/m^2 convert tot ug/kg  (ng/g)       
        ENDDO
        
 !        PTMDL(MAPTM(L))= PTMDL(MAPTM(L))+ TOXA_TMP/KC*1.476
 !       PTMDLS(MAPTM(L))=PTMDLS(MAPTM(L))+TOXB_TMP*1.56
 !       TOXAV(L)=TOXAV(L)+TOXA_TMP/KC*1.476
 !       TOXAVS(L)=TOXAVS(L)+TOXB_TMP*1.56
       ENDDO
       
      ENDIF
      
!
C**********************************************************************C
	ENDDO
      IF(ISCRAY.EQ.0) THEN
        THDMT=THDMT+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        THDMT=T2TMP-TTMP
        WTHDMT=(WT2TMP-WTTMP)*0.001
      END IF

      WRITE(6,6995)THDMT/60.0,TCONG
      
      open(1,file='toxrestart.out')
      DO NT=1,NPCB
      DO L=2,LA
      WRITE(1,100)L,(TOX(L,K,NT),K=1,KC),
     &  (TOXB(L,K,NT),K=1,2),(PORBEDP(L,K),K=1,2)
      ENDDO
      ENDDO
      CLOSE(1)
 100  format(I8,50E12.5)
 
       open(2,file='pcbttmdl.txt')
       write(2,*)'L,X,Y,P1,P2,P3,P4,P5,P6,P7,P8,P9,SED,tPCB'
  !     write(2,*)TIME
       DO L=2,LA
        DO K=1,KC
        TOX_TMP(K)=0
        ENDDO
        TOXB_TMP=0.
        TOXA_TMP=0.
        DO NT=1,NPCB
        DO K=1,KC
        TOX_TMP(K)=TOX_TMP(K)+TOX(L,K,NT)
        TOXA_TMP=TOXA_TMP+TOX(L,K,NT)
        ENDDO
        TOXB_TMP=TOXB_TMP+TOXB(L,KBP,NT)*TOXPFTB(L,KBP,NT)/HBEDP(L,KBP)/
     &    (2500.0*(1-PORBEDP(L,KBP)))   ! model unit ug/m^2 convert tot ug/kg  (ng/g)       
        ENDDO
        TOXAV(L)=TOXAV(L)/NTMDLAVG
        TOXAVS(L)=TOXAVS(L)/NTMDLAVG    
        write(2,202)L,DLON(L),DLAT(L),(TOX_TMP(K)*1.476,K=1,KC),
     &   TOXAVS(L),TOXAV(L)
 !    *  TOXB_TMP*1.56,TOXA_TMP/KC*1.476

       ENDDO
       close(2)
202   format(I8,',',F12.3,',',F12.3,',',10(F9.3,','),F9.3)

      open(1,file='TMDLsum.out')
      write(1,*)'NTMDLAVG ',NTMDLAVG
      DO K=1,NGTMDL
       PTMDL(K)=PTMDL(K)/NPTMDL(K)/NTMDLAVG
       PTMDLS(K)=PTMDLS(K)/NPTMDL(K)/NTMDLAVG
      write(1,*)K,PTMDL(K),PTMDLS(K)
      ENDDO
      CLOSE(1)
      
      return
	endif

      FOURDPI=4./PI
C
C**********************************************************************C
C
C **  INITIALIZE COURNT NUMBER DIAGNOSTICS
C
      DO K=1,KC
      DO L=2,LA
       CFLUUU(L,K)=0.
       CFLUUU(L,K)=0.
       CFLWWW(L,K)=0.
      END DO
      END DO
C
C**********************************************************************C
C
      ILOGC=0
C
C**********************************************************************C
C
C **  CALCULATE U AT V AND V AT U USING ENERGY CONSERVING WEIGHTING
C **  CALCULATE VELOCITY GRADIENTS
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      LNW=LNWC(L)
      LSE=LSEC(L)
      LSW=LSWC(L)
      H1C(L)=0.25*(H1P(L)+H1P(L-1)+H1P(LS)+H1P(LSW))
      HMC(L)=0.25*(HMP(L)+HMP(L-1)+HMP(LS)+HMP(LSW))
      UV(L)=0.25*(HP(LS)*(U(LSE,1)+U(LS,1))
     $         +HP(L)*(U(L+1,1)+U(L,1)))*HVI(L)
      U1V(L)=0.25*(H1P(LS)*(U1(LSE,1)+U1(LS,1))
     $          +H1P(L)*(U1(L+1,1)+U1(L,1)))*H1VI(L)
      VU(L)=0.25*(HP(L-1)*(V(LNW,1)+V(L-1,1))
     $         +HP(L)*(V(LN,1)+V(L,1)))*HUI(L)
      V1U(L)=0.25*(H1P(L-1)*(V1(LNW,1)+V1(L-1,1))
     $          +H1P(L)*(V1(LN,1)+V1(L,1)))*H1UI(L)
      END DO
C
C**********************************************************************C
C
C **  CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
C
      IF (ISWAVE.EQ.1) CALL WAVEBL
      IF (ISWAVE.EQ.2) CALL WAVESXY
      IF(ISWAN.eq.1.and.istran(6).ge.1) then
      N=1           ! the first time, N has no value yet
      CALL wavetau  !Ji@wave, used SWAN's output,3/6/2000
      endif
C
C**********************************************************************C
C
C **  FIRST CALL TO INITIALIZE BOTTOM STRESS COEFFICINETS
C
      CALL CALTBXY
C
C**********************************************************************C
C
C **  CALCULATE HORIZONTAL VISCOSITY AND DIFFUSIVE MOMENTUM FLUXES
C
      IF (ISHDMF.GE.1) CALL CALHDMF
C
C**********************************************************************C
C
C **  CALCULATE BOTTOM AND SURFACE STRESS AT TIME LEVEL (N-1) AND N
C
C----------------------------------------------------------------------C
C
      N=-1
      CALL CALTSXY
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TBX1(L)=(AVCON1*H1UI(L)+STBX(L)*SQRT(V1U(L)*V1U(L)
     $        +U1(L,1)*U1(L,1)))*U1(L,1)
       TBY1(L)=(AVCON1*H1VI(L)+STBY(L)*SQRT(U1V(L)*U1V(L)
     $        +V1(L,1)*V1(L,1)))*V1(L,1)
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
       END DO
      END DO
C
c      IF(ISVEG.GE.1.AND.KC.GT.1)
c      DO ND=1,NDM
c       LF=2+(ND-1)*LDM
c       LL=LF+LDM-1
c       DO L=LF,LL
c       TBX1(L)=TBX1(L)+0.5*FXVEG(L,1)*SQRT(V1U(L)*V1U(L)
c     $        +U1(L,1)*U1(L,1))*U1(L,1)
c       TBY1(L)=TBY1(L)+0.5*FYVEG(L,1)*SQRT(U1V(L)*U1V(L)
c     $        +V1(L,1)*V1(L,1))*V1(L,1)
c       END DO
c      END DO
c      END IF
C
C**********************************************************************C
C
C **  SECOND CALL TO INITIALIZE BOTTOM STRESS COEFFICINETS
C
      CALL CALTBXY
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE STRESSES
C
C----------------------------------------------------------------------C
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TBX(L)=(AVCON1*HUI(L)+STBX(L)*SQRT(VU(L)*VU(L)
     $       +U(L,1)*U(L,1)))*U(L,1)
       TBY(L)=(AVCON1*HVI(L)+STBY(L)*SQRT(UV(L)*UV(L)
     $       +V(L,1)*V(L,1)))*V(L,1)
       END DO
      END DO
C
c      IF(ISVEG.GE.1.AND.KC.GT.1)
c      DO ND=1,NDM
c       LF=2+(ND-1)*LDM
c       LL=LF+LDM-1
c       DO L=LF,LL
c       TBX(L)=TBX(L)+0.5*FXVEG(L,1)*SQRT(VU(L)*VU(L)
c     $        +U(L,1)*U(L,1))*U(L,1)
c       TBY(L)=TBY(L)+0.5*FYVEG(L,1)*SQRT(U1V(L)*UV(L)
c     $        +V(L,1)*V(L,1))*V(L,1)
c       END DO
c      END DO
c      END IF
C
      N=0
      CALL CALTSXY
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
C
C----------------------------------------------------------------------C
C
C     IF (KC.GT.1.OR.ISTRAN(4).GE.1) THEN
C
      IF (ISWAVE.EQ.0) THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY1(LNC(L))
       TVAR3W(L)=TSX1(L+1)
       TVAR3E(L)=TBX1(L+1   )
       TVAR3N(L)=TBY1(LNC(L))
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       QQ1(L,0 )=0.5*CTURB2*SQRT((TVAR3E(L)+TBX1(L))**2
     $                        +(TVAR3N(L)+TBY1(L))**2)
       QQ1(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX1(L))**2
     $                         +(TVAR3S(L)+TSY1(L))**2)
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       QQ(L,0 )=0.5*CTURB2*SQRT((TVAR3E(L)+TBX(L))**2
     $                        +(TVAR3N(L)+TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2
     $                         +(TVAR3S(L)+TSY(L))**2)
       END DO
      END DO
C
      END IF
C
C     END IF
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
C
C----------------------------------------------------------------------C
C
C     IF (KC.GT.1.OR.ISTRAN(4).GE.1) THEN
C
      IF (ISWAVE.GE.1) THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY1(LNC(L))
       TVAR3W(L)=TSX1(L+1)
       TVAR3E(L)=TBX1(L+1   )
       TVAR3N(L)=TBY1(LNC(L))
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
         TAUBC2=0.25*( (TVAR3E(L)+TBX1(L))**2
     $                        +(TVAR3N(L)+TBY1(L))**2 )
         TAUBC=SQRT(TAUBC2)
         UTMP=0.5*STCUV(L)*(U1(L+1,1)+U1(L,1))+1.E-12
         VTMP=0.5*STCUV(L)*(V1(LN,1)+V1(L,1))
         CURANG=ATAN2(VTMP,UTMP)
         TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))
     $        +FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))
         TAUB2=MAX(TAUB2,0.)
         QQ1(L,0 )=CTURB2*SQRT(TAUB2)
         QQ1(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX1(L))**2
     $                         +(TVAR3S(L)+TSY1(L))**2)
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
         TAUBC2=0.25*( (TVAR3E(L)+TBX(L))**2
     $                        +(TVAR3N(L)+TBY(L))**2 )
         TAUBC=SQRT(TAUBC2)
         UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
         VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))
         CURANG=ATAN2(VTMP,UTMP)
         TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))
     $        +FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))
         TAUB2=MAX(TAUB2,0.)
         QQ(L,0 )=CTURB2*SQRT(TAUB2)
         QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2
     $                         +(TVAR3S(L)+TSY(L))**2)
       END DO
      END DO
C
      END IF
C
C     END IF
C
C**********************************************************************C
C
C **   SET SWITCHES FOR THREE TIME LEVEL STEP
C
       ISTL=3
       DELT=DT2
       DELTD2=DT
       DZDDELT=DZ/DELT
       ROLD=0.
       RNEW=1.
C
C**********************************************************************C
C**********************************************************************C
C
C **  BEGIN TIME LOOP FOR FULL HYDRODYNAMIC AND MASS TRANSPORT
C **  CALCULATION
C
C **  SET CYCLE COUNTER AND CALL TIMER
C
      NTIMER=0
      N=0
C
C **  dtime and FLUSH are supported on SUN systems, but may not be
C **  supported on other systems.
C
      CALL TIMELOG(N)
c     CALL dtime (tarray)
      WRITE(9,200)N, tarray(1),tarray(2)
c     CALL FLUSH(9)
  200 FORMAT(1X,'N=',I10,5X,'USER TIME=',E12.4,5X,'SYSTEM TIME=',E12.4)
      NTIMER=1
C
C----------------------------------------------------------------------C
! Test toxic model transport using salinity. Set ini. condition
!      DO LL=1,NPCB
!       DO K=1,KC
!        DO L=1,LA
!        TOX(L,K,LL)=SAL(L,K)
!        TOX1(L,K,LL)=SAL(L,K)
!        ENDDO
!       ENDDO
!      ENDDO
C
      DO 1000 N=1,NTS   

      if(NCSTEP.GT.0)then
       SECDLAST=SECDLAST+DT        ! % time in second
       DAYLAST=SECDLAST/86400.0    ! % time in days  J.S. 12/24/2010
       TIME=SECDLAST/TCON+TBEGIN   ! % Time in day
      endif 
      
! Longterm option
! =========================================================================
       
	if(ilong.eq.1) then
	   iyr_tm=365
         if(mod(iyear1,4).eq.0)iyr_tm=366
! % change timestep       
      IF(NCSTEP.GT.0)then
          if ( TIME-(tplus+iyr_tm).GT.0.00001 ) then
            iyear1=iyear1+1        
		    write(iyear,'(I4)')iyear1
		    write(6,*)'YEAR ',iyear
	      tplus=tplus+365
	      if(mod(iyear1-1,4).eq.0)tplus=tplus+1
!	      if(NWSER.LT.100) NWSER=8
!	      if(iyear1.ge.2005) then
!		      if(NWSER.LT.100)  NWSER=10
!	      endif
	      call newinput  !(iyear, tplus)        
          endif
      ELSE    
         if (N.eq.int((tplus+iyr_tm-TBEGIN)*NTSPTC) )  then  ! at day 365
            iyear1=iyear1+1
		  write(iyear,'(I4)')iyear1
		  write(6,*)'YEAR ',iyear
!	      iyear='2004'
	      tplus=tplus+365
	      if(mod(iyear1-1,4).eq.0)tplus=tplus+1
!	      if(iyear1.ge.2005) then
!		      if(NWSER.LT.100)  NWSER=10
!		    endif
	      call newinput  !(iyear, tplus)
	   endif
	ENDIF  
	endif
!Longterm =========================================================================
C 


C
      NLOGTMP=2*NLOGTMP
      IF (ILOGC.EQ.NTSMMT) THEN
        CLOSE(8,STATUS='DELETE')
        OPEN(8,FILE='efdc.log',STATUS='UNKNOWN')
        IF (ISDRY.GT.0) THEN
          OPEN(1,FILE='drywet.log',STATUS='UNKNOWN')
          CLOSE(1,STATUS='DELETE')
        END IF
        IF (ISCFL.EQ.1) THEN
          OPEN(1,FILE='cfl.out',STATUS='UNKNOWN')
          CLOSE(1,STATUS='DELETE')
        END IF
        ILOGC=0
      END IF
      ILOGC=ILOGC+1
c	  WRITE(8,6626)ILOGC
 6626 FORMAT(' ILOGC = ',I10)
C
C
      IF (N.LE.NLTS) SNLT=0.
      IF (N.GT.NLTS.AND.N.LE.NTTS) THEN
         NTMP1=N-NLTS
         NTMP2=NTTS-NLTS+1
         SNLT=FLOAT(NTMP1)/FLOAT(NTMP2)
      END IF
      IF (N.GT.NTTS) SNLT=1.
C
      IF (N.LE.NTSVB) THEN
       GP=GPO*(FLOAT(N)/FLOAT(NTSVB))
      ELSE
       GP=GPO
      END IF
C
C----------------------------------------------------------------------C
C
C **  INITIALIZE VOLUME, MASS, MOMENTUM, AND ENERGY BALANCE
C
      IF (NCTBC.NE.NTSTBC.AND.ISBAL.GE.1) THEN
         CALL CALBAL1
         NTMP=MOD(N,2)
         IF(NTMP.EQ.0) THEN
           CALL CBALEV1
          ELSE
           CALL CBALOD1
         END IF
       END IF
c
c  ** initialize sediment budget calculation   (dlk 10/15)
c
      IF (NCTBC.NE.NTSTBC.AND.ISSBAL.GE.1) THEN
         CALL BUDGET1
C         NTMP=MOD(N,2)
C         IF(NTMP.EQ.0) THEN
C           CALL BUDGEV1
C          ELSE
C           CALL BUDGOD1
C         END IF
       END IF
C
C----------------------------------------------------------------------C
C
C **  REENTER HERE FOR TWO TIME LEVEL CORRECTION
C
  500 CONTINUE
C
C**********************************************************************C
C
C **  CALCULATE VERTICAL VISCOSITY AND DIFFUSIVITY AT TIME LEVEL (N)
C
      IF(ISCRAY.EQ.0) THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      END IF
      IF (KC.GT.1) THEN
         IF(N.EQ.1.OR.ISTOPT(0).EQ.0) CALL CALAVB (ISTL)
         IF(N.GT.1.AND.ISTOPT(0).GT.1) CALL CALAVB2 (ISTL)
      END IF
      IF(ISCRAY.EQ.0) THEN
        TAVB=TAVB+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TAVB=TAVB+T2TMP-T1TMP
        WTAVB=WTAVB+(WT2TMP-WT1TMP)*0.001
      END IF
C
C**********************************************************************C
C
C **  CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
C
      IF(ISTL.EQ.3) THEN
        IF (ISWAVE.EQ.1) CALL WAVEBL
        IF (ISWAVE.EQ.2) CALL WAVESXY
      IF(ISWAN.eq.1.and.istran(6).ge.1) CALL wavetau  !Ji@wave, used SWAN's output,3/6/2000
      END IF
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT MOMENTUM EQUATION TERMS
C
      IF(ISCRAY.EQ.0) THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      END IF
cx     IF (KC.EQ.1.AND.NDM.GE.2) THEN
cx        CALL CALEXP2 (ISTL) ! ORIGINAL OF CALEXP2 MOVED TO CALEXP1
cx       ELSE    ISCDMA.EQ.5 or 6 !CALEXP2 COPY CALEXP
        IF(ISCDMA.LE.2) CALL CALEXP (ISTL)
        IF(ISCDMA.EQ.3) CALL CALEXP3 (ISTL)
        IF(ISCDMA.EQ.4) CALL CALEXP3 (ISTL)
        IF(ISCDMA.EQ.5) CALL CALEXP2 (ISTL)
        IF(ISCDMA.EQ.6) CALL CALEXP2 (ISTL)
        IF(ISCDMA.EQ.9) CALL CALEXP9 (ISTL)
cx      END IF
      IF(ISCRAY.EQ.0) THEN
        TCEXP=TCEXP+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TCEXP=TCEXP+T2TMP-T1TMP
        WTCEXP=WTCEXP+(WT2TMP-WT1TMP)*0.001
      END IF
C
C**********************************************************************C
C
C **  UPDATE TIME VARIABLE VOLUME SOURCES AND SINKS, CONCENTRATIONS
C **  AND SURFACE ELEVATIONS
C
      CALL CALCSER (ISTL)
      CALL CALQVS (ISTL)
      PSERT(0)=0.
      IF (NPSER.GE.1) CALL CALPSER (ISTL)
C
C**********************************************************************C
C
C **  ADVANCE TIME VARIABLE SURFACE WIND STRESS AND LOAD INTO INTERNAL
C **  MODE FORCING (S&M SOLUTION ONLY)
C
C----------------------------------------------------------------------C
C
      IF (ISCDMA.GE.3.AND.ISCDMA.LE.8) THEN
      IF (ISTL.EQ.3) THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
       END DO
      END DO
C
      CALL CALTSXY
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
       DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
       END DO
      END DO
C
      DO L=2,LA
        FXE(L)=FXE(L)+DT*SUB(L)*DYU(L)*(TSX(L)-TSX1(L))
        FYE(L)=FYE(L)+DT*SVB(L)*DXV(L)*(TSY(L)-TSY1(L))
      END DO
C
      END IF
      END IF
C
C**********************************************************************C
C
C **  SOLVE EXTERNAL MODE EQUATIONS FOR P, UHDYE, AND VHDXE
C
      IF(ISCRAY.EQ.0) THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      END IF
C
      IF (ISDRY.EQ.-1) THEN
         CALL CALPUV(ISTL)
         GO TO 5555
      END IF
C
      IF (ISEVER.GE.1) THEN
        CALL CALPUVE(ISTL)
        GO TO 5555
      END IF
C
      IF(ISCDMA.LE.2.OR.ISCDMA.GE.9) THEN
        IF (ISDRY.EQ.0.AND.IRVEC.EQ.0) CALL CALPUV2(ISTL)
        IF (ISDRY.EQ.0.AND.IRVEC.EQ.2) CALL CALPUV5(ISTL)
        IF (ISDRY.EQ.0.AND.IRVEC.EQ.3) CALL CALPUV2(ISTL)
        IF (ISDRY.EQ.0.AND.IRVEC.EQ.4) CALL CALPUV5(ISTL)
        IF (ISDRY.EQ.0.AND.IRVEC.EQ.9) CALL CALPUV9(ISTL)
        IF (ISDRY.EQ.0.AND.IRVEC.EQ.99) CALL CALPUV9(ISTL)
        IF (ISDRY.EQ.0.AND.IRVEC.EQ.10) CALL CALPUVA(ISTL)
        IF (ISDRY.EQ.0.AND.IRVEC.EQ.11) CALL CALPUVB(ISTL)
        IF (ISDRY.EQ.0.AND.IRVEC.EQ.12) CALL CALPUVC(ISTL)
        IF (ISDRY.EQ.0.AND.IRVEC.EQ.13) CALL CALPUVD(ISTL)
        IF (ISDRY.EQ.0.AND.IRVEC.EQ.14) CALL CALPUVF(ISTL)
        IF (ISDRY.EQ.1.OR.ISDRY.EQ.2) CALL CALPUV5(ISTL)
        IF (ISDRY.EQ.11.OR.ISDRY.EQ.12) CALL CALPUV5(ISTL)
        IF (ISDRY.EQ.3.OR.ISDRY.EQ.4) CALL CALPUV5(ISTL)
        IF (ISDRY.EQ.99) CALL CALPUV5(ISTL)
        GO TO 5555
      END IF
C
      IF(ISCDMA.GT.2.AND.ISCDMA.LT.5) THEN
        IF (ISDRY.EQ.0.AND.IRVEC.EQ.0) CALL CALPUV3(ISTL)
        IF (ISDRY.EQ.0.AND.IRVEC.EQ.3) CALL CALPUV3(ISTL)
        GO TO 5555
      END IF
C
      IF (ISCDMA.GE.5.AND.ISCDMA.LE.8) THEN
        IF (ISDRY.EQ.1.OR.ISDRY.EQ.2) CALL CALPUV6(ISTL)
        IF (ISDRY.EQ.11.OR.ISDRY.EQ.12) CALL CALPUV6(ISTL)
        IF (ISDRY.EQ.3.OR.ISDRY.EQ.4) CALL CALPUV6(ISTL)
        IF (ISDRY.EQ.99) CALL CALPUV7(ISTL)
CX      IF (ISDRY.EQ.3.OR.ISDRY.EQ.4) CALL CALPUV6(ISTL) !7 MOVED TO 8
C       IF (ISDRY.EQ.3.OR.ISDRY.EQ.4) CALL CALPUV8(ISTL)
      END IF
C
 5555 CONTINUE
C
      IF(ISCRAY.EQ.0) THEN
        TPUV=TPUV+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TPUV=TPUV+T2TMP-T1TMP
        WTPUV=WTPUV+(WT2TMP-WT1TMP)*0.001
      END IF
C
C**********************************************************************C
C
C **  WRITE DIAGNOSTICS
C
C----------------------------------------------------------------------C
C
C **  dtime and FLUSH are supported on SUN systems, but may not be
C **  supported on other systems.
C
      IF (ISLOG.GE.1) THEN
      WRITE(8,17)N,ITER,RSQ,CFMAX,AVMAX,ABMIN,ABMAX,ABMIN
c     CALL FLUSH(8)
      END IF
C
C  17 FORMAT(1X,'N,ITER,AVMA,AVMI,ABMA,ABMI',2I5,4(1X,F8.4))
   17 FORMAT(1X,'N,ITER,RSQ,CFMAX,AVMA,AVMI,ABMA,ABMI',
     $        I7,I5,2E12.4,4(1X,F8.4))
C
      ERRMAX=MAX(ERRMAX,ERR)
      ERRMIN=MIN(ERRMIN,ERR)
      ITRMAX=MAX(ITRMAX,ITER)
      IRRMIN=MIN(ITRMIN,ITER)
C
C**********************************************************************C
C
C **  ADVANCE INTERNAL VARIABLES FOR THREE TIME LEVEL STEP
C
C----------------------------------------------------------------------C
C
      IF (ISTL.EQ.3) THEN
C
      DO K=1,KC
      DO L=2,LA
      UHDY2(L,K)=UHDY1(L,K)
      UHDY1(L,K)=UHDY(L,K)
      VHDX2(L,K)=VHDX1(L,K)
      VHDX1(L,K)=VHDX(L,K)
      U2(L,K)=U1(L,K)
      V2(L,K)=V1(L,K)
      U1(L,K)=U(L,K)
      V1(L,K)=V(L,K)
      W2(L,K)=W1(L,K)
      W1(L,K)=W(L,K)
      END DO
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  ADVANCE TIME VARIABLE SURFACE WIND STRESS AND LOAD INTO INTERNAL
C **  MODE FORCING
C
C----------------------------------------------------------------------C
C
      IF (ISCDMA.LE.2.OR.ISCDMA.GE.9) THEN
      IF (ISTL.EQ.3) THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TSX1(L)=TSX(L)
       TSY1(L)=TSY(L)
       END DO
      END DO
C
      CALL CALTSXY
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
       DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
       END DO
      END DO
C
      END IF
      END IF
C
C**********************************************************************C
C
C **  SOLVE INTERNAL SHEAR MODE EQUATIONS FOR U, UHDY, V, VHDX, AND W
C
C----------------------------------------------------------------------C
C
      IF(ISCRAY.EQ.0) THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      END IF
      IF (KC.GT.1) THEN
        CALL CALUVW (ISTL)
       ELSE
        DO ND=1,NDM
         LF=2+(ND-1)*LDM
         LL=LF+LDM-1
         DO L=LF,LL
          UHDY(L,1)=UHDYE(L)
          U(L,1)=UHDYE(L)*HUI(L)*DYIU(L)
          VHDX(L,1)=VHDXE(L)
          V(L,1)=VHDXE(L)*HVI(L)*DXIV(L)
          W(L,1)=0.
         END DO
        END DO
        CALL CALUVW (ISTL)
      END IF
      IF(ISCRAY.EQ.0) THEN
        TUVW=TUVW+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TUVW=TUVW+T2TMP-T1TMP
        WTUVW=WTUVW+(WT2TMP-WT1TMP)*0.001
      END IF
C
C**********************************************************************C
C
C **  CALCULATE SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATIONS
C **  AT TIME LEVEL (N+1)
C
C----------------------------------------------------------------------C
C
      CALL CALCONC (ISTL)
      
C
C----------------------------------------------------------------------C
C
      DO K=1,KB
      DO L=1,LC
       SEDBT(L,K)=0.
       SNDBT(L,K)=0.
      END DO
      END DO
C
      DO NS=1,NSED
       DO K=1,KB
       DO L=1,LC
        SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
       END DO
       END DO
      END DO
C
      DO NS=1,NSND
       DO K=1,KB
       DO L=1,LC
        SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)
       END DO
       END DO
      END DO
C
      DO K=1,KC
       DO L=1,LC
        SEDT(L,K)=0.
        SNDT(L,K)=0.
       END DO
      END DO
C
      DO NS=1,NSED
       DO K=1,KC
        DO L=1,LC
         SEDT(L,K)=SEDT(L,K)+SED(L,K,NS)
        END DO
       END DO
      END DO
C
      DO NS=1,NSND
       DO K=1,KC
        DO L=1,LC
         SNDT(L,K)=SNDT(L,K)+SND(L,K,NS)
        END DO
       END DO
      END DO
C
cjh5/13/97      DO NT=1,NTOX
cjh5/13/97       DO K=1,KC
cjh5/13/97        DO L=1,LC
cjh5/13/97         TOXPFT(L,K,NT)=0.
cjh5/13/97        END DO
cjh5/13/97       END DO
cjh5/13/97      END DO
C
cjh5/13/97      DO NT=1,NTOX
cjh5/13/97       DO NS=1,NSED+NSND
cjh5/13/97        DO K=1,KC
cjh5/13/97         DO L=1,LC
cjh5/13/97          TOXPFT(L,K,NT)=TOXPFT(L,K,NT)+TOXPF(L,K,NS,NT)
cjh5/13/97         END DO
cjh5/13/97        END DO
cjh5/13/97       END DO
cjh5/13/97      END DO
C
C----------------------------------------------------------------------C
C
C **  CHECK RANGE OF SALINITY AND DYE CONCENTRATION
C
      IF (ISMMC.EQ.1) THEN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF (SAL(L,K).GT.SALMAX) THEN
       SALMAX=SAL(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      END IF
      IF (SAL(L,K).LT.SALMIN) THEN
       SALMIN=SAL(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      END IF
      END DO
      END DO
C
      WRITE(6,6001)N
      WRITE(6,6002)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6003)SALMIN,IMIN,JMIN,KMIN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF (DYE(L,K).GT.SALMAX) THEN
       SALMAX=DYE(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      END IF
      IF (DYE(L,K).LT.SALMIN) THEN
       SALMIN=DYE(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      END IF
      END DO
      END DO
C
      WRITE(6,6004)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6005)SALMIN,IMIN,JMIN,KMIN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF (SFL(L,K).GT.SALMAX) THEN
       SALMAX=SFL(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      END IF
      IF (SFL(L,K).LT.SALMIN) THEN
       SALMIN=SFL(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      END IF
      END DO
      END DO
C
      WRITE(6,6006)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6007)SALMIN,IMIN,JMIN,KMIN
C
      END IF
C
C
      IF (ISMMC.EQ.2) THEN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF (TEM(L,K).GT.SALMAX) THEN
       SALMAX=TEM(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      END IF
      IF (TEM(L,K).LT.SALMIN) THEN
       SALMIN=TEM(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      END IF
      END DO
      END DO
C
      WRITE(6,6001)N
      WRITE(6,6008)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6009)SALMIN,IMIN,JMIN,KMIN
C
      END IF
C
 6001 FORMAT(1X,'N=',I10)
 6002 FORMAT(1X,'SALMAX=',F14.4,5X,'I,J,K=',(3I10))
 6003 FORMAT(1X,'SALMIN=',F14.4,5X,'I,J,K=',(3I10))
 6004 FORMAT(1X,'DYEMAX=',F14.4,5X,'I,J,K=',(3I10))
 6005 FORMAT(1X,'DYEMIN=',F14.4,5X,'I,J,K=',(3I10))
 6006 FORMAT(1X,'SFLMAX=',F14.4,5X,'I,J,K=',(3I10))
 6007 FORMAT(1X,'SFLMIN=',F14.4,5X,'I,J,K=',(3I10))
 6008 FORMAT(1X,'TEMMAX=',F14.4,5X,'I,J,K=',(3I10))
 6009 FORMAT(1X,'TEMMIN=',F14.4,5X,'I,J,K=',(3I10))
C
C**********************************************************************C
C
C **  CALCULATE SHELL FISH LARVAE AND/OR WATER QUALITY CONSTITUENT
C **  CONCENTRATIONS AT TIME LEVEL (N+1) AFTER SETTING DOULBE TIME
C **  STEP TRANSPORT FIELD
C
C----------------------------------------------------------------------C
C
      IF(ISTRAN(5).GE.1.OR.ISTRAN(8).GE.1.OR.ifed_inc.GT.0) THEN
      NTMP=MOD(N,2)
      IF (NTMP.EQ.0.AND.ISTL.EQ.3) THEN
C
C **  CALCULATE CONSERVATION OF VOLUME FOR THE WATER QUALITY ADVECTION
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        UHDY2E(L)=0.
        VHDX2E(L)=0.
        HWQ(L)=0.25*(H2P(L)+2.*H1P(L)+HP(L))
       END DO
      END DO
C
      DO K=1,KC
       DO L=2,LA
        UHDYWQ(L,K)=0.5*(UHDY1(L,K)+UHDY2(L,K))
        VHDXWQ(L,K)=0.5*(VHDX1(L,K)+VHDX2(L,K))
        UHDY2E(L)=UHDY2E(L)+UHDYWQ(L,K)*DZC(K)
        VHDX2E(L)=VHDX2E(L)+VHDXWQ(L,K)*DZC(K)
        UWQ(L,K)=0.5*(U1(L,K)+U2(L,K))
        VWQ(L,K)=0.5*(V1(L,K)+V2(L,K))
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3E(L)=UHDY2E(L+1)
       TVAR3N(L)=VHDX2E(LNC(L))
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
        TVAR2E(L,K)=UHDYWQ(L+1   ,K)
        TVAR2N(L,K)=VHDXWQ(LNC(L),K)
        END DO
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KS
        DO L=LF,LL
        WWQ(L,K)=SWB(L)*(WWQ(L,K-1)
     $       -DZC(K)*(TVAR2E(L,K)-UHDYWQ(L,K)-TVAR3E(L)+UHDY2E(L)
     $       +TVAR2N(L,K)-VHDXWQ(L,K)-TVAR3N(L)+VHDX2E(L))*DXYIP(L))
     $        +SWB(L)*( QSUM(L,K)-DZC(K)*QSUME(L) )*DXYIP(L)
        END DO
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        HPPTMP=H2WQ(L)+DT2*DXYIP(L)*(QSUME(L)
     $      -(TVAR3E(L)-UHDY2E(L)+TVAR3N(L)-VHDX2E(L)))
        HWQ(L)=SPB(L)*HPPTMP+(1.-SPB(L))*HWQ(L)
       END DO
      END DO
c
c     add channel interactions
C

      IF (MDCHH.GE.1) THEN
        DO NMD=1,MDCHH
        IF (MDCHTYP(NMD).EQ.1) THEN
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))
     $    +DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))
          HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD))
     $    -DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))
        END IF
        IF (MDCHTYP(NMD).EQ.2) THEN
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))
     $    +DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
          HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD))
     $    -DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))
        END IF
        IF (MDCHTYP(NMD).EQ.3) THEN
          HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))
     $    +DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))
     $    +DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
          HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD))
     $    -DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))
          HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD))
     $    -DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))
        END IF
        END DO
      END IF

	if(ifed_inc.GT.0) call flow_out
c
c     end add channel interactions
C
C      IF(ISTRAN(6).GE.1) CALL CALWQC(2)
      IF(ISTRAN(8).GE.1) 
     & CALL WQ3D(ISTL,N,DT,TCON,TBEGIN,TIDALP,NTSPTC,IWQDT,
     & LA,KC,IC,JC,IWQS,IOP_SAVE,iyear1,NCSTEP,SECDLAST)
     
    !    if(IOP_SAVE.GT.0) Call zwqfield
        
   !       Call zwqfield_bin
     &  

      IF(ISTRAN(4).GE.1) CALL CALSFT(2)        
C--------------------------------------------------------------
C Calculate PCB here 
C---------------------------------------------------------------
      if(ISTRAN(5).GT.0.and.NPCB.GT.0) then
 
       CALL CALPCB(3)                ! transport  2-timestep    
!	 CALL SSEDTOX_PCB(3,CORDT)     ! PCB dynamics

       if(idiaout.gt.0) then
       ACV=86400.0
       ACV1=1/2500.0
       if(mod(N,240).EQ.0) then
        Open(1,file='outpcb.txt',err=9771)
        Open(2,file='outsd.txt',err=9771)
        Open(3,file='outcon.txt',err=9771)
	  write(1,*)'L,X,Y,SEDPOC,TOXb,TOXs,TOXFc, TOXFa,SEDFP'
	  write(2,*)'L,X,Y,TOXB2,TOXB1,PRO2,PRO1,SEDFP1,SEDFP2'
	  write(3,*)'L,X,Y,PBPC,DBPCB,PCBb,PFb,PCBs,PFs,PFsed'
	  DO L=2,LA
!	write(1,11)L,DLON(L),DLAT(L),HP(L)+BELV(L),SAL(L,1),SAL(L,KC)
	 write(1,111)L,DLON(L),DLAT(L),SEDPOC(L,2),TOX(L,1,1),TOX(L,KC,1),
     &	 TOXF(L,0,1)*ACV,TOXF(L,0,2)*ACV,(SEDFP(L,0,1)+SEDFP(L,0,2))*ACV
       write(2,112)L,DLON(L),DLAT(L),TOXB(L,2,1),TOXB(L,1,1),
     &   PORBEDP(L,2),PORBEDP(L,1),SEDFP(L,0,1)*ACV,SEDFP(L,0,2)*ACV

	 FA_1=TOXPFTB(L,2,1)/HBEDP(L,2)/2500.0/(1-PORBEDP(L,2))
	 FA_2=1/HBEDP(L,2)/2500.0/(1-PORBEDP(L,2))
	 write(3,113)L,DLON(L),DLAT(L),TOXB(L,2,1)*FA_2,            ! PPCB: ug/kg  DPCB:ug/m^3
     &   HBEDP(L,2), TOX(L,1,1), TOXPFTW(L,1,1),
     &   TOX(L,KC,1), TOXPFTW(L,KC,1),TOXPFTB(L,2,1)                            ! PCB : ug/m^3
       ENDDO
	 close(1)
	 close(2)
	 close(3)
9771	 continue
	 endif     
	endif
	endif
      
C-----------------------------------------------------------------------
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        H2WQ(L)=HWQ(L)
       END DO
      END DO
C 
      END IF
      END IF
C
C**********************************************************************C
C
C **  UPDATE BUOYANCY AND CALCULATE NEW BUOYANCY USING
C **  AN EQUATION OF STATE
C
      IF (ISTL.EQ.3) THEN
       DO K=1,KC
       DO L=2,LA
       B1(L,K)=B(L,K)
       END DO
       END DO
      END IF
C
      CALL CALBUOY
C
      IF (NCTBC.NE.NTSTBC.AND.ISBAL.GE.1) THEN
         CALL CALBAL4
         NTMP=MOD(N,2)
         IF(NTMP.EQ.0) THEN
           CALL CBALEV4
          ELSE
           CALL CBALOD4
         END IF
      END IF
C
C
C**********************************************************************C
C
C **  CALCULATE U AT V AND V AT U AT TIME LEVEL (N+1)
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      LNW=LNWC(L)
      LSE=LSEC(L)
      LSW=LSWC(L)
      H1C(L)=0.25*(H1P(L)+H1P(L-1)+H1P(LS)+H1P(LSW))
      UV(L)=0.25*(HP(LS)*(U(LSE,1)+U(LS,1))
     $         +HP(L)*(U(L+1,1)+U(L,1)))*HVI(L)
      VU(L)=0.25*(HP(L-1)*(V(LNW,1)+V(L-1,1))
     $         +HP(L)*(V(LN,1)+V(L,1)))*HUI(L)
      END DO
C
C**********************************************************************C
C
C **  CALCULATE HORIZONTAL VISCOSITY AND MOMENTUM DIFFUSION FLUXES
C **  AT TIME LEVEL (N)
C
      IF (ISTL.NE.2.AND.ISHDMF.GE.1) CALL CALHDMF
C
C**********************************************************************C
C
C **  UPDATE BOTTOM STRESSES AND SURFACE AND BOTTOM TURBULENT
C **  INTENSITIES
C
C----------------------------------------------------------------------C
C
C      IF (ISTL.EQ.2) THEN
C
C      DO K=1,KC
C       DO L=2,LA
C        QQ(L,K)=SQRT(QQ(L,K)*QQ1(L,K))
C       END DO
C      END DO
C
C      END IF
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3) THEN
      IF(ISCDMA.EQ.2) THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TBX1(L)=TBX(L)
        TBY1(L)=TBY(L)
        QQ2(L,0)=QQ(L,0)+QQ1(L,0)
        QQ2(L,KC)=QQ(L,KC)+QQ1(L,KC)
        QQ1(L,0)=QQ(L,0)
        QQ1(L,KC)=QQ(L,KC)
       END DO
      END DO
C
      ELSE
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TBX1(L)=TBX(L)
        TBY1(L)=TBY(L)
        QQ2(L,0)=QQ1(L,0)+QQ1(L,0)
        QQ2(L,KC)=QQ1(L,KC)+QQ1(L,KC)
        QQ1(L,0)=QQ(L,0)
        QQ1(L,KC)=QQ(L,KC)
       END DO
      END DO
C
      END IF
      END IF
C
C**********************************************************************C
C
C **  CALCULATE BOTTOM STRESS AT LEVEL (N+1)
C
      CALL CALTBXY
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TBX(L)=(AVCON1*HUI(L)+STBX(L)*SQRT(VU(L)*VU(L)
     $       +U(L,1)*U(L,1)))*U(L,1)
        TBY(L)=(AVCON1*HVI(L)+STBY(L)*SQRT(UV(L)*UV(L)
     $       +V(L,1)*V(L,1)))*V(L,1)
       END DO
      END DO
C
c      IF(ISVEG.GE.1.AND.KC.GT.1)
c      DO ND=1,NDM
c       LF=2+(ND-1)*LDM
c       LL=LF+LDM-1
c       DO L=LF,LL
c       TBX(L)=TBX(L)+0.5*FXVEG(L,1)*SQRT(VU(L)*VU(L)
c     $        +U(L,1)*U(L,1))*U(L,1)
c       TBY(L)=TBY(L)+0.5*FYVEG(L,1)*SQRT(U1V(L)*UV(L)
c     $        +V(L,1)*V(L,1))*V(L,1)
c       END DO
c      END DO
c      END IF
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
C
C----------------------------------------------------------------------C
C
c     IF (KC.GT.1.OR.ISTRAN(4).GE.1) THEN
C
      IF(ISWAVE.EQ.0) THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       QQ(L,0 )=0.5*CTURB2*SQRT((TVAR3E(L)+TBX(L))**2
     $                        +(TVAR3N(L)+TBY(L))**2)
       QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2
     $                         +(TVAR3S(L)+TSY(L))**2)
       END DO
      END DO
C
      END IF
C
c     END IF
C
C**********************************************************************C
C
C **  SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
C
C----------------------------------------------------------------------C
C
c     IF (KC.GT.1.OR.ISTRAN(4).GE.1) THEN
C
      IF(ISWAVE.GE.1) THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       TVAR3S(L)=TSY(LNC(L))
       TVAR3W(L)=TSX(L+1)
       TVAR3E(L)=TBX(L+1   )
       TVAR3N(L)=TBY(LNC(L))
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
         TAUBC2=0.25*( (TVAR3E(L)+TBX(L))**2
     $                        +(TVAR3N(L)+TBY(L))**2 )
         TAUBC=SQRT(TAUBC2)
         UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
         VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))
         CURANG=ATAN2(VTMP,UTMP)
         TAUB2=TAUBC*TAUBC+0.5*(QQWV1(L)*QQWV1(L))
     $        +FOURDPI*TAUBC*QQWV1(L)*COS(CURANG-WACCWE(L))
         TAUB2=MAX(TAUB2,0.)
         QQ(L,0 )=CTURB2*SQRT(TAUB2)
         QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2
     $                         +(TVAR3S(L)+TSY(L))**2)
       END DO
      END DO
C
      END IF
C
c     END IF
C
C**********************************************************************C
C
C **  CALCULATE TURBULENT INTENSITY SQUARED
C
      IF(ISCRAY.EQ.0) THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      END IF
      IF (KC.GT.1) THEN
       IF (ISQQ.EQ.1) CALL CALQQ1 (ISTL)
       IF (ISQQ.EQ.3) CALL CALQQ1a (ISTL)
       IF (ISQQ.EQ.2) CALL CALQQ2 (ISTL)
      END IF
      IF(ISCRAY.EQ.0) THEN
        TQQQ=TQQQ+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TQQQ=TQQQ+T2TMP-T1TMP
        WTQQQ=WTQQQ+(WT2TMP-WT1TMP)*0.001
      END IF
C
C**********************************************************************C
C
C **  HYDRODYNAMIC CALCULATIONS FOR THIS TIME STEP ARE COMPLETED
C **  IF NCTBC EQ NTSTBC APPLY TRAPEZOIDAL CORRECTION
C
C----------------------------------------------------------------------C
C
      IF (NCTBC.EQ.NTSTBC) THEN
       NCTBC=0
       ISTL=2
       DELT=DT
       DELTD2=0.5*DT
       DZDDELT=DZ/DELT
       ROLD=0.5
       RNEW=0.5
       GO TO 500
      ELSE
       NCTBC=NCTBC+1
       ISTL=3
       DELT=DT2
       DELTD2=DT
       DZDDELT=DZ/DELT
       ROLD=0.
       RNEW=1.
      END IF
C
C**********************************************************************C
C
C **  WRITE TO TIME SERIES FILES
C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014
C
      IF (ISTMSR.GE.1) THEN
!        IF (N.GE.NBTMSR.AND.N.LE.NSTMSR) THEN
          IF (NCTMSR.EQ.NWTMSR) THEN
!            if(IGrADS.ne.1) CALL TMSR     ! ji, 10/16/98
           CALL TMSR
            ICALLTP=1
            NCTMSR=1
           ELSE
            NCTMSR=NCTMSR+1
          END IF
!        END IF
      END IF
C
C**********************************************************************C
C
C **  WRITE TO DUMP FILES
C
      IF (ISDUMP.GE.1) THEN
        IF (TIME.GE.TSDUMP.AND.TIME.LE.TEDUMP) THEN
          IF (NCDUMP.EQ.NSDUMP) THEN
!            if(IGrADS.ne.1) CALL DUMP     ! ji, 10/16/98                ! NOT USED
c           CALL DUMP
            ICALLTP=1
            NCDUMP=1
           ELSE
            NCDUMP=NCDUMP+1
          END IF
        END IF
      END IF
C **  WRITE TO hyts.bin, and hy3d.bin for GrADS, jeff ji, 10/16/98
C
      IF (IGrADS.GE.1) CALL DUMP2
C
C**********************************************************************C
C
C ** WRITE BOTTOM VELOCITIES TO FILE FOR AJ MINE SEDIMENT BED MODEL
C     (DLK - ADDED 3/25/97)
      IF (IHYDOUT.GE.1) THEN
        IF (N.GE.NBTMSR.AND.N.LE.NSTMSR) THEN
          IF (NCHYDOUT.EQ.NWHYDOUT) THEN
            CALL HYDOUT
            NCHYDOUT=1
           ELSE
            NCHYDOUT=NCHYDOUT+1
          END IF
        END IF
      END IF
C
C**********************************************************************C
C
C **  OUTPUT ZERO DIMENSION VOLUME BALANCE
C
C----------------------------------------------------------------------C
C
      IF (ISDRY.GE.1.AND.ICALLTP.EQ.1) THEN
        OPEN(1,FILE='zvolbal.out',ACCESS='APPEND',STATUS='UNKNOWN')
        DO LS=1,LORMAX
        IF (VOLZERD.GE.VOLSEL(LS).AND.VOLZERD.LT.VOLSEL(LS+1)) THEN
           WTM=VOLSEL(LS+1)-VOLZERD
           WTMP=VOLZERD-VOLSEL(LS)
           DELVOL=VOLSEL(LS+1)-VOLSEL(LS)
           WTM=WTM/DELVOL
           WTMP=WTMP/DELVOL
           SELZERD=WTM*BELSURF(LS)+WTMP*BELSURF(LS+1)
           ASFZERD=WTM*ASURFEL(LS)+WTMP*ASURFEL(LS+1)
        END IF
        END DO
        TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
       IF(NCSTEP.GT.0) TIME=(SECDLAST+TCON*TBEGIN)/TCTMSR  !% J.S. 1/31/2014 
        WRITE(1,5304) TIME,SELZERD,ASFZERD,VOLZERD,VETZERD
        CLOSE(1)
      END IF
      ICALLTP=0
C
 5304 FORMAT(2X,F10.4,2X,F10.5,3(2X,E12.4))
C
C**********************************************************************C
C
C **  WRITE VERTICAL SCALAR FIELD PROFILES
C
      IF (ISVSFP.EQ.1) THEN
        IF (N.GE.NBVSFP.AND.N.LE.NSVSFP) THEN
          CALL VSFP
        END IF
      END IF
C
C**********************************************************************C
C
C **  CALCULATE MEAN MASS TRANSPORT FIELD
C
      IF(ISSSMMT.NE.2) CALL CALMMT
C
C**********************************************************************C
C
C **  ADVANCE NEUTRALLY BUOYANT PARTICLE DRIFTER TRAJECTORIES
C
c      IF (ISPD.EQ.1) THEN                     ! Ji, use drifter2.for, 4/1/2001
c        IF (N.GE.NPDRT) CALL DRIFTER
c      END IF
      IF (ISPD.EQ.1.and.N.GE.NPDRT) CALL DRIFTER
      IF (ISPD.EQ.2.and.N.GE.NPDRT) CALL DRIFTER2   ! Ji, 4/1/2001
c
      IF (ISLRPD.GE.1) THEN
        IF(ISCRAY.EQ.0) THEN
          T1TMP=SECNDS(0.0)
         ELSE
          T1TMP=SECOND( )
          CALL TIMEF(WT1TMP)
        END IF
        IF (ISLRPD.LE.2) THEN
          IF (N.GE.NLRPDRT(1)) CALL LAGRES
        END IF
        IF (ISLRPD.GE.3) THEN
          IF (N.GE.NLRPDRT(1)) CALL GLMRES
        END IF
        IF(ISCRAY.EQ.0) THEN
          TLRPD=TLRPD+SECNDS(T1TMP)
         ELSE
          T2TMP=SECOND( )
          CALL TIMEF(WT2TMP)
          TLRPD=TLRPD+T2TMP-T1TMP
          WTLRPD=WTLRPD+(WT2TMP-WT1TMP)*0.001
        END IF
      END IF
C
C**********************************************************************C
C
C **  CALCULATE VOLUME MASS, MOMENTUM AND ENERGY BALANCES
C
      IF (ISBAL.GE.1) THEN
         CALL CALBAL5
         NTMP=MOD(N,2)
         IF(NTMP.EQ.0) THEN
           CALL CBALEV5
          ELSE
           CALL CBALOD5
         END IF
       END IF
C
C   SEDIMENT BUDGET CALCULATION     (DLK 10/15)
C
c       IF (ISSBAL.GE.1) THEN  !Ji, 10/21/00
       IF (ISSBAL.GE.1.and.(ISTRAN(6)+ISTRAN(7)).ge.1) THEN
       CALL BUDGET5
       END IF
C       NTMP=MOD(N,2)
C       IF(NTMP.EQ.0) THEN
C         CALL BUDGEV5
C        ELSE
C         CALL BUDGOD5
C       END IF
C
C**********************************************************************C
C
C **  PERFORM AN M2 TIDE HARMONIC ANALYSIS EVERY 2 M2 PERIODS
C
      IF (ISHTA.EQ.1) CALL CALHTA
C
C**********************************************************************C
C
C **  CALCULATE DISPERSION COEFFICIENTS
C
C     IF(N.GE.NDISP) THEN
      IF(N.GE.NDISP.AND.NCTBC.EQ.1) THEN
       IF (ISDISP.EQ.2) CALL CALDISP2
       IF (ISDISP.EQ.3) CALL CALDISP3
      END IF
C
C**********************************************************************C
C
C **  PERFORM LEAST SQUARES HARMONIC ANALYSIS AT SELECTED LOCATIONS
C
      IF (ISLSHA.EQ.1.AND.N.EQ.NCLSHA) THEN
       CALL LSQHARM
       NCLSHA=NCLSHA+(NTSPTC/24)
      END IF
C
C**********************************************************************C
C
C **  PRINT INTERMEDIATE RESULTS
C
C----------------------------------------------------------------------C
C
      IF (NPRINT .EQ. NTSPP) THEN
       NPRINT=1
       CALL OUTPUT1
      ELSE
       NPRINT=NPRINT+1
      END IF
C
C**********************************************************************C
C
C **  WRITE TO TIME VARYING GRAPHICS FILES
C
C----------------------------------------------------------------------C
C
      IF (N.EQ.NCPPH.AND.ISPPH.EQ.1) THEN
       CALL SURFPLT
       NCPPH=NCPPH+(NTSPTC/NPPPH)
      END IF
C
C----------------------------------------------------------------------C
C
      IF (N.EQ.NCVPH.AND.ISVPH.GE.1) THEN
       CALL VELPLTH
       NCVPH=NCVPH+(NTSPTC/NPVPH)
      END IF
C
C----------------------------------------------------------------------C
C
      IF (N.EQ.NCVPV.AND.ISVPV.GE.1) THEN
       CALL VELPLTV
       NCVPV=NCVPV+(NTSPTC/NPVPV)
      END IF
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=1,LC
        TVAR1S(L,K)=TOX(L,K,1)
       END DO
      END DO
C
      IF (N.EQ.NCSPH(1).AND.ISSPH(1).EQ.1) THEN
       IF (ISTRAN(1).GE.1) CALL SALPLTH (1,SAL)
       NCSPH(1)=NCSPH(1)+(NTSPTC/NPSPH(1))
      END IF
C
      IF (N.EQ.NCSPH(2).AND.ISSPH(2).EQ.1) THEN
       IF (ISTRAN(2).GE.1) CALL SALPLTH (2,TEM)
       NCSPH(2)=NCSPH(2)+(NTSPTC/NPSPH(2))
      END IF
C
      IF (N.EQ.NCSPH(3).AND.ISSPH(3).EQ.1) THEN
       IF (ISTRAN(3).GE.1) CALL SALPLTH (3,DYE)
       NCSPH(3)=NCSPH(3)+(NTSPTC/NPSPH(3))
      END IF
C
      IF (N.EQ.NCSPH(4).AND.ISSPH(4).EQ.1) THEN
       IF (ISTRAN(4).GE.1) CALL SALPLTH (4,SFL)
       NCSPH(4)=NCSPH(4)+(NTSPTC/NPSPH(4))
      END IF
C
      IF (N.EQ.NCSPH(5).AND.ISSPH(5).EQ.1) THEN
       IF (ISTRAN(5).GE.1) CALL SALPLTH (5,TVAR1S)
       NCSPH(5)=NCSPH(5)+(NTSPTC/NPSPH(5))
      END IF
C
      IF (N.EQ.NCSPH(6).AND.ISSPH(6).EQ.1) THEN
       IF (ISTRAN(6).GE.1) CALL SALPLTH (6,SEDT)
       NCSPH(6)=NCSPH(6)+(NTSPTC/NPSPH(6))
      END IF
C
      IF (N.EQ.NCSPH(7).AND.ISSPH(7).EQ.1) THEN
       IF (ISTRAN(7).GE.1) CALL SALPLTH (7,SNDT)
       NCSPH(7)=NCSPH(7)+(NTSPTC/NPSPH(7))
      END IF
C
C----------------------------------------------------------------------C
C
      DO ITMP=1,7
      IF (N.EQ.NCSPV(ITMP).AND.ISSPV(ITMP).GE.1) THEN
       CALL SALPLTV(ITMP)
       NCSPV(ITMP)=NCSPV(ITMP)+(NTSPTC/NPSPV(ITMP))
      END IF
      END DO
C
C**********************************************************************C
C
C **  WRITE TO TIME VARYING 3D HDF GRAPHICS FILES
C
C----------------------------------------------------------------------C
C
      IF (N.EQ.NC3DO.AND.IS3DO.EQ.1) THEN
       CALL OUT3D
       NC3DO=NC3DO+(NTSPTC/NP3DO)
      END IF
C
C**********************************************************************C
C
C **  WRITE RESTART FILE EVERY ISRESTO M2 TIDAL CYCLES
C
      IF (ISRESTO.GE.1) THEN
        NRESTO=ISRESTO*NTSPTC
        ISSREST=MOD(N,NRESTO)
        IF (ISSREST.EQ.0) THEN
c          CALL RESTOUT(0)  ! Ji, 10/31/00
          CALL RESTOUT9
          IF (ISTRAN(8).GE.1) THEN
!            IF (IWQRST.EQ.1) CALL WWQRST
!            IF (IWQBEN.EQ.1 .AND. ISMRST.EQ.1) CALL WSMRST
          END IF
        END IF
      END IF
C
C**********************************************************************C
C
C **  RECORD TIME
C
C **  dtime and FLUSH are supported on SUN systems, but may not be
C **  supported on other systems.
C
      IF (NTIMER.EQ.NTSPTC) THEN
      CALL TIMELOG(N)
c     CALL dtime (tarray)
      WRITE(9,200)N, tarray(1),tarray(2)
c     CALL FLUSH(9)
      NTIMER=1
      ELSE
      NTIMER=NTIMER+1
      END IF
C
C**********************************************************************C
C
      IF (ISHOW.EQ.1) CALL SHOWVAL1
      IF (ISHOW.EQ.2) CALL SHOWVAL2
c
c  % Change timestep here J.S. 12/24/2010
c     
c   NCDAY(NLAST)=last timestep, it is set in efdc.for
C
      IF(NCSTEP.GT.0) THEN
       IF(N.EQ.NCDAY(NLAST))THEN
        NTSPTC=NTSPTCC(NLAST)         ! change to new timestep
        DT=TIDALP*FLOAT(NFLTMT)/FLOAT(NTSPTCC(NLAST))
        DTI=1./DT
        DT2=2.*DT
        NSTDAY=N             ! record the last N before changing timestep
        DTD = DT/86400.0     ! wq step in day
        nrefx=NTSPTC*nref
        ISDRY=ISDRYC(NLAST)  !Reset wet-dry
        IF(ISDRY.GE.1) then
         IRVEC=2             ! this need double check for Baylong J.S. 1/31/2014
        ELSE
         IRVEC=9
         DO L=2,LA
          ISCDRY(L)=0
         ENDDO
        ENDIF
        NLAST=NLAST+1
        Iaverage=NTSPTC*aver3d/24.0+0.5
        write(841,*) "2  ",i,Iaverage
        write(*,*)'Change timestep ',NSTDAY,Iaverage,nrefx,NLAST-1
        write(*,*)IAVGBIN ,NTSPTC,DTWQ,SMHSED(1)
!
! Average results CALMMT   ! Fix 2 chose i.e hourly or daily
! Note to ensure to get correct average results, only change timestep by day
!
        if(NTSMMT.LT.NTSPTCC(NLAST-1))THEN
         NTSMMT=NTSPTC/24
        else
         NTSMMT=NTSPTC       
        endif

! Reset TMSR
        Itimes=NTSPTC/24  ! % force output in every hour *savets/24.0+0.5
        NWTMSR=ITIMES
        ICALLTP=1
        NCTMSR=1
! Reset DUMP   ! not used when using efdcwin.inp
!       NSDUMP 
        ICALLTP=1
        NCDUMP=1 
! Reset DUMP2 and DUMP3
!   DUMP2 & DUMP3 use new control (MOD((N-NSTDAY),ITIMES)),
!   NSTDAY is the last N before changing timestep, ITIMES is updated   
c
! WQ   
        IF(ISTRAN(8).GE.1) THEN
         DTWQ = DTD*REAL(IWQDT)  !( IWQDT1) ! this will     
         DTWQO2 = DTWQ*0.5
         IWQTSDT=Itimes  
         WQTSDT=24.0*ITIMES/NTSPTC
         ISMTSDT=NTSPTC !/24*IAVGBIN 
         IAVGBIN=NTSPTC/24*IAVGBIN0 
         IF(IWQBEN.GT.0)then 
          ISMTDMBS1=ISMTDMBS
          ISMTCMBS1=ISMTCMBS
          IWQTSDT = NINT(WQTSDT*3600.0/DT) 
          ISMTSDT = NINT(savets*3600.0/DT)
          ISMTDMBS = NINT(SMTDMBS/DTWQ)
          ISMTCMBS = NINT(SMTCMBS/DTWQ)
          SM1OKMDP = 1.0/SMKMDP
          SMBST1 = 1.0 / (1.0 + SMKBST*DTWQ)
          R1=ISMTDMBS/ISMTDMBS1
          R2=ISMTCMBS/ISMTCMBS1
          DO LL=2,LA
           IF(ISMHYPD(L).GT.0)ISMHYPD(L)= ISMHYPD(L)*R1
          ENDDO
          write(*,*)'ISMH ',ISMTDMBS1,ISMTDMBS,ISMTCMBS1,ISMTCMBS,
     &      ISMHYPD(20)
        
          DO I=1,ISMZ
           SMDTOH(I) = DTWQ/SMHSED(I)
           SMHODT(I) = SMHSED(I)/DTWQ
           SM1DIFT(I) = SMDIFT * SMDTOH(I)/(SMHSED(I)+ 1.E-18)
           SM2DIFT(I) = 1.0 / (1.0 + SM1DIFT(I))
           SMW2DTOH(I) = 1.0 + SMW2(I)*SMDTOH(I)
           SMW2PHODT(I) = SMW2(I) + SMHODT(I)
          ENDDO
         ENDIF
        ENDIF   ! end if WQ parameter
c        
       ENDIF
      ENDIF
C==========
C
      if(mod(N,ISHPRT).EQ.0) then
      Open(1,file='out.txt',err=1000)
!	write(1,*)'L,X,Y,Z,BB,SS'
	DO L=2,LA
!	write(1,11)L,DLON(L),DLAT(L),HP(L)+BELV(L),SAL(L,1),SAL(L,KC)
      IF(NSHTYPE.LE.4) THEN
	write(1,110)L,IL(L),JL(L),SAL(L,1),SAL(L,KC),TIME 
	ELSEIF(NSHTYPE.EQ.5) THEN
	write(1,110)L,IL(L),JL(L),TEM(L,1),TEM(L,KC),TIME 
	ELSE
	write(1,110)L,IL(L),JL(L),WQVO(L,1,NSHTYPE-5),WQVO(L,KC,NSHTYPE-5),TIME 
	ENDIF
	enddo
	close(1)
	endif
11    format(I8,',',F12.3,',',F12.3,',',F9.2,',',F9.2,',',F9.2)
110   format(I8,',',I6,',',I6,',',F9.2,',',F12.4,',',F12.5)
c
C**********************************************************************C
C
 1000 CONTINUE
C
C**********************************************************************C
C
C **  TIME LOOP COMPLETED
C
      IF(ISCRAY.EQ.0) THEN
        THDMT=THDMT+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        THDMT=T2TMP-TTMP
        WTHDMT=(WT2TMP-WTTMP)*0.001
      END IF
C
      WRITE(6,6995)THDMT,TCONG
      WRITE(6,6996)TPUV,TCGRS
      WRITE(6,6997)TCEXP,TAVB
      WRITE(6,6998)TUVW,TQQQ
      WRITE(6,6999)TVDIF,TSADV
      WRITE(6,6994)TLRPD
      WRITE(8,6995)THDMT,TCONG
      WRITE(8,6996)TPUV,TCGRS
      WRITE(8,6997)TCEXP,TAVB
      WRITE(8,6998)TUVW,TQQQ
      WRITE(8,6999)TVDIF,TSADV
      WRITE(8,6994)TLRPD
      WRITE(9,6995)THDMT,TCONG
      WRITE(9,6996)TPUV,TCGRS
      WRITE(9,6997)TCEXP,TAVB
      WRITE(9,6998)TUVW,TQQQ
      WRITE(9,6999)TVDIF,TSADV
      WRITE(9,6994)TLRPD
      WRITE(9,6993)TWQDIF,TWQADV
      WRITE(9,6992)TWQKIN,TWQSED
      WRITE(6,7995)WTHDMT,WTCONG
      WRITE(6,7996)WTPUV,WTCGRS
      WRITE(6,7997)WTCEXP,WTAVB
      WRITE(6,7998)WTUVW,WTQQQ
      WRITE(6,7999)WTVDIF,WTSADV
      WRITE(6,7994)WTLRPD
      WRITE(8,7995)WTHDMT,WTCONG
      WRITE(8,7996)WTPUV,WTCGRS
      WRITE(8,7997)WTCEXP,WTAVB
      WRITE(8,7998)WTUVW,WTQQQ
      WRITE(8,7999)WTVDIF,WTSADV
      WRITE(8,7994)WTLRPD
      WRITE(9,7995)WTHDMT,WTCONG
      WRITE(9,7996)WTPUV,WTCGRS
      WRITE(9,7997)WTCEXP,WTAVB
      WRITE(9,7998)WTUVW,WTQQQ
      WRITE(9,7999)WTVDIF,WTSADV
      WRITE(9,7994)WTLRPD
      WRITE(9,7993)WTWQDIF,WTWQADV
      WRITE(9,7992)WTWQKIN,WTWQSED
C
 6992 FORMAT(' TWQKIN = ',F14.4,'  TWQSED = ',F14.4/)
 6993 FORMAT(' TWQDIF = ',F14.4,'  TWQADV = ',F14.4/)
 6994 FORMAT(' TLRPD = ',F14.4/)
 6995 FORMAT(' THDMT = ',F14.4,'  TCONG = ',F14.4/)
 6996 FORMAT(' TPUV  = ',F14.4,'  TCGRS = ',F14.4/)
 6997 FORMAT(' TCEXP = ',F14.4,'  TAVB  = ',F14.4/)
 6998 FORMAT(' TUVW  = ',F14.4,'  TQQQ  = ',F14.4/)
 6999 FORMAT(' TVDIF = ',F14.4,'  TSADV = ',F14.4/)
C
 7992 FORMAT(' WTWQKIN = ',F14.4,'  WTWQSED = ',F14.4/)
 7993 FORMAT(' WTWQDIF = ',F14.4,'  WTWQADV = ',F14.4/)
 7994 FORMAT(' WTLRPD = ',F14.4/)
 7995 FORMAT(' WTHDMT = ',F14.4,'  WTCONG = ',F14.4/)
 7996 FORMAT(' WTPUV  = ',F14.4,'  WTCGRS = ',F14.4/)
 7997 FORMAT(' WTCEXP = ',F14.4,'  WTAVB  = ',F14.4/)
 7998 FORMAT(' WTUVW  = ',F14.4,'  WTQQQ  = ',F14.4/)
 7999 FORMAT(' WTVDIF = ',F14.4,'  WTSADV = ',F14.4/)
C
C**********************************************************************C
C**********************************************************************C
C
C **  CALCULATE VECTOR POTENTIAL AND VECTOR POTENTIAL TRANSPORTS
C **  USING RESULTS OF THE HARMONIC ANALYSIS
C
c     IF (ISVPTHA.NE.1) GO TO 2000
C
C----------------------------------------------------------------------C
C
c     DO K=1,KC
c     DO L=2,LA
c     LS=LSC(L)
c     VPZ(L,K)=TCVP*SUB(L)*SUB(LS)*SVB(L)*SVB(L-1)*HMC(L)*
c    $         ((AMSU(L,K)+AMSU(LS,K))*(AMCV(L,K)+AMCV(L-1,K))
c    $         -(AMCU(L,K)+AMCU(LS,K))*(AMSV(L,K)+AMSV(L-1,K)))
c     END DO
c     END DO
c
c     DO K=1,KS
c     DO L=2,LA
c     LS=LSC(L)
c     VPX(L,K)=TCVP*((AMSV(L,K)+AMSV(L,K+1))*(AMCW(L,K)+AMCW(LS,K))
c    $              -(AMCV(L,K)+AMCV(L,K+1))*(AMSW(L,K)+AMSW(LS,K)))
c     VPY(L,K)=TCVP*((AMSW(L,K)+AMSW(L-1,K))*(AMCU(L,K)+AMCU(L,K+1))
c    $              -(AMCW(L,K)+AMCW(L-1,K))*(AMSU(L,K)+AMSU(L,K+1)))
c     END DO
c     END DO
C
C----------------------------------------------------------------------C
C
c     DO K=1,KC
c     DO L=2,LA
c     LS=LSC(L)
c     LN=LNC(L)
c     UVPT(L,K)=(VPZ(LN,K)-VPZ(L,K))/DYU(L)-DZI*(VPY(L,K)-VPY(L,K-1))
c     VVPT(L,K)=DZI*(VPX(L,K)-VPX(L,K-1))-(VPZ(L+1,K)-VPZ(L,K))/DXV(L)
c     END DO
c     END DO
C
c     DO K=1,KS
c     DO L=2,LA
c     LS=LSC(L)
c     LN=LNC(L)
c     WVPT(L,K)=(VPY(L+1,K)-VPY(L,K))/DXP(L)
c    $         -(VPX(LN,K)-VPX(L,K))/DYP(L)
c     END DO
c     END DO
C
C----------------------------------------------------------------------C
C
 2000 CONTINUE
C
C**********************************************************************C
C
C **  PRINT FINAL RESULTS
C
      CALL OUTPUT2
C
C**********************************************************************C
C
C **  WRITE RESTART FILE
C
      IF (ISRESTO.EQ.-1.OR.ISRESTO.EQ.-11) THEN
        CALL RESTOUT(0)
        IF (ISTRAN(8).GE.1) THEN
!          IF (IWQRST.EQ.1) CALL WWQRST
!          IF (IWQBEN.EQ.1 .AND. ISMRST.EQ.1) CALL WSMRST
        END IF
      END IF
      IF (ISRESTO.EQ.-2) THEN
        CALL RESTMOD
      END IF
C
C**********************************************************************C
C
C **  COMPLETE LEAST SQUARES HARMONIC ANALYSIS
C
      LSLSHA=1
      IF (ISLSHA.EQ.1) CALL LSQHARM
C
C**********************************************************************C
C
C **  OUTPUT COURANT NUMBER DIAGNOSTICS
C
      OPEN(1,FILE='cflmax.out')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='cflmax.out')
C
      DO L=2,LA
       WRITE(1,1991)IL(L),JL(L),(CFLUUU(L,K),K=1,KC)
       WRITE(1,1992)(CFLVVV(L,K),K=1,KC)
       WRITE(1,1992)(CFLWWW(L,K),K=1,KC)
      END DO
C
      CLOSE(1)
C
 1991 FORMAT(2I5,12F7.2)
 1992 FORMAT(10X,12F7.2)
C
C**********************************************************************C
C
      close(550)
      close(551)
      close(552)      
      close(553) 
      close(554)
      close(555)
              DO NW=1,NTRNVA 
              IDDD=700+NW
              close(IDDD)
              ENDDO
      RETURN
      END
