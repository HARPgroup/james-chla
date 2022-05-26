      SUBROUTINE zwqfield_bin
      INCLUDE 'wq.par' 
      INCLUDE 'efdc.par'
      INCLUDE 'wqcom.cmn'
      INCLUDE 'efdc.cmn'
      LOGICAL, SAVE :: initwq=.true. 
      INTEGER, SAVE :: nwq_re = 0,iopen=0,hrsave=1
	REAL, SAVE :: T_1,T_2,Dbeg ,RCDAY0   
	REAL :: lmd_1,lmd_2  
      LOGICAL :: op
 
      character*4 tday
      character*60 pathwq
            pathwq=PATHpcb
	TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON     
      TIMEhr=DT*FLOAT(N)/TCON*24                 ! get hour
      
      IF(NCSTEP.GT.0) THEN
      TIME=SECDLAST/TCON+TBEGIN                  !  J.S. 1/31/2014  
      TIMEhr=SECDLAST/TCON*24                    !  get hour
      ENDIF
      
      IF(initwq)THEN
      initwq=.false.
      write(tday,'(i4.4)')irecdwq !nfin
      write(*,*)'Read/write WQ inary files !',tday
      
      INQUIRE(261,opened=op)
      IF(.not.op)OPEN(261,FILE=trim(pathwq)//'WQAG'//tday//'.outb',
     &        	STATUS='UNKNOWN', form='unformatted')

      INQUIRE(262,opened=op)
      IF(.not.op)OPEN(262,FILE=trim(pathwq)//'WQPOC'//tday//'.outb',
     &        	STATUS='UNKNOWN', form='unformatted')

c	INQUIRE(263,opened=op)     
c      IF(.not.op)OPEN(263,FILE=trim(pathwq)//'WQSEDPOC'//tday//'.outb',
c     &        	STATUS='UNKNOWN', form='unformatted')
      INQUIRE(263,opened=op)
      IF(.not.op)OPEN(263,FILE=trim(pathwq)//'WQTP'//tday//'.outb',
     &        	STATUS='UNKNOWN', form='unformatted')
     
      INQUIRE(264,opened=op)
      IF(.not.op)OPEN(264,FILE=trim(pathwq)//'WQTN'//tday//'.outb',
     &        	STATUS='UNKNOWN', form='unformatted')
      
      ENDIF
     
      IF(IOP_SAVE.GT.0) THEN  ! Save results hourly
  !    IF(mod(N,IOP_SAVE).EQ.0.OR.N.EQ.2)THEN
      IF(abs(TIMEhr-hrsave).LT.0.00001) THEN
       if(mod(hrsave,24).eq.0) write(*,*)'Save WQ ',TIMEhr, hrsave,TIME     ! save hourly
       hrsave=hrsave+1
       nwq_re=nwq_re+1
       DO K=1,KC
       DO L=2,LA
        PCBAGD(L,K)=WQV(L,K,2)            ! diaton
        PCBAGG(L,K)=WQV(L,K,3)+WQV(L,K,1) ! green
       END DO
       END DO
       WRITE(261)TIME,N,PCBDOC,PCBAGD,PCBAGG         ! Algae
       DO K=1,KC
       DO L=2,LA
        PCBLPC(L,K)=WQV(L,K,4)+ WQV(L,K,5)           ! particulate carbon
        PCBDOC(L,K)=WQV(L,K,6)                       ! dissolved carbon
       END DO
       END DO     
       DO L=1,LA
         SE_POC(L,1)=SMPOC(L,1)  ! sediment particula C  12/6/2017 JS
	   SE_POC(L,2)=SMPOC(L,2)
         SE_POC(L,3)=SMPOC(L,3)
	 ENDDO 
       WRITE(262)TIME,N,PCBLPC,PCBDOC,SE_POC                      ! Sediment POC 
       DO K=1,KC
       DO L=2,LA
        PCBLPC(L,K)=WQV(L,K,7)+ WQV(L,K,8)+ WQV(L,K,9)           ! OP
        PCBDOC(L,K)=WQV(L,K,10)                                  ! PO4
       END DO
       END DO 
       DO L=1,LA 
        SE_POC(L,1)=SMPOP(L,1)                                    ! sediment PON
	  SE_POC(L,2)=SMPOP(L,2)
        SE_POC(L,3)=SMPOP(L,3)
       ENDDO 
       WRITE(263)TIME,N,PCBLPC,PCBDOC         
       DO K=1,KC
       DO L=2,LA
        PCBLPC(L,K)=WQV(L,K,11)+WQV(L,K,12)+WQV(L,K,13)         ! ON
        PCBDOC(L,K)=WQV(L,K,14)+WQV(L,K,15)                      ! dissolved carbon
       END DO
       END DO 
       DO L=1,LA 
        SE_POC(L,1)=SMPON(L,1)                                    ! sediment PON
	  SE_POC(L,2)=SMPON(L,2)
        SE_POC(L,3)=SMPON(L,3)
       ENDDO
       WRITE(264)TIME,N,PCBLPC,PCBDOC,SE_POC                
      ENDIF
      ENDIF
C==================================================================================
      IF(IOP_SAVE.LT.0) THEN
      if(iopen.eq.0) then
      iopen=1
      read(261)T_2,N_2,PCBAGD,PCBAGG   
      read(262)T_2,N_2,PCBLPC,PCBDOC,SE_POC 
      PCBDOC1=PCBDOC
      PCBAGG1=PCBDOC 
      PCBAGG1=PCBAGG
      PCBLPC1=PCBLPC
      PCBDOC1=PCBDOC
      SE_POC1=SE_POC 
      T_1=T_2-1.0/24.0   
      nwq_re=nwq_re+1      
      endif
      
      IF(TIME.GT.T_2) THEN              
       DO K=1,KC
        DO L=2,LA
        PCBAGD1(L,K)=PCBAGD(L,K)          ! particulate carbon
        PCBAGG1(L,K)=PCBAGG(L,K) 
        PCBLPC1(L,K)=PCBLPC(L,K)
        PCBDOC1(L,K)=PCBDOC(L,K)	                 
        END DO
       END DO     
       DO K=1,3
	  DO L=2,LA
         SE_POC1(L,K)=SE_POC(L,K)
	  ENDDO
	 ENDDO
!	 write(*,*)'nwq_re ==',nwq_re
        nwq_re=nwq_re+1 !inc_re
	  T_1=T_2
        read(261)T_2,N_2,PCBAGD,PCBAGG        
        read(262)T_2,N_2,PCBLPC,PCBDOC,SE_POC          	
!        T_2=T_2+ Dbeg
       ENDIF
        if(mod(N,1440).eq.0) write(*,*)'Input wq record at ',T_1, T_2
       
        t_dff=T_2-T_1 
        lmd_1=(TIME-T_1)/t_dff
	  lmd_2=(T_2-TIME)/t_dff
        if(lmd_1.LT.0) then
	   lmd_1=0.0
	   lmd_2=1
	  endif
C
C     Get results using linear interpolation
C
        DO  K=1,KC
         DO L=2,LA
         WQV(L,K,1)=0
         WQV(L,K,3)=PCBAGG1(L,K)*lmd_2+PCBAGG(L,K)*lmd_1 
         WQV(L,K,2)=PCBAGD1(L,K)*lmd_2+PCBAGD(L,K)*lmd_1 
         WQV(L,K,4)=0	    
!         WQV(L,K,5)=PCBLPC1(L,K)*lmd_2+PCBLPC(L,K)*lmd_1 
!         WQV(L,K,6)=PCBDOC1(L,K)*lmd_2+PCBDOC(L,K)*lmd_1
         WQV(L,K,5)=PCBLPC1(L,K)*lmd_2+PCBLPC(L,K)*lmd_1 
         WQV(L,K,6)=PCBDOC1(L,K)*lmd_2+PCBDOC(L,K)*lmd_1         
!         WQ_SET_C(L,K)=PCBCSET1(L,K)*lmd_2+PCBCSET(L,K)*lmd_1
!         WQ_SET_G(L,K)=PCBGSET1(L,K)*lmd_2+PCBGSET(L,K)*lmd_1
         ENDDO
        ENDDO
        
	  DO  K=1,3
         DO L=2,LA
		 SMPOC(L,K)=SE_POC1(L,K)*lmd_2+SE_POC(L,K)*lmd_1  
!		 SMPOC(L,K)=max(SMPOC(L,K),20.0)
         SMPOC(L,K)=max(SMPOC(L,K),1.0)
	   ENDDO
	  ENDDO
c
c   Output diagnostics
c   
       if(idiaout.eq.2) then
        OPEN(266,FILE='wqdia.out',STATUS='UNKNOWN',ACCESS='APPEND')
        DO M=1,IWQTS                              ! #200
         DO K=1,KC                                ! #100
          LL=LWQTS(M)	  
          write(266,2660)Time, LL,K,PCBLPC(LL,K),WQV(LL,K,5),
     &    PCBDOC(LL,K),WQV(LL,K,6), PCBAGG(LL,K),WQV(LL,K,3),
     &    SE_POC(LL,1),SE_POC(LL,2),SE_POC(LL,3)
          ENDDO
        ENDDO
       close(266)
       endif  
      
      ENDIF
 2660  format(F8.2,2I7,20F10.4)
C==================================================================================
       
      if(nwq_re .eq. 24*30)then            ! save hourly for 30 days
	  !  open new file next timestep
	  do i=261,264
	  close(i)
	  enddo
	  initwq=.true.
	  nwq_re = 0
	  irecdwq=irecdwq+1
	  write(*,*)'Open new files ',irecdwq
	 endif
      Return
      END
C======================================================================
      SUBROUTINE zwqfield
!
!  Save 3D flow fileld every FLO_INT
! 
!      IOP_SAVE > 0 save resuls
!               < 0 read resuls
!      inc_re   = hour (=1 save one hour   
!      irec     = number of saved record which will not be used for the current run 
!                 Note that it need to figure out based on save frequency and start of the model run
!      isWQFELD save resutls when mod(N,isWQFELD).EQ.0
!
      INCLUDE 'wq.par' 
      INCLUDE 'efdc.par'
      INCLUDE 'wqcom.cmn'
      INCLUDE 'efdc.cmn'
      INTEGER, SAVE :: nwq_re = 0,iopen=0
	REAL, SAVE :: T_1,T_2,Dbeg ,RCDAY0   
	REAL :: lmd_1,lmd_2  
      INTEGER :: oute=11
      LOGICAL :: op
      LOGICAL, SAVE :: init=.true.
      
      nrecs=0  !24*TBEGIN/real(inc_re)
	TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON
      NM=LCMWQ*KCWM*2

      IF(init)THEN
      init=.false.

!      if(IOP_SAVE.GT.0)then
!      write(*,*)'Write WQ carbon !'
!	 open(250,file='wqfield.bin',form='unformatted') 
!      write(250)LCMWQ,ICMWQ,JCMWQ,KCWM
!      write(250)LIJW,ILW,JLW
!      write(250)TBEGIN,DT,isWQDF_inc
!      close(250)
!      endif
      
!      OPEN(266,FILE='wqdia.out',STATUS='UNKNOWN')
!      CLOSE(266,STATUS='DELETE')
!      OPEN(266,FILE='wqdia.out',STATUS='UNKNOWN')
!      WRITe(266,*)'Time LL K PCBLPC,WQV(LL,K,5)'
!      WRITE(266,*)'PCBDOC,WQV(6) PCBAGG WQV(3)'
!      WRITE(266,*)'SE_POC(LL,1) SE_POC(LL,2) SE_POC(LL,3)'
!      close(266)
      
      INQUIRE(261,opened=op)
      IF(.not.op)OPEN(261,FILE='WQAG.outb',
     &        	ACCESS='direct',RECL=4*(2+NM))

      INQUIRE(262,opened=op)
      IF(.not.op)OPEN(262,FILE='WQPOC.outb', 
     &           ACCESS='direct',RECL=4*(2+NM))

	INQUIRE(263,opened=op)
      IF(.not.op)OPEN(263,FILE='WQSEDPOC.outb', 
     &           ACCESS='direct',RECL=4*(2+LCMWQ*3))
      
      ENDIF
      
      if(iopen.eq.0) then
       iopen=1
        RCDAY0=RCDAY
        Dbeg=0  
	 IF(IOP_SAVE.LE.0.and.nwq_re.eq.0) then   ! Read results first time
                
        nwq_re = nrecs + inc_re                       ! count for record, forced to increase hourly
        read(261,REC=nwq_re)T_1,N_1,PCBAGD1,PCBAGG1        
        read(262,REC=nwq_re)T_1,N_1,PCBLPC1,PCBDOC1     
	  read(263,REC=nwq_re)T_1,N_1,SEDPOC1    
!       read(263,REC=n_re)T_1,N_,PCBCSET1,PCBGSET1      
        nwq_re=nwq_re+inc_re
        read(261,REC=nwq_re)T_2,N_2,PCBAGD,PCBAGG        
        read(262,REC=nwq_re)T_2,N_2,PCBLPC,PCBDOC         
	  read(263,REC=nwq_re)T_2,N_2,SE_POC   
!      read(263,REC=n_re)T_2,N_2,PCBCSET,PCBGSET 
       write(*,*)'Input wq data at the 1st time ',T_1,T_2,nwq_re
       END IF
      endif  
             
      IF(IOP_SAVE.GT.0) THEN
       IF(mod(N,IOP_SAVE).EQ.0.OR.N.EQ.2)THEN
         nwq_re=nwq_re+1
        DO K=1,KC
        DO L=2,LA
         PCBAGD(L,K)=WQV(L,K,2)   ! diaton
         PCBAGG(L,K)=WQV(L,K,3)   ! green
        END DO
        END DO
        WRITE(261,REC=nwq_re)TIME,N,PCBAGD,PCBAGG         ! Algae
        DO K=1,KC
        DO L=2,LA
         PCBLPC(L,K)=WQV(L,K,4)+ WQV(L,K,5)           ! particulate carbon
         PCBDOC(L,K)=WQV(L,K,6)                       ! dissolved carbon
        END DO
        END DO    
        WRITE(262,REC=nwq_re)TIME,N,PCBLPC,PCBDOC         ! PC,DOC      
        DO L=1,LA
          SE_POC(L,1)=SMPOC(L,1)
	    SE_POC(L,2)=SMPOC(L,2)
          SE_POC(L,3)=SMPOC(L,3)
	  ENDDO	 
         WRITE(263,REC=nwq_re)TIME,N,SE_POC           ! Sediment POC 
        if(mod(nwq_re,1440).eq.0)write(*,*)'Out wq ',TIME,N,nwq_re
       ELSE
         Return
       ENDIF
      ENDIF
      
      IF(IOP_SAVE.LT.0) THEN
      
       IF(TIME.GT.T_2) THEN              
        DO K=1,KC
         DO L=2,LA
         PCBAGD1(L,K)=PCBAGD(L,K)          ! particulate carbon
         PCBAGG1(L,K)=PCBAGG(L,K) 
         PCBLPC1(L,K)=PCBLPC(L,K)
         PCBDOC1(L,K)=PCBDOC(L,K)	     
!         PCBCSET1(L,K)= PCBCSET(L,K)
!         PCBGSET1(L,K)= PCBGSET(L,K)             
         END DO
        END DO     
        DO K=1,3
	   DO L=2,LA
	    SE_POC1(L,K)=SE_POC(L,K)
	   ENDDO
	  ENDDO
        nwq_re=nwq_re+inc_re
	  T_1=T_2
        read(261,REC=nwq_re)T_2,N_2,PCBAGD,PCBAGG        
        read(262,REC=nwq_re)T_2,N_2,PCBLPC,PCBDOC         
        read(263,REC=nwq_re)T_2,N_2,SE_POC	
        T_2=T_2+ Dbeg
        if(mod(nwq_re,240).eq.0)write(*,*)'WQ input T_2 = ',T_2,N
       ENDIF
        if(mod(N,1440).eq.0) write(*,*)'Input wq record at ',T_2,nwq_re       
       
        t_dff=T_2-T_1 
        lmd_1=(TIME-T_1)/t_dff
	  lmd_2=(T_2-TIME)/t_dff
        if(lmd_1.LT.0) then
	   lmd_1=0.0
	   lmd_2=1
	  endif

C
C     Get results using linear interpolation
C
        DO  K=1,KC
         DO L=2,LA
         WQV(L,K,1)=0
         WQV(L,K,3)=PCBAGG1(L,K)*lmd_2+PCBAGG(L,K)*lmd_1 
         WQV(L,K,2)=PCBAGD1(L,K)*lmd_2+PCBAGD(L,K)*lmd_1 
         WQV(L,K,4)=0	    
!         WQV(L,K,5)=PCBLPC1(L,K)*lmd_2+PCBLPC(L,K)*lmd_1 
!         WQV(L,K,6)=PCBDOC1(L,K)*lmd_2+PCBDOC(L,K)*lmd_1
         WQV(L,K,5)=PCBLPC1(L,K)*lmd_2+PCBLPC(L,K)*lmd_1 
         WQV(L,K,6)=PCBDOC1(L,K)*lmd_2+PCBDOC(L,K)*lmd_1         
!         WQ_SET_C(L,K)=PCBCSET1(L,K)*lmd_2+PCBCSET(L,K)*lmd_1
!         WQ_SET_G(L,K)=PCBGSET1(L,K)*lmd_2+PCBGSET(L,K)*lmd_1
         ENDDO
        ENDDO
        
	  DO  K=1,3
         DO L=2,LA
		 SMPOC(L,K)=SE_POC1(L,K)*lmd_2+SE_POC(L,K)*lmd_1  
		 SMPOC(L,K)=max(SMPOC(L,K),20.0)
	   ENDDO
	  ENDDO
c
c   Output diagnostics
c   
       if(idiaout.eq.1) then
        OPEN(266,FILE='wqdia.out',STATUS='UNKNOWN',ACCESS='APPEND')
        DO M=1,IWQTS                              ! #200
         DO K=1,KC                                ! #100
          LL=LWQTS(M)	  
          write(266,2660)Time, LL,K,PCBLPC(LL,K),WQV(LL,K,5),
     &    PCBDOC(LL,K),WQV(LL,K,6), PCBAGG(LL,K),WQV(LL,K,3),
     &    SE_POC(LL,1),SE_POC(LL,2),SE_POC(LL,3)
          ENDDO
        ENDDO
       close(266)
       endif
      ENDIF
      
!     if(mod(N,14400).eq.0)then
!	  do i=261,263
!	  close(i)
!	  enddo
!	  init=.true.
!	endif
      if(IOP_SAVE.LE.0.and.RCDAY.GT.1)then
	if((TIME-TBEGIN).GE.RCDAY0-0.00001)then
	  do i=261,264
	  close(i)
	  enddo
	  init=.true.
        nwq_re=nrecs
	  Dbeg=Dbeg+RCDAY
	  RCDAY0=RCDAY0+RCDAY
	endif 
	endif 
 2660  format(F8.2,2I7,20F10.4)

	RETURN
      END  !SUBROUTINE flow_out
