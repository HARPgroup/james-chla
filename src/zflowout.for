      SUBROUTINE flow_out
!
!  Save 3D flow fileld every FLO_INT
! 
      INCLUDE 'efdc.par'  
      INCLUDE 'efdc.cmn'
      INTEGER, SAVE :: n_re = 0, nofile=0, nvar=0
      INTEGER :: oute=11
      LOGICAL :: op
      LOGICAL, SAVE :: init=.true.
      character*4 tday
!     ifed_inc  = timestep to save the resutls

	TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON
      NM=LCM*KCM*2
	NSD=LCM*KCM*NSCM
	if(n<=2) nofile=int(TBEGIN)/30
	   
      IF(init)THEN
      init=.false.
      write(tday,'(i4.4)')nofile
      write(*,*)'Write binary flow field files !',tday,TBEGIN,nofile
      
      if(nvar.eq.0) then
	open(250,file='dynamicinf.bin',form='unformatted') 
      write(250)LCM,KCM,KSM,NSCM,LA,KC,ICM,JCM
      write(250)LIJ,IL,JL
      write(250)TBEGIN,DT,ifed_inc
      close(250)
      nvar=1
    
      INQUIRE(259,opened=op)
      IF(.not.op)OPEN(259,FILE='infor.out') 
       CLOSE(259,STATUS='DELETE')   
      endif
      
      INQUIRE(251,opened=op)
      IF(.not.op)OPEN(251,FILE='HWQ'//tday//'.outb',
     &        	STATUS='UNKNOWN', form='unformatted')

      INQUIRE(252,opened=op)
      IF(.not.op)OPEN(252,FILE='UHDXYWQ'//tday//'.outb', 
     &           	STATUS='UNKNOWN', form='unformatted')
      
	  INQUIRE(253,opened=op)
      IF(.not.op)OPEN(253,FILE='UVWQ'//tday//'.outb', 
     &         	STATUS='UNKNOWN', form='unformatted')

	  INQUIRE(254,opened=op)
      IF(.not.op)OPEN(254,FILE='WWQ'//tday//'.outb', 
     &    	STATUS='UNKNOWN', form='unformatted')

	  INQUIRE(255,opened=op)
      IF(.not.op)OPEN(255,FILE='AB'//tday//'.outb',
     &   	STATUS='UNKNOWN', form='unformatted')

	  INQUIRE(256,opened=op)
      IF(.not.op)OPEN(256,FILE='TEM'//tday//'.outb', 
     &     	STATUS='UNKNOWN', form='unformatted')

	  INQUIRE(257,opened=op)
      IF(.not.op)OPEN(257,FILE='SED'//tday//'.outb', 
     &    	STATUS='UNKNOWN', form='unformatted')

	  INQUIRE(258,opened=op)
      IF(.not.op)OPEN(258,FILE='SAL'//tday//'.outb', 
     &    	STATUS='UNKNOWN', form='unformatted')

	  INQUIRE(259,opened=op)
      IF(.not.op)OPEN(259,FILE='infor.out',
     &   ACCESS='APPEND',STATUS='UNKNOWN') 
        write(259,*)iyear1,nofile,n_re,TIME
        CLOSE(259)

      END IF

       IF(mod(N,ifed_inc).NE.0) RETURN
       n_re=n_re+1
       write(251)TIME,N,HWQ,HU,HV             ! Elevtion
       WRITE(252)TIME,N,VHDXWQ,UHDYWQ   ! Horizontal transport
       WRITE(253)TIME,N,UWQ,VWQ         ! U,V velocity
       WRITE(254)TIME,N,WWQ             ! W velocity
       WRITE(255)TIME,N,AB              ! Edd difusivity
       WRITE(256)TIME,N,TEM             ! Temperature
       WRITE(257)TIME,N,SEDT            ! Sendiment
       WRITE(258)TIME,N,SAL             ! salt 
!	 WRITE(259,REC=n_re)TIME,N,UUU,VVV
!	 WRITE(260,REC=n_re)TIME,N,WWW
	
	 if(mod(N,NTSPTC*30).eq.0)then   ! 2880 timestep, 30 day 
       nofile=nofile+1
	  do i=251,258
	  close(i)
	  enddo
	  init=.true.
	  n_re=0
	 endif
	
	RETURN
      END  !SUBROUTINE flow_out
C
C======================================================================================================== 
C
      SUBROUTINE flow_out_r
!
!  Save 3D flow fileld every FLO_INT
! 
      INCLUDE 'efdc.par'  
      INCLUDE 'efdc.cmn'
      INTEGER, SAVE :: n_re = 0, nofile=0, nvar=0
      INTEGER :: oute=11
      LOGICAL :: op
      LOGICAL, SAVE :: init=.true.
      character*4 tday

	TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON
      NM=LCM*KCM*2
	NSD=LCM*KCM*NSCM
      
      IF(init)THEN
      init=.false.
      write(tday,'(i4.4)')nofile
      write(*,*)'Write binary flow field files !',tday
      
      if(nvar.eq.0) then
	open(250,file='dynamicinf.bin',form='unformatted') 
      write(250)LCM,KCM,KSM,NSCM,LA,KC,ICM,JCM
      write(250)LIJ,IL,JL
      write(250)TBEGIN,DT,ifed_inc
      close(250)
      nvar=1
      endif
      
      INQUIRE(251,opened=op)
      IF(.not.op)OPEN(251,FILE='HWQ'//tday//'.outb',
     &        	ACCESS='direct',RECL=4*(2+LCM*3))

      INQUIRE(252,opened=op)
      IF(.not.op)OPEN(252,FILE='UHDXYWQ'//tday//'.outb', 
     &           ACCESS='direct',RECL=4*(2+NM))
      
	  INQUIRE(253,opened=op)
      IF(.not.op)OPEN(253,FILE='UVWQ'//tday//'.outb', 
     &           ACCESS='direct',RECL=4*(2+NM))

	  INQUIRE(254,opened=op)
      IF(.not.op)OPEN(254,FILE='WWQ'//tday//'.outb', 
     &        ACCESS='direct',RECL=4*(2+LCM*(KCM+1)))

	  INQUIRE(255,opened=op)
      IF(.not.op)OPEN(255,FILE='AB'//tday//'.outb',
     &        ACCESS='direct',RECL=4*(2+LCM*KSM))

	  INQUIRE(256,opened=op)
      IF(.not.op)OPEN(256,FILE='TEM'//tday//'.outb', 
     &        ACCESS='direct',RECL=4*(2+LCM*KCM))

	  INQUIRE(257,opened=op)
      IF(.not.op)OPEN(257,FILE='SED'//tday//'.outb', 
     &        ACCESS='direct',RECL=4*(2+LCM*KCM))

	  INQUIRE(258,opened=op)
      IF(.not.op)OPEN(258,FILE='SAL'//tday//'.outb', 
     &        ACCESS='direct',RECL=4*(2+LCM*KCM))


      END IF

       IF(mod(N,ifed_inc).NE.0) RETURN
       n_re=n_re+1
       write(251,REC=n_re)TIME,N,HWQ,HU,HV             ! Elevtion
       WRITE(252,REC=n_re)TIME,N,VHDXWQ,UHDYWQ   ! Horizontal transport
       WRITE(253,REC=n_re)TIME,N,UWQ,VWQ         ! U,V velocity
       WRITE(254,REC=n_re)TIME,N,WWQ             ! W velocity
       WRITE(255,REC=n_re)TIME,N,AB              ! Edd difusivity
       WRITE(256,REC=n_re)TIME,N,TEM             ! Temperature
       WRITE(257,REC=n_re)TIME,N,SEDT             ! Sendiment
       WRITE(258,REC=n_re)TIME,N,SAL             ! salt 

	
	 if(mod(N,NTSPTC*30).eq.0)then   ! 2880 timestep, 30 day 
        nofile=nofile+1
	  do i=251,258
	  close(i)
	  enddo
	  init=.true.
	  n_re=0
	 endif
	
	RETURN
      END  !SUBROUTINE flow_out
C
C==========================================================
C
      SUBROUTINE flow_out_a
!
!  Save 3D flow fileld every FLO_INT
! 
      INCLUDE 'efdc.par'  
      INCLUDE 'efdc.cmn'
      INTEGER, SAVE :: n_re = 0, nofile=0,nvar=0
      INTEGER :: oute=11
      LOGICAL :: op
      LOGICAL, SAVE :: init=.true.
      character*4 tday

	TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON
      NM=LCM*KCM*2
	NSD=LCM*KCM*NSCM
      IF(init)THEN
      init=.false.
      write(tday,'(i4.4)')nofile
      write(*,*)'Write binary flow field files !',tday
      if(nvar.eq.0) then
	open(250,file='dynamicinf.bin',form='unformatted') 
      write(250)LCM,KCM,KSM,NSCM,LA,KC,ICM,JCM
      write(250)LIJ,IL,JL
      write(250)TBEGIN,DT,ifed_inc
      close(250)
      nvar=1
      endif

      INQUIRE(251,opened=op)
      IF(.not.op)OPEN(251,FILE='HWQ'//tday//'.outb',
     &        	 FORM='UNFORMATTED')

      INQUIRE(252,opened=op)
      IF(.not.op)OPEN(252,FILE='UHDXYWQ'//tday//'.outb', 
     &           FORM='UNFORMATTED')
      
	  INQUIRE(253,opened=op)
      IF(.not.op)OPEN(253,FILE='UVWQ'//tday//'.outb', 
     &           FORM='UNFORMATTED')

	  INQUIRE(254,opened=op)
      IF(.not.op)OPEN(254,FILE='WWQ'//tday//'.outb', 
     &         FORM='UNFORMATTED')

	  INQUIRE(255,opened=op)
      IF(.not.op)OPEN(255,FILE='AB'//tday//'.outb',
     &         FORM='UNFORMATTED')

	  INQUIRE(256,opened=op)
      IF(.not.op)OPEN(256,FILE='TEM'//tday//'.outb', 
     &         FORM='UNFORMATTED')

	  INQUIRE(257,opened=op)
      IF(.not.op)OPEN(257,FILE='SED'//tday//'.outb', 
     &        FORM='UNFORMATTED')

	  INQUIRE(258,opened=op)
      IF(.not.op)OPEN(258,FILE='SAL'//tday//'.outb', 
     &         FORM='UNFORMATTED')

!	  INQUIRE(259,opened=op)
!      IF(.not.op)OPEN(259,FILE='UUU.outb', 
!    &        ACCESS='direct',RECL=4*(2+NM))

!	  INQUIRE(260,opened=op)
!      IF(.not.op)OPEN(260,FILE='WWW.outb', 
!     &        ACCESS='direct',RECL=4*(2+LCM*(KCM+1)))

      END IF

       IF(mod(N,ifed_inc).NE.0) RETURN
       n_re=n_re+1
       write(251)TIME,N,HWQ,HU,HV             ! Elevtion
       WRITE(252)TIME,N,VHDXWQ,UHDYWQ   ! Horizontal transport
       WRITE(253)TIME,N,UWQ,VWQ         ! U,V velocity
       WRITE(254)TIME,N,WWQ             ! W velocity
       WRITE(255)TIME,N,AB              ! Edd difusivity
       WRITE(256)TIME,N,TEM             ! Temperature
       WRITE(257)TIME,N,SED             ! Sendiment
       WRITE(258)TIME,N,SAL             ! salt 
!	 WRITE(259,REC=n_re)TIME,N,UUU,VVV
!	 WRITE(260,REC=n_re)TIME,N,WWW
	
	 if(mod(N,2880*30).eq.0)then   ! 2880 timestep, 30 day 
       nofile=nofile+1
	  do i=251,258
	  close(i)
	  enddo
	  init=.true.
	  n_re=0
	 endif
	
	RETURN
      END  !SUBROUTINE flow_out
