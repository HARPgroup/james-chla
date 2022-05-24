      SUBROUTINE flow_in
!
!  Save 3D flow fileld every FLO_INT
!  Last Version 2/24/2014 by J.S.
!
      INCLUDE 'efdc.par'  
      INCLUDE 'efdc.cmn'
      include 'wq.par'
      INCLUDE 'wqcom.cmn'
      INTEGER, SAVE :: n_re = 0,N_1,N_2,nopen=0,nfin=0
	REAL, SAVE :: t_1,t_2,Dbeg,RCDAY0,TBg_save
	REAL :: lmd_1,lmd_2,t_t11,t_12,t_13,t_14
      INTEGER :: oute=11
      LOGICAL :: op
      LOGICAL, SAVE :: init=.true.
      character*4 tday
!      character*10 path
!      path='..\Hydro1\'
!      path='   Hydro1\' 
C       
      iskp=1 
      nrecs=24*TBEGIN*real(inc_re)
      if(mod(int(TBEGIN),30).eq.0)nrecs=0
	TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014  
       
      NM=LCM*KCM*2
	NSD=LCM*KCM*NSCM

!      if(nopen.eq.0)irecd=TBEGIN/30  Read from input Card 92

      IF(init)THEN
      init=.false.
      
      write(tday,'(i4.4)')irecd !nfin
      write(*,*)'Read binary flow field files !',tday

      INQUIRE(251,opened=op)
      IF(.not.op)OPEN(251,FILE=trim(path)//'HWQ'//tday//'.outb',
     &       	STATUS='UNKNOWN', form='binary')

      INQUIRE(252,opened=op)
      IF(.not.op)OPEN(252,FILE=trim(path)//'UHDXYWQ'//tday//'.outb', 
     &          	STATUS='UNKNOWN', form='binary')
      
	  INQUIRE(253,opened=op)
      IF(.not.op)OPEN(253,FILE=trim(path)//'UVWQ'//tday//'.outb', 
     &       STATUS='UNKNOWN', form='binary')

	  INQUIRE(254,opened=op)
      IF(.not.op)OPEN(254,FILE=trim(path)//'WWQ'//tday//'.outb', 
     &    	STATUS='UNKNOWN', form='binary')
     
	  INQUIRE(255,opened=op)
      IF(.not.op)OPEN(255,FILE=trim(path)//'AB'//tday//'.outb',
     &       	STATUS='UNKNOWN', form='binary')
     
      if(ISTRAN(2).eq.1)then
	  INQUIRE(256,opened=op)
      IF(.not.op)OPEN(256,FILE=trim(path1)//'TEM'//tday//'.outb' 
     &      	,STATUS='UNKNOWN', form='binary')
      endif
      
      if(ISTRAN(6).eq.1) then
	  INQUIRE(257,opened=op)
      IF(.not.op)OPEN(257,FILE=trim(path1)//'SED'//tday//'.outb' 
     &       	,STATUS='UNKNOWN', form='binary')
      endif
      
      if(ISTRAN(1).eq.1)then
	  INQUIRE(258,opened=op)
      IF(.not.op)OPEN(258,FILE=trim(path)//'SAL'//tday//'.outb', 
     &      	STATUS='UNKNOWN', form='binary')
      endif
      
       if(nopen.eq.0)then
C !  Read 2 record for starting
        nopen=1
        RCDAY0=RCDAY
        Dbeg=0 
        n_re = n_re +1  ! Save half hour            

!       write(251)TIME,N,HWQ,HU,HV             ! Elevtion
!       WRITE(252)TIME,N,VHDXWQ,UHDYWQ   ! Horizontal transport
!       WRITE(253)TIME,N,UWQ,VWQ         ! U,V velocity
!       WRITE(254)TIME,N,WWQ             ! W velocity
!       WRITE(255)TIME,N,AB              ! Edd difusivity
!       WRITE(256)TIME,N,TEM             ! Temperature
!       WRITE(257)TIME,N,SEDT            ! Sendiment
!       WRITE(258)TIME,N,SAL             ! salt 
       
       read(251)t_1,N_1,H2P,U1V,V1U    ! Elevtion
       read(252)t_1,N_1,VHDX2,UHDY2    ! Horizontal transport
       read(253)t_1,N_1,U2,V2          ! U,V velocity
       read(254)t_1,N_1,W2             ! W velocity
       read(255)t_1,N_1,ABLPF             ! Edd difusivity
       if(ISTRAN(2).eq.1)read(256)t_1,N_1,TEM1           ! Temperature
       if(ISTRAN(6).eq.1)read(257)t_1,N_1,SEDTLPF         ! Sendiment
       if(ISTRAN(1).eq.1)read(258)t_1,N_1,SAL1  
       TBg_save=t_1 
       IF(IRELSOP.EQ.5)t_1=NTC+30-t_1
111    n_re=n_re+1
       read(251)t_2,N_2,H1P,H1U,H1V             ! Elevtion
       read(252)t_2,N_2,VHDX1,UHDY1    ! Horizontal transport
       read(253)t_2,N_2,U1,V1          ! U,V velocity
       read(254)t_2,N_2,W1             ! W velocity
       read(255)t_2,N_2,ABEFF          ! Edd difusivity
       if(ISTRAN(2).eq.1)read(256)t_2,N_2,TEMLPF         ! Temperature
       if(ISTRAN(6).eq.1)read(257)t_2,N_2,SNDTLPF        ! Sendiment
       if(ISTRAN(1).eq.1)read(258)t_2,N_2,SALLPF   

        IF(IRELSOP.EQ.5)t_2=NTC+30-t_2

       DO L=1,LA
        HWQ(L)=H2P(L)       
       ENDDO
       
       if(TIME.GT.t_2) THEN   ! Take care restart to fine current recored. 
        write(*,*)'Start hydro enter ',TIME, t_2,n_re
       	DO L=1,LA
        H2P(L)=max(H1P(L),0.3)
        U1V(L)=H1U(L)
	  V1U(L)=H1V(L) 
	  if(abs(H2P(L)- U1V(L)).GT.0.5.and.abs(V1U(L)-H2P(L)).Gt.0.5)then
!	  write(*,*)'Check HUV  ',L,H2P(L),U1V(L),V1U(L)
	  endif
!	  HWQ(L)=max(H1P(L),0.3)
	  ENDDO

	  DO K=1,KC 
	  DO L=1,LA 
        VHDX2(L,K)=VHDX1(L,K)
	  UHDY2(L,K)=UHDY1(L,K)
!	  ABLPF(L,K)=ABEFF(L,K)
	  TEM1(L,K)=TEMLPF(L,K)
	  SEDTLPF(L,K)=SNDTLPF(L,K)
	  SAL1(L,K)=SALLPF(L,K)
	  U2(L,K)=U1(L,K)
	  V2(L,K)=V1(L,K)
	  ENDDO
	  ENDDO

	  DO K=0,KC
	  DO L=1,LA
        W2(L,K)=W1(L,K)
        ENDDO
	  ENDDO

	  DO K=1,KC-1 
	  DO L=1,LA 
	  ABLPF(L,K)=ABEFF(L,K)
	  ENDDO
	  ENDDO
  	  t_1=t_2
  	  goto 111  	  
  	 endif
      endif

      END IF ! end IF(init)

      IF(TIME.GT.t_2) THEN
  !      write(*,*)'Time == ',TIME
	  DO L=2,LA
        H2P(L)=H1P(L)
        U1V(L)=H1U(L)
	  V1U(L)=H1V(L) 
!	  DO II=1,IWQPS
!	  if(L.EQ.LIJ(ICPSL(II),JCPSL(II))) goto 1278
!	  ENDDO
!	  if(JL(L).NE.94) then
 !       aa_u=(H2P(L)+H2P(L-1))*0.5
!        aa_v=(H2P(L)+H2P(LSC(L)))*0.5 
!	  if(abs(aa_u-U1V(L)).GT.0.2.or.abs(V1U(L)-aa_v).Gt.0.2 )then
!	  write(*,*)'Check HUV  ',L,IL(L),JL(L)
!	  write(*,*)H2P(L),U1V(L),V1U(L),aa_u,aa_v
!	  endif
!	  endif
!1278    continue	  	  
	  ENDDO

	  DO K=1,KC 
	  DO L=1,LA 
        VHDX2(L,K)=VHDX1(L,K)
	  UHDY2(L,K)=UHDY1(L,K)
	  TEM1(L,K)=TEMLPF(L,K)
	  SEDTLPF(L,K)=SEDTLPF(L,K)
	  SAL1(L,K)=SALLPF(L,K)
	  U2(L,K)=U1(L,K)
	  V2(L,K)=V1(L,K)
	  ENDDO
	  ENDDO

	  DO K=0,KC
	  DO L=1,LA
        W2(L,K)=W1(L,K)
        ENDDO
	  ENDDO

	  DO K=1,KC-1 
	  DO L=1,LA 
	  ABLPF(L,K)=ABEFF(L,K)
	  ENDDO
	  ENDDO

  	  t_1=t_2
  	  
  	  do LLL=1,iskp
        n_re=n_re+1       
        read(251)t_11,N_2,H1P,H1U,H1V   ! Elevtion  	  
        read(252)t_11,N_2,VHDX1,UHDY1   ! Horizontal transport
        read(253)t_12,N_2,U1,V1         ! U,V velocity
        read(254)t_13,N_2,W1            ! W velocity
        read(255)t_2,N_2,ABEFF         ! Edd difusivity      
        if(ISTRAN(2).eq.1)read(256)t_12,N_2,TEMLPF        ! Temperature
        if(ISTRAN(6).eq.1)read(257)t_13,N_2,SNDTLPF       ! Sendiment
        if(ISTRAN(1).eq.1)read(258)t_2,N_2,SALLPF  
        
        IF(IRELSOP.EQ.5)t_2=NTC+30-t_2
                
        enddo
        if(n_re+iskp .eq. 24*30*inc_re)then 
 !        if(mod(N+NTSPTC/(24*inc_re),NTSPTC*30).eq.0)then 
	  !  open new file next timestep
	   do i=251,258 
	    close(i)
	   enddo
	   init=.true.
	   n_re = 0
	   irecd=irecd+1
	   nfin=nfin+1
	   write(*,*)'Open new files ',irecd
	  endif
C	  
	ENDIF

       t_dff=t_2-t_1
	 lmd_1=(TIME-t_1)/t_dff
	 lmd_2=(t_2-TIME)/t_dff
       if(lmd_1.LT.0) then
	  lmd_1=0.0 
	  lmd_2=1
	 endif
       if(lmd_2.LT.0) then
	  lmd_1=1.0 
	  lmd_2=0
	 endif
!
! Check
!	 
        if(mod(N,1440).eq.0)then
 !       if(mod(n_re,24).eq.0) then
        write(*,*)
     $    'Input at t1 t2 t ',t_1,t_2,time,n_re
 !       write(*,*)' Lamd 1, 2 ',lmd_1,lmd_2
 !       write(*,*)'Data input time ',t_2
        endif
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)        
       DO L=1,LA
	  H2WQ(L)=HWQ(L)
        HU(L)=U1V(L)*lmd_2+H1U(L)*lmd_1
        HV(L)=V1U(L)*lmd_2+H1V(L)*lmd_1
        HWQ(L)=H2P(L)*lmd_2+H1P(L)*lmd_1
        HWQ(L)=max(HWQ(L),0.3)
        HPWQ(L)=HWQ(L)
        HP(L)=HWQ(L)
        P(L)=HWQ(L)*9.8
	 ENDDO

	 DO K=1,KC
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
	 DO L=1,LA
        VHDXWQ(L,K)=VHDX2(L,K)*lmd_2+VHDX1(L,K)*lmd_1
        UHDYWQ(L,K)=UHDY2(L,K)*lmd_2+UHDY1(L,K)*lmd_1
        UWQ(L,K)=U2(L,K)*lmd_2+U1(L,K)*lmd_1
	  VWQ(L,K)=V2(L,K)*lmd_2+V1(L,K)*lmd_1
	  UWQS(L)=max(UWQ(L,KC), UWQ(L+1,KC))
	  VWQS(L)=max(VWQ(L,KC), VWQ(LNC(L),KC)) 	  
	  SALWQ(L,K)=SAL1(L,K)*lmd_2+SALLPF(L,K)*lmd_1
	  TEMWQ(L,K)=TEM1(L,K)*lmd_2+TEMLPF(L,K)*lmd_1
	  TEMWQ(L,K)=min(TEMWQ(L,K),40.0)
	  TEMWQ(L,K)=max(TEMWQ(L,K),-4.0)
       ENDDO
	 ENDDO
	 
	  if(ISTRAN(6).eq.0) then
 	  DO K=1,KC
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
	  DO L=2,LA
	   SEDTWQ(L,K)=15
        ENDDO
	  ENDDO
	  else
 	  DO K=1,KC
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
	  DO L=2,LA	  
	   SEDTWQ(L,K)=SEDTLPF(L,K)*lmd_2+SNDTLPF(L,K)*lmd_1	       
        ENDDO
	  ENDDO
        endif

	  
!	 if(idiaout.LE.0)then
	 DO K=1,KC
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
	 DO L=1,LA
	  SAL(L,K)=SALWQ(L,K)
	  TEM(L,K)=TEMWQ(L,K)	
	  U(L,K)=UWQ(L,K)
	  V(L,K)=VWQ(L,K)
	  SEDT(L,K)=SEDTWQ(L,K)
	 ENDDO
	 ENDDO 
!	 endif

	 DO K=0,KC
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
	 DO L=1,LA
        WWQ(L,K)=(W2(L,K)*lmd_2+W1(L,K)*lmd_1) !*FLXPND(L,4)
       ENDDO
	 ENDDO
 
 	 DO K=1,KC-1 
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
	 DO L=1,LA 
	  AB(L,K)=(ABLPF(L,K)*lmd_2+ABEFF(L,K)*lmd_1) !*FLXPND(L,5)
	 ENDDO
	 ENDDO

       if(mod(N,20).EQ.0) THEN
	  DTCFL=1.0e10
	  DO K=2,KS
        DO L=2,LA
        LN=LNC(L)
        UWTMP=ABS(DXIU(L  )*UWQ(L  ,K))
        UETMP=ABS(DXIU(L+1)*UWQ(L+1,K))
        VSTMP=ABS(DYIV(L  )*VWQ(L  ,K))
        VNTMP=ABS(DYIV(LN )*VWQ(LN ,K)) 
        WBTMP=ABS(DZIC(K)*WWQ(L,K-1)/HWQ(L))
        WTTMP=ABS(DZIC(K)*WWQ(L,K  )/HWQ(L))
        DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)
     $         +1.0E-12
        DTMAX=0.5/DTMAXI
        IF (DTMAX.LT.DTCFL) THEN
          DTCFL=DTMAX
          ICFL=IL(L)
          JCFL=JL(L)
          KCFL=K
        END IF
        END DO
        END DO
 !       write(610,'(3I8,2F15.3)')ICFL,JCFL,K,TIME,DTCFL    ! J.S. 6/2016
        ENDIF
!	  if(n.lt.5) write(*,*)w_ts,ab_ts
	  
c	  if(mod(N,NTSPTC*30-2).eq.0)then
c	  !  open new file next timestep
c	  do i=251,258
c	  close(i)
c	  enddo
c	  init=.true.
c	  n_re = nrecs
c	  nfin=nfin+1
c	  write(*,*)'Ready to read new files ',nfin
c	  endif
c
      if((TIME-(TBEGIN)).GE.RCDAY0-0.000001) then
 	  do i=251,258
	  close(i)
	  enddo
	  init=.true.  
	  Dbeg=Dbeg+RCDAY 
	  n_re=nrecs 
	  RCDAY0=RCDAY0+RCDAY
	  
! Reinitical time series

      DO NS=1,NQSER
       MQTLAST(NS)=1
      END DO

      NTOX1=NTOX
      IF(NPCB.GT.0)NTOX1=NPCB
      
      NTMP=4+NSED+NSND+NTOX1
      DO NC=1,NTMP
      DO NN=1,NCSER(NC)
      MCTLAST(NN,NC)=1
      END DO
      END DO
      
      DO NS=1,NWSER
        MWTLAST(NS)=1
      END DO
           
      DO NS=1,NQSER
       DO M=1,MQSER(NS)
        TQSER(M,NS)=TQSER(M,NS)+RCDAY
       END DO
      END DO
      
      NC=MSVTOX(1)
      DO NS=1,NCSER(NC)     
       DO M=1,MCSER(NS,NC)
        TCSER(M,NS,NC)=TCSER(M,NS,NC)+RCDAY     
       ENDDO
      ENDDO
 
      DO N=1,NWSER
       DO M=1,MWSER(N)
       TWSER(M,N)=TWSER(M,N)+RCDAY  
       END DO
      ENDDO
C     
      endif 	  
1000	RETURN
      END  !SUBROUTINE flowin
C==============================================
      SUBROUTINE flow_in_r
!
!  Save 3D flow fileld every FLO_INT
! 
      INCLUDE 'efdc.par'  
      INCLUDE 'efdc.cmn'
      include 'wq.par'
      INCLUDE 'wqcom.cmn'
      INTEGER, SAVE :: n_re = 0,N_1,N_2,nopen=0,nfin=0
	REAL, SAVE :: t_1,t_2,Dbeg,RCDAY0
	REAL :: lmd_1,lmd_2,t_t11,t_12,t_13,t_14
      INTEGER :: oute=11
      LOGICAL :: op
      LOGICAL, SAVE :: init=.true.
      character*4 tday
      
      nrecs=24*TBEGIN*real(inc_re)
      if(mod(int(TBEGIN),30).eq.0)nrecs=0
	TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON

      NM=LCM*KCM*2
	NSD=LCM*KCM*NSCM

      IF(init)THEN
      init=.false.
      write(tday,'(i4.4)')irecd !nfin
      write(*,*)'Read binary flow field files !',tday

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
     &        ACCESS='direct',RECL=4*(2+NSD))

	  INQUIRE(258,opened=op)
      IF(.not.op)OPEN(258,FILE='SAL'//tday//'.outb', 
     &        ACCESS='direct',RECL=4*(2+LCM*KCM))

       if(nopen.eq.0)then
C !  Read 2 record for starting
        nopen=1
        RCDAY0=RCDAY
        Dbeg=0 
        n_re = nrecs +1  ! Save half hour            
C             Start ftom beging, inc_re=0. Start from day 10, 
C        inc_re=1 save hourly, =2 save half hourly, 
C        nrecs record will not read
C
       read(251,REC=n_re)t_1,N_1,H2P,U1V,V1U    ! Elevtion
       read(252,REC=n_re)t_1,N_1,VHDX2,UHDY2    ! Horizontal transport
       read(253,REC=n_re)t_1,N_1,U2,V2          ! U,V velocity
       read(254,REC=n_re)t_1,N_1,W2             ! W velocity
       read(255,REC=n_re)t_1,N_1,ABLPF             ! Edd difusivity
       read(256,REC=n_re)t_1,N_1,TEM1           ! Temperature
       read(257,REC=n_re)t_1,N_1,SED1           ! Sendiment
       read(258,REC=n_re)t_1,N_1,SAL1   
       n_re=n_re+1
       read(251,REC=n_re)t_2,N_2,H1P,H1U,H1V             ! Elevtion
       read(252,REC=n_re)t_2,N_2,VHDX1,UHDY1    ! Horizontal transport
       read(253,REC=n_re)t_2,N_2,U1,V1          ! U,V velocity
       read(254,REC=n_re)t_2,N_2,W1             ! W velocity
       read(255,REC=n_re)t_2,N_2,ABEFF             ! Edd difusivity
       read(256,REC=n_re)t_2,N_2,TEM            ! Temperature
       read(257,REC=n_re)t_2,N_2,SED            ! Sendiment
       read(258,REC=n_re)t_2,N_2,SAL   
       
 !      write(*,*)'model time and hydro time at start ',Time,t_1
 !      write(*,*)'current record and data record not use',n_re,nrecs
 !      write(*,*)'Timestep of first 2 record, & time ', N_1,N_2,t_1,t_2
  
       n_re = nrecs
       t_2=TBEGIN
      endif
! Test transport
!	  DO LL=1,NPCB
!	  DO K=1,KC 
!	  DO L=1,LA 
!	  TOX(L,K,LL)=SAL(L,K)
!	  TOX1(L,K,LL)=SAL(L,K)
!	  ENDDO
!	  ENDDO
!        ENDDO

      END IF

      IF(TIME.GT.t_2) THEN
      
	  DO L=1,LA
        H2P(L)=H1P(L)
        U1V(L)=H1U(L)
	  V1U(L)=H1V(L) 
	  ENDDO

	  DO K=1,KC 
	  DO L=1,LA 
        VHDX2(L,K)=VHDX1(L,K)
	  UHDY2(L,K)=UHDY1(L,K)
!	  ABLPF(L,K)=ABEFF(L,K)
	  TEM1(L,K)=TEM(L,K)
	  SED1(L,K,1)=SED(L,K,1)
	  SAL1(L,K)=SAL(L,K)
	  U2(L,K)=U1(L,K)
	  V2(L,K)=V1(L,K)
	  ENDDO
	  ENDDO

	  DO K=0,KC
	  DO L=1,LA
        W2(L,K)=W1(L,K)
        ENDDO
	  ENDDO

	  DO K=1,KC-1 
	  DO L=1,LA 
	  ABLPF(L,K)=ABEFF(L,K)
	  ENDDO
	  ENDDO

  !      n_re=n_re+inc_re
        n_re=n_re+1
	  t_1=t_2
        read(251,REC=n_re)t_11,N_2,H1P,H1U,H1V   ! Elevtion  	  
        read(252,REC=n_re)t_11,N_2,VHDX1,UHDY1   ! Horizontal transport
        read(253,REC=n_re)t_12,N_2,U1,V1         ! U,V velocity
        read(254,REC=n_re)t_13,N_2,W1            ! W velocity
        read(255,REC=n_re)t_14,N_2,ABEFF         ! Edd difusivity      
        read(256,REC=n_re)t_12,N_2,TEM           ! Temperature
        read(257,REC=n_re)t_13,N_2,SED            ! Sendiment
        read(258,REC=n_re)t_2,N_2,SAL  
        t_2=t_2+Dbeg
        N00=500
 !       write(*,*)'last record n_re, t1, t2 = ',n_re, t_1,t_2,time      
 !       write(*,*)SAL(N00,KC),TEM(N00,KC),H1P(N00),
 !    &           UHDY1(N00,KC),VHDX1(N00,KC)
 !       write(*,*)U1(N00,KC),U1(N00,1),V1(N00,KC),V1(N00,1)
 !       write(*,*)ABEFF(N00,2),W1(N00,2)
        
        if(n_re .eq. 24*30*inc_re)then
	  !  open new file next timestep
	  do i=251,258
	  close(i)
	  enddo
	  init=.true.
	  n_re = 0
	  irecd=irecd+1
	  nfin=nfin+1
	  write(*,*)'Open new files ',nfin
	  endif
	  
	ENDIF
        if(mod(N,1440).eq.0)then
        write(*,*)
     $    'Input data at ',t_2,time,n_re,n_re-nrecs
 !       write(*,*)'Data input time ',t_2
        endif
       t_dff=t_2-t_1
	 lmd_1=(TIME-t_1)/t_dff
	 lmd_2=(t_2-TIME)/t_dff
       if(lmd_1.LT.0) then
	  lmd_1=0.0 
	  lmd_2=1
	 endif

       DO L=1,LA
	  H2WQ(L)=HWQ(L)
        HU(L)=U1V(L)*lmd_2+H1U(L)*lmd_1
        HV(L)=V1U(L)*lmd_2+H1V(L)*lmd_1
        HWQ(L)=H2P(L)*lmd_2+H1P(L)*lmd_1
        HPWQ(L)=HWQ(L)
        HP(L)=HWQ(L)
        P(L)=HWQ(L)*9.8
	 ENDDO

	 DO K=1,KC
	 DO L=1,LA
	  SALWQ(L,K)=SAL1(L,K)*lmd_2+SAL(L,K)*lmd_1
	  TEMWQ(L,K)=TEM1(L,K)*lmd_2+TEM(L,K)*lmd_1
	  TEMWQ(L,K)=min(TEMWQ(L,K),40.0)
	  TEMWQ(L,K)=max(TEMWQ(L,K),-4.0)
	  SEDTWQ(L,K)=SED1(L,K,1)*lmd_2+SED(L,K,1)*lmd_1	 
!       if( idiaout.GT.0) then	 
        VHDXWQ(L,K)=VHDX2(L,K)*lmd_2+VHDX1(L,K)*lmd_1
        UHDYWQ(L,K)=UHDY2(L,K)*lmd_2+UHDY1(L,K)*lmd_1
        UWQ(L,K)=U2(L,K)*lmd_2+U1(L,K)*lmd_1
	  VWQ(L,K)=V2(L,K)*lmd_2+V1(L,K)*lmd_1
	  UWQS(L)=max(UWQ(L,KC), UWQ(L+1,KC))
	  VWQS(L)=max(VWQ(L,KC), VWQ(LNC(L),KC)) 
	  SEDT(L,K)=SEDTWQ(L,K)	
	  TEM(L,K)=TEMWQ(L,K)	    	  
!	 endif
       ENDDO
	 ENDDO
	 
!	 if(idiaout.eq.2)then
!	 DO K=1,KC
!	 DO L=1,LA
!	  SAL(L,K)=SALWQ(L,K)
!	  TEM(L,K)=TEMWQ(L,K)	
!	  SED1(L,K,1)= SEDTWQ(L,K)
!	 ENDDO
!	 ENDDO 
!	 endif

	 DO K=0,KC
	 DO L=1,LA
	  if(HWQ(L).GE.5) then
        WWQ(L,K)=(W2(L,K)*lmd_2+W1(L,K)*lmd_1) !*FLXPND(L,4)
        else
        WWQ(L,K)=(W2(L,K)*lmd_2+W1(L,K)*lmd_1) !*FLXPND(L,4)       
        endif
       ENDDO
	 ENDDO
 
 	 DO K=1,KC-1 
	 DO L=1,LA 
	  AB(L,K)=(ABLPF(L,K)*lmd_2+ABEFF(L,K)*lmd_1)*FLXPND(L,5)
	 ENDDO
	 ENDDO
	  if(n.lt.5) write(*,*)w_ts,ab_ts
	  
        return
        
c	  if(mod(N,NTSPTC*30-2).eq.0)then
c	  !  open new file next timestep
c	  do i=251,258
c	  close(i)
c	  enddo
c	  init=.true.
c	  n_re = nrecs
c	  nfin=nfin+1
c	  write(*,*)'Ready to read new files ',nfin
c	  endif
c
      if((TIME-(TBEGIN)).GE.RCDAY0-0.000001) then
 	  do i=251,258
	  close(i)
	  enddo
	  init=.true.  
	  Dbeg=Dbeg+RCDAY 
	  n_re=nrecs 
	  RCDAY0=RCDAY0+RCDAY
	  
! Reinitical time series

      DO NS=1,NQSER
       MQTLAST(NS)=1
      END DO

      NTOX1=NTOX
      IF(NPCB.GT.0)NTOX1=NPCB
      
      NTMP=4+NSED+NSND+NTOX1
      DO NC=1,NTMP
      DO NN=1,NCSER(NC)
      MCTLAST(NN,NC)=1
      END DO
      END DO
      
      DO NS=1,NWSER
        MWTLAST(NS)=1
      END DO
           
      DO NS=1,NQSER
       DO M=1,MQSER(NS)
        TQSER(M,NS)=TQSER(M,NS)+RCDAY
       END DO
      END DO
      
      NC=MSVTOX(1)
      DO NS=1,NCSER(NC)     
       DO M=1,MCSER(NS,NC)
        TCSER(M,NS,NC)=TCSER(M,NS,NC)+RCDAY     
       ENDDO
      ENDDO
 
      DO N=1,NWSER
       DO M=1,MWSER(N)
       TWSER(M,N)=TWSER(M,N)+RCDAY  
       END DO
      ENDDO
C     
      endif 	  
1000	RETURN
      END  !SUBROUTINE flowin