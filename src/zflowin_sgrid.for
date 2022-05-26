      SUBROUTINE flow_in_sub
!
!  Save 3D flow fileld every FLO_INT
!  Last Version 2/24/2014 by J.S.
!  Use Jm_make_sbugrid_lij_mapping to make new newlijmap.inp mapping file
!
      INCLUDE 'efdc.par'  
      INCLUDE 'efdc.cmn'
      include 'wq.par'
      INCLUDE 'wqcom.cmn'
      INTEGER, SAVE :: n_re = 0,N_1,N_2,nopen=0,nfin=0,LAOLD,LANEW
	REAL, SAVE :: t_1,t_2,Dbeg,RCDAY0
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
      open(1,file='newlijmap.inp')
      read(1,*)LANEW,LAOLD
      do i=2,LANEW
      read(1,*)NSUBLMP(i,1),NSUBLMP(i,2),NSUBLMP(i,3),
     &     NSUBLMP(i,4)
      enddo
      NSUBLMP(LANEW+1,1)=i
      NSUBLMP(LANEW+1,2)=1
      NSUBLMP(1,2)=1
      close(1)
      
      write(tday,'(i4.4)')irecd !nfin
      write(*,*)'Read binary flow field files !',tday

      INQUIRE(251,opened=op)
      IF(.not.op)OPEN(251,FILE=trim(path)//'HWQ'//tday//'.outb',
     &       	STATUS='UNKNOWN', form='unformatted')

      INQUIRE(252,opened=op)
      IF(.not.op)OPEN(252,FILE=trim(path)//'UHDXYWQ'//tday//'.outb', 
     &          	STATUS='UNKNOWN', form='unformatted')
      
	  INQUIRE(253,opened=op)
      IF(.not.op)OPEN(253,FILE=trim(path)//'UVWQ'//tday//'.outb', 
     &       STATUS='UNKNOWN', form='unformatted')

	  INQUIRE(254,opened=op)
      IF(.not.op)OPEN(254,FILE=trim(path)//'WWQ'//tday//'.outb', 
     &    	STATUS='UNKNOWN', form='unformatted')

	  INQUIRE(255,opened=op)
      IF(.not.op)OPEN(255,FILE=trim(path)//'AB'//tday//'.outb',
     &       	STATUS='UNKNOWN', form='unformatted')
 
      if(ISTRAN(2).eq.1)then
	  INQUIRE(256,opened=op)
      IF(.not.op)OPEN(256,FILE=trim(path1)//'TEM'//tday//'.outb' 
     &      	,STATUS='UNKNOWN', form='unformatted')
      endif
      
      if(ISTRAN(6).eq.1)then
	  INQUIRE(257,opened=op)
      IF(.not.op)OPEN(257,FILE=trim(path1)//'SED'//tday//'.outb' 
     &       	,STATUS='UNKNOWN', form='unformatted')
      endif
      
      if(ISTRAN(1).eq.1)then
	  INQUIRE(258,opened=op)
      IF(.not.op)OPEN(258,FILE=trim(path)//'SAL'//tday//'.outb', 
     &      	STATUS='UNKNOWN', form='unformatted')
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
       
       IF(IRELSOP.EQ.5)t_1=NTC+30-t_1
       
111    n_re=n_re+1
       read(251)t_22,N_2,H1P,H1U,H1V             ! Elevtion
       read(252)t_2,N_2,VHDX1,UHDY1    ! Horizontal transport
       read(253)t_2,N_2,U1,V1          ! U,V velocity
       read(254)t_2,N_2,W1             ! W velocity
       read(255)t_2,N_2,ABEFF          ! Edd difusivity
       if(ISTRAN(2).eq.1)read(256)t_2,N_2,TEMLPF         ! Temperature
       if(ISTRAN(6).eq.1)read(257)t_2,N_2,SNDTLPF        ! Sendiment
       if(ISTRAN(1).eq.1)read(258)t_2,N_2,SALLPF   
       
       IF(IRELSOP.EQ.5)t_2=NTC+30-t_2
       write(*,*)'First read new file ',tday,TIME,t_2,t_22

       DO L=1,LAOLD+1
        HWQ(L)=H2P(L)       
       ENDDO
       
       if(TIME.GT.t_2) THEN   ! Take care restart to fine current recored. 
        write(*,*)'Start hydro enter ',TIME, t_2,n_re
       	DO L=1,LAOLD+1
        H2P(L)=H1P(L)
        U1V(L)=H1U(L)
	  V1U(L)=H1V(L) 
	  HWQ(L)=H1P(L)
	  ENDDO

	  DO K=1,KC 
	  DO L=1,LAOLD+1 
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
	  DO L=1,LAOLD+2
        W2(L,K)=W1(L,K)
        ENDDO
	  ENDDO

	  DO K=1,KC-1 
	  DO L=1,LAOLD+1 
	  ABLPF(L,K)=ABEFF(L,K)
	  ENDDO
	  ENDDO
  	  t_1=t_2
  	  goto 111  	  
  	 endif
      endif

      END IF ! end IF(init)
 !     write(*,*)TIME,t_1,t_2
      IF(TIME.GT.t_2) THEN
      
	  DO L=1,LAOLD+1
        H2P(L)=H1P(L)
        U1V(L)=H1U(L)
	  V1U(L)=H1V(L) 
	  ENDDO

	  DO K=1,KC 
	  DO L=1,LAOLD+1 
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
	  DO L=1,LAOLD+1
        W2(L,K)=W1(L,K)
        ENDDO
	  ENDDO

	  DO K=1,KC-1 
	  DO L=1,LAOLD+2 
	  ABLPF(L,K)=ABEFF(L,K)
	  ENDDO
	  ENDDO

  	  t_1=t_2
  	  
  	  do LLL=1,iskp
        n_re=n_re+1       
        read(251)t_11,N_2,H1P,H1U,H1V   ! Elevtion  	  
        read(252)t_11,N_2,VHDX1,UHDY1   ! Horizontal transport
        read(253)t_12,N_2,U1,V1         ! U,V velocity
        read(254)t_2,N_2,W1            ! W velocity
        read(255)t_2,N_2,ABEFF         ! Edd difusivity      
        if(ISTRAN(2).eq.1)read(256)t_12,N_2,TEMLPF        ! Temperature
        if(ISTRAN(6).eq.1)read(257)t_13,N_2,SNDTLPF       ! Sendiment
        if(ISTRAN(1).eq.1)read(258)t_2,N_2,SALLPF  
        
        IF(IRELSOP.EQ.5)t_2=NTC+30-t_2
 !        write(*,*)'Update t1 t2',time,t_1,t_2     
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
        write(*,*)
     $    'Input at t1 t2 t ',t_1,t_2,time,n_re
 !       write(*,*)' Lamd 1, 2 ',lmd_1,lmd_2
 !       write(*,*)'Data input time ',t_2
        endif
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)        
       DO LL=2,LA
        L=NSUBLMP(LL,2)
	  H2WQ(LL)=HWQ(LL)
        HU(LL)=U1V(L)*lmd_2+H1U(L)*lmd_1
        HV(LL)=V1U(L)*lmd_2+H1V(L)*lmd_1
        HWQ(LL)=H2P(L)*lmd_2+H1P(L)*lmd_1
        HWQ(LL)=max(HWQ(LL),0.5)
        HPWQ(LL)=HWQ(LL)
        HP(LL)=HWQ(LL)
        P(LL)=HWQ(LL)*9.8
	 ENDDO

	 DO K=1,KC
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
	 DO LL=2,LA
	  L=NSUBLMP(LL,2)
        VHDXWQ(LL,K)=VHDX2(L,K)*lmd_2+VHDX1(L,K)*lmd_1
        UHDYWQ(LL,K)=UHDY2(L,K)*lmd_2+UHDY1(L,K)*lmd_1
        UWQ(LL,K)=U2(L,K)*lmd_2+U1(L,K)*lmd_1
	  VWQ(LL,K)=V2(L,K)*lmd_2+V1(L,K)*lmd_1
	  UWQS(LL)=max(UWQ(LL,KC), UWQ(LL+1,KC))
	  VWQS(LL)=max(VWQ(LL,KC), VWQ(LNC(LL),KC)) 	  
	  SALWQ(LL,K)=SAL1(L,K)*lmd_2+SALLPF(L,K)*lmd_1
	  TEMWQ(LL,K)=TEM1(L,K)*lmd_2+TEMLPF(L,K)*lmd_1
	  TEMWQ(LL,K)=min(TEMWQ(LL,K),40.0)
	  TEMWQ(LL,K)=max(TEMWQ(LL,K),-4.0)
       ENDDO
	 ENDDO
	 
	  if(ISTRAN(6).eq.0) then
 	  DO K=1,KC
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
	  DO LL=2,LA
	   L=NSUBLMP(LL,2)
	   SEDTWQ(LL,K)=15
        ENDDO
	  ENDDO
	  else
 	  DO K=1,KC
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
	  DO LL=2,LA
	 	 L=NSUBLMP(LL,2) 	  
	   SEDTWQ(LL,K)=SEDTLPF(L,K)*lmd_2+SNDTLPF(L,K)*lmd_1	       
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
	 DO LL=2,LA
	 	 L=NSUBLMP(LL,2)
        WWQ(LL,K)=(W2(L,K)*lmd_2+W1(L,K)*lmd_1) !*FLXPND(L,4)
       ENDDO
	 ENDDO
 
 	 DO K=1,KC-1 
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
	 DO LL=2,LA 
	 	L=NSUBLMP(LL,2)
	  AB(LL,K)=(ABLPF(L,K)*lmd_2+ABEFF(L,K)*lmd_1) !*FLXPND(L,5)
	 ENDDO
	 ENDDO
	 
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
