c inputw.for, read efdcwin.inp
c      program inputw
      SUBROUTINE INPUTw
C 
C**********************************************************************C
C
C **  SUBROUTINE INPUTw READS efdcwin.inp
c
c     jeff Ji, 10/16/987
C
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
      parameter (ncread=120,ntides=MTM) !
      dimension TCP9(ntides)
      character aaa(ncread),xsymbol*5,symbol9(ntides)*5
C
C**********************************************************************C
C
      G=9.81
      PI=3.14159265358987
      PI2=2.*PI
c tidal periods, James Martin, 1999, P546 & \toolkit\...\tides.gif
      symbol9(1)='M2'
      TCP9(1)   =12.42*3600
      symbol9(2)='S2'
      TCP9(2)   =12.00*3600
      symbol9(3)='N2'
      TCP9(3)   =12.66*3600
      symbol9(4)='K2'
      TCP9(4)   =11.97*3600
      symbol9(5)='K1'
      TCP9(5)   =23.93*3600
      symbol9(6)='O1'
      TCP9(6)   =25.82*3600
      symbol9(7)='P1'
      TCP9(7)   =24.07*3600
      symbol9(8)='Q1'
      TCP9(8)   =26.87*3600
      symbol9(9)='MF'
      TCP9(9)   =327.9*3600
      symbol9(10)='MM'
      TCP9(10)   =661.30*3600
c     symbol9(11)='SSA'
c     TCP9(11)   =4383.30*3600
c     symbol9(12)='M4'
c     TCP9(12)   =6.21*3600
c     symbol9(13)='MS4'
c     TCP9(13)   =6.10*3600

c
C**********************************************************************C
C**********************************************************************C
C
C --- READ MAIN INPUT FILE efdcwin.inp  ------------------------------------
       nefdcwin=211
       open(nefdcwin,file='efdcwin.inp',status='old',err=977)
c
       nefdcout=221
       open(nefdcout,file='efdcwin.out')
c      go to 2000
c-----------------------------------------------
C1** Submodel swiches
C
c
       ncard=1
       call PasteC(nefdcwin,nefdcout)
c
       read(nefdcwin,*,IOSTAT=ISO)
     * ISTRAN(1),ISTRAN(2),ISTRAN(3),ISTRAN(8),
c        sal       tem       dye      wq
c    * ISTRAN(6),ISTRAN(7),ISTRAN(5)
     * NSED9,      NSND9,     NTOX9,IWQS
c        sed       snd        tx
       if(ISO.NE.0) go to 987
       write(nefdcout,7201)
     * ISTRAN(1),ISTRAN(2),ISTRAN(3),ISTRAN(8),
c    * ISTRAN(6),ISTRAN(7),ISTRAN(5)
     * NSED9,      NSND9,     NTOX9
 7201  format(999i8)
c
      ISTRAN(6)=0
      ISTRAN(7)=0
      ISTRAN(5)=0
      if(NSED9.gt.0) ISTRAN(6)=1
      if(NSND9.gt.0) ISTRAN(7)=1
      if(NTOX9.gt.0) ISTRAN(5)=1
c
      NSED=MAX(1,NSED9)
      NSND=MAX(1,NSND9)
      NTOX=MAX(1,NTOX9)
c
      MTMP=4             ! efdc.inp, Card 22
      DO N=1,NTOX
       MTMP=MTMP+1
       MSVTOX(N)=MTMP
      END DO
      DO N=1,NSED
       MTMP=MTMP+1
       MSVSED(N)=MTMP
      END DO
      DO N=1,NSND
       MTMP=MTMP+1
       MSVSND(N)=MTMP
      END DO
C
      NCSER(4)=0    ! shell fish
c     NQSER=0
c     ncser(1)=0              ! serach for real max later
c     ncser(2)=0
c     ncser(3)=0
c     NTOXSER=0
c     NSEDSER=0
c     NSNDSER=0
c     NPSER=0
c
c     do i=0,8
c     write(6,*) "istran ", i,istran(i)   ! debug
c     enddo
C
C2** Major Parameters
C
       ncard=2
       call PasteC(nefdcwin,nefdcout)
c
       read(nefdcwin,*,IOSTAT=ISO) isresti,NTC,NTSPTC,TBEGIN,ZBRADJ
     & ,Alati,TEMO,RKDYE
       if(ISO.NE.0) go to 987
       write(nefdcout,7202) isresti,NTC,NTSPTC,TBEGIN,ZBRADJ
     & ,Alati,TEMO,RKDYE
 7202  format(3i8,4f11.3,E10.3)
c check for SWAN, JI@WAVE, 3/5/00
	if((NSED9+NSND9).eq.0) ISWAN=0   ! when no SED => SWAN not needed
      if(ISWAN.gt.0) then
      if(dayswan1.gt.TBEGIN) then
      write(6,*) "hsdir.inp error: dayswan1,tbegin=",dayswan1,tbegin
      stop
      endif
      tend=tbegin+NTC
      if(dayswan2.lt.tend) then
      write(6,*) "hsdir.inp error: dayswan2,tend=",dayswan2,tend
      stop
      endif
      endif
c-
c check
      if(NTSPTC*2.gt.NTSM   ) then   !
      write(6,*) " Error:  NTSPTC*2>  NTSM   ", NTSPTC, NTSM
      stop
      endif
c
c       if(isresti.gt.0) then
       do i=1,8
       isci(i)=0
       if(ISTRAN(I).ge.1) isci(i)=1
       isco(i)=isci(i)
c	write(501,*) i,isci(i),isco(i),istran(i)
       enddo
c       endif
       CF=2*7.29e-5*sin(3.14*alati/180.0)
C
c      go to 2000
C3** Layer Thickness in vertical
 8003 continue
C
       ncard=3
       call PasteC(nefdcwin,nefdcout)
c
       do ncountx=1,99999
       read(nefdcwin,*,IOSTAT=ISO) KDUM,DZC(KDUM)
       if(ISO.NE.0) exit
       write(nefdcout,7203) KDUM,DZC(KDUM)
 7203  format(i10,f16.2)
       enddo
c
       backspace (nefdcwin)
       read(nefdcwin,997) aaa
       if(aaa(1).ne.'C'.and.aaa(1).ne.'c') go to 987
       write(nefdcout,997) aaa
c
c      go to 2000
       KC=ncountx-1
       write(6,*) "KC=  ",KC
       if(kc.gt.kcm) then                    ! ensure KC <= KCM
       write(6,*) "Error in kc,kcm: ", kc,kcm
       go to 987
       endif
       xxsum=0.0                             ! ensure sum of DZC = 1.0
       do k=1,kc
       xxsum=xxsum+dzc(k)
       enddo
       if(abs(xxsum-1.0).ge.1.0e-5) then
       write(6,*) "Error in KC or DZC(k):  ",xxsum
       go to 987
       endif
C
c      go to 2000
c     IF (KC.GE.2) ISITB=0           ! YES! impact efdcwin.inp logic, ji, 10/16/987
      IF (KC.GE.2.AND.ISVEG.EQ.0) ISITB=0  !Ji, The test case, 1/19/00,
c      go to 2000
C4 Harmonic Tidal FORCING SURF ELEV OR PRESSURE BOUNDARY COND. FORCINGS
C
       ncard=4
       call PasteC(nefdcwin,nefdcout)
c
       M=0
       NP=0
       do N=1,99999
c      read(nefdcwin,*,IOSTAT=ISO) NDUM,xsymbol,PFAM(NP,M),PFPH(NP,M)
       read(nefdcwin,*,IOSTAT=ISO) NP,xsymbol,xPFAM,xPFPH
       if(ISO.NE.0) exit
       if(N.eq.1) NP2=NP  ! first column
c match tidal components
       if(NP.EQ.NP2) then
       M=M+1
       ELSE
       M=1
       ENDIF
       NP2=NP
c
       PFAM(NP,M)=xpfam
       PFPH(NP,M)=xpfph
       SYMBOL(M)=xsymbol
c
       write(nefdcout,7205) NP,xsymbol,PFAM(NP,M),PFPH(NP,M)
 7205  format(i6,4x,a5,2e14.5)
       enddo
c
       MTIDE=M
       write(6,*) "MTIDE=  ",MTIDE
       if(MTIDE.gt.MTM) then
       write(6,*) "Error in MTIDE,MTM: ",MTIDE,MTM
       go to 987
       endif
c
       npfor=np
       write(6,*) "NPFOR=  ",NPFOR
       if(NPFOR.GT.NPFORM) THEN
       write(6,*) "Error in NPFOR,NPFORM: ",NPFOR,NPFORM
       go to 987
       endif
c match tidal periods
       do M=1,MTIDE
       nxx=0
       do i=1,ntides
       if(symbol(m).eq.symbol9(i)) then
       TCP(m)=TCP9(i)
       write(nefdcout,7904) M, symbol(m),tcp(m)
7904   format(i6,2x,a5,f14.2)
       nxx=1
       exit
       endif 
       enddo
       if(nxx.eq.0) then      ! tidal period is not found!
       write(6,*) "Error: Tidal period is not found at Card ", ncard
       stop
       endif
       enddo
c
       backspace (nefdcwin)
       read(nefdcwin,997) aaa
       if(aaa(1).ne.'C'.and.aaa(1).ne.'c') go to 987
c      if(m.ne.1) then
c      write(6,*) " m: ", m
c      go to 987        ! ensure record is complete, and no incomplete record
c      endif
       write(nefdcout,997) aaa
c
c      go to 2000
C5** Point source discharge locations, concentrations, and loadings
C
       ncard=5
       call PasteC(nefdcwin,nefdcout)
c zeroes to remove values from efdc.inp, added on 10/26/02
      NQSER=0
      NCSER(1)=0
      NCSER(2)=0
      NCSER(3)=0
      NTOXSER=0
      NSEDSER=0
      NSNDSER=0
c
	npstmsr2=0
       do L=1,99999
       read(nefdcwin,*,IOSTAT=ISO)
     $  IQS(L),JQS(L),QSSE,
     $          NQSERQ(L),NCSERQ(L,1),NCSERQ(L,2),NCSERQ(L,3),
     $                      NTOXSRQ,NSEDSRQ,NSNDSRQ,ITMP
       if(ISO.NE.0) exit
       write(nefdcout,7271)
     $  IQS(L),JQS(L),QSSE,
     $          NQSERQ(L),NCSERQ(L,1),NCSERQ(L,2),NCSERQ(L,3),
     $                      NTOXSRQ,NSEDSRQ,NSNDSRQ,ITMP
7271  format(2i8,f12.3,999i5)
c
      if(ISTRAN(8).gt.0)  then
!!      MVPSL2(L)=ITMP
      NPSTMSR2=MAX(NPSTMSR2,ITMP)
      endif
c search
      NQSER=MAX(NQSER,NQSERQ(L))
      NCSER(1)=MAX(NCSER(1),NCSERQ(L,1))
      NCSER(2)=MAX(NCSER(2),NCSERQ(L,2))
      NCSER(3)=MAX(NCSER(3),NCSERQ(L,3))
      NTOXSER=MAX(NTOXSER,NTOXSRQ)
      NSEDSER=MAX(NSEDSER,NSEDSRQ)
      NSNDSER=MAX(NSNDSER,NSNDSRQ)
c
        DO K=1,KC
        QSS(K,L)=QSSE*DZC(K)
        END DO
       DO N=1,NTOX
        M=MSVTOX(N)
        NCSERQ(L,M)=NTOXSRQ
       END DO
       DO N=1,NSED
        M=MSVSED(N)
        NCSERQ(L,M)=NSEDSRQ
       END DO
       DO N=1,NSND
        M=MSVSND(N)
        NCSERQ(L,M)=NSNDSRQ
       END DO
       enddo
c
       backspace (nefdcwin)
       read(nefdcwin,997) aaa
       if(aaa(1).ne.'C'.and.aaa(1).ne.'c') go to 987
       write(nefdcout,997) aaa
c
       NQJPIJ=0           ! Jet Plume
       NQSIJ=L-1
       write(6,*) "NQSij=  ",NQSij
       if(NQSij.GT.NQSijM) THEN
       write(6,*) "Error in NQSij,NQSijM: ",NQSij,NQSijM
       go to 987
       endif
c
      MMAX=4+NTOX                     ! C25*  READ TIME CONSTANT VOLUMETRIC SOURCE INFLOW CONCENTRAT
      DO L=1,NQSIJ
      do m=1,mmax
      CQSE(M)=0.0
      enddo
       DO MS=1,MMAX
        DO K=1,KC
        CQS(K,L,MS)=CQSE(MS)
        END DO
       END DO
      END DO
C
      MMIN=MMAX+1                      ! C26*  READ TIME CONSTANT VOLUMETRIC SOURCE INFLOW CONCENTRA
      MMAX=MMAX+NSED+NSND              ! C     SED(1 TO NSED),SND(1 TO NSND)
      DO L=1,NQSIJ
      do M=mmin,mmax
      CQSE(M)=0.0
      enddo
       DO MS=MMIN,MMAX
        DO K=1,KC
        CQS(K,L,MS)=CQSE(MS)
        END DO
       END DO
      END DO
C
C6**  READ SURF ELEV OR PRESS DEPENDENT FLOW CONTROL STRUCTURE INFO
C
      NCARD=6
       call PasteC(nefdcwin,nefdcout)
c
       do L=1,99999
      READ (nefdcwin,*,IOSTAT=ISO)IQCTLU(L),JQCTLU(L),IQCTLD(L),
     &                  JQCTLD(L), NQCTYP(L),NQCTLQ(L)
       if(ISO.NE.0) exit
      write(nefdcout,7201)        IQCTLU(L),JQCTLU(L),IQCTLD(L),
     &                  JQCTLD(L),NQCTYP(L),NQCTLQ(L)
c search
      NQCTLT=MAX(NQCTLT,NQCTLQ(L))
       enddo
c
       backspace (nefdcwin)
       read(nefdcwin,997) aaa
       if(aaa(1).ne.'C'.and.aaa(1).ne.'c') go to 987
       write(nefdcout,997) aaa
c
       NQctl=L-1
       write(6,*) "NQctl=  ",NQctl
       if(NQctl.GT.NQctlM) THEN
       write(6,*) "Error in NQctl,NQctlM: ",NQctl,NQctlM
       go to 987
       endif
c
      DO L=1,NQCTL
      NQCMUL(L)=0
      NQCMFU(L)=0
      NQCMFD(L)=0
      IQCAX(L)=0
      JQCAX(L)=0
       DO K=1,KC
         QCTLTO(K,L)=0.
         QCTLT(K,L)=0.
       END DO
      END DO
C
C7**  READ FLOW WITHDRAWAL, HEAT OR MATERIAL ADDITION, FLOW RETURN DATA
C
      NCARD=7
       call PasteC(nefdcwin,nefdcout)
c
      MMAX=4+NTOX
      MMIN=MMAX+1
      MMAX2=MMAX+NSED+NSND
      DO L=1,99999
      READ (nefdcwin,*,IOSTAT=ISO)IQWRU(L),JQWRU(L),KQWRU(L),
     $                     IQWRD(L),JQWRD(L),KQWRD(L),
     &          QWR(L),NQWRSERQ(L),
     * (CQWR(L,MS),MS=1,2),(CQWR(L,MS),MS=4,MMAX) ! skip CQWR(L,3)=SFL
     *,(CQWR(L,MS),MS=MMIN,MMAX2)  ! SED & SND
       if(ISO.NE.0) exit
      write(nefdcout,7401)        IQWRU(L),JQWRU(L),KQWRU(L),
     $                     IQWRD(L),JQWRD(L),KQWRD(L),
     &          QWR(L),NQWRSERQ(L),
     * (CQWR(L,MS),MS=1,2),(CQWR(L,MS),MS=4,MMAX) ! skip CQWR(L,3)=SFL
     *,(CQWR(L,MS),MS=MMIN,MMAX2)  ! SED & SND
7401  format(6I5,f7.2,I5,999f7.2)
       CQWR(L,3)=0.0                              ! SFL
c search
      NQWRSR=MAX(NQWRSR,NQWRSERQ(L))
       enddo
c
       backspace (nefdcwin)
       read(nefdcwin,997) aaa
       if(aaa(1).ne.'C'.and.aaa(1).ne.'c') go to 987
       write(nefdcout,997) aaa
c
       NQwr =L-1
       write(6,*) "NQwr =  ",NQwr
       if(NQwr .GT.NQwrm ) THEN
       write(6,*) "Error in NQwr, NQwrm : ",NQwr ,NQwrm
       go to 987
       endif
C
C8** COHESIVE SUSPENDED SEDIMENT SOURCE/SINK PARAMETERS REPEAT DATA LINE       |
C
       ncard=8
       call PasteC(nefdcwin,nefdcout)
c
      NSEDTMP=MAX(NSED,1)
c     DO N=1,NSEDTMP
      DO N=1,NSED9
       READ (nefdcwin,*,IOSTAT=ISO) ISEDINT,SEDO(N),SEDBO(N),
     $                  WSEDO(N),WRSPO(N),TAUR(N),RTAUD
      if(ISO.NE.0) go to 987
      write(nefdcout,9667)  ISEDINT,SEDO(N),SEDBO(N),
     $                  WSEDO(N),WRSPO(N),TAUR(N),RTAUD
9667  format(i10,999e9.2)
c
	TAUD(N)=RTAUD*TAUR(N)
c
      TAUN(N)=TAUR(N)           ! = TAUR(N)
      if(n.gt.1) then           ! bug !!, ji, 10/27/98
      SDEN(N)=SDEN(1)           ! use efdc.inp default value
      SSG(N)=SSG(1)
      SEDN(N)=SEDN(1)
      SEXP(N)=SEXP(1)
      TEXP(N)=TEXP(1)
      SDBLV(N)=SDBLV(1)
      endif
      enddo
c!@ read spacially varing sediment property, jeff ji, 2/7/00
      xWSEDo(1)=wsedo(1)
      xtaud(1)=taud(1)
      xwrspo(1)=wrspo(1)
      xtaur(1)=taur(1)
      do k=1,kbm
      taurc(k)=taur(1)              ! only ONE class of sediment is considered !!
      enddo
c
      if(nsed9.gt.0) then    ! #9
c     do N=2,4                        ! hardwired for test, Ji
c     READ (nefdcwin,*,IOSTAT=ISO)xWSEDO(N),xTAUD(N),xWRSPO(N),xTAUR(N)
c     if(ISO.NE.0) go to 987
c     write(nefdcout,9687) xWSEDO(N),xTAUD(N),xWRSPO(N),xTAUR(N)
c9687  format(28x,999e9.2)
c     enddo
c   spread to entire domain in input.for, after cellmap is called and LIJ(I,J) is ready  !@
c
c bed critical resuspension stress, 10/21/00
c
       do ncountx=1,99999
       read(nefdcwin,*,IOSTAT=ISO) KDUM,TAURC(KDUM)
       if(ISO.NE.0) exit
       write(nefdcout,7273) KDUM,TAURC(KDUM)
7273  format(38x,i6,999e12.4)
       enddo
       KB=ncountx
	write(6,*) "KB  =", kb
       kb=max(kb,1)
       if(kb.ne.kbm.and.kb.gt.1) then
       write(6,*) "?? KB.ne.KBM  ", kb,kbm
       stop
       endif
       backspace (nefdcwin)
c
      endif                ! #9
c
C
C9*** NONCOHESIVE SUSPENDED SEDIMENT SOURCE/SINK PARAMETERS REPEAT DATA LINE    |
C
       ncard=9
       call PasteC(nefdcwin,nefdcout)
c
      NSNDTMP=MAX(NSND,1)
c     write(6,*) "NSEDTMP,NSNDTMP", NSEDTMP,NSNDTMP
c     DO NX=1,NSNDTMP
      DO NX=1,NSND9
       N=NX+NSEDTMP
       READ (nefdcwin,*,IOSTAT=ISO) SEDO(N),SEDBO(N),
     $                  WSEDO(N),SEDN(N),SEXP(N),TAUD(N),WRSPO(N),
     $                  TAUR(N),        TEXP(N)
      if(ISO.NE.0) go to 987
      write(nefdcout,967)  SEDO(N),SEDBO(N),
     $                  WSEDO(N),SEDN(N),SEXP(N),TAUD(N),WRSPO(N),
     $                  TAUR(N),        TEXP(N)
      TAUN(N)=TAUR(N)            ! = TAUR(N)
      if(nx.gt.1) then           ! bug!!, ji, 10/27/98
      N1=1+NSEDTMP               ! take NX=1
      SDEN(N)=SDEN(N1)           ! use efdc.inp default value
      SSG(N)=SSG(N1)
      SDBLV(N)=SDBLV(N1)
      endif
      enddo
C
C10**TOXIC CONTAMINANT INITIAL CONDITIONS AND PARAMETERS                       |
C
       ncard=10
       call PasteC(nefdcwin,nefdcout)
c
      NTOXT=MAX(NTOX,1)
c     DO NT=1,NTOXT
      DO NT=1,NTOX9
        READ (nefdcwin,*,IOSTAT=ISO)NDUM,ITXINT(NT),ITXBDUT(NT)
     $     ,TOXINTW(NT), TOXINTB(NT),DIFTOX(NT)
      if(ISO.NE.0) go to 987
       write (nefdcout,7211)        NDUM,ITXINT(NT),ITXBDUT(NT)
     $     ,TOXINTW(NT), TOXINTB(NT),DIFTOX(NT)
 7211 format(3i6,3E12.4)
      if(NT.gt.1) then
      RKTOXW(NT)=RKTOXW(1)
      TKTOXW(NT)=TKTOXW(1)
      RKTOXB(NT)=RKTOXB(1)
      TRTOXB(NT)=TRTOXB(1)
c
      VOLTOX(NT)=VOLTOX(1)
      RMOLTX(NT)=RMOLTX(1)
      RKTOXP(NT)=RKTOXP(1)
      SKTOXP(NT)=SKTOXP(1)
      endif
      END DO
C
C11  TOXIC CONTAMINANT SEDIMENT INTERACTION PARAMETERS                         |
C
       ncard=11
       call PasteC(nefdcwin,nefdcout)
c
      NSEDTMP=MAX(NSED,1)
      NSNDTMP=MAX(NSND,1)
      NTOXT=MAX(NTOX,1)
c     write(6,*) "NTOXT,NSEDTMP",NTOXT,NSEDTMP
c     DO NT=1,NTOXT
c     DO N=1,NSEDTMP
      DO NT=1,NTOX9
      DO N=1,NSED9
        READ (nefdcwin,*,IOSTAT=ISO)
     $         ITXPARW(N,NT),TOXPARW(N,NT),CONPARW(N,NT),
     $         ITXPARB(N,NT),TOXPARB(N,NT),CONPARB(N,NT)
      if(ISO.NE.0) go to 987
        write(nefdcout,7212)
     $         ITXPARW(N,NT),TOXPARW(N,NT),CONPARW(N,NT),
     $         ITXPARB(N,NT),TOXPARB(N,NT),CONPARB(N,NT)
 7212 format( i4,2f9.4,i6,2f9.4)
      END DO
c     DO NX=1,NSNDTMP
      DO NX=1,NSND9
        N=NX+NSEDTMP
        READ (nefdcwin,*,IOSTAT=ISO)
     $         ITXPARW(N,NT),TOXPARW(N,NT),CONPARW(N,NT),
     $         ITXPARB(N,NT),TOXPARB(N,NT),CONPARB(N,NT)
      if(ISO.NE.0) go to 987
        write(nefdcout,7212)
     $         ITXPARW(N,NT),TOXPARW(N,NT),CONPARW(N,NT),
     $         ITXPARB(N,NT),TOXPARB(N,NT),CONPARB(N,NT)
      END DO
      END DO
C
C12 LOCATION OF PERIODIC FORCING (TIDAL) SURF ELEV OR PRESSURE, and CONC BC'S ON
C     SOUTH BOUNDARIES
C                                         ! be careful! C13-C16 , ji, 10/16/987
       ncard=12
       call PasteC(nefdcwin,nefdcout)
	NPSER=0
c
      DO L=1,99999
c     write(99,*) L
      READ (nefdcwin,*,IOSTAT=ISO) IXX99,JPBS(L),  ! to avoid IPBS(L) exceed bound
     & ISPBS(L),NPFORS,NPSERS(L),                  ! & all similar statements should be corrected!, 3/5/99
     $   NCSERS(L,1),NCSERS(L,2),NCSERS(L,3),      ! & yes, I did.
     $   NTOXSRC,NSEDSRC,NSNDSRC,NTSCRS(L)
       if(ISO.NE.0) exit
      IPBS(L)=IXX99
      ICBs(L)=IPBs(L)
      JCBs(L)=JPBs(L)
      NCSERS(L,4)=0      ! shell fish
      write(nefdcout,7213)         IPBS(L),JPBS(L),
     * ISPBS(L),NPFORS,NPSERS(L)
     &                                      ,
     $   NCSERS(L,1),NCSERS(L,2),NCSERS(L,3),
     $   NTOXSRC,NSEDSRC,NSNDSRC,NTSCRS(L)
 7213 format(999i6)
c search
      IF(NPFORs.GT.NPFOR) go to 987
      NPSER=MAX(NPSER,NPSERs(L))
      NCSER(1)=MAX(NCSER(1),NCSERs(L,1))
      NCSER(2)=MAX(NCSER(2),NCSERs(L,2))
      NCSER(3)=MAX(NCSER(3),NCSERs(L,3))
      NTOXSER=MAX(NTOXSER,NTOXSRC)
      NSEDSER=MAX(NSEDSER,NSEDSRC)
      NSNDSER=MAX(NSNDSER,NSNDSRC)
c tides
      if(NPFORS.gt.0) then
      DO M=1,MTIDE
      RAD=PI2*PFPH(NPFORS,M)/TCP(M)
      AMP=G*PFAM(NPFORS,M)
      PCBS(L,M)=AMP*COS(RAD)
      PSBS(L,M)=AMP*SIN(RAD)
      END DO
      endif
c conc.
       DO N=1,NTOX
        M=MSVTOX(N)
        NCSERS(L,M)=NTOXSRC
       END DO
       DO N=1,NSED
        M=MSVSED(N)
        NCSERS(L,M)=NSEDSRC
       END DO
       DO N=1,NSND
        M=MSVSND(N)
        NCSERS(L,M)=NSNDSRC
       END DO
c
       enddo   ! DO L
c
       backspace (nefdcwin)
       read(nefdcwin,997) aaa
       if(aaa(1).ne.'C'.and.aaa(1).ne.'c') go to 987
       write(nefdcout,997) aaa
c
       NPBS=L-1
       NCBS=NPBS
       write(6,*) "NPBS=NCBS=  ",NPBS
       if(NPBS.GT.NPBSM) THEN
       write(6,*) "Error in NPBS,MPBSM: ",NPBS,NPBSM
       go to 987
       endif
C     CONSTANT bottom CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
      MMAX=4+NTOX
      DO L=1,NCBS
      DO M=1,MMAX
      CBS(L,1,M)=0.0
      END DO
      enddo
C     CONSTANT bottom CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBS
      do m=Mmin,mmax
      CBS(L,1,M)=0.0
      END DO
      enddo
C     CONSTANT surface CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
C
      MMAX=4+NTOX
      DO L=1,NCBS
      do M=1,MMAX
      CBS(L,2,M)=0.0
      END DO
      enddo
C     CONSTANT surface CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBS
      do M=MMIN,MMAX
      CBS(L,2,M)=0.0
      END DO
      enddo
C
C13 LOCATION OF PERIODIC FORCING (TIDAL) SURF ELEV OR PRESSURE, and CONC BC'S ON
C     West  BOUNDARIES
C
       ncard=13
       call PasteC(nefdcwin,nefdcout)
c
      DO L=1,99999
      READ (nefdcwin,*,IOSTAT=ISO) IXX99,JPBw(L),ISPBw(L),
     & NPFORw,NPSERw(L),
     $   NCSERw(L,1),NCSERw(L,2),NCSERw(L,3),
     $   NTOXSRC,NSEDSRC,NSNDSRC,NTSCRw(L)
       if(ISO.NE.0) exit
      IPBw(L)=IXX99
      ICBw(L)=IPBw(L)
      JCBw(L)=JPBw(L)
      NCSERw(L,4)=0
      write(nefdcout,7213)         IPBw(L),JPBw(L),ISPBw(L),
     & NPFORw,NPSERw(L)
     &                                      ,
     $   NCSERw(L,1),NCSERw(L,2),NCSERw(L,3),
     $   NTOXSRC,NSEDSRC,NSNDSRC,NTSCRw(L)
c search
      IF(NPFORw.GT.NPFOR) go to 987
      NPSER=MAX(NPSER,NPSERw(L))
      NCSER(1)=MAX(NCSER(1),NCSERw(L,1))
      NCSER(2)=MAX(NCSER(2),NCSERw(L,2))
      NCSER(3)=MAX(NCSER(3),NCSERw(L,3))
      NTOXSER=MAX(NTOXSER,NTOXSRC)
      NSEDSER=MAX(NSEDSER,NSEDSRC)
      NSNDSER=MAX(NSNDSER,NSNDSRC)
c tides
      if(NPFORw.gt.0) then
      DO M=1,MTIDE
      RAD=PI2*PFPH(NPFORw,M)/TCP(M)
      AMP=G*PFAM(NPFORw,M)
      PCBw(L,M)=AMP*COS(RAD)
      PSBw(L,M)=AMP*SIN(RAD)
      END DO
      ENDIF
c conc.
       DO N=1,NTOX
        M=MSVTOX(N)
        NCSERw(L,M)=NTOXSRC
       END DO
       DO N=1,NSED
        M=MSVSED(N)
        NCSERw(L,M)=NSEDSRC
       END DO
       DO N=1,NSND
        M=MSVSND(N)
        NCSERw(L,M)=NSNDSRC
       END DO
       enddo
c check
       backspace (nefdcwin)
       read(nefdcwin,997) aaa
       if(aaa(1).ne.'C'.and.aaa(1).ne.'c') go to 987
       write(nefdcout,997) aaa
c process
       NPBw=L-1
       NCBw=NPBw
       write(6,*) "NPBw=NCBw=  ",NPBw
       if(NPBw.GT.NPBwM) THEN
       write(6,*) "Error in NPBw,MPBwM: ",NPBw,NPBwM
       go to 987
       endif
C     CONSTANT bottom CONCENTRATION ON west  CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
      MMAX=4+NTOX
      DO L=1,NCBw
      DO M=1,MMAX
      CBw(L,1,M)=0.0
      END DO
      enddo
C     CONSTANT bottom CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBw
      do m=Mmin,mmax
      CBw(L,1,M)=0.0
      END DO
      enddo
C     CONSTANT surface CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
C
      MMAX=4+NTOX
      DO L=1,NCBw
      do M=1,MMAX
      CBw(L,2,M)=0.0
      END DO
      enddo
C     CONSTANT surface CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBw
      do M=MMIN,MMAX
      CBw(L,2,M)=0.0
      END DO
      enddo
C
C14 LOCATION OF PERIODIC FORCING (TIDAL) SURF ELEV OR PRESSURE, and CONC BC'S ON
C     East  BOUNDARIES
C
       ncard=14
       call PasteC(nefdcwin,nefdcout)
c
      DO L=1,99999
      READ (nefdcwin,*,IOSTAT=ISO) IXX99,JPBe(L),
     & ISPBe(L),NPFORe,NPSERe(L)
     &                                     ,
     $   NCSERe(L,1),NCSERe(L,2),NCSERe(L,3),
     $   NTOXSRC,NSEDSRC,NSNDSRC,NTSCRe(L)
       if(ISO.NE.0) exit
      IPBe(L)=IXX99
      ICBe(L)=IPBe(L)
      JCBe(L)=JPBe(L)
      NCSERe(L,4)=0
      write(nefdcout,7213)         IPBe(L),JPBe(L),
     & ISPBe(L),NPFORe,NPSERe(L)
     &                                      ,
     $   NCSERe(L,1),NCSERe(L,2),NCSERe(L,3),
     $   NTOXSRC,NSEDSRC,NSNDSRC ,NTSCRe(L)
c search
      IF(NPFORe.GT.NPFOR) go to 987
      NPSER=MAX(NPSER,NPSERe(L))
      NCSER(1)=MAX(NCSER(1),NCSERe(L,1))
      NCSER(2)=MAX(NCSER(2),NCSERe(L,2))
      NCSER(3)=MAX(NCSER(3),NCSERe(L,3))
      NTOXSER=MAX(NTOXSER,NTOXSRC)
      NSEDSER=MAX(NSEDSER,NSEDSRC)
      NSNDSER=MAX(NSNDSER,NSNDSRC)
c tides
      if(NPFORe.gt.0) then
      DO M=1,MTIDE
      RAD=PI2*PFPH(NPFORe,M)/TCP(M)
      AMP=G*PFAM(NPFORe,M)
      PCBe(L,M)=AMP*COS(RAD)
      PSBe(L,M)=AMP*SIN(RAD)
      END DO
      endif
c conc.
       DO N=1,NTOX
        M=MSVTOX(N)
        NCSERe(L,M)=NTOXSRC
       END DO
       DO N=1,NSED
        M=MSVSED(N)
        NCSERe(L,M)=NSEDSRC
       END DO
       DO N=1,NSND
        M=MSVSND(N)
        NCSERe(L,M)=NSNDSRC
       END DO
       enddo
c check
       backspace (nefdcwin)
       read(nefdcwin,997) aaa
       if(aaa(1).ne.'C'.and.aaa(1).ne.'c') go to 987
       write(nefdcout,997) aaa
c process
       NPBe=L-1
       NCBe=NPBe
       write(6,*) "NPBe=NCBe=  ",NPBe
       if(NPBe.GT.NPBeM) THEN
       write(6,*) "Error in NPBe,MPBeM: ",NPBe,NPBeM
       go to 987
       endif
C     CONSTANT bottom CONCENTRATION ON west  CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
      MMAX=4+NTOX
      DO L=1,NCBe
      DO M=1,MMAX
      CBe(L,1,M)=0.0
      END DO
      enddo
C     CONSTANT bottom CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBe
      do m=Mmin,mmax
      CBe(L,1,M)=0.0
      END DO
      enddo
C     CONSTANT surface CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
C
      MMAX=4+NTOX
      DO L=1,NCBe
      do M=1,MMAX
      CBe(L,2,M)=0.0
      END DO
      enddo
C     CONSTANT surface CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBe
      do M=MMIN,MMAX
      CBe(L,2,M)=0.0
      END DO
      enddo
C
C15 LOCATION OF PERIODIC FORCING (TIDAL) SURF ELEV OR PRESSURE, and CONC BC'S ON
C     North BOUNDARIES
C
       ncard=15
       call PasteC(nefdcwin,nefdcout)
c
      DO L=1,99999
      READ (nefdcwin,*,IOSTAT=ISO) IXX99,JPBn(L),
     & ISPBn(L),NPFORn,NPSERn(L),
     $   NCSERn(L,1),NCSERn(L,2),NCSERn(L,3),
     $   NTOXSRC,NSEDSRC,NSNDSRC,NTSCRn(L)
       if(ISO.NE.0) exit
      IPBn(L)=IXX99
      ICBn(L)=IPBn(L)
      JCBn(L)=JPBn(L)
      NCSERn(L,4)=0
      write(nefdcout,7213)         IPBn(L),JPBn(L),
     & ISPBn(L),NPFORn,NPSERn(L),
     $   NCSERn(L,1),NCSERn(L,2),NCSERn(L,3),
     $   NTOXSRC,NSEDSRC,NSNDSRC ,NTSCRn(L)
c search
      IF(NPFORn.GT.NPFOR) then
      write(6,*) "error: NPFORn.GT.NPFOR  ", NPFORn,NPFOR
      go to 987
      endif
      NPSER=MAX(NPSER,NPSERn(L))
      NCSER(1)=MAX(NCSER(1),NCSERn(L,1))
      NCSER(2)=MAX(NCSER(2),NCSERn(L,2))
      NCSER(3)=MAX(NCSER(3),NCSERn(L,3))
      NTOXSER=MAX(NTOXSER,NTOXSRC)
      NSEDSER=MAX(NSEDSER,NSEDSRC)
      NSNDSER=MAX(NSNDSER,NSNDSRC)

c tides
      if(NPFORn.gt.0) then
      DO M=1,MTIDE
      RAD=PI2*PFPH(NPFORn,M)/TCP(M)
      AMP=G*PFAM(NPFORn,M)
      PCBn(L,M)=AMP*COS(RAD)
      PSBn(L,M)=AMP*SIN(RAD)
      END DO
      endif
c conc.
       DO N=1,NTOX
        M=MSVTOX(N)
        NCSERn(L,M)=NTOXSRC
       END DO
       DO N=1,NSED
        M=MSVSED(N)
        NCSERn(L,M)=NSEDSRC
       END DO
       DO N=1,NSND
        M=MSVSND(N)
        NCSERn(L,M)=NSNDSRC
       END DO
c
       enddo   ! DO L
c
       backspace (nefdcwin)
       read(nefdcwin,997) aaa
       write(nefdcout,997) aaa
       if(aaa(1).ne.'C'.and.aaa(1).ne.'c') go to 987
c
       NPBn=L-1
       NCBn=NPBn
       write(6,*) "NPBn=NCBn=  ",NPBn
       if(NPBn.GT.NPBnM) THEN
       write(6,*) "Error in NPBn,MPBnM: ",NPBn,NPBnM
       go to 987
       endif
c modify & write -----------------
7801  format(" NQSER,NSSER,NTSER,NDSER,NTXSER,NDSER,NNSER,NPSER", /
     &,999i6)
      write(nefdcout,7801) NQSER,
     & NCSER(1),NCSER(2),NCSER(3),NTOXSER,NSEDSER,NSNDSER,NPSER
c     write(6       ,7801) NQSER,
c    & NCSER(1),NCSER(2),NCSER(3),NTOXSER,NSEDSER,NSNDSER,NPSER
      if(ISTRAN(1).EQ.0) NCSER(1)=0
      if(ISTRAN(2).EQ.0) NCSER(2)=0
      if(ISTRAN(3).EQ.0) NCSER(3)=0
      if(ISTRAN(5).EQ.0) NTOXSER=0
      if(ISTRAN(6).EQ.0) NSEDSER=0
      if(ISTRAN(7).EQ.0) NSNDSER=0
      write(6       ,7801) NQSER,
     & NCSER(1),NCSER(2),NCSER(3),NTOXSER,NSEDSER,NSNDSER,NPSER
c
      MTMP=4
      DO N=1,NTOX
       MTMP=MTMP+1
       MSVTOX(N)=MTMP
      END DO
      DO N=1,NSED
       MTMP=MTMP+1
       MSVSED(N)=MTMP
c         write(6,*) "???? MSVSeD  ",n,mtmp,msvsed(N)
      END DO
      DO N=1,NSND
       MTMP=MTMP+1
       MSVSND(N)=MTMP
      END DO
      DO N=1,NTOX
       M=MSVTOX(N)
       NCSER(M)=NTOXSER
      END DO
      DO N=1,NSED
       M=MSVSED(N)
       NCSER(M)=NSEDSER
      END DO
      DO N=1,NSND
       M=MSVSND(N)
       NCSER(M)=NSNDSER
      END DO
c----------------
C     CONSTANT bottom CONCENTRATION ON west  CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
      MMAX=4+NTOX
      DO L=1,NCBn
      DO M=1,MMAX
      CBn(L,1,M)=0.0
      END DO
      enddo
C     CONSTANT bottom CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBn
      do m=Mmin,mmax
      CBn(L,1,M)=0.0
      END DO
      enddo
C     CONSTANT surface CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
C
      MMAX=4+NTOX
      DO L=1,NCBn
      do M=1,MMAX
      CBn(L,2,M)=0.0
      END DO
      enddo
C     CONSTANT surface CONCENTRATION ON SOUTH CONC BOUNDARIES
C     SED(1 TO NSED),SND(1,NSND)
C
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBn
      do M=MMIN,MMAX
      CBn(L,2,M)=0.0
      END DO
      enddo
C
C16 INITIAL DRIFTER POSITIONS (FOR USE WITH SUB DRIFTER2)
C
       ncard=16
       call PasteC(nefdcwin,nefdcout)
c
       do ncountx=1,99999
       read(nefdcwin,*,IOSTAT=ISO) RI(ncountx),RJ(ncountx)
       if(ISO.NE.0) exit
       write(nefdcout,7416) RI(ncountx),RJ(ncountx)
 7416  format(3x,2f7.1)
       enddo
c
       backspace (nefdcwin)
       read(nefdcwin,997) aaa
       if(aaa(1).ne.'C'.and.aaa(1).ne.'c') go to 987
       write(nefdcout,997) aaa
c
c      go to 2000
       NPD=max(ncountx-1,0)
       if(NPD.gt.0) ISPD=2                     ! call drifter2
       if(NPD.gt.NPDM) then                    ! ensure NPD <= NPDM
       write(6,*) "Error in NPD,NPDm: ", NPD,NPDM
       go to 987
       endif
C
      NPDRT=max(1,NPDRT)                       ! efdc.inp C67
c     NWPD is given in C19
C
C17**  Decide what variables to be saved
c
      NCARD=17
      call PasteC(nefdcwin,nefdcout)
C
      DO L=1,99999
      READ(nefdcwin,*,IOSTAT=ISO) IGrADS,
     1 IDUMPe   ,IDUMPu   ,IDUMPv   ,IDUMPw   ,IDUMPs
     1,IDUMPt   ,IDUMPd   ,IDUMPc   ,IDUMPn   ,IDUMPx
     1,IDUMPb   ,IDUMPfx   ,IDUMPfy
       if(ISO.NE.0) exit
      write(nefdcout,7216) IGrADS,
     1 IDUMPe   ,IDUMPu   ,IDUMPv   ,IDUMPw   ,IDUMPs
     1,IDUMPt   ,IDUMPd   ,IDUMPc   ,IDUMPn   ,IDUMPx
     1,IDUMPb   ,IDUMPfx   ,IDUMPfy
 7216 format(8x,i2,i6,999i2)
c
      IF(ISTRAN(1).EQ.0) IDUMPs   =0           ! ensure consistency
      IF(ISTRAN(2).EQ.0) IDUMPt   =0
      IF(ISTRAN(3).EQ.0) IDUMPd   =0
c     IF(ISTRAN(4).EQ.0) IDUMPs(L)=0           ! shell fish
      IF(ISTRAN(5).EQ.0) IDUMPx   =0
      IF(ISTRAN(6).EQ.0) IDUMPc   =0
      IF(ISTRAN(7).EQ.0) IDUMPn   =0
c     IF(ISTRAN(8).EQ.0) IDUMPw(L)=0           ! WQ
c check IDUMPb       tox       sed       snd        WQ
      ixx98=max(ISTRAN(5),ISTRAN(6),ISTRAN(7),ISTRAN(8))
      if(ixx98.eq.0) IDUMPb=0
c
       IDUMP2=IGrADS
       enddo
c
       backspace (nefdcwin)
       read(nefdcwin,997) aaa
       if(aaa(1).ne.'C'.and.aaa(1).ne.'c') go to 987
       write(nefdcout,997) aaa
C
C18** Specifying locations FOR WRITING TO TIME SERIES FILES (hyts.bin)
c
      NCARD=18
      call PasteC(nefdcwin,nefdcout)
C
      DO L=1,99999
      READ(nefdcwin,*,IOSTAT=ISO) ILTMSR(L),JLTMSR(L),CLTMSR(L)
       if(ISO.NE.0) exit
      write(nefdcout,7246)        ILTMSR(L),JLTMSR(L),CLTMSR(L)
 7246 format(2i8,5x,a20)
       enddo
c
       backspace (nefdcwin)
       read(nefdcwin,997) aaa
       if(aaa(1).ne.'C'.and.aaa(1).ne.'c') go to 987
       write(nefdcout,997) aaa
c
       MLTMSR=L-1
       write(6,*) "MLTMSR=  ",MLTMSR
       if(MLTMSR.gt.MLTMSRM) then
       write(6,*) "Error in MLTMSR,MLTMSRM: ", MLTMSR,MLTMSRM
       go to 987
       endif
c make it consistent with efdc.inp C79, & C82
       do M=1,MLTMSR
       NTSSSS(M)=1
       MTMSRP(M)=IDUMPe
       MTMSRC(M)=max(IDUMPs,IDUMPt,IDUMPd,IDUMPx,IDUMPc,IDUMPn)
!       MTMSRA(M)=0
       MTMSRUE(M)=0
       MTMSRUT(M)=0
       MTMSRU(M)=max(IDUMPu   ,IDUMPv   )
       MTMSRQE(M)=0
       MTMSRQ(M)=0
c      MTMSRQE(M)=1  ! Hardware for Wister, JI, 7/30/99
c      MTMSRQ(M)=1   ! Hardware for Wister, Ji, 7/30/99
       enddo

C
C19** Specify averaging and saving for time series (hyts.bin) and
C     3D data (hy3d.bin) for GrADS
c
      NCARD=19
       call PasteC(nefdcwin,nefdcout)
c
      READ(nefdcwin,*,IOSTAT=ISO) SaveTs,StartDay,EndDay,Aver3D, Save3D
       if(ISO.NE.0) go to 987
      write(nefdcout,7217) SaveTs,StartDay,EndDay,Aver3D, Save3D
 7217  format( 999f10.2)

      Itimes=NTSPTC*savets/24.0+0.5
      Iaverage=NTSPTC*aver3d/24.0+0.5
      isave=NTSPTC*save3d/24.0+0.5
      if(Itimes.le.0) itimes=1
      if(iaverage.le.0.or.save3d.le.0.0) iaverage=1
      if(isave.le.0) isave=1
c change efdc.inp's
      NWTMSR=ITIMES
      NTSMMT=ITIMES  ! efdc.inp, C7, for budget5.for, sediment flux, 9/10/00
      NWPD  =ITIMES  ! efdc.inp, C67, for drifter, 4/1/2001
c
c Notes: Itimes must be even number to save WQ time series in the same way as hydro,
c because WQ is calculated at every even steps!!, Jeff Ji, 8/11/99
c
      if(ISTRAN(8).GT.0.and.MOD(ITIMES,2).ne.0) then
      write(6,771) NCARD
771   format( "efdcwin.inp, ITIMES should be an even number for",
     & " WQ model, Card ", i8)
      stop
      endif
c checking
      if(iaverage.gt.isave) then
      write(6,*) "Error: Iaverage>isave ", Iaverage, isave
      stop
      endif
c
c calculate 3D array saving number
      dayx=EndDay-StartDay
      dayx2=ntc
      JHM=amin1(dayx,dayx2)*NTSPTC/ISAVE+0.5
      if(JHM.ge.1000) then
      write(6,7417) ITIMES,StartDay,EndDay,IAVERAGE,ISAVE,JHM
7417  format("Huge hy3d.bin might be produced! ", i8,2f10.2,2i8,i8)
      endif
c
c make ihist(7000,2)
      ISKIP=startday*NTSPTC+0.5  ! skip startday days from beginning
      ncount=0
      nsteps=NTC*NTSPTC
      do 700 n=1,nsteps
      if(n.le.ISKIP) go to 700
      if(mod(n,isave).eq.0) then
      ncount=ncount+1
c     write(199,*) ncount  ! ji
      ihist(ncount,2)=n                ! save and average upper bound
      IHIST(ncount,1)=n-IAVerage+1     ! average lower bound
      endif
 700  continue
c write
      write(nefdcout,917) ncount,IAVERAGE
	do i=1,ncount
      write(nefdcout,917) i,ihist(i,1),ihist(i,2)
	enddo
 917  format(10i8)
C
      GO TO 2000
C
c---------------------------------------
c
 987    write(6,*) "Read error: EFDCWIN.INP, Card ", Ncard
       stop
c
 977   write(6,*) "EFDCWIN.INP does not exist and is not used"
       go to 2001
c
 2000  continue
       NCSER(1)=NCSER(1)+NLCDA            !J.S.
       close(nefdcwin)
       close(nefdcout)
       write(6,*) "EFDCWIN.INP is read successfully."
 2001  continue
 997   format(120a1)
 967    format(999e9.2)
       return
       end
