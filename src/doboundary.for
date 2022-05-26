
      subroutine doboundary(Os,Asss,dotem,dosal,dowindst,UM, VM,
	1  HP,DZC,iop,IWQKA_1)
c      INCLUDE 'efdc.par'
c      INCLUDE 'efdc.cmn'
c      DIMENSION OS(LCMWQ)
	real Os,Asss,dotem,dosal,dowindst,UM,VM,HP,DZC
	real TWQ,SWQ,windrea
	integer IWQKA_1
c
c          TIME=DT*FLOAT(N)+TCON*TBEGIN
c          TIME=TIME/TCON
C	if K.eq.KC  THEN ! surface

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C  The following modification to the D.O. saturation calculation made
C  by J.M. Hamrick / M.R. Morton on 03/08/97.  See Chapra (1997) pg. 361-364.
C
C        TWQ(L)=TEM(L,K)
C	  SWQ(L)=MAX(SAL(L,K), 0.0)
        if(iop.eq.1) then
        TWQ=dotem
	  SWQ=MAX(dosal,0.0)
C
        TVAL1=1./(TWQ+273.15)
        TVAL2=TVAL1*TVAL1
        TVAL3=TVAL1*TVAL2
        TVAL4=TVAL2*TVAL2
        RLNSAT1=-139.3441+(1.575701E+5*TVAL1)-(6.642308E+7*TVAL2)
     $                   +(1.2438E+10*TVAL3)-(8.621949E+11*TVAL4)
        RLNSAT2=RLNSAT1-SWQ*( 1.7674E-2-(1.0754E+1*TVAL1)
     $                           +(2.1407E+3*TVAL2) )
c        WQDOS(L) = EXP(RLNSAT2)
        Os = EXP(RLNSAT2)
        endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Do not allow wind speeds above 11 m/sec in the following equation:
        windrea = dowindst
        if (dowindst .gt. 11.0) windrea = 11.0
        WQWREA = 0.728*SQRT(windrea) + (0.0372*windrea - 0.317)*windrea
c        WQWREA = 0.728*SQRT(WINDST(L))
c     $           +(0.0372*WINDST(L) - 0.317)*WINDST(L)
CTt       WQWREA = 0.0
chong-------------------------------------------------------
!       IWQKA_1=2
        WQKRO=3.933
chong-------------------------------------------------------
c
c MRM 04/29/1999  User specifies constant reaeration WQKRO:
c
          IF (IWQKA_1 .EQ. 0) THEN
            WQVREA = WQKRO
            WQWREA = 0.0
          END IF
c
c MRM 04/12/1999  Constant reaeration due to water velocity,
c                 wind velocity computed above:
c
          IF (IWQKA_1 .EQ. 1) THEN
            WQVREA = WQKRO
          END IF
c
c MRM 03/06/1999  O'Connor-Dobbins (1958) equation for reaeration is:
c    WQKRO = 3.933 typically
c
          IF (IWQKA_1 .EQ. 2) THEN
           WQKRO = 3.933
           umrm=UM
	     VMRM=VM   
            xmrm = SQRT(umrm*umrm + vmrm*vmrm)
C            WQVREA = WQKRO * xmrm**0.5 / HP(L)**1.5
            WQVREA = WQKRO * xmrm**0.5 / HP**1.5
          END IF
c
c MRM 04/12/1999  Owens and Gibbs (1964) reaeration equation:
c    WQKRO = 5.32 typically
c
          IF (IWQKA_1 .EQ. 3) THEN
            WQKRO = 5.32
            umrm = UM
            vmrm = VM
            xmrm = SQRT(umrm*umrm + vmrm*vmrm)
            WQVREA = WQKRO * xmrm**0.67 / HP**1.85
          END IF

c
c now combine reaeration due to water velocity and wind stress:
c
c          WQVREA = WQVREA * REAC(IWQZMAP(L,K))
c          WQWREA = WQWREA * REAC(IWQZMAP(L,K))

C-----------------------------------------

        DZWQ = 1.0 / (DZC*HP)
	  TT20 = dotem-20.0
	  theta=1.06
chong	  WQKTR = 1.06  ! 1.06~1.08
chong   WQTDKR(M) = WQKTR**TT20

	Asss= (WQVREA + WQWREA) * DZWQ* theta**TT20 
chong          WQP19(L) = - (WQVREA + WQWREA) * DZWQ(L)* WQTDKR(IWQT(L))

      END

