C 
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SMMBE(LA,KC,DT,N,TCON,TBEGIN,NTSPTC,NCSTEP,SECDLAST)
C
C**********************************************************************C
C
C  CONTROL SUBROUTINE FOR SEDIMENT COMPONENT OF WATER QUALITY MODEL
C
C  ORGINALLY CODED BY K.-Y. PARK
C  OPTIMIZED AND MODIFIED BY J. M. HAMRICK
C  LAST MODIFIED BY M.R. MORTON  2 JUNE 1999
C  06/02/99 MRM: added BFO2sum, BFNH4sum, BFNO3sum, BFPO4sum, BFCODsum,
C   and BFSADsum arrays to keep track of benthic fluxes at all cells
C   for later storage in binary file (WQSDTS.bin).
C
C**********************************************************************C
C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
C
      COMMON/TSM/SMW12(LCMWQ),SMKL12(LCMWQ),XSMO20(LCMWQ),SMTMP(LCMWQ),
     $           ISMT(LCMWQ)
      REAL*8 SECDLAST  
      NDMWQ=1
	LDM=LA
C
C  Sed temp., & find an index for look-up table for temperature
C  dependency
C: SM1DIFT(IZ)=SMDIFT*SMDTOH/SMHSED,SM2DIFT=1/(1+SM1DIFT)
C
      DO L=2,LA
        IZ = ISMZMAP(L)                                                  ! ISMZMAP=1, usually
        SMT(L) = (SMT(L) + SM1DIFT(IZ)*TEMWQ(L,1)) * SM2DIFT(IZ)           ! get new sediment model temperature
        ISMT(L) = 10.0*SMT(L) + 51    !must be consistent with the look-up table:-5C-> 50C in smrin1 & wqwin, Ji, 8/29/02
        IF (ISMT(L).LT.1 .OR. ISMT(L).GT.NWQTD) THEN
c MRM +++++++++ added by M. Morton 07/24/98
          timtmp = (dt*float(n) + tcon*tbegin)/86400.0
          OPEN(1,FILE='error.log',ACCESS='APPEND',STATUS='UNKNOWN')
          write(1,911) timtmp, L, ILW(L), JLW(L), temwq(L,1), SMT(L)
911       format(/,'ERROR!! invalid sediment temperature',/,
     +     'TIME, L, I, J, TEMwq(L,1), SMT(L) = ', f10.5, 3i4, 2f10.4)
          close(1)
c MRM +++++++++ added by M. Morton 07/24/98
          PRINT*, 'L, TEMwq(L,1), SMT(L) = ', L,TEMwq(L,1),SMT(L)
c
c MRM: the following STOP was commented out by M. Morton 07/24/98 and
c ISMT(L) was set equal to the bounds if it exceeded the bounds, thus
c the model is now allowed to continue to run.  The user should check
c the error.log file for sediment temperatures out of range.
c          STOP 'ERROR!! invalid sediment temperature'

          if (ismt(L) .lt. 1) ismt(L)=1
          if (ismt(L) .gt. nwqtd) ismt(L) = nwqtd
        END IF
      END DO
C
C1 Depositional flux from water column (SMDFN,SMDFP,SMDFC)
C: SMHODT(IZ)=SMHSED/DTWQ,SMDTOH(IZ)=DTWQ/SMHSED
C
      DO M=1,NSMG                                                       ! =3 = G1, G2, G3
        DO L=2,LA
          SMDFNA = SMFNBC(M)*WQANCC*WQDFBC(L)                           ! (4-3), last term
     *      + SMFNBD(M)*WQANCD*WQDFBD(L) + SMFNBG(M)*WQANCG*WQDFBG(L)
          SMDFPA = ( SMFPBC(M)*WQDFBC(L) + SMFPBD(M)*WQDFBD(L)          ! (4-4), last term
     *      + SMFPBG(M)*WQDFBG(L) ) * WQAPC(L)
          SMDFCA = SMFCBC(M)*WQDFBC(L) + SMFCBD(M)*WQDFBD(L)            ! (4-2), last term
     *      + SMFCBG(M)*WQDFBG(L)
          SMDFN(L,M) =SCBWQ(L)*( SMDFNA + SMFNR(ISMZMAP(L),M)*WQDFRN(L)) ! (4-3), second term
          SMDFP(L,M) =SCBWQ(L)*( SMDFPA + SMFPR(ISMZMAP(L),M)*WQDFRP(L))
          SMDFC(L,M) =SCBWQ(L)*( SMDFCA + SMFCR(ISMZMAP(L),M)*WQDFRC(L))
        END DO
      END DO

      do L=2,la
        SMDFN(L,1) = scbwq(L)*( SMDFN(L,1) + WQDFLN(L) )                  ! (4-3), first term. No fraction needed,
        SMDFP(L,1) = scbwq(L)*( SMDFP(L,1) + WQDFLP(L) )                  ! since all labile ones in WC go to labile ones in SM
        SMDFC(L,1) = scbwq(L)*( SMDFC(L,1) + WQDFLC(L) )
      end do
C
C2 Diagenesis Eq and diagenesis flux (SMDGFN,SMDGFP,SMDGFC)
C: SMW2 in m/d,SMW2DTOH(IZ)=1.0+SMW2*SMDTOH
C
c Integrate (4-6)
c
      DO M=1,NSMG
        DO L=2,LA                                                         !SMDTOH=dt/H 
        SMPON(L,M)=SCBWQ(L)*(SMPON(L,M) + SMDFN(L,M)*SMDTOH(ISMZMAP(L)))  !PON=particulate organic nitrogen, class Gm
     *      / (SMW2DTOH(ISMZMAP(L)) + SMTDND(ISMT(L),M)*DTWQ+ 1.E-18)     ! (4-6)
        SMPOP(L,M)=SCBWQ(L)*(SMPOP(L,M) + SMDFP(L,M)*SMDTOH(ISMZMAP(L)))
     *      / (SMW2DTOH(ISMZMAP(L)) + SMTDPD(ISMT(L),M)*DTWQ+ 1.E-18)
        SMPOC(L,M)=SCBWQ(L)*(SMPOC(L,M) + SMDFC(L,M)*SMDTOH(ISMZMAP(L)))
     *      / (SMW2DTOH(ISMZMAP(L)) + SMTDCD(ISMT(L),M)*DTWQ+ 1.E-18)
        END DO
      END DO
C
!      DO ND=1,NDMWQ
!      LF=2+(ND-1)*LDMWQ
!      LL=LF+LDM-1
      DO L=2,LA !LF,LL
        SMDGFN(L) = SCBWQ(L)*SMHSED(ISMZMAP(L)) *
     *    (SMTDND(ISMT(L),1)*SMPON(L,1) + SMTDND(ISMT(L),2)*SMPON(L,2))   ! (4-6)
        SMDGFP(L) = SCBWQ(L)*SMHSED(ISMZMAP(L)) *
     *    (SMTDPD(ISMT(L),1)*SMPOP(L,1) + SMTDPD(ISMT(L),2)*SMPOP(L,2))
        SMDGFC(L) = SCBWQ(L)*SMHSED(ISMZMAP(L)) *
     *    (SMTDCD(ISMT(L),1)*SMPOC(L,1) + SMTDCD(ISMT(L),2)*SMPOC(L,2))
C
C3 Mass-balance eq's for fluxes
C
C  Common parameters: SMBST1=1/(1+SMKBST*DTWQ),SM1OKMDP=1/SMKMDP
C: ISMTDMBS=NINT(SMTDMBS/DTWQ),ISMTCMBS=NINT(SMTCMBS/DTWQ)
C: use SMTMP(L) to store old SMBST(L)
C
c       XSMO20(L) = MAX( WQV(L,1,19), 0.001 ) ! bug, Ji, 9/24/99. Too small -> large P flux.
        XSMO20(L) = MAX( WQV(L,1,19), 0.1)  ! Changed for Tenkiller, 11/30/99

        SMTMP(L) = SCBWQ(L)*SMBST(L)
        IF (XSMO20(L).LT.SMKMDP) THEN
          SMBST(L) =SCBWQ(L)*(SMTMP(L)
     $                        +DTWQ*(1.0-XSMO20(L)*SM1OKMDP)) * SMBST1
         ELSE
          SMBST(L) = SCBWQ(L)*SMBST(L)*SMBST1
        END IF
      END DO
!      END DO
C
      IF (ISMHYST.EQ.1) THEN
!        DO ND=1,NDMWQ
!        LF=2+(ND-1)*LDMWQ
!        LL=LF+LDM-1
        DO L=2,LA !LF,LL
          IF(SCBWQ(L).GT.0.5) THEN
          IF (SMHYST(L)) THEN
            IF (XSMO20(L).GE.SMO2BS) ISMHYPD(L) = ISMHYPD(L) - 1
            IF (ISMHYPD(L).EQ.0) THEN
              SMHYST(L) = .FALSE.
              ISMHYPD(L) = 0
            END IF
            SMBST(L) = SMTMP(L)
           ELSE
            IF (XSMO20(L).LT.SMO2BS) ISMHYPD(L) = ISMHYPD(L) + 1
            IF (ISMHYPD(L).EQ.ISMTCMBS) THEN
              SMHYST(L) = .TRUE.
              ISMHYPD(L) = ISMTDMBS
            END IF
          END IF
          END IF
        END DO
!        END DO
      END IF
C
C: SMDP(IZ)=SMDP/(SMHSED*SMPOCR),SMDD(IZ)=SMDD/SMHSED
C: SMDPMIN(IZ)=SMDPMIN/SMHSED
C
!      DO ND=1,NDMWQ
!      LF=2+(ND-1)*LDMWQ
!      LL=LF+LDM-1
      DO L=2,LA ! LF,LL
        IZ = ISMZMAP(L)
        SMW12(L) = SMDP(IZ)*SMTDDP(ISMT(L)) * SMPOC(L,1) * XSMO20(L)
     *    * (1.0-SMKBST*SMBST(L)) / (SMKMDP+XSMO20(L)+ 1.E-18)
     *      + SMDPMIN(IZ)
        SMKL12(L) = SMDD(IZ)*SMTDDD(ISMT(L)) + SMRBIBT*SMW12(L)
        SMKL12(L) =SCBWQ(L)*SMKL12(L)
      END DO
 !     END DO
C
C  NH4, NO3: SMKMO2N=SMKMO2N*2,SMKNH4=SMKNH4^2*SMKMNH4,
C: SMW2PHODT(IZ)=SMW2+SMHODT,SMK1NO3=SMK1NO3^2
C
      DO L=2,LA
      IF(SCBWQ(L).GT.0.5) THEN
        IZ = ISMZMAP(L)
        SMO20 = XSMO20(L)
        SK1NH4SM = ( SMKNH4(IZ)*SMTDNH4(ISMT(L)) * SMO20 )
     *    / ( (SMKMO2N+SMO20+ 1.E-12) * (SMKMNH4+SM1NH4(L)) )
        A1NH4SM = SMKL12(L)*SMFD1NH4 + SMW12(L)*SMFP1NH4 + SMW2(IZ)
        A2NH4SM = SMKL12(L)*SMFD2NH4 + SMW12(L)*SMFP2NH4
        A22NH4SM = A2NH4SM + SMW2PHODT(IZ)
        B1NH4SM = WQV(L,1,14)
        B2NH4SM = SMDGFN(L) + SMHODT(IZ)*SM2NH4(L)

        SK1NO3SM = SMK1NO3(IZ)*SMTDNO3(ISMT(L))
        A1NO3SM = SMKL12(L) + SMW2(IZ)
        A2NO3SM = SMKL12(L)
        RK2NO3SM = SMK2NO3(IZ)*SMTDNO3(ISMT(L))
        A22NO3SM = A2NO3SM + SMW2PHODT(IZ) + RK2NO3SM
        B1NO3SM = WQV(L,1,15)
        B2NO3SM = SMHODT(IZ)*SM2NO3(L)
C
C  H2S/CH4
C: SMK1H2S=(SMKD1HS^2*SMFD1H2S+SMKP1HS^2*SMFP1H2S)/(2*SMKMH2S)*SMTHH2S**(T-20)
C: SMTD1CH4(IT)=SMTD1CH4*20
C
        SMO2JC = SMO2C*SMDGFC(L)
        SMSAL0 = SALWQ(L,1)
        IF (SMSAL0.GT.SMCSHSCH) THEN
          SK1H2SSM = SMK1H2S(ISMT(L)) * SMO20
          A1H2SSM = SMKL12(L)*SMFD1H2S + SMW12(L)*SMFP1H2S + SMW2(IZ)
          A2H2SSM = SMKL12(L)*SMFD2H2S + SMW12(L)*SMFP2H2S
          A22H2SSM = A2H2SSM + SMW2PHODT(IZ)
          B1H2SSM = 0.0
          B2H2SSM = SMHODT(IZ)*SM2H2S(L)
         ELSE
          SMCH4S = (10.0 + HPWQ(L) + SMHSED(IZ)) * SMTD1CH4(ISMT(L))
     *      * SMKL12(L)
          SMK1CH4 = SMTD2CH4(ISMT(L))
        END IF
C
C back substitution to get SMSS
C
        SMSOD = ZBRENT(ISMERR)
        IF (ISMERR.EQ.1) THEN
 !         OPEN(1,FILE='zbrent.log',STATUS='UNKNOWN',ACCESS='APPEND')
 !         WRITE(1,401) ITNWQ,L,ILW(L),JLW(L),SMSOD,
 !    *      ' Root must be bracketed for ZBRENT  '
 !         CLOSE(1)
         ELSE IF (ISMERR.EQ.2) THEN
          OPEN(1,FILE='zbrent.log',STATUS='UNKNOWN',ACCESS='APPEND')
          WRITE(1,401) ITNWQ,L,ILW(L),JLW(L),SMSOD,
     *      ' ZBRENT exceeding maximum iterations'
          CLOSE(1)
        END IF
C
        SMSS(L) = RSMSS
        SM1NH4(L) = RSM1NH4
        SM2NH4(L) = RSM2NH4
        SM1NO3(L) = RSM1NO3
        SM2NO3(L) = RSM2NO3
        SM1H2S(L) = RSM1H2S
        SM2H2S(L) = RSM2H2S
c        WQBFO2(L) = -SMSOD
c mrm  SODMULT is a spatially variable adjustment factor for SOD:
        WQBFO2(L) = -SMSOD * SODmult(iz)*FLXPND(L,3)

        SMCSOD(L) = -CSODSM
        SMNSOD(L) = -RNSODSM
        SMJNIT(L) = RJNITSM
        SMJDEN(L) = RJDENSM
        SMJAQH2S(L) = AQJH2SSM + AQJCH4SM
        SMJGCH4(L) = GJCH4SM
C
        WQBFNH4(L) = SMSS(L) * (SMFD1NH4*SM1NH4(L) - WQV(L,1,14))* BFHN
     &               * FLXPND(L,2)  ! positive = sediment -> WC
        WQBFNO3(L) = SMSS(L) * (SM1NO3(L) - WQV(L,1,15))* BFNO
     &               * FLXPND(L,2)
        WQBFCOD(L) = SMJAQH2S(L) - SMSS(L)*WQV(L,1,18)
C       WQBFCOD(L) = SMJAQH2S(L)
C       WQBFCOD(L) = AQJH2SSM-SMSS(L)*WQV(L,1,18) + AQJCH4SM
C
        END IF
      END DO
C
C PO4
C
      DO L=2,LA
      IF(SCBWQ(L).GT.0.5) THEN
        IF (XSMO20(L).LT.SMCO2PO4) THEN
          SMP1PO4 = SMP2PO4
     $         * SMDP1PO4(ISMZMAP(L))**(XSMO20(L)/(SMCO2PO4+ 1.E-18))
         ELSE
          SMP1PO4 = SMP2PO4 * SMDP1PO4(ISMZMAP(L))
        END IF
        SMFD1PO4 = 1.0 / (1.0 + SMM1*SMP1PO4)
        SMFP1PO4 = 1.0 - SMFD1PO4
        A1PO4SM = SMKL12(L)*SMFD1PO4 + SMW12(L)*SMFP1PO4
     *    + SMW2(ISMZMAP(L))
        A2PO4SM = SMKL12(L)*SMFD2PO4 + SMW12(L)*SMFP2PO4
        A11PO4SM = SMSS(L)*SMFD1PO4 + A1PO4SM
        A22PO4SM = A2PO4SM + SMW2PHODT(ISMZMAP(L))
        B11PO4SM = SMSS(L) * WQPO4D(L,1)
        B22PO4SM = SMDGFP(L) + SMHODT(ISMZMAP(L))*SM2PO4(L)
        CALL SOLVSMBE(RSM1PO4,RSM2PO4,A11PO4SM,A22PO4SM,A1PO4SM,A2PO4SM,
     *    B11PO4SM,B22PO4SM)
        SMD1PO4(L) = SMFD1PO4*RSM1PO4
        WQBFPO4D(L) = SMSS(L) * (SMD1PO4(L) - WQPO4D(L,1))* BFP*
     *   FLXPND(L,1)
        if(FLXPND(L,1).GT.1.01) WQBFPO4D(L) =max(WQBFPO4D(L),0.001)

4567	format(2i6,999e13.5)
        SM1PO4(L) = RSM1PO4
        SM2PO4(L) = RSM2PO4
      END IF
      END DO
C
C Si
C
      IF (IWQSI.EQ.1) THEN
!        DO ND=1,NDMWQ
!         LF=2+(ND-1)*LDMWQ
!         LL=LF+LDM-1
         DO L=2,LA !LF,LL
         IF(SCBWQ(L).GT.0.5) THEN
          SMDFSI(L) = (WQASCD*WQDFBD(L) + WQDFSI(L) + SMJDSI)
     *      * SMDTOH(ISMZMAP(L))
          WQTT = DTWQ * SMTDSI(ISMT(L)) * (SMSISAT-SMFD2SI*SM2SI(L))
     *      / (SMPSI(L)+SMKMPSI+ 1.E-18)
          SMPSI(L) = (SMPSI(L)+SMDFSI(L)) /
     *               (SMW2DTOH(ISMZMAP(L))+WQTT+ 1.E-18)
          IF (XSMO20(L).LT.SMCO2SI) THEN
            SMP1SI = SMP2SI * SMDP1SI**(XSMO20(L)/(SMCO2SI+ 1.E-18))
           ELSE
            SMP1SI = SMP2SI * SMDP1SI
          END IF
          SMFD1SI = 1.0 / (1.0 + SMM1*SMP1SI)
          SMFP1SI = 1.0 - SMFD1SI
          A1SISM = SMKL12(L)*SMFD1SI + SMW12(L)*SMFP1SI
     *      + SMW2(ISMZMAP(L))
          A2SISM = SMKL12(L)*SMFD2SI + SMW12(L)*SMFP2SI
          A11SISM = SMSS(L)*SMFD1SI + A1SISM
          WQTT = SMTDSI(ISMT(L)) * SMPSI(L) * SMHSED(ISMZMAP(L))
     *      / (SMPSI(L)+SMKMPSI+ 1.E-18)
          SMJ2SI = WQTT * SMSISAT
          A22SISM = A2SISM + SMW2PHODT(ISMZMAP(L)) + WQTT*SMFD2SI
          B11SISM = SMSS(L) * WQSAD(L,1)
          B22SISM = SMHODT(ISMZMAP(L))*SM2SI(L) + SMJ2SI
          CALL SOLVSMBE(RSM1SI,RSM2SI,A11SISM,A22SISM,A1SISM,A2SISM,
     *      B11SISM,B22SISM)
          SMD1SI(L) = SMFD1SI*RSM1SI
          WQBFSAD(L) = SMSS(L) * (SMD1SI(L) - WQSAD(L,1))
          SM1SI(L)  = RSM1SI
          SM2SI(L)  = RSM2SI
         END IF
         END DO
!        END DO
      END IF
C
C Keep track of benthic flux rates here for later save to binary file:
C
      DO L=2,LA
        BFO2sum(L)  = BFO2sum(L)  + WQBFO2(L)
        BFNH4sum(L) = BFNH4sum(L) + WQBFNH4(L)
        BFNO3sum(L) = BFNO3sum(L) + WQBFNO3(L)
        BFPO4sum(L) = BFPO4sum(L) + WQBFPO4D(L)
        BFSADsum(L) = BFSADsum(L) + WQBFSAD(L)
        BFCODsum(L) = BFCODsum(L) + WQBFCOD(L)
      END DO
      nbfcnt = nbfcnt + 1
 !     TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
      TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/86400.0
      IF(NCSTEP.GT.0) TIMTMP=SECDLAST/TCON+TBEGIN 
      timebf = timebf + TIMTMP
C
  401 FORMAT(I8,3I5,E12.3,A36)
C
      RETURN
      END
