C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WWQTSbin(LA,KC,DT,N,TCON,TBEGIN,TIDALP,NTSPTC,
     &  NCSTEP,SECDLAST)
c
c   This sub might not be needed, since we have GrADS savings 
c   using dump3.for, jeff ji, 8/26/02
c 
C
C**********************************************************************C
C Write time-series output: WQCHLx=1/WQCHLx to BINARY file
C**********************************************************************C
C M. Morton  28 Jun 1998
C Averages WQ variables over IWQTSDT time steps (e.g., daily averages).
C Also, determines maximum and minimum over IWQTSDT time steps.
C Averaged diurnal DO variables over IWQDIUDT time steps
C
C Revisions:
C  23 Jul 98 MRM: added binary file write to dump all cells to output
C     file.  See SUBROUTINE INITbin for header record information.
C**********************************************************************C
C **  SHEN'S MODIFICATION TO OUTPUT MACROALGAE
C**********************************************************************C
C
c These are the index numbers of the parameters in WQVo array:
c WQVo(LL,K, 1) = cyanobacteria as C  WQVo(LL,K,12) = labile PON
c WQVo(LL,K, 2) = diatoms as C        WQVo(LL,K,13) = diss. org. nitrogen
c WQVo(LL,K, 3) = green algae as C    WQVo(LL,K,14) = ammonia nitrogen
c WQVo(LL,K, 4) = refractory POC      WQVo(LL,K,15) = nitrate nitrogen
c WQVo(LL,K, 5) = labile POC          WQVo(LL,K,16) = part. biogenic silica
c WQVo(LL,K, 6) = diss. org. carbon   WQVo(LL,K,17) = available silica
c WQVo(LL,K, 7) = refractory POP      WQVo(LL,K,18) = COD
c WQVo(LL,K, 8) = labile POP          WQVo(LL,K,19) = dissolved oxygen
c WQVo(LL,K, 9) = diss. org. phos.    WQVo(LL,K,20) = total active metal
c WQVo(LL,K,10) = tot. inorg. phos.   WQVo(LL,K,21) = fecal coliform bacteria
c WQVo(LL,K,11) = refractory PON      WQVo(LL,K,22) = macroalgae
c
      INCLUDE 'wq.par' 
      INCLUDE 'wqcom.cmn'
      REAL*8 SECDLAST 
      TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/86400
      IF(NCSTEP.GT.0) TIMTMP=SECDLAST/TCON+TBEGIN 
          OPEN(UNIT=82, FILE='WQDOCOMP.bin',ACCESS='append',
     +       FORM='unformatted',STATUS='unknown')
          write(82)TIMTMP
          DO LL=2,LA
            DO K=1,KC
              xLimNc(LL,K) = xLimNc(LL,K) / NLIM
              xLimNd(LL,K) = xLimNd(LL,K) / NLIM
              xLimNg(LL,K) = xLimNg(LL,K) / NLIM
              xLimNm(LL,K) = xLimNm(LL,K) / NLIM
              xLimPc(LL,K) = xLimPc(LL,K) / NLIM
              xLimPd(LL,K) = xLimPd(LL,K) / NLIM
              xLimPg(LL,K) = xLimPg(LL,K) / NLIM
              xLimPm(LL,K) = xLimPm(LL,K) / NLIM
              xLimIc(LL,K) = xLimIc(LL,K) / NLIM
              xLimId(LL,K) = xLimId(LL,K) / NLIM
              xLimIg(LL,K) = xLimIg(LL,K) / NLIM
              xLimIm(LL,K) = xLimIm(LL,K) / NLIM
              xLimTc(LL,K) = xLimTc(LL,K) / NLIM
              xLimTd(LL,K) = xLimTd(LL,K) / NLIM
              xLimTg(LL,K) = xLimTg(LL,K) / NLIM
              xLimTm(LL,K) = xLimTm(LL,K) / NLIM
              xDOdz(LL,K) = xDOdz(LL,K) / NLIM

              xmrm = 1.0 / xDOdz(LL,K)
              xDOsat(LL,K) = xDOsat(LL,K) *xmrm
              xDOdef(LL,K) = xDOdef(LL,K) *xmrm
              xDOall(LL,K) = xDOall(LL,K) *xmrm
              xnumer = xDOsat(LL,K) - xDOdef(LL,K)
              xnumer = max (xnumer, 0.0)
              xdenom = abs (xDOall(LL,K))
c              xdenom = max (xdenom, 0.1)
c              xratio = xnumer / (DTWQ*(1.0e-9+xdenom))
              xratio = 1.0 
              xDOdz(LL,K) = abs (xDOall(LL,K)) *xratio

              xDOpsl(LL,K) = xDOpsl(LL,K) *xmrm *xratio
              xDOsod(LL,K) = xDOsod(LL,K) *xmrm *xratio
              xDOkar(LL,K) = xDOkar(LL,K) *xmrm *xratio
              xDOdoc(LL,K) = xDOdoc(LL,K) *xmrm *xratio
              xDOnit(LL,K) = xDOnit(LL,K) *xmrm *xratio
              xDOcod(LL,K) = xDOcod(LL,K) *xmrm *xratio
              xDOppB(LL,K) = xDOppB(LL,K) *xmrm *xratio
              xDOrrB(LL,K) = xDOrrB(LL,K) *xmrm *xratio
              xDOppM(LL,K) = xDOppM(LL,K) *xmrm *xratio
              xDOrrM(LL,K) = xDOrrM(LL,K) *xmrm *xratio
              xDOtrn(LL,K) = xDOtrn(LL,K) *xmrm *xratio
              xDOacoef(LL,K)=xDOacoef(LL,K)/ NLIM
            END DO
          END DO
c
c write to BINARY file WQDOCOMP.bin here:
c
           goto 800    ! Temperol shut down 2010 for Paul
           
            WRITE(82) xLimNc
			WRITE(82) xLimNd
			WRITE(82) xLimNg
            WRITE(82) xLimNm
			WRITE(82) xLimPc
			WRITE(82) xLimPd
			WRITE(82) xLimPg
			WRITE(82) xLimPm
			WRITE(82) xLimIc
			WRITE(82) xLimId
			WRITE(82) xLimIg
			WRITE(82) xLimIm
			WRITE(82) xLimTc
			WRITE(82) xLimTd
			WRITE(82) xLimTg
			WRITE(82) xLimTm
			WRITE(82) xDOsat
			WRITE(82) xDOpsl
			WRITE(82) xDOsod
			WRITE(82) xDOkar
			WRITE(82) xDOdoc
			WRITE(82) xDOnit
			WRITE(82) xDOcod
			WRITE(82) xDOppB
			WRITE(82) xDOrrB
			WRITE(82) xDOppM
			WRITE(82) xDOrrM
			WRITE(82) xDOdef
			WRITE(82) xDOtrn
			WRITE(82) xDOall
			WRITE(82) xDOdz 
	      CALL FLUSH(82) 
	      
  800       WRITE(82)xDOacoef      ! Output for Paul
            WRITE(82)xDOkar
 	      CALL FLUSH(82) 
            CLOSE(82)
          
           OPEN(UNIT=83, FILE='WQDOCOMP.out',ACCESS='append',
     +       STATUS='unknown')
          DO M=1,IWQTS                              
             LL=LWQTS(M)
           DO k=1,KC
           write(83,820)ILW(LL),JLW(LL),K,TIMTMP
     +      ,xLimNc(LL,k),xLimNd(LL,k),xLimNg(LL,k)
     +		,xLimNm(LL,k),xLimPc(LL,k),xLimPd(LL,k),xLimPg(LL,k)
     +      ,xLimPm(LL,k),xLimIc(LL,k),xLimId(LL,k),xLimIg(LL,k)
     +	    ,xLimIm(LL,k),xLimTc(LL,k),xLimTd(LL,k),xLimTg(LL,k)
     +      ,xLimTm(LL,k),xDOsat(LL,k),xDOpsl(LL,k),xDOsod(LL,k)
     +      ,xDOkar(LL,k),xDOdoc(LL,k),xDOnit(LL,k),xDOcod(LL,k)
     +      ,xDOppB(LL,k),xDOrrB(LL,k),xDOppM(LL,k),xDOrrM(LL,k)
     +      ,xDOdef(LL,k),xDOtrn(LL,k),xDOall(LL,k),xDOdz(LL,k)
     +      ,SEDTWQ(LL,K),TEMWQ(LL,K),UWQS(LL),VWQS(LL),SALWQ(LL,K)
     +      ,xKe(LL)    
           ENDDO
         ENDDO
          CLOSE(83)
          call WQzero3(LA,KC)
 820     format(3I6,44F9.3)
c
      RETURN
      END
