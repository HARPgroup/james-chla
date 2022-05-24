!**********************************************************************C
! 
      SUBROUTINE RWQPSL(DT,N,TBEGIN,TCON,LA,KC,IWQS,NTSPTC,
     &  NCSTEP,SECDLAST)
!
!**********************************************************************C
!
! **  NEW VERSION BY J. M. HAMRICK  7 APRIL 1997
!
! **  LAST MODIFIED BY JOHN HAMRICK ON  7 APRIL 1997
!     
! **  Last modified by AEE on 3/21/2007
!
!**********************************************************************C
!
! **  ALL LOADS ARE IN KG/DAY EXPECT COLIFORM IN MPN/DAY
! **  FCB: MPN/L*10*(xPSQ m^3/s)*(86400 s/d) -> output will be MPN/100mL
! **  INTERNAL CONVERSION TO GM/DAY FOR FIRST 20 STATE VARIABLES
! **  TAM is not used
!
!**********************************************************************C
!
      PARAMETER (NACP=12,NACP1=9)
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
      REAL*8 SECDLAST 
      DIMENSION RLDTMP(21),IACP(NACP),IACP1(NACP1)
      DATA   IACP/3,5,6,8,9,10,12,13,14,15,18,19/
      DATA   IACP1/3,6,9,10,13,14,15,18,19/
C
C**********************************************************************C
C
      IF(ITNWQ.GT.0) GO TO 1000
C
C**********************************************************************C
C
C **  READ IN LOADING SERIES FROM FILE 'wqpsl.inp'
C
C----------------------------------------------------------------------C
C
      IF( NPSTMSR.GE.1) THEN
        if(IWQS.LE.0) then
         OPEN(1,FILE='wqpsl.inp',STATUS='UNKNOWN')
        elseif(IWQS.EQ.1) then
         OPEN(1,FILE='wqpsl12.inp',STATUS='UNKNOWN') 
        elseif(IWQS.EQ.2)then
         OPEN(1,FILE='wqpsl9.inp',STATUS='UNKNOWN')        
        endif      
       DO NW=1,20
         RLDTMP(NW)=0.0 
       END DO
C
C **  SKIP OVER TITLE AND HEADER LINES
C
       DO IS=1,13
         READ(1,1)
       END DO
C
       DO NS=1,NPSTMSR
        MWQPTLT(NS)=1
        READ(1,*,IOSTAT=ISO)MWQPSR(NS),TCWQPSR(NS),
     $                  TAWQPSR(NS),RMULADJ,ADDADJ
        IF(ISO.GT.0) GO TO 900
c
        READ(1,*,IOSTAT=ISO) (WKWQ(NS,K),K=1,KC) ! similar to qser.inp, Ji, 9/16/99
        IF(ISO.GT.0) GO TO 900
c 
        if(MWQPSR(NS).gt.NDWQPSR) then   ! check array dimensions
         write(6,*) " Error:  MWQPSR(NS)>NDWQPSR ", MWQPSR(NS),NDWQPSR
         stop
        endif
c
         RMULADJ=1000.*RMULADJ   ! kg -> g, except FC
         ADDADJ=ADDADJ
        if(IWQS.LE.0)THEN
         DO M=1,MWQPSR(NS)
          READ(1,*,IOSTAT=ISO)TWQPSER(M,NS),(RLDTMP(NW),NW=1,7)
          IF(ISO.GT.0) GO TO 9001
          READ(1,*,IOSTAT=ISO)(RLDTMP(NW),NW=8,14)
          IF(ISO.GT.0) GO TO 9002
          READ(1,*,IOSTAT=ISO)(RLDTMP(NW),NW=15,21)
          IF(ISO.GT.0) GO TO 900
          DO NW=11,15
          RLDTMP(NW)=RLDTMP(NW)*(1-RN_R)
          ENDDO
          DO NW=4,6
          RLDTMP(NW)=RLDTMP(NW)*(1-RC_R)
          ENDDO
          DO NW=7,10
          RLDTMP(NW)=RLDTMP(NW)*(1-RP_R)
          ENDDO
          TWQPSER(M,NS)=TWQPSER(M,NS)+TAWQPSR(NS)
          DO NW=1,20
            WQPSSER(M,NW,NS)=RMULADJ*RLDTMP(NW)    ! no unit convert for FC 
          END DO
         END DO    
        elseif(IWQS.EQ.1)THEN
          DO M=1,MWQPSR(NS)
          READ(1,*,IOSTAT=ISO)TWQPSER(M,NS),(RLDTMP(IACP(NW)),NW=1,NACP)
          IF(ISO.GT.0) GO TO 9001
          DO NW=11,15
          RLDTMP(NW)=RLDTMP(NW)*(1-RN_R/100.0)
          ENDDO
          DO NW=4,6
          RLDTMP(NW)=RLDTMP(NW)*(1-RC_R/100.0)
          ENDDO
          DO NW=7,10
          RLDTMP(NW)=RLDTMP(NW)*(1-RP_R/100.0)
          ENDDO
           TWQPSER(M,NS)=TWQPSER(M,NS)+TAWQPSR(NS)
          DO NW=1,20
           WQPSSER(M,NW,NS)=RMULADJ*RLDTMP(NW)     
          END DO
         END DO       
        else
          DO M=1,MWQPSR(NS)
          READ(1,*,IOSTAT=ISO)TWQPSER(M,NS),
     *         (RLDTMP(IACP1(NW)),NW=1,NACP1)
          IF(ISO.GT.0) GO TO 9001
          DO NW=11,15
          RLDTMP(NW)=RLDTMP(NW)*(1-RN_R/100.0)
          ENDDO
          DO NW=4,6
          RLDTMP(NW)=RLDTMP(NW)*(1-RC_R/100.0)
          ENDDO
          DO NW=7,10
          RLDTMP(NW)=RLDTMP(NW)*(1-RP_R/100.0)
          ENDDO
          TWQPSER(M,NS)=TWQPSER(M,NS)+TAWQPSR(NS)
          DO NW=1,20
            WQPSSER(M,NW,NS)=RMULADJ*RLDTMP(NW)    
          END DO
         END DO    
        endif          
        
        imov=0
!   moving avearge    J.S. 12/24/2010
        IF(imov.gt.0) THEN
         DO NW=1,20         
          DO M=MWQPSR(NS),imov,-1
           a_tmp=0.
           DO J=M,M-imov+1,-1
            a_tmp=a_tmp+WQPSSER(M,NW,NS)
           ENDDO
            a_tmp=a_tmp/real(imov)
            WQPSSER(M,NW,NS)=a_tmp
          ENDDO                
         ENDDO
        ENDIF
!                    
       ENDDO
       CLOSE(1)
      END IF
706	format(i4,f12.4,3i4,99f12.4)

      GO TO 901

 9001 CONTINUE
      WRITE(6,*)" wqpsl.inp, 1 ",NS,M
      STOP
C
 9002 CONTINUE
      WRITE(6,*)" wqpsl.inp 2 ",NS,M
      write(6,9561) TWQPSER(M,NS)
9561  format(7e12.4)
      STOP
C
  900 CONTINUE
      WRITE(6,*)" wqpsl.inp 3 ",NS,M
      STOP
C
  901 CONTINUE
C
    1 FORMAT(120X)
  601 FORMAT(' READ ERROR WQPS TIME SERIES, NSER,MDATA = ',2I5)
  602 FORMAT(' READ OF FILE wqpsl.inp SUCCESSFUL'/)
C
C**********************************************************************C
C
 1000 CONTINUE
C
C**********************************************************************C
C
C **  INITIALIZE NULL SERIES LOADING TO ZERO
C
      DO NW=1,NWQV
       WQPSSRT(NW,0)=0.
      END DO
C
C**********************************************************************C
C
C **  LOADING SERIES INTERPOLTATION
C
    !  write(*,*)'Total nonpint source ',NPSTMSR
      DO NS=1,NPSTMSR

       TIME=DT*FLOAT(N-1)/TCWQPSR(NS)
     &           +TBEGIN*(TCON/TCWQPSR(NS))
     
      IF(NCSTEP.GT.0) TIME=(SECDLAST-DT)/TCON+TBEGIN  !% J.S. 6/16/2014  
C
       M1=MWQPTLT(NS)
  100  CONTINUE
       M2=M1+1
       IF (TIME.GT.TWQPSER(M2,NS)) THEN
         M1=M2
         GO TO 100
        ELSE
         MWQPTLT(NS)=M1
       END IF
C
       TDIFF=TWQPSER(M2,NS)-TWQPSER(M1,NS)
       WTM1=(TWQPSER(M2,NS)-TIME)/TDIFF
       WTM2=(TIME-TWQPSER(M1,NS))/TDIFF
        DO NW=1,NWQV
         WQPSSRT(NW,NS)=WTM1*WQPSSER(M1,NW,NS)
     &                   +WTM2*WQPSSER(M2,NW,NS)
        END DO
C
      END DO
C
C
!      IF(ITNWQ.EQ.0) THEN
C
!       OPEN(1,FILE='wqpslt.dia',STATUS='UNKNOWN')
!       CLOSE(1,STATUS='DELETE')
!       OPEN(1,FILE='wqpslt.dia',STATUS='UNKNOWN')
C
!       WRITE(1,112)N,TIME
C
!       DO NS=1,NPSTMSR
!            WRITE(1,111)NS,(WQPSSRT(NW,NS),NW=1,NWQV)
!       END DO
C
!      CLOSE(1)
C
!      END IF
C
C**********************************************************************C
!
! Initialize combined PSL array to zero:
!
      DO NW=1,NWQV
       DO K=1,KC
        DO L=2,LA
         WQWPSL(L,K,NW) = 0.0
        END DO
       END DO
      END DO
!
!  WQWPSLC:  is for constant point source input. Not supported (AEE, 3/20/2007)
!  WQWPSL:   final point loads go to model
!  WQPSSRT : input point time series
!

!       WRITE(1,112)N,TIME

       IF(ITNWQ.EQ.0) THEN
       OPEN(10,FILE='wqpsl.dia',STATUS='UNKNOWN')
       CLOSE(10,STATUS='DELETE')
       ENDIF

      IF(ITNWQ.LE.15) THEN
       OPEN(10,FILE='wqpsl.dia',ACCESS='APPEND')
       WRITE(10,112)N,TIME
      ENDIF

      DO NS=1,IWQPS
        L = LIJW(ICPSL(NS), JCPSL(NS))
        ITMP = MVPSL(NS)
        do k=1,kc
        DO NW=1,NWQV
	    WQWPSL(L,K,NW) = WQWPSL(L,K,NW)
     +         + WQWPSLC(NS,NW) + WQPSSRT(NW,ITMP)*WkWQ(ITMP,K)
        END DO

         IF(ITNWQ.LE.15) THEN
          WRITE(10,110)ILW(L),JLW(L),K,(WQWPSL(L,K,NW),NW=1,NWQV)
	   ENDIF
        enddo
      END DO
705   format(i4, f12.4,7i4,99f12.4)



!       DO K=1,KC  
!         DO L=2,LA
!          ITMP=IWQPSC(L,1)
!          IF(ITMP.GT.0) THEN
!             WRITE(10,110)ILW(L),JLW(L),K,(WQWPSL(L,K,NW),NW=1,NWQV)
!          END IF
!         END DO
!       END DO

      IF(ITNWQ.LE.15) CLOSE(10)

!      END IF

  110 FORMAT(1X,3I4,2X,7E12.4,/,15X,7E12.4,/,15X,7E12.4)
  111 FORMAT(1X,I4,2X,7E12.4,/,7X,7E12.4,/,7X,7E12.4)
  112 FORMAT(' N, TIME = ', I10, F12.5/)
!
!**********************************************************************C
!
      RETURN
      END
