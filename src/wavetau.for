c wavetau.for, calculate bottom shear stress due to wind wave & current, jeff ji@wave, 3/2/00
c     significant wave height & period come from SWAN (readswan3.for)
c
      subroutine wavetau
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      dimension curang2(LCM),qbar(LCM),ubout(LCM)
     & ,abmout(LCM),phiout(LCM)
     *  ,HS(LCM),DIR(LCM),TM(LCM)     ! not in efdc.cmn, due to dir is used somewhere else
c---------------------------------------------------------
c     iwcm=0              ! =0 = use wcm.inp made by wcm.for, =1 = call wcm2000
c     iwcm=1
c
      nub= 50             ! based on wcm.for
      nab= 40
      nur= 30
      nzr= 41
      nphicw=  6
      ddeg=90.0/(nphicw-1)
c------
      deg2rad=3.14159/180.0
      rad2deg=180.0/3.14159
c
      L1=135   ! (38,12) A, hardwired, Ji, 10/16/00
      L2=1686  ! (26,46) B
      L3=857   ! (39,29) C
      L4=677   ! (23,26) E
c     L1=127 -3  !  27   12  '       '   ! 5, check 1
c     L2=610 -3  !   7   25  'str    '   ! 6        2
c     L3=941 -3  !  10   31  'str    '   ! 7        3
c     L4=1425-3  !  11   40  'str    '   ! 8        4
c     L4=857   ! (39,29) C
c
      TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
c     write(301,391) N,DT,TCON,TBEGIN,TCTMSR, TIME
391   format(i4,999f12.4)
c update SWAN data
c
699   continue
c     write(301,*) "days1,days2,time  ", days1,days2,time
      if(time.gt.days2) then   ! #1
      read(480,692) dayswan,NX
692   format(f14.4,4x,i9)
c      write(6,*) "hsdir.inp N= ",NX,dayswan
      read(480,691,IOSTAT=ISO) (ihs(L),L=1,LSWAN),(IDIR(L),L=1,LSWAN),
     +   (ITM(L),L=1,LSWAN)
691   format(20i5)
      IF(ISO.GT.0) then
      write(6,*) "hsdir.inp READ error: ", dayswan,NX
      stop
      endif
c update wind wave
      days1=days2
      do L=2,LA
      HS1(L)=HS2(L)
      DIR1(L)=DIR2(L)
      TM1(L)=TM2(L)
      ENDDO
      days2=dayswan
      do L=1,LSWAN
      L2=L+1                 ! shift back to original EFDC L sequence
      HS2(L2) =IHS(L)*0.001  ! -> M
      DIR2(L2)=Idir(L)*0.1   ! -> deg
      if(idir(L).eq.-999) dir2(L2)=-999.9
      TM2(L2) =ITM(L)*0.001  ! -> second
      if(ITM(L).eq.-999) TM2(L2)=-999.9
      ENDDO
c
      go to 699
      endif                    ! #1

c interpolate
      TDIFF=days2-days1
      WTM1=(days2-TIME)/TDIFF
      WTM2=(TIME-days1)/TDIFF
C
      do L=2,LA
      HS(L) =WTM1*HS1(L)+ WTM2*HS2(L)
c     DIR(L)=WTM1*DIR1(L)+WTM2*DIR2(L)   ! direction can not be interpolated this way
      if(dir1(L).lt.-900.0.or.dir2(L).lt.-900.0) then
      dir(L)=-999.9
      else
      xrad1=dir1(L)*deg2rad
      xrad2=dir2(L)*deg2rad
      xu1=sin(xrad1)           ! unit vector, w.r.t North, cw
      xv1=cos(xrad1)
      xu2=sin(xrad2)           ! unit vector
      xv2=cos(xrad2)
      xu=WTM1*xu1+WTM2*xu2
      xv=WTM1*xv1+WTM2*xv2
      xrad=atan2(xu,(xv+1.0e-12))
      dir(L)=xrad*rad2deg
      if(dir(L).lt.0.0) dir(L)=dir(L)+360.0
      endif
      TM(L) =WTM1*TM1(L)+ WTM2*TM2(L)
      if(TM1(L).le.0.0 .or.TM2(L).le.0.0 ) TM(L)=-999.9
      ENDDO
c     write(311,381) N,time,WTM1,WTM2,HS1(L1),HS2(L1),HS(L1),
c    +                      dir1(L1),dir2(L1),dir(L1),
c    +                      tm1(L1),tm2(L1),tm(L1)
381   format(i4,f12.4,999f9.3)
c-----------------------------------------------------
C--- modified from s.wavetau.for ------
c     SUBROUTINE WAVETAU
C
C  CALC. BOTTOM SHEAR STRESS DUE TO WIND WAVES AND CURRENTS
C
C  USES GRANT/MADSEN/GLENN WAVE MODEL
C
C  REVISION DATE:  MAY 13, 1996
c
c Modified by Jeff Ji,  23/2/2000
C
C**************************************************************
C
C
      PI=3.14159
      PI2=2.*PI
      PIHALF=PI/2.
      GRAV=9.806
      EPS=1.0E-4      ! for Newton's method, test its sensitivity!
c
c=============================================================
c
c  This section calulate significant wave height, and is not used here
c  Significant wave height and wave direction come from SWAM model
c  in this application.
c  Jeff Ji, 3/2/2000
c
C  CALM WIND ( < 1 MPH = 0.45 m/s), PURE CURRENTS
c    !!Ji, We migtht need to deal with the case of small hs, 3/5/00
c
c     IF (WINDSP.LT.0.45) GOTO 99
C
C  CALC. SIGNIFICANT WAVE HEIGHT AND PERIOD
C  USE SHALLOW WATER SMB THEORY
C
C  HSIG = SIG. WAVE HEIGHT (m)
C  TSIG = SIG. WAVE PERIOD (s)
C
C  WINDSP = WIND SPEED (m/s)
C  DIR(L) = WIND DIRECTION (degrees)
C  HMEAN(N) = MEAN DEPTH ALONG FETCH FOR DIRECTION n (m)
C  FETCH(N) = FETCH FOR DIRECTION n (m)
C
c     NDIR=NINT(DIR(L)/10.)
C
c     IF (NDIR.EQ.0) NDIR=36
C
c     WRAT0=GRAV/(WINDSP*WINDSP)
C
c     HRAT0=0.283/WRAT0
C
c     TRAT0=2.4*PI*WINDSP/GRAV
C
c     DO 10 J=2,JM-1
c       DO 10 I=2,IM-2
c         IF (FSM(L,1).GT.0.0) THEN
c           WRAT1=WRAT0*HMEAN(NDIR,I,J)
c           WRAT2=WRAT0*FETCH(NDIR,I,J)
c
c           HA1=TANH(0.53*(WRAT1**0.75))
c           TA1=TANH(0.833*(WRAT1**0.375))
C
c           hs(L)=HRAT0*HA1*TANH(0.0125*(WRAT2**0.42)/HA1)
C
C  MODIFIED BY C.K.Z. ON 5/13/96
C
C  LIMIT WAVE HEIGHT TO BREAKING WAVE HEIGHT (APPROX.)
C  REF:  USCOE SHORE PROTECTION MANUAL, p. 2-121
C
c           HBREAK=0.78*HP(L)
c           hs(L)=AMIN1(HBREAK,hs(L))
C
c           tm(L)=TRAT0*TA1*TANH(0.077*(WRAT2**0.25)/TA1)
c         ENDIF
c10   CONTINUE
C
C  CALC. bottom CURRENT VELOCITY AND ANGLE
C
      DO L=2,LA
      Umean=0.5*STCUV(L)*(U(L+1,1)+U(L,1))
      Vmean=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
C  CALC. CURRENT ANGLE W.R.T. E1 AXIS (ccw)
C  (CURANG IS IN radians)
C
            IF (UMEAN.NE.0.0) THEN
c              CURANG2(L)=ATAN(VMEAN/UMEAN)   ! SIZ bug, atan-> (-90,90), not (1080,180), Ji, 3/7/00
              CURANG2(L)=ATAN2(VMEAN,UMEAN)
            ELSE
              IF (VMEAN.GT.0.0) THEN
                CURANG2(L)=PI/2.
              ELSE
                IF (VMEAN.EQ.0.0) THEN
                  CURANG2(L)=0.0
                ELSE
                  CURANG2(L)=1.5*PI
                ENDIF
              ENDIF
            ENDIF
C  CONVERT CURANG SO THAT IT IS W.R.T. (with respect to) EAST AXIS (ccw) (counter-clock wise)
C
            CURANG2(L)=CURANG2(L)+ANGleM(L)
C
C  CONVERT CURANG SO THAT IT IS W.R.T. NORTH AXIS (cw) (clock wise)
C
            CURANG2(L)=PIHALF-CURANG2(L)
c
            if(curang2(L).lt.0) curang2(L)=curang2(L)+pi2
c
      qbar(L)=SQRT(umean*umean+vmean*vmean)
c
c     if(L.eq.L1) then
c     write(302,371) N,time,U(L,1),U(L+1,1),V(L,1),V(LNC(L),1),
c    + umean,vmean,qbar(L),curang2(L)*rad2deg
371   format(i5,f9.4,999f8.3)
c     endif
      enddo
c
C -----------------------------------------------------------------
C  CALC. TOTAL BOTTOM SHEAR STRESS
C
      TAU(1)=0.0
      TAU(LC)=0.0
      DO 30 L=2,LA
C
C  CALC. WAVE NUMBER (1/m)
C  ASSUME OMEGAR=OMEGAA
C
C  CALC. SHALLOW WATER WAVE SPEED (m/s)
C
            C0=SQRT(GRAV*HP(L))
C
            IF (tm(L).le.0.0.or.dir(L).lt.-900.0) THEN
              UB=0.00001
              ABM=0.00001
              PHIBAR=0.0
              GOTO 110
            ENDIF
C
C  ESTIMATE WAVE NO.
C
            WAVEK0=PI2/(C0*tm(L))
C
C  ITERATE FOR WAVE NO.
C
            TRAT=PI2*PI2/(GRAV*tm(L)*tm(L))
C
c----------------------------------------------------------------------
c---- copied from DV Fortran exmaple on FUNCTION
c The following example uses the Newton-Raphson iteration method (F(X) = cosh(X) + cos(X) - A = 0) to get the root of the function:
c
c FUNCTION ROOT(A)
c   X  = 1.0
c   DO
c     EX = EXP(X)
c     EMINX = 1./EX
c     ROOT  = X - ((EX+EMINX)*.5+COS(X)-A)/((EX-EMINX)*.5-SIN(X))
c     IF (ABS((X-ROOT)/ROOT) .LT. 1E-6) RETURN
c     X  = ROOT
c   END DO
c END
c-----------------------------------------------
C  USE NEWTON'S METHOD
C
c           DO 40 NX=1,20
            DO 40 NX=1,100
              TANHKH=TANH(WAVEK0*HP(L))
C
              FKH=WAVEK0*TANHKH-TRAT
C
              FPKH=TANHKH+WAVEK0*HP(L)*(1.-TANHKH*TANHKH)
C
              WAVEK=WAVEK0-FKH/FPKH
C
C  CHECK FOR CONVERGENCE
C
              DELTA=ABS(WAVEK/WAVEK0-1.)
              IF (DELTA.LE.EPS) GOTO 45 ! test sensitivity !!
              WAVEK0=WAVEK
 40         CONTINUE
            if(delta.gt.eps) then
            write(6,*) "Newton, wavek", NX, wavek,delta
            stop
            endif
C
C  CALC. BOTTOM ORBITAL AMPLITUDE (m)
C
 45         ABM=hs(L)/(2.*SINH(WAVEK*HP(L))) + 0.00001
c
c           xx=sinh(wavek*HP(L))
c           xwaveL=pi2/wavek
c           if(L.eq.L1) then
c           write(303,361) N,time,hs(L),abm,wavek,xwaveL,xx
361   format(i5,999f9.4)
c           endif
C
C  CALC. WAVE FREQUENCY (1/s)
C
            OMEGAA=PI2/tm(L)
C
C  CALC. BOTTOM ORBITAL VELOCITY (m/s)
C
            UB=ABM*OMEGAA + 0.00001
C
C  CALC. ANGLE BETWEEN WAVE AND CURRENT  (0 < PHIBAR < PI/2 radians)
C
            PHIBAR=ABS((DIR(L))*pi/180.0-CURANG2(L))
            PHIBAR=XMOD(PHIBAR,PIHALF)
C
C
C  CALC. SHEAR VELOCITY DUE TO WAVES & CURRENTS
C
C  NEAR BOTTOM CURRENT VELOCITY
C
 110        UR=QBAR(L)
             ubout(L)=ub
             abmout(L)=abm
             phiout(L)=180.*phibar/pi
C
C  HEIGHT OF VELOCITY ABOVE BOTTOM
C
            ZR=HP(L)*DZC(1)/2.
C
C  BOTTOM ROUGHNESS
C
C  MODIFIED BY C.K.Z. ON 3/22/96
C
C  PER JIM LEWIS, Z0 TO WAVE MODEL (KB IN WCM93) IS 30 * Z0B
C
c           IF (IBMSK(L).EQ.0) THEN
c             Z0=30.*Z0BCOH
c           ELSE
c     z0b=amax1(ZBR(L),0.0025)    ! bottom roughness, check this!
c             Z0=30.*Z0B
c             Z0= 1.*Z0B  ! test
c           ENDIF
c-----
c     if(UB.ge.0.0  ) then ! always use wcm2000
c     if(UB.gt.0.001) then ! when orbital velocity is large enough
c     if(UB.gt.0.010) then ! when orbital velocity is large enough, wcm93 has problem when ZBRADJ is large & water is shallow, alog(1)=large
c     if(UB.gt.0.0001) then ! when orbital velocity is large enough  ! tests indicate both give similar results, jeff, 12/15/00
C
c           CALL WCM93(UB,ABM,UR,PHIBAR,ZR,Z0,USC,USCW,fcw)
C-- Input parameters for wcm2000 ----------------------------
c     sg_ub=40
c     sg_ab=80
c     sg_ur=60
c     sg_zr=100
c     sg_phicw=0
      sg_ub=ub*100.0      ! m/s -> cm/s
      sg_ab=abm*100.0     ! m   -> cm
      sg_ur=ur*100.0      ! m/s -> cm/s
      sg_zr=zr*100.0      ! m   -> cm
      sg_phicw=phibar     ! radian
c     if(L.eq.L3) write(603,1401) N,time,sg_ub,sg_ab,sg_ur,
c    +    sg_zr,sg_phicw
c---------
      if(iwcm.eq.0) then                 ! #3, based on wcm.for
      iub=sg_ub+0.5
      iub=max(1,iub)
      N1 =min(nub,iub)
c
      iab=sg_ab+0.5
      iab=max(1,iab)
      N2 =min(nab,iab)
c
      iur=sg_ur+0.5
      iur=max(1,iur)
      N3 =min(nur,iur)
c
      izr=sg_zr-13.5
      izr=max(1,izr)
      N4 =min(nzr,izr)
c
      iphicw=sg_phicw*rad2deg/ddeg+1.5
      iphicw=max(1,iphicw)
      N5    =min(nphicw,iphicw)
c
      Irec=n5+Nphicw*(N4-1)+Nphicw*NZR*(N3-1)+nphicw*NZR*NUR*(N2-1)+
     *   nphicw*NZR*NUR*NAB*(N1-1)
c
c     write(701,751) Irec,N1,N2,N3,N4,N5,L,N
751   format(999i8)
      read(481,rec=Irec) xtau ! in dynes/cm**2
      tau(L)=xtau*1.0e-4     ! dynes/cm**2 -> (m/s)**2 used in EFDC
c
      else                               ! #3
            call wcm2000(sg_ub,sg_ab,sg_ur,sg_zr,sg_phicw,     ! Input
     +                      ustarc,   ustarwm,   ustarcw)      ! Output
      uscw=ustarcw*0.01   ! cm/s -> m/s
c     if(L.eq.L3) write(503,1401) N,time,ustarc,   ustarwm,   ustarcw
c     if(L.eq.L4) write(504,1401) N,time,ustarc,   ustarwm,   ustarcw
c
      TAU(L)=USCW*USCW      ! in (m/s)**2, , I figured it out on 12/14/00
      endif                              ! #3
c
 30   CONTINUE
C
C  FOR WIND WAVE OUTPUT
C
      NXX=NTSPTC/24  ! hourly
c     NOUT=MOD(N,1)
      NOUT=MOD(N,NXX)      ! save hourly
      IF (NOUT.EQ.0) THEN
c look at A, B, C, E
      write(401,1401) N,time,HP(L1),hs(L1),dir(L1),tm(L1),
     +   qbar(L1),phiout(L1)    ,ubout(L1),tau(L1)*1.0E4                !  uscw**2 (m/s)**2 -> real stress in dynes/cm**2
      write(402,1401) N,time,HP(L2),hs(L2),dir(L2),tm(L2),
     +   qbar(L2),phiout(L2)    ,ubout(L2),tau(L2)*1.0E4
      write(403,1401) N,time,HP(L3),hs(L3),dir(L3),tm(L3),
     +   qbar(L3),phiout(L3)    ,ubout(L3),tau(L3)*1.0E4
      write(404,1401) N,time,HP(L4),hs(L4),dir(L4),tm(L4),
     +   qbar(L4),phiout(L4)    ,ubout(L4),tau(L4)*1.0E4
1401  FORMAT (I8,f9.4,7f10.4,f10.3)
      ENDIF
c
C
      RETURN
      END
C
C*****************************************************************
C
      FUNCTION XMOD(A,B)
C
C  MODIFICATION OF 'AMOD' FUNCTION IN FORTRAN 77
C
C   AMOD(A,B) = A - (INT(A/B))*B
C   (THIS CAN GIVE WRONG ANSWERS FOR THIS SUBROUTINE)
C
C  REVISION DATE:  SEPT. 28, 1995
C
C****************************************************************
C

      A=A-(NINT(A/B))*B
C
c     A=ABS(A) ! bug, xmod has no value !, Ji, 3/5/00
      xmod=abs(A) ! Ji
C
      RETURN
      END
