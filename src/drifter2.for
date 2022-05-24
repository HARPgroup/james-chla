C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE DRIFTER2
C
c a) Particle tracking subroutine based on Euler method, i.e., forward marching.
c Since the tracking calculation
c is done with hydrodynamic calculation simultanously, dt should be very small
c already and other heigher order accuracy methods, such 4th Runge-Kutta,
c might not be neceesary.
c b) The idea is to track particles based on their (I,J) locations directly.
c It makes finding (u,v) much easier.
c c) Go to efdc.inp, Card 67 - 69 for parameter input definitions
c d) See notes for more details
c
C Jeff Ji, 4/1/2001
C
C**********************************************************************C
C
C **  SUBROUTINE DRIFTER CALCULATES THREE DIMENSIONAL TRAJECTORIES
C **  OF NEUTRALLY BUOYANT PARTICLES
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      TIME=DT*FLOAT(N)/tcon+TBEGIN        ! in days
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
c
      IF(JSPD.NE.1) GO TO 5
C
      JSPD=0                              ! first time
      OPEN(96,FILE='drifter2.out',STATUS='UNKNOWN')
      CLOSE(96,STATUS='DELETE')
      OPEN(96,FILE='drifter2.out',STATUS='UNKNOWN')
C
      write(96,209) (NP,NP=1,NPD)
209   format(8x,999(1x,i12))
c     WRITE(96,207)N,time,(RI(NP),RJ(NP),NP=1,NPD)
      WRITE(96,208) time,(RI(NP),RJ(NP),NP=1,NPD)
208   format(f8.3,999(1x,2f6.2))
207   format(i8,f12.4,999(1x,2f6.2))
      CLOSE(96)
C
    5 CONTINUE
C----------------------------------------------------------------------C
C
C **  LOOP OVER THE NUMBER OF PARTICLES
C
      kx=kc          ! surface layer
c     kx=kc/2.0+0.6  ! middle  layer
c     kx=1           ! bottom  layer
      DO 100 NP=1,NPD
c
      i=RI(NP)
      j=RJ(NP)
      xi=RI(NP)-i
      xj=RJ(NP)-j
      L=LIJ(I,J)
      LN=LNC(L)
c
c I-direction,
      xu=u1(L,Kx)*(1.0-xi)+u1(L+1,Kx)*xi
c     xu=0.10                               ! put theoretical values here for testing
      RI1   =RI(NP)+dt*xu/dxp(L)            ! Euler method & convert to (I,J)
c J-direction
      xv=v1(L,Kx)*(1.0-xj)+v1(LN ,Kx)*xj
c     xv=0.05                               ! put theoretical values here for testing
      RJ1   =RJ(NP)+dt*xV/dyp(L)            ! Euler method & convert to (I,J)
c
c check if it is still in the domain
      i=RI1
      j=RJ1
      if(LIJ(I,J).ne.1) then ! based on cellmap.for, it is not land
      RI(NP)=RI1
      RJ(NP)=RJ1
      else
      write(772,905) NP,N,time,RI1,RJ1      ! if (RI1,RJ1) is land, keep the old (RI,RJ)
905   format(2i8,9f12.4)
      endif

100   continue
C
C**********************************************************************C
C
      IF (NCPD.EQ.NWPD) THEN
C
      NCPD=1
      OPEN(96,FILE='drifter2.out',ACCESS='APPEND',STATUS='UNKNOWN')
c     WRITE(96,207)N,time,(RI(NP),RJ(NP),NP=1,NPD)
      WRITE(96,208) time,(RI(NP),RJ(NP),NP=1,NPD)
      CLOSE(96)
c
      ELSE
      NCPD=NCPD+1
      END IF
C
C**********************************************************************C
C
      RETURN
      END
