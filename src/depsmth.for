C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE DEPSMTH
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      IF (WSMH.GT.0.) GO TO 1000
C
C **  SMOOTH BOTTOM ELEVATION
C
      WSMH=-WSMH
C
      DO NSM=1,NSHMAX
C
      DO L=2,LA
       IF(LCT(L).GT.0.AND.LCT(L).LT.9) THEN
        I=IL(L)
        J=JL(L)
        HTN=BELV(LNC(L))
        HTS=BELV(LSC(L))
        IF (IJCT(I,J+1).EQ.9) HTN=BELV(L)
        IF (IJCT(I,J-1).EQ.9) HTS=BELV(L)
        HTMP(L)=(1.-WSMH)*BELV(L)+0.5*WSMH*(HTN+HTS)
       END IF
      END DO
C
      DO L=2,LA
       IF(LCT(L).GT.0.AND.LCT(L).LT.9) THEN
        I=IL(L)
        J=JL(L)
        HTE=HTMP(L+1)
        HTW=HTMP(L-1)
        IF (IJCT(I+1,J).EQ.9) HTE=HTMP(L)
        IF (IJCT(I-1,J).EQ.9) HTW=HTMP(L)
        BELV(L)=(1.-WSMH)*HTMP(L)+0.5*WSMH*(HTE+HTW)
       END IF
      END DO
C
      END DO
C
C **  SMOOTH DEPTH
C
      DO NSM=1,NSHMAX
C
      DO L=2,LA
       IF(LCT(L).GT.0.AND.LCT(L).LT.9) THEN
        I=IL(L)
        J=JL(L)
        HTN=HMP(LNC(L))
        HTS=HMP(LSC(L))
        IF (IJCT(I,J+1).EQ.9) HTN=HMP(L)
        IF (IJCT(I,J-1).EQ.9) HTS=HMP(L)
        HTMP(L)=(1.-WSMH)*HMP(L)+0.5*WSMH*(HTN+HTS)
       END IF
      END DO
C
      DO L=2,LA
       IF(LCT(L).GT.0.AND.LCT(L).LT.9) THEN
        I=IL(L)
        J=JL(L)
        HTE=HTMP(L+1)
        HTW=HTMP(L-1)
        IF (IJCT(I+1,J).EQ.9) HTE=HTMP(L)
        IF (IJCT(I-1,J).EQ.9) HTW=HTMP(L)
        HMP(L)=(1.-WSMH)*HTMP(L)+0.5*WSMH*(HTE+HTW)
       END IF
      END DO
C
      END DO
C
      GO TO 2000
C
C**********************************************************************C
C
 1000 CONTINUE
C
C **  SMOOTH BOTTOM ELEVATION
C
      DO L=2,LA
       HTMP(L)=BELV(L)
      END DO
C
      DO NSM=1,NSHMAX
C
      DO L=2,LA
       IF(LCT(L).GT.0.AND.LCT(L).LT.9) THEN
        I=IL(L)
        J=JL(L)
        HTN=HTMP(LNC(L))
        HTS=HTMP(LSC(L))
        HTE=HTMP(L+1)
        HTW=HTMP(L-1)
        IF (IJCT(I  ,J+1).EQ.9) HTN=HTMP(L)
        IF (IJCT(I  ,J-1).EQ.9) HTS=HTMP(L)
        IF (IJCT(I+1,J  ).EQ.9) HTE=HTMP(L)
        IF (IJCT(I-1,J  ).EQ.9) HTW=HTMP(L)
        BELV(L)=(1.-WSMH)*HTMP(L)+0.25*WSMH*(HTN+HTS+HTE+HTW)
       END IF
      END DO
C
      DO L=2,LA
       HTMP(L)=BELV(L)
      END DO
C
      END DO
C
C **  SMOOTH DEPTH
C
      DO L=2,LA
       HTMP(L)=HMP(L)
      END DO
C
      DO NSM=1,NSHMAX
C
      DO L=2,LA
       IF(LCT(L).GT.0.AND.LCT(L).LT.9) THEN
        I=IL(L)
        J=JL(L)
        HTN=HTMP(LNC(L))
        HTS=HTMP(LSC(L))
        HTE=HTMP(L+1)
        HTW=HTMP(L-1)
        IF (IJCT(I  ,J+1).EQ.9) HTN=HTMP(L)
        IF (IJCT(I  ,J-1).EQ.9) HTS=HTMP(L)
        IF (IJCT(I+1,J  ).EQ.9) HTE=HTMP(L)
        IF (IJCT(I-1,J  ).EQ.9) HTW=HTMP(L)
        HMP(L)=(1.-WSMH)*HTMP(L)+0.25*WSMH*(HTN+HTS+HTE+HTW)
       END IF
      END DO
C
      DO L=2,LA
       HTMP(L)=HMP(L)
      END DO
C
      END DO 
C
C**********************************************************************C
C
 2000 CONTINUE
C
       OPEN(1,FILE='newdxdy.tmp',STATUS='UNKNOWN')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='newdxdy.tmp',STATUS='UNKNOWN')
       DO L=2,LA
        WRITE(1,339)IL(L),JL(L),DXP(L),DYP(L),HMP(L),BELV(L),
     $              ZBR(L)
       END DO
       CLOSE(1)
c hong=-------------------
!           IF (NQCTL .GT.0 .OR. NQCTLT .GT. 0) THEN 
!	! for add dam, keep wanted waterdepth 
!       do j=2,30  !28
!	 L=LIJ(126,J)
!	HMP(L)=1.0
!	BELV(L)=-1.0
!	enddo
!	    ENDIF
C-HONG-----------------------
  339 FORMAT(1X,I5,2X,I5,2X,E12.4,2X,E12.4,2X,E12.4,2X,E12.4,2X,E12.4)
C
      RETURN
      END
