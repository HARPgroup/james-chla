C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      PROGRAM EFDC
C
C     NON-CRAY VERSION 10 april 1999 BETA RELEASE
C
C**********************************************************************C
C
C **  WELCOME TO THE ENVIRONMENTAL FLUID DYNAMICS COMPUTER CODE SERIES
C **  DEVELOPED BY JOHN M. HAMRICK.  THE EFDC CODE WAS ORGINALLY
C **  DEVELOPED AT VIRGINIA INSTITUTE OF MARINE SCIENCE
C **  /SCHOOL OF MARINE SCIENCE, THE COLLEGE OF
C **  WILLIAM AND MARY, GLOUCESTER POINT, VA 23062
C
C **  THIS SOURCE FILE IS A DIRECT RELEASE BY THE DEVELOPER
C **  AND DIFFERS SIGNIFICANTLY FROM PRE 1 MARCH 96 VIMS RELEASES OF
C **  EFDC AND POST 1 MARCH 96 VIMS RELEASES OF HEM3D (THE NAMED
C **  CURRENTLY USED BY VIMS FOR THE VERSION OF EFDC EXISTING AT
C **  THE TIME OF MY DEPARTURE) WITH RESPECT TO ERROR FIXES AND
C **  APPLICATION CAPABILITIES
C
C **  ENVIRONMENTAL FLUID DYNAMICS CODE AND EFDC ARE
C **  TRADEMARKS OF JOHN M. HAMRICK, PH.D., P.E. 
C
C **  EFDC SOLVES THE 3D REYNOLDS AVERAGED NAVIER-STOKES
C **  EQUATIONS (WITH HYDROSTATIC AND BOUSINESSQ APPROXIMATIONS) AND
C **  TRANSPORT EQUATIONS FOR TURBULENT INTENSITY, TURBULENT
C **  INTENSITYxLENGHT SCALE, SALINITY (OR WATER VAPOR CONTENT),
C **  TEMPERATURE, AN INERT TRACER (CALLED DYE), A DYNAMICALLY ACTIVE
C **  SUSPENDED SETTLING PARTICLE FIELD (CALLED SEDIMENT).  A FREE
C **  SURFACE OR RIGID LID IS PRESENT ON THE VERTICAL BOUNDARY Z=1
C **  IN THE SIGMA STRETCHED VERTICAL COORDINATE.  THE HORIZONTAL
C **  COORDINATE SYSTEM IS CURVILINEAR AND ORTHOGONAL.
C
C **  THE NUMERICAL SOLUTION SCHEME IS ON A SPATIALLY STAGGERED MAC
C **  OR C GRID AND THE TIME INTEGRATION USES A THREE TIME LEVEL
C **  LEAPFROG INTEGRATION WITH PERIODIC TRAPEZOIDAL CORRECTIONS TO
C **  SUPPRESS THE COMPUTATIONAL MODE AND REDUCE NOISE.
C **  SPATIAL SOLUTION OF THE EXTERNAL MODE FOR THE FREE SURFACE
C **  ELEVATION OR KINEMATIC PRESSURE UNDER THE RIGID LID IS BY
C **  RED-BLACK SUCCESSIVE OVER RELAXATION (RB SOR) OR CONJUGATE
C **  GRADIENT SOLUTION OF A PSEUDO-HEMHOLTZ EQUATION.  THE INTERNAL
C **  SOLUTION IS IMPLICIT FOR THE VERTICAL SHEAR OR VELOCITY STRUCTURE.
C **  A NUMBER OF OPTIONS ARE AVAILABLE FOR REPRESENTING THE ADVECTIVE
C **  TRANSPORT TERMS IN THE MOMENTUM AND SCALAR TRANSPORT EQUATIONS.
C
C **  FULL DOCUMENTATION OF THE NUNERICAL SCHEME IS FOUND IN:
C
C **     HAMRICK, J. M., (1992)  'A THREE-DIMENSIONAL ENVIRONMENTAL
C **       FLUID DYNAMICS COMPUTER CODE',  THE COLLEGE OF WILLIAM
C **       MARY, VIRGINIA INSTITUTE OF MARINE SCIENCE, SPECIAL REPORT
C **       IN APPLIED MARINE SCIENCE AND OCEAN ENGINEERING, NO. 317.
C
C **  CHANGES MADE TO THIS CODE BY UNAUTHORIZED PERSONS WILL BE
C **  SUPPORTED ON A COST REIMBURSED BASIS ONLY.  SUPPORT IS AVAILABLE
C **  FROM JOHN M. HAMRICK, 2520 WEST WHITTAKER CLOSE
C **  WILLIAMSBURG, VA, TEL. 804-258-0608, FAX 804-258-9698
C **  EMAIL: ham@visi.net
C
C **  THE AUTHOR ASSUMES NO LIABILITY FOR USE
C **  OF THIS CODE FOR ENVIRONMENTAL AND ENGINEERING STUDIES.
C
C **  THE FOLLOWING FILES ARE NECESSARY TO COMPILE THIS CODE:
C
C **      efdc.com
C **      efdc.par
C
C **  THIS CODE HAS BEEN COMPLIED ON SUN SPARC, HP 9000/700, AND SGI
C **  WORKSTATIONS, VAX VMS AND DEC ALPHA SYSTEMS,
C **  CRAY Y/MP AND C90 SYSTEMS, MACINTOSH SYSTEMS USING LSI AND
C **  ABSOFT FORTRAN, AND 486 AND PENTIUM PC'S USING LAHEY
C **  PROFESSIONAL FORTRAN AND MICROSOFT POWERSTATION
C **  LINES IN THE CODE BEGINNING WITH CDHP IMPLEMENT THE
C **  HP 9000/700 SERIES VECTOR LIBRARY AND MAY BE ACTIVATES BY
C **  REPLACING CDHP WITH 4 BLANK SPACES.  VAX EXTENSION TIMING
C **  UTILITIES USING THE CALL SECNDS MAY BE ACTIVATED BY UNCOMMENTING
C **  ALL LINES CONTAINING SECNDS.
C **  TO RUN ON CRAY SYSTEMS, REPLACE ACCESS='APPEND' WITH
C **  POSITION='APPEND' IN APPROPRIATE FILE OPEN STATEMENTS
C
C **  THE FOLLOWING FILES MAY BE NECESSARY TO RUN THIS CODE:
C
C **      efdc.inp
C **      dxdy.inp
C **      lxly.inp
C **      salt.inp or restart.inp or restran.inp
C **      temp.inp or restart.inp or restran.inp
C **      aser.inp
C **      qser.inp
C **      pser.inp
C **      sser.inp
C **      tser.inp
C **      dser.inp
C **      sdser.inp
C **      snser.inp
C **      txser.inp
C **      sfser.inp
C **      txser.inp
C **      qctl.inp
C **      mask.inp
C **      show.inp
C **      vege.inp
C **      moddxdy.inp
C **      modchan.inp
C **      gwater.inp
C
C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  Last modified by J. Shen on Oct. 2018
C
C     This code is revsied for James River Water quality project. The code was only used
C     for the James River and some parts of the code is hardwired, which are note
C     suitable for use it in different area. The model is customized for run
C     water quality model in decoupled mode.
C     To run water quality model, it needs to run hydrodynamc model first to save
C     model resutls. Inputs and outputs files use different format in order 
C     to run mutiple years.
C  
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      CHARACTER*80 TITLE
	character fname*10
C
C**********************************************************************C
C
C **  set ieee exception traps (SUN SYSTEMS ONLY)
C
C     CALL SETUP_IEEE
C
C**********************************************************************C
C
c     CALL WELCOME
C
C**********************************************************************C
C
C **  OPEN OUTPUT FILES
C
C----------------------------------------------------------------------C
C
      OPEN(7,FILE='efdc.out',STATUS='UNKNOWN')
      OPEN(8,FILE='efdc.log',STATUS='UNKNOWN')
      OPEN(9,FILE='time.log',STATUS='UNKNOWN')
      OPEN(10,FILE='drywet.log',STATUS='UNKNOWN')
      OPEN(1,FILE='vsfp.out',STATUS='UNKNOWN')
C
      CLOSE(7,STATUS='DELETE')
      CLOSE(8,STATUS='DELETE')
      CLOSE(9,STATUS='DELETE')
      CLOSE(10,STATUS='DELETE')
      CLOSE(1,STATUS='DELETE')
C
      OPEN(7,FILE='efdc.out',STATUS='UNKNOWN')
      OPEN(8,FILE='efdc.log',STATUS='UNKNOWN')
      OPEN(9,FILE='time.log',STATUS='UNKNOWN')
C
      OPEN(1,FILE='sediag.out',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
C
      OPEN(1,FILE='cfl.out',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      
      OPEN(1,FILE='assim.out',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='assim.out',STATUS='UNKNOWN')
      write(1,*)'Data Assimilation log'
      CLOSE(1)
      
            
C
c------------------------
c
      IC9=0
      JC9=0
      LDM9=0
      call grid(IC9,JC9,LDM9)             ! read grid.inp, ji, 9/10/99
c
c check & read hsdir.inp from SWAN model for the first time, Ji@wave, 3/5/00
      ISWAN=0
      IWCM=1        ! for calling wcm2000.for
      dayswan1=-1.0e10
      dayswan2=1.0e10
      lswan=-999
c
c      open(480,file='hsdir.inp',status='old',err=97) ! from readswan3.for
      open(480,file='..\hsdir.inp',status='old',err=97) ! from readswan3.for
      open(481,file='..\wcm.inp',status='old',form='unformatted',
     &                            access='DIRECT',recL=1,err=98)  ! from wcm.for
      iwcm=0         ! used wcm.inp from wcm.for
      write(6,*) "wcm.inp from wcm.for is available"
      ISWAN=1
      read(480,*) dayswan1,dayswan2,Lswan
      write(6,*) "hsdir.inp from SWAM model is available"
      read(480,692) dayswan,NX
692   format(f14.4,4x,i9)
      write(6,*) "hsdir.inp N= ",NX,dayswan
      read(480,691,IOSTAT=ISO) (ihs(L),L=1,LSWAN),(IDIR(L),L=1,LSWAN),
     +   (ITM(L),L=1,LSWAN)
691   format(20i5)
      IF(ISO.GT.0) then
      write(6,*) "hsdir.inp READ error: ", dayswan,NX
      stop
      endif
c
      days2=dayswan
      do L=1,LSWAN
      L2=L+1                 ! shift back to original EFDC L sequence
      HS2(L2) =IHS(L)*0.001  ! -> M
      DIR2(L2)=Idir(L)*0.1   ! -> deg
      if(idir(L).eq.-999) dir2(L2)=-999.9
      TM2(L2) =ITM(L)*0.001  ! -> second
      if(ITM(L).eq.-999) TM2(L2)=-999.9
      enddo
      go to 96
c
97    write(6,*) "hsdir.inp does not exist and SWAN is not used"
      go to 96
98    write(6,*) "wcm.inp is not used"
96    continue
c---------------------------
C
C**********************************************************************C
C
C **  CALL INPUT SUBROUTINE
C
      CALL INPUT(TITLE)
C
C**********************************************************************C
C
C **  CALL SUBROUTINE TO ADJUST, CONVERT AND SMOOTH DEPTH FIELD
C
      IF(NSHMAX.GE.1) CALL DEPSMTH
C
C**********************************************************************C
C
C **  SET TIME RELATED PARAMETERS
C
C **  THE PARAMETER NTC=NUMBER OF TIME CYCLES, CONTROLS
C **  THE LENGTH OF RUN (NUMBER OF TIME STEPS)
C
C----------------------------------------------------------------------C
C
      ISCRAY=0
C
      TCYCLE=0.0
      TLRPD=0.0
      THDMT=0.0
      TVDIF=0.0
      TCGRS=0.0
      TCONG=0.0
      TSADV=0.0
      TRELAXV=0.0
      TPUV=0.0
      TCEXP=0.0
      TAVB=0.0
      TUVW=0.0
      TQQQ=0.0
      TQCLT=0.0
      TWQADV=0.0
      TWQDIF=0.0
      TWQKIN=0.0
      TWQSED=0.0
      WTLRPD=0.0
      WTHDMT=0.0
      WTVDIF=0.0
      WTCGRS=0.0
      WTCONG=0.0
      WTSADV=0.0
      WTRELAXV=0.0
      WTPUV=0.0
      WTCEXP=0.0
      WTAVB=0.0
      WTUVW=0.0
      WTQQQ=0.0
      WTQCLT=0.0
      WTWQADV=0.0
      WTWQDIF=0.0
      WTWQKIN=0.0
      WTWQSED=0.0
      CFMAX=CF
      PI=3.1415926535898
      NBAN=49
      TPN=FLOAT(NTSPTC)
      NTS=NTC*NTSPTC/NFLTMT
      NLTS=NTSPTC*NLTC
      NTTS=NTSPTC*NTTC
      SNLT=0.
      NCTBC=1
      NPRINT=1
      NTSPP=NTCPP*NTSPTC/NFLTMT
      NTSNB=NTCNB*NTSPTC
      NTSVB=NTCVB*NTSPTC
      ITRMAX=0
      ITRMIN=1000
      ERRMAX=1E-9
      ERRMIN=1000.
      NMMT=1
      NBAL=1
      NBALE=1
      NBALO=1
      NBUD=1
      NHAR=1
      NTSPTC2=2*NTSPTC/NFLTMT
C     NDISP=NTS-NTSPTC+1
      NDISP=NTS-NTSPTC+2
      NSHOWR=0
      NSHOWC=0
C
      DO NS=1,NASER
        MATLAST(NS)=1
      END DO
C
      DO NS=1,NWSER
        MWTLAST(NS)=1
      END DO
C
      DO NS=1,NPSER
       MPTLAST(NS)=1
      END DO
C
      DO NS=1,NQSER
       MQTLAST(NS)=1
      END DO
C
      DO NS=1,NQWRSR
       MQWRTLST(NS)=1
      END DO
C
      NTOX1=NTOX
      IF(NPCB.GT.0)NTOX1=NPCB
      
      NTMP=4+NSED+NSND+NTOX1
      DO NC=1,NTMP
      DO NN=1,NCSER(NC)
      MCTLAST(NN,NC)=1
      END DO
      END DO
C
       MSFTLST=1
       nefdcwin=1
       SECDLAST=0
       NSTDAY=0
       NCSTEP=0
       DAYLAST=0
       NTS_S=0
       NLAST1=0
       NLAST=1
       SECDLAST4=0
       DO I=1,200
       CDAY(I)=5000
       NTSPTCC(I)=NTSPTC
       ENDDO
       open(nefdcwin,file='Timestep.inp',status='old',err=940)
       DO I=1,9
       READ(nefdcwin,*)
       ENDDO
       READ(nefdcwin,*)NCSTEP
       IF(NCSTEP.GT.0)THEN
        READ(nefdcwin,*)
        READ(nefdcwin,*)
        DO I=1,NCSTEP
         READ(nefdcwin,*)CDAY(I),NTSPTCC(I),ISDRYC(I)
        ENDDO
        NTSPTC=NTSPTCC(1)  ! Reset timestep for the model using 1st day input
!
! make sure change time step can only occur at correspoinding days for average resutls
! 
        open(1,file='Timestep_n.inp')
        Write(1,*)'Revised Timestep file'
        DO i=1,8
        Write(1,*)'C'
        ENDDO

        
        DO I=2,NCSTEP
         IF(CDAY(I).gt.TBEGIN)THEN
         NLAST=I
         GOTO 102
         ENDIF
        ENDDO
  102   continue ! NTSPTCC(NLAST-1)=NTSPTC     ! Set 1st day and timestep eq. EFDC.inp
 
        if(NTSMMT.LT.NTSPTCC(NLAST-1))THEN
           NTSMMT=NTSPTCC(NLAST-1)/24           ! force avg. in hours
        else
           NTSMMT=NTSPTCC(NLAST-1)              ! force avg. in days
        endif  
        
        CDAY(NLAST-1)=TBEGIN
        Write(1,'(F9.1,2I10)')CDAY(NLAST-1),NTSPTCC(NLAST-1),
     $    ISDRYC(NLAST-1)
c 
        IS_DAY=1 !int(aver3d/24)
        aver3d=real(NTSMMT)/NTSPTC*24.0  ! using one day
        DO I=NLAST,NCSTEP
        IF(CDAY(I).EQ.CDAY(I-1))CDAY(I)=CDAY(I)+1

        IF(mod(int(CDAY(I)-CDAY(I-1)),IS_DAY).NE.0) THEN
        LL=mod(int(CDAY(I)-CDAY(I-1)),IS_DAY)
        CDAY(I)=CDAY(I)-(LL)
        IF( (CDAY(I)-CDAY(I-1)).LT. IS_DAY) THEN
        CDAY(I)=CDAY(I-1)+IS_DAY
        ENDIF
        ENDIF
        write(1,'(F9.1,4I10)')CDAY(I),NTSPTCC(I),ISDRYC(I),
     *  int(CDAY(I)-CDAY(I-1)),int(CDAY(I)-CDAY(I-1))*NTSPTCC(I-1)      
        ENDDO      
        
        DO I=NLAST,NCSTEP
         IF(CDAY(I).gt.TBEGIN+NTC) THEN
          NLAST1=I
          CDAY(I)=TBEGIN+NTC
          GOTO 103
         ENDIF
        ENDDO
        IF(NLAST1.EQ.0) THEN
         write(*,*)'Check Timestep.inp file to match the starting time'
         stop
        ENDIF
        
  103   CDAY(NLAST1)=TBEGIN+NTC
        DO I=NLAST,NLAST1
         NTS_S=NTS_S+(CDAY(I)-CDAY(I-1))*NTSPTCC(I-1)
         NCDAY(I)=(CDAY(I)-CDAY(I-1))*NTSPTCC(I-1)+NCDAY(I-1) 
         write(1,*)'Change timestep = ',NCDAY(I)
         IF(CDAY(I)-CDAY(I-1).LT.0) THEN
         write(*,*)'Check Timestep.inp file at ', I
         stop
         ENDIF
        ENDDO
        NTS=NTS_S                ! Reset
      Close(1)
      ENDIF

  940 CONTINUE
C
C**********************************************************************C
C
C **  SET CONTROLS FOR WRITING TO INSTANTANEOUS 2D SCALAR CONTOURING
C **  AND 2D VELOCITY VECTOR PLOTTING FILES
C
C----------------------------------------------------------------------C
C
C **  SCALAR FIELD CONTOURING IN HORIZONTAL PLANES: SUBROUTINE SALPLTH
C
      DO N=1,7
      IF(ISSPH(N).EQ.1) THEN
        NCSPH(N)=NTS-(NTSPTC-(NTSPTC/NPSPH(N)))/NFLTMT
        JSSPH(N)=1
       ELSE
        NCSPH(N)=0
        JSSPH(N)=0
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
C **  FREE SURFACE ELEVATION OR PRESSURE CONTOURING IN HORIZONTAL
C **  PLANES: SUBROUTINE SURFPLT
C
C     ISPPH=0
C     NPPPH=12
      IF(ISPPH.EQ.1)  THEN
        NCPPH=NTS-(NTSPTC-(NTSPTC/NPPPH))/NFLTMT
        JSPPH=1
       ELSE
        NCPPH=0
        JSPPH=0
      END IF
C
C----------------------------------------------------------------------C
C
C **  VELOCITY VECTOR PLOTTING IN HORIZONTAL PLANES: SUBROUTINE VELPLTH
C
      IF(ISVPH.EQ.1) THEN
        NCVPH=NTS-(NTSPTC-(NTSPTC/NPVPH))/NFLTMT
        JSVPH=1
       ELSE
        NCVPH=0
        JSVPH=0
      END IF
C
C----------------------------------------------------------------------C
C
C **  SCALAR FIELD CONTOURING IN VERTICAL PLANES: SUBROUTINE SALPLTV
C
      DO N=1,7
      IF(ISSPV(N).GE.1) THEN
        NCSPV(N)=NTS-(NTSPTC-(NTSPTC/NPSPV(N)))/NFLTMT
        JSSPV(N)=1
         DO IS=1,ISECSPV
         CCTITLE(20+IS)=CCTITLE(10+IS)
         CCTITLE(30+IS)=CCTITLE(10+IS)
         CCTITLE(40+IS)=CCTITLE(10+IS)
         CCTITLE(50+IS)=CCTITLE(10+IS)
         END DO
       ELSE
        NCSPV(N)=0
        JSSPV(N)=0
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
C **  NORMAL VELOCITY CONTOURING AND TANGENTIAL VELOCITY VECTOR
C **  PLOTTING IN VERTICAL PALNES: SUBROUTINE VELPLTV
C
      IF(ISVPV.GE.1) THEN
        NCVPV=NTS-(NTSPTC-(NTSPTC/NPVPV))/NFLTMT
        JSVPV=1
         DO IS=1,ISECVPV
         CVTITLE(20+IS)=CVTITLE(10+IS)
         CVTITLE(30+IS)=CVTITLE(10+IS)
         CVTITLE(40+IS)=CVTITLE(10+IS)
         CVTITLE(50+IS)=CVTITLE(10+IS)
         CVTITLE(60+IS)=CVTITLE(10+IS)
         CVTITLE(70+IS)=CVTITLE(10+IS)
         CVTITLE(80+IS)=CVTITLE(10+IS)
         CVTITLE(90+IS)=CVTITLE(10+IS)
         END DO
       ELSE
        NCVPV=0
        JSVPV=0
      END IF
C
C----------------------------------------------------------------------C
C
C **  THREE-DIMENSIONAL HDF FORMAT GRAPHICS FILES: SUBROUTINE OUT3D
C
      IF(IS3DO.EQ.1) THEN
       NC3DO=NTS-(NTSPTC-(NTSPTC/NP3DO))/NFLTMT
      END IF
C
C**********************************************************************C
C
C **  SET CONTROLS FOR WRITING TO FILTERED, AVERAGED OR RESIDUAL
C **  2D SCALAR CONTOURING AND 2D VELOCITY VECTOR PLOTTING FILES
C
C----------------------------------------------------------------------C
C
C **  RESIDUAL SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATION
C **  CONTOURING IN HORIZONTAL: SUBROUTINE RSALPLTH
C
      DO N=1,7
      IF(ISRSPH(N).GE.1) JSRSPH(N)=1
      END DO
C
C----------------------------------------------------------------------C
C
C **  RESIDUAL VELOCITY VECTOR PLOTTING IN HORIZONTAL PLANES:
C **  SUBROUTINE RVELPLTH
C
      IF(ISRVPH.GE.1) JSRVPH=1
C
C----------------------------------------------------------------------C
C
C **  RESIDUAL SURFACE ELEVATION PLOTTING IN HORIZONTAL PLANES:
C **  SUBROUTINE RVELPLTH
C
      IF(ISRPPH.GE.1) JSRPPH=1
C
C----------------------------------------------------------------------C
C
C **  RESIDUAL SCALAR FIELD CONTOURING IN VERTICAL
C **  PLANES: SUBROUTINE RSALPLTV
C
      DO N=1,7
      IF(ISRSPV(N).GE.1) JSRSPV(N)=1
      END DO
C
C----------------------------------------------------------------------C
C
C **  RESIDUAL NORMAL AND TANGENTIAL VELOCITY CONTOURING AND AND
C **  TANGENTIAL VELOCITY VECTOR PLOTTING IN VERTICAL PLANES:
C **  SUBROUTINE RVELPLTV
C
      IF(ISRVPV.GE.1) JSRVPV=1
C
C**********************************************************************C
C
C **  SET CONTROLS FOR WRITING TO DRIFTER, HARMONIC ANALYSIS,
C **  RESIDUAL TRANSPORT, AND BLANCE OUTPUT FILES
C
C----------------------------------------------------------------------C
C
      JSPD=1
      NCPD=1
C
      JSLSHA=1
      IF (ISLSHA.EQ.1) THEN
       LSLSHA=0
       NCLSHA=NTS-NTCLSHA*NTSPTC
      END IF
C
      IF(ISRESTR.EQ.1) JSRESTR=1
      JSWASP=0
      IF(ISWASP.GE.1) JSWASP=1
C
      IF (ISBAL.GE.1) THEN
        JSBAL=1
        JSBALO=1
        JSBALE=1
      END IF
      JSSBAL=1
C
C**********************************************************************C
C
C **  SET CONTROL FOR CALCULATION OF LAGRANGIAN MEAN VELOCITIY FIELDS
C **  BY PARTICLE TRACKING
C
C----------------------------------------------------------------------C
C
      IF(ISLRPD.GE.1) THEN
        NLRPDRT(1)=NTS-NTSPTC-(NTSPTC-(NTSPTC/MLRPDRT))/NFLTMT
        DO M=2,MLRPDRT
        NLRPDRT(M)=NLRPDRT(M-1)+(NTSPTC/MLRPDRT)
        END DO
        JSLRPD=1
       ELSE
        NLRPDRT(1)=NTS+2
        JSLRPD=0
      END IF
C
C**********************************************************************C
C
C **  SET SOME CONSTANTS
C
C----------------------------------------------------------------------C
C
      JSTBXY=0
C
      CTURB2=CTURB**0.667
      CTURB3=CTURB**0.333
C
      KS=KC-1
      IF (KS.EQ.0) KS=1
      DZI=FLOAT(KC)
      DZ=1./DZI
      DZS=DZ*DZ
C
      DT=TIDALP*FLOAT(NFLTMT)/FLOAT(NTSPTC)
      DTI=1./DT
      DT2=2.*DT
C
      AVCON1=2.*(1.-AVCON)*DZI*AVO
      G=9.81
      GPO=G*BSC
      GI=1./G
      GID2=.5*GI
      PI=3.1415926535898
      PI2=2.*PI
C
      TCVP=0.0625*TIDALP/PI
C
C**********************************************************************C
C
C **  SET CONSTANTS FOR M2 TIDAL CYCLE HARMONIC ANALYSIS
C
C----------------------------------------------------------------------C
C
      AC=0.
      AS=0.
      ACS=0.
C
      DO N=1,NTSPTC
      TNT=FLOAT(N)
      NP=NTSPTC+N
      WC(N)=COS(2.*PI*TNT/TPN)
      WS(N)=SIN(2.*PI*TNT/TPN)
      WC(NP)=WC(N)
      WS(NP)=WS(N)
      AC=AC + 2.*WC(N)*WC(N)
      AS=AS + 2.*WS(N)*WS(N)
C     ACS=ACS +2.*WC(N)*WS(N)
      ACS=0.
      WC2(N)=COS(4.*PI*TNT/TPN)
      WS2(N)=SIN(4.*PI*TNT/TPN)
      WC2(NP)=WC2(N)
      WS2(NP)=WC2(N)
      AC2=AC2 + 2.*WC2(N)*WC2(N)
      AS2=AS2 + 2.*WS2(N)*WS2(N)
C     ACS2=ACS2 +2.*WC2(N)*WS2(N)
      ACS2=0.
      END DO
C
      DET=AC*AS-ACS*ACS
      AS=AS/DET
      AC=AC/DET
      ACS=ACS/DET
      DET=AC2*AS2-ACS2*ACS2
      AS2=AS2/DET
      AC2=AC2/DET
      ACS2=ACS2/DET
C
C**********************************************************************C
C
C **  SET WEIGHTS FOR SALINITY AND TEMPERATURE BOUNDARY INTERPOLATION
C
C----------------------------------------------------------------------C
C
      IF (KC.GT.1) THEN
      DO K=1,KC
      WTCI(K,1)=FLOAT(K-KC)/FLOAT(1-KC)
      WTCI(K,2)=FLOAT(K-1)/FLOAT(KC-1)
      END DO
      ELSE
      WTCI(1,1)=0.5
      WTCI(1,2)=0.5
      END IF
C
C**********************************************************************C
C
C **  INITIALIZE ARRAYS
C
      CALL AINIT
C
C**********************************************************************C
C
C **  READ IN XLON AND YLAT OR UTME AND UTMN OF CELL CENTERS OF
C **  CURVILINEAR PORTION OF THE  GRID
C
      IF (ISCLO.EQ.1) THEN
C
      OPEN(1,FILE='lxly.inp',STATUS='UNKNOWN')
C
      DO NS=1,4
      READ(1,1111)
      END DO
 1111 FORMAT(80X)
C
      IF(ISCORV.EQ.1) THEN
       DO LL=1,LVC
       READ(1,*,ERR=3000)I,J,XLNUTME,YLTUTMN,CCUE,CCVE,CCUN,CCVN,TMPCOR
       L=LIJ(I,J)
       DLON(L)=XLNUTME
       DLAT(L)=YLTUTMN
       CUE(L)=CCUE
       CVE(L)=CCVE
       CUN(L)=CCUN
       CVN(L)=CCVN
       CVN(L)=CCVN
       FCORC(L)=TMPCOR
       DETTMP=1./( CUE(L)*CVN(L)-CUN(L)*CVE(L) )
       IF(DETTMP.EQ.0.0) THEN
         WRITE(6,6262)
         WRITE(6,6263)IL(L),JL(L)
         STOP
       END IF
       END DO
      ELSE
       DO LL=1,LVC
       READ(1,*,ERR=3000)I,J,XLNUTME,YLTUTMN,CCUE,CCVE,CCUN,CCVN
       L=LIJ(I,J)
       DLON(L)=XLNUTME
       DLAT(L)=YLTUTMN
       CUE(L)=CCUE
       CVE(L)=CCVE
       CUN(L)=CCUN
       CVN(L)=CCVN
       CVN(L)=CCVN
       FCORC(L)=CF
       DETTMP=1./( CUE(L)*CVN(L)-CUN(L)*CVE(L) )
       IF(DETTMP.EQ.0.0) THEN
         WRITE(6,6262)
         WRITE(6,6263)IL(L),JL(L)
         STOP
       END IF
       END DO
      END IF
C
 6262 FORMAT('  SINGULAR INVERSE TRANSFORM FROM E,N TO CURV X,Y')
 6263 FORMAT('  I,J =',2I10/)
C
      CLOSE (1)
C
      END IF
c calculate curvelinear coordinate angle, JI@wave, 3/5/00
      anglem(1)=0.0
      anglem(LC)=0.0
      do L=2,LA
      angle1=acos(cue(L))               ! acos -> (0,PI)
c     angle2=-asin(cve(L))
c     angle3=asin(cun(L))
      angle4=acos(cvn(L))
c     anglem(l)=(angle1+angle2+angle3+angle4)*0.25  ! in radians, not in deg!
      anglem(l)=(angle1+              angle4)*0.50  ! in radians, not in deg!
      enddo
c-----------------
C
      FCORC(1)=FCORC(2)
      FCORC(LC)=FCORC(LA)
C
      GO TO 3002
 3000 WRITE(6,3001)
 3001 FORMAT(1X,'READ ERROR FOR FILE lxly.inp ')
      STOP
 3002 CONTINUE
C
c     COSTMP=COSD(15.)
c     SINTMP=SIND(15.)
c     SINNEG=-SINTMP
c      OPEN(1,FILE='newlxly.inp',STATUS='UNKNOWN')  !Ji, 11/3/00
      OPEN(1,FILE='newlxly.out',STATUS='UNKNOWN')
      IF(ISCORV.EQ.1) THEN
       DO L=2,LA
        WRITE(1,1112)IL(L),JL(L),DLON(L),DLAT(L),CUE(L),CVE(L),CUN(L),
     $             CVN(L),FCORC(L)
       END DO
      ELSE
       DO L=2,LA
        WRITE(1,1112)IL(L),JL(L),DLON(L),DLAT(L),CUE(L),CVE(L),CUN(L),
     $             CVN(L)
       END DO
      END IF
      CLOSE(1)
C
      ZERO=0.
      OPEN(1,FILE='lijmap.out',STATUS='UNKNOWN')
      DO L=2,LA
      WRITE(1,1113)L,IL(L),JL(L),ZERO
      END DO
      CLOSE(1)
C
 1112 FORMAT (2I5,2F12.4,4F12.7)
 1113 FORMAT (3I5,F10.2)
C
C**********************************************************************C
C
C **  READ VEGETATION DATA
C
      IF (ISVEG.GE.1) THEN
C
      OPEN(1,FILE='vege.inp',STATUS='UNKNOWN')
C
      DO NS=1,9
      READ(1,1111)
      END DO
C
      READ(1,*)MVEGTYP,MVEGOW,UVEGSCL
      DO M=1,MVEGTYP
      READ(1,*,ERR=3120)IDUM,RDLPSQ(M),BPVEG(M),HPVEG(M),ALPVEG(M),
     $         BETVEG(M),GAMVEG(M),SCVEG(M)
      BDLTMP=BPVEG(M)*BPVEG(M)*RDLPSQ(M)
      PVEGX(M)=1.-BETVEG(M)*BDLTMP
      PVEGY(M)=1.-BETVEG(M)*BDLTMP
      PVEGZ(M)=1.-ALPVEG(M)*BDLTMP
      BDLPSQ(M)=BPVEG(M)*RDLPSQ(M)
      END DO
C
      CLOSE(1)
C
      END IF
C
      GO TO 3122
 3120 WRITE(6,3121)
 3121 FORMAT(1X,'READ ERROR FOR FILE vege.inp ')
      STOP
 3122 CONTINUE
C
C**********************************************************************C
C
C **  READ IN COUNTER CLOCKWISE ANGLE FROM EAST SPECIFYING
C **  PRINCIPAL FLOOD FLOW DIRECTION
C
      IF(ISTRAN(4).GE.1.AND.ISSFLFE.GE.1) THEN
C
      OPEN(1,FILE='fldang.inp',STATUS='UNKNOWN')
C
      DO LL=2,LA
      READ(1,*,ERR=3130)I,J,ANGTMP1,ANGTMP2
      L=LIJ(I,J)
      ACCWFLD(L,1)=0.0174533*ANGTMP1
      ACCWFLD(L,2)=0.0174533*ANGTMP2
      END DO
C
      CLOSE(1)
C
      END IF
C
      GO TO 3132
 3130 WRITE(6,3131)
 3131 FORMAT(1X,'READ ERROR FOR FILE fldang.inp ')
      STOP
 3132 CONTINUE
C
C**********************************************************************C
C
C **  SET BOUNDARY CONDITION SWITCHES
C
      CALL SETBCS
C
      IF (ISRESTI.EQ.0.AND.ISTRAN(1).GE.1.AND.ifed_inc.GE.0) THEN
        OPEN(1,FILE='salt.inp',STATUS='UNKNOWN')
C **    SKIP OVER TITLE AND AND HEADER LINES
        DO IS=1,4
         READ(1,*)
        END DO
        READ(1,*)ISALTYP
        IF(ISALTYP.EQ.0) THEN
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO) (SALINIT(L,K),K=1,KC)
           IF(ISO.GT.0) GO TO 840
          END DO
         ELSE
          DO L=2,LC-1
           READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,
     &       (SALINIT(LIJ(IDUM,JDUM),K),K=1,KC)
           IF(ISO.GT.0) GO TO 840
          END DO
        END IF
       DO K=1,KC
       DO L=2,LA
       SAL(L,K)=SALINIT(L,K)
       SAL1(L,K)=SALINIT(L,K)
       END DO
       END DO     
        CLOSE(1)
      END IF
 840  continue   
C
C
C  J.S  4/8/2011
C   Read tidelkb.inp to modify resuspension
C
       DO L=2,LA
       TAU_MD(L)=1.0
       ENDDO
       TT0=-999
       IF(iMDTau.GE.1.and.irunpcb.EQ.1) then
        open(1,file='tidelkb.inp') 
        DO I=1,4
        read(1,*)
        ENDDO
        IF(iMDTau.EQ.1.or.iMDTau.EQ.2)THEN       
         DO L=2,LA
          read(1,*)IT_I,IT_J,TP_X,TP_Y,TP_1,TP_2,TP_3
          LL=LIJ(IT_I,IT_J)
          TAU_MD(LL)=abs(TP_1)*TP_2  ! normlized shear stress with modification
          TAU_MD1(LL)=TP_3           ! modification of suspension rate
          IF(TP_1.GT.TT0)TT0=TP_1        
         ENDDO
         DO L=2,LA
       !   TAU_MD(L)=TAU_MD(L)/TT0  !divid by max. stress
      !    TAU_MD(L)=TAU_MD(L)**Tauexp       JS. 2017 don't do this
   !       TAU_MD(L)=min(1.0/TAU_MD(L),5.0)
         TAU_MD(L)=min(TAU_MD(L),5.0)
         ENDDO 
         IF(iMDTau.EQ.2)THEN       
         DO L=2,LA
          TAU_MD(L)=TAU_MD(L)*TAU_MD1(L)
         ENDDO 
         ENDIF 
         open(2,file='tidelkb.dia')  
         DO L=2,LA
          write(2,*)IL(L),JL(L),TAU_MD(L),TAU_MD1(L)       
         ENDDO                 
         CLOSE(2) 
        ENDIF

        IF(iMDTau.GE.3) then
         DO L=2,LA
          read(1,*)IT_I,IT_J,TP_1,TP_2
          TAU_MD(LL)=TP_2      
         ENDDO         
        ENDIF
        CLOSE(1)
       ENDIF        
 
C**********************************************************************C
C
C **  READ RESTART CONDITIONS OR INITIALIZE SCALAR FIELDS
C
C     ISRESTI.EQ.10 READS AND OLD RESTART FILE GENERATED BY
C     PRE SEPTEMBER 8, 1992 VERSIONS OF EFDC.FOR
C
C----------------------------------------------------------------------C
C
      IF (ISLTMT.EQ.0) THEN
      IF (ISRESTI.GE.1) THEN
c       IF (ISRESTI.EQ.1) CALL RESTIN1  ! Ji, 10/31/00
       IF (ISRESTI.EQ.1) CALL RESTIN9
       IF (ISRESTI.EQ.2) CALL RESTIN2
       IF (ISRESTI.EQ.10) CALL RESTIN10
      END IF
      IF (ISRESTI.EQ.-1) CALL RESTIN1
      END IF
C
C----------------------------------------------------------------------C
C
C **  INTIALIZE SALINITY FIELD IF NOT READ IN FROM RESTART FILE
C
      IF (ISRESTI.EQ.0.AND.ISTRAN(1).GE.1) THEN
      IF (ISLTMT.EQ.0.AND.ISTOPT(1).GE.1) THEN
C
       NREST=0
C
       DO K=1,KC
       DO L=2,LA
       SAL(L,K)=SALINIT(L,K)
       SAL1(L,K)=SALINIT(L,K)
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,1)=SALINIT(L,K)
       NLOS(LL,K,1)=0
       IF(NCSERS(LL,1).EQ.0) THEN
        SAL(L,K)=WTCI(K,1)*CBS(LL,1,1)+WTCI(K,2)*CBS(LL,2,1)
        SAL1(L,K)=SAL(L,K)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,1)=SALINIT(L,K)
       NLOW(LL,K,1)=0
       IF(NCSERW(LL,1).EQ.0) THEN
        SAL(L,K)=WTCI(K,1)*CBW(LL,1,1)+WTCI(K,2)*CBW(LL,2,1)
        SAL1(L,K)=SAL(L,K)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,1)=SALINIT(L,K)
       NLOE(LL,K,1)=0
       IF(NCSERE(LL,1).EQ.0) THEN
        SAL(L,K)=WTCI(K,1)*CBE(LL,1,1)+WTCI(K,2)*CBE(LL,2,1)
        SAL1(L,K)=SAL(L,K)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,1)=SALINIT(L,K)
       NLON(LL,K,1)=0
       IF(NCSERN(LL,1).EQ.0) THEN
        SAL(L,K)=WTCI(K,1)*CBN(LL,1,1)+WTCI(K,2)*CBN(LL,2,1)
        SAL1(L,K)=SAL(L,K)
       END IF
       END DO
       END DO
C
       OPEN(1,FILE='newsalt.out',STATUS='UNKNOWN')
       IONE=1
       WRITE(1,9101)IONE
       DO L=2,LC-1
       WRITE(1,9102)L,IL(L),JL(L),(SAL(L,K),K=1,KC)
       END DO
       CLOSE(1)
C
      END IF
      END IF
C
 9101 FORMAT(I5)
 9102 FORMAT(3I5,12F8.2)
C

C
C **  INTIALIZE TEMP FIELD IF NOT READ IN FROM RESTART FILE
C
      IF (ISRESTI.EQ.0.AND.ISTRAN(2).GE.1) THEN
      IF (ISLTMT.EQ.0.AND.ISTOPT(2).EQ.1) THEN
C
       NREST=0
C
       DO K=1,KC
       DO L=2,LA
       TEM(L,K)=TEMINIT(L,K)
       TEM1(L,K)=TEM(L,K)
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,2)=TEMINIT(L,K)
       NLOS(LL,K,2)=0
       IF(NCSERS(LL,2).EQ.0) THEN
        TEM(L,K)=WTCI(K,1)*CBS(LL,1,2)+WTCI(K,2)*CBS(LL,2,2)
        TEM1(L,K)=TEM(L,K)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,2)=TEMINIT(L,K)
       NLOW(LL,K,2)=0
       IF(NCSERW(LL,2).EQ.0) THEN
        TEM(L,K)=WTCI(K,1)*CBW(LL,1,2)+WTCI(K,2)*CBW(LL,2,2)
        TEM1(L,K)=TEM(L,K)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,2)=TEMINIT(L,K)
       NLOE(LL,K,2)=0
       IF(NCSERE(LL,2).EQ.0) THEN
        TEM(L,K)=WTCI(K,1)*CBE(LL,1,2)+WTCI(K,2)*CBE(LL,2,2)
        TEM1(L,K)=TEM(L,K)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,2)=TEMINIT(L,K)
       NLON(LL,K,2)=0
       IF(NCSERN(LL,2).EQ.0) THEN
        TEM(L,K)=WTCI(K,1)*CBN(LL,1,2)+WTCI(K,2)*CBN(LL,2,2)
        TEM1(L,K)=TEM(L,K)
       END IF
       END DO
       END DO
C
      END IF
      END IF
C
C **  INTIALIZE TEMPERATURE BC IF NOT READ IN FROM RESTART FILE
C     AND CONSTANT INTIAL CONDITION IS USED
C
      IF (ISRESTI.EQ.0.AND.ISTRAN(2).GE.1) THEN
      IF (ISLTMT.EQ.0.AND.ISTOPT(2).GE.2) THEN
       M=2
C
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       NSID=NCSERS(LL,M)
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
       CLOS(LL,K,M)=TEMO
       NLOS(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       NSID=NCSERW(LL,M)
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       CLOW(LL,K,M)=TEMO
       NLOW(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       NSID=NCSERE(LL,M)
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       CLOE(LL,K,M)=TEMO
       NLOE(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       NSID=NCSERN(LL,M)
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
       CLON(LL,K,M)=TEMO
       NLON(LL,K,M)=0
       END DO
       END DO
C
      END IF
      END IF
C
C **  INTIALIZE DYE FIELD IF NOT READ IN FROM RESTART FILE
C
      IF (ISRESTI.EQ.0.AND.ISTRAN(3).GE.1) THEN
      IF (ISLTMT.EQ.0.AND.ISTOPT(3).GE.1) THEN
C
       NREST=0
C
       DO K=1,KC
       DO L=2,LA
       DYE(L,K)=DYEINIT(L,K)
       DYE1(L,K)=DYE(L,K)
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,3)=DYEINIT(L,K)
       NLOS(LL,K,3)=0
       IF(NCSERS(LL,3).EQ.0) THEN
        DYE(L,K)=WTCI(K,1)*CBS(LL,1,3)+WTCI(K,2)*CBS(LL,2,3)
        DYE1(L,K)=DYE(L,K)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,3)=DYEINIT(L,K)
       NLOW(LL,K,3)=0
       IF(NCSERW(LL,3).EQ.0) THEN
        DYE(L,K)=WTCI(K,1)*CBW(LL,1,3)+WTCI(K,2)*CBW(LL,2,3)
        DYE1(L,K)=DYE(L,K)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,3)=DYEINIT(L,K)
       NLOE(LL,K,3)=0
       IF(NCSERE(LL,3).EQ.0) THEN
        DYE(L,K)=WTCI(K,1)*CBE(LL,1,3)+WTCI(K,2)*CBE(LL,2,3)
        DYE1(L,K)=DYE(L,K)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,3)=DYEINIT(L,K)
       NLON(LL,K,3)=0
       IF(NCSERN(LL,3).EQ.0) THEN
        DYE(L,K)=WTCI(K,1)*CBN(LL,1,3)+WTCI(K,2)*CBN(LL,2,3)
        DYE1(L,K)=DYE(L,K)
       END IF
       END DO
       END DO
C
      END IF
      END IF
C
C **  INTIALIZE DYE BC IF NOT READ IN FROM RESTART FILE
C **  AND CONSTANT INITIAL CONDITIONS ARE USED
C
      IF (ISRESTI.EQ.0.AND.ISTRAN(3).GE.1) THEN
      IF (ISLTMT.EQ.0.AND.ISTOPT(3).EQ.0) THEN
C
       M=3
C
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       NSID=NCSERS(LL,M)
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
       CLOS(LL,K,M)=0.
       NLOS(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       NSID=NCSERW(LL,M)
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       CLOW(LL,K,M)=0.
       NLOW(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       NSID=NCSERE(LL,M)
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       CLOE(LL,K,M)=0.
       NLOE(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       NSID=NCSERN(LL,M)
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
       CLON(LL,K,M)=0.
       NLON(LL,K,M)=0
       END DO
       END DO
C
      END IF
      END IF
C
C **  INTIALIZE TOX AND BC IF NOT READ IN FROM RESTART FILE
C **  AND VARIABLE INITIAL CONDITIONS ARE USED
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(5).EQ.0) IISTMP=0
C
      IF (IISTMP.EQ.0.AND.ISTRAN(5).EQ.1) THEN
       DO NT=1,NTOX1
       IF(ITXINT(NT).EQ.1.OR.ITXINT(NT).EQ.3) THEN
       M=4+NT

       DO K=1,KC
       DO L=2,LA
       TOX(L,K,NT)=TOXINIT(L,K,NT)
       TOX1(L,K,NT)=TOX(L,K,NT)
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,M)=TOXINIT(L,K,NT)
       NLOS(LL,K,M)=0
       IF(NCSERS(LL,M).EQ.0) THEN
        TOX(L,K,NT)=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
        TOX1(L,K,NT)=TOX(L,K,NT)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,M)=TOXINIT(L,K,NT)
       NLOW(LL,K,M)=0
       IF(NCSERW(LL,M).EQ.0) THEN
        TOX(L,K,NT)=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
        TOX1(L,K,NT)=TOX(L,K,NT)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,M)=TOXINIT(L,K,NT)
       NLOE(LL,K,M)=0
       IF(NCSERE(LL,3).EQ.0) THEN
        TOX(L,K,NT)=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
        TOX1(L,K,NT)=TOX(L,K,NT)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,M)=TOXINIT(L,K,NT)
       NLON(LL,K,M)=0
       IF(NCSERN(LL,M).EQ.0) THEN
        TOX(L,K,NT)=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
        TOX1(L,K,NT)=TOX(L,K,NT)
       END IF
       END DO
       END DO
C
      END IF
      END DO
      END IF
C
C **  INTIALIZE TOX BC IF NOT READ IN FROM RESTART FILE
C **  AND CONSTANT INITIAL CONDITIONS ARE USED
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(5).EQ.0) IISTMP=0
C
      IF (IISTMP.EQ.0.AND.ISTRAN(5).EQ.1) THEN
       DO NT=1,NTOX1
       IF(ITXINT(NT).EQ.0.OR.ITXINT(NT).EQ.2) THEN
       M=4+NT
C
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       NSID=NCSERS(LL,M)
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
       CLOS(LL,K,M)=TOXINTW(NT)
       NLOS(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       NSID=NCSERW(LL,M)
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       CLOW(LL,K,M)=TOXINTW(NT)
       NLOW(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       NSID=NCSERE(LL,M)
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       CLOE(LL,K,M)=TOXINTW(NT)
       NLOE(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       NSID=NCSERN(LL,M)
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
       CLON(LL,K,M)=TOXINTW(NT)
       NLON(LL,K,M)=0
       END DO
       END DO
C
      END IF
      END DO
      END IF
C
C **  INTIALIZE TOX BED IF NOT READ IN FROM RESTART FILE
C **  AND VARIABLE INITIAL CONDITIONS ARE USED
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(5).EQ.0) IISTMP=0
C
      IF (IISTMP.EQ.0.AND.ISTRAN(5).EQ.1) THEN
       DO NT=1,NTOX1
       IF(ITXINT(NT).EQ.2.OR.ITXINT(NT).EQ.3) THEN
C
       DO K=1,KB
       DO L=2,LA
       TOXB(L,K,NT)=TOXBINIT(L,K,NT)
       TOXB1(L,K,NT)=TOXB(L,K,NT)
       END DO
       END DO
C
       END IF
       END DO
      END IF
           
C
C **  INTIALIZE SED AND BC IF NOT READ IN FROM RESTART FILE
C **  AND VARIABLE INITIAL CONDITIONS ARE USED
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(6).EQ.0) IISTMP=0
C
      IF (IISTMP.EQ.0.AND.ISTRAN(6).EQ.1) THEN
       IF(ISEDINT.EQ.1.OR.ISEDINT.EQ.3) THEN
       DO NS=1,NSED
       M=4+NTOX1+NS

       DO K=1,KC
       DO L=2,LA
       SED(L,K,NS)=SEDINIT(L,K,NS)
       SED1(L,K,NS)=SED(L,K,NS)
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,M)=SEDINIT(L,K,NS)
       NLOS(LL,K,M)=0
       IF(NCSERS(LL,M).EQ.0) THEN
        SED(L,K,NS)=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
        SED1(L,K,NS)=SED(L,K,NS)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,M)=SEDINIT(L,K,NS)
       NLOW(LL,K,M)=0
       IF(NCSERW(LL,M).EQ.0) THEN
        SED(L,K,NS)=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
        SED1(L,K,NS)=SED(L,K,NS)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,M)=SEDINIT(L,K,NS)
       NLOE(LL,K,M)=0
       IF(NCSERE(LL,3).EQ.0) THEN
        SED(L,K,NS)=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
        SED1(L,K,NS)=SED(L,K,NS)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,M)=SEDINIT(L,K,NS)
       NLON(LL,K,M)=0
       IF(NCSERN(LL,M).EQ.0) THEN
        SED(L,K,NS)=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
        SED1(L,K,NS)=SED(L,K,NS)
       END IF
       END DO
       END DO
C
      END DO
      END IF
      END IF
C
C **  INTIALIZE SED BC IF NOT READ IN FROM RESTART FILE AND
C **  CONSTANT INITIAL CONDITIONS ARE USED
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(6).EQ.0) IISTMP=0
C
      IF (IISTMP.EQ.0.AND.ISTRAN(6).EQ.1) THEN
      IF (ISEDINT.EQ.0.OR.ISEDINT.EQ.2) THEN
      DO NS=1,NSED
       M=4+NTOX1+NS
C
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       NSID=NCSERS(LL,M)
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
       CLOS(LL,K,M)=SEDO(NS)
       NLOS(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       NSID=NCSERW(LL,M)
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       CLOW(LL,K,M)=SEDO(NS)
       NLOW(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       NSID=NCSERE(LL,M)
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       CLOE(LL,K,M)=SEDO(NS)
       NLOE(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       NSID=NCSERN(LL,M)
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
       CLON(LL,K,M)=SEDO(NS)
       NLON(LL,K,M)=0
       END DO
       END DO
C
      END DO
      END IF
      END IF
C
C **  INTIALIZE SED BED IF NOT READ IN FROM RESTART FILE
C **  AND VARIABLE INITIAL CONDITIONS ARE USED
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(6).EQ.0) IISTMP=0
C
      IF (IISTMP.EQ.0.AND.ISTRAN(6).EQ.1) THEN
       IF(ISEDINT.EQ.2.OR.ISEDINT.EQ.3) THEN
C
       DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
       SEDB(L,K,NS)=SEDBINIT(L,K,NS)
       SEDB1(L,K,NS)=SEDB(L,K,NS)
       END DO
       END DO
       END DO
C
       END IF
      END IF
C
C **  INTIALIZE SND AND BC IF NOT READ IN FROM RESTART FILE
C **  AND VARIABLE INITIAL CONDITIONS ARE USED
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(7).EQ.0) IISTMP=0
C
      IF (IISTMP.EQ.0.AND.ISTRAN(7).EQ.1) THEN
       IF(ISEDINT.EQ.1.OR.ISEDINT.EQ.3) THEN
       DO NS=1,NSND
       M=4+NTOX1+NSED+NS

       DO K=1,KC
       DO L=2,LA
       SND(L,K,NS)=SNDINIT(L,K,NS)
       SND1(L,K,NS)=SND(L,K,NS)
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,M)=SNDINIT(L,K,NS)
       NLOS(LL,K,M)=0
       IF(NCSERS(LL,M).EQ.0) THEN
        SND(L,K,NS)=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
        SND1(L,K,NS)=SND(L,K,NS)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,M)=SNDINIT(L,K,NS)
       NLOW(LL,K,M)=0
       IF(NCSERW(LL,M).EQ.0) THEN
        SND(L,K,NS)=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
        SND1(L,K,NS)=SND(L,K,NS)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,M)=SNDINIT(L,K,NS)
       NLOE(LL,K,M)=0
       IF(NCSERE(LL,3).EQ.0) THEN
        SND(L,K,NS)=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
        SND1(L,K,NS)=SND(L,K,NS)
       END IF
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,M)=SNDINIT(L,K,NS)
       NLON(LL,K,M)=0
       IF(NCSERN(LL,M).EQ.0) THEN
        SND(L,K,NS)=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
        SND1(L,K,NS)=SED(L,K,NS)
       END IF
       END DO
       END DO
C
      END DO
      END IF
      END IF
C
C **  INTIALIZE SND BC IF NOT READ IN FROM RESTART FILE AND
C **  CONSTANT INITIAL CONDITIONS ARE USED
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(7).EQ.0) IISTMP=0
C
      IF (IISTMP.EQ.0.AND.ISTRAN(7).EQ.1) THEN
      IF (ISEDINT.EQ.0.OR.ISEDINT.EQ.2) THEN
       DO NX=1,NSND
       NS=NSED+NX
       M=4+NTOX1+NSED+NX
C
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       NSID=NCSERS(LL,M)
       CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)
       CLOS(LL,K,M)=SEDO(NS)
       NLOS(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       NSID=NCSERW(LL,M)
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       CLOW(LL,K,M)=SEDO(NS)
       NLOW(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       NSID=NCSERE(LL,M)
       CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       CLOE(LL,K,M)=SEDO(NS)
       NLOE(LL,K,M)=0
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       NSID=NCSERN(LL,M)
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
       CLON(LL,K,M)=SEDO(NS)
       NLON(LL,K,M)=0
       END DO
       END DO
C
      END DO
      END IF
      END IF
C
C **  INTIALIZE SND BED IF NOT READ IN FROM RESTART FILE
C **  AND VARIABLE INITIAL CONDITIONS ARE USED
C
      IISTMP=1
      IF(ISRESTI.EQ.0) IISTMP=0
      IF(ISRESTI.GE.1.AND.ISCI(7).EQ.0) IISTMP=0
C
      IF (IISTMP.EQ.0.AND.ISTRAN(7).EQ.1) THEN
       IF(ISEDINT.EQ.2.OR.ISEDINT.EQ.3) THEN
C
       DO NX=1,NSND
       DO K=1,KB
       DO L=2,LA
       SNDB(L,K,NX)=SNDBINIT(L,K,NX)
       SNDB1(L,K,NX)=SNDB(L,K,NX)
       END DO
       END DO
       END DO
C
       END IF
      END IF
C
C**********************************************************************C
C
C **  INITIALIZE SEDIMENT BED
C
      IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1) CALL BEDINIT
c     IF(NSED.GE.1.OR.NSND.GE.1) CALL BEDINIT
C
C**********************************************************************C
C
C **  INITIALIZE BUOYANCY AND EQUATION OF STATE
C
      CALL CALBUOY
C
C**********************************************************************C
C
C **  INITIALIZE SFL IF (ISRESTI.EQ.0.AND ISTRAN(4).GE.1)
C
      IF (ISRESTI.EQ.0.AND.ISTRAN(4).GE.1) THEN
      IF (ISLTMT.EQ.0.AND.ISTOPT(4).EQ.11) THEN
C
       DO K=1,KC
       DO L=1,LC
       SFL(L,K)=SAL(L,K)
       SFL2(L,K)=SAL(L,K)
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBS
       L=LCBS(LL)
       CLOS(LL,K,5)=SALINIT(L,K)
       NLOS(LL,K,5)=0
       SFL(L,K)=WTCI(K,1)*CBS(LL,1,5)+WTCI(K,2)*CBS(LL,2,5)
       SFL2(L,K)=SFL(L,K)
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBW
       L=LCBW(LL)
       CLOW(LL,K,5)=SALINIT(L,K)
       NLOW(LL,K,5)=0
       SFL(L,K)=WTCI(K,1)*CBW(LL,1,5)+WTCI(K,2)*CBW(LL,2,5)
       SFL2(L,K)=SFL(L,K)
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBE
       L=LCBE(LL)
       CLOE(LL,K,5)=SALINIT(L,K)
       NLOE(LL,K,5)=0
       SFL(L,K)=WTCI(K,1)*CBE(LL,1,5)+WTCI(K,2)*CBE(LL,2,5)
       SFL2(L,K)=SFL(L,K)
       END DO
       END DO
C
       DO K=1,KC
       DO LL=1,NCBN
       L=LCBN(LL)
       CLON(LL,K,5)=SALINIT(L,K)
       NLON(LL,K,5)=0
       SFL(L,K)=WTCI(K,1)*CBN(LL,1,5)+WTCI(K,2)*CBN(LL,2,5)
       SFL2(L,K)=SFL(L,K)
       END DO
       END DO
C
      END IF
      END IF
C
C**********************************************************************C
C
C **  ACTIVATE DYE TRACER CONTINUITY CHECK
C
C----------------------------------------------------------------------C
C
      IF (ISMMC.EQ.1) THEN
C
      DO K=1,KC
      DO L=1,LC
      DYE(L,K)=1.
      DYE1(L,K)=1.
      END DO
      END DO
C
      DO K=1,KC
C
      DO LL=1,NCBS
      CLOS(LL,K,3)=1.
      NLOS(LL,K,3)=0
      END DO
C
      DO LL=1,NCBW
      CLOW(LL,K,3)=1.
      NLOW(LL,K,3)=0
      END DO
C
      DO LL=1,NCBE
      CLOE(LL,K,3)=1.
      NLOE(LL,K,3)=0
      END DO
C
      DO LL=1,NCBN
      CLON(LL,K,3)=1.
      NLON(LL,K,3)=0
      END DO
C
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  MASK CELLS TO BE CONVERTED FROM WATER TO LAND
C
	n=0
	fname='mask.inp'
      IF (ISMASK.EQ.1) CALL CELLMASK(fname)
C
C**********************************************************************C
C
C **  SET VERTICAL GRID DEPENDENT ARRAYS AND HARDWIRE DIMENSIONLESS
C **  MIXING LENGTH
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
C     DZC(K)=DZ
      DZIC(K)=1./DZC(K)
      END DO
C
      DZIG(0)=0.
      DZIG(KC)=0.
      DO K=1,KS
      DZG(K)=0.5*(DZC(K)+DZC(K+1))
      DZIG(K)=1./DZG(K)
      DZIGSD4(K)=0.25*DZIG(K)*DZIG(K)
      CDZU(K)=-DZC(K)/(DZC(K)+DZC(K+1))
      CDZL(K)=-DZC(K+1)/(DZC(K)+DZC(K+1))
      CDZF(K)=DZC(K)*DZC(K+1)/(DZC(K)+DZC(K+1))
      CDZM(K)=0.5*DZC(K)*DZC(K+1)
      END DO
C
      CDZR(1)=DZC(1)-1.
      CDZD(1)=DZC(1)
      DO K=2,KS
      CDZR(K)=DZC(K)+CDZR(K-1)
      CDZD(K)=DZC(K)+CDZD(K-1)
      END DO
C
      DO K=1,KS
      CDZR(K)=CDZR(K)*DZG(K)*CDZL(1)
      END DO
C
      CDZKMK(1)=0.
      DO K=2,KC
      CDZKMK(K)=DZIG(K-1)*DZIC(K)
      END DO
C
      DO K=1,KS
      CDZKK(K)=DZIC(K)*DZIG(K)
      CDZKKP(K)=DZIG(K)*DZIC(K+1)
      END DO
      CDZKK(KC)=0.
C
      Z(0)=0.
      IF(KC.GT.1) THEN
      DO K=1,KS
      Z(K)=Z(K-1)+DZC(K)
      ZZ(K)=Z(K)-0.5*DZC(K)
C     FPROX(K)=1./(VKC*Z(K)*(1.-Z(K)))**2
      FPROX(K)=(1./(VKC*Z(K))**2)+0.25*(1./(VKC*(1.-Z(K)))**2)/CTE2
      END DO
      END IF
C
      Z(KC)=Z(KS)+DZC(KC)
      ZZ(KC)=Z(KC)-0.5*DZC(KC)
C
      IF (ISRESTI.EQ.0) THEN
       DO K=0,KC
       DO L=1,LC
       DML(L,K)=VKC*Z(K)*(1.-Z(K))
       END DO
       END DO
      END IF
C
C**********************************************************************C
C
C **  INITIALIZE UNSTRETCHING PROCEDURE
C
C----------------------------------------------------------------------C
C
      DZPC=(SELVMAX-BELVMIN)/FLOAT(KPC)
C
      ZP(0)=BELVMIN
      DO KP=1,KPC
      ZP(KP)=ZP(KP-1)+DZPC
      END DO
C
      DO KP=1,KPC
      ZZP(KP)=0.5*(ZP(KP)+ZP(KP-1))
      END DO
C
      DO L=2,LA
      TMP=(BELV(L)-BELVMIN)/DZPC
      KPB(L)=NINT(0.5+TMP)
      END DO
C
C**********************************************************************C
C
C **  CALCULATE CONSTANT HORIZONTAL SPATIAL ARRAYS
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      DXYU(L)=DXU(L)*DYU(L)
      DXYV(L)=DXV(L)*DYV(L)
      DXYP(L)=STCAP(L)*DXP(L)*DYP(L)
      DXIU(L)=1./DXU(L)
      DYIU(L)=1./DYU(L)
      DXIV(L)=1./DXV(L)
      DYIV(L)=1./DYV(L)
      DXYIP(L)=1./(STCAP(L)*DXP(L)*DYP(L))
      DXYIU(L)=1./(DXU(L)*DYU(L))
      DXYIV(L)=1./(DXV(L)*DYV(L))
      HRU(L)=SUB(L)*HMU(L)*DYU(L)*DXIU(L)
      HRV(L)=SVB(L)*HMV(L)*DXV(L)*DYIV(L)
      HRUO(L)=SUB(L)*DYU(L)*DXIU(L)
      HRVO(L)=SVB(L)*DXV(L)*DYIV(L)
      SBX(L)=0.5*SBX(L)*DYU(L)
      SBY(L)=0.5*SBY(L)*DXV(L)
      SBXO(L)=SBX(L)
      SBYO(L)=SBY(L)
      SNLPX(L)=GID2*SNLPX(L)*DYU(L)
      SNLPY(L)=GID2*SNLPY(L)*DXV(L)
      END DO
C
C**********************************************************************C
C
C **  DEACTIVATE DRY CELLS
C
C     IF (ISDRY.GE.1.AND.ISDRY.LE.3) THEN
C       OPEN(1,FILE='drywet.log',ACCESS='APPEND',STATUS='UNKNOWN')
C       DO L=2,LA
C       IF (HP(L).LE.HDRY) THEN
C         LN=LNC(L)
C         NTMP=0
C         WRITE(1,6902)NTMP,IL(L),JL(L),HP(L),H1P(L),H2P(L)
C         WRITE(6,6902)NTMP,IL(L),JL(L),HP(L),H1P(L),H2P(L)
C         WRITE(8,6902)NTMP,IL(L),JL(L),HP(L),H1P(L),H2P(L)
C         NATDRY(L)=0
C         ISCDRY(L)=2
C         SUB(L)=0.
C         SUB(L+1)=0.
C         SVB(L)=0.
C         SVB(LN)=0.
C         SBX(L)=0.
C         SBX(L+1)=0.
C         SBY(L)=0.
C         SBY(LN)=0.
C       END IF
C       END DO
C       CLOSE(1)
C     END IF
C
 6902 FORMAT(1X,'DRYING AT N,I,J =',I10,2I6,'  H,H1,H2 =',3(2X,E12.4))
C
C**********************************************************************C
C
C **  CHECK FOR DRYING AND SET SWITCHES ON RESTART
C
C----------------------------------------------------------------------C
C
      IF (ISDRY.GE.1.AND.ISRESTI.EQ.1) THEN
      OPEN(1,FILE='drywet.log',ACCESS='APPEND',STATUS='UNKNOWN')
C
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      HUTMPP=0.5*(HP(L)+HP(L-1))
      IF (HUTMPP.LE.HUWET(L)) THEN
        SUB(L)=0.
        SBX(L)=0.
      END IF
      HUTMPP=0.5*(HP(L)+HP(L+1))
      IF (HUTMPP.LE.HUWET(L+1)) THEN
        SUB(L+1)=0.
        SBX(L+1)=0.
      END IF
      HVTMPP=0.5*(HP(L)+HP(LS))
      IF (HVTMPP.LE.HVWET(L)) THEN
        SVB(L)=0.
        SBY(L)=0.
      END IF
      HVTMPP=0.5*(HP(L)+HP(LN))
      IF (HVTMPP.LE.HVWET(LN)) THEN
        SVB(LN)=0.
        SBY(LN)=0.
      END IF
      IF (HP(L).LE.HDRY) THEN
        SUB(L)=0.
        SVB(L)=0.
        SUB(L+1)=0.
        SVB(LN)=0.
        SBX(L)=0.
        SBY(L)=0.
        SBX(L+1)=0.
        SBY(LN)=0.
      END IF
      END DO
C
      CLOSE(1)
      END IF
C
C**********************************************************************C
C
C **  INITIALIZE ZERO DIMENSION VOLUME BALANCE
C
C----------------------------------------------------------------------C
C
      IF (ISDRY.GE.1) THEN
        OPEN(1,FILE='zvolbal.out',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='avsel.out',STATUS='UNKNOWN')
        LPBTMP=0
        DO L=2,LA
        ISLUSED(L)=0
        IF (SPB(L).EQ.0) THEN
          LPBTMP=LPBTMP+1
          ISLUSED(L)=1
        END IF
        LORDER(L)=0
        END DO
        ISLUSED(1)=1
        ISLUSED(LC)=1
        LORMAX=LC-2-LPBTMP
        DO LS=1,LORMAX
        BELMIN=100000.
          DO L=2,LA
           IF (SPB(L).NE.0.AND.ISLUSED(L).NE.1) THEN
             IF (BELV(L).LT.BELMIN) THEN
               LBELMIN=L
               BELMIN=BELV(L)
             END IF
           END IF
           END DO
         LORDER(LS)=LBELMIN
         ISLUSED(LBELMIN)=1
        END DO
        WRITE(1,5300)
        LS=1
        L=LORDER(LS)
        BELSURF(LS)=BELV(L)
        ASURFEL(LS)=DXYP(L)
        VOLSEL(LS)=0.
        WRITE(1,5301)LS,BELSURF(LS),ASURFEL(LS),VOLSEL(LS)
        DO LS=2,LORMAX
          L=LORDER(LS)
          BELSURF(LS)=BELV(L)
          ASURFEL(LS)=ASURFEL(LS-1)+DXYP(L)
          VOLSEL(LS)=VOLSEL(LS-1)+0.5*(BELSURF(LS)-BELSURF(LS-1))*
     $                                (ASURFEL(LS)+ASURFEL(LS-1))
        WRITE(1,5301)LS,BELSURF(LS),ASURFEL(LS),VOLSEL(LS)
        END DO
        LS=LORMAX+1
        BELSURF(LS)=BELV(L)+10.0
        ASURFEL(LS)=ASURFEL(LS-1)
        VOLSEL(LS)=VOLSEL(LS-1)+0.5*(BELSURF(LS)-BELSURF(LS-1))*
     $                                (ASURFEL(LS)+ASURFEL(LS-1))
        WRITE(1,5301)LS,BELSURF(LS),ASURFEL(LS),VOLSEL(LS)
        VOLZERD=0.
        VOLLDRY=0.
        DO L=2,LA
        IF (SPB(L).NE.0) THEN
          VOLZERD=VOLZERD+DXYP(L)*HP(L)
          IF (HP(L).GT.HDRY) VOLLDRY=VOLLDRY+DXYP(L)*HP(L)
        END IF
        END DO
        DO LS=1,LORMAX
        IF (VOLZERD.GE.VOLSEL(LS).AND.VOLZERD.LT.VOLSEL(LS+1)) THEN
           WTM=VOLSEL(LS+1)-VOLZERD
           WTMP=VOLZERD-VOLSEL(LS)
           DELVOL=VOLSEL(LS+1)-VOLSEL(LS)
           WTM=WTM/DELVOL
           WTMP=WTMP/DELVOL
           SELZERD=WTM*BELSURF(LS)+WTMP*BELSURF(LS+1)
           ASFZERD=WTM*ASURFEL(LS)+WTMP*ASURFEL(LS+1)
        END IF
        END DO
        VETZERD=VOLZERD
        WRITE(1,5302)
        WRITE(1,5303) SELZERD,ASFZERD,VOLZERD,VOLLDRY
        CLOSE(1)
      END IF
C
 5300 FORMAT('   M    BELSURF     ASURFEL     ',
     $       '   VOLSEL',/)
 5301 FORMAT(1X,I3,2X,F10.5,2X,E12.4,2X,E12.4)
 5302 FORMAT(/)
 5303 FORMAT(2X,F10.5,3(2X,E12.4))
C
C**********************************************************************C
C
C **  INITIALIZE ELEVATION OF ACTIVE GROUNDWATER ZONE FOR COLD START
C
C----------------------------------------------------------------------C
C
c       WRITE(6,5301) ISGWIE,DAGWZ,RNPOR,RIFTRM
C
        IF (ISGWIE.GE.1.AND.ISRESTI.EQ.0) THEN
C
        DO L=2,LA
          IF (HP(L).GT.HDRY) THEN
            AGWELV(L)=BELV(L)
           ELSE
            IF (BELAGW(L).LT.SELZERD) THEN
              AGWELV(L)=SELZERD
              AGWELV(L)=MIN(AGWELV(L),BELV(L))
             ELSE
              AGWELV(L)=BELAGW(L)
            END IF
          END IF
        END DO
C
        DO L=2,LA
        AGWELV1(L)=AGWELV(L)
        AGWELV2(L)=AGWELV(L)
        END DO
C
        OPEN(1,FILE='gwelv.out',STATUS='UNKNOWN')
        WRITE(1,5400)
        WRITE(1,5402)
        DO L=2,LA
        WRITE(1,5401)IL(L),JL(L),BELV(L),BELAGW(L),AGWELV(L)
        END DO
        CLOSE(1)
C
        END IF
C
C
 5400 FORMAT('   I   J    BELELV      BELAGW     ',
     $       '   AGWELV',/)
 5401 FORMAT(1X,2I5,2X,F10.5,2X,F10.5,2X,F10.5)
 5402 FORMAT(/)
C
C**********************************************************************C
C
C **  CALCULATE CONSTANT C ARRAYS FOR EXTERNAL P SOLUTION
C **  HRU=SUB*HMU*DYU/DXU & HRV=SVB*HMV*DXV/DYV
C **  DXYIP=1/(DXP*DYP)
C
C----------------------------------------------------------------------C
C
      IF (ISLTMT.EQ.0) THEN
C
C----------------------------------------------------------------------C
C
      IF(IRVEC.NE.9) THEN
C
      DO L=2,LA
      CC(L)=1.
      CCC(L)=1.
      END DO
C
      IF(ISRLID.EQ.1) THEN
       DO L=2,LA
       CC(L)=0.
       CCC(L)=0.
       IF(SPB(L).EQ.0.) CC(L)=1.
       IF(SPB(L).EQ.0.) CCC(L)=1.
       END DO
      END IF
C
      DO L=2,LA
      LN=LNC(L)
      C1=-G*DT*DT*SPB(L)*DXYIP(L)
      CS(L)=C1*HRV(L)
      CW(L)=C1*HRU(L)
      CE(L)=C1*HRU(L+1)
      CN(L)=C1*HRV(LN)
      CC(L)=CC(L)-CS(L)-CW(L)-CE(L)-CN(L)
      CCI(L)=1./CC(L)
      CCS(L)=0.25*CS(L)
      CCW(L)=0.25*CW(L)
      CCE(L)=0.25*CE(L)
      CCN(L)=0.25*CN(L)
      CCC(L)=CCC(L)-CCS(L)-CCW(L)-CCE(L)-CCN(L)
      CCCI(L)=1./CCC(L)
      END DO
C
      DO LR=1,NRC
      L=LRC(LR)
      CCSR(LR)=CCS(L)*CCCI(L)
      CCWR(LR)=CCW(L)*CCCI(L)
      CCER(LR)=CCE(L)*CCCI(L)
      CCNR(LR)=CCN(L)*CCCI(L)
      CSR(LR)=CS(L)*CCI(L)
      CWR(LR)=CW(L)*CCI(L)
      CER(LR)=CE(L)*CCI(L)
      CNR(LR)=CN(L)*CCI(L)
      END DO
C
      DO LB=1,NBC
      L=LBC(LB)
      CCSB(LB)=CCS(L)*CCCI(L)
      CCWB(LB)=CCW(L)*CCCI(L)
      CCEB(LB)=CCE(L)*CCCI(L)
      CCNB(LB)=CCN(L)*CCCI(L)
      CSB(LB)=CS(L)*CCI(L)
      CWB(LB)=CW(L)*CCI(L)
      CEB(LB)=CE(L)*CCI(L)
      CNB(LB)=CN(L)*CCI(L)
      END DO
C
      DO LL=1,NPBW
      L=LPBW(LL)
      IF(ISPBW(LL).EQ.0)THEN
       IF (ISRED(L).EQ.1) THEN
        LR=LLRC(L)
        CER(LR)=0.
        CCER(LR)=0.
       ELSE
        LB=LLBC(L)
        CEB(LB)=0.
        CCEB(LB)=0.
       END IF
      END IF
      END DO
C
      DO LL=1,NPBE
      L=LPBE(LL)
      IF(ISPBE(LL).EQ.0)THEN
       IF (ISRED(L).EQ.1) THEN
        LR=LLRC(L)
        CWR(LR)=0.
        CCWR(LR)=0.
       ELSE
        LB=LLBC(L)
        CWB(LB)=0.
        CCWB(LB)=0.
       END IF
      END IF
      END DO
C
      DO LL=1,NPBS
      L=LPBS(LL)
      IF(ISPBS(LL).EQ.0)THEN
       IF (ISRED(L).EQ.1) THEN
        LR=LLRC(L)
        CNR(LR)=0.
        CCNR(LR)=0.
       ELSE
        LB=LLBC(L)
        CNB(LB)=0.
        CCNB(LB)=0.
       END IF
      END IF
      END DO
C
      DO LL=1,NPBN
      L=LPBN(LL)
      IF(ISPBN(LL).EQ.0)THEN
       IF (ISRED(L).EQ.1) THEN
        LR=LLRC(L)
        CSR(LR)=0.
        CCSR(LR)=0.
       ELSE
        LB=LLBC(L)
        CSB(LB)=0.
        CCSB(LB)=0.
       END IF
      END IF
      END DO
C
      END IF
C
C----------------------------------------------------------------------C
C
      ELSE
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LE=L+1
      LN=LNC(L)
      C1=-SPB(L)/(HRV(L)+HRU(L)+HRU(LE)+HRV(LN))
      CS(L)=C1*HRV(L)
      CW(L)=C1*HRU(L)
      CE(L)=C1*HRU(LE)
      CN(L)=C1*HRV(LN)
      CC(L)=C1
      END DO
C
      DO LL=1,NPBE
      L=LPBE(LL)-1
      LE=L+1
      LN=LNC(L)
      C1=(HRV(L)+HRU(L)+HRV(LN))
      IF (C1.EQ.0.) C1=1.
      C1=-SPB(L)/C1
      CS(L)=C1*HRV(L)
      CW(L)=C1*HRU(L)
      CE(L)=0.
      CN(L)=C1*HRV(LN)
      CC(L)=C1
      END DO
C
      DO LL=1,NPBW
      L=LPBW(LL)+1
      LE=L+1
      LN=LNC(L)
      C1=(HRV(L)+HRU(LE)+HRV(LN))
      IF (C1.EQ.0.) C1=1.
      C1=-SPB(L)/C1
      CS(L)=C1*HRV(L)
      CW(L)=0.
      CE(L)=C1*HRU(LE)
      CN(L)=C1*HRV(LN)
      CC(L)=C1
      END DO
C
      DO LL=1,NPBN
      L=LPBN(LL)
      L=LSC(L)
      LE=L+1
      LN=LNC(L)
      C1=(HRV(L)+HRU(L)+HRU(LE))
      IF (C1.EQ.0.) C1=1.
      C1=-SPB(L)/C1
      CS(L)=C1*HRV(L)
      CW(L)=C1*HRU(L)
      CE(L)=C1*HRU(LE)
      CN(L)=0.
      CC(L)=C1
      END DO
C
      DO LL=1,NPBS
      L=LPBS(LL)
      L=LNC(L)
      LE=L+1
      LN=LNC(L)
      C1=(HRU(L)+HRU(LE)+HRV(LN))
      IF (C1.EQ.0.) C1=1.
      C1=-SPB(L)/C1
      CS(L)=0.
      CW(L)=C1*HRU(L)
      CE(L)=C1*HRU(LE)
      CN(L)=C1*HRV(LN)
      CC(L)=C1
      END DO
C
C----------------------------------------------------------------------C
C
      END IF
C
C**********************************************************************C
C
C **  SMOOTH INITIAL SALINITY
C
C----------------------------------------------------------------------C
C
      IF (NSBMAX.GE.1) THEN
      CALL SALTSMTH
      END IF
C
C**********************************************************************C
C
C **  OUTPUT INITIAL DEPTH AND SALINITY FIELDS
C
C----------------------------------------------------------------------C
C
C **  PLOT SMOOTHED CELL CENTER STATIC DEPTHS
C
      DO L=2,LA
      PAM(L)=HMP(L)
      END DO
      WRITE (7,16)
      CALL PPLOT (2)
C
C----------------------------------------------------------------------C
C
      CALL DEPPLT
C
C----------------------------------------------------------------------C
C
C **  PLOT INITIAL SALINITY IN SURFACE AND BOTTOM LAYERS
C
      DO KK=1,KC,KS
      DO L=2,LA
      PAM(L)=SAL(L,KK)
      END DO
      WRITE (7,316) KK
      CALL PPLOT (1)
      END DO
C
C----------------------------------------------------------------------C
C
   16 FORMAT (1H1,' CELL CENTER STATIC DEPTHS',//)
  316 FORMAT (1H1,'INITIAL SALINITY IN LAYER',I5,//)
C
C**********************************************************************C
C
C **  INITIALIZE WATER QUALITY MODEL AND READ INPUT
C
!      IF(ISTRAN(8).GE.1) CALL WQ3DINP(LA,KC,IC,JC,IWQDT,DT,NTSPTC)
C
C**********************************************************************C
C
C **  SELECT FULL HYDRODYNAMIC AND MASS TRANSPORT CALCULATION OR
C **  LONG-TERM MASS TRANSPORT CALCULATION
C
      IF (ISLTMT.EQ.0) THEN
         IF(IS1DCHAN.EQ.0) CALL HDMT
         IF(IS1DCHAN.GE.1) CALL HDMT1D
      END IF
      IF (ISLTMT.GE.1) CALL LTMT
C
C**********************************************************************C
C
C **  CLOSE OUTPUT  FILES
C
C----------------------------------------------------------------------C
C
      CLOSE(7)
      CLOSE(8)
      CLOSE(9)
C
C**********************************************************************C
C
      STOP
      END
