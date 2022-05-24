C
C*********************************************************
C*********************************************************
C*********************************************************
C
      subroutine INITbin3(TBEGIN, DT,LA,KC)
c
c M.R. Morton    29 Apr 1999
c
c Initializes binary file for EFDC output.  Places control
c parameters for post-processor in header section of binary
c file WQDOCOMP.BIN for D.O. component analysis.
c
c-------------------------------------------------------------------
c
      include 'wq.par'
      include 'wqcom.cmn'

      PARAMETER(mxparm=31)
      real tend
	integer nparm, ncells
      logical fexist, is1open, is2open
      character*20 wqname(mxparm)
      character*10 wqunits(mxparm)
      character*3  wqcode(mxparm)

c
c The following parameters are specified in EFDC.INP and WQ3DWC.INP:
c KC       = number of vertical layers
c IWQTSDT  = number of time steps per data dump
c DT       = time step of EFDC model in seconds
c LA       = number of active cells + 1 in model
c TBEGIN   = beginning time of run in days
c
c The parameter NPARM must be changed if the output data
c is changed in SUBROUTINE WWQTSbin:
c nparm   = number of parameters written to binary file
c
c nrec3   = number of records written to binary file (one record
c           is a complete data dump for time interval IWQDIUDT)
c
      nparm = 31
      ncells = LA-1
      nrec3 = 0
      tend = TBEGIN
c
c The following water quality names, units, and 3-character codes
c should be modified to match the parameters written to the binary
c file in SUBROUTINE WWQTSbin.  The character strings must be
c exactly the length specified below in order for the post-processor
c to work correctly.
c
c Be sure WQNAME strings are exactly 20-characters long:
c------------------'         1         2'
c------------------'12345678901234567890'
      wqname( 1) = 'Nitrogen_Limit_CYA  '
      wqname( 2) = 'Nitrogen_Limit_DIA  '
      wqname( 3) = 'Nitrogen_Limit_GRN  '
      wqname( 4) = 'Nitrogen_Limit_MAC  '
      wqname( 5) = 'Phosphorus_Limit_CYA'
      wqname( 6) = 'Phosphorus_Limit_DIA'
      wqname( 7) = 'Phosphorus_Limit_GRN'
      wqname( 8) = 'Phosphorus_Limit_MAC'
      wqname( 9) = 'Light_Limit_CYA     '
      wqname(10) = 'Light_Limit_DIA     '
      wqname(11) = 'Light_Limit_GRN     '
      wqname(12) = 'Light_Limit_MAC     '
      wqname(13) = 'Temp_Limit_CYA      '
      wqname(14) = 'Temp_Limit_DIA      '
      wqname(15) = 'Temp_Limit_GRN      '
      wqname(16) = 'Temp_Limit_MAC      '
      wqname(17) = 'DO_saturation       '
      wqname(18) = 'DO_point_sources    '
      wqname(19) = 'DO_sed_oxygen_demand'
      wqname(20) = 'DO_reaeration       '
      wqname(21) = 'DO_DOC_decay        '
      wqname(22) = 'DO_NH4_nitrification'
      wqname(23) = 'DO_COD_oxidation    '
      wqname(24) = 'DO_photosynth_CHL   '
      wqname(25) = 'DO_respiration_CHL  '
      wqname(26) = 'DO_photosynth_MAC   '
      wqname(27) = 'DO_respiration_MAC  '
      wqname(28) = 'DO_deficit          '
      wqname(29) = 'DO_transport        '
      wqname(30) = 'DO_ALL_COMPONENTS   '
      wqname(31) = 'Layer_Thickness     '
c
c Be sure WQUNITS strings are exactly 10-characters long:
c-------------------'         1'
c-------------------'1234567890'
      wqunits( 1) = 'unitless  '
      wqunits( 2) = 'unitless  '
      wqunits( 3) = 'unitless  '
      wqunits( 4) = 'unitless  '
      wqunits( 5) = 'unitless  '
      wqunits( 6) = 'unitless  '
      wqunits( 7) = 'unitless  '
      wqunits( 8) = 'unitless  '
      wqunits( 9) = 'unitless  '
      wqunits(10) = 'unitless  '
      wqunits(11) = 'unitless  '
      wqunits(12) = 'unitless  '
      wqunits(13) = 'unitless  '
      wqunits(14) = 'unitless  '
      wqunits(15) = 'unitless  '
      wqunits(16) = 'unitless  '
      wqunits(17) = 'mg/L/day  '
      wqunits(18) = 'mg/L/day  '
      wqunits(19) = 'mg/L/day  '
      wqunits(20) = 'mg/L/day  '
      wqunits(21) = 'mg/L/day  '
      wqunits(22) = 'mg/L/day  '
      wqunits(23) = 'mg/L/day  '
      wqunits(24) = 'mg/L/day  '
      wqunits(25) = 'mg/L/day  '
      wqunits(26) = 'mg/L/day  '
      wqunits(27) = 'mg/L/day  '
      wqunits(28) = 'mg/L/day  '
      wqunits(29) = 'mg/L/day  '
      wqunits(30) = 'mg/L/day  '
      wqunits(31) = 'meters    '
c
c Be sure WQCODE strings are exactly 3-characters long:
c
c------------------'123'
      wqcode( 1) = 'NLC'
      wqcode( 2) = 'NLD'
      wqcode( 3) = 'NLG'
      wqcode( 4) = 'NLM'
      wqcode( 5) = 'PLC'
      wqcode( 6) = 'PLD'
      wqcode( 7) = 'PLG'
      wqcode( 8) = 'PLM'
      wqcode( 9) = 'LLC'
      wqcode(10) = 'LLD'
      wqcode(11) = 'LLG'
      wqcode(12) = 'LLM'
      wqcode(13) = 'TLC'
      wqcode(14) = 'TLD'
      wqcode(15) = 'TLG'
      wqcode(16) = 'TLM'
      wqcode(17) = 'DCS'
      wqcode(18) = 'DPS'
      wqcode(19) = 'DSO'
      wqcode(20) = 'DKA'
      wqcode(21) = 'DCA'
      wqcode(22) = 'DNH'
      wqcode(23) = 'DCO'
      wqcode(24) = 'DPC'
      wqcode(25) = 'DRC'
      wqcode(26) = 'DPM'
      wqcode(27) = 'DRM'
      wqcode(28) = 'DEF'
      wqcode(29) = 'DTR'
      wqcode(30) = 'DAL'
      wqcode(31) = 'DZZ'
c
c-------------------------------------------------------------------
c
c If WQDOCOMP.BIN already exists, delete it here.
c
      if (isCOMP .gt. 0) then
          OPEN(84, FILE='WQDOCOMP.out',STATUS='unknown')
          close(84, STATUS='DELETE') 
          
        io = 1
10      io = io+1
        if (io .gt. 99) then
          write(0,*) ' No available IO units ... io > 99'
          stop ' EFDC Halted in Subroutine INITbin3'
        end if
        INQUIRE(UNIT=io, OPENED=is2open)
        if (is2open) go to 10
        INQUIRE(FILE='WQDOCOMP.bin', EXIST=fexist)
        if (fexist) then
          OPEN(UNIT=io, FILE='WQDOCOMP.bin', ACCESS='append',
     +      FORM='binary', STATUS='unknown')
          CLOSE(UNIT=io, STATUS='DELETE')
          
          write(0,*) 'Old file WQDOCOMP.BIN/out deleted...'
        end if
        io=82
        open(UNIT=io, FILE='WQDOCOMP.bin', ACCESS='append',
c         open(UNIT=io, FILE='WQDOCOMP.bin',
     +     FORM='binary', STATUS='unknown')
c--------------------------------------------------------------------
c write control parameters for post-processor to header
c section of the WQDOCOMP.bin binary file:
c
        write(io) LCMWQ,ICMWQ,JCMWQ,KCWM,TBEGIN, DT, IWQTSDT,KC 
        write(io) wqname
        write(io) wqunits
        write(io) wqcode
c
c write cell I,J mapping reference to header section of binary file:
c
        write(io) ILW
        write(io) JLW
        write(io) LIJW
	  call flush(io)
	  close(io)
	  
	  open(84, FILE='WQDOCOMP.out', ACCESS='append')
        do i=1,31
        write(84,'(I4,2x,A20,2x,A10,2x,A3)')i,wqname(i),wqunits(i),
     +  wqcode(i)
        enddo
        close(84)
        
      endif
      return
      end
