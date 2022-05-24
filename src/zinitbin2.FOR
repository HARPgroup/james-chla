C
C*********************************************************
C*********************************************************
C*********************************************************
C
      subroutine INITbin2
c
c M.R. Morton    01 Apr 1999
c
c Initializes binary file for EFDC output.  Places control
c parameters for post-processor in header section of binary
c file WQDIURDO.BIN for diurnal dissolved oxygen calculations.
c
c---------------------------------------------------------
c
      include 'wq.par'
      include 'wqcom.cmn'

      PARAMETER(mxparm=30)
      real tend
      integer nparm, ncells
      logical fexist
      character*20 wqname(mxparm)
      character*10 wqunits(mxparm)
      character*3  wqcode(mxparm)

c
c The following parameters are specified in EFDC.INP and WQ3DWC.INP:
c KC       = number of vertical layers
c IWQDIUDT = number of time steps per data dump
c DT       = time step of EFDC model in seconds
c LA       = number of active cells + 1 in model
c TBEGIN   = beginning time of run in days
c
c The parameter NPARM must be changed if the output data
c is changed in SUBROUTINE WWQTSbin:
c nparm   = number of parameters written to binary file
c
c nrec2   = number of records written to binary file (one record
c           is a complete data dump for time interval IWQDIUDT)
c
      nparm = 20
!      ncells = LA-1
      nrec2 = 0
!      tend = TBEGIN
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
      wqname( 1) = 'Salinity            '
      wqname( 2) = 'Temperature         '
      wqname( 3) = 'DO_Saturation       '
      wqname( 4) = 'Dissolved_Oxygen    '
      wqname( 5) = 'Sed_Oxygen_Demand   '
      wqname( 6) = 'Ka_reaeration       '
      wqname( 7) = 'Layer_Thickness     '
      wqname( 8) = 'Prim_Production_CYA '
      wqname( 9) = 'Respiration_CYA     '
      wqname(10) = 'Prim_Production_DIA '
      wqname(11) = 'Respiration_DIA     '
      wqname(12) = 'Prim_Production_GRN '
      wqname(13) = 'Respiration_GRN     '
      wqname(14) = 'Prim_Production_MAC '
      wqname(15) = 'Respiration_MAC     '
      wqname(16) = 'Algae_Cyanobacteria '
      wqname(17) = 'Algae_Diatoms       '
      wqname(18) = 'Algae_Greens        '
      wqname(19) = 'Total_Chlorophylla  '
      wqname(20) = 'Macroalgae          '
c
c Be sure WQUNITS strings are exactly 10-characters long:
c-------------------'         1'
c-------------------'1234567890'
      wqunits( 1) = 'g/L       '
      wqunits( 2) = 'degC      '
      wqunits( 3) = 'mg/L      '
      wqunits( 4) = 'mg/L      '
      wqunits( 5) = 'g/m2/day  '
      wqunits( 6) = 'per_day   '
      wqunits( 7) = 'meters    '
      wqunits( 8) = 'mgO2/L/day'
      wqunits( 9) = 'mgO2/L/day'
      wqunits(10) = 'mgO2/L/day'
      wqunits(11) = 'mgO2/L/day'
      wqunits(12) = 'mgO2/L/day'
      wqunits(13) = 'mgO2/L/day'
      wqunits(14) = 'mgO2/L/day'
      wqunits(15) = 'mgO2/L/day'
      wqunits(16) = 'ug/L      '
      wqunits(17) = 'ug/L      '
      wqunits(18) = 'ug/L      '
      wqunits(19) = 'ug/L      '
      wqunits(20) = 'ug/L      '
c
c Be sure WQCODE strings are exactly 3-characters long:
c
c------------------'123'
      wqcode( 1) = 'SAL'
      wqcode( 2) = 'TEM'
      wqcode( 3) = 'DOS'
      wqcode( 4) = 'DOO'
      wqcode( 5) = 'SOD'
      wqcode( 6) = 'RKA'
      wqcode( 7) = 'DEP'
      wqcode( 8) = 'PPC'
      wqcode( 9) = 'RRC'
      wqcode(10) = 'PPD'
      wqcode(11) = 'RRD'
      wqcode(12) = 'PPG'
      wqcode(13) = 'RRG'
      wqcode(14) = 'PPM'
      wqcode(15) = 'RRM'
      wqcode(16) = 'CYA'
      wqcode(17) = 'DIA'
      wqcode(18) = 'GRN'
      wqcode(19) = 'CHL'
      wqcode(20) = 'MAC'
c
c-------------------------------------------------------------------
c
c If WQDIURDO.BIN already exists, delete it here.
c
      if (isDIURDO .gt. 0) then
        INQUIRE(FILE='WQDIURDO.bin', EXIST=fexist)
        if (fexist) then
          OPEN(UNIT=2, FILE='WQDIURDO.bin', ACCESS='append',
     +      FORM='unformatted', STATUS='unknown')
          CLOSE(UNIT=2, STATUS='DELETE')
          write(0,*) 'Old file WQDIURDO.BIN deleted...'
        end if

        open(UNIT=2, FILE='WQDIURDO.bin', ACCESS='append',
     +     FORM='unformatted', STATUS='unknown')
c--------------------------------------------------------------------
c write control parameters for post-processor to header
c section of the WQWCavg.bin binary file:
c
!        write(2) tend, DT, IWQDIUDT, nparm, ncells, KC
        write(2) (wqname(i), i=1,nparm)
        write(2) (wqunits(i), i=1,nparm)
        write(2) (wqcode(i), i=1,nparm)
c
c write cell I,J mapping reference to header section of binary file:
c
          write(2) ILW
          write(2) JLW
       close(2)
      end if

      return
      end
