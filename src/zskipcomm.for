C
C*******************************************************************************
C*******************************************************************************
C*******************************************************************************
C
      subroutine SkipComm(iunit, cc)
c
c Skips over comment lines in input files
c
      integer iunit
      character cc*1, line*80
100   read(iunit, *, end=999) line
      if (line(1:1) .eq. cc) go to 100
      if (line(1:1) .eq. 'C') go to 100
      if (line(1:1) .eq. 'c') go to 100
      backspace(iunit)
999   return
      end
