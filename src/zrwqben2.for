C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RWQBEN2 (timtmp,LA)
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 10 april 1999
C
C**********************************************************************C
C Read in spatially and/or temporally varying parameters for BENTHIC
C FLUXES of PO4d, NH4, NO3, SAd, COD, O2
C**********************************************************************C
C **  LAST MODIFIED BY Mike Morton ON 30 Jan 1998
c M. Morton 01/30/98: changed code to allow for temporally
c   varying benthic fluxes in the BENFN file.  Previous version only
c   provided spatially varying flux (no provision for time varying).
C**********************************************************************C
c Format of BENFN file is:
c    Title 1
c    Title 2
c    Title 3
c  270.00000  <-- day at which following fluxes become active
c     1  FPO4  FNH4  FNO3  FSAD  FCOD  FSOD <-- fluxes for zone 1
c     2  FPO4  FNH4  FNO3  FSAD  FCOD  FSOD <-- fluxes for zone 2
c     etc. for IWQZ zones
c  350.00000  <-- day at which following fluxes become active
c     1  FPO4  FNH4  FNO3  FSAD  FCOD  FSOD <-- fluxes for zone 1
c     2  FPO4  FNH4  FNO3  FSAD  FCOD  FSOD <-- fluxes for zone 2
c     etc. for IWQZ zones
c 9999.99999 <-- enter large day at end of file
C**********************************************************************C
C
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
c
      DIMENSION XBFPO4D(NWQZM),XBFNH4(NWQZM),XBFNO3(NWQZM)
      DIMENSION XBFCOD(NWQZM),XBFO2(NWQZM),XBFSAD(NWQZM)
      dimension izone(nwqzm)
      CHARACTER TITLE(3)*79, dumrm*1
C
      OPEN(1,FILE=benfn,STATUS='UNKNOWN')
      OPEN(2,FILE='wq3d.out',STATUS='UNKNOWN',ACCESS='APPEND')
C
c Skip over three header records:
      READ(1,50) (TITLE(M),M=1,3)
      WRITE(2,999)
      WRITE(2,50) (TITLE(M),M=1,3)
c skip over 6 more comment cards
      do m = 1,6
        read(1,5) dumrm
      end do
5     format(a1)
      dumrm=dumrm

      read(1, *) ibenz
      write(2, 65) timtmp, ibenz
65    format(' * Benthic fluxes at     ', f10.5,' days of model run',/,
     +       ' Number of benthic flux zones = ', i4)
c
c Sequentially read through benthic flux file until the appropriate
c time is found:
c
c   bday   = current day at which benthic flux is in effect
c   benday = next day at which benthic flux changes (passed to main program)
c
10    read(1, *, end=15) benday
      if (benday .gt. timtmp) go to 20
      bday = benday
      DO I=1,ibenz
        READ(1,*,END=15) MM, XBFPO4D(mm), XBFNH4(mm), XBFNO3(mm),
     *    XBFSAD(mm), XBFCOD(mm), XBFO2(mm)
        izone(i) = mm
      end do
      go to 10

c Unexpected end-of-file encountered:
15    write(2,16) benfn
16    format(//,' ************* WARNING *************',/,
     +          ' End-of-file encountered in file: ', A20,/,/
     +          ' Benthic fluxes set to values corresponding ',
     +          ' to LAST DAY in file.',/)

20    continue

      write(2, 48) bday
48    format(/,' Day in Benthic Flux File: ',F10.5,/,
     + '    Zone    FPO4    FNH4    FNO3    FSAD    FCOD    FSOD')
      do i=1,ibenz
        mm = izone(i)
        WRITE(2,51) mm,XBFPO4D(mm),XBFNH4(mm),XBFNO3(mm),XBFSAD(mm),
     *    XBFCOD(mm),XBFO2(mm)
      END DO
C
c Determine benthic flux for each cell (L) by interpolating between
c the mud and sand fluxes.  XBENMUD(L) is the percent mud for each
c cell.
      DO L=2,LA
        izm = ibenmap(L,1)
        izs = ibenmap(L,2)
        xm = xbenmud(L)
        WQBFPO4D(L) = xm*XBFPO4D(izm) + (1.0-xm)*XBFPO4D(izs)
        WQBFNH4(L)  = xm*XBFNH4(izm)  + (1.0-xm)*XBFNH4(izs)
        WQBFNO3(L)  = xm*XBFNO3(izm)  + (1.0-xm)*XBFNO3(izs)
        WQBFSAD(L)  = xm*XBFSAD(izm)  + (1.0-xm)*XBFSAD(izs)
        WQBFCOD(L)  = xm*XBFCOD(izm)  + (1.0-xm)*XBFCOD(izs)
        WQBFO2(L)   = xm*XBFO2(izm)   + (1.0-xm)*XBFO2(izs)
      END DO
C
      close(1)
      CLOSE(2)
C
  999 FORMAT(1X)
   50 FORMAT(A79)
   51 FORMAT(I8, 10F8.3)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)
C
      RETURN
      END
