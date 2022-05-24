C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      REAL FUNCTION CSNDSET(SND,SDEN,IOPT)
C
C **  CALCULATES HINDERED SETTLING CORRECTION FOR CLASS NS NONCOHESIVE 
C **  SEDIMENT 
C
c      ROPT=FLOAT(IOTP)  !bug, Ji, 3/5/00
      ROPT=FLOAT(IOPT)
      CSNDSET=(1.-SDEN*SND)**ROPT
C
      RETURN
      END 
