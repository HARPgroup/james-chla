C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      FUNCTION SNRM(n,sx,itol)
      INCLUDE 'efdc.par'
      DIMENSION sx(LCM)
      if (itol.le.3)then
        SNRM=0.
        do 11 i=1,n
          SNRM=SNRM+sx(i)**2
11      continue
        SNRM=sqrt(SNRM)
      else
        isamax=1
        do 12 i=1,n
          if(abs(sx(i)).gt.abs(sx(isamax))) isamax=i
12      continue
        SNRM=abs(sx(isamax))
      endif
      return
      END
