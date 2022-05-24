C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      REAL FUNCTION VALKH(HFFDG)
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
      if(HFFDG.le.0.02) then
        VALKH=HFFDG*HFFDG
        return
      end if
c
      if(HFFDG.ge.10.) then
        VALKH=HFFDG
        return
      end if
c
      DO NTAB=2,1001
      FTMPM1=FUNKH(NTAB-1)
      FTMP  =FUNKH(NTAB  )
      IF (FTMPM1.LE.HFFDG.AND.HFFDG.LT.FTMP) THEN
        VALKH=RKHTAB(NTAB)
     $   -(RKHTAB(NTAB)-RKHTAB(NTAB-1))*(FTMP-HFFDG)/(FTMP-FTMPM1)
        RETURN
      END IF
      END DO
C
      IF (NTAB.EQ.1001) THEN
        WRITE(6,600) RKHTAB(1001)
        WRITE(8,600) RKHTAB(1001)
        STOP
      END IF
c
c **  initialize wave dispersion relation table
c
c      do n=1,501
c       rkh(n)=0.02*float(n-1)
c       frkh(n)=rkh(n)*tanh(rkh(n))
c      end do
c
c      wfrkh=wvfrq*wvfrq*dep(n)/9.8
c      if(wfrkh.le.0.02) wrkh=wfrkh*wfrkh
c      if(wfrkh.ge.10.) wrkh=wfrkh
c      if(wfrkh.gt.0.02.and.wfrkh.lt.10.) then
c        do m=1,500
c         if(wfrkh.gt.frkh(m).and.wfrkh.lt.frkh(m+1)) then
c           drdf=(rkh(m+1)-rkh(m))/(frkh(m+1)-frkh(m))
c           wrkh=drdf*(wfrkh-frkh(m))+rkh(m)
c           go to 200 
c         end if
c        end do
c      end if
c  200 continue
c
C
C
  600 FORMAT(' WAVE DISPERSION TABLE OUT OF BOUNDS KH = ',E12.4)
C 
      RETURN
      END
