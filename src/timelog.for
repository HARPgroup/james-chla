C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
c
      subroutine timelog(n)
c
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
C
c     character*8 mrmdate,mrmtime*11
	character*8 mrmdate*9,mrmtime*11  ! Ji, 7/23/99
c Write out Model Time Step and SUN/PC system clock time to time.log file:
c     call date(mrmdate)
      call timef(mrmtime) 
      write(9,100) n, mrmdate, mrmtime
  100 format(' ','N=',I10,5x,'Date = ',A8,5x,'Time = ',A11)
      call flush(9)
      return
      end
 
