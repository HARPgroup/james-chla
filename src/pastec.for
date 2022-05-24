c PasteC.for, jeff ji, 10/20/98
c
c----------------------------------------------------
       subroutine PasteC(nefdcwin,nefdcout)
c purpose: read lines starting with C or c from nefdcwin and paste to nefdcout
c          The current line should be after the second CXX
       parameter (ncread=120)
       character aaa(ncread),cnumber(10)
       nc=0
c
      cnumber(1)='1'
      cnumber(2)='2'
      cnumber(3)='3'
      cnumber(4)='4'
      cnumber(5)='5'
      cnumber(6)='6'
      cnumber(7)='7'
      cnumber(8)='8'
      cnumber(9)='9'
      cnumber(10)='0'
c ensure first reading is a comment line
       read(nefdcwin,99) aaa
       write(nefdcout,99) aaa
       if(aaa(1).ne."C".and.aaa(1).ne."c") then
       write(6,*) "read error in efdcwin.inp or wqwin.inp"
	 write(6,*) "open efdcwin.out or wqwin.out to check"
       stop
       endif
c--
       do nread=1,99999
       read(nefdcwin,99) aaa
 99    format(120a1)
       if(aaa(1).eq.'C'.or.aaa(1).eq.'c') then
       do i=1,10
       if(aaa(2).eq.cnumber(i)) nc=nc+1
       enddo
       write(nefdcout,99) aaa
       if(nc.ge.2) exit
       else
       exit
       endif
       enddo
       if(aaa(1).ne."C".and.aaa(1).ne."c") backspace (nefdcwin)
       return
       end
