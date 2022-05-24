c grid.for, read and modify grid.inp, and make dxdy.inp, lxly.inp cell.inp, jeff ji,  9/10/99
c
c     program grid2
      subroutine grid(IC9,JC9,LDM9)
      INCLUDE 'efdc.par'   ! efdc.cmn is not used for simplicity
c
      parameter (IL=ICM,JL=JCM,IM=IL,JM=JL)
      DIMENSION x(Il,JL),y(IL,JL),dx(IL,JL),dy(IL,JL)
     & ,bdepth(IL,JL),depth(iL,JL),alpha(iL,JL),zrough(IL,JL)
     & ,fsm(iL,JL),mask(IL,JL),cue(IL,JL),cve(IL,JL)
     & ,cun(IL,JL),cvn(IL,JL)
      character a(180)
c -------------------------------------------------------------------
c
c     return         ! disable grid.for
c     write(6,*) "grid.for called"
c
c     open(80,file='grid.inp')      ! from Grid Generator
      open(80,file='grid.inp',status='old',err=97)
c
c     open(21,file='grid2.txt')     !
      open(81,file='dxdy.inp')      ! EFDC grid files
      open(82,file='lxly.inp')
      open(83,file='cell.inp')
      hmin=0.1                      ! minimal water depth
      iveg=1
      idxdy=81
      ilxly=82
      icell=83
c--------------------------------------------------
      undef=-999.9
      do i=1,iL
      do j=1,jL
      x(i,j)=undef
      y(i,j)=undef
      depth(i,j)=undef
      bdepth(i,j)=undef
      fsm(i,j)=0.0
      alpha(i,j)=0.0
      zrough(i,j)=0.0
      enddo
      enddo
c read grid.inp file
      do i=1,22
      read(80,210) (a(j),j=1,180)
210   format(180a1)
c     write(21,210) (a(j),j=1,180)
      enddo
c
      IC9=0
      JC9=0
      LDM9=0
      do m=1,99999
      read(80,*,end=199,err=3009) i,j,x(i,j),y(i,j),dx(i,j),dy(i,j),
     + depth(i,j),bdepth(i,j),alpha(i,j),zrough(i,j)
c
c     write(21,120) i,j,x(i,j),y(i,j),dx(i,j),dy(i,j),depth(i,j),
c    + bdepth(i,j),alpha(i,j),zrough(i,j)
120   format(2i4,2E16.6,999f10.3)
c
      IC9=max(ic9,i)
      JC9=max(Jc9,J)
      LDM9=LDM9+1
      enddo
199   continue
      IC9=IC9+2
      JC9=JC9+2
c
c check array dimensions
      write(6,*) 'wet cell #= ', LDM9, '  (IC,JC) ', ic9,jc9
      if(LDM9.le.0) then
      write(6,*) "Error in grid.inp:  ", LDM9
      stop
      endif
c
      if(LDM9  .gt.(LCM  -2)) then   ! 1
      write(6,*) " Error:  LDM9  >(LCM  -2)", LDM9  , LCM
      stop
      endif
c
      if(IC9   .gt.ICM      ) then   ! 2
      write(6,*) " Error:  IC9   >ICM      ", IC9   , ICM
      stop
      endif
c
      if(JC9   .gt.JCM      ) then   ! 3
      write(6,*) " Error:  JC9   >JCM      ", JC9   , JCM
      stop
      endif
c check array dimensions, jeff Ji, 8/5/99
c make grid files --
c     call makeg (
c    & iL,jL,dx,dy,depth,bdepth,alpha,x,y,hmin,81,82,83,zrough,iveg, ! input arrays & parameters
c    & fsm,mask,cue,cve,cun,cvn)                                     ! working arrays
c-----
c999  continue
c-------------------------------------------------------------
c     stop
c     end
c-----------------------------------------------
c     subroutine makeg(im,jm,dx,dy,depth,bdepth,alpha,x,y,
c    & hmin,idxdy,ilxly,icell,zrough,iveg,                       ! input
c    & fsm,mask,cue,cve,cun,cvn)                                 ! working arrays
c
c
c   Writen by Jeff Ji, on December, 19, 1997. Latest update: 8/25/99
c
c   Attention: modify 39  format!!!
c
c -- Usage --
c  Mkae EFDC grid files,
c
c -- INPUT ---
c
c dx(im,jm),dy(im,jm)              = grid information
c depth(im,jm)                     = initial water depth
c bdepth(im,jm)                    = bottom water depth, it is actually what used in
c                                    calculation
c                                    1) when bottom slope=0 => bdepth=-depth, which means that
c                                    the datum is taken at initial water surface (=0)
c                                    2) when slope.NE.0, bdepth should be elevation referring
c                                    to certain datum, and there is not direct connection
c                                    between depth and bdepth, 1/2/98, talking to John
c alpha(im,jm)                     = angle defined in EFDC manual
c x(im,jm), y(im,jm)               = center grid coordinates
c hmin                             = minimum depth allowed, in meters
c iveg                             = 0, vegetation type
c idxdy,ilxly,icell                = output units for dxdy.inp, lxly.inp, cell.inp
c
c -- OUTPUT ---
c
c in the main program, xdxdy.inp, xlxly.inp, xcell.inp
c
c
c -- WORKING ARRAYS ---
c
c fsm(im,jm),mask(im,jm),cue(im,jm),cve(im,jm),cun(im,jm),cvn(im,jm)
c
c
c     dimension dx(im,jm),dy(im,jm),depth(im,jm),bdepth(im,jm)
c    &, alpha(im,jm),zrough(im,jm)
c    &, fsm(im,jm),mask(im,jm),x(im,jm),y(im,jm)
c    &, cue(im,jm),cve(im,jm),cun(im,jm),cvn(im,jm)
c--------------------------------------------------------------------
c--------------------------------------------------------------------
c     do j=1,jm
c     write(21,*) j,bdepth(5,j)            ! checking
c     enddo
c mask
      do i=1,im
      do j=1,jm
      mask(i,j)=0
      if(dx(i,j).gt.0.0) mask(i,j)=5
      enddo
      enddo
c
      do i=1,im
      do j=1,jm
      fsm(i,j)=mask(i,j)
      enddo
      enddo
c cos & sin
      do i=1,im
      do j=1,jm
      cue(i,j)=0.0
      cve(i,j)=0.0
      cun(i,j)=0.0
      cvn(i,j)=0.0
      if(mask(i,j).ne.0) then
      angle=alpha(i,j)*3.14159/180.0
      cue(i,j)=cos(angle)
      cve(i,j)=-sin(angle)
      cun(i,j)=sin(angle)
      cvn(i,j)=cos(angle)
      endif
      enddo
      enddo
c dxdy.inp
      write(idxdy,1410)
 1410 format('C dxdy.inp file, in free format across columns',/,
     &       'C',/,
     &       'C     I     J        DX            DY            DEPTH
     & BOTTOM ELEV      ZROUGH  VEG TYPE',/,
     &       'C')
C dxdy.inp file, in free format across columns
C
C     I     J        DX            DY            DEPTH     BOTTOM ELEV      ZROUGH  VEG TYPE
C
      do i=1,im
      do j=1,jm
      if(mask(i,j).gt.0.and.mask(i,j).lt.9) then
      if(depth(i,j).lt.hmin) depth(i,j)=hmin
c     Bdepth(i,j)=-depth(i,j)
      write(idxdy,1400) i,j,dx(i,j),dy(i,j),depth(i,j),Bdepth(i,j),
     & zrough(i,j),iveg
 1400 FORMAT(1X,I5,2X,I5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,
     $       2X,E12.5,2X,I3)
      endif
      enddo
      enddo
c lxly.inp
c lxly.inp file, in free format across line
C
C    I     J    XLNUTME       YLTUTMN        CCUE            CCVE          CCUN         CCVN
C
      write(ilxly,1610)
 1610 format('c lxly.inp file, in free format across line',/,
     &       'C',/,
     &       'C    I     J    XLNUTME       YLTUTMN        CCUE
     &       CCVE          CCUN         CCVN',/,
     &       'C')
      do i=1,im
      do j=1,jm
      if(mask(i,j).ne.0) then
      write(ilxly,1600) i,j,x(i,j),y(i,j),cue(i,j)
     &  ,cve(i,j),cun(i,j),cvn(i,j)
 1600 FORMAT(1X,I5,1X,I5,6(1X,E13.6))
c1600 FORMAT(1X,I5,1X,I5,2f13.3,4(1X,E13.6))
      endif
      enddo
      enddo
c cell.inp
      do i=2,im-1
      do j=2,jm-1
      xx=fsm (i+1,j)+fsm (i-1,j)+fsm (i,j+1)+fsm (i,j-1)
     &  +fsm(i-1,j-1)+fsm(i+1,j-1)+fsm(i-1,j+1)+fsm (i+1,j+1)
      if(xx.ge.1.0.and.mask(i,j).eq.0) mask(i,j)=9
      enddo
      enddo
c 4 sides
      do i=1,ic9
       if(mask(i,2).gt.0.and.mask(i,2).lt.9) mask(i,1)=9
       if(mask(i,jc9-1).gt.0.and.mask(i,jc9-1).lt.9) mask(i,jc9)=9
      mask(i,1)=9
      mask(i,jc9)=9
      enddo
      do j=1,jc9
      if(mask(2,j).gt.0.and.mask(2,j).lt.9) mask(1,j)=9
      if(mask(ic9-1,j).gt.0.and.mask(ic9-1,j).lt.9) mask(ic9,j)=9
      mask(1,j)=9
      mask(ic9,j)=9
      enddo
C cell.inp file, i columns and j rows, for Blackstone River, jeff ji, 1/5/98
C             1         2         3         4         5         6         7
C    123456789012345678901234567890123456789012345678901234567890123456789012
C
c
      write(icell,3910)
 3910 format('C cell.inp file, i columns and j rows',/,
     &       'C             1         2         3         4',9x,
     & '5         6         7',/,
     &       'C    12345678901234567890123456789012345678901234567890123
     &4567890123456789012',/,
     &       'C')
      do j=jc9,1,-1
      write(icell,39) j,(mask(i,j),i=1,ic9)  !,j
 39   format(i3,2x,999i1)
c39   format(i3,2x,10i1,1x,i3)
      enddo
c-----------------------------------
      close(80)
      close(81)
      close(82)
      close(83)
c
      GO TO 2000
 3009 WRITE(6,3001)
 3001 FORMAT(1X,'READ ERROR FOR FILE grid.inp ')
      STOP
c---------------------------------------
 97    write(6,*) "grid.inp does not exist and is not used"
       go to 2001
c
 2000  continue
       write(6,*) "grid.inp is read successfully and "
       write(6,*) "  new dxdy.inp, lxly.inp, and cell.inp are produced"
 2001  continue
c     stop
      return
      end
