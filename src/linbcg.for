C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE LINBCG (n,b,x,itol,tol,itmax,iter,err,irvec)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
C                LINBCG (LC,FPTMP,P,ITMP,RSQM,ITERM,ITER,RSQ,IRVEC)
C
      INCLUDE 'efdc.par'
C
      PARAMETER (EPS=1.d-14)
C
CU    USES atimes,asolve,snrm
C
      DIMENSION p(LCM),pp(LCM),r(LCM),rr(LCM),z(LCM),zz(LCM),
     $          b(LCM),x(LCM)
C
      iter=0
      call ATIMES(x,r)
      do 11 j=1,n
        r(j)=b(j)-r(j)
        rr(j)=r(j)
11    continue
      IF (IRVEC.GE.50) call ATIMES(r,rr)
      znrm=1.d0
      if(itol.eq.1) then
        bnrm=SNRM(n,b,itol)
      else if (itol.eq.2) then
        call ASOLVE(b,z)
        bnrm=SNRM(n,z,itol)
      else if (itol.eq.3.or.itol.eq.4) then
        call ASOLVE(b,z)
        bnrm=SNRM(n,z,itol)
        call ASOLVE(r,z)
        znrm=SNRM(n,z,itol)
      else
        pause 'illegal itol in linbcg'
      endif
      call ASOLVE(r,z)
100   if (iter.le.itmax) then
        iter=iter+1
        zm1nrm=znrm
        call ASOLVE(rr,zz)
        bknum=0.d0
        do 12 j=1,n
          bknum=bknum+z(j)*rr(j)
12      continue
        if(iter.eq.1) then
          do 13 j=1,n
            p(j)=z(j)
            pp(j)=zz(j)
13        continue
        else
          bk=bknum/bkden
          do 14 j=1,n
            p(j)=bk*p(j)+z(j)
            pp(j)=bk*pp(j)+zz(j)
14        continue
        endif
        bkden=bknum
        call ATIMES(p,z)
        akden=0.d0
        do 15 j=1,n
          akden=akden+z(j)*pp(j)
15      continue
        ak=bknum/akden
        call ATIMES(pp,zz)
        do 16 j=1,n
          x(j)=x(j)+ak*p(j)
          r(j)=r(j)-ak*z(j)
          rr(j)=rr(j)-ak*zz(j)
16      continue
        call ASOLVE(r,z)
        if(itol.eq.1.or.itol.eq.2)then
          znrm=1.d0
          err=SNRM(n,r,itol)/bnrm
        else if(itol.eq.3.or.itol.eq.4)then
          znrm=SNRM(n,z,itol)
          if(abs(zm1nrm-znrm).gt.EPS*znrm) then
            dxnrm=abs(ak)*SNRM(n,p,itol)
            err=znrm/abs(zm1nrm-znrm)*dxnrm
          else
            err=znrm/bnrm
            goto 100
          endif
          xnrm=SNRM(n,x,itol)
          if(err.le.0.5d0*xnrm) then
            err=err/xnrm
          else
            err=znrm/bnrm
            goto 100
          endif
        endif
c       write (*,*) ' iter=',iter,' err=',err
      if(err.gt.tol) goto 100
      endif
      return
      END
