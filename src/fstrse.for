C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
c
      FUNCTION FSTRSE(void)
c
      voidlog=log(void)
c
      fstrselog=-2.03249*(voidlog**6)+13.2918*(voidlog**5)
     $          -29.9227*(voidlog**4)+22.6251*(voidlog**3)
     $          +6.84764*(voidlog**2)-19.3663*voidlog
     $        + 19.8503
      FSTRSE=exp(fstrselog)
c
      return
      end
