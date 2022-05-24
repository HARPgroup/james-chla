c wcm2000, newer version of wcm93, modified by jeff ji for EFDC, 12/1/2000
c          based on from bblm00_cp1.for
c          Go the end of the file for comments and the original email from Rich Styles
c
      subroutine wcm2000(sg_ub,sg_ab,sg_ur,sg_zr,sg_phicw,     ! Input
     +                      ustarc,   ustarwm,   ustarcw)      ! Output
c
C INPUT: sg_xxxx
C ub      - bed wave orbital velocity (cm/s)
C ab      - bed wave excursion amplitude (cm)
C ur      - mean current at a known height above the bed, zr (cm/s)
C zr      - height of mean current (cm)
C phicw   - angle between the wave and current (radians)
c
C OUTPUT:
C ustarc  - time average bottom shear stress for current (cm/s)
C ustarwm - maximum bottom shear stress for the wave (cm/s)
C ustarcw - maximum combined bottom shear stress (cm/s), the ONE !!
c
c------------------------------------------------------------------------
c     program bblm00_cpl
c
C This is the main driver routine which can be scrapped for
C individual modeling needs.
c
      real sg_ub,sg_ab,sg_ur,sg_zr,sg_phicw,sg_z100,sg_z1p
      real sg_pi,sg_kappa,sg_alpha,sg_dd,sg_ss,sg_nu,sg_g,sg_zntflg
      real sg_tol,sg_ustarcdef,sg_znotdef,sg_znotcdef
c     real sg_shldcr,sg_delta,sg_znot,sg_eta,sg_lambda,sg_znotc
      real sg_shldcr,sg_delta,sg_znot,sg_eta,sg_lambda           ! bug, Ji, 12/1/00
      real sg_znotc,sg_ustarc,sg_ustarcw,sg_u100
      complex sg_mp
      integer sg_n

      common /sg_aa1/ sg_pi,sg_kappa,sg_alpha,sg_dd,sg_ss,sg_nu,
     Z                sg_g,sg_z1p,sg_mp,sg_tol,sg_n
      common /sg_aa2/ sg_ustarcdef,sg_znotdef,sg_znotcdef,sg_zntflg,
     Z                sg_z100

      common /sg_a3/ sg_znot,sg_eta,sg_lambda,sg_znotc,sg_ustarc,
     Z               sg_ustarcw,sg_u100,sg_ustarwm
C-- Input parameters for testing ---------------------
c     sg_ub=40
c     sg_ab=80
c     sg_ur=60
c     sg_zr=100
c     sg_phicw=0
c-------------------------
C set constants & parameters   ! nd-nondimensional
      sg_pi=4.*atan(1.)  ! pi (3.14...) (nd)
      sg_kappa=0.4       ! von Karman's constant (nd)
      sg_alpha=1.        ! non-dimensional closure constant (nd)
      sg_dd=0.04         ! particle diameter (cm/s)
      sg_ss=2.65         ! relative sediment density (nd)
      sg_nu=0.0119       ! viscosity of sea water at 15C (cm^2/s)
      sg_g=981.          ! acceleration due to gravity (cm/s^2)
      sg_z1p=sg_alpha    ! nd-height of the inner layer (nd)
      sg_delta=1./sqrt(2.*sg_z1p) ! (nd)
      sg_mp=cmplx(sg_delta,sg_delta)!  (nd)
      sg_tol=1.0e-4      ! tolerence (nd)
      sg_ustarcdef=2.0   ! default ustarc (cm/s)
      sg_znotdef=0.1     ! default znot (cm)
      sg_shldcr=0.04     ! critical shields parameter for initiation of
                         ! sediment motion (set to a default value) (nd)
      sg_znotcdef=.2     ! default znotc (cm)
      sg_zntflg=0        ! znot flag if = 0 then calculate znot with this
                         ! model, else zntflg = 1 and specify znot (nd)
      sg_z100=100.       ! z=100 cm (1 m)
      sg_n=20            ! max # of itterations for bisection method (nd)


c     print*,'angle',sg_phicw*180./sg_pi
c --- These are the output from above testing parameters, jeff ji, 12/1/00 ----
c D:\Test\Zothers\WCM>wcm
c  angle  0.0000000E+00
c  sg_shldcr  3.7313867E-02
c  zn znotc ustarc ustarcw eta lambda  0.7667559       2.510718       6.513620
c    11.36870       9.978918       82.02074
c
c
c  z1ozn z2ozn   11.86162       20.70296
c
c
c  ur u100   60.00000       60.00024
c --- These are the output from above testing parameters ----

      call sg_bblm99(sg_ub,sg_ab,sg_ur,sg_zr,sg_phicw,
     Z               sg_z1ozn,sg_z2ozn)

C write output for a test run

c     print *,'zn znotc ustarc ustarcw eta lambda',
c    Z   sg_znot,sg_znotc,sg_ustarc,sg_ustarcw,
c    Z   sg_eta,sg_lambda
c     print*,' '
c     print*,' '
c     print*,'z1ozn z2ozn',sg_z1ozn,sg_z2ozn
c     print*,' '
c     print*,' '
c     print*,'ur u100',sg_ur,sg_u100
c
c--
      ustarc=sg_ustarc
      ustarwm=sg_ustarwm
      ustarcw=sg_ustarcw
c     write(6,*) "USTARC,WM,CM",ustarc,ustarwm,ustarcw
c--
c     stop
      return
      end
c--------------
      subroutine sg_bblm99(sg_ub,sg_ab,sg_ur,sg_zr,sg_phicw,
     Z                   sg_z1ozn,sg_z2ozn)

C This subroutine calculates mu=ustarwm/ustarcw, epsilon=ustarc/ustarcw
C and sigma=ub/ustarcw based on an eddy viscosity profile
C for the wave that is linear increasing above znot and constant above
C z1.

C References:  Styles, R. and S. M. Glenn, "Modeling stratified wave and
C              current bottom bondary layers on the continental shelf",
C              JGR (In Press)  (as of Oct-00)
C Styles, R. and S. M. Glenn, "An optimized combined wave and
C            current model for arbitrary bed roughness",
C            Submitted to J. Atmos. and Oceanic Tech.


C nd      -nondimensional
C INPUT:
C omega   - wave radian frequency (s-1)
C ub      - bed wave orbital velocity (cm/s)
C ab      - bed wave excursion amplitude (cm)
C d       - sediment grain diameter (cm)  ! not used in this version
C dd      - representative mean grain diameter to calculate ripple
C           properties (cm)
C s       - relative sediment density = sed. density/fluid density
C           ! not used in this version (nd)
C ss      - representative relative sediment density to calculate
C           ripple properties = sed. density/fluid density (nd)
C ur      - mean current at a known height above the bed, zr (cm/s)
C zr      - height of mean current (cm)
C phicw   - angle between the wave and current (radians)


C OUTPUT:
C sigma   - ratio of bed wave orbital velocity, ub, over combined wave
C           and current shear stress, ustarcw
C mu      - ustarwm/ustarcw (nd)
C epsilon - ustarc/ustarcw (nd)
C ro      - internal friction Rossby number = ustarcw/(omega*znot) (nd)
C eta     - ripple height (cm)
C lambda  - ripple length (cm)
C znot    - hydraulic roughness (cm)
C znotc   - apparent hydraulic roughness (cm)
C u100    - mean current at 100 cm above the bottom (cm/s)
C ustarc  - time average bottom shear stress (cm/s)
C ustarwm - maximum bottom shear stress for the wave (cm/s)
C ustarcw - maximum combined bottom shear stress (cm/s)


      real sg_znotp,sg_z1p,sg_epsilon,sg_tol,sg_ro,sg_phicw
      real sg_alpha,sg_mu,sg_zntflg,sg_u100,sg_ubokur,sg_chi
      real sg_pi,sg_sigma,sg_znotdef,sg_z2,sg_znotc,sg_z1
      real sg_g,sg_nu,sg_eta,sg_lambda,sg_ss,sg_dd,sg_z100
      real sg_ub,sg_ab,sg_ur,sg_zr,sg_kappa,sg_abokb,sg_fwm
      real sg_ustarcdef,sg_ustarc,sg_ustarcw,sg_znotcdef,sg_ustarwm
      real sg_a1,sg_b1,sg_c1,sg_fofa,sg_fofb,sg_fofc,sg_star,sg_shldcr
      complex sg_mp
      integer sg_n,sg_idx

      common /sg_aa1/ sg_pi,sg_kappa,sg_alpha,sg_dd,sg_ss,sg_nu,
     Z                sg_g,sg_z1p,sg_mp,sg_tol,sg_n

      common /sg_aa2/ sg_ustarcdef,sg_znotdef,sg_znotcdef,sg_zntflg,
     Z                sg_z100


      common /sg_a3/ sg_znot,sg_eta,sg_lambda,sg_znotc,sg_ustarc,
     Z               sg_ustarcw,sg_u100,sg_ustarwm

C set constants & parameters

      sg_ubokur=sg_ub/(sg_kappa*sg_ur)

C Chekc #1 - if wave = 0 and ur = 0 then no calculation, set
C ustarc and znot to default values, then exit

      if (sg_ub .le. 0  .and. sg_ur .le. 0) then
          sg_ustarc=sg_ustarcdef
          sg_znot=sg_znotdef
          sg_znotc=sg_znot
          return
      endif

C Check #2 - if wave = 0 then pure current.  Set znot=default, then exit

      if (sg_ub .le. 0) then
        print *, 'ub = 0 ub ur',sg_ub,sg_ur
          sg_znot=sg_znotdef
          sg_znotc=sg_znotdef
          if (sg_zr/sg_znot .le. 1) then
              print *,'warning!!! zr/znot < 1',sg_zr/sg_znot
              sg_ustarc=sg_ustarcdef
              return
          endif
          sg_ustarc=sg_ur*sg_kappa/log(sg_zr/sg_znot)
         return
      endif

C Check znot flag to see if we calculate znot or implement the default.
C If sg_znotflg = 0, then compute znot based on ripple height,
C otherwise (sg_znotflg = 1), take the default znot.

      if (sg_zntflg .eq. 1) then
          sg_znot=sg_znotdef
          sg_chi=4.*sg_nu*sg_ub**2/(sg_dd*((sg_ss-1.)*
     Z           sg_g*sg_dd)**1.5)
          if (sg_chi .le. 2) then
              sg_eta=sg_ab*0.32*sg_chi**(-0.34)
              sg_lambda=sg_ab*2.04*sg_chi**(-0.23)
          else
              sg_eta=sg_ab*0.52*sg_chi**(-1.01)
              sg_lambda=sg_ab*2.7*sg_chi**(-0.78)
          endif
      else
C check initiation of sediment motion criteria, to see if we
C compute znot based on the wave-formed ripples.  If the skin
C friction calculation indicates that sediment is NOT in motion,
C the ripple model is invalid and take the default znot.
C Compute critical Shield's parameter based on grain diameter (dd).

         sg_star=sg_dd/(4.*sg_nu)*sqrt((sg_ss-1.)*sg_g*sg_dd)
         call sg_cshld(sg_star,sg_shldcr)
c        print*,'sg_shldcr',sg_shldcr


C Calculate skin friction shear stress based on Ole Madsen's (1994)
C empirical formula.

          sg_abokb=sg_ab/sg_dd
          if (sg_abokb .le. 100) then
              sg_fwm=exp(7.02*sg_abokb**-0.078-8.82)
          else
              sg_fwm=exp(5.61*sg_abokb**-0.109-7.30)
          endif
          sg_ustarwm=sqrt(sg_fwm/2)*sg_ub
          sg_shdnrm=(sg_ss-1)*sg_dd*sg_g
          sg_shld=sg_ustarwm**2/sg_shdnrm

          if (sg_shld/sg_shldcr .le. 1) then
              sg_znot=sg_znotdef
              sg_eta=0.
              sg_lambda=0.
          else
C Calculate ripple height and length and bottom roughness
              sg_chi=4.*sg_nu*sg_ub**2/(sg_dd*((sg_ss-1.)*
     Z               sg_g*sg_dd)**1.5)
              if (sg_chi .le. 2) then
                  sg_eta=sg_ab*0.30*sg_chi**(-0.39)
                  sg_lambda=sg_ab*1.96*sg_chi**(-0.28)
              else
                  sg_eta=sg_ab*0.45*sg_chi**(-0.99)
                  sg_lambda=sg_ab*2.71*sg_chi**(-0.75)
              endif
              sg_kbs=sg_ab*0.0655*
     Z               (sg_ub**2/((sg_ss-1.)*sg_g*sg_ab))**1.4
              sg_znot=(sg_dd+2.3*sg_eta+sg_kbs)/30
          endif
      endif

C Check #3 - if zr/znot is less than one, set ustarc
C to default value, return.

      sg_zrozn=sg_zr/sg_znot

      if (sg_zrozn .le. 1) then
          sg_ustarc=sg_ustarcdef
          sg_znotc=sg_znotdef
          sg_znot=sg_znotc
             print *,'warning!!! zr/znot < 1',sg_zr/sg_znot
          return
      endif

C Check #4 - if ur = 0 then pure wave.  Calculate wave stress
C based on Ole Madsen's (1994) empirical formula.  Set ustarc
C to a small value to keep model from blowing up and return.

      if (sg_ur .le. 0) then
          sg_abokb=sg_ab/(30*sg_znot)
          sg_row=sg_ab/sg_znot
          if (sg_abokb .le. 100) then
              sg_fwm=exp(7.02*sg_abokb**-0.078-8.82)
          else
              sg_fwm=exp(5.61*sg_abokb**-0.109-7.30)
          endif
          sg_ubouwm=sqrt(2/sg_fwm)
          call sg_purwv(sg_row,sg_ubouwm,sg_znotp,sg_ro)
          sg_ustarwm=sg_ub/sg_ubouwm
          sg_ustarc=sg_ustarcdef
          sg_znotc=sg_znotcdef
          return
      endif

C calculate bottom stresses based on ripple roughness

      sg_row=sg_ab/sg_znot
      sg_a1=1.0e-6
      call sg_bstrssn(sg_row,sg_zrozn,sg_a1,sg_phicw,sg_ubokur,
     Z                sg_mu,sg_epsilon,sg_ro,sg_fofa)
      sg_abokb=sg_ab/(30*sg_znot)
      if (sg_abokb .le. 100) then
          sg_fwm=exp(7.02*sg_abokb**-0.078-8.82)
      else
          sg_fwm=exp(5.61*sg_abokb**-0.109-7.30)
      endif
      sg_ubouwm=sqrt(2/sg_fwm)
      call sg_purwv(sg_row,sg_ubouwm,sg_znotp,sg_ro)
      sg_b1=sg_ubouwm
      sg_fofb=-1*sg_fofa
      sg_c1=0.5*(sg_a1+sg_b1)
      call sg_bstrssn(sg_row,sg_zrozn,sg_c1,sg_phicw,sg_ubokur,
     Z                sg_mu,sg_epsilon,sg_ro,sg_fofc)
      sg_idx=1
      do 10 jj=1,sg_n
         if (sg_idx .eq. 1) then
             sg_sgn=sg_fofb*sg_fofc
             if (sg_sgn .lt. 0) then
                 sg_a1=sg_c1
             else
                 sg_b1=sg_c1
             endif
             sg_c1=0.5*(sg_a1+sg_b1)
             call sg_bstrssn(sg_row,sg_zrozn,sg_c1,sg_phicw,
     Z                       sg_ubokur,sg_mu,sg_epsilon,sg_ro,
     Z                       sg_fofc)
             if (sg_b1-sg_c1 .lt. sg_tol) sg_idx=0
         endif
   10 continue
      sg_sigma=sg_c1

C calculate apparent hydraulic roughness & u100
      sg_lcw=sg_kappa*sg_ab/sg_sigma
      sg_z1=sg_alpha*sg_lcw
      sg_z2=sg_z1/sg_epsilon
      sg_ustarcw=sg_ub/sg_sigma
      sg_ustarc=sg_epsilon*sg_ustarcw
      sg_z1ozn=sg_z1/sg_znot
      sg_z2ozn=sg_z2/sg_znot

      sg_znotc=sg_z2*
     Z         exp(-1*(1-sg_epsilon+sg_epsilon*log(sg_z1ozn)))
      if (sg_z100 .gt. sg_z2) then
          sg_u100=sg_ustarc/sg_kappa*(log(sg_z100/sg_z2)+1-
     Z            sg_epsilon+sg_epsilon*log(sg_z1ozn))
      elseif (sg_z100 .le. sg_z2 .and. sg_zr .gt. sg_z1) then
          sg_u100=sg_ustarc/sg_kappa*sg_epsilon*(sg_z100/sg_z1-1+
     Z           log(sg_z1ozn))
      else
          sg_u100=sg_ustarc/sg_kappa*sg_epsilon*log(sg_z100/sg_znot)
      endif

      return
      end


      subroutine sg_bstrssn(sg_row,sg_zrozn,sg_sigma,sg_phicw,
     Z                      sg_ubokur,sg_mu,sg_epsilon,sg_ro,sg_fofx)

      real sg_x,sg_znotp,sg_z2p,sg_epsilon,sg_mu,sg_sigma,sg_ubokur
      real sg_alpha,sg_kappa,sg_pi,sg_z1p,sg_tol,sg_phi,sg_phicw
      real sg_ber,sg_bei,sg_ker,sg_kei,sg_berp,sg_beip,sg_kerp
      real sg_keip,sg_eps2,sg_fofx
      integer sg_n,sg_idx
      complex sg_bnot,sg_knot,sg_b1,sg_k1,sg_b1p,sg_k1p,sg_mp
      complex sg_bnotp,sg_knotp,sg_ll,sg_nn,sg_gammai,sg_argi

      common /sg_aa1/ sg_pi,sg_kappa,sg_alpha,sg_dd,sg_ss,sg_nu,
     Z                sg_g,sg_z1p,sg_mp,sg_tol,sg_n

      sg_idx=1
      do 20 ii=1,sg_n
         if (sg_idx .eq. 1) then
         sg_ro=sg_row/sg_sigma
         sg_znotp=1./(sg_kappa*sg_ro)
         if (sg_z1p/sg_znotp .gt. 1) then
             do 10 i=1,2
                if (i .eq. 1) sg_x=2.*sqrt(sg_znotp)
                if (i .eq. 2) sg_x=2.*sqrt(sg_z1p)
                if (sg_x .le. 8) then
                    CALL sg_KLVN2(sg_X,sg_BER,sg_BEI,sg_KER,
     Z                            sg_KEI,sg_BERP,sg_BEIP,sg_KERP,
     Z                            sg_KEIP)
                else
                    CALL sg_KEL2(sg_X,sg_KER,sg_KEI,sg_BER,
     Z                   sg_BEI,sg_KERP,sg_KEIP,sg_BERP,sg_BEIP)
                endif
                if (i .eq. 1) then
                    sg_bnot=cmplx(sg_ber,sg_bei)
                    sg_knot=cmplx(sg_ker,sg_kei)
                    sg_bnotp=cmplx(sg_berp,sg_beip)/sqrt(sg_znotp)
                    sg_knotp=cmplx(sg_kerp,sg_keip)/sqrt(sg_znotp)
                endif
                if (i .eq. 2) then
                    sg_b1=cmplx(sg_ber,sg_bei)
                    sg_k1=cmplx(sg_ker,sg_kei)
                    sg_b1p=cmplx(sg_berp,sg_beip)/sqrt(sg_z1p)
                    sg_k1p=cmplx(sg_kerp,sg_keip)/sqrt(sg_z1p)
                endif
   10        continue
             sg_ll=sg_mp*sg_b1+sg_b1p
             sg_nn=sg_mp*sg_k1+sg_k1p
             sg_argi=sg_bnotp*sg_nn/(sg_bnot*sg_nn-sg_knot*sg_ll)+
     Z            sg_knotp*sg_ll/(sg_knot*sg_ll-sg_bnot*sg_nn)
             sg_gammai=-sg_kappa*sg_znotp*sg_argi
             sg_phi=cabs(sg_gammai)
         else
             sg_gammai=-sg_kappa*sg_z1p*sg_mp
             sg_phi=cabs(sg_gammai)
         endif

         if (sg_sigma .gt. 1./sg_phi) then
             sg_sigma=1./sg_phi
         else
             sg_idx=0
         endif
       endif
   20 continue

      sg_mu=sqrt(sg_sigma*sg_phi)
      sg_eps2=-sg_mu**2*abs(cos(sg_phicw))+
     Z        sqrt(1.-sg_mu**4*(sin(sg_phicw))**2)
      sg_epsilon=sqrt(sg_eps2)
      if (sg_epsilon .eq. 0) return

      sg_z2p=sg_z1p/sg_epsilon
      sg_ror=sg_ro/sg_zrozn
      sg_zroz1=1./(sg_alpha*sg_kappa*sg_ror)
      sg_zroz2=sg_epsilon*sg_zroz1
      sg_z1ozn=sg_alpha*sg_kappa*sg_ro
      sg_z2ozn=sg_z1ozn/sg_epsilon

      if (sg_zroz2 .gt. 1 .and. sg_z1ozn .gt. 1) then
          sg_fofx=sg_ubokur*sg_epsilon*(log(sg_zroz2)+1-sg_epsilon+
     Z         sg_epsilon*log(sg_z1ozn))-sg_sigma

      elseif (sg_zroz2 .le. 1 .and. sg_zroz1 .gt. 1
     Z                 .and. sg_z1ozn .gt. 1) then
          sg_fofx=sg_ubokur*sg_epsilon**2*
     Z            (sg_zroz1-1+log(sg_z1ozn))-sg_sigma

      elseif (sg_zroz1 .le. 1 .and. sg_z1ozn .gt. 1) then
          sg_fofx=sg_ubokur*sg_epsilon**2*log(sg_zrozn)-sg_sigma

      elseif (sg_zroz2 .gt. 1 .and. sg_z1ozn .le. 1
     Z                 .and. sg_z2ozn .gt. 1) then
          sg_fofx=sg_ubokur*sg_epsilon*(log(sg_zroz2)+1.
     Z         -1./sg_z2ozn)-sg_sigma

      elseif (sg_zroz2 .le. 1 .and. sg_zroz1 .gt. 1 .and.
     Z    sg_z1ozn .le. 1 .and. sg_z2ozn .gt. 1) then
          sg_fofx=sg_ubokur*sg_epsilon**2*(sg_zroz1-1./sg_z1ozn)-
     Z            sg_sigma

      elseif (sg_zroz2 .gt. 1 .and. sg_z2ozn .le. 1) then
          sg_fofx=sg_ubokur*sg_epsilon*log(sg_zrozn)-sg_sigma
      endif

      return
      end




      subroutine  sg_purwv(sg_row,sg_ubouwm,sg_znotp,sg_ro)

      real sg_x,sg_znotp,sg_ubouwm,sg_ubouwmn
      real sg_alpha,sg_kappa,sg_pi,sg_z1p,sg_tol
      real sg_ber,sg_bei,sg_ker,sg_kei,sg_berp
      real sg_beip,sg_kerp,sg_keip
      integer sg_n
      complex sg_bnot,sg_knot,sg_b1,sg_k1,sg_b1p
      complex sg_k1p,sg_mp,sg_bnotp,sg_knotp
      complex sg_ll,sg_nn,sg_gammai,sg_argi

      common /sg_aa1/ sg_pi,sg_kappa,sg_alpha,sg_dd,sg_ss,sg_nu,
     Z                sg_g,sg_z1p,sg_mp,sg_tol,sg_n

      do 20 ii=1,sg_n
         sg_ro=sg_row/sg_ubouwm
         sg_znotp=1./(sg_kappa*sg_ro)
         if (sg_z1p/sg_znotp .gt. 1) then
             do 10 i=1,2
                if (i .eq. 1) sg_x=2.*sqrt(sg_znotp)
                if (i .eq. 2) sg_x=2.*sqrt(sg_z1p)
                if (sg_x .le. 8) then
                    CALL SG_KLVN2(SG_X,SG_BER,SG_BEI,SG_KER,
     Z                            SG_KEI,SG_BERP,SG_BEIP,SG_KERP,
     Z                            SG_KEIP)
                endif
                if (sg_x .gt. 8) then
                    CALL SG_KEL2(SG_X,SG_KER,SG_KEI,SG_BER,
     Z                   SG_BEI,SG_KERP,SG_KEIP,SG_BERP,SG_BEIP)
                endif
                if (i .eq. 1) then
                    sg_bnot=cmplx(sg_ber,sg_bei)
                    sg_knot=cmplx(sg_ker,sg_kei)
                    sg_bnotp=cmplx(sg_berp,sg_beip)/sqrt(sg_znotp)
                    sg_knotp=cmplx(sg_kerp,sg_keip)/sqrt(sg_znotp)
                endif
                if (i .eq. 2) then
                    sg_b1=cmplx(sg_ber,sg_bei)
                    sg_k1=cmplx(sg_ker,sg_kei)
                    sg_b1p=cmplx(sg_berp,sg_beip)/sqrt(sg_z1p)
                    sg_k1p=cmplx(sg_kerp,sg_keip)/sqrt(sg_z1p)
                endif
   10        continue
             sg_ll=sg_mp*sg_b1+sg_b1p
             sg_nn=sg_mp*sg_k1+sg_k1p
             sg_argi=sg_bnotp*sg_nn/(sg_bnot*sg_nn-sg_knot*sg_ll)+
     Z            sg_knotp*sg_ll/(sg_knot*sg_ll-sg_bnot*sg_nn)
             sg_gammai=-sg_kappa*sg_znotp*sg_argi
             sg_phi=cabs(sg_gammai)
         else
             sg_gammai=-sg_kappa*sg_z1p*sg_mp
             sg_phi=cabs(sg_gammai)
         endif
         sg_ubouwmn=1./sg_phi
         if (abs((sg_ubouwmn-sg_ubouwm)/sg_ubouwmn) .le. sg_tol)then
             sg_ubouwm=sg_ubouwmn
             return
         else
             sg_ubouwm=sg_ubouwmn
         endif
   20 continue

      return
      end


      SUBROUTINE SG_KLVN2(SG_X,SG_BER,SG_BEI,SG_KER,SG_KEI,SG_BERP,
     Z                    SG_BEIP,SG_KERP,SG_KEIP)
C THIS SUB-PROGRAM CALCULATES THE KELVIN FUNCTIONS FOR ARGUMENTS
C LESS THAN 8.
      REAL SG_KER,SG_KEI,SG_KERP,SG_KEIP,SG_BER,SG_BEI,SG_BERP
      REAL SG_BEIP,SG_X,SG_X2,SG_X8
      REAL SG_PI,SG_KAPPA,SG_ALPHA,SG_DD,SG_SS,SG_NU,SG_G
      REAL SG_Z1P,SG_TOL
      INTEGER SG_N
      COMPLEX SG_MP
      common /sg_aa1/ sg_pi,sg_kappa,sg_alpha,sg_dd,sg_ss,sg_nu,
     Z                sg_g,sg_z1p,sg_mp,sg_tol,sg_n

      SG_X2=SG_X/2.
      SG_X8=SG_X/8.

      SG_BER=1.-64.*(SG_X8**4)+113.77777774*(SG_X8**8)-
     Z32.36345652*(SG_X8**12)+2.64191397*(SG_X8**16)-
     Z.08349609*(SG_X8**20)+.00122552*(SG_X8**24)-
     Z.00000901*(SG_X8**28)

      SG_BEI=16.*(SG_X8**2)-113.77777774*(SG_X8**6)+
     Z72.81777742*(SG_X8**10)-10.56765779*(SG_X8**14)+
     Z.52185615*(SG_X8**18)-.01103667*(SG_X8**22)+
     Z.00011346*(SG_X8**26)

      SG_KER=-SG_BER*LOG(SG_X2)+SG_PI/4.*SG_BEI-
     Z.57721566-59.05819744*(SG_X8**4)
     Z+171.36272133*(SG_X8**8)-60.60977451*(SG_X8**12)
     Z+5.65539121*(SG_X8**16)-.19636347*(SG_X8**20)
     Z+.00309699*(SG_X8**24)-.00002458*(SG_X8**28)

      SG_KEI=-SG_BEI*LOG(SG_X2)-SG_PI/4.*SG_BER+6.76454936*(SG_X8**2)
     Z-142.91827687*(SG_X8**6)+124.23569650*(SG_X8**10)
     Z-21.30060904*(SG_X8**14)+1.17509064*(SG_X8**18)
     Z-.02695875*(SG_X8**22)+.00029532*(SG_X8**26)

      SG_BERP=SG_X*(-4.*(SG_X8**2)+14.22222222*(SG_X8**6)-
     Z6.06814810*(SG_X8**10)+.66047849*(SG_X8**14)
     Z-.02609253*(SG_X8**18)+.00045957*(SG_X8**22)
     Z-.00000394*(SG_X8**26))

      SG_BEIP=SG_X*(.5-10.66666666*(SG_X8**4)+11.37777772*(SG_X8**8)
     Z-2.31167514*(SG_X8**12)+.14677204*(SG_X8**16)
     Z-.00379386*(SG_X8**20)+.00004609*(SG_X8**24))

      SG_KERP=-SG_BERP*LOG(SG_X2)-SG_BER/SG_X+SG_PI/4*SG_BEIP
     Z+SG_X*(-3.69113734*(SG_X8**2)+21.42034017*(SG_X8**6)
     Z-11.36433272*(SG_X8**10)+1.41384780*(SG_X8**14)
     Z-.06136358*(SG_X8**18)+.00116137*(SG_X8**22)
     Z-.00001075*(SG_X8**26))

      SG_KEIP=-SG_BEIP*LOG(SG_X2)-SG_BEI/SG_X-SG_PI/4.*SG_BERP
     Z+SG_X*(.21139217-13.39858846*(SG_X8**4)+19.41182758*(SG_X8**8)
     Z-4.65950823*(SG_X8**12)+.33049424*(SG_X8**16)
     Z-.00926707*(SG_X8**20)+.00011997*(SG_X8**24))

      RETURN
      END


      SUBROUTINE SG_KEL2(SG_X,SG_KER,SG_KEI,SG_BER,SG_BEI,SG_KERP,
     Z                   SG_KEIP,SG_BERP,SG_BEIP)

C THIS SUB-PROGRAM CALCULATES THE KELVIN FUNCTIONS FOR ARGUMENTS
C GREATER THAN 8.
      REAL SG_X,SG_X8,SG_X8M,SG_BER,SG_BEI,SG_KER,SG_KEI
      REAL SG_BERP,SG_BEIP,SG_KERP,SG_KEIP
      REAL SG_PI,SG_KAPPA,SG_ALPHA,SG_DD,SG_SS,SG_NU,SG_G
      REAL SG_Z1P,SG_TOL
      INTEGER SG_N
      COMPLEX SG_THETA1,SG_THETAM,SG_PHI1,SG_PHIM,SG_ARGP,SG_ARGM
      COMPLEX SG_FOFX,SG_GOFX,SG_MP
      COMMON /SG_AA1/ SG_PI,SG_KAPPA,SG_ALPHA,SG_DD,SG_SS,SG_NU,
     Z                SG_G,SG_Z1P,SG_MP,SG_TOL,SG_N

      SG_X8=8./SG_X
      SG_X8M=-8./SG_X

      SG_THETA1=CMPLX(0.0,-0.3926991)+
     ZCMPLX(0.0110486,-0.0110485)*SG_X8+
     ZCMPLX(0.0,-0.0009765)*SG_X8**2+
     ZCMPLX(-0.0000906,-0.0000901)*SG_X8**3+
     ZCMPLX(-0.0000252,0.0)*SG_X8**4+
     ZCMPLX(-0.0000034,0.0000051)*SG_X8**5+
     ZCMPLX(0.0000006,0.0000019)*SG_X8**6

      SG_THETAM=CMPLX(0.0,-0.3926991)+
     ZCMPLX(0.0110486,-0.0110485)*SG_X8M+
     ZCMPLX(0.0,-0.0009765)*SG_X8M**2+
     ZCMPLX(-0.0000906,-0.0000901)*SG_X8M**3+
     ZCMPLX(-0.0000252,0.0)*SG_X8M**4+
     ZCMPLX(-0.0000034,0.0000051)*SG_X8M**5+
     ZCMPLX(0.0000006,0.0000019)*SG_X8M**6

      SG_PHI1=CMPLX(0.7071068,0.7071068)+
     ZCMPLX(-0.0625001,-0.0000001)*SG_X8+
     ZCMPLX(-0.0013813,0.0013811)*SG_X8**2+
     ZCMPLX(0.0000005,0.0002452)*SG_X8**3+
     ZCMPLX(0.0000346,0.0000338)*SG_X8**4+
     ZCMPLX(0.0000117,-0.0000024)*SG_X8**5+
     ZCMPLX(0.0000016,-0.0000032)*SG_X8**6

      SG_PHIM=CMPLX(0.7071068,0.7071068)+
     ZCMPLX(-0.0625001,-0.0000001)*SG_X8M+
     ZCMPLX(-0.0013813,0.0013811)*SG_X8M**2+
     ZCMPLX(0.0000005,0.0002452)*SG_X8M**3+
     ZCMPLX(0.0000346,0.0000338)*SG_X8M**4+
     ZCMPLX(0.0000117,-0.0000024)*SG_X8M**5+
     ZCMPLX(0.0000016,-0.0000032)*SG_X8M**6

      SG_ARGM=-1./SQRT(2.)*SG_X*CMPLX(1.,1.)+SG_THETAM
      SG_FOFX=SQRT(SG_PI/(2.*SG_X))*EXP(SG_ARGM)
      SG_KER=REAL(SG_FOFX)
      SG_KEI=IMAG(SG_FOFX)

      SG_ARGP=1./SQRT(2.)*SG_X*CMPLX(1.,1.)+SG_THETA1
      SG_GOFX=1./SQRT(2.*SG_PI*SG_X)*EXP(SG_ARGP)
      SG_BER=REAL(SG_GOFX)-SG_KEI/SG_PI
      SG_BEI=IMAG(SG_GOFX)+SG_KER/SG_PI

      SG_KERP=REAL(-SG_FOFX*SG_PHIM)
      SG_KEIP=IMAG(-SG_FOFX*SG_PHIM)

      SG_BERP=REAL(SG_GOFX*SG_PHI1)-SG_KEIP/SG_PI
      SG_BEIP=IMAG(SG_GOFX*SG_PHI1)+SG_KERP/SG_PI

      RETURN
      END


      SUBROUTINE SG_CSHLD(SG_SSTAR,SG_SHLDC)
C CALCULATE CRITICAL SHIELDS PARAMETER
C FOR INITIATION OF SEDIMENT MOTION
C FROM SHIELDS DIAGRAM.
C SCF IS CORRECTION FACTOR.
      SG_SCF=1.
      IF(SG_SSTAR .GT. 1.5) GO TO 10
         SG_SHLDC=SG_SCF*0.0932*SG_SSTAR**(-.707)
      RETURN
   10 IF(SG_SSTAR .GT. 4.0) GO TO 20
         SG_SHLDC=SG_SCF*0.0848*SG_SSTAR**(-.473)
      RETURN
   20 IF(SG_SSTAR .GT. 10.0) GO TO 30
         SG_SHLDC=SG_SCF*0.0680*SG_SSTAR**(-.314)
      RETURN
   30 IF(SG_SSTAR .GT. 34.0) GO TO 40
         SG_SHLDC=SG_SCF*0.033
      RETURN
   40 IF(SG_SSTAR .GT. 270.) GO TO 50
         SG_SHLDC=SG_SCF*0.0134*SG_SSTAR**.255
      RETURN
   50 SG_SHLDC=SG_SCF*0.056
      RETURN
      END
c
c--Email from Rich Styles -------------------------------------------------
c
c From: Richard Styles [styles@delaware.rutgers.edu]
c Sent: Wednesday, October 11, 2000 9:55 AM
c To: styles@arctic.rutgers.edu; Jeff.Ji@mms.gov
c Cc: bricker@stanford.edu
c Subject: RE: Your Wind Current Model
c
c Jeff and Jeremy,
c
c I have placed the fortran code in our ftp directory:
c
c ftp into arctic.rutgers.edu (anonymous, email for password, the usual stuff)
c goto ftp/pub/styles
c the file is called: bblm00_cpl.f
c
c Below is some info for users that may be helpful.  Anyway, could both
c of you send me a quick reply including the following information:
c the name(s) of all model users, their affiliation, and for what purpose
c they are using the model.  We use this information in our progress reports
c and for writing proposals for further model development.  The more
c people using the model, the more the need for further refinement/development,
c and the more likely to get funding.  You guys understand
c
c thanks,
c
c -rich
c------------
c
c Users,
c
c I have put a version of the bblm code in our ftp directory.
c The bblm consists of a driver routine (bblm00_cpl) which
c calls the subroutine sg_bblm99.  This is where the 'meat and
c potatos' are. A few general comments to help make the transitioning
c of this model a little easier.
c
c 1) if you want to run test cases, the easiest way to check if
c the model is working properly is to set the input current
c height = to 100 cm (default).  The model already predicts what the
c current should be at 100 cm.  If the calculated current does
c not match the input current (less small roundoff error), then something is wrong.
c If this happens, try these things first:
c reduce the convergence tolerence (tol),
c make sure that z1/znot and z2/znot are both greater than
c one.  If z1/znot is less than one, then the equation
c used to compute u100 is wrong and you will never get
c the right answer.
c
c 2) For coupling to shelf circulation models.
c In very shallow water under high roughness values, we found
c that the circulation model (ROMS) sometimes blew up.  This
c occured because the lowest model grid level was on the order
c of znot, the roughness height.  If you are using an 's' type
c vertical grid in shallow water (! 10 m or less) this
c may be a problem.  What we did to correct this, was to set
c a minimum value of
c zr = 100 cm and interpolate the current on a log
c scale between the two grid levels that bracket 100 cm.
c In deeper water, the lowest grid level is usually greater
c than 100 cm, so in that case we used the bottom grid
c level to compute zr.

c 3) the model does not include code to produce concentration
c or current profiles.  Given
c that the model produces the shear stress components, it is
c easy to reconstruct these profiles.  The model also produces
c skin friction components so that the reference concentration
c can be computed.
c
c 4) You are some of the first people other than Hernan to use
c this code, so there may be some bugs.  I think I have
c done a good extermination job, but you never know
c until a lot of people start using the code.  If you
c have problems, let me know.
c
c 5) there are a 3 common blocks.  The common block sg_a3
c contains output variables, and may be eliminated.  The other
