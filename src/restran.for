C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RESTRAN
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
C **  SUBROUTINE RESTRAN READS A RESIDUAL TRANSPORT FILE 
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      IF(NTSMMT.LT.NTSPTC) THEN
       DO L=2,LA
       READ(99,907)HMP(L),HLPF(L),QSUMELPF(L)
       READ(99,907)(UHLPF(L,K),K=1,KC)
       READ(99,907)(VHLPF(L,K),K=1,KC)
       READ(99,907)(AHULPF(L,K),K=1,KC)
       READ(99,907)(AHVLPF(L,K),K=1,KC)
       READ(99,907)(SALLPF(L,K),K=1,KC)
       READ(99,907)(ABLPF(L,K),K=1,KS)
c      READ(99,907)(ABEFF(L,K),K=1,KS)
       END DO
      ELSE
       DO L=2,LA
       READ(99,907)HMP(L),HLPF(L),QSUMELPF(L)
       READ(99,907)(UHLPF(L,K),K=1,KC)
       READ(99,907)(VHLPF(L,K),K=1,KC)
       READ(99,907)(VPZ(L,K),K=1,KC)
       READ(99,907)(AHULPF(L,K),K=1,KC)
       READ(99,907)(AHVLPF(L,K),K=1,KC)
       READ(99,907)(SALLPF(L,K),K=1,KC)
       READ(99,907)(VPX(L,K),K=1,KS)
       READ(99,907)(VPY(L,K),K=1,KS)
       READ(99,907)(ABLPF(L,K),K=1,KS)
c      READ(99,907)(ABEFF(L,K),K=1,KS)
       END DO
      END IF
c
      DO K=1,KC
      DO L=2,LA
      AHULPF(L,K)=AHULPF(L,K)+AHO
      AHVLPF(L,K)=AHVLPF(L,K)+AHO
      END DO
      END DO
c
      DO K=1,KC
      DO L=2,LA
      AH(L,K)=0.25*(AHULPF(L,K)+AHULPF(L+1   ,K)
     $             +AHVLPF(L,K)+AHVLPF(LNC(L),K))
      END DO
      END DO
C
      IF(NTSMMT.LT.NTSPTC.OR.ISLTMT.EQ.2) THEN
       DO K=1,KC
       DO L=2,LA
       VPZ(L,K)=0.
       END DO
       END DO
       DO K=1,KS
       DO L=2,LA
       VPX(L,K)=0.
       VPY(L,K)=0.
       END DO
       END DO
      END IF
C
  907 FORMAT(12E12.4)
C
C**********************************************************************C
C
      RETURN
      END
