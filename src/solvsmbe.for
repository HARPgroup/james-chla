C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SOLVSMBE(SMV1,SMV2,SMA11,SMA22,SMA1,SMA2,SMB11,SMB22)
C
C**********************************************************************C
C
C Solve 2x2 matrix
C
C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 10 april 1999
C
C**********************************************************************C
C
      SMA12 = -SMA2
      SMA21 = -SMA1
      SMDET = SMA11*SMA22 - SMA12*SMA21
      IF (SMDET.EQ.0.0) THEN
        PRINT*, 'Singular matrix: A11, A12, A21, A22, B11, B22'
        PRINT*, SMA11,SMA12,SMA21,SMA22,SMB11,SMB22
        STOP
      END IF
C
      SMDET = 1.0 / SMDET
      SMV1 = (SMB11*SMA22 - SMB22*SMA12) * SMDET
      SMV2 = (SMB22*SMA11 - SMB11*SMA21) * SMDET
C
      RETURN
      END
