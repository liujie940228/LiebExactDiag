! ********************************************************************
!       
! For Lieb model, call the function DSYEV() to calculate eigenvalues
! and eigenvectors
!
! ********************************************************************
      

! --------------------------------------------------------------------
! Input:
!
! IErr	error code
!----------------------------------------------------------------------

SUBROUTINE Input(IErr)

  USE MyNumbers
  
  USE CConstants
  USE IConstants
  
  USE IPara
  USE DPara
  
  USE IChannels

  
  INTEGER IErr, ILine
!  REAL(KIND=RKIND) RIter
  
  !	PRINT*,"DBG: Input()"
  
  IErr = 0
  ILine= 0
  
  OPEN(UNIT= IChInp, ERR=120, FILE= "Exactdiag.inp",STATUS= 'OLD')

  ILine= ILine+1
  READ(IChInp,10,ERR=20) ISeed
  !PRINT*,"ISeed        = ",ISeed

  ILine= ILine+1
  READ(IChInp,10,ERR=20) NSeed
  !PRINT*,"NSeed        = ",NSeed

  ILine= ILine+1
  READ(IChInp,10,ERR=20) dm
  !PRINT*,"Dim         = ",Dim

  ILine= ILine+1
  READ(IChInp,10,ERR=20) nu
  !PRINT*,"Nx          = ",Nx
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IBCFlag
  !PRINT*,"IBCFlag      = ", IBCFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IRNGFlag
  !PRINT*,"IRNGFlag     = ", IRNGFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IKeepFlag
  !PRINT*,"IKeepFlag    = ", IKeepFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) IWriteFlag
  !PRINT*,"IWriteFlag   = ", IWriteFlag
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) Width0
  !PRINT*,"Width0       = ",Width0
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) Width1
  !PRINT*,"Width1       = ", Width1
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20) dWidth
  !PRINT*,"dWidth       = ", dWidth
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) DiagDis0
  !PRINT*,"DiagDis0     = ", DiagDis0
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) DiagDis1
  !PRINT*,"DiagDis1     = ", DiagDis1
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20) dDiagDis
  !PRINT*,"dDiagDis     = ", dDiagDis

  ILine= ILine+1
  READ(IChInp,15,ERR=20) RimDiagDis
  !PRINT*,"dDiagDis     = ", dDiagDis

  
10 FORMAT(16X,I15.1)
  ! 10	FORMAT("IMAXIteration= ",I15.1)
15 FORMAT(16X,F18.9)
  ! 15     FORMAT("IMAXIteration= ",F18.9)

  ! check the parameters for validity
  
  IF(IWriteFlag.GE.2) THEN
     PRINT*,"ISeed        = ", ISeed
     PRINT*,"NSeed        = ", NSeed
     PRINT*,"dm           = ", dm
     PRINT*,"nu           = ", nu
     PRINT*,"IBCFlag      = ", IBCFlag
     PRINT*,"IRNGFlag     = ", IRNGFlag
     PRINT*,"IKeepFlag    = ", IKeepFlag
     PRINT*,"IWriteFlag   = ", IWriteFlag
     PRINT*,"Width0       = ", Width0
     PRINT*,"Width1       = ", Width1
     PRINT*,"dWidth       = ", dWidth
     PRINT*,"DiagDis0     = ", DiagDis0
     PRINT*,"DiagDis1     = ", DiagDis1
     PRINT*,"dDiagDis     = ", dDiagDis
     Print*,"RimDiagDis   = ", RimDiagDis
    
  ENDIF

  CLOSE(IChInp)
  RETURN

120 PRINT*,"Input(): ERR in OPEN()"  
  IErr= 1
  RETURN

20 PRINT*,"Input(): ERR in WRITE()"  
  IErr= 1
  RETURN
  
END SUBROUTINE Input

! --------------------------------------------------------------------
! Input:
!
! IErr	error code
!----------------------------------------------------------------------

FUNCTION GetFileName(fnamestr,vdata,n)
!write(*,*) trim(GetFileName('/Output/QH_L(I4)_NL(I1)_NS(I1)_B(F5.3)_C(F4.2)_S(I5).dat',&
!     (/L_INPUT,NLevels,NSpins,BField,Coul,Seed/),6))


  IMPLICIT NONE

  INTEGER,INTENT(in)         :: n
  REAL(8),INTENT(in)         :: vdata(n)
  CHARACTER(len=*),INTENT(in) :: fnamestr

  CHARACTER(len=*)   :: GetFileName
  CHARACTER(len=500) :: str,fstr,vstr
  
  INTEGER           :: i

  str = fnamestr
  
  DO i=1,n
     
     fstr = str(INDEX(str,'('):INDEX(str,')'))
     IF(fstr(2:2).EQ.'I')THEN
        WRITE(vstr,fstr) INT(vdata(i))
     ELSE
        WRITE(vstr,fstr) vdata(i)
     END IF
  
     str = str(1:INDEX(str,'(')-1) // TRIM(ADJUSTL(vstr)) // str(INDEX(str,')')+1:LEN(str)) 
     
  END DO
  
  GetFileName = TRIM(ADJUSTL(str))

  RETURN
END FUNCTION GetFileName


!--------------------------------------------------------------------
! CheckOutput:
!
! IErr	error code

SUBROUTINE CheckOutput( IWidth, Energy, DiagDis, ISeed, IErr )

  USE MyNumbers 
  USE IChannels
!  USE DPara
!  USE IPara

  
  INTEGER(KIND=IKIND) IWidth, IErr, ERR, Iseed
  REAL(KIND=RKIND) DiagDis, Energy
  
  CHARACTER*28 FileName
  
  !PRINT*,"DBG: CheckOutput()"
  
  IErr= 0
  
  !   WRITE out the input parameter
  IF( Energy .GE. 0.0D0) THEN
     WRITE(FileName, '(A4,I4.4,A1,I4.4,A1,I4.4,A1,I4.4,A4)') &
          "Eval",IWidth, "-", &
          NINT(100.0D0*ABS(Energy)),"-", &
          NINT(100.0D0*ABS(DiagDis)), "-", &
          !Number, "-", &
          ISeed, ".raw"
  ELSE
     WRITE(FileName, '(A4,I4.4,A2,I4.4,A1,I4.4,A1,I4.4,A4)') &
          "Eval",IWidth, "-m", &
          NINT(100.0D0*ABS(Energy)),"-", &
          NINT(100.0D0*ABS(DiagDis)), "-", &
          !Number, "-", &
          ISeed, ".raw"
  ENDIF
  
  OPEN(UNIT= IChOut, ERR= 10, STATUS= 'NEW', FILE=FileName)
  
  IErr= 0
  
20 CLOSE(UNIT= IChOut, ERR= 100)
  
  RETURN
  
10 WRITE(*,15) FileName
15 FORMAT(" CheckOutput(): ", A28,&
        " exists -- skipped!")
  
  IErr= 2
  GOTO 20
  
  !  ERR in CLOSE detected
100 &
  PRINT*,"CheckOutput(): ERR in CLOSE()"
  IErr= 1
  RETURN
  
END SUBROUTINE CheckOutput


!--------------------------------------------------------------------
! WriteOutputEVal:
!
! IErr	error code

SUBROUTINE WriteOutputEVal(NEVals, EIGS, IWidth, Energy, DiagDis, ISeed, IErr)


  USE MyNumbers
  USE IChannels
!  USE DPara
!  USE IPara

  
  INTEGER(KIND=IKIND) IWidth, IErr, ERR, NEVals, Iseed, i
  REAL(KIND=RKIND) DiagDis, Energy
  REAL(KIND=RKIND) EIGS(NEVals)
  
  CHARACTER*28 FileName
  
  IErr= 0
  
  !   WRITE out the input parameter
  IF(Energy.GE.0.0D0) THEN
     WRITE(FileName, '(A4,I4.4,A1,I4.4,A1,I4.4,A1,I4.4,A4)') &
          "Eval",IWidth, "-", &
          NINT(100.0D0*ABS(Energy)),"-", &
          NINT(100.0D0*ABS(DiagDis)), "-", &
          !Number, "-", &
          ISeed, ".raw"
  ELSE
     WRITE(FileName, '(A4,I4.4,A2,I4.4,A1,I4.4,A1,I4.4,A4)') &
          "Eval",IWidth, "-m", &
          NINT(100.0D0*ABS(Energy)),"-", &
          NINT(100.0D0*ABS(DiagDis)), "-", &
          !Number, "-", &
          ISeed, ".raw"
  ENDIF
  
  !        IF(IWriteFlag .GE. 2) THEN
  PRINT*, "eVAL file    ", FileName
  !        ENDIF
  
  OPEN(UNIT= IChEVal, ERR= 10, STATUS='UNKNOWN', FILE=FileName)
!!$  DO i=1,NEVals
!!$     WRITE(IChEVal, FMT=15, ERR=20) EIGS(i)
!!$  ENDDO
!!$  CLOSE(UNIT=IChEVal, ERR= 30)
  
  RETURN
  
15 FORMAT(f30.20)

  !	error in OPEN detected
10 PRINT*,"WriteOutputEVals(): ERR in OPEN()"
  IErr= 1
  RETURN
  
  !	error in WRITE detected
20 PRINT*,"WriteOutputEVals(): ERR in WRITE()"
  IErr= 1
  RETURN

  ! ERR in CLOSE detected
30 PRINT*,"OutputEVals(): ERR in CLOSE()"
  IErr= 1
  RETURN
  
END SUBROUTINE WriteOutputEVal


!--------------------------------------------------------------------
! WriteOutputEVec:
!
! IErr	error code

SUBROUTINE WriteOutputEVec( Inum, NEVals, Lsize, VECS, VECS_size, &
                   IWidth, Energy, DiagDis, ISeed, IErr)

  USE MyNumbers
  USE IChannels
!  USE DPara
!  USE IPara

  INTEGER(KIND=IKIND) Inum, Iseed, IWidth, IErr, ERR, Lsize, VECS_size, NEVals, i
  REAL(KIND=RKIND) DiagDis, Energy

  REAL(KIND=RKIND) VECS(VECS_size)

  CHARACTER*40 FileName

  IErr= 0
  
  !   WRITE out the input parameter
  IF(Energy.GE.0.0D0) THEN
     WRITE(FileName, '(A4,I4.4,A1,I4.4,A1,I4.4,A1,I4.4,A1,I4.4,A4)') &
          "Evec",IWidth, "-", &
          NINT(100.0D0*ABS(Energy)),"-", &
          NINT(100.0D0*ABS(DiagDis)), "-", &
          Inum, "-", &
          ISeed, ".raw"
  ELSE
     WRITE(FileName, '(A4,I4.4,A1,I4.4,A1,I4.4,A1,I4.4,A1,I4.4,A4)') &
          "Evec",IWidth, "-m", &
          NINT(100.0D0*ABS(Energy)),"-", &
          NINT(100.0D0*ABS(DiagDis)), "-", &
          Inum, "-", &
          ISeed, ".raw"
  ENDIF

  print*,'evector file ',FileName

  OPEN(UNIT= IChEVec, ERR= 40, STATUS= 'UNKNOWN', FILE=FileName)

  DO i= 1+( Lsize*( Inum -1) ), Lsize*Inum
     WRITE(UNIT=IChEVec, FMT=45, ERR=50) VECS(i)
  ENDDO

  CLOSE(UNIT= IChEVec, ERR= 60)

  RETURN

45 FORMAT(f30.20)

  !	error in OPEN detected
40 PRINT*,"WriteOutputEVec(): ERR in OPEN()"
  IErr= 1
  RETURN

  !	error in WRITE detected
50 PRINT*,"WriteOutputEVec(): ERR in WRITE()"
  IErr= 1
  RETURN

  ! ERR in CLOSE detected
60 PRINT*,"OutputEVec(): ERR in CLOSE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteOutputEVec


