!
!
!   Call function DSYEV() to calculate eigenvalues and eigenvectors
!   for Lieb matrix and its extendsions(2D and 3D)
!
!
!--------------------------------------------------------------------------------------

PROGRAM Lieb

!!$  use, intrinsic :: iso_c_binding
  USE MyNumbers  
  USE CConstants
  USE IConstants
  USE IPara
  USE DPara
  USE IChannels
  USE RNG
  
  IMPLICIT NONE

  !----------------------------------------
  ! Variable declaration
  !----------------------------------------

  ! paramters for Lieb matrix
  INTEGER(KIND=IKIND) IWidth
  
  INTEGER(KIND=IKIND) &
       ucl, &  ! the number of atoms in a unit cell
       n_uc,& ! the number of unit cell
       nt  ! the whole number of atoms in system

  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE:: matr

  ! Parameters for call function DSYEV()
  
  INTEGER(KIND=IKIND) LDA, LWMAX, INFO, LWORK

  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE:: W, WORK

  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE:: matr_W

  INTRINSIC        INT, MIN
  EXTERNAL         DSYEV

  ! Parameters for eigenverctor, participation numbers
  
  INTEGER(KIND=IKIND) i, j, IErr, sumIErr
  REAL(KIND=RKIND),ALLOCATABLE :: norm(:), part_nr(:)

  Character*38 Matrixname


  ! ----------------------------------------------------------
  ! start of main code
  ! ----------------------------------------------------------
  
  ! ----------------------------------------------------------
  ! protocol feature
  ! ----------------------------------------------------------
  
  PRINT*,"LiebExactDia ", RStr, DStr, AStr 

  ! ----------------------------------------------------------
  ! inout handling
  ! ----------------------------------------------------------
  
  CALL  Input(IErr)
  IF(IErr.NE.0) THEN
     PRINT*,"main: Input() finds IErr=", IErr
     STOP
  ENDIF

  ! ----------------------------------------------------------
  ! start of main IWidth loop
  ! ----------------------------------------------------------
  
  DO IWidth= Width0, Width1, dWidth

     ! ----------------------------------------------------------
     IF(IWriteFlag.GE.0) THEN
        PRINT*, "START@ IWidth=", IWidth
     ENDIF


     !--------------------------------------------------------------------------
     ! Setting the parameters size passed to function of generating Lieb Matrix
     !--------------------------------------------------------------------------
     
     ucl = (Dim * Nx) + 1
     n_uc = IWidth**Dim
     nt = ucl * n_uc

     !--------------------------------------------------------------------------
     ! Setting the parameters size passed to function DSYEV
     !--------------------------------------------------------------------------

     LDA = nt
     LWMAX = 100000
     
     ! ----------------------------------------------------------
     ! ALLOCATing memory
     ! ----------------------------------------------------------

     sumIErr= 0
     ALLOCATE ( matr(nt, nt), STAT=IErr ); sumIErr=sumIErr+IErr
     ALLOCATE ( w( nt ), STAT=IErr ); sumIErr=sumIErr+IErr
     ALLOCATE ( WORK( LWMAX ), STAT=IErr ); sumIErr=sumIErr+IErr
     ALLOCATE ( matr_W( nt, nt ), STAT=IErr ); sumIErr=sumIErr+IErr
     ALLOCATE ( norm(nt), STAT=IErr ); sumIErr=sumIErr+IErr
     ALLOCATE ( part_nr(nt), STAT=IErr ); sumIErr=sumIErr+IErr

     IF(IErr.NE.0) THEN
        PRINT*,"main: Input() finds sumIErr=", sumIErr
        STOP
     ENDIF
 

     matr(:,:) = 0.0D0
     matr_W(:,:) = 0.0D0

     CALL MakeLiebMatrixStructrue(Dim, Nx, IWidth, ucl, n_uc, nt, matr)

!!$
!!$     END DO
!!$     DO i= 1, nt
!!$        write(*,108) i
!!$        Do j= 1, nt 
!!$           IF( matr(i,j).ne.(0.0) )Then
!!$              Write(*,108) j
!!$           END IF
!!$           
!!$        END DO
!!$        write(*,*)" "
!!$     END DO
!!$108  format(1x,1I3\)
     
     norm(:) = 0d0
     part_nr(:) = 0d0

     DO DiagDis= DiagDis0,DiagDis1,dDiagDis
        
        DO Seed= ISeed, ISeed + NSeed -1
           
           CALL SRANDOM(Seed)

           matr_W(:,:) = matr(:,:)

           ! Give the Lieb matrix different onsite potensial
           DO i=1, n_uc

              matr_W( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = DiagDis*(DRANDOM(Seed) - 0.5D0)

              DO j=1, ucl-1

                 matr_W((i-1)*ucl + j + 1 , (i-1)*ucl + j + 1) = RimDiagDis*(DRANDOM(Seed) - 0.5D0)

              END DO

           END DO

           !-------------------------------------------------------------------------
           ! Generate the file containing full matrix to comparing with sparse matrix
           !-------------------------------------------------------------------------
           
!!$           WRITE(Matrixname, '(A9,I1,I1,A2,I4.4,A6,I4,A4)') &
!!$                "FullMat-L", Dim, Nx, "-M", IWidth,  &
!!$                "-Seed-", Seed, &
!!$                ".txt"
!!$
!!$
!!$           OPEN(UNIT= 10, FILE=Matrixname)   
!!$
!!$           DO i=1, nt
!!$              DO j=1, nt
!!$                 WRITE(10,102) matr_W(i,j)              
!!$              END DO
!!$              WRITE(10,*) " "
!!$           END DO
!!$
!!$102        FORMAT(1x,F15.6\)
!!$
!!$           CLOSE(10)

           !----------------------------------
           ! calling the lapack function
           !----------------------------------

           LWORK =  -1  !3*nt

           CALL DSYEV( 'V', 'Upper', nt, matr_W, nt, W, WORK, LWORK, INFO )

           LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

           CALL DSYEV( 'V', 'Upper', nt, matr_W, nt, W, WORK, LWORK, INFO )      


           CALL WriteEvals(Dim, Nx, IWidth, nt, DiagDis, RimDiagDis, W, matr_W, norm, part_nr, Seed, INFO)

           matr_W(:,:) = 0d0
           norm(:) = 0d0
           part_nr(:) = 0d0

        END DO ! Seed cycle

     END DO  ! Disorder cycle

     DEALLOCATE ( matr, w, WORK, matr_W, norm, part_nr )

  END DO ! IWidth cycle


END PROGRAM Lieb

