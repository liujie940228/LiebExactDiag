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
  
  INTEGER(KIND=IKIND) tt, i, j, IErr
  REAL(KIND=RKIND),ALLOCATABLE :: norm(:), part_nr(:)


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
        PRINT*, "START@ IWidth=", IWidth, " Seed=", ISeed
     ENDIF


     !--------------------------------------------------------------------------
     ! Setting the parameters size passed to function of generating Lieb Matrix
     !--------------------------------------------------------------------------
     
     ucl = (dm * nu) + 1
     n_uc = IWidth**dm
     nt = ucl * n_uc

     !--------------------------------------------------------------------------
     ! Setting the parameters size passed to function DSYEV
     !--------------------------------------------------------------------------

     LDA = nt
     LWMAX = 100000
     
     ! ----------------------------------------------------------
     ! ALLOCATing memory
     ! ----------------------------------------------------------
     
     ALLOCATE ( matr(nt, nt) )
     ALLOCATE ( w( nt ) )
     ALLOCATE ( WORK( LWMAX ) )
     ALLOCATE ( matr_W( nt, nt ) )
     ALLOCATE ( norm(nt) )
     ALLOCATE ( part_nr(nt) )
 

     matr(:,:) = 0.0D0
     matr_W(:,:) = 0.0D0

     CALL MakeLiebMatrixStructrue(dm, nu, IWidth, ucl, n_uc, nt, matr)

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

        DO tt=ISeed, ISeed+NSeed-1

           CALL SRANDOM(tt)

           matr_W(:,:) = matr(:,:)

           ! Give the Lieb matrix different onsite potensial
           DO i=1, n_uc

              matr_W( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = DiagDis*(DRANDOM(tt) - 0.5D0)

              DO j=1, ucl-1

                 matr_W((i-1)*ucl + j + 1 , (i-1)*ucl + j + 1) = RimDiagDis*(DRANDOM(tt) - 0.5D0)

              END DO

           END DO

           LWORK =  -1  !3*nt

           CALL DSYEV( 'V', 'Upper', nt, matr_W, nt, W, WORK, LWORK, INFO )

           LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

           CALL DSYEV( 'V', 'Upper', nt, matr_W, nt, W, WORK, LWORK, INFO )      


           CALL WriteEvals(dm, nu, IWidth, nt, DiagDis, RimDiagDis, W, matr_W, norm, part_nr, tt, INFO)

           matr_W(:,:) = 0d0
           norm(:) = 0d0
           part_nr(:) = 0d0

        END DO ! ISeed cycle

     END DO  ! Disorder cycle

     DEALLOCATE ( matr, w, WORK, matr_W, norm, part_nr )

  END DO ! IWidth cycle


END PROGRAM Lieb

 

  


!!$Function GenerateFileName(dm, nu, n, HudDiagDis, RimDiagDis)
!!$  character(len=100) :: fid,fid2,fid3,fidU,fidD
!!$
!!$  CHARACTER(len=*)   :: GenerateFileName
!!$
!!$  WRITE(fid, '(I5)') dm ; fid = ADJUSTL(fid) 
!!$  WRITE(fid2, '(I5)') nu ; fid2 = ADJUSTL(fid2) 
!!$  WRITE(fid3, '(I5)') n ; fid3 = ADJUSTL(fid3) 
!!$
!!$  WRITE(fidU, '(3ES26.16)') HubDiagDis
!!$  WRITE(fidU, '(F9.1)') HubDiagDis ; fidU = ADJUSTL(fidU)  
!!$
!!$  WRITE(fidD, '(3ES26.16)') RimDiagDis
!!$  WRITE(fidD, '(F9.2)') HubDiagDis ; fidD = ADJUSTL(fidD)  
!!$
!!$  open(unit=1, file='Part_nr_Lieb_'//TRIM(fid)//'_'//TRIM(fid2)//&
!!$       '_dis_'//TRIM(fidU)//'_'//TRIM(fidD)//'_n_'//TRIM(fid3)//'.dat',&
!!$       form='formatted', action='write',status='replace') 
!!$
!!$  RETURN
!!$
!!$END Function GenerateFileName










