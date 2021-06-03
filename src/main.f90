!
!
!   Call function DSYEV() to calculate eigenvalues and eigenvectors
!   for Lieb matrix and its extendsions(2D and 3D)
!
!
!--------------------------------------------------------------------------------------

PROGRAM Lieb

  use, intrinsic :: iso_c_binding
  
  IMPLICIT NONE

  INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(15,307)

  !----------------------------------------
  ! Variable declaration
  !----------------------------------------

  ! paramters for Lieb matrix
  
  INTEGER(KIND=IKIND) &
       dm, & ! the dimension
       n, &  ! the number of unit cell in each dimension
       nu, & ! the number of site in each edge
       ucl, &  ! the number of atoms in a unit cell
       n_uc,& ! the number of unit cell
       nt  ! the whole number of atoms in system

  PARAMETER ( dm = 3, n = 5, nu = 1)
  PARAMETER ( ucl = (dm * nu) + 1 )
  PARAMETER ( nt = ucl * n**dm )
  PARAMETER ( n_uc = n**dm )

  REAL(KIND=RKIND) HubDiagDis, RimDiagDis, x
  REAL(KIND=RKIND) matr( nt, nt )

  PARAMETER ( HubDiagDis = 2.0D0, RimDiagDis = 0.0D0 )

  
  ! Parameters for call function DSYEV()
  
  INTEGER(KIND=IKIND) LDA, LWMAX, INFO, LWORK
  
  PARAMETER        ( LDA = nt )
  PARAMETER        ( LWMAX = 100000 )
  REAL(KIND=RKIND) W( nt ), WORK( LWMAX ), matr_W( nt, nt )
  
  INTRINSIC        INT, MIN
  EXTERNAL         DSYEV

  
  ! Parameters for eigenverctor, participation numbers
  
  INTEGER(KIND=IKIND) Seed, ISeed, & ! the number of disorder realization
       tt, i, j
  PARAMETER( ISeed = 10 )
  REAL(KIND=RKIND),ALLOCATABLE :: norm(:), part_nr(:)
 
 
  ! Parameters for clock
  
  INTEGER(KIND=IKIND) t1,t2,clock_rate,clock_max
  REAL(KIND=RKIND) :: time_CPU_before
  
  call init_random_seed()
  call srand(seed) 

  ! Start timing
  CALL system_clock (t1, clock_rate, clock_max)  

  matr(:,:) = 0.0D0
  matr_W(:,:) = 0.0D0
  
  CALL MakeLiebMatrixStructrue(dm, nu, n, ucl, n_uc, nt, matr)

! EASY to check the Lieb matrix structure in screen
!!$ DO i= 1, nt
!!$    Print*, "---------", i, " line         -----------"
  !!$
!!$    Do j= 1, nt 
!!$       Print*, j, " ", matr(i,j)    
!!$    END DO
!!$
!!$    Print*, "-------------------------"
!!$    
!!$ END DO
 
  ALLOCATE ( norm(nt) )
  ALLOCATE ( part_nr(nt) )

  norm(:) = 0d0
  part_nr(:) = 0d0

  
  DO tt=1,ISeed

     matr_W(:,:) = matr(:,:) 

     DO i=1, n_uc
        CALL Random_number(x)        
        matr_W( (i-1)*ucl + 1 , (i-1)*ucl + 1 ) = HubDiagDis*(x - 0.5D0)
        
        DO j=1, ucl-1
           CALL Random_number(x) 
           matr_W((i-1)*ucl + j + 1 , (i-1)*ucl + j + 1) = RimDiagDis*(x - 0.5D0)
           
        END DO
     END DO
           
     LWORK =  -1  !3*nt

     CALL DSYEV( 'V', 'Upper', nt, matr_W, nt, W, WORK, LWORK, INFO )

     LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

     CALL DSYEV( 'V', 'Upper', nt, matr_W, nt, W, WORK, LWORK, INFO )      
     
     
     CALL WriteEvals(dm, nu, n, nt, HubDiagDis, RimDiagDis, W, matr_W, norm, part_nr, tt, INFO)
    
     matr_W(:,:) = 0d0
     norm(:) = 0d0
     part_nr(:) = 0d0
     
  END DO


  CALL cpu_time(time_CPU_before)



END PROGRAM Lieb

 
SUBROUTINE WriteEvals(dm, nu, n, nt, HubDiagDis, RimDiagDis, W, matr_W, norm, part_nr, ISample, INFO)

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(15,307)

  INTEGER(KIND=IKIND) dm, nu, n, nt, ISample, INFO
  INTEGER(KIND=IKIND) i, j
  REAL(KIND=RKIND) HubDiagDis, RimDiagDis
  REAL(KIND=RKIND) W( nt ), matr_W( nt, nt ), norm( nt ), part_nr( nt )
  
  CHARACTER*100 FileName

  WRITE(FileName, '(A5,A1,I1,I1,A2,I4.4,A1,A2,I4.4,A1,A2,I4.4,A1,I4.4,A4)')&
       "Eval-","L",dm,nu,&
       "-M",n, "-",&
       "WH", NINT(100.D0*ABS(HubDiagDis)),&
       "-","WR", NINT(100.D0*ABS(RimDiagDis)),&
       "-",ISample,".raw"

  Print*, "FileName: ", FileName

  OPEN(Unit=9, FILE=FileName)

  IF(INFO==0)THEN

     DO i=1,nt
        DO j=1,nt	
           norm(i) = norm(i) + matr_W(i,j)**2
           part_nr(i) = part_nr(i) + matr_W(i,j)**4
        END DO
     END DO

     DO i=1,nt
        write(9,'(2f30.20)') W(i) , (norm(i)**2) / part_nr(i)
     END DO

  ELSE
     PRINT*, "ERROR IN CALL DSYEV()" 
  END IF

  CLOSE(9)

  RETURN
       
END SUBROUTINE WRITEEVALS
  


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

 
Subroutine MakeLiebMatrixStructrue(dm, nu, n, ucl, n_uc, nt, matr)

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: IKIND = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(15,307)

  INTEGER(KIND=IKIND) &
       dm, & ! the dimension
       n, &  ! the number of unit cell in each dimension
       nu, & ! the number of site in each edge
       ucl, & ! the number of atoms in a unit cell
       nt, &  ! the whole number of atoms in system
       n_uc   ! the number of unit cell

  INTEGER(KIND=IKIND) i, j, k, ind

  LOGICAL(KIND=8) Flag

  INTEGER(KIND=IKIND), ALLOCATABLE :: ucl_d(:) 
  REAL(KIND=RKIND) matr(nt, nt), matr_W( nt, nt )
  
  matr(:,:) = 0.0D0

  IF(dm==2)THEN
     ALLOCATE(ucl_d(dm))
     ucl_d(1) = 2          ! The first Rim atoms
     ucl_d(2) = nu + 2
  ELSE IF(dm==3)THEN
     ALLOCATE(ucl_d(dm))
     ucl_d(1) = 2
     ucl_d(2) = nu + 2
     ucl_d(3) = 2 * nu + 2
  ELSE
     Print*, "We Only Finished the 2D and 3D cases for Lieb model"
     STOP
  END IF
    
!!!!!!! Start construct Lieb Matrix
  
  DO i=1,n_uc
     ! ----------------------- hub atoms ------------------------!
     ind = (i-1)*ucl + 1

     matr(ind,ind) = 0.0D0
     ! hopping terms for hub atoms 
     DO j=1,dm

        IF(j==1)THEN
           Flag = ( MOD(i,n)==1 )
        ELSE IF(j==2)THEN
           Flag = ( MOD(i,n**2)>0 .AND. MOD(i,n**2)<=n )
        ELSE
           Flag = ( i<=n**2 )
        END IF

        matr(ind, (i-1)*ucl + ucl_d(j)) = 1.0D0
        matr((i-1)*ucl + ucl_d(j), ind) = 1.0D0

        IF(Flag)THEN
           matr(ind, (i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1) + ucl*(n)**j) = 1.0D0
           matr((i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1) + ucl*(n)**j, ind) = 1.0D0
        ELSE
           matr(ind, (i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1)) = 1.0D0
           matr((i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1), ind) = 1.0D0
        END IF

     END DO
     ! ----------------------- Rim atoms ------------------------!

     ! ------- For rim atoms except close to hub atom of orther unit cell-----------

     IF(nu>1)THEN
        DO k=1, nu-1

           DO j=1, dm

              ind = (i-1)*ucl + ucl_d(j) + k - 1

              matr(ind,ind) = 0.0D0
              !Hopping term
              IF(k==1)THEN
                 matr(ind, ind - (j-1)*nu -1) = 1.0D0
                 matr(ind - (j-1)*nu -1, ind) = 1.0D0                
              ELSE
                 matr(ind, ind - 1) = 1.0D0
                 matr(ind - 1, ind) = 1.0D0
              END IF
              matr(ind, ind + 1) = 1.0D0
              matr(ind + 1, ind) = 1.0D0

           END DO

        END DO
     END IF


     ! For rim atoms close to hub atoms of other unit cells

     DO j=1,dm
        IF(j==1)THEN
           Flag = ( MOD(i,n**j)==0 )
        ELSE IF(j==2)THEN
           Flag = ( MOD(i,n**j)==0 .OR. (MOD(i,n**j).GT.(n**j-n)) )
        ELSE
           Flag = ( i.GT.(n**j-n**2) )
        END IF

        ind = (i-1)*ucl + ucl_d(j) + nu - 1

        matr(ind, ind) = 0.0D0 

        IF(nu==1) THEN
           matr(ind, ind - ucl_d(j) +1) = 1.0D0
           matr(ind - ucl_d(j) +1, ind) = 1.0D0
        ELSE
           matr(ind, ind - 1) = 1.0D0
           matr(ind - 1, ind) = 1.0D0
        END IF

        IF(Flag)THEN
           matr(ind, (i-1)*ucl + 1 - (n-1)*ucl*(n)**(j-1)) = 1.0D0
           matr((i-1)*ucl + 1 - (n-1)*ucl*(n)**(j-1), ind) = 1.0D0
        ELSE
           matr(ind, (i-1)*ucl + 1 + ucl*(n)**(j-1)) = 1.0D0
           matr((i-1)*ucl + 1 + ucl*(n)**(j-1), ind) = 1.0D0
        END IF
     END DO

  END DO
!!!!! End construct Lieb Matrix

  RETURN

END Subroutine MAKELIEBMATRIXSTRUCTRUE


!---random number
Subroutine init_random_seed()

  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)

END Subroutine init_random_seed












