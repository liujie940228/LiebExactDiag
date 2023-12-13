Subroutine MakeLiebMatrixStructrue(dm, nu, n, ucl, n_uc, nt, matr)

  USE MyNumbers
  USE IChannels
  
  IMPLICIT NONE

  INTEGER(KIND=IKIND) &
       dm, & ! the dimension
       n, &  ! the number of unit cell in each dimension
       nu, & ! the number of site in each edge
       ucl, & ! the number of atoms in a unit cell
       nt, &  ! the whole number of atoms in system
       n_uc   ! the number of unit cell

  INTEGER(KIND=IKIND) i, j, k, ind

  LOGICAL(KIND=8) Flag
  REAL(KIND=RKIND) randhopping

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
        CALL RANDOM_NUMBER(randhopping)
        randhopping = randhopping +0.5D0
        matr(ind, (i-1)*ucl + ucl_d(j)) = randhopping
        matr((i-1)*ucl + ucl_d(j), ind) = randhopping

        IF(Flag)THEN
           CALL RANDOM_NUMBER(randhopping)
           randhopping = randhopping +0.5D0
           matr(ind, (i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1) + ucl*(n)**j) = randhopping
           matr((i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1) + ucl*(n)**j, ind) = randhopping
        ELSE
           CALL RANDOM_NUMBER(randhopping)
           randhopping = randhopping +0.5D0
           matr(ind, (i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1)) = randhopping
           matr((i-1)*ucl + ucl_d(j) + (nu-1) - ucl*(n)**(j-1), ind) = randhopping
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
                 CALL RANDOM_NUMBER(randhopping)
                 randhopping = randhopping +0.5D0
                 matr(ind, ind - (j-1)*nu -1) = randhopping
                 matr(ind - (j-1)*nu -1, ind) = randhopping              
              ELSE
                 CALL RANDOM_NUMBER(randhopping)
                 randhopping = randhopping +0.5D0
                 matr(ind, ind - 1) = randhopping
                 matr(ind - 1, ind) = randhopping
              END IF
              CALL RANDOM_NUMBER(randhopping)
              randhopping = randhopping +0.5D0
              matr(ind, ind + 1) = randhopping
              matr(ind + 1, ind) = randhopping

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
           CALL RANDOM_NUMBER(randhopping)
           randhopping = randhopping +0.5D0
           matr(ind, ind - ucl_d(j) +1) = randhopping
           matr(ind - ucl_d(j) +1, ind) = randhopping
        ELSE
           CALL RANDOM_NUMBER(randhopping)
           randhopping = randhopping +0.5D0
           matr(ind, ind - 1) = randhopping
           matr(ind - 1, ind) = randhopping
        END IF

        IF(Flag)THEN
           CALL RANDOM_NUMBER(randhopping)
           randhopping = randhopping +0.5D0
           matr(ind, (i-1)*ucl + 1 - (n-1)*ucl*(n)**(j-1)) = randhopping
           matr((i-1)*ucl + 1 - (n-1)*ucl*(n)**(j-1), ind) = randhopping
        ELSE
           CALL RANDOM_NUMBER(randhopping)
           randhopping = randhopping +0.5D0
           matr(ind, (i-1)*ucl + 1 + ucl*(n)**(j-1)) = randhopping
           matr((i-1)*ucl + 1 + ucl*(n)**(j-1), ind) = randhopping
        END IF
     END DO

  END DO
!!!!! End construct Lieb Matrix

  RETURN

END Subroutine MAKELIEBMATRIXSTRUCTRUE





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


