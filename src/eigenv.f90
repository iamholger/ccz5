RECURSIVE SUBROUTINE PDEEigenvalues(Lambda,Q,nv)
  USE iso_c_binding
  IMPLICIT NONE
  INTEGER, PARAMETER             :: nVar = 59                           ! The number of variables of the PDE system 
  TYPE, bind(C) :: tEquations
      REAL(8)    :: CCZ4e=0.1, CCZ4ds=0.11, CCZ4GLMc = 0.75, CCZ4GLMd = 0.75 
  END TYPE tEquations
  TYPE(tEquations) :: EQN
  REAL :: Lambda(nVar), nv(3), Q(nVar), tempA, tempB
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: Lambda 
  ! Local Variables 
  REAL :: lam, mu, irho, VPR(3), cs,c0,uu,alpha
#if defined(CCZ4EINSTEIN)  
    alpha = MAX( 1.0, EXP(Q(17)) )*MAX( 1.0, EXP(Q(55)) )/MIN( SQRT(Q(1)), SQRT(Q(4)), SQRT(Q(6)) )     
#else
    alpha = 1.0 
#endif
    tempA = alpha*MAX(SQRT(2.0), EQN%CCZ4e, EQN%CCZ4ds, EQN%CCZ4GLMc/alpha, EQN%CCZ4GLMd/alpha )
    tempB = DOT_PRODUCT(Q(18:20),nv(:))
    Lambda = 0.0 
    Lambda(1) = -tempA-tempB 
    Lambda(2) =  tempA-tempB
    !Lambda(1) = -alpha*MAX(SQRT(2.0), EQN%CCZ4e, EQN%CCZ4ds, EQN%CCZ4GLMc/alpha, EQN%CCZ4GLMd/alpha ) - DOT_PRODUCT(Q(18:20),nv(:))   ! MAX( SQRT(2.0), EQN%CCZ4e, EQN%CCZ4ds ) + SQRT(SUM(Q(18:20)**2)) 
    !Lambda(2) = +alpha*MAX(SQRT(2.0), EQN%CCZ4e, EQN%CCZ4ds, EQN%CCZ4GLMc/alpha, EQN%CCZ4GLMd/alpha ) - DOT_PRODUCT(Q(18:20),nv(:))   ! MAX( SQRT(2.0), EQN%CCZ4e, EQN%CCZ4ds ) + SQRT(SUM(Q(18:20)**2)) 
    !
END SUBROUTINE PDEEigenvalues

program main
  implicit none
  REAL :: Q(59), nv(3), Lambda(59),dummy(3)
  integer::i, seed_size
  integer,allocatable :: seed(:)
  call random_seed() ! initialize with system generated seed
  call random_seed(size=seed_size) ! find out size of seed
  allocate(seed(seed_size))
  call random_seed(get=seed) ! get system generated seed
  seed=314159261
  call random_seed(put=seed) ! set current seed
  call random_seed(get=seed) ! get current seed
  deallocate(seed)           ! safe


  call random_number(Q)

  nv = (/1,0,0/)
  dummy=(/2,0,0/)
  !print*,Q

  ! The initial data for Q is chosen such that we do not divide by 0 in the inverse determinant

  !do i=1,1000000000
  call PDEEigenvalues(Lambda, Q, nv)
  !enddo

  print*,Lambda
  print*,nv
  print*, dot_product(nv, dummy)
  !do i=1,59
    !print*,i,Lambda(i)
  !enddo
  !call PDENCP(BgradQ, Q, gradQin)
  !print*, BgradQ
end program main

