RECURSIVE SUBROUTINE PDESource(S,Q) 
  IMPLICIT NONE
    INTEGER, PARAMETER :: nDim = 3                   ! The number of space dimensions
    INTEGER, PARAMETER :: nVar = 59                           ! The number of variables of the PDE system 
    INTEGER, PARAMETER :: nParam=0
    INTEGER, PARAMETER :: d=3
  !USE MainVariables, :ONLY:  nVar,nDim,EQN
  !USE iso_c_binding
  TYPE, bind(C) :: tEquations
      REAL(8)    :: gamma, Pi, c0, g = 9.81, friction = 1.0     
      REAL(8)    :: CCZ4k1=0.1, CCZ4k2=0.1, CCZ4k3=0.1, CCZ4eta=1.0, CCZ4itau=1.0, CCZ4f=0.1, CCZ4g, CCZ4xi=0.1, CCZ4e=1, CCZ4c=1.0, CCZ4mu=0.2, CCZ4ds=1.0, CCZ4sk=0.1, CCZ4bs=0.1  
      REAL(8)    :: CCZ4GLMc0 = 0.5, CCZ4GLMc = 0.75, CCZ4GLMd = 0.75, CCZ4GLMepsD = 1e-2, CCZ4GLMepsA = 1e-2, CCZ4GLMepsP = 1e-2, cs, alpha, beta, lambda, cv, rho0, p0, tau1, tau2, mu, kappa ,tau 
      INTEGER :: CCZ4LapseType=0, EinsteinAutoAux = 0, ReferenceDepth = 1.0    
      REAL(8)    :: DivCleaning_a = 1.0 
  END TYPE tEquations
  TYPE(tEquations) :: EQN
  ! --------------------------------------------
  ! Argument list declaration
  REAL(8),INTENT(IN)  :: Q(nvar) 
  REAL(8),INTENT(OUT) :: S(nvar)
  ! --------------------------------------------
  ! Local variables
  INTEGER :: i,j,k,l,m,n
  REAL(8) :: mu
  REAL(8) :: src(nVar)
  REAL(8) :: k1,k2,k3,fff,ggg,e,c,ds,xi,sk,sknl,bs,dgup(3,3,3),eta,itau
  REAL(8) :: fa,faa
  ! Q and their gradients
  REAL(8) :: g_cov(3,3),g_contr(3,3),det                                ! Q(1-6)
  REAL(8) :: Aex(3,3),Amix(3,3),Aup(3,3),Aupdown,traceA                 ! Q(7-12)
  REAL(8) :: Theta                                                      ! Q(13)
  REAL(8) :: Ghat(3)                                                    ! Q(14-16)
  REAL(8) :: alpha                                                      ! Q(17)
  REAL(8) :: beta(3)                                                    ! Q(18-20)
  REAL(8) :: b(3)                                                       ! Q(21-23)
  REAL(8) :: AA(3)                                                      ! Q(24-26)
  REAL(8) :: BB(3,3),traceB                                             ! Q(27-35)
  REAL(8) :: DD(3,3,3)                                                  ! Q(36-53)
  REAL(8) :: traceK                                                     ! Q(54)
  REAL(8) :: phi                                                        ! Q(55)
  REAL(8) :: PP(3)                                                      ! Q(56-58)
  REAL(8) :: K0                                                         ! Q(59)
  ! time derivatives of Q
  REAL(8) :: dtgamma(3,3)                                               ! Q(1-7)
  REAL(8) :: dtK(3,3)                                                   ! Q(7-12)
  REAL(8) :: dtTheta                                                    ! Q(13)
  REAL(8) :: dtGhat(3)                                                  ! Q(14-16)
  REAL(8) :: dtalpha                                                    ! Q(17)
  REAL(8) :: dtbeta(3)                                                  ! Q(18-20)
  REAL(8) :: dtbb(3)                                                    ! Q(21-23)
  REAL(8) :: dtA(3)                                                     ! Q(24-26)
  REAL(8) :: dtB(3,3)                                                   ! Q(27-35)
  REAL(8) :: dtD(3,3,3)                                                 ! Q(36-53)
  REAL(8) :: dtTraceK                                                   ! Q(54)
  REAL(8) :: dtphi                                                      ! Q(55)
  REAL(8) :: dtP(3)                                                     ! Q(56-58)
  ! intermediate quantities
  REAL(8) :: Gtilde(3)
  REAL(8) :: Kex(3,3),Kmix(3,3),Kup(3,3)
  REAL(8) :: Christoffel(3,3,3),Christoffel_tilde(3,3,3)
  REAL(8) :: dChristoffelSrc(3,3,3,3)
  REAL(8) :: dChristoffel_tildeSrc(3,3,3,3)
  REAL(8) :: RiemannSrc(3,3,3,3)
  REAL(8) :: RicciSrc(3,3)
  REAL(8) :: RSrc
  REAL(8) :: dGtildeSrc(3,3)
  REAL(8) :: Z(3),Zup(3),dZ(3,3),dZSrc(3,3),nablaZSrc(3,3)
  REAL(8) :: RicciPlusNablaZSrc(3,3)
  REAL(8) :: RPlusTwoNablaZSrc
  REAL(8) :: nablaijalphaSrc(3,3)
  REAL(8) :: nablanablaalphaSrc
  REAL(8) :: SecondOrderTermsSrc(3,3),traceSrc
  REAL(8) :: ov(3)

#if defined(CCZ4EINSTEIN) || defined(CCZ4GRHD) || defined(CCZ4GRMHD) || defined(CCZ4GRGPR)
  !
  k1   = EQN%CCZ4k1                             ! kappa_1
  k2   = EQN%CCZ4k2                             ! kappa_2
  k3   = EQN%CCZ4k3                             ! kappa_3
  fff  = EQN%CCZ4f                              ! multiplied to \partial_k(b^i) in the evolution eqn for BB
  ggg  = EQN%CCZ4g                              !
  eta  = EQN%CCZ4eta                            ! eta
  itau = EQN%CCZ4itau                           ! tau^-1
  e    = EQN%CCZ4e                              ! e
  c    = EQN%CCZ4c                              ! c
  mu   = EQN%CCZ4mu                             ! mu
  ds   = EQN%CCZ4ds                             ! only multiplied to Z_i and \nabla_i(Z_j)
  bs   = EQN%CCZ4bs                             ! only used in dtbb and dtB
  xi   = EQN%CCZ4xi                             ! only used in dtbb
  sk   = EQN%CCZ4sk                             ! s multiplied to \partial_i(B_j^k) and other places
  !
  ! these are the tilde quantities, so be careful!
  g_cov(1,1) = Q(1)
  g_cov(1,2) = Q(2)
  g_cov(1,3) = Q(3)
  g_cov(2,1) = Q(2)
  g_cov(2,2) = Q(4)
  g_cov(2,3) = Q(5)
  g_cov(3,1) = Q(3)
  g_cov(3,2) = Q(5)
  g_cov(3,3) = Q(6)
  ! this determinant should be unity, since we use the conformal decomposition
  det = (Q(1)*Q(4)*Q(6)-Q(1)*Q(5)**2-Q(2)**2*Q(6)+2*Q(2)*Q(3)*Q(5)-Q(3)**2*Q(4))
  !
  g_contr(1,1) =  ( g_cov(2,2)*g_cov(3,3)-g_cov(2,3)*g_cov(3,2))/det
  g_contr(1,2) = -( g_cov(1,2)*g_cov(3,3)-g_cov(1,3)*g_cov(3,2))/det
  g_contr(1,3) = -(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/det
  g_contr(2,1) = -( g_cov(2,1)*g_cov(3,3)-g_cov(2,3)*g_cov(3,1))/det
  g_contr(2,2) =  ( g_cov(1,1)*g_cov(3,3)-g_cov(1,3)*g_cov(3,1))/det
  g_contr(2,3) = -( g_cov(1,1)*g_cov(2,3)-g_cov(1,3)*g_cov(2,1))/det
  g_contr(3,1) = -(-g_cov(2,1)*g_cov(3,2)+g_cov(2,2)*g_cov(3,1))/det
  g_contr(3,2) = -( g_cov(1,1)*g_cov(3,2)-g_cov(1,2)*g_cov(3,1))/det
  g_contr(3,3) =  ( g_cov(1,1)*g_cov(2,2)-g_cov(1,2)*g_cov(2,1))/det
  !
  alpha = DEXP(DMAX1(-20.D0,DMIN1(20.D0,Q(17))))
  SELECT CASE(EQN%CCZ4LapseType)
  CASE(0)  ! harmonic
    fa  = 1.0D0
    faa = 0.0D0
  CASE DEFAULT  ! 1 + log
    fa  = 2.0D0/alpha
    faa = -2.0D0/alpha**2
  END SELECT
  !
  K0  = Q(59)
  !
  Aex(1,1) = Q(7 )
  Aex(1,2) = Q(8 )
  Aex(1,3) = Q(9 )
  Aex(2,1) = Q(8 )
  Aex(2,2) = Q(10)
  Aex(2,3) = Q(11)
  Aex(3,1) = Q(9 )
  Aex(3,2) = Q(11)
  Aex(3,3) = Q(12)
  !
  traceA = SUM(g_contr*Aex)
  Aex = Aex-1.0D0/3.0D0*g_cov*traceA
  !
  Amix = MATMUL(g_contr,          Aex  )
  Aup  = MATMUL(g_contr,TRANSPOSE(Amix))
  !
  Theta  = Q(13)
  !
  Ghat = (/Q(14),Q(15),Q(16)/)
  !
  b = Q(21:23)
  !
  AA = (/Q(24),Q(25),Q(26)/)
  !
  traceK = Q(54)
  !
  phi = DEXP(DMAX1(-20.0D0,DMIN1(20.0D0,Q(55))))
  !
  PP  = Q(56:58)
  !
  beta = (/Q(18),Q(19),Q(20)/)
  BB(1,1) = Q(27)
  BB(2,1) = Q(28)
  BB(3,1) = Q(29)
  BB(1,2) = Q(30)
  BB(2,2) = Q(31)
  BB(3,2) = Q(32)
  BB(1,3) = Q(33)
  BB(2,3) = Q(34)
  BB(3,3) = Q(35)
  !
  DD(1,1,1)=Q(36)
  DD(1,1,2)=Q(37)
  DD(1,1,3)=Q(38)
  DD(1,2,1)=Q(37)
  DD(1,2,2)=Q(39)
  DD(1,2,3)=Q(40)
  DD(1,3,1)=Q(38)
  DD(1,3,2)=Q(40)
  DD(1,3,3)=Q(41)
  !
  DD(2,1,1)=Q(42)
  DD(2,1,2)=Q(43)
  DD(2,1,3)=Q(44)
  DD(2,2,1)=Q(43)
  DD(2,2,2)=Q(45)
  DD(2,2,3)=Q(46)
  DD(2,3,1)=Q(44)
  DD(2,3,2)=Q(46)
  DD(2,3,3)=Q(47)
  !
  DD(3,1,1)=Q(48)
  DD(3,1,2)=Q(49)
  DD(3,1,3)=Q(50)
  DD(3,2,1)=Q(49)
  DD(3,2,2)=Q(51)
  DD(3,2,3)=Q(52)
  DD(3,3,1)=Q(50)
  DD(3,3,2)=Q(52)
  DD(3,3,3)=Q(53)
  !
  dgup = 0.0D0
  DO n=1,3
    DO j=1,3
      DO l=1,3
        DO m=1,3
          DO k=1,3
            dgup(k,m,l) = dgup(k,m,l)-2.0D0*g_contr(m,n)*g_contr(j,l)*DD(k,n,j)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  Kex  = Aex/phi**2+1.0D0/3.0D0*traceK*g_cov/phi**2
  Kmix = MATMUL(phi**2*g_contr,          Kex  )
  Kup  = MATMUL(phi**2*g_contr,TRANSPOSE(Kmix))
  !
  Christoffel_tilde = 0.0D0
  Christoffel       = 0.0D0
  Gtilde = 0.0D0
  !
  DO k=1,3
    DO j=1,3
      DO i=1,3
        DO l=1,3
          Christoffel_tilde(i,j,k) = Christoffel_tilde(i,j,k)+g_contr(k,l)*(DD(i,j,l)+DD(j,i,l)-DD(l,i,j))
          Christoffel(i,j,k)       = Christoffel(i,j,k)      +g_contr(k,l)*(DD(i,j,l)+DD(j,i,l)-DD(l,i,j)) &
                  &                         -g_contr(k,l)*(g_cov(j,l)*PP(i)+g_cov(i,l)*PP(j)-g_cov(i,j)*PP(l))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  DO l=1,3
    DO j=1,3
      DO i=1,3
        Gtilde(i) = Gtilde(i)+g_contr(j,l)*Christoffel_tilde(j,l,i)
      ENDDO
    ENDDO
  ENDDO
  !
  Z   = 0.5D0*ds*MATMUL(g_cov,Ghat-Gtilde)
  Zup = MATMUL(phi**2*g_contr,Z)
  !
  dChristoffelSrc = 0.0D0
  dChristoffel_tildeSrc = 0.0D0
  DO l=1,3
    DO m=1,3
      DO j=1,3
        DO i=1,3
          DO k=1,3
            dChristoffelSrc(k,i,j,m) = dChristoffelSrc(k,i,j,m)+dgup(k,m,l)*(DD(i,j,l)+DD(j,i,l)-DD(l,i,j))                      &
                    &                         -dgup(k,m,l)*(g_cov(j,l)*PP(i)+g_cov(i,l)*PP(j)-g_cov(i,j)*PP(l)) &
                    &                         -2.0D0*g_contr(m,l)*(DD(k,j,l)*PP(i)+DD(k,i,l)*PP(j)-DD(k,i,j)*PP(l))
            !
            dChristoffel_tildeSrc(k,i,j,m) = dChristoffel_tildeSrc(k,i,j,m)+dgup(k,m,l)*(DD(i,j,l)+DD(j,i,l)-DD(l,i,j))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  RiemannSrc = 0.0D0
  DO m=1,3
    DO j=1,3
      DO k=1,3
        DO i=1,3
          RiemannSrc(i,k,j,m) = dChristoffelSrc(k,i,j,m)-dChristoffelSrc(j,i,k,m)
          DO l=1,3
            RiemannSrc(i,k,j,m) = RiemannSrc(i,k,j,m)+Christoffel(i,j,l)*Christoffel(l,k,m)-Christoffel(i,k,l)*Christoffel(l,j,m)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  RicciSrc = 0.0D0
  DO l=1,3
    DO n=1,3
      DO m=1,3
        RicciSrc(m,n) = RicciSrc(m,n)+RiemannSrc(m,l,n,l)
      ENDDO
    ENDDO
  ENDDO
  !
  RSrc = phi**2*SUM(g_contr*RicciSrc)
  !
  ! Here we directly compute the derivative of Gtilde from its original definition as contracted Christoffel symbol,
  ! without assuming unit determinant of the conformal metric. Back to the roots, and as few assumptions as possible...
  !
  dGtildeSrc = 0.0D0
  DO l=1,3
    DO j=1,3
      DO i=1,3
        DO k=1,3
          dGtildeSrc(k,i) = dGtildeSrc(k,i)+dgup(k,j,l)*Christoffel_tilde(j,l,i)+g_contr(j,l)*dChristoffel_tildeSrc(k,j,l,i)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  dZSrc = 0.0D0
  DO j=1,3
    DO i=1,3
      DO k=1,3
        dZSrc(k,i) = dZSrc(k,i)+ds*(DD(k,i,j)*(Ghat(j)-Gtilde(j))+0.5D0*g_cov(i,j)*(-dGtildeSrc(k,j)))
      ENDDO
    ENDDO
  ENDDO
  !
  nablaZSrc = 0.0D0
  DO j=1,3
    DO i=1,3
      nablaZSrc(i,j) = dZSrc(i,j)
      DO k=1,3
        nablaZSrc(i,j) = nablaZSrc(i,j)-Christoffel(i,j,k)*Z(k)
      ENDDO
    ENDDO
  ENDDO
  !
  RicciPlusNablaZSrc = RicciSrc+(nablaZSrc+TRANSPOSE(nablaZSrc))
  !
  RPlusTwoNablaZSrc = phi**2*SUM(g_contr*RicciPlusNablaZSrc)
  !
  nablaijalphaSrc = 0.0D0
  DO j=1,3
    DO i=1,3
      nablaijalphaSrc(i,j) =      alpha*AA(j)*AA(i)
      DO k=1,3
        nablaijalphaSrc(i,j) = nablaijalphaSrc(i,j)-alpha*Christoffel(i,j,k)*AA(k)
      ENDDO
    ENDDO
  ENDDO
  nablanablaalphaSrc = phi**2*SUM(g_contr*nablaijalphaSrc)
  !
  SecondOrderTermsSrc = -nablaijalphaSrc+alpha*RicciPlusNablaZSrc
  traceSrc = SUM(g_contr*SecondOrderTermsSrc)
  SecondOrderTermsSrc = SecondOrderTermsSrc-1.0D0/3.0D0*g_cov*traceSrc
  !
  ! now assemble all this terrible stuff...
  !
  dtgamma = -2.0D0*alpha*Aex-itau*(det-1.0D0)*g_cov
  DO k=1,3
    DO j=1,3
      DO i=1,3
        dtgamma(i,j) = dtgamma(i,j)+g_cov(k,i)*BB(j,k)+g_cov(k,j)*BB(i,k)-2.0D0/3.0D0*g_cov(i,j)*BB(k,k)+2.0D0*beta(k)*DD(k,i,j)
      ENDDO
    ENDDO
  ENDDO
  !
  ! Main variables of the CCZ4 system
  ! extrinsic curvature
  dtK = phi**2*SecondOrderTermsSrc + alpha*Aex*(traceK-2.0D0*Theta) - 2.0D0*alpha*MATMUL(Aex,Amix)-itau*g_cov*traceA
  DO j=1,3
    DO i=1,3
      DO k=1,3
        dtK(i,j) = dtK(i,j)+Aex(k,i)*BB(j,k)+Aex(k,j)*BB(i,k)-2.0D0/3.0D0*Aex(i,j)*BB(k,k)
      ENDDO
    ENDDO
  ENDDO
  !
  ! dtTraceK = -nablanablaalphaNCP-nablanablaalphaSrc+alpha*(RPlusTwoNablaZNCP+RPlusTwoNablaZSrc+traceK**2-2.0D0*Theta*traceK)-3.0D0*alpha*k1*(1.0D0+k2)*Theta+SUM(beta(:)*dtraceK(:))   ! Baojiu
  dtTraceK = -nablanablaalphaSrc+alpha*(RPlusTwoNablaZSrc+traceK**2-2.0D0*c*Theta*traceK)-3.0D0*alpha*k1*(1.0D0+k2)*Theta ! Baojiu
  !
  traceB  = BB(1,1)+BB(2,2)+BB(3,3)
  dtphi   = beta(1)*PP(1)+beta(2)*PP(2)+beta(3)*PP(3)+1.0D0/3.0D0*alpha*traceK-1.0D0/3.0D0*traceB
  dtalpha = -alpha*fa*(traceK-K0-2.0D0*c*Theta)+beta(1)*AA(1)+beta(2)*AA(2)+beta(3)*AA(3)
  !
  Aupdown = SUM(Aex*Aup)
  ! *** original
  dtTheta = 0.5D0*alpha*e**2*RplusTwoNablaZSrc+    &            ! temporal Z
          !         & 0.5D0*alpha*e**2*(-Aupdown+2.0D0/3.0D0*traceK**2)-alpha*Theta*traceK-SUM(Zup*alpha*AA)-alpha*k1*(2.0D0+k2)*Theta    ! Baojiu
          & 0.5D0*alpha*e**2*(-Aupdown+2.0D0/3.0D0*traceK**2)-c*alpha*Theta*traceK-alpha*SUM(Zup*AA)-alpha*k1*(2.0D0+k2)*Theta  ! Baojiu
  !
  dtGhat = 0.0D0
  DO i=1,3
    dtGhat(i) = dtGhat(i)  &
            & +2.0D0*alpha*(SUM(Christoffel_tilde(:,:,i)*Aup(:,:))-3.0D0*SUM(Aup(i,:)*PP(:)))    &
            & +2.0D0*alpha*SUM(g_contr(:,i)*(-Theta*AA(:)-2.0D0/3.0D0*traceK*Z(:)))  &
            & -2.0D0*alpha*SUM(Aup(i,:)*AA(:))-2.0D0*alpha*k1*SUM(g_contr(i,:)*Z(:))-SUM(Gtilde(:)*BB(:,i))   &
            & +2.0D0/3.0D0*Gtilde(i)*traceB
    DO l=1,3
      DO k=1,3
        dtGhat(i) = dtGhat(i)+ &
                &  2.0D0*k3*(2.0D0/3.0D0*g_contr(i,l)*Z(l)*BB(k,k)-g_contr(l,k)*Z(l)*BB(k,i))
      ENDDO
    ENDDO
  ENDDO
  !
  DO k=1,3
    ov(k) = 2.0D0*alpha*SUM(dgup(k,:,:)*Aex(:,:))    ! here we can use the constraint that trace A tilde = 0.
  ENDDO
  !
  dtGhat = dtGhat+sk*MATMUL(g_contr,ov)                                               ! the above ordering constraint is "down", so we must raise the index via g_contr.
  !
  dtbb = xi*dtGhat-eta*b                                                            !  <= be careful, this damping term -eta*b may be dangerous for the gamma driver, since it may kill the waves that you want !
  dtbb = sk*dtbb
  !
  dtbeta  = fff*b
  ! Add the following term if you want to have shift convection in the PDE for beta^i
  ! Do not add it if you want a real Lie derivative for beta. In this case, the advection term cancels out.
  dtbeta = dtbeta+bs*(beta(1)*BB(1,:)+beta(2)*BB(2,:)+beta(3)*BB(3,:))
  dtbeta = sk*dtbeta
  !
  ! Auxiliary variables
  dtA = -alpha*AA*(fa+alpha*faa)*(traceK-K0-2.0D0*c*Theta)+MATMUL(BB,AA)
  DO k=1,3
    dtA(k) = dtA(k)-sk*alpha*fa*SUM(dgup(k,:,:)*Aex(:,:))  ! here we can use the constraint that trace A tilde = 0.
  ENDDO
  !
  dtB = sk*MATMUL(BB,BB)     ! Baojiu
  !
  dtD = 0.0D0
  DO m=1,3
    DO j=1,3
      DO i=1,3
        DO k=1,3
          DO n=1,3
            dtD(k,i,j) = dtD(k,i,j)+1.0D0/3.0D0*alpha*g_cov(i,j)*dgup(k,n,m)*Aex(n,m)   ! explicitly remove the trace of tilde A again
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  DO j=1,3
    DO i=1,3
      DO k=1,3
        dtD(k,i,j) = dtD(k,i,j)-alpha*AA(k)*Aex(i,j) !trace removing missing here
        DO l=1,3
          dtD(k,i,j) = dtD(k,i,j)+BB(k,l)*DD(l,i,j)+DD(k,l,i)*BB(j,l)+DD(k,l,j)*BB(i,l)-2.0D0/3.0D0*DD(k,i,j)*BB(l,l)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  dtP = MATMUL(BB,PP)

  
  DO k=1,3
    dtP(k) = dtP(k)+1.0D0/3.0D0*alpha*AA(k)*traceK+sk*1.0D0/3.0D0*alpha*SUM(dgup(k,:,:)*Aex(:,:))
  ENDDO

  src(1:6)   = (/dtgamma(1,1),dtgamma(1,2),dtgamma(1,3),dtgamma(2,2),dtgamma(2,3),dtgamma(3,3)/)       ! \tilde \gamma_ij
  src(7:12)  = (/dtK(1,1),dtK(1,2),dtK(1,3),dtK(2,2),dtK(2,3),dtK(3,3)/)                               ! \tilde A_ij
  src(13)    = dtTheta                                                                                 ! Theta
  src(14:16) = dtGhat(1:3)                                                                             ! \hat \Gamma^i
  src(17)    = dtalpha                                                                                 ! log alpha
  src(18:20) = dtbeta                                                                                  ! beta^i
  src(21:23) = dtbb                                                                                    ! b^i
  src(24:26) = dtA(1:3)                                                                                ! A_k
  src(27:35) = (/dtB(1,1),dtB(2,1),dtB(3,1),dtB(1,2),dtB(2,2),dtB(3,2),dtB(1,3),dtB(2,3),dtB(3,3)/)    ! B_k^i
  src(36:41) = (/dtD(1,1,1),dtD(1,1,2),dtD(1,1,3),dtD(1,2,2),dtD(1,2,3),dtD(1,3,3)/)                   ! D_kij
  src(42:47) = (/dtD(2,1,1),dtD(2,1,2),dtD(2,1,3),dtD(2,2,2),dtD(2,2,3),dtD(2,3,3)/)                   ! D_kij
  src(48:53) = (/dtD(3,1,1),dtD(3,1,2),dtD(3,1,3),dtD(3,2,2),dtD(3,2,3),dtD(3,3,3)/)                   ! D_kij
  src(54)    = dtTraceK                                                                                ! traceK
  src(55)    = dtphi                                                                                   ! log phi
  src(56:58) = dtP                                                                                     ! P_k
  src(59)    = 0.0D0
  !
  S = src    ! here, we do not have to change sign, since we work on the right hand side in the fused subroutine
  !
  RETURN
#endif
end SUBROUTINE PDESource

program main
  implicit none
  REAL(8) ::  Q(59), S(59), Qtest(59)
  integer :: i,j,k,l,m

  Qtest = (/1.03876392436862219348e+00,1.13327632951146378931e-17,1.13328125850752231240e-17,9.81163922513405450943e-01,7.28991960382758263436e-20,9.81163922513405450943e-01,1.60702820377219834924e-01,-2.51854521837995827255e-14,-2.51859943549822510890e-14,-7.58958825491163335819e-02,-1.66695209611466999131e-16,-7.58958825491163335819e-02,5.73443499779581738335e-04,3.01713236980058585601e-01,-5.73062863551328161271e-12,-5.73061639353499262517e-12,2.85203957292370774423e-02,0.00000000000000000000e+00,0.00000000000000000000e+00,0.00000000000000000000e+00,0.00000000000000000000e+00,0.00000000000000000000e+00,0.00000000000000000000e+00,2.33452764998758000026e-01,-4.17807322717328905611e-12,-4.17806121758503744557e-12,0.00000000000000000000e+00,0.00000000000000000000e+00,0.00000000000000000000e+00,0.00000000000000000000e+00,0.00000000000000000000e+00,0.00000000000000000000e+00,0.00000000000000000000e+00,0.00000000000000000000e+00,0.00000000000000000000e+00,1.66029261165760388952e-01,3.52110122269561003271e-16,3.52110623323521089469e-16,-7.84114259822868420180e-02,-1.00892874234056400197e-18,-7.84114259822868420180e-02,-2.87874071850877923045e-12,-3.11299149429065573305e-16,4.67768108082361629111e-18,1.35940313205684060638e-12,-2.07569443609550208564e-18,1.35970991086859628736e-12,-2.87874063909529684826e-12,1.43318865709343599115e-18,-3.14548568405474184678e-16,1.35970993559255365097e-12,-2.04322863271445272127e-18,1.35940303232291987362e-12,2.34302828568053289615e-01,-9.50679857641236088217e-03,-7.99171325660463643947e-02,1.44427998332769005893e-12,1.44427998345849223420e-12,0.00000000000000000000e+00/)

  ! The initial data for Q is chosen such that we do not divide by 0 in the inverse determinant
  Q=1
  Q(1)=2
  Q(2)=3

  S=0
  !do i=1,1000000
  !call PDESource(S,Q)
  call PDESource(S,Qtest)
  !enddo

  !do i=1,1000000
    !call PDENCP(BgradQ, Q, gradQin)
  !enddo
  do i=1,59
    print*,i,S(i)
  enddo
  !call PDENCP(BgradQ, Q, gradQin)
  !print*, BgradQ
end program main
