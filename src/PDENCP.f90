subroutine PDENCP(BgradQ, Q, gradQin) 
    !USE ISO_C_BINDING
!$omp declare target
   !use matrix3
   IMPLICIT NONE
   ! 11. Oct 21:40: This was a matrix BGradQ(nVar, nDim) but is a vector in spaceTimePredictorNonlinear
#if defined(Dim3)
    INTEGER, PARAMETER :: nDim = 3                   ! The number of space dimensions
#elif defined(Dim2)
    INTEGER, PARAMETER :: nDim = 2                   ! The number of space dimensions
#endif
    INTEGER, PARAMETER :: nVar = 59                           ! The number of variables of the PDE system 
    INTEGER, PARAMETER :: nParam=0
    INTEGER, PARAMETER :: d=3

   REAL, intent(out) :: BgradQ(nVar)
   REAL, intent(in) :: gradQin(nVar, 3)
   REAL, intent(in) :: Q(nVar)
  TYPE, bind(C) :: tEquations
      REAL(8)    :: gamma, Pi, c0, g = 9.81, friction = 1.0     
      REAL(8)    :: CCZ4k1, CCZ4k2, CCZ4k3, CCZ4eta, CCZ4itau, CCZ4f=0.0, CCZ4g, CCZ4xi=0.0, CCZ4e=1, CCZ4c=1.0, CCZ4mu=0.2, CCZ4ds=1.0, CCZ4sk=0.1, CCZ4bs=0.0  
      REAL(8)    :: CCZ4GLMc0 = 0.5, CCZ4GLMc = 0.75, CCZ4GLMd = 0.75, CCZ4GLMepsD = 1e-2, CCZ4GLMepsA = 1e-2, CCZ4GLMepsP = 1e-2, cs, alpha, beta, lambda, cv, rho0, p0, tau1, tau2, mu, kappa ,tau 
      INTEGER :: CCZ4LapseType=0, EinsteinAutoAux = 0, ReferenceDepth = 1.0    
      REAL(8)    :: DivCleaning_a = 1.0 
  END TYPE tEquations
  TYPE(tEquations) :: EQN

  ! Set parameters here according to Properties, case GaugeWave
  !EQN%CCZ4GLMc0   = 1.5   ! 0.1      
  !EQN%CCZ4GLMc    = 1.2   ! 2.0    
  !EQN%CCZ4GLMd    = 2.0   ! 1.0     
  !EQN%CCZ4GLMepsA = 1.0   ! 5. 
  !EQN%CCZ4GLMepsP = 1.0   ! 5.  
  !EQN%CCZ4GLMepsD = 1.0   ! 0.1 
  !!
  !EQN%CCZ4itau  = 1.0 

  !EQN%CCZ4k1  = 0.0  !modified according to the version in ExaHyPE 1
  !EQN%CCZ4k2  = 0.0 
  !EQN%CCZ4k3  = 0.0 
  !EQN%CCZ4eta = 0.0 
  !EQN%CCZ4f   = 0.0 
  !EQN%CCZ4g   = 0.0 
  !EQN%CCZ4xi  = 0.0 
  !EQN%CCZ4e   = 1.0 
  !EQN%CCZ4c   = 1.0 
  !EQN%CCZ4mu  = 0.2 
  !EQN%CCZ4ds  = 1.0 
  !EQN%CCZ4sk  = 0.0
  !EQN%CCZ4bs   = 0.0      ! set bs=1 if you want to activate the shift convection for beta, b and B (standard CCZ4 formulation). set it to bs=0 to switch off shift convection for those quantities 
  !EQN%CCZ4LapseType   = 0 ! harmonic lapse 
  !EQN%EinsteinAutoAux = 0 




   
   REAL  :: gradQ(nVar, 3)
    ! Argument list 
    REAL :: par(nParam)  
    ! Local variables 
    INTEGER :: i,j,k,l,m,n,ip,iq,ii,jj,kk,ll,iErr,count    
    REAL :: p, irho, lam, mu 
    REAL :: Qx(nVar), Qy(nVar), Qz(nVar) 
    REAL :: k1,k2,k3,fff,ggg,e,c,ds,xi,sk,sknl,bs,g_cov(3,3),g_contr(3,3),dgup(3,3,3)
    REAL :: det, alpha, fa, k0, dk0(3), beta0(3), b0(3), u(3), ialpha 
    REAL :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, eta, itau    
    REAL :: AQx(nVar), BQy(nVar), CQz(nVar)  
    REAL :: lapse, shift(3), gammaij(6), delta(3,3), bv(3), vxb(3), vxb_contr(3), psi, qv_contr(3), qb_contr(3), bv_contr(3) 
    REAL :: v2,vf(3),uem,b2,e2,gp,gm,vc(nVar),lf,w,ww,gamma1,rho,vf_cov(3), s_contr(3), w_ij, wim
    REAL :: v_cov(3) 
#ifdef VECTOR    
#ifdef AVX512 
  INTEGER, PARAMETER :: nVarGRMHD = 24                           ! The number of variables of the PDE system 
#else   
  INTEGER, PARAMETER :: nVarGRMHD = 20                           ! The number of variables of the PDE system 
#endif 
#else
  INTEGER, PARAMETER :: nVarGRMHD = 19                           ! The number of variables of the PDE system 
#endif
    REAL :: QGRMHD(nVarGRMHD), gradQGRMHD(nVarGRMHD,d), BgradQGRMHD(nVarGRMHD), ov(3) 
#ifdef VECTOR    
      INTEGER, PARAMETER :: nVarGRGPR = 32                           ! The number of variables of the PDE system 
#else
      INTEGER, PARAMETER :: nVarGRGPR = 30                           ! The number of variables of the PDE system 
#endif
    REAL :: QGRGPR(nVarGRGPR), gradQGRGPR(nVarGRGPR,d), BgradQGRGPR(nVarGRGPR)
    !
    REAL :: Christoffel(3,3,3), RiemannNCP(3,3,3,3), RiemannSrc(3,3,3,3), dChristoffelNCP(3,3,3,3), dChristoffel_tildeNCP(3,3,3,3), dChristoffelSrc(3,3,3,3), DD(3,3,3), dDD(3,3,3,3)  
    REAL :: AA(3), dAA(3,3), BB(3,3), dBB(3,3,3), beta(3), Kex(3,3), Kmix(3,3), Kup(3,3), Z(3), dZ(3,3), nablaZNCP(3,3), nablaZSrc(3,3), RplusNablaZNCP, RplusNablaZSrc  
    REAL :: Theta, dTheta(3), nablaijalphaNCP(3,3), nablaijalphaSrc(3,3), Ricci(3,3), RicciNCP(3,3), RicciSrc(3,3), dtraceK(3), dtraceKNCP(3), dKtempNCP(3), dZNCP(3,3), dZSrc(3,3)  
    REAL :: dtgamma(3,3), dtK(3,3), dK(3,3,3), dtTheta, dtZ(3), dtalpha, dtGhat(3), dtbeta(3), dtbb(3), dxbb(3,3), dtA(3), dtB(3,3), dtD(3,3,3)  
    REAL :: Aupdown, Aex(3,3), dAex(3,3,3), Amix(3,3), Aup(3,3), Ghat(3), Gtilde(3), dGhat(3,3), traceK, Kupdown, phi, phi2, PP(3), dPP(3,3), TwoNablaZNCP, TwoNablaZSrc, dKex(3,3,3)      
    REAL :: dGtildeSrc(3,3), dGtildeNCP(3,3), RiccitildeNCP(3,3), RicciphiNCP(3,3), RiccitildeSrc(3,3), RicciphiSrc(3,3), Mom(3), Ham, Pup(3), DDcontr(3)   
    REAL :: Christoffel_tilde(3,3,3), Christoffel_kind1(3,3,3), Zup(3), RicciPlusNablaZNCP(3,3), RicciPlusNablaZSrc(3,3), traceA, traceB, QG(3), b(3), faa, temp   
    REAL :: SecondOrderTermsNCP(3,3), SecondOrderTermsSrc(3,3), traceNCP, traceSrc, dtphi, dtTraceK, dtP(3), dtX(3), XX(3), dXX(3,3), nablaXNCP(3,3)     
    REAL :: RPlusTwoNablaZNCP, RNCP, RSrc, RPlusTwoNablaZSrc, nablanablaalpha, nablanablaalphaNCP, nablanablaalphaSrc, Riemann(3,3,3,3), dChristoffel(3,3,3,3), divAupNCP(3), divAupSrc(3) 
    !
#ifdef GLMROT     
    REAL(8) :: dpsiA(3,3), dpsiP(3,3), dpsiB(3,3,3), dpsiD(3,3,3,3), dphiA(3), dphiP(3), dphiB(3,3), dphiD(3,3,3)   
    REAL(8) :: dtpsiA(3), dtpsiP(3), dtpsiB(3,3), dtpsiD(3,3,3), dtphiA, dtphiP, dtphiB(3), dtphiD(3,3)    
#endif         
    !
    BgradQ = 0.0

#if defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) || defined(CCZ4GRHD) || defined(CCZ4GRGPR) 
    !
    Qx = gradQin(:,1) 
    Qy = gradQin(:,2)
#if defined(Dim2)
    Qz = 0.0
    gradQ(:,1:2)=gradQin(:,1:2)
    gradQ(:,3)=0.0
#elif defined(Dim3)
    Qz = gradQin(:,3)
    gradQ=gradQin
#endif


    k1   = EQN%CCZ4k1  
    k2   = EQN%CCZ4k2  
    k3   = EQN%CCZ4k3   
    fff  = EQN%CCZ4f 
    ggg  = EQN%CCZ4g 
    e    = EQN%CCZ4e 
    itau = EQN%CCZ4itau  
    eta  = EQN%CCZ4eta  
    c    = EQN%CCZ4c 
    mu   = EQN%CCZ4mu 
    ds   = EQN%CCZ4ds 
    bs   = EQN%CCZ4bs 
    xi   = EQN%CCZ4xi 
    sk   = EQN%CCZ4sk
    !
    ! These are the tilde quantities, so be careful !    
    g_cov(1,1) = Q(1)
    g_cov(1,2) = Q(2)
    g_cov(1,3) = Q(3)
    g_cov(2,1) = Q(2)
    g_cov(2,2) = Q(4)
    g_cov(2,3) = Q(5)
    g_cov(3,1) = Q(3)
    g_cov(3,2) = Q(5)
    g_cov(3,3) = Q(6)
    ! This determinant should be close to unity, since we use the conformal decomposition 
    ! NEED to ensure that this is not 0
    det = 1./(Q(1)*Q(4)*Q(6)-Q(1)*Q(5)**2-Q(2)**2*Q(6)+2*Q(2)*Q(3)*Q(5)-Q(3)**2*Q(4))  ! HS inverse det really
    !print*,"DET",det
    
    g_contr(1,1) =  ( g_cov(2,2)*g_cov(3,3)-g_cov(2,3)*g_cov(3,2)) 
    g_contr(1,2) = -( g_cov(1,2)*g_cov(3,3)-g_cov(1,3)*g_cov(3,2))
    g_contr(1,3) = -(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2)) 
    g_contr(2,1) = -( g_cov(2,1)*g_cov(3,3)-g_cov(2,3)*g_cov(3,1)) 
    g_contr(2,2) =  ( g_cov(1,1)*g_cov(3,3)-g_cov(1,3)*g_cov(3,1)) 
    g_contr(2,3) = -( g_cov(1,1)*g_cov(2,3)-g_cov(1,3)*g_cov(2,1)) 
    g_contr(3,1) = -(-g_cov(2,1)*g_cov(3,2)+g_cov(2,2)*g_cov(3,1)) 
    g_contr(3,2) = -( g_cov(1,1)*g_cov(3,2)-g_cov(1,2)*g_cov(3,1)) 
    g_contr(3,3) =  ( g_cov(1,1)*g_cov(2,2)-g_cov(1,2)*g_cov(2,1)) 
     
    g_contr = g_contr * det

    alpha = EXP(MAX(-20.,MIN(20.,Q(17))))  
    !print*,"alpha2",alpha*alpha
    SELECT CASE(EQN%CCZ4LapseType) 
    CASE(0)  ! harmonic 
        fa = 1.0 
        faa = 0.0 
    CASE DEFAULT  ! 1 + log 
        fa = 2.0/alpha
        faa = -2.0/alpha**2   
    END SELECT 
    ! 
    K0    = Q(59)
    dK0   = 0.0 ! sk*gradQ(59,:) 
    !  
    Aex(1,1) = Q(7) 
    Aex(1,2) = Q(8) 
    Aex(1,3) = Q(9) 
    Aex(2,1) = Q(8) 
    Aex(2,2) = Q(10) 
    Aex(2,3) = Q(11) 
    Aex(3,1) = Q(9) 
    Aex(3,2) = Q(11) 
    Aex(3,3) = Q(12) 
    !print*,"Aex",Aex
    !
    traceA = SUM(g_contr*Aex)
    !print*,"traceA",traceA
    Aex = Aex - 1./3.*g_cov*traceA 
    !
    dAex(:,1,1) = gradQ(7,:) 
    dAex(:,1,2) = gradQ(8,:) 
    dAex(:,1,3) = gradQ(9,:) 
    dAex(:,2,1) = gradQ(8,:) 
    dAex(:,2,2) = gradQ(10,:) 
    dAex(:,2,3) = gradQ(11,:) 
    dAex(:,3,1) = gradQ(9,:) 
    dAex(:,3,2) = gradQ(11,:) 
    dAex(:,3,3) = gradQ(12,:) 
    !
    Amix = matmul(g_contr, Aex)
    Aup  = matmul(g_contr, transpose(Amix)) 
    !
    Theta = Q(13)
    dTheta = gradQ(13,:) 
    ! 
    Ghat = (/ Q(14), Q(15), Q(16) /)
    dGhat(:,1) = gradQ(14,:)
    dGhat(:,2) = gradQ(15,:)
    dGhat(:,3) = gradQ(16,:)
    !
    b = Q(21:23) 
    !
    AA    = (/ Q(24), Q(25), Q(26) /) 
    dAA(:,1) = gradQ(24,:) 
    dAA(:,2) = gradQ(25,:) 
    dAA(:,3) = gradQ(26,:) 
    !print*,"dAA",dAA
    !
    traceK = Q(54) 
    dtraceK = gradQ(54,:) 
    !
    phi   = EXP(MAX(-20.,MIN(20.,Q(55))))
    phi2 =  phi*phi
    !print*,"phi2",phi2

    PP    = Q(56:58) 
    dPP(:,1) = gradQ(56,:) 
    dPP(:,2) = gradQ(57,:) 
    dPP(:,3) = gradQ(58,:) 
    !
    beta = (/ Q(18), Q(19), Q(20) /) 
    !print*,"beta",beta
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
    dBB(:,1,1) = sk*gradQ(27,:) 
    dBB(:,2,1) = sk*gradQ(28,:) 
    dBB(:,3,1) = sk*gradQ(29,:) 
    dBB(:,1,2) = sk*gradQ(30,:) 
    dBB(:,2,2) = sk*gradQ(31,:) 
    dBB(:,3,2) = sk*gradQ(32,:) 
    dBB(:,1,3) = sk*gradQ(33,:) 
    dBB(:,2,3) = sk*gradQ(34,:) 
    dBB(:,3,3) = sk*gradQ(35,:) 
    !
    !dBB = dBB*sk    
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
    dDD(:,1,1,1)=gradQ(36,:) 
    dDD(:,1,1,2)=gradQ(37,:) 
    dDD(:,1,1,3)=gradQ(38,:) 
    dDD(:,1,2,1)=gradQ(37,:) 
    dDD(:,1,2,2)=gradQ(39,:) 
    dDD(:,1,2,3)=gradQ(40,:)
    dDD(:,1,3,1)=gradQ(38,:) 
    dDD(:,1,3,2)=gradQ(40,:) 
    dDD(:,1,3,3)=gradQ(41,:)
    dDD(:,2,1,1)=gradQ(42,:) 
    dDD(:,2,1,2)=gradQ(43,:) 
    dDD(:,2,1,3)=gradQ(44,:) 
    dDD(:,2,2,1)=gradQ(43,:) 
    dDD(:,2,2,2)=gradQ(45,:) 
    dDD(:,2,2,3)=gradQ(46,:)
    dDD(:,2,3,1)=gradQ(44,:) 
    dDD(:,2,3,2)=gradQ(46,:) 
    dDD(:,2,3,3)=gradQ(47,:) 
    dDD(:,3,1,1)=gradQ(48,:) 
    dDD(:,3,1,2)=gradQ(49,:) 
    dDD(:,3,1,3)=gradQ(50,:) 
    dDD(:,3,2,1)=gradQ(49,:) 
    dDD(:,3,2,2)=gradQ(51,:) 
    dDD(:,3,2,3)=gradQ(52,:)
    dDD(:,3,3,1)=gradQ(50,:) 
    dDD(:,3,3,2)=gradQ(52,:) 
    dDD(:,3,3,3)=gradQ(53,:)



    !
    dgup = 0.0 
    DO k = 1, 3 
     DO m = 1, 3 
      DO l = 1, 3 
       DO n = 1, 3
        DO j = 1, 3 
           dgup(k,m,l) = dgup(k,m,l)-g_contr(m,n)*g_contr(j,l)*2*DD(k,n,j) 
        ENDDO
       ENDDO 
      ENDDO 
     ENDDO 
    ENDDO         
    !
    Kex  = Aex/phi2 + 1./3.*traceK*g_cov/phi2 
    Kmix = matmul( phi2*g_contr, Kex  ) 
    Kup  = matmul( phi2*g_contr, transpose(Kmix)) 
    !
    Christoffel_tilde = 0.0  
    Christoffel       = 0.0 
    Gtilde = 0.0 
    !
    DO j = 1, 3
     DO i = 1, 3
      DO k = 1, 3
       Christoffel_kind1(i,j,k) = DD(k,i,j)+DD(j,i,k)-DD(i,j,k)      ! this definition seems to work ! 
       DO l = 1, 3
          Christoffel_tilde(i,j,k) = Christoffel_tilde(i,j,k) + g_contr(k,l)*( DD(i,j,l)+DD(j,i,l)-DD(l,i,j) ) 
          Christoffel(i,j,k)       = Christoffel(i,j,k)       + g_contr(k,l)*( DD(i,j,l)+DD(j,i,l)-DD(l,i,j) ) -g_contr(k,l)*( g_cov(j,l)*PP(i)+g_cov(i,l)*PP(j)-g_cov(i,j)*PP(l) ) 
        ENDDO 
      ENDDO
     ENDDO
    ENDDO
    DO l = 1, 3
     DO j = 1, 3
      DO i = 1, 3
          Gtilde(i) = Gtilde(i) + g_contr(j,l)*Christoffel_tilde(j,l,i) 
      ENDDO
     ENDDO     
    ENDDO   
    Z   = 0.5*matmul( g_cov, Ghat - Gtilde ) 
    Zup = matmul(phi2*g_contr, Z)
    !
    !print*,"SUMDD",sum(dDD) 
    !
    DO i = 1, 3 
     DO ip = 1, 3 
      DO m = 1, 3 
       DO k = 1, 3
          dChristoffelNCP(k,i,ip,m) = 0 
          dChristoffel_tildeNCP(k,i,ip,m) = 0 
          DO l = 1, 3 
            ! 
            dChristoffelNCP(k,i,ip,m) = dChristoffelNCP(k,i,ip,m) + g_contr(m,l)*( 0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)) )         & 
                                                                  - g_contr(m,l)*( g_cov(ip,l)*0.5*(dPP(k,i)+dPP(i,k))+g_cov(i,l)*0.5*(dPP(k,ip)+dPP(ip,k))-g_cov(i,ip)*0.5*(dPP(k,l)+dPP(l,k)) ) 
            ! 
            dChristoffel_tildeNCP(k,i,ip,m) = dChristoffel_tildeNCP(k,i,ip,m) + g_contr(m,l)*( 0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)) )           
            ! 
          ENDDO 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    ! 
    RiemannNCP = 0.0 
    DO i = 1, 3 
     DO ip = 1, 3 
      DO m = 1, 3 
       DO k = 1, 3
          RiemannNCP(i,k,ip,m) = dChristoffelNCP(k,i,ip,m)-dChristoffelNCP(ip,i,k,m)
       ENDDO
      ENDDO
     ENDDO
    ENDDO    
    ! 
    RicciNCP = 0.0 
    DO m = 1, 3 
     DO n = 1, 3
      DO l = 1, 3    
         RicciNCP(m,n) = RicciNCP(m,n) + RiemannNCP(m,l,n,l)  
      ENDDO
     ENDDO
    ENDDO
    !print*,"RicciNCP",RicciNCP
    !DO i = 1, 3 
     !DO k = 1, 3
     !print*,i,k,RicciNCP(i,k)
     !enddo
     !enddo
    !
    RNCP = phi2*SUM(g_contr*RicciNCP)
    !
    dGtildeNCP = 0.0
    !
    ! Here we directly compute the derivative of Gtilde from its original definition as contracted Christoffel symbol,
    ! without assuming unit determinant of the conformal metric. Back to the roots, and as few assumptions as possible...      
    !
    DO i = 1, 3 
     DO k = 1, 3
      DO j = 1, 3
       DO l = 1, 3
           dGtildeNCP(k,i) = dGtildeNCP(k,i) + g_contr(j,l)*dChristoffel_tildeNCP(k,j,l,i) 
       ENDDO
      ENDDO
     ENDDO
    ENDDO

    !DO i = 1, 3 
    !DO j = 1, 3 
    !print*,"dGhat",i,j,dGhat(i,j)
     !ENDDO
    !ENDDO
    !DO i = 1, 3 
    !DO j = 1, 3 
    !print*,"dGtildeNCP",i,j,dGtildeNCP(i,j)
     !ENDDO
    !ENDDO
    !DO i = 1, 3 
    !DO j = 1, 3 
    !print*,"g_cov",i,j,g_cov(i,j)
     !ENDDO
    !ENDDO
    !
    dZNCP = 0.0
    DO j = 1, 3 
     DO i = 1, 3
      DO k = 1, 3    
        dZNCP(k,i) = dZNCP(k,i) + ds*0.5*g_cov(i,j)*(dGhat(k,j)-dGtildeNCP(k,j))  
       ENDDO 
      ENDDO 
    ENDDO
    !print*,"dZNCP",dZNCP    
    !DO i = 1, 3 
    !DO j = 1, 3 
    !print*,"dZNCP",i,j,dZNCP(i,j)
     !ENDDO
    !ENDDO
    !
    DO j = 1, 3 
     DO i = 1, 3 
      nablaZNCP(i,j) = dZNCP(i,j)
     ENDDO
    ENDDO    
    !
    RicciPlusNablaZNCP = RicciNCP + ( nablaZNCP + transpose(nablaZNCP) ) 
    !DO i = 1, 3 
    !DO j = 1, 3 
    !print*,"RicciPlusNableZNCP",i,j,RicciPlusNablaZNCP(i,j)
     !ENDDO
    !ENDDO
    !print*,"RicciPlusNablaZNCP",RicciPlusNablaZNCP
    !
    RPlusTwoNablaZNCP = phi2*SUM(g_contr*RicciPlusNablaZNCP) 
    !
    nablaijalphaNCP = 0.0
    DO j = 1, 3 
     DO i = 1, 3 
       nablaijalphaNCP(i,j) = alpha*0.5*( dAA(i,j)+dAA(j,i) ) 
     ENDDO
    ENDDO 
    nablanablaalphaNCP = phi2*SUM( g_contr*nablaijalphaNCP )
    !!print*,"nablanablaalphaNCP",nablanablaalphaNCP
    !print*,"nablaijalphaNCP",nablaijalphaNCP
    !DO i = 1, 3 
    !DO j = 1, 3 
    !print*,"nablaijalphaNCP",i,j,nablaijalphaNCP(i,j)
     !ENDDO
    !ENDDO
    !
    SecondOrderTermsNCP = -nablaijalphaNCP + alpha*RicciPlusNablaZNCP 
    !print*,"SecondOrderTermsNCP",SecondOrderTermsNCP
    traceNCP = SUM( g_contr*SecondOrderTermsNCP )
    !print*,"traceNCP",traceNCP
    SecondOrderTermsNCP = SecondOrderTermsNCP - 1./3.*g_cov*traceNCP 
    !
    ! Now assemble all this terrible stuff... 
    !
    dtgamma = 0.0 
    !
    ! Main variables of the CCZ4 system 
    dtK = phi2*SecondOrderTermsNCP + beta(1)*dAex(1,:,:) + beta(2)*dAex(2,:,:) + beta(3)*dAex(3,:,:)      ! extrinsic curvature
    !
    dtTraceK = -nablanablaalphaNCP + alpha*RPlusTwoNablaZNCP + SUM(beta(:)*dtraceK(:)) 
    !!print*,"dtTraceK",dtTraceK
    !
    traceB = BB(1,1) + BB(2,2) + BB(3,3) 
    !print*,"traceB",traceB
    dtphi   = 0.0 
    dtalpha = 0.0 

    !print*,"dTheta",dTheta
    !print*,"beta",beta
    !print*,"alpha",alpha,"e",e

    Aupdown = SUM(Aex*Aup) 
    dtTheta = 0.5*alpha*e**2*( RplusTwoNablaZNCP ) + beta(1)*dTheta(1) + beta(2)*dTheta(2) + beta(3)*dTheta(3)        ! *** original cleaning *** 
    !print*,"dtTheta",dtTheta
    !
    divAupNCP = 0.0
    DO i = 1, 3
        DO j = 1, 3
         DO l = 1, 3
          DO k = 1, 3    
            divAupNCP(i) = divAupNCP(i) + g_contr(i,l)*g_contr(j,k)*dAex(j,l,k) 
          ENDDO
         ENDDO
        ENDDO        
    ENDDO     
    DO i = 1, 3 
        Mom(i) = - 2./3.*SUM(g_contr(i,:)*dtraceK(:)) + divAupNCP(i)  
    ENDDO 
    !
    !
    DO i = 1, 3
        dtGhat(i) = - 4./3.*alpha*SUM(g_contr(i,:)*dtraceK(:))     &  
                    + 2.0*alpha*SUM( g_contr(:,i)*( dTheta(:)  ) ) &                    
                    + beta(1)*dGhat(1,i) + beta(2)*dGhat(2,i) + beta(3)*dGhat(3,i) 
        DO l = 1, 3
         DO k = 1, 3
             dtGhat(i) = dtGhat(i) + g_contr(k,l)*0.5*(dBB(k,l,i)+dBB(l,k,i)) + 1./3*g_contr(i,k)*0.5*(dBB(k,l,l)+dBB(l,k,l)) 
         ENDDO
        ENDDO         
    ENDDO
    
    DO k = 1, 3 
        ov(k) = 2*alpha*( SUM(g_contr(:,:)*dAex(k,:,:)) )           ! here we can use the constraint that trace A tilde = 0. 
    ENDDO
    !print*,"ov",ov
    dtGhat = dtGhat + sk*matmul(g_contr,ov)                         ! Ghat is an "up" vector, so we need to multiply with g_contr 
    !print*,"dtGhat",dtGhat
    !
    dtbb = xi*dtGhat + bs*( beta(1)*gradQ(21:23,1) + beta(2)*gradQ(21:23,2) + beta(3)*gradQ(21:23,3) - beta(1)*gradQ(14:16,1) - beta(2)*gradQ(14:16,2) - beta(3)*gradQ(14:16,3) ) 
    dtbb = sk*dtbb  
    !
    dtbeta  = 0.0    
    !
    ! Auxiliary variables
    !print*,"dtK",dtraceK,"dK0",dK0,"dTh",dTheta
    dtA = -alpha*fa*( dtraceK(:) -dK0(:) - c*2*dTheta(:) ) + beta(1)*dAA(1,:) + beta(2)*dAA(2,:) + beta(3)*dAA(3,:)
    !print*,"AAAA",dtA
    DO k = 1, 3 
        dtA(k) = dtA(k) - sk*alpha*fa*( SUM(g_contr(:,:)*dAex(k,:,:)) )   ! here we can use the constraint that trace A tilde = 0. 
    ENDDO
    !print*,"AAAA'",dtA
    !
    ! We have removed the conservative fluxes for CCZ4, so put all the stuff into the NCP and FusedNCP 
    dtB(:,1) = fff*gradQ(21,:)  
    dtB(:,2) = fff*gradQ(22,:)  
    dtB(:,3) = fff*gradQ(23,:)  
    ! #ordB1#     
    ! for the ordering constraints, we have to check whether they should be multiplied by alpha**2, or not... 
    DO i = 1, 3 
     DO k = 1, 3 
      DO j = 1, 3 
         dtB(k,i) = dtB(k,i) + mu*alpha**2*g_contr(i,j)*( (dPP(k,j)-dPP(j,k)) )  
         DO n = 1, 3 
          DO l = 1, 3 
            dtB(k,i) = dtB(k,i) - mu*alpha**2*g_contr(i,j)*g_contr(n,l)*( dDD(k,l,j,n)-dDD(l,k,j,n) )   
          ENDDO 
         ENDDO
      ENDDO 
      !
      ENDDO 
    ENDDO 
    dtB = dtB + bs*( beta(1)*dBB(1,:,:) + beta(2)*dBB(2,:,:) + beta(3)*dBB(3,:,:) ) 
    ! New stuff 1, which makes appear a Lie derivative and a second order ordering constraint 
    dtB = dtB*sk 
    !
    dtD = -alpha*dAex
    DO i = 1, 3
      DO j = 1, 3
       DO k = 1, 3 
        DO m = 1, 3
            !dtD(k,i,j) = dtD(k,i,j) + ( 0.5*(g_cov(m,i)*0.5*(dBB(k,j,m)+dBB(j,k,m))+g_cov(m,j)*0.5*(dBB(k,i,m)+dBB(i,k,m)) ) - 1./3.*g_cov(i,j)*0.5*(dBB(k,m,m)+dBB(m,k,m)) ) 
            dtD(k,i,j) = dtD(k,i,j) + 0.25*( g_cov(m,i)*(dBB(k,j,m)+dBB(j,k,m)) + g_cov(m,j)*(dBB(k,i,m)+dBB(i,k,m)) ) - 1./6.*g_cov(i,j)*(dBB(k,m,m)+dBB(m,k,m)) 
            DO n = 1, 3
                dtD(k,i,j) = dtD(k,i,j) + 1./3*alpha*g_cov(i,j)*g_contr(n,m)*dAex(k,n,m)     ! explicitly remove the trace of tilde A again 
            ENDDO 
        ENDDO        
       ENDDO
      ENDDO
    ENDDO 


    dtD = dtD + beta(1)*dDD(1,:,:,:) + beta(2)*dDD(2,:,:,:) + beta(3)*dDD(3,:,:,:)
    !DO i = 1, 3
      !DO j = 1, 3
       !DO k = 1, 3 
       !print*,i,j,k,dtD(i,j,k)
       !enddo
       !enddo
       !enddo
    !
    dtP = + beta(1)*dPP(1,:) + beta(2)*dPP(2,:) + beta(3)*dPP(3,:)    
    DO k = 1, 3 
     dtP(k) = dtP(k) + 1./3.*alpha*dtraceK(k) + sk*1./3.*alpha*SUM(g_contr(:,:)*dAex(k,:,:))  ! use the fact that trace A tilde = 0 
     DO i = 1, 3 
          dtP(k) = dtP(k) - 1./3.*0.5*(dBB(k,i,i)+dBB(i,k,i))  
     ENDDO
    ENDDO 
    !
    BgradQ(1:6)    = (/ dtgamma(1,1), dtgamma(1,2), dtgamma(1,3), dtgamma(2,2), dtgamma(2,3), dtgamma(3,3) /)          ! \tilde \gamma_ij 
    BgradQ(7:12)   = (/ dtK(1,1), dtK(1,2), dtK(1,3), dtK(2,2), dtK(2,3), dtK(3,3) /)                                  ! \tilde A_ij 
    BgradQ(13)     = dtTheta                                                                                           ! Theta       
    BgradQ(14:16)  = dtGhat(1:3)                                                                                       ! \hat \Gamma^i           
    BgradQ(17)     = dtalpha                                                                                           ! log alpha 
    BgradQ(18:20)  = dtbeta                                                                                            ! beta^i 
    BgradQ(21:23)  = dtbb                                                                                              ! b^i 
    BgradQ(24:26)  = dtA(1:3)                                                                                          ! A_k       
    BgradQ(27:35)  = (/ dtB(1,1), dtB(2,1), dtB(3,1), dtB(1,2), dtB(2,2), dtB(3,2), dtB(1,3), dtB(2,3), dtB(3,3) /)    ! B_k^i 
    BgradQ(36:41)  = (/ dtD(1,1,1), dtD(1,1,2), dtD(1,1,3), dtD(1,2,2), dtD(1,2,3), dtD(1,3,3) /)                      ! D_kij 
    BgradQ(42:47)  = (/ dtD(2,1,1), dtD(2,1,2), dtD(2,1,3), dtD(2,2,2), dtD(2,2,3), dtD(2,3,3) /)                      ! D_kij 
    BgradQ(48:53)  = (/ dtD(3,1,1), dtD(3,1,2), dtD(3,1,3), dtD(3,2,2), dtD(3,2,3), dtD(3,3,3) /)                      ! D_kij 
    BgradQ(54)     = dtTraceK                                                                                          ! traceK 
    BgradQ(55)     = dtphi                                                                                             ! log phi 
    BgradQ(56:58)  = dtP                                                                                               ! P_k 
    BgradQ(59:nVar) = 0.0  
    !
    BgradQ = -BgradQ ! change sign, since we work on the left hand side in PDENCP 
#endif
END subroutine PDENCP

program main
  implicit none
  REAL :: BgradQ(59), gradQin(59,3), Q(59), dummy(4,2)
  integer :: i,j,k,l,m
  BgradQ = 1

  gradQin =0
  gradQin(:,1)=1

  ! The initial data for Q is chosen such that we do not divide by 0 in the inverse determinant
  Q=1
  Q(1)=2
  Q(2)=3

  !do i=1,1000000
    call PDENCP(BgradQ, Q, gradQin)
  !enddo
  do i=1,59
    print*,i,BgradQ(i)
  enddo
  !call PDENCP(BgradQ, Q, gradQin)
  !print*, BgradQ
end program main
