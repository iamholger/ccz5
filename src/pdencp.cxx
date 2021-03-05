#include <cmath>
#include <iostream>



#pragma omp declare target
void pdencp_(double* BgradQ, const double* const Q, const double* const gradQSerialised, const int normal)
{
#if defined(Dim3)
    constexpr int nDim = 3;
#elif defined(Dim2)
    constexpr int nDim = 2;
#endif
    constexpr int nVar(59), nParam(0), d(3);

    // First model parameter ds here (ds == CCZ4ds)
    const double ds = 1.0; // NOTE this param seems to always be 1
    // Second model parameter CCZ4e
    const double e = 1.0; 
    const double e2 = e*e;
    // Parameter CCZ4sk here
    // TODO test for values that from 0 to check logic
    const double sk=0.1;
    // Input and output variables
    // Params CCZ4xi, CCZ4bs
    const double xi=0.0;
    const double bs=0.0;
    // CCZ4 parameeter -- harmonic lapse CCZ4c
    const double c   = 1.0;
    // Param CCZ4fff, CCZ4mu
    const double fff = 0;
    const double mu  = 0.2;
    const double fa  = 1.0;
    const double faa = 1.0;
    double gradQin[59][3] ={0};

    // De-serialise input data and fill static array
    // FIXME the use of 2D arrays can be avoided: all terms not in the normal are 0
    for (int i=0; i<nVar; i++) gradQin[i][normal] = gradQSerialised[i+normal*nVar];


    // Note g_cov is symmetric
    const double g_cov[3][3] = { {Q[0], Q[1], Q[2]}, {Q[1], Q[3], Q[4]}, {Q[2], Q[4], Q[5]} };
    const double invdet = 1./( Q[0]*Q[3]*Q[5] - Q[0]*Q[4]*Q[4] - Q[1]*Q[1]*Q[5] + 2*Q[1]*Q[2]*Q[4] -Q[2]*Q[2]*Q[3]);
    //std::cout << "INVDET " << invdet << "\n";

    const double g_contr[3][3] = {
        { ( Q[3]*Q[5]-Q[4]*Q[4])*invdet, -( Q[1]*Q[5]-Q[2]*Q[4])*invdet, -(-Q[1]*Q[4]+Q[2]*Q[3])*invdet},
        {-( Q[1]*Q[5]-Q[4]*Q[2])*invdet,  ( Q[0]*Q[5]-Q[2]*Q[2])*invdet, -( Q[0]*Q[4]-Q[2]*Q[1])*invdet},
        {-(-Q[1]*Q[4]+Q[3]*Q[2])*invdet, -( Q[0]*Q[4]-Q[1]*Q[2])*invdet,  ( Q[0]*Q[3]-Q[1]*Q[1])*invdet}
    };
    //std::cout << "g_contr\n";
    //for (int i=0;i<3;i++)
    //{
      //for (int j=0;j<3;j++) std::cout << " " << g_contr[i][j];
      //std::cout << "\n";
    //}
    //std::cout << "\n";


    // NOTE Aex is symmetric
    double Aex[3][3] = { {Q[6], Q[7], Q[8]}, {Q[7], Q[9], Q[10]}, {Q[8], Q[10], Q[11]} };
    //std::cout << "Aex\n";
    //for (int i=0;i<3;i++)
    //{
      //for (int j=0;j<3;j++) std::cout << " " << Aex[i][j];
      //std::cout << "\n";
    //}
    //std::cout << "\n";

    double traceA = 0;
    for (int i=0;i<3;i++)
    for (int j=0;j<3;j++) traceA+=g_contr[i][j]*Aex[i][j];
    traceA *= 1./3;
    //std::cout << "traceA " <<3*traceA << "\n";
    
    for (int i=0;i<3;i++)
    for (int j=0;j<3;j++) Aex[i][j] -= traceA * g_cov[i][j];

    // Matrix multiplications Amix = matmul(g_contr, Aex) Aup  = matmul(g_contr, mytranspose(Amix))
    double Amix[3][3]={0};
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    for (int u = 0; u < 3; u++) Amix[i][j] += g_contr[i][u] * Aex[u][j];
    double Aup[3][3]={0};
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    for (int u = 0; u < 3; u++) Aup[i][j] += g_contr[i][u] * Amix[j][u]; // Note the transposition is in the indices

    const double DD[3][3][3] = { 
        {{Q[35], Q[36], Q[37]}, {Q[36], Q[38], Q[39]}, {Q[37], Q[39], Q[40]}},
        {{Q[41], Q[42], Q[43]}, {Q[42], Q[44], Q[45]}, {Q[43], Q[45], Q[46]}},
        {{Q[47], Q[48], Q[49]}, {Q[48], Q[50], Q[51]}, {Q[49], Q[51], Q[52]}}
    };
    double dgup[3][3][3] {0};
    for (int k = 0; k < 3; k++)
    for (int m = 0; m < 3; m++)
    for (int l = 0; l < 3; l++) 
    for (int n = 0; n < 3; n++) 
    for (int j = 0; j < 3; j++) dgup[k][m][l] -= g_contr[m][n]*g_contr[j][l]*2*DD[k][n][j];

    const double PP[3] = {Q[55], Q[56], Q[57]};

    double Christoffel[3][3][3]       = {0};
    double Christoffel_tilde[3][3][3] = {0};
    double Christoffel_kind1[3][3][3] = {0};
    
    for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
    for (int k = 0; k < 3; k++)
    {
        Christoffel_kind1[i][j][k] = DD[k][i][j] + DD[j][i][k] - DD[i][j][k];
        for (int l = 0; l < 3; l++)
        {
            Christoffel_tilde[i][j][k] += g_contr[k][l] * ( DD[i][j][l] + DD[j][i][l] - DD[l][i][j] );
            Christoffel[i][j][k]       += g_contr[k][l] * ( DD[i][j][l] + DD[j][i][l] - DD[l][i][j] ) - g_contr[k][l] * ( g_cov[j][l] * PP[i] + g_cov[i][l] * PP[j] - g_cov[i][j] * PP[l] );
        }
    }

    double Gtilde[3] = {0};
    for (int l = 0; l < 3; l++) 
    for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++) Gtilde[i] += g_contr[j][l] * Christoffel_tilde[j][l][i];

    const double Ghat[3] = {Q[13], Q[14], Q[15]};

    const double phi = std::exp(std::fmax(-20., std::fmin(20.,Q[54])));
    const double phi2 = phi*phi;

    //std::cout << "phi2 " << phi2 << "\n";

    double Z[3] = {0};
    for (int i=0;i<3;i++)
    for (int j=0;j<3;j++) Z[i] += ( g_cov[i][j]* (Ghat[j] - Gtilde[j]));
    double Zup[3] = {0};
    for (int i=0;i<3;i++)
    for (int j=0;j<3;j++) Zup[i] += phi2 * g_contr[i][j] * Z[j];


    const double dDD[3][3][3][3] = { 
        {
            {
                {gradQin[35][0],gradQin[36][0],gradQin[37][0]}, {gradQin[36][0],gradQin[38][0],gradQin[39][0]}, {gradQin[37][0],gradQin[39][0],gradQin[40][0]},
            },
            {
                {gradQin[41][0],gradQin[42][0],gradQin[43][0]}, {gradQin[42][0],gradQin[44][0],gradQin[45][0]}, {gradQin[43][0],gradQin[45][0],gradQin[46][0]},
            },
            {
                {gradQin[47][0],gradQin[48][0],gradQin[49][0]}, {gradQin[48][0],gradQin[50][0],gradQin[51][0]}, {gradQin[49][0],gradQin[51][0],gradQin[52][0]}
            }
        },
        {
            {
                {gradQin[35][1],gradQin[36][1],gradQin[37][1]}, {gradQin[36][1],gradQin[38][1],gradQin[39][1]}, {gradQin[37][1],gradQin[39][1],gradQin[40][1]},
            },
            {
                {gradQin[41][1],gradQin[42][1],gradQin[43][1]}, {gradQin[42][1],gradQin[44][1],gradQin[45][1]}, {gradQin[43][1],gradQin[45][1],gradQin[46][1]},
            },
            {
                {gradQin[47][1],gradQin[48][1],gradQin[49][1]}, {gradQin[48][1],gradQin[50][1],gradQin[51][1]}, {gradQin[49][1],gradQin[51][1],gradQin[52][1]}
            }
        },
        {
            {
                {gradQin[35][2],gradQin[36][2],gradQin[37][2]}, {gradQin[36][2],gradQin[38][2],gradQin[39][2]}, {gradQin[37][2],gradQin[39][2],gradQin[40][2]},
            },
            {
                {gradQin[41][2],gradQin[42][2],gradQin[43][2]}, {gradQin[42][2],gradQin[44][2],gradQin[45][2]}, {gradQin[43][2],gradQin[45][2],gradQin[46][2]},
            },
            {
                {gradQin[47][2],gradQin[48][2],gradQin[49][2]}, {gradQin[48][2],gradQin[50][2],gradQin[51][2]}, {gradQin[49][2],gradQin[51][2],gradQin[52][2]}
            }
        }
    };


    const double dPP[3][3] = {
        {gradQin[55][0],gradQin[56][0],gradQin[57][0]},
        {gradQin[55][1],gradQin[56][1],gradQin[57][1]},
        {gradQin[55][2],gradQin[56][2],gradQin[57][2]}
    };



    double dChristoffelNCP[3][3][3][3] = {0};
    double dChristoffel_tildeNCP[3][3][3][3] = {0};
    for (int i = 0; i < 3; i++)
    for (int ip = 0; ip < 3; ip++)
    for (int m = 0; m < 3; m++)
    for (int k = 0; k < 3; k++)
    {
        dChristoffelNCP[k][i][ip][m] = 0;
        dChristoffel_tildeNCP[k][i][ip][m] = 0;
        for (int l = 0; l < 3; l++)
        {
            dChristoffelNCP[k][i][ip][m] +=  0.5*g_contr[m][l] * ( 
                    dDD[k][i][ip][l] + dDD[i][k][ip][l] + dDD[k][ip][i][l] + dDD[ip][k][i][l] - dDD[k][l][i][ip] + dDD[l][k][i][ip] 
                    - g_cov[ip][l]*(dPP[k][i] + dPP[i][k]) - g_cov[i][l]*(dPP[k][ip]+dPP[ip][k]) +  g_cov[i][ip]*(dPP[k][l]+dPP[l][k]) );



            dChristoffel_tildeNCP[k][i][ip][m] += 0.5*g_contr[m][l]*(dDD[k][i][ip][l] + dDD[i][k][ip][l] + dDD[k][ip][i][l] + dDD[ip][k][i][l] - dDD[k][l][i][ip] - dDD[l][k][i][ip]);
            
        }
    }


    double RiemannNCP[3][3][3][3] = {0};
    for (int i = 0; i < 3; i++)
    for (int ip = 0; ip < 3; ip++)
    for (int m = 0; m < 3; m++)
    for (int k = 0; k < 3; k++) RiemannNCP[i][k][ip][m] = dChristoffelNCP[k][i][ip][m] - dChristoffelNCP[ip][i][k][m];

    double RicciNCP[3][3] = {0};
    for (int m = 0; m < 3; m++)
    for (int n = 0; n < 3; n++)
    for (int l = 0; l < 3; l++) RicciNCP[m][n] += RiemannNCP[m][l][n][l];  

    double dGtildeNCP[3][3] = {0};
    for (int i = 0; i < 3; i++)
    for (int k = 0; k < 3; k++)
    for (int j = 0; j < 3; j++)
    for (int l = 0; l < 3; l++) dGtildeNCP[k][i] += g_contr[j][l]*dChristoffel_tildeNCP[k][j][l][i];

    //std::cout << "dGtildeNCP\n"; // Checked!
    //for (int i=0;i<3;i++)
    //{
      //for (int j=0;j<3;j++) std::cout << i+1 << " " << j+1 << " " << dGtildeNCP[i][j] << "\n";
      ////std::cout << "\n";
    //}
    //std::cout << "g_cov\n"; // Checked!
    //for (int i=0;i<3;i++)
    //{
      //for (int j=0;j<3;j++) std::cout << i+1 << " " << j+1 << " " << g_cov[i][j] << "\n";
      ////std::cout << "\n";
    //}
    
    const double dGhat[3][3] = {
        {gradQin[13][0],gradQin[14][0],gradQin[15][0]},
        {gradQin[13][1],gradQin[14][1],gradQin[15][1]},
        {gradQin[13][2],gradQin[14][2],gradQin[15][2]}
    };
    //std::cout << "dGhat\n"; // Checked!
    //for (int i=0;i<3;i++)
    //{
      //for (int j=0;j<3;j++) std::cout << i+1 << " " << j+1 << " " << dGhat[i][j] << "\n";
      ////std::cout << "\n";
    //}

    //std::cout << "\n";
    //std::cout << "g_cov\n"; // Checked!
    //for (int i=0;i<3;i++)
    //{
      //for (int j=0;j<3;j++) std::cout << " " << g_cov[i][j];
      //std::cout << "\n";
    //}
    //std::cout << "\n";

    //std::cout << "RicciNCP\n"; // Checked!
    //for (int i=0;i<3;i++)
    //{
      //for (int j=0;j<3;j++) std::cout << i+1 << " " << j+1 << " " << RicciNCP[i][j] << "\n";
    //}

    //BUG?
    double dZNCP[3][3] = {0};
    for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
    for (int k = 0; k < 3; k++) dZNCP[k][i] += ds*0.5*g_cov[i][j]*(dGhat[k][j]-dGtildeNCP[k][j]);  
    
    //std::cout << "dZNCP\n";
    //for (int i=0;i<3;i++)
    //{
      //for (int j=0;j<3;j++) std::cout << i+1 << " " << j+1 << " " << dZNCP[i][j] << "\n";
    //}

    // !BUG!
    double RicciPlusNablaZNCP[3][3];
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) RicciPlusNablaZNCP[i][j] = RicciNCP[i][j] + dZNCP[i][j] + dZNCP[j][i];
    //std::cout << "RicciPlusNablaZNCP\n"; // Checked!
    //for (int i=0;i<3;i++)
    //{
      //for (int j=0;j<3;j++) std::cout << " " << i+1 << " " << j+1 << " " << RicciPlusNablaZNCP[i][j] <<"\n";
    //}

    double RPlusTwoNablaZNCP = 0; 
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) RPlusTwoNablaZNCP += g_contr[i][j]*RicciPlusNablaZNCP[i][j]; // TODO fuse these steps
    RPlusTwoNablaZNCP*=phi2;

    const double alpha = std::exp(std::fmax(-20., std::fmin(20.,Q[16])));
    const double alpha2 = alpha*alpha;


    const double AA[3] = {Q[23], Q[24], Q[25]};
    const double dAA[3][3] = {
        {gradQin[23][0],gradQin[24][0],gradQin[25][0]},
        {gradQin[23][1],gradQin[24][1],gradQin[25][1]},
        {gradQin[23][2],gradQin[24][2],gradQin[25][2]}
    };
    //std::cout << "dAA\n";
    //for (int i=0;i<3;i++)
    //{
      //for (int j=0;j<3;j++) std::cout << " " << dAA[i][j];
      //std::cout << "\n";
    //}
    //std::cout << "\n";

    double nablaijalphaNCP[3][3];
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) nablaijalphaNCP[i][j] = alpha*0.5*( dAA[i][j] + dAA[j][i] ); 
    
    double nablanablaalphaNCP = 0;
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) nablanablaalphaNCP += g_contr[i][j]*nablaijalphaNCP[i][j]; 
    nablanablaalphaNCP*=phi2;
    //std::cout << "nablanablaalphaNCP " << nablanablaalphaNCP << "\n";
   
    double SecondOrderTermsNCP[3][3];
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)  SecondOrderTermsNCP[i][j] = -nablaijalphaNCP[i][j] + alpha*RicciPlusNablaZNCP[i][j];
    
    //std::cout << "nablaijalphaNCP\n";
    //for (int i=0;i<3;i++)
    //{
      //for (int j=0;j<3;j++) std::cout << i+1 << " " <<j+1<<" " << nablaijalphaNCP[i][j] << "\n";
    //}
    
    double traceNCP = 0;
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) traceNCP += g_contr[i][j]*SecondOrderTermsNCP[i][j];
    //std::cout << "traceNCP " << traceNCP << "\n";
    
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)  SecondOrderTermsNCP[i][j] -= 1./3 * traceNCP * g_cov[i][j];
        

    const double beta[3] = {Q[17], Q[18], Q[19]};
    //std::cerr << "beta " << beta[0] << " " << beta[1] << " " << beta[2] << "\n"; 


    // checked!
    const double dAex[3][3][3] = { 
        {{gradQin[6][0],gradQin[7][0],gradQin[8][0]}, {gradQin[7][0], gradQin[9][0], gradQin[10][0]},  {gradQin[8][0], gradQin[10][0], gradQin[11][0]}},
        {{gradQin[6][1],gradQin[7][1],gradQin[8][1]}, {gradQin[7][1], gradQin[9][1], gradQin[10][1]},  {gradQin[8][1], gradQin[10][1], gradQin[11][1]}},
        {{gradQin[6][2],gradQin[7][2],gradQin[8][2]}, {gradQin[7][2], gradQin[9][2], gradQin[10][2]},  {gradQin[8][2], gradQin[10][2], gradQin[11][2]}}
    };
    //! Now assemble all this terrible stuff... 
    //!
    //! Main variables of the CCZ4 system 
    double dtK[3][3];
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)  dtK[i][j] = phi2*SecondOrderTermsNCP[i][j] + beta[0] * dAex[0][i][j] + beta[1] * dAex[1][i][j] + beta[2] * dAex[2][i][j]; // extrinsic curvature
   
    const double dtraceK[3] = {gradQin[53][0], gradQin[53][1], gradQin[53][1]};

    double dtTraceK = -nablanablaalphaNCP + alpha*RPlusTwoNablaZNCP + beta[0]*dtraceK[0] + beta[1]*dtraceK[1] + beta[2]*dtraceK[2];

    //std::cout << "dtTraceK: " << dtTraceK << "\n";



    const double BB[3][3] = {
        {Q[26], Q[27], Q[28]}, {Q[29], Q[30], Q[31]}, {Q[32], Q[33], Q[34]}
    };

    double traceB = BB[0][0] + BB[1][1] + BB[2][2]; // TODO direct from Q! 

    //std::cout << "traceB " << traceB << "\n";

    double Aupdown = 0;
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) Aupdown += Aex[i][j]*Aup[i][j];  

  
    const double dTheta[3] = {gradQin[13][0],gradQin[13][1],gradQin[13][2]};
    //std::cout << "dTheta: " << dTheta[0] << " " << dTheta[1] << " " << dTheta[2] << "\n"; 
    //std::cout << "beta: " << beta[0] << " " << beta[1] << " " << beta[2] << "\n"; 
    const double dtTheta = 0.5*alpha*e2*( RPlusTwoNablaZNCP ) + beta[0]*dTheta[0] + beta[1]*dTheta[1] + beta[2]*dTheta[2]; // *** original cleaning *** 
  
    //std::cout << "dtTheta: " << dtTheta << "\n"; 

    double divAupNCP[3] = {0};

    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    for (int l = 0; l < 3; l++)
    for (int k = 0; k < 3; k++) divAupNCP[i] += g_contr[i][l]*g_contr[j][k]*dAex[j][l][k];


    //double Mom[3];
    //for (int i = 0; i < 3; i++)
    //{
        //Mom[i] = divAupNCP[i];
        //double temp=0;
        //for (int j = 0; j < 3; j++) temp +=g_contr[i][j]*dtraceK[j];
        //Mom[i] -= 2./3*temp;
    //}




    const double dBB[3][3][3] = {
        {
            {sk*gradQin[26][0],sk*gradQin[27][0],sk*gradQin[28][0]}, {sk*gradQin[29][0],sk*gradQin[30][0],sk*gradQin[31][0]}, {sk*gradQin[32][0],sk*gradQin[33][0],sk*gradQin[34][0]}
        },
        {
            {sk*gradQin[26][1],sk*gradQin[27][1],sk*gradQin[28][1]}, {sk*gradQin[29][1],sk*gradQin[30][1],sk*gradQin[31][1]}, {sk*gradQin[32][1],sk*gradQin[33][1],sk*gradQin[34][1]}
        },
        {
            {sk*gradQin[26][2],sk*gradQin[27][2],sk*gradQin[28][2]}, {sk*gradQin[29][2],sk*gradQin[30][2],sk*gradQin[31][2]}, {sk*gradQin[32][2],sk*gradQin[33][2],sk*gradQin[34][2]}
        }
    };

    //DO i = 1, 3
        //dtGhat(i) = - 4./3.*alpha*SUM(g_contr(i,:)*dtraceK(:))     &  
                    //+ 2.0*alpha*SUM( g_contr(:,i)*( dTheta(:)  ) ) &                    
                    //+ beta(1)*dGhat(1,i) + beta(2)*dGhat(2,i) + beta(3)*dGhat(3,i) 
        //DO l = 1, 3
         //DO k = 1, 3
             //dtGhat(i) = dtGhat(i) + g_contr(k,l)*0.5*(dBB(k,l,i)+dBB(l,k,i)) + 1./3*g_contr(i,k)*0.5*(dBB(k,l,l)+dBB(l,k,l)) 
         //ENDDO
        //ENDDO         
    //ENDDO

    double dtGhat[3];
    for (int i = 0; i < 3; i++)
    {
        double temp=0, temp2=0;
        for (int j = 0; j < 3; j++)
        {
            temp  +=g_contr[i][j]*dtraceK[j];
            temp2 +=g_contr[j][i]*dTheta[j];
        }
        dtGhat[i] = -4./3.*alpha*temp + 2*alpha*temp2 + beta[0]*dGhat[0][i] + beta[1]*dGhat[1][i] + beta[2]*dGhat[2][i];
        for (int l = 0; l < 3; l++)
        for (int k = 0; k < 3; k++) dtGhat[i] += g_contr[k][l]*0.5*(dBB[k][l][i] + dBB[l][k][i]) + 1./3.*g_contr[i][k]*0.5*(dBB[k][l][l] + dBB[l][k][l]);
    }


    double ov[3];
    for (int k = 0; k < 3; k++)
    {
        double temp=0;
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) temp += g_contr[i][j]*dAex[k][i][j];
        ov[k] = 2*alpha*temp;
    }
    //std::cout << "ov: " << ov[0] << " " <<ov[1] << " " << ov[2]<<"\n";

    // TODO check this for sk non-zero
    // matrix vector multiplication in a loop and add result to existing vector
    // checked for sk=0
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) dtGhat[i] += sk*g_contr[i][j]*ov[j];
    
    //for (int i = 0; i < 3; i++) std::cout << "dtGhat " << i+1 << " " << dtGhat[i] << "\n";

    double dtbb[3];
    for (int i = 0; i < 3; i++)
    {
        dtbb[i] = xi*dtGhat[i] + bs * ( beta[0]*gradQin[20+i][0] + beta[1]*gradQin[20+i][1] + beta[2]*gradQin[20+i][2] - beta[0]*gradQin[13+i][0] - beta[1]*gradQin[13+i][1] - beta[2]*gradQin[13+i][2]);
        dtbb[i] *= sk;

    }


    // Auxiliary variables 
    double dtA[3];
    double dK0[3] = {0}; // FIXME we just add 0 all the time???
    for (int i = 0; i < 3; i++)
    {
        dtA[i] = -alpha*fa*(dtraceK[i] - dK0[i] - c*2*dTheta[i]) + beta[0]*dAA[0][i] + beta[1]*dAA[1][i] + beta[2]*dAA[2][i];
    }
    //for (int i = 0; i < 3; i++) std::cout << "dtA " << i+1 << " " << dtA[i] << "\n";

    //dtA = -alpha*fa*( dtraceK(:) -dK0(:) - c*2*dTheta(:) ) + beta(1)*dAA(1,:) + beta(2)*dAA(2,:) + beta(3)*dAA(3,:)

    // TODO check correctness for non zero params (CCZ4%sk)
    for (int k = 0; k < 3; k++)
    {
        double temp=0;
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) temp += g_contr[i][j]*dAex[k][i][j]; // TODO we computed this quantity  few lines earlier alrady
        dtA[k] -= sk*alpha*fa*temp;
    }

    double dtB[3][3] = {
        {fff*gradQin[20][0],fff*gradQin[20][1],fff*gradQin[20][2]},
        {fff*gradQin[21][0],fff*gradQin[21][1],fff*gradQin[21][2]},
        {fff*gradQin[22][0],fff*gradQin[22][1],fff*gradQin[22][2]}
    };


    for (int i = 0; i < 3; i++)
    for (int k = 0; k < 3; k++)
    for (int j = 0; j < 3; j++)
    {
        dtB[k][i] += mu*alpha2 * g_contr[i][j]*( dPP[k][j] - dPP[j][k]);
        for (int n = 0; n < 3; n++)
        for (int l = 0; l < 3; l++) dtB[k][i] -= mu*alpha2 * g_contr[i][j]*g_contr[n][l]*( dDD[k][l][j][n] - dDD[l][k][j][n]);
    }

    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) dtB[i][j] += bs*(beta[0]*dBB[0][i][j] + beta[1]*dBB[1][i][j] + beta[2]*dBB[2][i][j]);

    // NOTE 0 value param
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) dtB[i][j] *= sk;

    double dtD[3][3][3];
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    for (int k = 0; k < 3; k++) dtD[i][j][k] = -alpha*dAex[i][j][k];
    
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    for (int k = 0; k < 3; k++)
    for (int m = 0; m < 3; m++)
    {
        dtD[k][i][j] += ( 0.25*(g_cov[m][i]*(dBB[k][j][m] + dBB[j][k][m]) + g_cov[m][j]*(dBB[k][i][m]+dBB[i][k][m])) - 1./6.*g_cov[i][j]*(dBB[k][m][m]+dBB[m][k][m]) );
        for (int n = 0; n < 3; n++)
            dtD[k][i][j] += 1./3*alpha*g_cov[i][j]*g_contr[n][m]*dAex[k][n][m]; // explicitly remove the trace of tilde A again  
    }

    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    for (int k = 0; k < 3; k++) dtD[i][j][k] += beta[0]*dDD[0][i][j][k] + beta[1]*dDD[1][i][j][k] + beta[2]*dDD[2][i][j][k];
   
    // Checked!

/*    for (int i = 0; i < 3; i++)*/
    //for (int j = 0; j < 3; j++)
    //for (int k = 0; k < 3; k++) 
    //{
      //std::cout << "dtD " << i+1 << " " << j+1 << " "  << k+1 << " " << dtD[i][j][k] << "\n";
    //}


    double dtP[3];
    for (int i = 0; i < 3; i++) dtP[i] = beta[0]*dPP[0][i] + beta[1]*dPP[1][i] + beta[2]*dPP[2][i];

    // NOTE test for non zero params
    for (int k = 0; k < 3; k++)
    {
        double temp=0;
        for (int m = 0; m < 3; m++)
        for (int n = 0; n < 3; n++) temp += g_contr[m][n]*dAex[k][m][n]; // TODO we computed this quantity  few lines earlier alrady
        dtP[k] += 1./3*alpha*(dtraceK[k] + sk*temp);
        for (int i = 0; i < 3; i++)
            dtP[k] -= 1./6*(dBB[k][i][i] + dBB[i][k][i]);

    }

    double dtgamma[3][3] = {0};
    double dtalpha = 0;
    double dtbeta[3] = {0};
    double dtphi = 0;

    BgradQ[0]  = -dtgamma[0][0];
    BgradQ[1]  = -dtgamma[0][1];
    BgradQ[2]  = -dtgamma[0][2];
    BgradQ[3]  = -dtgamma[1][1];
    BgradQ[4]  = -dtgamma[1][2];
    BgradQ[5]  = -dtgamma[2][2];
    BgradQ[6]  = -dtK[0][0];
    BgradQ[7]  = -dtK[0][1];
    BgradQ[8]  = -dtK[0][2];
    BgradQ[9]  = -dtK[1][1];
    BgradQ[10] = -dtK[1][2];
    BgradQ[11] = -dtK[2][2]; // ok
    BgradQ[12] = -dtTheta;   // ok
    for (int i = 0; i < 3; i++) BgradQ[13+i] = -dtGhat[i]; // buggy
    BgradQ[16] = -dtalpha;
    for (int i = 0; i < 3; i++) BgradQ[17+i] = -dtbeta[i];
    for (int i = 0; i < 3; i++) BgradQ[20+i] = -dtbb[i];
    for (int i = 0; i < 3; i++) BgradQ[23+i] = -dtA[i];
    BgradQ[26] = -dtB[0][0]; // note: thes are all 0 for default sk=0
    BgradQ[27] = -dtB[1][0];
    BgradQ[28] = -dtB[2][0];
    BgradQ[29] = -dtB[0][1];
    BgradQ[30] = -dtB[1][1];
    BgradQ[31] = -dtB[2][1];
    BgradQ[32] = -dtB[0][2];
    BgradQ[33] = -dtB[1][2];
    BgradQ[34] = -dtB[2][2];
    BgradQ[35] = -dtD[0][0][0];
    BgradQ[36] = -dtD[0][0][1];
    BgradQ[37] = -dtD[0][0][2];
    BgradQ[38] = -dtD[0][1][1];
    BgradQ[39] = -dtD[0][1][2];
    BgradQ[40] = -dtD[0][2][2];
    BgradQ[41] = -dtD[1][0][0];
    BgradQ[42] = -dtD[1][0][1];
    BgradQ[43] = -dtD[1][0][2];
    BgradQ[44] = -dtD[1][1][1];
    BgradQ[45] = -dtD[1][1][2];
    BgradQ[46] = -dtD[1][2][2];
    BgradQ[47] = -dtD[2][0][0];
    BgradQ[48] = -dtD[2][0][1];
    BgradQ[49] = -dtD[2][0][2];
    BgradQ[50] = -dtD[2][1][1];
    BgradQ[51] = -dtD[2][1][2];
    BgradQ[52] = -dtD[2][2][2];
    BgradQ[53] = -dtTraceK;
    BgradQ[54] = -dtphi;
    for (int i = 0; i < 3; i++) BgradQ[55+i] = -dtP[i];
    BgradQ[58] = 0;

    // Output data
    
}
#pragma omp end declare target

int main()
{
  const int nVar(59);
  double gradQSerialised[nVar*3], BgradQ[nVar]={2}, Q[nVar]={4};

  // Set up initial test data
  for (int i=0;i<nVar;i++)
  {
    Q[i] = 1;
    BgradQ[i] = 1;
    gradQSerialised[i]=1;
    gradQSerialised[i +nVar]=1;
    gradQSerialised[i +nVar + nVar]=1;
  }
   //The initial data for Q is chosen such that we do not divide by 0 in the inverse determinant
  Q[0]=2;
  Q[1]=3;

  /*for (int i=0;i<1000000;i++) */pdencp_(BgradQ, Q, gradQSerialised, 0);
  //pdencp_(BgradQ, Q, gradQSerialised, 0);

  for (int i=0;i<nVar;i++) std::cout << i+1 << " "  << BgradQ[i] << "\n";
  return 0;

}
