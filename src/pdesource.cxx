#include <cmath>
#include <iostream>
#include "Constants.h"

using namespace examples::exahype2::ccz4;

#pragma omp declare target
void pdesource_(double* S, const double* const Q)
{
#if defined(Dim3)
    constexpr int nDim = 3;
#elif defined(Dim2)
    constexpr int nDim = 2;
#endif
    constexpr int nVar(59), nParam(0), d(3);

    // Input and output variables
    // CCZ4 parameeter -- harmonic lapse CCZ4c
    // Param CCZ4fff, CCZ4mu
    const double fa  = 1.0;
    const double faa = 0.0;

    // De-serialise input data and fill static array
    // FIXME the use of 2D arrays can be avoided: all terms not in the normal are 0


    // Note g_cov is symmetric
    const double g_cov[3][3] = { {Q[0], Q[1], Q[2]}, {Q[1], Q[3], Q[4]}, {Q[2], Q[4], Q[5]} };
    const double det = Q[0]*Q[3]*Q[5] - Q[0]*Q[4]*Q[4] - Q[1]*Q[1]*Q[5] + 2*Q[1]*Q[2]*Q[4] -Q[2]*Q[2]*Q[3];
    const double invdet = 1./det;

    const double g_contr[3][3] = {
        { ( Q[3]*Q[5]-Q[4]*Q[4])*invdet, -( Q[1]*Q[5]-Q[2]*Q[4])*invdet, -(-Q[1]*Q[4]+Q[2]*Q[3])*invdet},
        {-( Q[1]*Q[5]-Q[4]*Q[2])*invdet,  ( Q[0]*Q[5]-Q[2]*Q[2])*invdet, -( Q[0]*Q[4]-Q[2]*Q[1])*invdet},
        {-(-Q[1]*Q[4]+Q[3]*Q[2])*invdet, -( Q[0]*Q[4]-Q[1]*Q[2])*invdet,  ( Q[0]*Q[3]-Q[1]*Q[1])*invdet}
    };


    // NOTE Aex is symmetric
    double Aex[3][3] = { {Q[6], Q[7], Q[8]}, {Q[7], Q[9], Q[10]}, {Q[8], Q[10], Q[11]} };

    double traceA = 0;
    for (int i=0;i<3;i++)
    for (int j=0;j<3;j++) traceA += g_contr[i][j]*Aex[i][j];
    
    for (int i=0;i<3;i++)
    for (int j=0;j<3;j++) Aex[i][j] -= 1./3. * traceA * g_cov[i][j];

    // Matrix multiplications Amix = matmul(g_contr, Aex) Aup  = matmul(g_contr, transpose(Amix))
    double Amix[3][3]={0};
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    for (int u = 0; u < 3; u++) Amix[i][j] += g_contr[i][u] * Aex[u][j];
    double Aup[3][3]={0};
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    for (int u = 0; u < 3; u++) Aup[i][j] += g_contr[i][u] * Amix[j][u]; // Note the transposition is in the indices

    const double Theta = Q[12]; // needed later

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
    for (int j = 0; j < 3; j++) dgup[k][m][l] -= 2.0*g_contr[m][n]*g_contr[j][l]*DD[k][n][j];

    const double PP[3] = {Q[55], Q[56], Q[57]};

    double Christoffel[3][3][3]       = {0};
    double Christoffel_tilde[3][3][3] = {0};
    
    for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
    for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
    {
        Christoffel_tilde[i][j][k] += g_contr[k][l] * ( DD[i][j][l] + DD[j][i][l] - DD[l][i][j] );
        Christoffel[i][j][k]       += g_contr[k][l] * ( DD[i][j][l] + DD[j][i][l] - DD[l][i][j] ) - g_contr[k][l] * ( g_cov[j][l] * PP[i] + g_cov[i][l] * PP[j] - g_cov[i][j] * PP[l] );
    }

    double Gtilde[3] = {0};
    for (int l = 0; l < 3; l++) 
    for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++) Gtilde[i] += g_contr[j][l] * Christoffel_tilde[j][l][i];

    const double Ghat[3] = {Q[13], Q[14], Q[15]};

    const double phi = std::exp(std::fmax(-20., std::fmin(20.,Q[54])));
    const double phi2 = phi*phi;

    double Z[3] = {0}; // Matrix vector multiplications
    for (int i=0;i<3;i++)
    for (int j=0;j<3;j++) Z[i] += 0.5*CCZ4ds*( g_cov[i][j]* (Ghat[j] - Gtilde[j]));
    double Zup[3] = {0};
    for (int i=0;i<3;i++)
    for (int j=0;j<3;j++) Zup[i] += phi2 * g_contr[i][j] * Z[j];

    double dChristoffelSrc[3][3][3][3] = {0};
    double dChristoffel_tildeSrc[3][3][3][3] = {0};
    for (int l = 0; l < 3; l++)
    for (int m = 0; m < 3; m++)
    for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
    for (int k = 0; k < 3; k++)
    {
        dChristoffelSrc[k][i][j][m] +=     dgup[k][m][l]*(      DD[i][j][l]+      DD[j][i][l]-DD[l][i][j])
                                      -    dgup[k][m][l]*(g_cov[j][l]*PP[i]+g_cov[i][l]*PP[j]-g_cov[i][j]*PP[l])
                                      -2.0*g_contr[m][l]*(DD[k][j][l]*PP[i]+DD[k][i][l]*PP[j]-DD[k][i][j]*PP[l]);
          
        dChristoffel_tildeSrc[k][i][j][m] += dgup[k][m][l]*(DD[i][j][l]+DD[j][i][l]-DD[l][i][j]);
    }


    double RiemannSrc[3][3][3][3] = {0};
    for (int m = 0; m < 3; m++)
    for (int j = 0; j < 3; j++)
    for (int k = 0; k < 3; k++)
    for (int i = 0; i < 3; i++)
    {
      RiemannSrc[i][k][j][m] = dChristoffelSrc[k][i][j][m] - dChristoffelSrc[j][i][k][m];
      for (int l = 0; l < 3; l++) RiemannSrc[i][k][j][m] += Christoffel[i][j][l]*Christoffel[l][k][m] - Christoffel[i][k][l]*Christoffel[l][j][m];
    }

    double RicciSrc[3][3] = {0};
    for (int l = 0; l < 3; l++)
    for (int n = 0; n < 3; n++)
    for (int m = 0; m < 3; m++) RicciSrc[m][n] += RiemannSrc[m][l][n][l];

    // Here we directly compute the derivative of Gtilde from its original definition as contracted Christoffel symbol,
    // without assuming unit determinant of the conformal metric. Back to the roots, and as few assumptions as possible...
    double dGtildeSrc[3][3] = {0};
    for (int l = 0; l < 3; l++) 
    for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
    for (int k = 0; k < 3; k++) dGtildeSrc[k][i] += dgup[k][j][l]*Christoffel_tilde[j][l][i] + g_contr[j][l]*dChristoffel_tildeSrc[k][j][l][i];

    double dZSrc[3][3] = {0};
    for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
    for (int k = 0; k < 3; k++) dZSrc[k][i] += CCZ4ds*(DD[k][i][j]*(Ghat[j]-Gtilde[j]) - 0.5*g_cov[i][j]*dGtildeSrc[k][j]);
    
    double nablaZSrc[3][3] = {0};
    for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
    {
      nablaZSrc[i][j] = dZSrc[i][j];
      for (int k = 0; k < 3; k++) nablaZSrc[i][j] -= Christoffel[i][j][k]*Z[k];
    }
    
    double RicciPlusNablaZSrc[3][3];
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) RicciPlusNablaZSrc[i][j] = RicciSrc[i][j] + nablaZSrc[i][j] + nablaZSrc[j][i];

    double RPlusTwoNablaZSrc = 0; 
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) RPlusTwoNablaZSrc += g_contr[i][j]*RicciPlusNablaZSrc[i][j]; // TODO fuse these steps
    RPlusTwoNablaZSrc*=phi2;


    const double alpha = std::exp(std::fmax(-20., std::fmin(20.,Q[16])));

    const double AA[3] = {Q[23], Q[24], Q[25]};

    double nablaijalphaSrc[3][3];
    for (int j = 0; j < 3; j++) 
    for (int i = 0; i < 3; i++)
    {
      nablaijalphaSrc[i][j] = alpha*AA[j]*AA[i]; 
      for (int k = 0; k < 3; k++) nablaijalphaSrc[i][j] -= alpha*Christoffel[i][j][k]*AA[k];
    }
    
    double nablanablaalphaSrc=0;
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) nablanablaalphaSrc += g_contr[i][j]*nablaijalphaSrc[i][j];
    nablanablaalphaSrc *= phi2;

   
    double SecondOrderTermsSrc[3][3];
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) SecondOrderTermsSrc[i][j] = -nablaijalphaSrc[i][j] + alpha*RicciPlusNablaZSrc[i][j];
    
    double traceSrc = 0;
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) traceSrc += g_contr[i][j]*SecondOrderTermsSrc[i][j];
    
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) SecondOrderTermsSrc[i][j] -= 1./3 * traceSrc * g_cov[i][j];
       
    double dtgamma[3][3];
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) dtgamma[i][j] = -2.0 * alpha * Aex[i][j] - CCZ4itau*(det -1.0)*g_cov[i][j];

    const double BB[3][3] = {
        {Q[26], Q[27], Q[28]}, {Q[29], Q[30], Q[31]}, {Q[32], Q[33], Q[34]}
    };

    double traceB = BB[0][0] + BB[1][1] + BB[2][2];
    const double beta[3] = {Q[17], Q[18], Q[19]};
    
    for (int k = 0; k < 3; k++) 
    for (int j = 0; j < 3; j++) 
    for (int i = 0; i < 3; i++) dtgamma[i][j] += g_cov[k][i]*BB[j][k] + g_cov[k][j]*BB[i][k] - 2./3. *g_cov[i][j]*BB[k][k] + 2.*beta[k]*DD[k][i][j];


    // MATMUL(Aex,Amix)
    double Atemp[3][3]={0};
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    for (int u = 0; u < 3; u++) Atemp[i][j] += Aex[i][u] * Amix[u][j];

    const double traceK = Q[53];

    //! Main variables of the CCZ4 system 
    double dtK[3][3];
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) dtK[i][j] = phi2*SecondOrderTermsSrc[i][j] + alpha*Aex[i][j]*(traceK-2.*Theta) - 2.*alpha*Atemp[i][j] - CCZ4itau*g_cov[i][j]*traceA;
    
    for (int j = 0; j < 3; j++) 
    for (int i = 0; i < 3; i++)
    for (int k = 0; k < 3; k++) dtK[i][j] += Aex[k][i]*BB[j][k] + Aex[k][j]*BB[i][k] - 2./3.*Aex[i][j]*BB[k][k];
   
    const double K0 = Q[58];
    
    double dtTraceK = -nablanablaalphaSrc + alpha*(RPlusTwoNablaZSrc + traceK*traceK - 2.0*CCZ4c*Theta*traceK) -3.0*alpha*CCZ4k1*(1.+CCZ4k2)*Theta;
    double dtphi = beta[0]*PP[0] + beta[1]*PP[1] + beta[2]*PP[2] + 1./3.*(alpha*traceK-traceB);
    double dtalpha = -alpha*fa*(traceK-K0-2.*CCZ4c*Theta)+beta[0]*AA[0]+beta[1]*AA[1]+beta[2]*AA[2];


    double Aupdown = 0;
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) Aupdown += Aex[i][j]*Aup[i][j];

  
    double sumzupaa = 0.0;
    for (int i = 0; i < 3; i++) sumzupaa += Zup[i]*AA[i];
    const double dtTheta = 0.5*alpha*CCZ4e*CCZ4e*(RPlusTwoNablaZSrc - Aupdown + 2./3.*traceK*traceK) - alpha*(CCZ4c*Theta*traceK + sumzupaa + CCZ4k1*(2.+CCZ4k2)*Theta);  // Baojiu


    double dtGhat[3];
    for (int i = 0; i < 3; i++)
    {
        double temp1=0, temp2=0, temp3=0, temp4=0, temp5=0, temp6=0;
        for (int m = 0; m < 3; m++)
        {
          temp1 += Aup[i][m]*PP[m];
          temp3 += Aup[i][m]*AA[m];
          temp2 += g_contr[m][i]*(-Theta*AA[m] -2./3.*traceK*Z[m]);
          temp4 += g_contr[i][m]*Z[m];
          temp5 += Gtilde[m]*BB[m][i];
          for (int n = 0; n < 3; n++) temp6  += Christoffel_tilde[m][n][i]*Aup[m][n];
        }
        dtGhat[i] += 2.*alpha*(temp6 - 3.*temp1 + temp2 - temp3 - CCZ4k1*temp4) - temp5 + 2./3.*Gtilde[i]*traceB;

        for (int l = 0; l < 3; l++)
        for (int k = 0; k < 3; k++) dtGhat[i] += 2.*CCZ4k3*(2./3.*g_contr[i][l]*Z[l]*BB[k][k] - g_contr[l][k]*Z[l]*BB[k][i]);
    }

    double ov[3];
    for (int k = 0; k < 3; k++)
    {
        double temp=0;
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) temp += dgup[k][i][j]*Aex[i][j];
        ov[k] = 2*alpha*temp;
    }

    // matrix vector multiplication in a loop and add result to existing vector
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) dtGhat[i] += CCZ4sk*g_contr[i][j]*ov[j];
    
    const double myb[3] = {Q[20], Q[21], Q[22]};
    
    double dtbb[3];
    for (int i = 0; i < 3; i++) dtbb[i] = CCZ4sk*(CCZ4xi*dtGhat[i] - CCZ4eta*myb[i]);

    double dtbeta[3];
    for (int i = 0; i < 3; i++) dtbeta[i] = CCZ4f*myb[i];
    for (int i = 0; i < 3; i++) dtbeta[i] += CCZ4bs*(beta[0]*BB[0][i] + beta[1]*BB[1][i] + beta[2]*BB[2][i]);
    for (int i = 0; i < 3; i++) dtbeta[i] *= CCZ4sk;


    // Auxiliary variables 
    double dtA[3];
    for (int i = 0; i < 3; i++)
    {
      dtA[i] = -alpha*AA[i]*(fa+alpha*faa)*(traceK - K0 - 2.*CCZ4c*Theta);
      for (int j = 0; j < 3; j++) dtA[i] += BB[i][j]*AA[j];
    }

    for (int k = 0; k < 3; k++)
    {
      double temp = 0;
      for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) temp+= dgup[k][i][j]*Aex[i][j];
      dtA[k] += -CCZ4sk*alpha*fa*temp;
    }

    double dtB[3][3] ={0};
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    for (int u = 0; u < 3; u++) dtB[i][j] += CCZ4sk*(BB[i][u] * BB[u][j]);

    double dtD[3][3][3];
    for (int m = 0; m < 3; m++)
    for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
    for (int k = 0; k < 3; k++)
    for (int n = 0; n < 3; n++) dtD[k][i][j] += 1./3*alpha*g_cov[i][j]*dgup[k][n][m]*Aex[n][m]; // explicitly remove the trace of tilde A again  

    for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
    for (int k = 0; k < 3; k++)
    {
      dtD[k][i][j] -= alpha*AA[k]*Aex[i][j];
      for (int l = 0; l < 3; l++) dtD[k][i][j] += BB[k][l]*DD[l][i][j] + BB[j][l]*DD[k][l][i] + BB[i][l]*DD[k][l][j] - 2./3.*BB[l][l]*DD[k][i][j];
    }

    double dtP[3] = {0};
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) dtP[i] += BB[i][j]*PP[j];

    for (int k = 0; k < 3; k++)
    {
      double temp=0;
      for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) temp += dgup[k][i][j]*Aex[i][j];
      dtP[k] += 1./3.*alpha*(AA[k]*traceK + CCZ4sk*temp);
    }


    S[0]  = dtgamma[0][0];
    S[1]  = dtgamma[0][1];
    S[2]  = dtgamma[0][2];
    S[3]  = dtgamma[1][1];
    S[4]  = dtgamma[1][2];
    S[5]  = dtgamma[2][2];
    S[6]  = dtK[0][0];
    S[7]  = dtK[0][1];
    S[8]  = dtK[0][2];
    S[9]  = dtK[1][1];
    S[10] = dtK[1][2];
    S[11] = dtK[2][2]; 
    S[12] = dtTheta;   
    for (int i = 0; i < 3; i++) S[13+i] = dtGhat[i];
    S[16] = dtalpha;
    for (int i = 0; i < 3; i++) S[17+i] = dtbeta[i];
    for (int i = 0; i < 3; i++) S[20+i] = dtbb[i];
    for (int i = 0; i < 3; i++) S[23+i] = dtA[i];
    S[26] = dtB[0][0];
    S[27] = dtB[1][0];
    S[28] = dtB[2][0];
    S[29] = dtB[0][1];
    S[30] = dtB[1][1];
    S[31] = dtB[2][1];
    S[32] = dtB[0][2];
    S[33] = dtB[1][2];
    S[34] = dtB[2][2];
    S[35] = dtD[0][0][0];
    S[36] = dtD[0][0][1];
    S[37] = dtD[0][0][2];
    S[38] = dtD[0][1][1];
    S[39] = dtD[0][1][2];
    S[40] = dtD[0][2][2];
    S[41] = dtD[1][0][0];
    S[42] = dtD[1][0][1];
    S[43] = dtD[1][0][2];
    S[44] = dtD[1][1][1];
    S[45] = dtD[1][1][2];
    S[46] = dtD[1][2][2];
    S[47] = dtD[2][0][0];
    S[48] = dtD[2][0][1];
    S[49] = dtD[2][0][2];
    S[50] = dtD[2][1][1];
    S[51] = dtD[2][1][2];
    S[52] = dtD[2][2][2];
    S[53] = dtTraceK;
    S[54] = dtphi;
    for (int i = 0; i < 3; i++) S[55+i] = dtP[i];
    S[58] = 0;
}
#pragma omp end declare target

int main()
{
#pragma omp target teams
  {
  const int nVar(59);
  double Q[nVar], S[nVar];

  // Set up initial test data
  for (int i=0;i<nVar;i++)
  {
    Q[i] = 1;
    S[i] = 0;
  }
   //The initial data for Q is chosen such that we do not divide by 0 in the inverse determinant
  Q[0]=2;
  Q[1]=3;

//#pragma omp distribute parallel for
  //for (int i=0;i<10;i++)
    pdesource_(S,Q);

  //printf("%f\n", temppp);

  for (int i=0;i<59;i++) printf("\t%d\t%.30f\n", i+1, S[i]);
  }
  return 0;
}
