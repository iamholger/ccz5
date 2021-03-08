#include <cmath>
#include <iostream>

void pdeeigenvalues_(double* lambda, const double* const Q, int normal)
{
  constexpr double sqrtwo = 1.4142135623730951;
  const double qmin = std::min({Q[0],Q[3],Q[5]});
  const double test = Q[16]+Q[54];
#if defined(CCZ4EINSTEIN)  
    //const double alpha = std::fmax( 1.0, std::exp(Q[16]) ) * std::fmax(1.0, std::exp(Q[54])) / std::sqrt(qmin);     
    double alpha =std::exp(fmax( 0, Q[16]+ Q[54]))/ std::sqrt(qmin);
    //double alpha =  1. / std::sqrt(qmin);
    //if (test>0) alpha *= std::exp(test);
#else
    const double alpha = 1.0;
#endif
    const double CCZ4e=0.1, CCZ4ds=0.11, CCZ4GLMc = 0.75, CCZ4GLMd = 0.75;
    const double tempA = alpha * std::max({sqrtwo, CCZ4e, CCZ4ds, CCZ4GLMc/alpha, CCZ4GLMd/alpha});
    const double tempB = Q[17+normal];//DOT_PRODUCT(Q(18:20),nv(:))
    lambda[1] = -tempA-tempB; 
    lambda[2] =  tempA-tempB;
}

int main()
{
  double lambda[3] = {0,0,0};
  double Q[59] = {0.144127786,0.603628755,0.273857057,0.476581693,0.392365575,0.673205435,0.361241937,1.21842027E-02,0.368709803,0.219898820,0.686351418,3.99175286E-02,0.290666640,0.234226882,0.509774327,0.604918957,0.916497707,0.992717385,0.445318341,0.647655010,0.107183695,6.34464025E-02,0.788047969,0.525901496,0.229169130,0.867747128,9.73913670E-02,0.746887863,5.47619462E-02,0.830878258,0.493397176,0.887370348,0.567643166,0.277938724,0.891554296,0.273388088,0.411233306,0.577602625,0.853790164,0.622290909,0.996963739,0.783919871,0.103851914,0.829477191,1.94624662E-02,0.860917091,0.317981422,0.244459093,6.99177980E-02,0.904407263,0.450430274,0.204707682,0.810942411,0.191997290,0.257384479,0.758928418,0.499860823,0.198119998,9.39652920E-02};

  for (int i=0;i<1000000000;i++)
    pdeeigenvalues_(lambda,Q, 0);
  std::cout << lambda[0] << " " << lambda[1] << " "  << lambda[2] << "\n";

  return 0;
}
