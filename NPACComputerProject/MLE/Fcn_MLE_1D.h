#ifndef Fcn_MLE_1D
#define Fcn_MLE_1D

#include <TMath.h>
#include <cmath>

const int N=100;
float data[N]={1.4,-0.0348,1.18,0.37,1.22,0.343,-0.501,0.325,0.408,-0.0108,1.79,-0.585,0.351,-1.04,-0.661,-0.988,1.17,-0.336,0.98,0.0179,2.46,-0.835,1.57,-0.0542,0.265,-0.0996,0.218,2.24,0.157,2.4,0.405,-0.0222,1.94,0.495,1.92,1.62,0.264,0.201,0.106,0.282,-0.76,0.412,0.422,0.329,-0.922,2.14,0.0676,0.153,-0.78,-0.79,-1.61,1.74,0.888,0.696,0.644,-0.692,-0.199,-0.388,-2.18,-0.816,-0.843,0.162,2.2,-0.0689,-1.84,0.101,2.18,-0.486,0.656,0.0315,-0.316,0.0132,0.099,0.273,1.25,0.855,0.369,1.6,0.759,-0.159,1.33,1.79,1.38,-0.495,0.678,-0.353,-0.217,-0.899,1.39,0.703,-1.09,1.74,0.62,0.222,0.0681,0.768,-0.0495,1.91,0.358,0.482};

double gaussian(double x, double mu, double sigma)
{
    return (1./sqrt(2*TMath::Pi())/sigma) * exp(- (x-mu)*(x-mu) / 2/sigma/sigma);
}

void fcn_to_minimize(int& npar, double* deriv, double& fvalue, double par[], int flag);

#endif

