#include "LineShape.hpp"

#include <cmath>
using namespace std;
/*
//===========================================
// Relativistic Breit-Wigner
//===========================================
double BWrel(double *xx, double *par)
{


}

//===========================================
// Gaussian
//===========================================
double Gaussian(double *x, double *par)
{

}
*/
//===========================================
// Breit-Wigner convoluted with a Gaussian
//===========================================
double GBWrel(double xx, double beta_ij, double norm, double mZ, double sigma, double gamma)
{
    return norm * TMath::Voigt(xx / (1 + beta_ij / 2) - mZ, sigma, gamma);
}
