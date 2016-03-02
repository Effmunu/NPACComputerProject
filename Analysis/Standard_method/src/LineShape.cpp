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
double GBWrel(double xx, double alpha_i, double alpha_j, double norm, double mZ, double sigma, double gamma)
{
    return norm * TMath::Voigt(xx / (1 + (alpha_i + alpha_j) / 2) - mZ, sigma, gamma);
}
