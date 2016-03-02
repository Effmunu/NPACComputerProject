//===============================================
//
//             Master NPAC
// Functions describing Z and J/psi line shape
//
//===============================================


#ifndef __ZLINESHAPE__
#define __ZLINESHAPE__

#include <TMath.h>
/*// Relativistic Breit-Wigner
double BWrel(double *xx, double *par);

// Gaussian
double Gaussian(double *x, double *par);*/

// Breit-Wigner convoluted with a Gaussian
double GBWrel(double xx, double beta_ij, double norm, double mZ, double sigma, double gamma);

#endif
