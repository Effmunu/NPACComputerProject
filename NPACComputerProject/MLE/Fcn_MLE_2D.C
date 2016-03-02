#include "Fcn_MLE_2D.h"

using namespace std;

void fcn_to_minimize(int& npar, double* deriv, double& fvalue, double par[], int flag)
{
    fvalue = 0; // initialization
    for(int i=0; i<N; i++) {
        fvalue += -2 * log(gaussian(data[i], par[0], par[1]));
    }
}

