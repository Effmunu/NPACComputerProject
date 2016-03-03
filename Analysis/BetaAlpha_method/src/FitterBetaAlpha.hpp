#ifndef __FITTERBETAALPHA__
#define __FITTERBETAALPHA__

#include "FitterBase.hpp"


class FitterBetaAlpha: public FitterBase
{

protected :
    int m_nbin;
    int m_npar;
    int** m_mapToLinear;

    double m_norm;
    double m_mZ;
    double m_sigma;
    double m_gamma;

    std::vector<double> m_Beta;
    std::vector<double> m_BetaEr;

public :

    FitterBetaAlpha(MappingTool* map=NULL,std::string name = "TOTO");
    void SetParameters(double norm, double mZ, double sigma, double gamma);

    static int GetLinearIndex(int i, int j); // Get linear index of triangular superior matrix from matrix indices (i,j)
    static int GetMatrixI(int k); // inverse method to get i from k
    static int GetMatrixJ(int k);// inverse method to get j from k

    //get betas
    std::vector<double> GetBetas() {return m_Beta;}

    //get errors
    std::vector<double> GetBetaErs() {return m_BetaEr;}

    //Virtual method, to be implemented in daughter class
    virtual void Execute();
    virtual void FcnForMinuit(int& npar2, double* deriv, double& f, double par[], int flag);
};


#endif
