#ifndef __FITTERSTANDARD__
#define __FITTERSTANDARD__

#include "FitterBase.hpp"


class FitterStandard: public FitterBase
{

protected :
    double m_norm;
    double m_mZ;
    double m_sigma;
    double m_gamma;

public :

    FitterStandard(MappingTool* map=NULL,std::string name = "TOTO");
    void SetParameters(double norm, double mZ, double sigma, double gamma);

    //Virtual method, to be implemented in daughter class
    virtual void Execute();
    virtual void FcnForMinuit(int& npar, double* deriv, double& f, double par[], int flag);
};


#endif
