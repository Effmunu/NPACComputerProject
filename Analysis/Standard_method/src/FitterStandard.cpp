//===============================================
//
//             Master NPAC
// Implementation of the standard method
// Inherit from FitterBase class
//
//===============================================


#include "FitterStandard.hpp"


#include "TMinuit.h"

FitterStandard::FitterStandard(MappingTool* map,std::string name):FitterBase(map,name){}

void FitterStandard::SetParameters(double norm, double mZ, double sigma, double gamma)
{
    m_norm = norm;
    m_mZ = mZ;
    m_sigma = sigma;
    m_gamma = gamma;
}

static FitterStandard* static_pointer_standard = NULL;

void fcn_wrapper_standard(int &npar, double *gin, double &f, double *par, int iflag) {
    if (static_pointer_standard==NULL) {
        std::cout << "static_pointer_standard is NULL!!!!!" << std::endl;
        //    exit(-1);
    }
    static_pointer_standard->FcnForMinuit(npar, gin, f, par, iflag);
}



void FitterStandard::Execute()
{
    //==========================================
    //Init
    //==========================================
    static_pointer_standard   = this;
    int npar = m_map->getNbOfBins();


    //==========================================
    //TMinuit settings
    //==========================================
    TMinuit minuit(npar);
    minuit.SetPrintLevel(-1);
    minuit.SetFCN(fcn_wrapper_standard);


    double par[npar];               // the start values
    double stepSize[npar];          // step sizes
    double minVal[npar];            // minimum bound on parameter
    double maxVal[npar];            // maximum bound on parameter
    std::vector<std::string> parName;
    parName.resize(npar);

    for(int i=0; i<npar;i++)
    {
        char tmp[100];
        par[i]= 0;
        stepSize[i] = 0.001;
        minVal[i] = -0.3;
        maxVal[i] = +0.3;
        sprintf(tmp,"%s%d","alpha",i);
        parName[i] = tmp;
    }

    for (int i=0; i<npar; i++)
    {
        minuit.DefineParameter(i, parName[i].c_str(),
                par[i], stepSize[i], minVal[i], maxVal[i]);
    }

    //==========================================
    //Do the minimization
    //==========================================
    minuit.SetPrintLevel(0);
    minuit.Migrad();

    //==========================================
    //Get results
    //==========================================

    for (int i=0; i<npar; i++)
    {
        minuit.GetParameter(i,m_Alpha[i],m_AlphaEr[i]);
//        cout<<"alpha"<<i<<"  "<<m_Alpha[i] <<"  "<<m_AlphaEr[i]<<endl;
    }
}




void FitterStandard::FcnForMinuit(int& npar, double* deriv, double& f, double par[], int flag)
{
    //===================================
    //Write here the function to minimize
    //===================================
    f = 0;
    for (int entry = 0; entry<m_InfoVector->size(); entry++) {
    //for (int entry = 0; entry<10; entry++) {
        /*cout << "mass: " << (*m_InfoVector)[entry].mass
             << "\t region0: " << par[(*m_InfoVector)[entry].region0]
             << "\t region1: " << par[(*m_InfoVector)[entry].region1]
             << "\t sigma: " << m_sigma
             << "\t gamma: " << m_gamma << endl;*/

        f += -2 * log(  GBWrel((*m_InfoVector)[entry].mass,
                        par[(*m_InfoVector)[entry].region0],
                        par[(*m_InfoVector)[entry].region1],
                        m_norm,
                        m_mZ,
                        m_sigma,
                        m_gamma));
    }
}
