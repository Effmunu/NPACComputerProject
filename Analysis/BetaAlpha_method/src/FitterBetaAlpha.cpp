//===============================================
//
//             Master NPAC
// Implementation of the BetaAlpha method
// Inherit from FitterBase class
//
//===============================================


#include "FitterBetaAlpha.hpp"


#include "TMinuit.h"

FitterBetaAlpha::FitterBetaAlpha(MappingTool* map,
    std::string name):FitterBase(map,name),
    m_nbin(m_map->getNbOfBins()),
    m_npar(m_nbin*(m_nbin+1)/2)
{
    if(map!=NULL) {
        m_Beta.resize(m_npar);
        m_BetaEr.resize(m_npar);
        for(int i = 0; i<m_Beta.size(); i++) {
            m_Beta[i]     = -999;
            m_BetaEr[i]   = -999;
        }
    }
    // Initialize map of triangular sup matrix indices to linear indices
    m_mapToLinear = new int*[m_nbin];
    for (int i = 0; i < m_nbin; i++) {
        m_mapToLinear[i] = new int[m_nbin];
    }
    // Fill the map
    for (int i = 0; i < m_nbin; i++) {
        for (int j = 0; j < m_nbin; j++) {
            m_mapToLinear[i][j] = GetLinearIndex(i, j);
        }
    }
}

void FitterBetaAlpha::SetParameters(double norm, double mZ, double sigma, double gamma)
{
    m_norm = norm;
    m_mZ = mZ;
    m_sigma = sigma;
    m_gamma = gamma;
}

int FitterBetaAlpha::GetLinearIndex(int i, int j)
{
    // We only use the triangular superior part of the beta_ij matrix.
    // So here is the calculation of the linear index of this triangular sup part from (i,j).
    if (i<=j)
        return m_nbin*(m_nbin-1)/2 - (m_nbin-i)*(m_nbin-i-1)/2 + j; // normal case
    else
        return m_nbin*(m_nbin-1)/2 - (m_nbin-j)*(m_nbin-j-1)/2 + i; // i,j inverted to get triangular sup.
}

int FitterBetaAlpha::GetMatrixI(int k)
{
    // We only use the triangular superior part of the beta_ij matrix.
    // So here is a way to get the first matrix index from the linear index of this triangular sup part.
    int i = 0;
    while(k+1 > (i+1)*m_nbin - i*(i+1)/2)
        i++;
    return i;
}

int FitterBetaAlpha::GetMatrixJ(int k)
{
    // We only use the triangular superior part of the beta_ij matrix.
    // So here is a way to get the second matrix index from the linear index of this triangular sup part.

    // Just invert the linear index formula by taking the i index with GetMatrixI
    int i = GetMatrixI(k);
    return k - m_nbin*(m_nbin-1)/2 + (m_nbin-i)*(m_nbin-i-1)/2;
}

static FitterBetaAlpha* static_pointer_betaAlpha = NULL;

void fcn_wrapper_betaAlpha(int &npar, double *gin, double &f, double *par, int iflag) {
    if (static_pointer_betaAlpha==NULL) {
        std::cout << "static_pointer_betaAlpha is NULL!!!!!" << std::endl;
        //    exit(-1);
    }
    static_pointer_betaAlpha->FcnForMinuit(npar, gin, f, par, iflag);
}



void FitterBetaAlpha::Execute()
{
    //==========================================
    //Init
    //==========================================
    static_pointer_betaAlpha   = this;

    //==========================================
    //TMinuit settings
    //==========================================
    TMinuit minuit(m_npar);
    minuit.SetPrintLevel(-1);
    minuit.SetFCN(fcn_wrapper_betaAlpha);


    double par[m_npar];               // the start values
    double stepSize[m_npar];          // step sizes
    double minVal[m_npar];            // minimum bound on parameter
    double maxVal[m_npar];            // maximum bound on parameter
    std::vector<std::string> parName;
    parName.resize(m_npar);

    for(int k=0; k<m_npar;k++)
    {
        char tmp[100];
        par[k]= 0;
        stepSize[k] = 0.001;
        minVal[k] = -0.5;
        maxVal[k] = +0.5;
        sprintf(tmp,"%s%d_%d","beta",GetMatrixI(k),GetMatrixJ(k));
        parName[k] = tmp;
    }

    for (int k=0; k<m_npar; k++)
    {
        minuit.DefineParameter(k, parName[k].c_str(),
                par[k], stepSize[k], minVal[k], maxVal[k]);
    }

    //==========================================
    //Do the minimization
    //==========================================
    minuit.SetPrintLevel(0);
    minuit.Migrad();

    //==========================================
    //Get results
    //==========================================

    for (int k=0; k<m_npar; k++)
    {
        minuit.GetParameter(k, m_Beta[k], m_BetaEr[k]);
        //cout << "beta" << GetMatrixI(k) << "_" << GetMatrixJ(k) << "  " << m_Beta[k] << "  " << m_BetaEr[k] << endl;
    }
}




void FitterBetaAlpha::FcnForMinuit(int& npar, double* deriv, double& f, double par[], int flag)
{
    //===================================
    //Write here the function to minimize
    //===================================
    f = 0;
    for (int entry = 0; entry<m_InfoVector->size(); entry++) {
    /*    if( ((*m_InfoVector)[entry].region0 == 2 && (*m_InfoVector)[entry].region1 == 11) || ((*m_InfoVector)[entry].region0 == 11 && (*m_InfoVector)[entry].region1 == 2) ){
            cout << "entry: " << entry << endl;
        }*/
    /*for (int entry = 0; entry<10; entry++) {
        cout << "mass: " << (*m_InfoVector)[entry].mass
             << "\t ij: " << (*m_InfoVector)[entry].region0 << ", " << (*m_InfoVector)[entry].region1
             << "\t n*i+j: " << npar * (*m_InfoVector)[entry].region0 +
             (*m_InfoVector)[entry].region1
             << "\t beta_ij: " << par[npar * (*m_InfoVector)[entry].region0 +
             (*m_InfoVector)[entry].region1]
             << endl;*/

        f += -2 * log(  GBWrel((*m_InfoVector)[entry].mass,
                        par[m_mapToLinear[(*m_InfoVector)[entry].region0][(*m_InfoVector)[entry].region1]],
                        m_norm,
                        m_mZ,
                        m_sigma,
                        m_gamma));
    }
    //cout << "Likelihood: " << f << endl;
}
