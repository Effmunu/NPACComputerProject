#include "MappingTool.hpp"

//ClassImp(MappingTool)

MappingTool::MappingTool():TObject()
{
    SetEtaBins(24, -2.4, 2.4);
    SetPhiBins(32, 0, TMath::Pi());
}

MappingTool::MappingTool(std::vector<double> etaBorders,
                            std::vector<double> phiBorders):TObject()
{
    SetEtaBins(etaBorders);
    SetPhiBins(phiBorders);
}

MappingTool::MappingTool(int Neta, double etamin, double etamax,
                            int Nphi, double phimin, double phimax):TObject()
{
    SetEtaBins(Neta, etamin, etamax);
    SetPhiBins(Nphi, phimin, phimax);
}



void MappingTool::SetEtaBins(std::vector<double> etaBorders)
{
    m_EtaBorders.clear();
    m_EtaBorders=etaBorders;
}

void MappingTool::SetPhiBins(std::vector<double> phiBorders)
{
    m_PhiBorders.clear();
    m_PhiBorders=phiBorders;
}

void MappingTool::SetEtaBins(int Neta, double etamin, double etamax)
{
    m_EtaBorders.clear();
    m_EtaBorders.resize(Neta+1);

    double step = (etamax-etamin) / Neta;
    for(int i=0; i<Neta+1; i++)
        m_EtaBorders[i] = etamin + i*step;
}

void MappingTool::SetPhiBins(int Nphi, double phimin, double phimax)
{
    m_PhiBorders.clear();
    m_PhiBorders.resize(Nphi+1);

    double step = (phimax-phimin) / Nphi;
    for(int i=0; i<Nphi+1; i++)
        m_PhiBorders[i] = phimin + i*step;
}


void MappingTool::MyPrint()
{
    cout<<"Eta bins= ";
    for (int i=0; i<m_EtaBorders.size(); i++)
        cout << m_EtaBorders[i] << " ";
    cout<<endl;

    cout<<"Phi bins= ";
    for (int i=0; i<m_PhiBorders.size(); i++)
        cout << m_PhiBorders[i] << " ";
    cout<<endl;

}


int MappingTool::getIndex(const double &eta, const double &phi)
{
    return getEtaIndex(eta) * (m_PhiBorders.size()-1) + getPhiIndex(phi);
}

int MappingTool::getEtaIndex(const double &eta)
{
    int iEta = -1;
    for(int ieta = 0; ieta < m_EtaBorders.size()-1; ieta++) {
        if(m_EtaBorders[ieta] <= eta && eta < m_EtaBorders[ieta+1]) {
            iEta = ieta;
            break;
        }
    }
    if(iEta < 0 || iEta >= getNbEtaBins())
        cout << "PRB  y=" << eta << endl;

    return iEta;
}

int MappingTool::getPhiIndex(const double &phi)
{
    int iPhi = -1;
    for(int iphi = 0; iphi < m_PhiBorders.size()-1; iphi++) {
        if(m_PhiBorders[iphi]<=phi  && phi < m_PhiBorders[iphi+1]) {
            iPhi = iphi;
            break;
        }
    }
    if(iPhi < 0 || iPhi >= getNbPhiBins())
        cout << "PRB  y=" << phi << endl;

    return iPhi;
}


void MappingTool::getEtaPhi(int index, double& etaCent, double& phiCent)
{
    etaCent = m_EtaBorders[index] +
                (m_EtaBorders[index+1] - m_PhiBorders[index])/2;
    phiCent = m_PhiBorders[index] +
                (m_PhiBorders[index+1] - m_PhiBorders[index])/2;
}



double  MappingTool::getEtaBinSize(int index)
{
    return m_EtaBorders[index+1] - m_EtaBorders[index];
}

double  MappingTool::getPhiBinSize(int index)
{
    return m_PhiBorders[index+1] - m_PhiBorders[index];
}
