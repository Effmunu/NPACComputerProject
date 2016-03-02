#include "MappingTool.hpp"
#include <cmath>

//ClassImp(MappingTool)

MappingTool::MappingTool():TObject()
{
  SetEtaBins(24, -2.4,2.4);
}

MappingTool::MappingTool(  std::vector<double> etaBorders):TObject()
{
  SetEtaBins(etaBorders);
}


MappingTool::MappingTool(  int Neta, double etamin, double etamax):TObject()
{
   SetEtaBins(Neta,etamin,etamax);
}



void MappingTool::SetEtaBins(std::vector<double> etaBorders)
{
  m_EtaBorders.clear();
  m_EtaBorders=etaBorders;
}

void MappingTool::SetEtaBins(int Neta, double etamin, double etamax)
{
  m_EtaBorders.clear();
  m_EtaBorders.resize(Neta+1);

  double step=(etamax-etamin)/ Neta;
  for(int i=0;i<Neta+1;i++)
    {
      m_EtaBorders[i]=etamin+i*step;
    }
}


void MappingTool::MyPrint()
{ 
  
  cout<<"Eta bins= ";
  for (int i=0; i<m_EtaBorders.size();i++)
    {
      cout<<m_EtaBorders[i]<<" ";
    }
  cout<<endl;
  
}

int MappingTool::getIndex(const double &eta)
{ 
  int iEta=-1;
  for(int ieta = 0;ieta<m_EtaBorders.size()-1;ieta++)
    {
      if(m_EtaBorders[ieta]<=eta  && eta<m_EtaBorders[ieta+1])
	{
	  iEta=ieta;
	  break;
	}
    }
  if(iEta<0 || iEta>= getNbOfBins()) cout<<"PRB  y="<<eta<<endl;
  
  return iEta;
}


void MappingTool::getEta(int index,double& etaCent)
{
  etaCent= m_EtaBorders[index]+(m_EtaBorders[index+1]-m_EtaBorders[index])/2;
}



double  MappingTool::getEtaBinSize(int index)
{
  return m_EtaBorders[index+1]-m_EtaBorders[index];
}

