#include "FitterBase.hpp"

#include <iomanip>
using namespace std;

FitterBase::FitterBase(MappingTool* map,std::string name)
{
  m_map=map;
  m_name=name;
  m_InfoVector=NULL;
  if(map!=NULL)
    {
      m_Alpha.resize(m_map->getNbOfBins());
      m_AlphaEr.resize(m_map->getNbOfBins());
      for(int i = 0;i<m_Alpha.size();i++)
	{
	  m_Alpha[i]     = -999; 
	  m_AlphaEr[i]   = -999;
	}
    }
}
  


void FitterBase::SetData(std::vector< InfoForFitter > *InfoVector)
{
  m_InfoVector=InfoVector;
}

void FitterBase::AlphaToFile(string filename)
{

 std::ofstream *outputFile = new ofstream(filename.c_str());
  if(outputFile==0)
    { 
      cerr<<"*********************************************************************"<<endl;
      cerr<<"FitterBase::WriteAlphaInAFile(string filename) : " <<endl;
      cerr<<"  Cannot create : " <<filename  <<endl; 
      cerr<<"*********************************************************************"<<endl;
      //      assert(0);
    }
  
  for(int i=0;i<m_Alpha.size();i++)
    {
	  *outputFile  <<i
		       <<" " 
		       <<setw(10)<<m_Alpha[i]<<"  "
		       <<setw(10)<<m_AlphaEr[i]<<"  "
		       <<endl;
    }
  outputFile->close();
  delete  outputFile;


}
