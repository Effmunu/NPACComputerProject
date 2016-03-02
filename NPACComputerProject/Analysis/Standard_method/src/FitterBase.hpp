//===============================================
//
//             Master NPAC
// Virtual class common to all fitting methods
// Each fitter inherits from this class
//
// the calibration coefficient are stored in a std::vector
//
//===============================================


#ifndef __FITTERBASE__
#define __FITTERBASE__

//root
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TGraph.h>

#include "MappingTool.hpp"
#include "InfoForFitter.h"
#include "LineShape.hpp"

#include <string>
#include <fstream>


class FitterBase
{
  
protected :

  //name associated to the instance of the class
  std::string m_name;
  
  // Mapping Tool 
  MappingTool* m_map;  
  
  //fitted b
  std::vector<double> m_Alpha;    //calibration factors
  std::vector<double> m_AlphaEr;  // and its errors
  
  //data
  std::vector< InfoForFitter >* m_InfoVector;
  

public:
  //constructors
  FitterBase(MappingTool* map=NULL,std::string name = "TOTO");

  //destructor
  virtual ~FitterBase(){};

  //setter for data
  void SetData(std::vector< InfoForFitter >* InfoVector);
  
  //get alphas
  std::vector<double> GetAlphas() {return m_Alpha;}

  //get errors 
  std::vector<double> GetAlphaErs() {return m_AlphaEr;}

  //Virtual method where the fit is performed, to be implemented in a daughter class
  virtual void Execute()=0;

  //Fcn to minimize by minuit
  virtual void FcnForMinuit(int& npar, double* deriv, double& f, double par[], int flag)=0;

  //Write the alphas into a text file
  void AlphaToFile(string name);



};




#endif
