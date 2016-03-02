//===============================================
//
//                 Master NPAC
//
// For each eta bins, an integer is associated and vice-versa
//        (eta)    <--->  integer
//
//
//  Constructor 1:  
//      MappingTool map(  24, -2.4,2.4);
//  Constructor 2:
//      vector<double> table;
//      MappingTool map2(table);
//      this constructor allows bins with different size
//
//===============================================

#ifndef __MAPPINGTOOL__
#define __MAPPINGTOOL__

#include <iostream>
#include <cmath>
using namespace std;

#include <TObject.h>

class MappingTool: public TObject
{
 protected:

  //bin borders
  std::vector<double> m_EtaBorders;


 public:

 /** constructor */
  MappingTool();

 /** constructor */
  MappingTool(  std::vector<double> etaBorders);

  /** constructor */
  MappingTool(  int Neta, double etamax, double etamin);
  
  /** Set Eta bins */
  void SetEtaBins(std::vector<double> etaBorders);

 /** Set Eta bins */
  void SetEtaBins(int Neta, double etamax, double etamin);



  /** get index from (eta) */
  int getIndex(const double &eta);
  
  /** get (eta) from index */
  void getEta(int index,double& etaCent);
  
  /** get the total number of regions */
  unsigned int  getNbOfBins(){ return (m_EtaBorders.size()-1);};
  
  /** get the size of the regions in eta*/
  double  getEtaBinSize(int);

  /** clear the vectors*/
  void clear(){m_EtaBorders.clear();}

  /** print the vectors*/
  void MyPrint();
 
  /** get bins*/
  std::vector<double> getEtaBorders() {return m_EtaBorders;}; 







  // technical (ignore)
  //  ClassDef(MappingTool,1)  //Mapping Tool
};

#endif
