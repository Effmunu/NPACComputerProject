//===============================================
//
//                 Master NPAC
//
// For each eta bins, an integer is associated and vice-versa
//        (eta)    <--->  integer
//
//
//  Constructor 1:
//      MappingTool map(24, -2.4, 2.4, 32, 0, TMath::Pi());
//  Constructor 2:
//      vector<double> tableEta, tablePhi;
//      MappingTool map2(tableEta, tablePhi);
//      this constructor allows bins with different size
//
//===============================================

#ifndef __MAPPINGTOOL__
#define __MAPPINGTOOL__

#include <iostream>
#include <cmath>
using namespace std;

#include <TObject.h>
#include <TMath.h>

class MappingTool: public TObject
{
protected:

    //bin borders
    std::vector<double> m_EtaBorders;
    std::vector<double> m_PhiBorders;


public:

    /** constructor */
    MappingTool();

    /** constructor */
    MappingTool(std::vector<double> etaBorders, std::vector<double> phiBorders);

    /** constructor */
    MappingTool(int Neta, double etamax, double etamin,
                int Nphi, double phimin, double phimax);

    /** Set Eta bins */
    void SetEtaBins(std::vector<double> etaBorders);
    void SetPhiBins(std::vector<double> phiBorders);

    /** Set Eta bins */
    void SetEtaBins(int Neta, double etamin, double etamax);
    void SetPhiBins(int Nphi, double phimin, double phimax);

    /** get index from (eta) */
    int getIndex(const double &eta, const double &phi);
    int getEtaIndex(const double &eta);
    int getPhiIndex(const double &phi);

    /** get (eta) from index */
    void getEtaPhi(int index, double& etaCent, double& phiCent);

    /** get the total number of regions */
    unsigned int  getNbEtaBins(){ return (m_EtaBorders.size()-1);};
    unsigned int  getNbPhiBins(){ return (m_PhiBorders.size()-1);};

    /** get the size of the regions in eta*/
    double  getEtaBinSize(int);
    double  getPhiBinSize(int);

    /** clear the vectors*/
    void clear()
    {
        m_EtaBorders.clear();
        m_PhiBorders.clear();
    }

    /** print the vectors*/
    void MyPrint();

    /** get bins*/
    std::vector<double> getEtaBorders() {return m_EtaBorders;};
    std::vector<double> getPhiBorders() {return m_PhiBorders;};

    // technical (ignore)
    //  ClassDef(MappingTool,1)  //Mapping Tool
};

#endif
