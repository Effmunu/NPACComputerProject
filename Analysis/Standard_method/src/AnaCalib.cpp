#define AnaCalib_cxx
#include "AnaCalib.hpp"

//root
#include <TMath.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TString.h>

//c++
#include <time.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//my files
#include "MappingTool.hpp"
#include "FitterStandard.hpp"


void AnaCalib::Loop(string& type, string& categ, string& nbEvents, int binning, int stained)
{
    // Measure the running time
    const clock_t startingTime = clock();

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();

    // Prefix for file output
    TString outputPrefix = Form("%s_%s_%s_%d_%s",
        type.c_str(), categ.c_str(), nbEvents.c_str(), binning,
        stained ? "stained" : "unstained");

    //data vector
    vector<InfoForFitter> infoVector;

    //Instanciate the mapping tool
    MappingTool map(binning, -2.4,2.4);
    vector<double> table;
    MappingTool map2(table);

    // Histogram of invariant masses
    TH1F* histInvMass = new TH1F("histInvMass",
        "Invariant Mass " + outputPrefix, 100, 80, 100);

    //===============================================
    // Loop over the events
    //===============================================
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        //========================================
        //Z events preselection
        //========================================
        if( ele0_et<20) continue;
        if( ele1_et<20) continue;
        if( fabs(ele0_eta)>2.4) continue;
        if( fabs(ele1_eta)>2.4) continue;

        double lambda0 = 0;
        double lambda1 = 0;
        double st_ele0_et = 0;
        double st_ele1_et = 0;

        // STAIN THE ENERGY (lambda_i = 0,01 * (-1)^i)
        if (stained) {
            lambda0 = map.getIndex(ele0_eta) % 2 == 0 ? 0.01 : -0.01;
            lambda1 = map.getIndex(ele1_eta) % 2 == 0 ? 0.01 : -0.01;
        }
        st_ele0_et = ele0_et * (1 + lambda0);
        st_ele1_et = ele1_et * (1 + lambda1);

        //========================================
        //Compute the invariant mass
        //========================================
        TLorentzVector el0_LV = TLorentzVector();
        TLorentzVector el1_LV = TLorentzVector();
        el0_LV.SetPtEtaPhiM(st_ele0_et, ele0_eta, ele0_phi, 511e-6);
        el1_LV.SetPtEtaPhiM(st_ele1_et, ele1_eta, ele1_phi, 511e-6);

        double mass = (el0_LV + el1_LV).M();

        // Cuts on invariant mass
        if(mass < 80) continue;
        if(mass > 100) continue;

        histInvMass->Fill(mass);

        //========================================
        //Fill InfoForFitter
        //========================================
        InfoForFitter info;
        info.eta0=ele0_eta;
        info.phi0=ele0_phi;
        info.pt0=st_ele0_et;
        info.region0=map.getIndex(ele0_eta);
        info.eta1=ele1_eta;
        info.phi1=ele1_phi;
        info.pt1=st_ele1_et;
        info.region1=map.getIndex(ele1_eta);
        info.mass=mass;
        infoVector.push_back(info);


    }//end of the event loop

    //========================================
    //Fit
    //========================================
    // Function to get the sigma of the histogram
    TF1* myVoigt = new TF1("myVoigt", "[0] * TMath::Voigt((x - [1]), [2], [3])", 80, 100);
    // Initialization of the parameters
    myVoigt->SetParameter(0, 5000);
    myVoigt->SetParameter(1, 91);
    myVoigt->SetParameter(2, 0.5);
    myVoigt->SetParameter(3, 2);

    // Fit and get the sigma of the histogram
    histInvMass->Fit("myVoigt", "Q"); // Q: quiet mode
    double norm = myVoigt->GetParameter(0);
    double mZ = myVoigt->GetParameter(1);
    double sigma = myVoigt->GetParameter(2);
    double gamma = myVoigt->GetParameter(3);

    //initialize the fitter tool
    FitterStandard fit_standard(&map,"std");
    fit_standard.SetParameters(norm, mZ, sigma, gamma);
    fit_standard.SetData(&infoVector);
    // do the fit
    fit_standard.Execute();
    //get the results
    vector<double > alpha_standard= fit_standard.GetAlphas();
    vector<double > alphaer_standard= fit_standard.GetAlphaErs();

    // End of measure of the running time
    fstream outputStream("timeResults_standard.txt", ios::app);
    outputStream << type << "_" << categ << "_" << nbEvents << "_"
        << binning << "_" << (stained ? "stained" : "unstained") << " - "
        << "Running time: " << double(clock() - startingTime) * 1000. / CLOCKS_PER_SEC << " ms"<< endl;
    outputStream.close();

/*    ///////////
    // Display
    ///////////
    // Plot diff_i = alpha_i - lambda_i and diff_i_over_sigma_i (alpha_i - lambda_i) / sigma_i
    int npar = map.getNbOfBins();
    double x[100], diff_i[100], diff_i_over_sigma_i[100];
    for (int i=0; i<npar; i++) {
        x[i] = i;
        diff_i[i] = alpha_standard[i] -
            (stained ? (i % 2 == 0 ? 0.01 : -0.01) : 0);
        diff_i_over_sigma_i[i] = diff_i[i] / alphaer_standard[i];
    }

    TGraph* graph = new TGraph(npar, x, diff_i);
    TCanvas* canv = new TCanvas("canv", outputPrefix + "_diff_i", 800, 600);
    graph->Draw("AB");
    canv->SaveAs("fig/" + outputPrefix + "_diff_i.png", "Q");

    TGraph* graph2 = new TGraph(npar, x, diff_i_over_sigma_i);
    TCanvas* canv2 = new TCanvas("canv2",
        outputPrefix + "_diff_i_over_sigma_i", 800, 600);
    graph2->Draw("AB");
    canv2->SaveAs("fig/" + outputPrefix + "_diff_i_over_sigma_i.png", "Q");
    // Display
    TCanvas* canv3 = new TCanvas("canv3",
        outputPrefix + "_InvMass", 800, 600);
    histInvMass->Draw();
    canv3->SaveAs("fig/" + outputPrefix + "_InvMass.png", "Q");*/
}
