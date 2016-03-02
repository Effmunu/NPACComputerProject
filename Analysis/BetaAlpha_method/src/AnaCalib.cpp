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
#include <TPaveText.h>
#include <TLegend.h>

//c++
#include <time.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//my files
#include "MappingTool.hpp"
#include "FitterBetaAlpha.hpp"


void AnaCalib::Loop(string& type, string& categ, string& nbEvents, int binning, int stained)
{
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetTitleOffset(1.4, "y");

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
    double lambda0 = 0;
    double lambda1 = 0;
    double st_ele0_et = 0;
    double st_ele1_et = 0;

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
    FitterBetaAlpha fit_betaAlpha(&map,"std");
    fit_betaAlpha.SetParameters(norm, mZ, sigma, gamma);
    fit_betaAlpha.SetData(&infoVector);
    // do the fit
    fit_betaAlpha.Execute();
    //get the results
    vector<double > beta_betaAlpha= fit_betaAlpha.GetBetas();
    vector<double > betaer_betaAlpha= fit_betaAlpha.GetBetaErs();

    // Export results to file (because it's looong to run...)
    fstream stream_out( (outputPrefix + "_betas.txt").Data(), ios::out);
    for (int i=0; i<beta_betaAlpha.size(); i++) {
        stream_out << beta_betaAlpha[i] << " " << betaer_betaAlpha[i] << endl;
    }
    stream_out.close();

    // End of measure of the running time
    fstream outputStream("timeResults_betaAlpha.txt", ios::app);
    outputStream << type << " " << categ << " " << nbEvents << " "
        << binning << "\t" << (stained ? "stained" : "unstained") << "\t"
        << double(clock() - startingTime) / CLOCKS_PER_SEC << " s"<< endl;
    outputStream.close();

    // Display
    TPaveText* info_text = new TPaveText(0.63, 0.65, 0.88, 0.8, "ndc");
    info_text->SetBorderSize(0);
    info_text->SetTextSize(0.04);
    info_text->SetFillColor(kWhite);
    info_text->AddText(Form("%s %s %s, %s", type.c_str(), categ.c_str(), nbEvents.c_str(), (stained ? "stained" : "unstained")));
    info_text->AddText(Form("Nb bins: %d", binning));

    TLegend* leg= new TLegend(0.20, 0.65, 0.45, 0.8);
    leg->SetTextSize(0.04);
    leg->AddEntry(histInvMass, "Data", "F");
    leg->AddEntry(myVoigt, "Voigtian fit", "L");

    histInvMass->GetXaxis()->SetTitle("m_{ee} [GeV]");
    histInvMass->GetYaxis()->SetTitle("Number of event / 0.2 GeV");
    TCanvas* canv = new TCanvas("canv",
        outputPrefix + "_InvMass", 800, 600);
    histInvMass->Draw();
    info_text->Draw();
    leg->Draw();
    canv->SaveAs("fig/" + outputPrefix + "_InvMass.png", "Q");
}
