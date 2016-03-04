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
#include "FitterStandard.hpp"

#define TIME_STUDY 1

void AnaCalib::Loop(string& type, string& categ, string& nbEvents, int binning, int stained, Long64_t nbEntriesToRead)
{
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.08);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetTitleOffset(1.4, "y");

    const double lowerBound = categ == "JPsi" ? 2.6 : 80;
    const double higherBound = categ == "JPsi" ? 3.6 : 100;

#if TIME_STUDY
    // Measure the running time
    const clock_t startingTime = clock();
#endif

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();

    // Prefix for file output
    TString outputPrefix = Form("%s_%s_%s_%d_%s_%lld",
        type.c_str(), categ.c_str(), nbEvents.c_str(), binning,
        stained ? "stained" : "unstained", nbEntriesToRead);

    //data vector
    vector<InfoForFitter> infoVector;

    //Instanciate the mapping tool
    MappingTool map(binning, -2.4,2.4);
    vector<double> table;
    MappingTool map2(table);

    // Histogram of invariant masses
    TH1F* histInvMass = new TH1F("histInvMass",
        "Invariant Mass " + outputPrefix, 100, lowerBound, higherBound);

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
        if (categ == "Z") {
            if( ele0_et<20) continue; // Z
            if( ele1_et<20) continue; // Z
        }
        else {
            if( ele0_et<7) continue; // JPSI
            if( ele1_et<7) continue; // JPSI
        }
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
        if(mass < lowerBound) continue; // Z
        if(mass > higherBound) continue; // Z

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
//    TF1* myVoigt = new TF1("myVoigt", "[0] * TMath::Voigt((x - [1]), [2], [3])", 80, 100); // Z
    TF1* myVoigt = new TF1("myVoigt", "[0] * TMath::Voigt((x - [1]), [2], [3])",
        lowerBound, higherBound); // JPSI

    // Initialization of the parameters
    myVoigt->SetParName(0, "Amp");
    myVoigt->SetParName(1, "Mean");
    myVoigt->SetParName(2, "Sigma");
    myVoigt->SetParName(3, "Gamma");

    // (we leave them all free to get a meaningful estimation of sigma
    // before any kind of correction)
    if (categ == "JPsi") {
        myVoigt->SetParameter(0, 2000);
        myVoigt->SetParameter(1, 3.);
        myVoigt->SetParameter(2, 1);
        myVoigt->FixParameter(3, 0); // Fit with a gaussian for J/Psi, not a GBW

        myVoigt->SetParLimits(1, 2.6, 3.6);
//        myVoigt->SetParLimits(3, 0, 1e-3);
    }
    else { // Z
        myVoigt->SetParameter(0, 5000);
        myVoigt->SetParameter(1, 91);
        myVoigt->SetParameter(2, 1);
        myVoigt->SetParameter(3, 2);
    }



    // Fit and get the amplitude and sigma of the histogram
    histInvMass->Fit("myVoigt"); // Q: quiet mode
    double norm = myVoigt->GetParameter(0);
    double mZ = categ == "JPsi" ? 3.096 : 91.1876; // PDG value Z
    double sigma = myVoigt->GetParameter(2);
    double gamma = categ == "JPsi" ? 92.9e-6 : 2.4952; // PDG value Z

    // Export fit parameters to file, to be used in other programs
    fstream fitparam_out( (outputPrefix + "_fitparam.txt").Data(), ios::out);
    fitparam_out << myVoigt->GetParameter(0) << endl;
    fitparam_out << myVoigt->GetParameter(1) << endl;
    fitparam_out << myVoigt->GetParameter(2) << endl;
    fitparam_out << myVoigt->GetParameter(3) << endl;
    fitparam_out.close();

    //initialize the fitter tool
    FitterStandard fit_standard(&map,"std");
    fit_standard.SetParameters(norm, mZ, sigma, gamma);
    // /!\ Remark: For the determination of the alphas, the mean and gamma of the distribution are fixed to the PDG values for Z!
    // The norm and sigma are fixed from the fit just above on the original distribution, leaving all parameters free to have a meaningful fit.
    fit_standard.SetData(&infoVector);
    // do the fit
    fit_standard.Execute();

    //get the results
    vector<double > alpha_standard= fit_standard.GetAlphas();
    vector<double > alphaer_standard= fit_standard.GetAlphaErs();

    // Export results to file (because it's looong to run...)
    fstream alphas_out( (outputPrefix + "_alphas.txt").Data(), ios::out);
    for (int i=0; i<alpha_standard.size(); i++) {
        alphas_out << alpha_standard[i] << " " << alphaer_standard[i] << endl;
    }
    alphas_out.close();

#if TIME_STUDY
    // End of measure of the running time
    fstream outputStream("timeResults_standard.txt", ios::app);
    outputStream << type << " " << categ << " " << nbEvents << " "
        << binning << "\t" << (stained ? "stained" : "unstained") << "\t"
        << nbEntriesToRead << "\t"
        << double(clock() - startingTime) / CLOCKS_PER_SEC << " s" << endl;
    outputStream.close();
#endif

    // Display
    TPaveText* info_text = new TPaveText(0.70, 0.52, 0.92, 0.82, "ndc");
    info_text->SetBorderSize(0);
    info_text->SetTextSize(0.04);
    info_text->SetTextAlign(32);
    info_text->SetFillColorAlpha(kWhite, 0.5);
    info_text->AddText(Form("%s %s %s, %s", type.c_str(), categ.c_str(),
                    nbEvents.c_str(), (stained ? "stained" : "unstained")));
    info_text->AddText(Form("Nb bins: %d", binning));
    info_text->AddText(Form("Nb bins %d", binning));
    info_text->AddText(Form("%s: %.2f", myVoigt->GetParName(0),
                                myVoigt->GetParameter(0)));
    info_text->AddText(Form("%s: %.2f", myVoigt->GetParName(1),
                                myVoigt->GetParameter(1)));
    info_text->AddText(Form("%s: %.2f", myVoigt->GetParName(2),
                                myVoigt->GetParameter(2)));
    info_text->AddText(Form("%s: %.2f", myVoigt->GetParName(3),
                                myVoigt->GetParameter(3)));

    TLegend* leg= new TLegend(0.15, 0.73, 0.35, 0.83);
    leg->SetFillColorAlpha(kWhite, 0.5);
    leg->SetTextSize(0.04);
    leg->AddEntry(histInvMass, "Data", "F");
    leg->AddEntry(myVoigt, "Voigtian fit", "L");

    histInvMass->GetXaxis()->SetTitle("m_{ee} [GeV]");
    histInvMass->GetYaxis()->SetTitle(categ == "JPsi" ?
        "Number of event / 0.01 GeV" : "Number of event / 0.2 GeV");
    TCanvas* canv = new TCanvas("canv",
        outputPrefix + "_InvMass", 800, 600);
    histInvMass->Draw();
    info_text->Draw();
    leg->Draw();
    canv->SaveAs("fig/" + outputPrefix + "_InvMass.png", "Q");
}
