#include "src/MappingTool.cpp"

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>

using namespace std;

void fudge(string& categ, string& nbEvents, int binning)
{
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetTitleOffset(1.4, "y");

    const double lowerBound = categ == "JPsi" ? 2.6 : 80;
    const double higherBound = categ == "JPsi" ? 3.6 : 100;

    TString inputPrefix = Form("data_%s_%s", categ.c_str(), nbEvents.c_str()); // Z_50000
    TFile* inputFile = TFile::Open((inputPrefix + ".root").Data());
    TTree* inputTree = (TTree*)inputFile->Get("tuple");
    Long64_t nbEntries = inputTree->GetEntries();
    cout << "Tree loaded." << endl;

    // read the alphas
    TString inputAplhas = inputPrefix + Form("_%d_unstained_alphas.txt", binning);
    cout << "Will read alphas from file: " << inputAplhas.Data() << endl;
    fstream stream_in( inputAplhas.Data(), ios::in);

    double* alphas = new double[binning];
    double dummy;
    for (int i = 0; i < binning; i++) {
        stream_in >> alphas[i] >> dummy; // dummy corresponds to alphaErs, not useful here.
        cout << alphas[i] << endl;
    }
    stream_in.close();
    cout << "Alphas loaded." << endl;

    // read MC fit parameters
    double amp_mc, mZ_mc, sigma_mc, gamma_mc;
    TString input_mc_fitparam_string;
    if (categ == "Z")
        input_mc_fitparam_string = Form("mc_%s_10000_%d_unstained_fitparam.txt", categ.c_str(), binning); // Z
    else
        input_mc_fitparam_string = Form("mc_%s_50000_%d_unstained_fitparam.txt", categ.c_str(), binning); // JPSI
    cout << "Will read mc fit param from file: " << input_mc_fitparam_string.Data() << endl;
    fstream input_mc_fitparam( input_mc_fitparam_string.Data(), ios::in);
    input_mc_fitparam >> amp_mc >> mZ_mc >> sigma_mc >> gamma_mc;
    input_mc_fitparam.close();
    cout << amp_mc << " " << mZ_mc << " " << sigma_mc << " " << gamma_mc << endl;

    // read data fit parameters
    double amp_data, mZ_data, sigma_data, gamma_data;
    TString input_data_fitparam_string = Form("data_%s_%s_%d_unstained_fitparam.txt", categ.c_str(), nbEvents.c_str(), binning);
    cout << "Will read data fit param from file: " << input_data_fitparam_string.Data() << endl;
    fstream input_data_fitparam( input_data_fitparam_string.Data(), ios::in);
    input_data_fitparam >> amp_data >> mZ_data >> sigma_data >> gamma_data;
    input_data_fitparam.close();
    cout << amp_data << " " << mZ_data << " " << sigma_data << " " << gamma_data << endl;

    MappingTool map(binning, -2.4,2.4);

    // Histogram of invariant masses
    TH1F* histInvMass = new TH1F("histInvMass",
        "Invariant Mass " + inputPrefix, 100, lowerBound, higherBound);
    histInvMass->SetLineColor(kBlue);

    TH1F* histInvMass_cor = new TH1F("histInvMass_cor",
        "Corrected Invariant Mass " + inputPrefix, 100, lowerBound, higherBound);
    histInvMass->SetLineColor(kRed);

    // SetBranchAddresses
    Float_t ele0_et;
    Float_t ele0_eta;
    Float_t ele0_phi;
    Float_t ele1_et;
    Float_t ele1_eta;
    Float_t ele1_phi;

    inputTree->SetBranchAddress("ele0_et", &ele0_et);
    inputTree->SetBranchAddress("ele0_eta", &ele0_eta);
    inputTree->SetBranchAddress("ele0_phi", &ele0_phi);
    inputTree->SetBranchAddress("ele1_et", &ele1_et);
    inputTree->SetBranchAddress("ele1_eta", &ele1_eta);
    inputTree->SetBranchAddress("ele1_phi", &ele1_phi);

    cout << "Branches initialized." << endl;

    //===============================================
    // Loop over the events
    //===============================================

    int nbUncorFilled = 0;
    int nbCorFilled = 0;

    TLorentzVector* el0_LV = new TLorentzVector();
    TLorentzVector* el1_LV = new TLorentzVector();
    TLorentzVector* el0_LV_cor = new TLorentzVector();
    TLorentzVector* el1_LV_cor = new TLorentzVector();

    for (Long64_t entry=0; entry<nbEntries; entry++) {
        if (entry % 1000 == 0)
            cout << "Entry: " << entry << endl;

        inputTree->GetEntry(entry);

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

        // UNCORRECTED DATA
        el0_LV->SetPtEtaPhiM(ele0_et, ele0_eta, ele0_phi, 511e-6);
        el1_LV->SetPtEtaPhiM(ele1_et, ele1_eta, ele1_phi, 511e-6);

        double mass = ((*el0_LV) + (*el1_LV)).M();

        // CORRECTED DATA (with the alphas read from file)
        el0_LV_cor->SetPtEtaPhiM(ele0_et / (1. + alphas[map.getIndex(ele0_eta)]),
            ele0_eta, ele0_phi, 511e-6);
        el1_LV_cor->SetPtEtaPhiM(ele1_et / (1. + alphas[map.getIndex(ele1_eta)]),
            ele1_eta, ele1_phi, 511e-6);

        double mass_cor = ((*el0_LV_cor) + (*el1_LV_cor)).M();

        // Cuts on invariant mass
        if(mass > lowerBound && mass < higherBound) {
            nbUncorFilled++;
            histInvMass->Fill(mass);
        }

        if(mass_cor > lowerBound && mass_cor < higherBound) {
            nbCorFilled++;
            histInvMass_cor->Fill(mass_cor);
        }

    }
    cout << "Loop over events finished." << endl;
    cout << "uncor filled: " << nbUncorFilled << endl;
    cout << "cor filled: " << nbCorFilled << endl;

    //========================================
    //Fit
    //========================================
    // Voigt fitting function
    TF1* myVoigt = new TF1("myVoigt", "[0] * TMath::Voigt((x - [1]), [2], [3])", 80, 100);
    // Initialization of the parameters
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

    // Fit
    histInvMass_cor->Fit("myVoigt");

    fstream fitparam_out( (inputPrefix + "_cor_fitparam.txt").Data(), ios::out);
    fitparam_out << myVoigt->GetParameter(0) << endl;
    fitparam_out << myVoigt->GetParameter(1) << endl;
    fitparam_out << myVoigt->GetParameter(2) << endl;
    fitparam_out << myVoigt->GetParameter(3) << endl;
    fitparam_out.close();

//    double norm = myVoigt->GetParameter(0);
    double mZ = myVoigt->GetParameter(1);
    double sigma = myVoigt->GetParameter(2);
//    double gamma = myVoigt->GetParameter(3);

    // Display
    TCanvas* canv = new TCanvas("canv",
        inputPrefix + "_InvMass", 800, 600);

//    TH1F* frame = canv->DrawFrame(80., 0., 100., 1., "frameInvMass");
    TH1F* frame = canv->DrawFrame(lowerBound, 0., higherBound, 1., "frameInvMass");
    frame->GetXaxis()->SetTitle("m_{ee} [GeV]");
    frame->GetYaxis()->SetTitle("Number of event / 0.2 GeV");

    frame->SetMaximum(1.05 * TMath::Max(histInvMass->GetMaximum(),
        histInvMass_cor->GetMaximum()));

    TPaveText* info_text = new TPaveText(0.63, 0.65, 0.88, 0.8, "ndc");
    info_text->SetBorderSize(0);
    info_text->SetTextSize(0.04);
    info_text->SetFillColor(kWhite);
    info_text->AddText(Form("data %s %s", categ.c_str(), nbEvents.c_str()));
    info_text->AddText(Form("Nb bins: %d", binning));

    TLegend* leg = new TLegend(0.20, 0.65, 0.45, 0.8);
    leg->SetTextSize(0.04);
    leg->AddEntry(histInvMass, "Data", "F");
    leg->AddEntry(histInvMass_cor, "Cor. Data", "F");
    leg->AddEntry(myVoigt, "Voigtian fit on Cor. Data", "L");

    histInvMass->Draw("same");
    histInvMass_cor->Draw("same");
    info_text->Draw();
    leg->Draw();
    canv->SaveAs("fig/" + inputPrefix + "_InvMassCor.png", "Q");

    // Print c factors before and after
    cout << sigma_data << " " << mZ_data << " " << sigma_mc << " " << mZ_mc << " " << sigma << " " << mZ << endl;

    cout << "Constant term before correction: "
        << sqrt(2 * (sigma_data * sigma_data / mZ_data / mZ_data
                    - sigma_mc * sigma_mc / mZ_mc / mZ_mc)) << endl;
    cout << "Constant term after correction: "
        << sqrt(2 * (sigma * sigma / mZ / mZ
                    - sigma_mc * sigma_mc / mZ_mc / mZ_mc)) << endl;

}
