#include "src/MappingTool.cpp"

#include <fstream>
#include <cstdlib>
#include <iostream>

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

    TString inputPrefix = Form("data_%s_%s", categ.c_str(), nbEvents.c_str()); // Z_50000
    TFile* inputFile = TFile::Open((inputPrefix + ".root").Data());
    TTree* inputTree = (TTree*)inputFile->Get("tuple");
    Long64_t nbEntries = inputTree->GetEntries();
    cout << "Tree loaded." << endl;

    TString inputAplhas = inputPrefix + Form("_%d_unstained_alphas.txt", binning);
    cout << "Will read alphas from file: " << inputAplhas.Data() << endl;
    fstream stream_in( inputAplhas.Data(), ios::in);

    double* alphas = new double[binning];
    double dummy;
    for (int i = 0; i < binning; i++) {
        stream_in >> alphas[i] >> dummy; // dummy corresponds to alphaErs, not useful here.
        cout << alphas[i] << endl;
    }
    cout << "Alphas loaded." << endl;

    MappingTool map(binning, -2.4,2.4);

    // Histogram of invariant masses
    TH1F* histInvMass = new TH1F("histInvMass",
        "Invariant Mass " + inputPrefix, 100, 80, 100);
    histInvMass->SetLineColor(kBlue);

    TH1F* histInvMass_cor = new TH1F("histInvMass_cor",
        "Corrected Invariant Mass " + inputPrefix, 100, 80, 100);
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
//cout << "before" << endl;
        inputTree->GetEntry(entry);
//cout << "after" << endl;
        //========================================
        //Z events preselection
        //========================================
        if( ele0_et<20) continue;
        if( ele1_et<20) continue;
        if( fabs(ele0_eta)>2.4) continue;
        if( fabs(ele1_eta)>2.4) continue;

//        cout << "after cuts" << endl;

        // UNCORRECTED
        el0_LV->SetPtEtaPhiM(ele0_et, ele0_eta, ele0_phi, 511e-6);
        el1_LV->SetPtEtaPhiM(ele1_et, ele1_eta, ele1_phi, 511e-6);

        double mass = ((*el0_LV) + (*el1_LV)).M();

//        cout << "after uncor computation" << endl;

        // CORRECTED
        el0_LV_cor->SetPtEtaPhiM(ele0_et / (1. + alphas[map.getIndex(ele0_eta)]),
            ele0_eta, ele0_phi, 511e-6);
        el1_LV_cor->SetPtEtaPhiM(ele1_et / (1. + alphas[map.getIndex(ele1_eta)]),
            ele1_eta, ele1_phi, 511e-6);

        double mass_cor = ((*el0_LV_cor) + (*el1_LV_cor)).M();
//        cout << "after cor computaiton" << endl;

        // Cuts on invariant mass
        if(mass > 80 && mass < 100) {
            nbUncorFilled++;
            histInvMass->Fill(mass);
        }

//        cout << "after uncor fill" << endl;

        if(mass_cor > 80 && mass_cor < 100) {
            nbCorFilled++;
            histInvMass_cor->Fill(mass_cor);
        }

//        cout << "after cor fill" << endl;
    }
    cout << "Loop over events finished." << endl;
    cout << "uncor filled: " << nbUncorFilled << endl;
    cout << "cor filled: " << nbCorFilled << endl;

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
    histInvMass_cor->Fit("myVoigt"); // Q: quiet mode
    double norm = myVoigt->GetParameter(0);
    double mZ = myVoigt->GetParameter(1);
    double sigma = myVoigt->GetParameter(2);
    double gamma = myVoigt->GetParameter(3);

    // Display
    TCanvas* canv = new TCanvas("canv",
        inputPrefix + "_InvMass", 800, 600);

    TH1F* frame = canv->DrawFrame(80., 0., 100., 1., "frameInvMass");
    frame->GetYaxis()->SetTitle("m_{ee} [GeV]");
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
    const double ampMC = 1.75828e+03;
    const double meanMC = 9.11448e+01;
    const double sigmaMC = 8.94999e-01;
    const double gammaMC = 2.36396e+00;

    const double ampDataRaw = 8.83974e+03;
    const double meanDataRaw = 8.83741e+01;
    const double sigmaDataRaw = 2.15181e+00;
    const double gammaDataRaw = 2.33554e+00;

    cout << "Constant term before correction: "
        << sqrt(2 * (sigmaDataRaw * sigmaDataRaw / meanDataRaw / meanDataRaw
                    - sigmaMC * sigmaMC / meanMC / meanMC)) << endl;
    cout << "Constant term after correction: "
        << sqrt(2 * (sigma * sigma / mZ / mZ
                    - sigmaMC * sigmaMC / meanMC / meanMC)) << endl;

}
