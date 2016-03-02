#include <TH1F.h>
#include <TMath.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TFile.h>

void fit()
{
    // Fit the truth_mass by a Breit Wigner distribution,
    // and the mass (reconstructed mass) by a Voigt function.

    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);

    // Retrieve histograms from file
    TFile *histFile = TFile::Open("histograms.root");

    TH1F * truth_mass = (TH1F*)histFile->Get("truth_mass"); // True Z mass
    truth_mass->GetXaxis()->SetTitle("m_{ee} [GeV]");
    truth_mass->GetYaxis()->SetTitle("Number of events");
    TH1F * mass = (TH1F*)histFile->Get("mass"); // Reconstructed Z mass
    mass->GetXaxis()->SetTitle("m_{ee} [GeV]");
    mass->GetYaxis()->SetTitle("Number of events");

    // Breit-Wigner function
    TF1 *myBW = new TF1("myBW", "[0] * ([2]*[2]*x*x/[1]/[1]) / ((x*x - [1]*[1])*(x*x - [1]*[1]) + [2]*[2]*x*x*x*x/[1]/[1])", truth_mass->GetXaxis()->GetXmin(), truth_mass->GetYaxis()->GetXmax());
    // Breit-Wigner distribution convoluted with Gaussian (Voigt function)
    TF1 *myVoigt = new TF1("myVoigt", "[0] * TMath::Voigt((x - [1]), [2], [3])", truth_mass->GetXaxis()->GetXmin(),
    truth_mass->GetYaxis()->GetXmax());

    // Name the parameters
    myBW->SetParName(0, "Amplitude");
    myBW->SetParName(1, "Mean");
    myBW->SetParName(2, "Width");
    myVoigt->SetParName(0, "Amplitude");
    myVoigt->SetParName(1, "Mean");
    myVoigt->SetParName(2, "Gaussian Sigma");
    myVoigt->SetParName(3, "Lorentzian Width");

    // Set initial parameter values
    myBW->SetParameter(0, 11e6);
    myBW->SetParameter(1, 90);
    myBW->SetParameter(2, 5);
    myVoigt->SetParameter(0, 5000);
//    myVoigt->SetParameter(1, 91); // Fixed later
    myVoigt->SetParameter(2, 0.5);
//    myVoigt->SetParameter(3, 2.5); // Fixed later

    // Canvas
    TCanvas *canv = new TCanvas("canv", "Fit of the true Z mass", 1600, 600);
    canv->Divide(2, 1); // 1 row, 2 columns

    // Plot to first pad
    canv->cd(1);
    // Fit
    truth_mass->Fit("myBW");
    // Legend
    TLegend* leg_truth = new TLegend(0.20, 0.7, 0.48, 0.85);
    leg_truth->SetTextSize(0.04);
    leg_truth->AddEntry(truth_mass, "Data", "F");
    leg_truth->AddEntry(myBW, "Fit", "L");
    leg_truth->Draw();

    TPaveText* text_truth = new TPaveText(0.6, 0.7, 0.9, 0.85, "ndc");
    text_truth->SetBorderSize(0);
    text_truth->SetTextSize(0.04);
    text_truth->SetFillColor(kWhite);
    text_truth->AddText("Fit of the true Z mass");
    text_truth->Draw();

    // Constrain the parameters (set limits): here, we fix Mass_Z (param. 1) and Gamma_Z (param. 3)
    myVoigt->FixParameter(1, myBW->GetParameter(1));
    myVoigt->FixParameter(3, myBW->GetParameter(2));

    // Plot to second pad
    canv->cd(2);
    // Fit
    mass->Fit("myVoigt");
    // Legend
    TLegend* leg_reco = new TLegend(0.20, 0.7, 0.48, 0.85);
    leg_reco->SetTextSize(0.04);
    leg_reco->AddEntry(truth_mass, "Data", "F");
    leg_reco->AddEntry(myBW, "Fit", "L");
    leg_reco->Draw();

    TPaveText* text_reco = new TPaveText(0.6, 0.7, 0.85, 0.85, "ndc");
    text_reco->SetBorderSize(0);
    text_reco->SetTextSize(0.04);
    text_reco->SetFillColor(kWhite);
    text_reco->AddText("Fit of the reconstructed");
    text_reco->AddText("Z mass");
    text_reco->Draw();

    // Save image
    canv->SaveAs("Fit of the true Z mass.png");
}
