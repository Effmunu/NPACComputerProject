#include <TH1F.h>
#include <TMath.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TFile.h>

void fit_sameHist()
{
    // Fit the truth_mass by a Breit Wigner distribution,
    // and the mass (reconstructed mass) by a Voigt function.

    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadRightMargin(0.08);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetTitleOffset(1.4, "y");

    // Retrieve histograms from file
    TFile *histFile = TFile::Open("histograms.root");

    TH1F * truth_mass = (TH1F*)histFile->Get("truth_mass"); // True Z mass
    TH1F * mass = (TH1F*)histFile->Get("mass"); // Reconstructed Z mass


    // Breit-Wigner function
    TF1 *myBW = new TF1("myBW", "[0] * ([2]*[2]*x*x/[1]/[1]) / ((x*x - [1]*[1])*(x*x - [1]*[1]) + [2]*[2]*x*x*x*x/[1]/[1])", 
        truth_mass->GetXaxis()->GetXmin(), truth_mass->GetYaxis()->GetXmax());
    // Breit-Wigner distribution convoluted with Gaussian (Voigt function)
    TF1 *myVoigt = new TF1("myVoigt", "[0] * TMath::Voigt((x - [1]), [2], [3])", 
        truth_mass->GetXaxis()->GetXmin(), truth_mass->GetYaxis()->GetXmax());

    // Name the parameters
    myBW->SetParName(0, "Amp.");
    myBW->SetParName(1, "Mean");
    myBW->SetParName(2, "Gamma");
    myVoigt->SetParName(0, "Amp.");
    myVoigt->SetParName(1, "Mean");
    myVoigt->SetParName(2, "Sigma");
    myVoigt->SetParName(3, "Gamma");

    // Set initial parameter values
    myBW->SetParameter(0, 11e6);
    myBW->SetParameter(1, 90);
    myBW->SetParameter(2, 5);
    myVoigt->SetParameter(0, 5000);
    myVoigt->SetParameter(1, 91);
    myVoigt->SetParameter(2, 0.5);
    myVoigt->SetParameter(3, 2.5);

    TCanvas* canv = new TCanvas("canv", "Zmass_fit", 800, 600);
    canv->cd();

    // Fit
    truth_mass->Fit("myBW");
    // Constrain the parameters (set limits): here, we fix Mass_Z (param. 1) and Gamma_Z (param. 3)
//    myVoigt->FixParameter(1, myBW->GetParameter(1));
//    myVoigt->FixParameter(3, myBW->GetParameter(2));
    mass->Fit("myVoigt");

    // Canvas
    TH1F* frame = canv->DrawFrame(80, 0., 100, 1., "frameInvMass");
    frame->GetXaxis()->SetTitle("m_{ee} [GeV]");
    frame->GetYaxis()->SetTitle("Number of events / 0.2 GeV");
    frame->SetMaximum(1.05 * TMath::Max(mass->GetMaximum(), 
                        truth_mass->GetMaximum()));

    truth_mass->SetLineColor(kGreen);
    truth_mass->SetFillColor(kGreen);
    myBW->SetLineColor(kBlack);
    truth_mass->GetFunction("myBW")->SetLineColor(kBlack);
    mass->SetLineColor(kBlue);
    myVoigt->SetLineColor(kRed);

    truth_mass->Draw("same");
    mass->Draw("same");
//    myBW->Draw("same");
//    myVoigt->Draw("same");

    // Legend
    TLegend* leg = new TLegend(0.70, 0.7, 0.90, 0.90);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.04);
    leg->AddEntry(truth_mass, "Truth mass", "F");
    leg->AddEntry(truth_mass->GetFunction("myBW"), "BW fit", "L");
    leg->AddEntry(mass, "Reco mass", "F");
    leg->AddEntry(mass->GetFunction("myVoigt"), "Voigtian fit", "L");
    leg->Draw();

    TPaveText* text_truth = new TPaveText(0.15, 0.75, 0.40, 0.90, "ndc");
    text_truth->SetBorderSize(0);
    text_truth->SetTextAlign(12);
    text_truth->SetTextSize(0.04);
    text_truth->SetFillColor(kWhite);
    text_truth->AddText("True Z mass fit:");
    text_truth->AddText(Form("%s: %.2f", myBW->GetParName(1),
                                myBW->GetParameter(1)));
    text_truth->AddText(Form("%s: %.2f", myBW->GetParName(2),
                                myBW->GetParameter(2)));
    text_truth->Draw();

    TPaveText* text_reco = new TPaveText(0.15, 0.50, 0.40, 0.70, "ndc");
    text_reco->SetBorderSize(0);
    text_reco->SetTextAlign(12);
    text_reco->SetTextSize(0.04);
    text_reco->SetFillColor(kWhite);
    text_reco->AddText("Reco Z mass fit:");
    text_reco->AddText(Form("%s: %.2f", myVoigt->GetParName(1),
                                myVoigt->GetParameter(1)));
    text_reco->AddText(Form("%s: %.2f", myVoigt->GetParName(2),
                                myVoigt->GetParameter(2)));
    text_reco->AddText(Form("%s: %.2f", myVoigt->GetParName(3),
                                myVoigt->GetParameter(3)));
    text_reco->Draw();

    // Save image
    canv->SaveAs("DoubleFit_Z.png");
}
