#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

#include <TMatrix.h>
#include <TVector.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TAxis.h>

using namespace std;

void plot_alphas(string& type, string& categ, string& nbEvents,
                    int binningEta, int binningPhi, int stained)
{
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.18);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.10);
    gStyle->SetTitleOffset(1.4, "z");

    // Prefix for file input and output
    TString inputPrefix = Form("%s_%s_%s_%d_%d_%s",
        type.c_str(), categ.c_str(), nbEvents.c_str(), binningEta, binningPhi,
        stained ? "stained" : "unstained");

    // Get the data back
    string data_file_path = (inputPrefix + "_alphas.txt").Data();
    fstream stream_in(data_file_path.c_str(), ios::in);
    if (!stream_in.is_open()) {
        cout << "Bad alpha file ! EXIT" << endl;
        exit(1);
    }

    ///////////////////////
    // Plot of the results
    ///////////////////////

    // Plot diff_i = alpha_i - lambda_i and
    // diff_i_over_sigma_i (alpha_i - lambda_i) / sigma_i
    TH2F* hist_diff = new TH2F("hist_diff", "alphas",
                                binningEta, -2.4, 2.4,
                                binningPhi, -TMath::Pi(),  TMath::Pi());
    hist_diff->GetXaxis()->SetTitle("#eta");
    hist_diff->GetYaxis()->SetTitle("#phi");
    hist_diff->GetZaxis()->SetTitle("#alpha_{i}");

    double** diff_i = new double*[binningEta];
    double** diff_i_over_sigma_i = new double*[binningEta];
    double** zerr = new double*[binningEta];
    double col1, col2;
    for (int i = 0; i < binningEta; i++) {
        diff_i[i] = new double[binningPhi];
        diff_i_over_sigma_i[i] = new double[binningPhi];
        zerr[i] = new double[binningPhi];
        for (int j = 0; j < binningPhi; j++) {
            stream_in >> col1 >> col2; // Elements from col1 are the alphas, and from col2 the alphaErs.
            diff_i[i][j] = col1 -
                (stained ? ( (i+j) % 2 == 0 ? 0.01 : -0.01) : 0);
            zerr[i][j] = col2;
            diff_i_over_sigma_i[i][j] = diff_i[i][j] / zerr[i][j];

            hist_diff->SetBinContent(i+1, j+1, col1); // ROOT convention for binning
        }
    }


    TPaveText* info_text = new TPaveText(0.5, 0.73, 0.75, 0.83, "ndc");
    info_text->SetBorderSize(0);
    info_text->SetTextSize(0.04);
    info_text->SetFillColorAlpha(kWhite, 0.5);
    info_text->AddText(Form("%s %s %s, %s", type.c_str(), categ.c_str(),
                nbEvents.c_str(), (stained ? "stained" : "unstained")));
    info_text->AddText(Form("Nb bins #eta x #phi: %dx%d",
                            binningEta, binningPhi));

    TCanvas* canv_diff = new TCanvas("canv_diff",
        inputPrefix + "_diff_i", 800, 600);
    hist_diff->Draw("COLZ");
//    hist_diff->Draw("LEGO2");
    info_text->Draw();
    canv_diff->SaveAs("fig/" + inputPrefix + "_diff_i.png", "Q");
}
