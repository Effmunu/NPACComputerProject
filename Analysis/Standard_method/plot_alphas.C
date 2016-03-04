#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

#include <TMatrix.h>
#include <TVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLine.h>

using namespace std;

void plot_alphas(string& type, string& categ, string& nbEvents, int binning, int stained, Long64_t nbEntriesToRead)
{
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetTitleOffset(1.4, "y");

    // Prefix for file input and output
    TString inputPrefix = Form("%s_%s_%s_%d_%s_%lld",
        type.c_str(), categ.c_str(), nbEvents.c_str(), binning,
        stained ? "stained" : "unstained", nbEntriesToRead);

    // Get the data back
    string data_file_path = (inputPrefix + "_alphas.txt").Data();
    fstream stream_in(data_file_path.c_str(), ios::in);

    ///////////////////////
    // Plot of the results
    ///////////////////////

    // Plot diff_i = alpha_i - lambda_i and diff_i_over_sigma_i (alpha_i - lambda_i) / sigma_i
    double* x = new double[binning];
    double* xerr = new double[binning];
    double* yerr = new double[binning];
    double* diff_i = new double[binning];
    double* diff_i_over_sigma_i = new double[binning];
    double col1, col2;
    for (int i = 0; i < binning; i++) {
        stream_in >> col1 >> col2; // Elements from col1 are the alphas, and from col2 the alphaErs.
        x[i] = i; // eta bins indices
        xerr[i] = 0.5; // bin number +/- 0.5
        yerr[i] = col2;
        diff_i[i] = col1 -
            (stained ? (i % 2 == 0 ? 0.01 : -0.01) : 0);
        diff_i_over_sigma_i[i] = diff_i[i] / yerr[i];
        //cout << "x_i: " << x[i] << "; diff_i: " << diff_i[i] << endl;
    }

    TPaveText* info_text = new TPaveText(0.6, 0.73, 0.85, 0.88, "ndc");
    info_text->SetBorderSize(0);
    info_text->SetTextSize(0.04);
    info_text->SetFillColor(kWhite);
    info_text->AddText(Form("%s %s %s, %s", type.c_str(), categ.c_str(), nbEvents.c_str(), (stained ? "stained" : "unstained")));
    info_text->AddText(Form("Nb bins: %d", binning));

    TGraphErrors* graph_diff = new TGraphErrors(binning, x, diff_i, xerr, yerr);
    graph_diff->GetXaxis()->SetTitle("#eta bin Number");
    graph_diff->GetYaxis()->SetTitle(stained ? "#alpha_{i} - #lambda_{i}" : "#alpha_{i}");
    graph_diff->SetMaximum(graph_diff->GetYaxis()->GetXmax()*1.3);

    TLine* line_diff = new TLine(-0.5, 0, binning - 0.5, 0);

    TCanvas* canv_diff = new TCanvas("canv_diff",
        inputPrefix + "_diff_i", 800, 600);
    graph_diff->Draw("AP");
    info_text->Draw();
    line_diff->Draw();
    canv_diff->SaveAs("fig/" + inputPrefix + "_diff_i.png", "Q");

    TGraph* graph_diff_over = new TGraph(binning, x, diff_i_over_sigma_i);
    graph_diff_over->GetXaxis()->SetTitle("#eta bin Number");
    graph_diff_over->GetYaxis()->SetTitle(stained ? "(#alpha_{i} - #lambda_{i}) / #sigma_{i}" : "#alpha_{i} / #sigma_{i}");
    graph_diff_over->SetMaximum(graph_diff_over->GetYaxis()->GetXmax()*1.3);

    TCanvas* canv_diff_over = new TCanvas("canv_diff_over",
        inputPrefix + "_diff_i_over_sigma_i", 800, 600);
    graph_diff_over->Draw("AB");
    info_text->Draw();
    canv_diff_over->SaveAs("fig/" + inputPrefix + "_diff_i_over_sigma_i.png", "Q");
}
