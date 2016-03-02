#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

#include <TMatrix.h>
#include <TVector.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPaveText.h>

using namespace std;

void plot_alphas(string& type, string& categ, string& nbEvents, int binning, int stained)
{
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);

    // Prefix for file input and output
    TString inputPrefix = Form("%s_%s_%s_%d_%s",
        type.c_str(), categ.c_str(), nbEvents.c_str(), binning,
        stained ? "stained" : "unstained");

    // Get the data back
    string data_file_path = (inputPrefix + "_alphas.txt").Data();
    fstream stream_in(data_file_path.c_str(), ios::in);
    vector<double> alphas;
    vector<double> alphaErs;

    double col1, col2;
    while(stream_in >> col1 >> col2) {
        alphas.push_back(col1);
        alphaErs.push_back(col2);
//        cout << col1 << " " << col2 << endl;
    }

    ///////////////////////
    // Plot of the results
    ///////////////////////

    // Plot diff_i = alpha_i - lambda_i and diff_i_over_sigma_i (alpha_i - lambda_i) / sigma_i
    const int n_alpha = alphas.size();
    double* x = new double[n_alpha];
    double* diff_i = new double[n_alpha];
    double* diff_i_over_sigma_i = new double[n_alpha];
    for (int i = 0; i < n_alpha; i++) {
        x[i] = i; // eta bins indices
        diff_i[i] = alphas[i] -
            (stained ? (i % 2 == 0 ? 0.01 : -0.01) : 0);
        diff_i_over_sigma_i[i] = diff_i[i] / alphaErs[i];
        //cout << "x_i: " << x[i] << "; diff_i: " << diff_i[i] << endl;
    }

    TPaveText* text_diff = new TPaveText(0.6, 0.7, 0.85, 0.85, "ndc");
    text_diff->SetBorderSize(0);
    text_diff->SetTextSize(0.04);
    text_diff->SetFillColor(kWhite);
    text_diff->AddText(stained ? "#alpha_{i} - #lambda_{i}" : "#alpha_{i}");

    TPaveText* text_diff_over = new TPaveText(0.6, 0.7, 0.85, 0.85, "ndc");
    text_diff_over->SetBorderSize(0);
    text_diff_over->SetTextSize(0.04);
    text_diff_over->SetFillColor(kWhite);
    text_diff_over->AddText(stained ? "#frac{#alpha_{i} - #lambda_{i}}{#sigma_{i}}" : "#frac{#alpha_{i}}{#sigma_{i}}");

    TGraph* graph_diff = new TGraph(n_alpha, x, diff_i);
    graph_diff->GetXaxis()->SetTitle("Eta bin Number");
    graph_diff->GetYaxis()->SetTitle("alpha_i");
    graph_diff->SetMaximum(graph_diff->GetYaxis()->GetXmax()*1.3);

    TCanvas* canv_diff = new TCanvas("canv_diff",
        inputPrefix + "_diff_i", 800, 600);
    graph_diff->Draw("AB");
    text_diff->Draw();
    canv_diff->SaveAs("fig/" + inputPrefix + "_diff_i.png", "Q");

    TGraph* graph_diff_over = new TGraph(n_alpha, x, diff_i_over_sigma_i);
    graph_diff_over->GetXaxis()->SetTitle("Eta bin Number");
    graph_diff_over->GetYaxis()->SetTitle("alpha_i");
    graph_diff_over->SetMaximum(graph_diff_over->GetYaxis()->GetXmax()*1.3);
    TCanvas* canv_diff_over = new TCanvas("canv_diff_over",
        inputPrefix + "_diff_i_over_sigma_i", 800, 600);
    graph_diff_over->Draw("AB");
    text_diff_over->Draw();
    canv_diff_over->SaveAs("fig/" + inputPrefix + "_diff_i_over_sigma_i.png", "Q");
}
