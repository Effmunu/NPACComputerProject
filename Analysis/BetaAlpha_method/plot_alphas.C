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

int getLinearIndex(int nbin, int i, int j)
{
    // We only use the triangular superior part of the beta_ij matrix.
    // So here is the calculation of the linear index of this triangular sup part from (i,j).
    if (i<=j)
        return nbin*(nbin-1)/2 - (nbin-i)*(nbin-i-1)/2 + j; // normal case
    else
        return nbin*(nbin-1)/2 - (nbin-j)*(nbin-j-1)/2 + i; // i,j inverted to get triangular sup.
}

int getMatrixI(int nbin, int k)
{
    // We only use the triangular superior part of the beta_ij matrix.
    // So here is a way to get the first matrix index from the linear index of this triangular sup part.
    int i = 0;
    while(k+1 > (i+1)*nbin - i*(i+1)/2)
        i++;
    return i;
}

int getMatrixJ(int nbin, int k)
{
    // We only use the triangular superior part of the beta_ij matrix.
    // So here is a way to get the second matrix index from the linear index of this triangular sup part.

    // Just invert the linear index formula by taking the i index with GetMatrixI
    int i = getMatrixI(nbin, k);
    return k - nbin*(nbin-1)/2 + (nbin-i)*(nbin-i-1)/2;
}

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
    string data_file_path = (inputPrefix + "_betas.txt").Data();
    fstream stream_in(data_file_path.c_str(), ios::in);
    vector<double> rawBetas;
    vector<double> rawBetaErs;

    double col1, col2;
    while(stream_in >> col1 >> col2) {
        rawBetas.push_back(col1);
        rawBetaErs.push_back(col2);
//        cout << col1 << " " << col2 << endl;
    }
    const int n_alpha = (-1 + sqrt(1. + 8*rawBetas.size())) / 2;

    ////////////////
    // Computation
    ////////////////

    //Initialize vectors and matrix
    TVectorT<double> alphas(n_alpha);
    TVectorT<double> alphaErs(n_alpha);
    TVectorT<double> betas(n_alpha);
    TMatrixT<double> alphasToBetas(n_alpha, n_alpha);

    // betas vector:
    for (int i=0; i<n_alpha; i++) {
        double b_i = 0;
        for (int j=0; j<n_alpha; j++) {
            b_i += rawBetas[getLinearIndex(n_alpha, i, j)] / rawBetaErs[getLinearIndex(n_alpha, i, j)] / rawBetaErs[getLinearIndex(n_alpha, i, j)];
        }
        betas[i] = b_i;
    }

    // alphasToBetas matrix:
    for (int i=0; i<n_alpha; i++) {
        for (int j=0; j<n_alpha; j++) {
            double u_ij = 0;
            if (fabs(rawBetas[getLinearIndex(n_alpha, i, j)]) > 0) { // If beta_ij is non zero
                if (i != j) {
                    alphasToBetas[i][j] = 1. / rawBetaErs[getLinearIndex(n_alpha, i, j)] / rawBetaErs[getLinearIndex(n_alpha, i, j)];
                }
                else {
                    for (int k=0; k<n_alpha; k++) {
                        u_ij += 1. / rawBetaErs[getLinearIndex(n_alpha, i, k)] / rawBetaErs[getLinearIndex(n_alpha, i, k)];
                    } // end for k
                    alphasToBetas[i][j] = u_ij;
                }
            } // end if beta_ij is non zero
            else { // if beta_ij is 0
                alphasToBetas[i][j] = 0;
            } // end else (beta_ij is 0)
        } // end for j
    } //end for i
    // alphasToBetas matrix initialized

    // Computation of alphas and alphaErs vectors
    TMatrixT<double> betasToAlphas = alphasToBetas.Invert();

    alphas = betasToAlphas * betas;

    for (int i = 0; i < n_alpha; i++)
        alphaErs[i] = sqrt(betasToAlphas[i][i]);

/*    for (int i = 0; i < n_alpha; i++) {
        cout << "alphas_" << i << ": " << alphas[i]
         << "\t" << "alphaErs_" << i << ": " << alphaErs[i]
         << "\t betasToAlphas[i][i]: " <<  betasToAlphas[i][i] << endl;
    }*/

    ///////////////////////
    // Plot of the results
    ///////////////////////
    // Plot diff_i = alpha_i - lambda_i and diff_i_over_sigma_i (alpha_i - lambda_i) / sigma_i
    double* x = new double[n_alpha];
    double* diff_i = new double[n_alpha];
    double* diff_i_over_sigma_i = new double[n_alpha];
    for (int i = 0; i < n_alpha; i++) {
        x[i] = i; // eta bins indices
        diff_i[i] = alphas[i] - (i % 2 == 0 ? 0.01 : -0.01);
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
    graph_diff->SetMaximum(graph_diff->GetMaximum()*1.3);
    TCanvas* canv_diff = new TCanvas("canv_diff",
        inputPrefix + "_diff_i", 800, 600);
    graph_diff->Draw("AB");
    text_diff->Draw();
    canv_diff->SaveAs("fig/" + inputPrefix + "_diff_i.png", "Q");

    TGraph* graph_diff_over = new TGraph(n_alpha, x, diff_i_over_sigma_i);
    graph_diff_over->GetXaxis()->SetTitle("Eta bin Number");
    graph_diff_over->GetYaxis()->SetTitle("alpha_i");
    graph_diff_over->SetMaximum(graph_diff_over->GetMaximum()*1.3);
    TCanvas* canv_diff_over = new TCanvas("canv_diff_over",
        inputPrefix + "_diff_i_over_sigma_i", 800, 600);
    graph_diff_over->Draw("AB");
    text_diff_over->Draw();
    canv_diff_over->SaveAs("fig/" + inputPrefix + "_diff_i_over_sigma_i.png", "Q");
}
