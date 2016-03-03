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
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetTitleOffset(1.4, "y");

    // Prefix for file input and output
    TString inputPrefix = Form("%s_%s_%s_%d_%s",
        type.c_str(), categ.c_str(), nbEvents.c_str(), binning,
        stained ? "stained" : "unstained");

    // Get the data back
    string data_file_path = (inputPrefix + "_betas.txt").Data();
    fstream stream_in(data_file_path.c_str(), ios::in);
    n_betas = binning * (binning+1) / 2;
    double* rawBetas = new double[n_betas];
    double* rawBetaErs = new double[n_betas];

    double col1, col2;
    for (int i = 0; i < n_betas; i++) {
        stream_in >> rawBetas[i] >> rawBetaErs[i];
    }

//    const int n_alpha = (-1 + sqrt(1. + 8*rawBetas.size())) / 2;

    ////////////////
    // Computation
    ////////////////

    //Initialize vectors and matrix
    TVectorT<double> alphas(binning);
    TVectorT<double> alphaErs(binning);
    TVectorT<double> betas(binning);
    TMatrixT<double> alphasToBetas(binning, binning);

    // betas vector:
    for (int i=0; i<binning; i++) {
        double b_i = 0;
        for (int j=0; j<binning; j++) {
            b_i += rawBetas[getLinearIndex(binning, i, j)] / rawBetaErs[getLinearIndex(binning, i, j)] / rawBetaErs[getLinearIndex(binning, i, j)];
        }
        betas[i] = b_i;
    }

    // alphasToBetas matrix:
    for (int i=0; i<binning; i++) {
        for (int j=0; j<binning; j++) {
            double u_ij = 0;
            if (fabs(rawBetas[getLinearIndex(binning, i, j)]) > 0) { // If beta_ij is non zero
                if (i != j) {
                    alphasToBetas[i][j] = 1. / rawBetaErs[getLinearIndex(binning, i, j)] / rawBetaErs[getLinearIndex(binning, i, j)];
                }
                else {
                    for (int k=0; k<binning; k++) {
                        u_ij += 1. / rawBetaErs[getLinearIndex(binning, i, k)] / rawBetaErs[getLinearIndex(binning, i, k)];
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

    for (int i = 0; i < binning; i++)
        alphaErs[i] = sqrt(betasToAlphas[i][i]);

/*    for (int i = 0; i < binning; i++) {
        cout << "alphas_" << i << ": " << alphas[i]
         << "\t" << "alphaErs_" << i << ": " << alphaErs[i]
         << "\t betasToAlphas[i][i]: " <<  betasToAlphas[i][i] << endl;
    }*/

    ///////////////////////
    // Plot of the results
    ///////////////////////
    // Plot diff_i = alpha_i - lambda_i and diff_i_over_sigma_i (alpha_i - lambda_i) / sigma_i
    double* x = new double[binning];
    double* xerr = new double[binning];
    double* yerr = new double[binning];
    double* diff_i = new double[binning];
    double* diff_i_over_sigma_i = new double[binning];
    for (int i = 0; i < binning; i++) {
        x[i] = i; // eta bins indices
        xerr[i] = 0.5; // bin number +/- 0.5
        yerr[i] = alphaErs[i];
        diff_i[i] = alphas[i] -
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
