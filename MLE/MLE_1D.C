#include "Fcn_MLE_1D.h"
#include <TMinuit.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include <TAxis.h>

double logL(double mu)
{
    double sigma = 1;
    double value = 0; // initialization
    for(int i=0; i<N; i++) {
        value += -2 * log(gaussian(data[i], mu, sigma));
    }
    return value;
}   

void MLE_1D()
{
    // root style
    gStyle->SetPalette(1);  
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetTitleOffset(1.4, "y");

    //load and compile fcn
    gROOT->ProcessLine(".L Fcn_MLE_1D.C");


    //plot -2ln(L) as a function of mu.
    TF1* tf1_logL = new TF1("tf1_logL", "logL(x)", -5, 5); 
    double muEstimator = tf1_logL->GetMinimumX(-5, 5);

    TCanvas* canv = new TCanvas("canv", "logL(mu)", 800, 600);

    TGraph* gr = new TGraph(tf1_logL);
    gr->GetXaxis()->SetTitle("mu");
    gr->GetYaxis()->SetTitle("logL(mu)");
    TPaveText* text = new TPaveText(0.4, 0.8, 0.6, 0.85, "ndc");
    text->AddText(Form("mu estimator: %f", muEstimator));
    text->SetFillColor(kWhite);

    gr->Draw();
    text->Draw();

    canv->SaveAs("logL(mu).png");

    //nb of parameters
    const int npar = 1;


    //==========================================
    //TMinuit settings
    //==========================================
    TMinuit minuit(npar);
    //  minuit.SetPrintLevel(-1);    // no print
    minuit.SetFCN(fcn_to_minimize);  //set the fcn to minimize 


    double par[npar];               // the start values
    double stepSize[npar];          // step sizes 
    double minVal[npar];            // minimum bound on parameter 
    double maxVal[npar];            // maximum bound on parameter
    string parName[npar];           // parameter name


    par[0]= 0;
    stepSize[0] = 0.0001;   
    minVal[0] = -1;
    maxVal[0] = +1;
    parName[0] = "mu";


    for (int i=0; i<npar; i++)
    {
        minuit.DefineParameter(i, parName[i].c_str(), 
                par[i], stepSize[i], minVal[i], maxVal[i]);
    }

    //==========================================
    //Do the minimization
    //==========================================
    minuit.Migrad();

    //==========================================
    //Get results
    //==========================================
    double alpha=0;
    double alphaer=0;
    minuit.GetParameter(0,alpha,alphaer);
}

