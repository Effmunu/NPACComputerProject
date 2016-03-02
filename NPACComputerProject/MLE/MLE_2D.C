#include "Fcn_MLE_2D.h"
#include <TMinuit.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TF2.h>
#include <TH2F.h>

double logL(double mu, double sigma)
{
    double value = 0; // initialization
    for(int i=0; i<N; i++) {
        value += -2 * log(gaussian(data[i], mu, sigma));
    }
    return value;
}   

void MLE_2D()
{
    // root style
    gStyle->SetPalette(1);  
    gStyle->SetOptTitle(0);

    //load and compile fcn
    gROOT->ProcessLine(".L Fcn_MLE_2D.C");

    //plot -2ln(L) as a function of mu and sigma.
    TF2* tf2_logL = new TF2("tf2_logL", "logL(x, y)", -1, 0, 2.0, 2.4); 
    double muEstimated = 0;
    double sigmaEstimated = 0;
    tf2_logL->GetMinimumXY(muEstimated, sigmaEstimated);
    // Transform function into histogram for 2D plot
    TH2F* hist = (TH2F*)tf2_logL->GetHistogram();

    TCanvas* canv = new TCanvas("canv", "logL(mu)", 0, 0, 800, 600);

    TGraph2D* gr = new TGraph2D(hist);
    TPaveText* text = new TPaveText(0.4, 0.7, 0.6, 0.85, "ndc");
    text->AddText(Form("mu estimated: %f", muEstimated));
    text->AddText(Form("sigma estimated: %f", sigmaEstimated));
    text->SetFillColor(kWhite);

    gr->Draw("surf1"); // Plot 2D
    text->Draw();

    canv->SaveAs("logL(mu, sigma).png");


    //nb of parameters
    const int npar = 2;


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
    par[1]= 2;
    stepSize[1] = 0.0001;   
    minVal[1] = 0;
    maxVal[1] = +5;
    parName[1] = "sigma";


    for (int i=0; i<npar; i++) {
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
    double muEstimator=0;
    double muEstimator_er=0;
    double sigmaEstimator=0;
    double sigmaEstimator_er=0;

    minuit.GetParameter(0, muEstimator, muEstimator_er);
    minuit.GetParameter(1, sigmaEstimator, sigmaEstimator_er);
}

