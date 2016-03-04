//===============================================
//
//             Master NPAC
// Main program for the calibration project
//
//===============================================

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

#include "TROOT.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TChain.h"

#include "AnaCalib.hpp"




int main(int argc, char ** argv)
{
    // exemple call ./bin/main mc Z 10000 24 1
    // to use mc sample with 10000 Z-events and binning 24, stained energy

    //=============================================
    //parameters
    //=============================================

    //command line option
    // default
    string type = "mc";
    string categ = "Z";
    string nbEvents = "10000";
    string inputfilename  =  "input.list";
    string outputfilename =  "output.root";
    int binning = 24;
    int stained = 0;
    Long64_t nbEntriesToRead = std::atoi(nbEvents.c_str());

    if (argc > 5) {
        type = argv[1];
        categ = argv[2];
        nbEvents = argv[3];
        binning = atoi(argv[4]);
        stained = atoi(argv[5]);

        inputfilename = "input_" + type + "_" + categ + "_" + nbEvents + ".list";
        outputfilename = "output/output_" + type + "_" + categ + "_" + nbEvents
        + "_" + std::to_string(binning) + "_"
        + (stained ? "stained" : "unstained") + ".root";
    }
    if (argc > 6) {
        nbEntriesToRead = std::atoi(argv[6]);
        outputfilename = "output/output_" + type + "_" + categ + "_" + nbEvents
        + "_" + std::to_string(binning) + "_"
        + (stained ? "stained" : "unstained") + std::to_string(nbEntriesToRead)
        + ".root";
    }

    cout << "Will read file:  " << inputfilename << endl;
    cout << "Will write file: " << outputfilename << endl;
    cout << "Will read " << nbEntriesToRead << " entries" << endl;

    // tree name
    string TreeName       =  "tuple";


    //=============================================
    // some root voodoo
    //=============================================
    TROOT test("test","test of PlotAll");
    test.SetBatch();
    gROOT->SetStyle("Plain");
    TApplication App(argv[0],&argc,argv);

    //=============================================
    //loop over the input files and make a TChain
    //=============================================
    TChain *data   =new TChain(TreeName.c_str());
    cout << ">>>Load Chain from file: " << inputfilename << endl;
    std::ifstream fList(inputfilename.c_str());
    if (!fList){
        cerr << ">>Can't open file " << inputfilename << endl;
        return 1;
    }

    char lineFromFile[10000];
    while (fList.getline(lineFromFile,10000)){
        TString fileName = lineFromFile;
        if(data->Add(fileName.Data())) cout << ">>File '"<< fileName<< "' has been loaded" <<endl;
        else cerr << ">>Can't load file '" << fileName << "'" << endl;

    }
    cout << ">>Total number of entries: " << data->GetEntries() << endl<<endl;
    fList.close();

    //=============================================
    //open the output file
    //=============================================
    TFile * file = new TFile(outputfilename.c_str(),"RECREATE");


    //=============================================
    //do the jobs
    //=============================================
    AnaCalib ana(data);
    ana.Loop(type, categ, nbEvents, binning, stained, nbEntriesToRead);

    //=============================================
    //Write in file
    //=============================================
    file->cd();
    file->Write();
    file->Close();
    delete file;

    return 0;
}
