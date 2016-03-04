//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug  2 22:59:33 2011 by ROOT version 5.28/00b
// from TTree tuple/tuple
// found on file: output.root
//////////////////////////////////////////////////////////

#ifndef AnaCalib_h
#define AnaCalib_h

#include <cstdlib>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class AnaCalib {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         ele0_et;
   Float_t         ele0_eta;
   Float_t         ele0_phi;
   Float_t         ele0_charge;
   Float_t         ele1_et;
   Float_t         ele1_eta;
   Float_t         ele1_phi;
   Float_t         ele1_charge;

   // List of branches
   TBranch        *b_ele0_et;   //!
   TBranch        *b_ele0_eta;   //!
   TBranch        *b_ele0_phi;   //!
   TBranch        *b_ele0_charge;   //!
   TBranch        *b_ele1_et;   //!
   TBranch        *b_ele1_eta;   //!
   TBranch        *b_ele1_phi;   //!
   TBranch        *b_ele1_charge;   //!

   AnaCalib(TTree *tree=0);
   virtual ~AnaCalib();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string& type, std::string& categ, std::string& nbEvent, int binning=24, int stained=0, Long64_t nbEntriesToRead=10000);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnaCalib_cxx
AnaCalib::AnaCalib(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output.root");
      if (!f) {
         f = new TFile("output.root");
      }
      tree = (TTree*)gDirectory->Get("tuple");

   }
   Init(tree);
}

AnaCalib::~AnaCalib()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnaCalib::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnaCalib::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnaCalib::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ele0_et", &ele0_et, &b_ele0_et);
   fChain->SetBranchAddress("ele0_eta", &ele0_eta, &b_ele0_eta);
   fChain->SetBranchAddress("ele0_phi", &ele0_phi, &b_ele0_phi);
   fChain->SetBranchAddress("ele0_charge", &ele0_charge, &b_ele0_charge);
   fChain->SetBranchAddress("ele1_et", &ele1_et, &b_ele1_et);
   fChain->SetBranchAddress("ele1_eta", &ele1_eta, &b_ele1_eta);
   fChain->SetBranchAddress("ele1_phi", &ele1_phi, &b_ele1_phi);
   fChain->SetBranchAddress("ele1_charge", &ele1_charge, &b_ele1_charge);
   Notify();
}

Bool_t AnaCalib::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnaCalib::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnaCalib::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnaCalib_cxx
