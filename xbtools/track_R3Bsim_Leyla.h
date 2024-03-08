//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 22 18:30:47 2016 by ROOT version 5.34/18
// from TTree qfs/R3BRoot into land02 tree
// found on file: 12C_ppn_11C_addgammas_exs0_newr3bjes_nofra_g4inclxx_r3bsim_cutE10keV_land02_leila_resx1_s296.root
//////////////////////////////////////////////////////////

#ifndef track_R3Bsim_Leyla_h
#define track_R3Bsim_Leyla_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class track {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          Xbmul;
   UInt_t          Xbi[182];   //[Xbmul]
   Float_t         Xbt[182];   //[Xbmul]
   Float_t         Xbe[182];   //[Xbmul]
   Float_t         Xbpe[182];   //[Xbmul]

   // List of branches
   TBranch        *b_Xbmul;   //!
   TBranch        *b_Xbi;   //!
   TBranch        *b_Xbt;   //!
   TBranch        *b_Xbe;   //!
   TBranch        *b_Xbpe;   //!

	 track(TTree *tree=0);
	 virtual ~track();
	 virtual void     AddChain(TTree* tree);
	 virtual void     Analyse(Char_t *out_file);
	 virtual Int_t    Cut(Long64_t entry);
	 virtual Int_t    GetEntry(Long64_t entry);
	 virtual Long64_t LoadTree(Long64_t entry);
	 virtual void     Init(TTree *tree);
	 //virtual void     Loop();
	 virtual Bool_t   Notify();
	 virtual void     Show(Long64_t entry = -1);
};

#endif

