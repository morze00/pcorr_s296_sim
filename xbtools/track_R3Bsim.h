//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 11 11:57:13 2016 by ROOT version 5.34/18
// from TTree h509/New R3BROOT tree before add-back
// found on file: ppn_gs_TGeant4_addback.root
//////////////////////////////////////////////////////////

#ifndef track_r3bsim_h
#define track_r3bsim_h

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
   Float_t         X0;
   Float_t         Y0;
   Float_t         T0;
   Float_t         fra_beta;
   Float_t         P_p[2];
   Float_t         P_e[2];
   Float_t         P_theta[2];
   Float_t         P_phi[2];
   UInt_t          Xbn;
   UInt_t          Xbi[100];   //[Xbn]
   Float_t         Xbe[100];   //[Xbn]
   Float_t         Xbhe[100];   //[Xbn]
   Float_t         Xbt[100];   //[Xbn]
   UInt_t          Ss01mul;
   Float_t         Ss01_x[64];   //[Ss01mul]
   Float_t         Ss01_y[64];   //[Ss01mul]
   Float_t         Ss01_z[64];   //[Ss01mul]
   Float_t         Ss01_e[64];   //[Ss01mul]
   Float_t         Ss01_theta[64];   //[Ss01mul]
   Float_t         Ss01_phi[64];   //[Ss01mul]
   UInt_t          Ss02mul;
   Float_t         Ss02_x[14];   //[Ss02mul]
   Float_t         Ss02_y[14];   //[Ss02mul]
   Float_t         Ss02_z[14];   //[Ss02mul]
   Float_t         Ss02_e[14];   //[Ss02mul]
   Float_t         Ss02_theta[14];   //[Ss02mul]
   Float_t         Ss02_phi[14];   //[Ss02mul]
   UInt_t          Ss03mul;
   Float_t         Ss03_x[46];   //[Ss03mul]
   Float_t         Ss03_y[46];   //[Ss03mul]
   Float_t         Ss03_z[46];   //[Ss03mul]
   Float_t         Ss03_e[46];   //[Ss03mul]
   Float_t         Ss03_theta[46];   //[Ss03mul]
   Float_t         Ss03_phi[46];   //[Ss03mul]
   UInt_t          Ss04mul;
   Float_t         Ss04_x[62];   //[Ss04mul]
   Float_t         Ss04_y[62];   //[Ss04mul]
   Float_t         Ss04_z[62];   //[Ss04mul]
   Float_t         Ss04_e[62];   //[Ss04mul]
   Float_t         Ss04_theta[62];   //[Ss04mul]
   Float_t         Ss04_phi[62];   //[Ss04mul]
   UInt_t          Ss05mul;
   Float_t         Ss05_x[64];   //[Ss05mul]
   Float_t         Ss05_y[64];   //[Ss05mul]
   Float_t         Ss05_z[64];   //[Ss05mul]
   Float_t         Ss05_e[64];   //[Ss05mul]
   Float_t         Ss05_theta[64];   //[Ss05mul]
   Float_t         Ss05_phi[64];   //[Ss05mul]
   UInt_t          Ss06mul;
   Float_t         Ss06_x[12];   //[Ss06mul]
   Float_t         Ss06_y[12];   //[Ss06mul]
   Float_t         Ss06_z[12];   //[Ss06mul]
   Float_t         Ss06_e[12];   //[Ss06mul]
   Float_t         Ss06_theta[12];   //[Ss06mul]
   Float_t         Ss06_phi[12];   //[Ss06mul]
   UInt_t          Ss07mul;
   Float_t         Ss07_x[1];   //[Ss07mul]
   Float_t         Ss07_y[1];   //[Ss07mul]
   Float_t         Ss07_z[1];   //[Ss07mul]
   Float_t         Ss07_e[1];   //[Ss07mul]
   Float_t         Ss07_theta[1];   //[Ss07mul]
   Float_t         Ss07_phi[1];   //[Ss07mul]
   UInt_t          Ss08mul;
   Float_t         Ss08_x[1];   //[Ss08mul]
   Float_t         Ss08_y[1];   //[Ss08mul]
   Float_t         Ss08_z[1];   //[Ss08mul]
   Float_t         Ss08_e[1];   //[Ss08mul]
   Float_t         Ss08_theta[1];   //[Ss08mul]
   Float_t         Ss08_phi[1];   //[Ss08mul]
   UInt_t          boxmul;
   Float_t         box_e[23];   //[boxmul]
   Float_t         box_theta[23];   //[boxmul]
   Float_t         box_phi[23];   //[boxmul]

   // List of branches
   TBranch        *b_X0;   //!
   TBranch        *b_Y0;   //!
   TBranch        *b_T0;   //!
   TBranch        *b_fra_beta;   //!
   TBranch        *b_P_p;   //!
   TBranch        *b_P_e;   //!
   TBranch        *b_P_theta;   //!
   TBranch        *b_P_phi;   //!
   TBranch        *b_Xbn;   //!
   TBranch        *b_Xbi;   //!
   TBranch        *b_Xbe;   //!
   TBranch        *b_Xbhe;   //!
   TBranch        *b_Xbt;   //!
   TBranch        *b_Ss01mul;   //!
   TBranch        *b_Ss01_x;   //!
   TBranch        *b_Ss01_y;   //!
   TBranch        *b_Ss01_z;   //!
   TBranch        *b_Ss01_e;   //!
   TBranch        *b_Ss01_theta;   //!
   TBranch        *b_Ss01_phi;   //!
   TBranch        *b_Ss02mul;   //!
   TBranch        *b_Ss02_x;   //!
   TBranch        *b_Ss02_y;   //!
   TBranch        *b_Ss02_z;   //!
   TBranch        *b_Ss02_e;   //!
   TBranch        *b_Ss02_theta;   //!
   TBranch        *b_Ss02_phi;   //!
   TBranch        *b_Ss03mul;   //!
   TBranch        *b_Ss03_x;   //!
   TBranch        *b_Ss03_y;   //!
   TBranch        *b_Ss03_z;   //!
   TBranch        *b_Ss03_e;   //!
   TBranch        *b_Ss03_theta;   //!
   TBranch        *b_Ss03_phi;   //!
   TBranch        *b_Ss04mul;   //!
   TBranch        *b_Ss04_x;   //!
   TBranch        *b_Ss04_y;   //!
   TBranch        *b_Ss04_z;   //!
   TBranch        *b_Ss04_e;   //!
   TBranch        *b_Ss04_theta;   //!
   TBranch        *b_Ss04_phi;   //!
   TBranch        *b_Ss05mul;   //!
   TBranch        *b_Ss05_x;   //!
   TBranch        *b_Ss05_y;   //!
   TBranch        *b_Ss05_z;   //!
   TBranch        *b_Ss05_e;   //!
   TBranch        *b_Ss05_theta;   //!
   TBranch        *b_Ss05_phi;   //!
   TBranch        *b_Ss06mul;   //!
   TBranch        *b_Ss06_x;   //!
   TBranch        *b_Ss06_y;   //!
   TBranch        *b_Ss06_z;   //!
   TBranch        *b_Ss06_e;   //!
   TBranch        *b_Ss06_theta;   //!
   TBranch        *b_Ss06_phi;   //!
   TBranch        *b_Ss07mul;   //!
   TBranch        *b_Ss07_x;   //!
   TBranch        *b_Ss07_y;   //!
   TBranch        *b_Ss07_z;   //!
   TBranch        *b_Ss07_e;   //!
   TBranch        *b_Ss07_theta;   //!
   TBranch        *b_Ss07_phi;   //!
   TBranch        *b_Ss08mul;   //!
   TBranch        *b_Ss08_x;   //!
   TBranch        *b_Ss08_y;   //!
   TBranch        *b_Ss08_z;   //!
   TBranch        *b_Ss08_e;   //!
   TBranch        *b_Ss08_theta;   //!
   TBranch        *b_Ss08_phi;   //!
   TBranch        *b_boxmul;   //!
   TBranch        *b_box_e;   //!
   TBranch        *b_box_theta;   //!
   TBranch        *b_box_phi;   //!


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

