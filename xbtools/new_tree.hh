#ifndef _NEW_TREE_HH_
#define _NEW_TREE_HH_

// CB variables(protons)
UInt_t          Xb_mul;    // multiplicity of all crystals in one event
Float_t         AoZ;    //AoZ (Brho)

UInt_t          Xbp_mul;    // multiplicity of proton clusters
Float_t         Xbpe_n[64];  //summed proton energies in clusters
UInt_t					Xbp_cmul[64]; // number of crystals in each cluster
UInt_t          Xbp_cr[64];  //crystal number
Float_t         Xbp_theta[64];
Float_t         Xbp_phi[64];
Float_t         Xbp_t[64];   //times of protons

//CB variables(gammas)
UInt_t          Xbg_mul;     // multiplicity of gamma clusters
Float_t         Xbge_n[165]; //summed gamma energies in clusters
UInt_t		Xbg_cmul[165]; // number of crystals in each cluster
UInt_t          Xbg_cr[165]; //crystal number
Float_t         Xbg_theta[165];
Float_t         Xbg_phi[165];
Float_t         Xbg_t[165];  //times of gammas

// SST HIT variables
Double_t       Ss01_theta[500];
Double_t       Ss01_phi[500];
Float_t        Ss01z[500];

Double_t       Ss02_theta[500];
Double_t       Ss02_phi[500];
Float_t        Ss02z[500];

Double_t       Ss03_theta[500];
Double_t       Ss03_phi[500];
Float_t        Ss03z[500];

Double_t       Ss04_theta[500];
Double_t       Ss04_phi[500];
Float_t        Ss04z[500];

Double_t       Ss05_theta[500];
Double_t       Ss05_phi[500];
Float_t        Ss05z[500];

Double_t       Ss06_theta[500];
Double_t       Ss06_phi[500];
Float_t        Ss06z[500];

new_tree->Branch("Xb_mul"   	,   &Xb_mul     , "Xb_mul/i" );
new_tree->Branch("AoZ"      	,   &AoZ	      , "AoZ/F" );

new_tree->Branch("Xbp_mul"   	,   &Xbp_mul     , "Xbp_mul/i" );
new_tree->Branch("Xbp_cmul" 	,   Xbp_cmul   , "Xbp_cmul[Xbp_mul]/i" );
new_tree->Branch("Xbpe_n"    	,   Xbpe_n       , "Xbpe_n[Xbp_mul]/F");
new_tree->Branch("Xbp_cr"    	,   Xbp_cr       , "Xbp_cr[Xbp_mul]/i");
new_tree->Branch("Xbp_theta" 	,   Xbp_theta    , "Xbp_theta[Xbp_mul]/F");
new_tree->Branch("Xbp_phi"   	,   Xbp_phi      , "Xbp_phi[Xbp_mul]/F");
new_tree->Branch("Xbp_t"     	,   Xbp_t        , "Xbp_t[Xbp_mul]/F");

new_tree->Branch("Xbg_mul"   	,   &Xbg_mul     , "Xbg_mul/i" );
new_tree->Branch("Xbg_cmul"   	,   Xbg_cmul     , "Xbg_cmul[Xbg_mul]/i" );
new_tree->Branch("Xbge_n"    	,   Xbge_n       , "Xbge_n[Xbg_mul]/F");
new_tree->Branch("Xbg_cr"    	,   Xbg_cr       , "Xbg_cr[Xbg_mul]/i");
new_tree->Branch("Xbg_theta" 	,   Xbg_theta    , "Xbg_theta[Xbg_mul]/F");
new_tree->Branch("Xbg_phi"   	,   Xbg_phi      , "Xbg_phi[Xbg_mul]/F");
new_tree->Branch("Xbg_t"   	,   Xbg_t        , "Xbg_t[Xbg_mul]/F");

//HIT variables
new_tree->Branch("Ss01_theta" 	, Ss01_theta  , "Ss01_theta[Ss01mul]/D");
new_tree->Branch("Ss01_phi"   	, Ss01_phi    , "Ss01_phi[Ss01mul]/D");
new_tree->Branch("Ss01z"     	, Ss01z      , "Ss01z[Ss01mul]/F");

new_tree->Branch("Ss02_theta" , Ss02_theta  , "Ss02_theta[Ss02mul]/D");
new_tree->Branch("Ss02_phi"   , Ss02_phi    , "Ss02_phi[Ss02mul]/D");
new_tree->Branch("Ss02z"     , Ss02z      , "Ss02z[Ss02mul]/F");

new_tree->Branch("Ss03_theta" , Ss03_theta  , "Ss03_theta[Ss03mul]/D");
new_tree->Branch("Ss03_phi"   , Ss03_phi    , "Ss03_phi[Ss03mul]/D");
new_tree->Branch("Ss03z"     , Ss03z      , "Ss03z[Ss03mul]/F");

new_tree->Branch("Ss04_theta" , Ss04_theta  , "Ss04_theta[Ss04mul]/D");
new_tree->Branch("Ss04_phi"   , Ss04_phi    , "Ss04_phi[Ss04mul]/D");
new_tree->Branch("Ss04z"     , Ss04z      , "Ss04z[Ss04mul]/F");

new_tree->Branch("Ss05_theta" , Ss05_theta  , "Ss05_theta[Ss05mul]/D");
new_tree->Branch("Ss05_phi"   , Ss05_phi    , "Ss05_phi[Ss05mul]/D");
new_tree->Branch("Ss05z"     , Ss05z      , "Ss05z[Ss05mul]/F");

new_tree->Branch("Ss06_theta" , Ss06_theta  , "Ss06_theta[Ss06mul]/D");
new_tree->Branch("Ss06_phi"   , Ss06_phi    , "Ss06_phi[Ss06mul]/D");
new_tree->Branch("Ss06z"     , Ss06z      , "Ss06z[Ss06mul]/F");

#endif
