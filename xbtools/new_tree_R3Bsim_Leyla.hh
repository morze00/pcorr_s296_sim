#ifndef _NEW_TREE_R3BSIM_HH_
#define _NEW_TREE_R3BSIM_HH_

// CB variables(protons)
UInt_t          Xb_mul;    // multiplicity of all crystals in one event
Float_t         AoZ;    //AoZ (Brho)

UInt_t          Xbh_mul;    // multiplicity of proton clusters
Float_t         Xbh_e[64];  //summed proton energies in clusters
UInt_t					Xbh_cmul[64]; // number of crystals in each cluster
UInt_t          Xbh_cr[64];  //crystal number
Float_t         Xbh_theta[64];
Float_t         Xbh_phi[64];
Float_t         Xbh_t[64];   //times of protons

//CB variables(gammas)
UInt_t          Xbg_mul;     // multiplicity of gamma clusters
Float_t         Xbg_elab[165]; //summed gamma energies in clusters
UInt_t					Xbg_cmul[165]; // number of crystals in each cluster
UInt_t          Xbg_cr[165]; //crystal number
Float_t         Xbg_theta[165];
Float_t         Xbg_phi[165];
Float_t         Xbg_t[165];  //times of gammas

new_tree->Branch("Xb_mul"   	,   &Xb_mul     , "Xb_mul/i" );

new_tree->Branch("Xbh_mul"   	,   &Xbh_mul     , "Xbh_mul/i" );
new_tree->Branch("Xbh_cmul" 	,   Xbh_cmul     , "Xbh_cmul[Xbh_mul]/i" );
new_tree->Branch("Xbh_e"    	,   Xbh_e        , "Xbh_e[Xbh_mul]/F");
new_tree->Branch("Xbh_cr"    	,   Xbh_cr       , "Xbh_cr[Xbh_mul]/i");
new_tree->Branch("Xbh_theta" 	,   Xbh_theta    , "Xbh_theta[Xbh_mul]/F");
new_tree->Branch("Xbh_phi"   	,   Xbh_phi      , "Xbh_phi[Xbh_mul]/F");
new_tree->Branch("Xbh_t"     	,   Xbh_t        , "Xbh_t[Xbh_mul]/F");

new_tree->Branch("Xbg_mul"   	,   &Xbg_mul     , "Xbg_mul/i" );
new_tree->Branch("Xbg_cmul"  	,   Xbg_cmul     , "Xbg_cmul[Xbg_mul]/i" );
new_tree->Branch("Xbg_elab"  	,   Xbg_elab     , "Xbg_elab[Xbg_mul]/F");
new_tree->Branch("Xbg_cr"    	,   Xbg_cr       , "Xbg_cr[Xbg_mul]/i");
new_tree->Branch("Xbg_theta" 	,   Xbg_theta    , "Xbg_theta[Xbg_mul]/F");
new_tree->Branch("Xbg_phi"   	,   Xbg_phi      , "Xbg_phi[Xbg_mul]/F");
new_tree->Branch("Xbg_t"   	  ,   Xbg_t        , "Xbg_t[Xbg_mul]/F");

#endif
