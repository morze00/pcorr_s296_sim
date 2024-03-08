#include "track_R3Bsim.h"
#include "libs.hh"
#include "xbtools.h"
#include "constants.hh"
#include "xb_list_energy.hh"
void track::Analyse(Char_t *out_file)
{
	if (fChain == 0) return;
	loadData();
	TFile* _file = new TFile(out_file,"Recreate");
	TTree *new_tree = fChain->CloneTree(0);	
#include"new_tree_R3Bsim.hh"
	//Containers of the crystals for the addback routine   
	std::list<xble> p_list; //output container of proton clusters 
	std::list<xble>::iterator p_iter;
	//Container of the crystals with doppler corrected gamma energies
	std::list<xble> g_doppler_list;
	std::list<xble>::iterator g_doppler_iter;
	xble cb_cryst; //Temporary dummy object needed for the XB loop	
	//11B beta at the half-thickness of the CH2 target(2.3mm) (from Alexandra's code)
	//double beta_atima = 0.711646; //half-thickness of the CH2 target(2.3mm)
	double beta_atima = 0.711398;//half-thickness of the CT(2.01mm)
	//double beta_atima = 0.7097;//half-thickness of the CT(2.01mm)

	/********************* EVENT LOOP ***************************/
	Long64_t nentries = fChain->GetEntries();
	Long64_t nbytes = 0, nb = 0;
	for(Long64_t jentry=0; jentry<nentries; jentry++){

		//std::cout << "\n\n ############  NEXT EVENT No: " << Evnt << "  ############"<<std::endl;
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if(jentry%1000==0) std::cout << jentry << " of " << nentries << " (" << (float)jentry/nentries*100 << "%)" << std::endl;
		//============ Sorting conditions

		//No overflows!
		bool bad_event = false;
		for (UInt_t cr = 0; cr<Xbn; cr++)
		{
					if(Xbe[cr]>GAMMA_MAX_THRESHOLD && isnan(Xbhe[cr]))
					bad_event = true;
		}
		if(bad_event) continue;

		memset(&Xb_mul,0,sizeof(Xb_mul));
		memset(&Xbh_mul,0,sizeof(Xbh_mul));
		memset(Xbh_cmul,0,sizeof(Xbh_cmul));
		memset(Xbh_e,0,sizeof(Xbh_e));
		memset(Xbh_cr,0,sizeof(Xbh_cr));
		memset(Xbh_theta,0,sizeof(Xbh_theta));
		memset(Xbh_phi,0,sizeof(Xbh_phi));
		memset(Xbh_t,0,sizeof(Xbh_t));

		memset(&Xbg_mul,0,sizeof(Xbg_mul));
		memset(Xbg_cmul,0,sizeof(Xbg_cmul));
		memset(Xbg_elab,0,sizeof(Xbg_elab));
		memset(Xbg_cr,0,sizeof(Xbg_cr));
		memset(Xbg_theta,0,sizeof(Xbg_theta));
		memset(Xbg_phi,0,sizeof(Xbg_phi));
		memset(Xbg_t,0,sizeof(Xbg_t));

		//***************** XB MULTIPLICITY LOOP ************************
		//empty containers from previously stored values
		std::list<xble> xb_global_list; //global container of all crystals
		std::list<xble>::iterator xb_iter;
		xb_global_list.clear(); 

		//~~~ Storing CB data into the global container ~~~
		UInt_t xbmul_count=0; 	int cl_mul=1;
		for (UInt_t cr = 0; cr<Xbn; cr++){
			cb_cryst.reset();

			//if(Xbe[cr]>GAMMA_MIN_THRESHOLD || //good gamma 
			if((Xbe[cr]>GAMMA_MIN_THRESHOLD && !isnan(Xbt[cr])) || //good gamma 
					(Xbe[cr]>GAMMA_MAX_THRESHOLD && isnan(Xbhe[cr])) || //overflow in backward hemisphere
					(Xbhe[cr]>PROTON_THRESHOLD && (Xbe[cr]>GAMMA_MIN_THRESHOLD || isinf(Xbe[cr])))) //good proton
			{
				xbmul_count++;
				cb_cryst.set(Xbi[cr],Xbhe[cr],Xbe[cr],Xbt[cr],cl_mul);
				xb_global_list.push_back(cb_cryst);
			}
		}

		//if(xbmul_count==0) continue;
		//only crystals with energy
		Xb_mul = xbmul_count;
		// sort the list by energy
		xb_global_list.sort();
		//std::cout << "\nSorted global list: \n";
		//print_all(xb_global_list);
		p_list.clear();
		g_doppler_list.clear();
		//container of gamma clusters
		std::list<xble> g_list;
		std::list<xble>::iterator g_iter;
		g_list.clear();

		//adbp(&xb_global_list, &p_list, &g_list);// add back routine
		//adbp_1(&xb_global_list, &p_list, &g_list);// add back routine
		//adbp_1t(&xb_global_list, &p_list, &g_list, T0, h_tdiff);// add back routine
		//adbp_2(&xb_global_list, &p_list, &g_list);// add back routine
		//adbp_3(&xb_global_list, &p_list, &g_list);// add back routine
		//adbp_4(&xb_global_list, &p_list, &g_list);// add back routine
		adbp_all(&xb_global_list, &p_list, &g_list);// add back routine

		//*************************************************************
		//number of "good" proton clusters found
		Xbh_mul = p_list.size();
		//if(Xbh_mul<1) continue;
		//if(Xbp_mul!=2) continue;
		//std::cout << "\n\n************** Proton list: \n";
		//print_all(p_list);
		//std::cout << "\n\n************** Gamma list: \n";
		//print_all(g_list);

		UInt_t it=0;

		for(p_iter = p_list.begin(); p_iter != p_list.end(); p_iter++)
		{
			Xbh_cr[it] = (*p_iter)._number;
			Xbh_e[it] = (*p_iter)._energy_p;
			Xbh_t[it] = (*p_iter)._time;
			Xbh_theta[it] = mycr[((*p_iter)._number)-1].theta;
			Xbh_phi[it] = mycr[((*p_iter)._number)-1].phi;
			Xbh_cmul[it] = (*p_iter)._mul;
			it++;
		}

		doppler(&g_list,beta_atima,&g_doppler_list);
		//doppler(&g_list,fra_beta,&g_doppler_list);
		it=0;//reset for the next loop
		Xbg_mul = g_doppler_list.size();

		//Sort gammas by energy
		g_doppler_list.sort();
		//Write into the output tree
		for(g_doppler_iter = g_doppler_list.begin(); g_doppler_iter != g_doppler_list.end(); g_doppler_iter++)
		{
			Xbg_cr[it] = (*g_doppler_iter)._number;
			Xbg_elab[it] = (*g_doppler_iter)._energy_g;
			Xbg_t[it] = (*g_doppler_iter)._time;
			Xbg_cmul[it] = (*g_doppler_iter)._mul;
			Xbg_theta[it] = mycr[((*g_doppler_iter)._number)-1].theta;
			Xbg_phi[it] = mycr[((*g_doppler_iter)._number)-1].phi;
			it++;
		}

		//std::cout << "\n\n************** Gamma doppler corrected list: \n";
		//print_all(g_doppler_list);

		new_tree->Fill();
	}//end_of_Eventloop

	//_file->cd();
	//new_tree->Write();
	//_file->Close();
	new_tree->Print();
	new_tree->AutoSave();
	//h_tdiff->Write();
	delete _file;

	//TFile f = TFile("hist_tdiff","Recreate");
	//TCanvas * c1 = new TCanvas("c1","c1",800,600);
	//c1->cd();

	//f->Close();


	return;
}//end_of_Analyse

//track::track(TTree *tree) : fChain(0) 
track::track(TTree *tree) 
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	//if (tree == 0) {
//		TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ppn_gs_TGeant4_addback.root");
//		if (!f || !f->IsOpen()) {
//			f = new TFile("ppn_gs_TGeant4_addback.root");
//		}
//		f->GetObject("h509",tree);

	if(tree!=0)
//}
	Init(tree);
}

void track::AddChain(TTree* tree)
{
	Init(tree);
}

track::~track()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}



Int_t track::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t track::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void track::Init(TTree *tree)
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
		
	 fChain->SetBranchAddress("X0", &X0, &b_X0);
   fChain->SetBranchAddress("Y0", &Y0, &b_Y0);
   fChain->SetBranchAddress("T0", &T0, &b_T0);
   fChain->SetBranchAddress("fra_beta", &fra_beta, &b_fra_beta);
   fChain->SetBranchAddress("P_p", P_p, &b_P_p);
   fChain->SetBranchAddress("P_e", P_e, &b_P_e);
   fChain->SetBranchAddress("P_theta", P_theta, &b_P_theta);
   fChain->SetBranchAddress("P_phi", P_phi, &b_P_phi);
   fChain->SetBranchAddress("Xbn", &Xbn, &b_Xbn);
   fChain->SetBranchAddress("Xbi", Xbi, &b_Xbi);
   fChain->SetBranchAddress("Xbe", Xbe, &b_Xbe);
   fChain->SetBranchAddress("Xbhe", Xbhe, &b_Xbhe);
   fChain->SetBranchAddress("Xbt", Xbt, &b_Xbt);
   fChain->SetBranchAddress("Ss01mul", &Ss01mul, &b_Ss01mul);
   fChain->SetBranchAddress("Ss01_x", Ss01_x, &b_Ss01_x);
   fChain->SetBranchAddress("Ss01_y", Ss01_y, &b_Ss01_y);
   fChain->SetBranchAddress("Ss01_z", Ss01_z, &b_Ss01_z);
   fChain->SetBranchAddress("Ss01_e", Ss01_e, &b_Ss01_e);
   fChain->SetBranchAddress("Ss01_theta", Ss01_theta, &b_Ss01_theta);
   fChain->SetBranchAddress("Ss01_phi", Ss01_phi, &b_Ss01_phi);
   fChain->SetBranchAddress("Ss02mul", &Ss02mul, &b_Ss02mul);
   fChain->SetBranchAddress("Ss02_x", Ss02_x, &b_Ss02_x);
   fChain->SetBranchAddress("Ss02_y", Ss02_y, &b_Ss02_y);
   fChain->SetBranchAddress("Ss02_z", Ss02_z, &b_Ss02_z);
   fChain->SetBranchAddress("Ss02_e", Ss02_e, &b_Ss02_e);
   fChain->SetBranchAddress("Ss02_theta", Ss02_theta, &b_Ss02_theta);
   fChain->SetBranchAddress("Ss02_phi", Ss02_phi, &b_Ss02_phi);
   fChain->SetBranchAddress("Ss03mul", &Ss03mul, &b_Ss03mul);
   fChain->SetBranchAddress("Ss03_x", Ss03_x, &b_Ss03_x);
   fChain->SetBranchAddress("Ss03_y", Ss03_y, &b_Ss03_y);
   fChain->SetBranchAddress("Ss03_z", Ss03_z, &b_Ss03_z);
   fChain->SetBranchAddress("Ss03_e", Ss03_e, &b_Ss03_e);
   fChain->SetBranchAddress("Ss03_theta", Ss03_theta, &b_Ss03_theta);
   fChain->SetBranchAddress("Ss03_phi", Ss03_phi, &b_Ss03_phi);
   fChain->SetBranchAddress("Ss04mul", &Ss04mul, &b_Ss04mul);
   fChain->SetBranchAddress("Ss04_x", Ss04_x, &b_Ss04_x);
   fChain->SetBranchAddress("Ss04_y", Ss04_y, &b_Ss04_y);
   fChain->SetBranchAddress("Ss04_z", Ss04_z, &b_Ss04_z);
   fChain->SetBranchAddress("Ss04_e", Ss04_e, &b_Ss04_e);
   fChain->SetBranchAddress("Ss04_theta", Ss04_theta, &b_Ss04_theta);
   fChain->SetBranchAddress("Ss04_phi", Ss04_phi, &b_Ss04_phi);
   fChain->SetBranchAddress("Ss05mul", &Ss05mul, &b_Ss05mul);
   fChain->SetBranchAddress("Ss05_x", Ss05_x, &b_Ss05_x);
   fChain->SetBranchAddress("Ss05_y", Ss05_y, &b_Ss05_y);
   fChain->SetBranchAddress("Ss05_z", Ss05_z, &b_Ss05_z);
   fChain->SetBranchAddress("Ss05_e", Ss05_e, &b_Ss05_e);
   fChain->SetBranchAddress("Ss05_theta", Ss05_theta, &b_Ss05_theta);
   fChain->SetBranchAddress("Ss05_phi", Ss05_phi, &b_Ss05_phi);
   fChain->SetBranchAddress("Ss06mul", &Ss06mul, &b_Ss06mul);
   fChain->SetBranchAddress("Ss06_x", Ss06_x, &b_Ss06_x);
   fChain->SetBranchAddress("Ss06_y", Ss06_y, &b_Ss06_y);
   fChain->SetBranchAddress("Ss06_z", Ss06_z, &b_Ss06_z);
   fChain->SetBranchAddress("Ss06_e", Ss06_e, &b_Ss06_e);
   fChain->SetBranchAddress("Ss06_theta", Ss06_theta, &b_Ss06_theta);
   fChain->SetBranchAddress("Ss06_phi", Ss06_phi, &b_Ss06_phi);
   fChain->SetBranchAddress("Ss07mul", &Ss07mul, &b_Ss07mul);
   fChain->SetBranchAddress("Ss07_x", &Ss07_x, &b_Ss07_x);
   fChain->SetBranchAddress("Ss07_y", &Ss07_y, &b_Ss07_y);
   fChain->SetBranchAddress("Ss07_z", &Ss07_z, &b_Ss07_z);
   fChain->SetBranchAddress("Ss07_e", &Ss07_e, &b_Ss07_e);
   fChain->SetBranchAddress("Ss07_theta", &Ss07_theta, &b_Ss07_theta);
   fChain->SetBranchAddress("Ss07_phi", &Ss07_phi, &b_Ss07_phi);
   fChain->SetBranchAddress("Ss08mul", &Ss08mul, &b_Ss08mul);
   fChain->SetBranchAddress("Ss08_x", &Ss08_x, &b_Ss08_x);
   fChain->SetBranchAddress("Ss08_y", &Ss08_y, &b_Ss08_y);
   fChain->SetBranchAddress("Ss08_z", &Ss08_z, &b_Ss08_z);
   fChain->SetBranchAddress("Ss08_e", &Ss08_e, &b_Ss08_e);
   fChain->SetBranchAddress("Ss08_theta", &Ss08_theta, &b_Ss08_theta);
   fChain->SetBranchAddress("Ss08_phi", &Ss08_phi, &b_Ss08_phi);
   fChain->SetBranchAddress("boxmul", &boxmul, &b_boxmul);
   fChain->SetBranchAddress("box_e", box_e, &b_box_e);
   fChain->SetBranchAddress("box_theta", box_theta, &b_box_theta);
   fChain->SetBranchAddress("box_phi", box_phi, &b_box_phi);



	Notify();
}

Bool_t track::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void track::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}

Int_t track::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}

void print_all(const std::list<xble> &rl)
{
	//std::cout << "The class room contains:" << std::endl;
	std::cout << "The list of crystal/cluster hits is:" << std::endl;
	std::list<xble>::const_iterator cii;
	for(cii=rl.begin() ; cii!=rl.end() ; cii++ )
		cii->print();
}

