#ifndef __MAIN_CXX__
#define __MAIN_CXX__

#include "libs.hh"
#include "track_R3Bsim.h"

Int_t main(Int_t argc, Char_t* argv[])
{
	Char_t *in_filename=0;
	Char_t *out_filename=0;

	Bool_t NeedHelp = kTRUE;

	if (argc > 1)
	{
		for (Int_t i = 0; i < argc; i++)
		{
			if (strncmp(argv[i],"--file=",7) == 0)
			{
				in_filename = argv[i]+7;
				NeedHelp = kFALSE;
			}

			else if (strncmp(argv[i],"--output=",9) == 0)
			{
				out_filename = argv[i]+9;
			}
		}
	}


	if (NeedHelp)
	{
		std::cout << "\nOptions:                                                             " << std::endl;
		std::cout << "  --file=/path/to/your/file/filename.root   : input file               " << std::endl;
		std::cout << "  --output=/path/to/your/file/filename.root : output file name         \n" << std::endl;
		return 0;
	}

	TFile *f_input = new TFile(in_filename,"Read");

	TTree *tr = (TTree*)f_input->Get("h509");
	
	track tra;


	if(tr!=NULL) tra.AddChain(tr);
	else {
		std::cout << "CANNOT FIND THE TREE!!!" << std::endl;
		return 0;
	}
	
	if(out_filename!=NULL) 
		tra.Analyse(out_filename);
	else {
		std::cout << "\nOUTPUT FILE IS NOT SPECIFIED!!!\nType: ./p_corr \n\n" << std::endl;
		return 0;
	}

	return 0;
}
#endif
