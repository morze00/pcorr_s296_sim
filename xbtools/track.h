//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb  4 15:49:08 2010 by ROOT version 5.22/00
// from TTree h509/track info
// found on file: /d/land/land/vpanin/Analysis/p_corr/tracked_allZ_CH2_1350skt.root
//////////////////////////////////////////////////////////

#ifndef track_h
#define track_h

#include "libs.hh"
//#include "LibPerso.h"

class track {
	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Declaration of leaf types
		Int_t           Evnt;
		Int_t           Trig;
		Int_t           Tpat;
		Float_t         Tprev;
		Float_t         Tprev2;
		Float_t         Tnext;
		Float_t         Inz;
		Float_t         Inaoverz;
		Float_t         Inbeta;
		Float_t         Indx;
		Float_t         Indy;
		Float_t         Indz;
		Float_t         X0;
		Float_t         Y0;
		Float_t         T0;
		Float_t         Gf1_01x;
		Float_t         Gf2_01x;
		Float_t         Po01t;
		Float_t         Po01e;
		Float_t         Ps01x;
		Float_t         Ps01y;
		Float_t         Ps01e;
		Float_t         Ps02x;
		Float_t         Ps02y;
		Float_t         Ps02e;
		Int_t           Pd1mul;
		Float_t         Pd1x[30];   //[Pd1mul]
		Float_t         Pd1y[30];   //[Pd1mul]
		Float_t         Pd1hq[30];   //[Pd1mul]
		Int_t           Pd2mul;
		Float_t         Pd2x[30];   //[Pd1mul]
		Float_t         Pd2y[30];   //[Pd1mul]
		Float_t         Pd2hq[30];   //[Pd1mul]
		Int_t           Ntfmul;
		Float_t         Ntfe[4];   //[Ntfmul]
		Float_t         Ntft[4];   //[Ntfmul]
		Float_t         Ntfx[4];   //[Ntfmul]
		Float_t         Ntfy[4];   //[Ntfmul]
		Float_t         Ntfz[4];   //[Ntfmul]
		Int_t           Tfmul;
		Float_t         Tfe[4];   //[Tfmul]
		Float_t         Tft[4];   //[Tfmul]
		Float_t         Tfx[4];   //[Tfmul]
		Float_t         Tfy[4];   //[Tfmul]
		Int_t           Nhmul;
		Float_t         Nhe[57];   //[Nhmul]
		Float_t         Nht[57];   //[Nhmul]
		Float_t         Nhx[57];   //[Nhmul]
		Float_t         Nhy[57];   //[Nhmul]
		Float_t         Nhz[57];   //[Nhmul]
		Int_t           Ntmul;
		Int_t           Nttag;
		Float_t         Nttype[57];   //[Ntmul]
		Float_t         Ntth[57];   //[Ntmul]
		Float_t         Ntphi[57];   //[Ntmul]
		Float_t         Ntrad[57];   //[Ntmul]
		Float_t         Ntvel[57];   //[Ntmul]
		Float_t         Ntfwz[57];   //[Ntmul]
		Float_t         Ntbkz[57];   //[Ntmul]
		Float_t         Ntqual[57];   //[Ntmul]
		Float_t         Ntnhit[57];   //[Ntmul]
		Float_t         Xbsume;
		Int_t           Xbmul;
		Int_t           Xbi[156];   //[Xbmul]
		Float_t         Xbt[156];   //[Xbmul]
		Float_t         Xbe[156];   //[Xbmul]
		Float_t         Xbpe[156];   //[Xbmul]
		Int_t           Ss01mul;
		Float_t         Ss01e[400];   //[Ss01mul]
		Float_t         Ss01x[400];   //[Ss01mul]
		Float_t         Ss01y[400];   //[Ss01mul]
		Int_t           Ss02mul;
		Float_t         Ss02e[400];   //[Ss02mul]
		Float_t         Ss02x[400];   //[Ss02mul]
		Float_t         Ss02y[400];   //[Ss02mul]
		Int_t           Ss03mul;
		Float_t         Ss03e[400];   //[Ss03mul]
		Float_t         Ss03x[400];   //[Ss03mul]
		Float_t         Ss03y[400];   //[Ss03mul]
		Int_t           Ss04mul;
		Float_t         Ss04e[400];   //[Ss04mul]
		Float_t         Ss04x[400];   //[Ss04mul]
		Float_t         Ss04y[400];   //[Ss04mul]
		Int_t           Ss05mul;
		Float_t         Ss05e[400];   //[Ss05mul]
		Float_t         Ss05x[400];   //[Ss05mul]
		Float_t         Ss05y[400];   //[Ss05mul]
		Int_t           Ss06mul;
		Float_t         Ss06e[400];   //[Ss06mul]
		Float_t         Ss06x[400];   //[Ss06mul]
		Float_t         Ss06y[400];   //[Ss06mul]
		Int_t           Ss01_smul;
		Float_t         Ss01_se[50];   //[Ss01_smul]
		Float_t         Ss01_spos[50];   //[Ss01_smul]
		Float_t         Ss01_sbw[50];   //[Ss01_smul]
		Float_t         Ss01_sarea[50];   //[Ss01_smul]
		Float_t         Ss01_seta[50];   //[Ss01_smul]
		Int_t           Ss01_kmul;
		Float_t         Ss01_ke[50];   //[Ss01_kmul]
		Float_t         Ss01_kpos[50];   //[Ss01_kmul]
		Float_t         Ss01_kbw[50];   //[Ss01_kmul]
		Float_t         Ss01_karea[50];   //[Ss01_kmul]
		Float_t         Ss01_keta[50];   //[Ss01_kmul]
		Int_t           Ss02_smul;
		Float_t         Ss02_se[50];   //[Ss02_smul]
		Float_t         Ss02_spos[50];   //[Ss02_smul]
		Float_t         Ss02_sbw[50];   //[Ss02_smul]
		Float_t         Ss02_sarea[50];   //[Ss02_smul]
		Float_t         Ss02_seta[50];   //[Ss02_smul]
		Int_t           Ss02_kmul;
		Float_t         Ss02_ke[50];   //[Ss02_kmul]
		Float_t         Ss02_kpos[50];   //[Ss02_kmul]
		Float_t         Ss02_kbw[50];   //[Ss02_kmul]
		Float_t         Ss02_karea[50];   //[Ss02_kmul]
		Float_t         Ss02_keta[50];   //[Ss02_kmul]
		Int_t           Ss03_smul;
		Float_t         Ss03_se[50];   //[Ss03_smul]
		Float_t         Ss03_spos[50];   //[Ss03_smul]
		Float_t         Ss03_sbw[50];   //[Ss03_smul]
		Float_t         Ss03_sarea[50];   //[Ss03_smul]
		Float_t         Ss03_seta[50];   //[Ss03_smul]
		Int_t           Ss03_kmul;
		Float_t         Ss03_ke[50];   //[Ss03_kmul]
		Float_t         Ss03_kpos[50];   //[Ss03_kmul]
		Float_t         Ss03_kbw[50];   //[Ss03_kmul]
		Float_t         Ss03_karea[50];   //[Ss03_kmul]
		Float_t         Ss03_keta[50];   //[Ss03_kmul]
		Int_t           Ss04_smul;
		Float_t         Ss04_se[50];   //[Ss04_smul]
		Float_t         Ss04_spos[50];   //[Ss04_smul]
		Float_t         Ss04_sbw[50];   //[Ss04_smul]
		Float_t         Ss04_sarea[50];   //[Ss04_smul]
		Float_t         Ss04_seta[50];   //[Ss04_smul]
		Int_t           Ss04_kmul;
		Float_t         Ss04_ke[50];   //[Ss04_kmul]
		Float_t         Ss04_kpos[50];   //[Ss04_kmul]
		Float_t         Ss04_kbw[50];   //[Ss04_kmul]
		Float_t         Ss04_karea[50];   //[Ss04_kmul]
		Float_t         Ss04_keta[50];   //[Ss04_kmul]
		Int_t           Ss05_smul;
		Float_t         Ss05_se[50];   //[Ss05_smul]
		Float_t         Ss05_spos[50];   //[Ss05_smul]
		Float_t         Ss05_sbw[50];   //[Ss05_smul]
		Float_t         Ss05_sarea[50];   //[Ss05_smul]
		Float_t         Ss05_seta[50];   //[Ss05_smul]
		Int_t           Ss05_kmul;
		Float_t         Ss05_ke[50];   //[Ss05_kmul]
		Float_t         Ss05_kpos[50];   //[Ss05_kmul]
		Float_t         Ss05_kbw[50];   //[Ss05_kmul]
		Float_t         Ss05_karea[50];   //[Ss05_kmul]
		Float_t         Ss05_keta[50];   //[Ss05_kmul]
		Int_t           Ss06_smul;
		Float_t         Ss06_se[50];   //[Ss06_smul]
		Float_t         Ss06_spos[50];   //[Ss06_smul]
		Float_t         Ss06_sbw[50];   //[Ss06_smul]
		Float_t         Ss06_sarea[50];   //[Ss06_smul]
		Float_t         Ss06_seta[50];   //[Ss06_smul]
		Int_t           Ss06_kmul;
		Float_t         Ss06_ke[50];   //[Ss06_kmul]
		Float_t         Ss06_kpos[50];   //[Ss06_kmul]
		Float_t         Ss06_kbw[50];   //[Ss06_kmul]
		Float_t         Ss06_karea[50];   //[Ss06_kmul]
		Float_t         Ss06_keta[50];   //[Ss06_kmul]
		Double_t        fra_chi2;
		Float_t         fra_Z;
		Double_t        fra_A;
		Double_t        fra_beta;
		Double_t        fra_p;
		Double_t        fra_dx;
		Double_t        fra_dy;
		Double_t        fra_dz;
		Double_t        fra_dx_sst;
		Double_t        fra_dy_sst;
		Double_t        fra_dz_sst;
		Double_t        fra_flightpath;
		Float_t         fres_x0;
		Float_t         fres_y0;
		Float_t         fres_sst1x;
		Float_t         fres_sst1y;
		Float_t         fres_sst2x;
		Float_t         fres_sst2y;
		Float_t         fres_ntfx;
		Float_t         fres_ntfy;
		Float_t         fres_ntft;
		Float_t         fres_gfi1x;
		Float_t         fres_gfi2x;
		Float_t         fres_psp1x;
		Float_t         fres_psp1y;
		Float_t         fres_psp2x;
		Float_t         fres_psp2y;
		Int_t           p_mul;
		Float_t         p_chi2[3];   //[p_mul]
		Float_t         p_x0[3];   //[p_mul]
		Float_t         p_y0[3];   //[p_mul]
		Float_t         p_z0[3];   //[p_mul]
		Float_t         p_dx[3];   //[p_mul]
		Float_t         p_dy[3];   //[p_mul]
		Float_t         p_dz[3];   //[p_mul]
		Float_t         p_beta[3];   //[p_mul]
		Float_t         p_beta_tof[3];   //[p_mul]
		Float_t         p_flightpath[3];   //[p_mul]
		Float_t         p_measured_tof[3];   //[p_mul]
		Float_t         p_ekin[3];   //[p_mul]
		Float_t         p_p[3];   //[p_mul]
		Int_t           p_pid[3];   //[p_mul]
		Int_t           p_ss1i[2];   //[p_mul]
		Int_t           p_ss2i[2];   //[p_mul]
		Int_t           p_pd1i[2];   //[p_mul]
		Int_t           p_pd2i[2];   //[p_mul]
		Int_t           p_tfi[3];   //[p_mul]
		Float_t         p_rx0[3];   //[p_mul]
		Float_t         p_ry0[3];   //[p_mul]
		Float_t         p_rss1x[3];   //[p_mul]
		Float_t         p_rss1y[3];   //[p_mul]
		Float_t         p_rss2x[3];   //[p_mul]
		Float_t         p_rss2y[3];   //[p_mul]
		Float_t         p_rpd1x[3];   //[p_mul]
		Float_t         p_rpd1y[3];   //[p_mul]
		Float_t         p_rpd2x[3];   //[p_mul]
		Float_t         p_rpd2y[3];   //[p_mul]
		Float_t         p_rtfx[3];   //[p_mul]
		Float_t         p_rtfy[3];   //[p_mul]
		Float_t         p_rtft[3];   //[p_mul]

		// List of branches
		TBranch        *b_Evnt;   //!
		TBranch        *b_Trig;   //!
		TBranch        *b_Tpat;   //!
		TBranch        *b_Tprev;   //!
		TBranch        *b_Tprev2;   //!
		TBranch        *b_Tnext;   //!
		TBranch        *b_Inz;   //!
		TBranch        *b_Inaover;   //!
		TBranch        *b_Inbeta;   //!
		TBranch        *b_Indx;   //!
		TBranch        *b_Indy;   //!
		TBranch        *b_Indz;   //!
		TBranch        *b_X0;   //!
		TBranch        *b_Y0;   //!
		TBranch        *b_T0;   //!
		TBranch        *b_Gf1_01x;   //!
		TBranch        *b_Gf2_01x;   //!
		TBranch        *b_Po01t;   //!
		TBranch        *b_Po01e;   //!
		TBranch        *b_Ps01x;   //!
		TBranch        *b_Ps01y;   //!
		TBranch        *b_Ps01e;   //!
		TBranch        *b_Ps02x;   //!
		TBranch        *b_Ps02y;   //!
		TBranch        *b_Ps02e;   //!
		TBranch        *b_Pd1mul;   //!
		TBranch        *b_Pd1x;   //!
		TBranch        *b_Pd1y;   //!
		TBranch        *b_Pd1hq;   //!
		TBranch        *b_Pd2mul;   //!
		TBranch        *b_Pd2x;   //!
		TBranch        *b_Pd2y;   //!
		TBranch        *b_Pd2hq;   //!
		TBranch        *b_Ntfmul;   //!
		TBranch        *b_Ntfe;   //!
		TBranch        *b_Ntft;   //!
		TBranch        *b_Ntfx;   //!
		TBranch        *b_Ntfy;   //!
		TBranch        *b_Ntfz;   //!
		TBranch        *b_Tfmul;   //!
		TBranch        *b_Tfe;   //!
		TBranch        *b_Tft;   //!
		TBranch        *b_Tfx;   //!
		TBranch        *b_Tfy;   //!
		TBranch        *b_Nhmul;   //!
		TBranch        *b_Nhe;   //!
		TBranch        *b_Nht;   //!
		TBranch        *b_Nhx;   //!
		TBranch        *b_Nhy;   //!
		TBranch        *b_Nhz;   //!
		TBranch        *b_Ntmul;   //!
		TBranch        *b_Nttag;   //!
		TBranch        *b_Nttype;   //!
		TBranch        *b_Ntth;   //!
		TBranch        *b_Ntphi;   //!
		TBranch        *b_Ntrad;   //!
		TBranch        *b_Ntvel;   //!
		TBranch        *b_Ntfwz;   //!
		TBranch        *b_Ntbkz;   //!
		TBranch        *b_Ntqual;   //!
		TBranch        *b_Ntnhit;   //!
		TBranch        *b_Xbsume;   //!
		TBranch        *b_Xbmul;   //!
		TBranch        *b_Xbi;   //!
		TBranch        *b_Xbt;   //!
		TBranch        *b_Xbe;   //!
		TBranch        *b_Xbpe;   //!
		TBranch        *b_Ss01mul;   //!
		TBranch        *b_Ss01e;   //!
		TBranch        *b_Ss01x;   //!
		TBranch        *b_Ss01y;   //!
		TBranch        *b_Ss02mul;   //!
		TBranch        *b_Ss02e;   //!
		TBranch        *b_Ss02x;   //!
		TBranch        *b_Ss02y;   //!
		TBranch        *b_Ss03mul;   //!
		TBranch        *b_Ss03e;   //!
		TBranch        *b_Ss03x;   //!
		TBranch        *b_Ss03y;   //!
		TBranch        *b_Ss04mul;   //!
		TBranch        *b_Ss04e;   //!
		TBranch        *b_Ss04x;   //!
		TBranch        *b_Ss04y;   //!
		TBranch        *b_Ss05mul;   //!
		TBranch        *b_Ss05e;   //!
		TBranch        *b_Ss05x;   //!
		TBranch        *b_Ss05y;   //!
		TBranch        *b_Ss06mul;   //!
		TBranch        *b_Ss06e;   //!
		TBranch        *b_Ss06x;   //!
		TBranch        *b_Ss06y;   //!
		TBranch        *b_Ss01_smul;   //!
		TBranch        *b_Ss01_se;   //!
		TBranch        *b_Ss01_spos;   //!
		TBranch        *b_Ss01_sbw;   //!
		TBranch        *b_Ss01_sarea;   //!
		TBranch        *b_Ss01_seta;   //!
		TBranch        *b_Ss01_kmul;   //!
		TBranch        *b_Ss01_ke;   //!
		TBranch        *b_Ss01_kpos;   //!
		TBranch        *b_Ss01_kbw;   //!
		TBranch        *b_Ss01_karea;   //!
		TBranch        *b_Ss01_keta;   //!
		TBranch        *b_Ss02_smul;   //!
		TBranch        *b_Ss02_se;   //!
		TBranch        *b_Ss02_spos;   //!
		TBranch        *b_Ss02_sbw;   //!
		TBranch        *b_Ss02_sarea;   //!
		TBranch        *b_Ss02_seta;   //!
		TBranch        *b_Ss02_kmul;   //!
		TBranch        *b_Ss02_ke;   //!
		TBranch        *b_Ss02_kpos;   //!
		TBranch        *b_Ss02_kbw;   //!
		TBranch        *b_Ss02_karea;   //!
		TBranch        *b_Ss02_keta;   //!
		TBranch        *b_Ss03_smul;   //!
		TBranch        *b_Ss03_se;   //!
		TBranch        *b_Ss03_spos;   //!
		TBranch        *b_Ss03_sbw;   //!
		TBranch        *b_Ss03_sarea;   //!
		TBranch        *b_Ss03_seta;   //!
		TBranch        *b_Ss03_kmul;   //!
		TBranch        *b_Ss03_ke;   //!
		TBranch        *b_Ss03_kpos;   //!
		TBranch        *b_Ss03_kbw;   //!
		TBranch        *b_Ss03_karea;   //!
		TBranch        *b_Ss03_keta;   //!
		TBranch        *b_Ss04_smul;   //!
		TBranch        *b_Ss04_se;   //!
		TBranch        *b_Ss04_spos;   //!
		TBranch        *b_Ss04_sbw;   //!
		TBranch        *b_Ss04_sarea;   //!
		TBranch        *b_Ss04_seta;   //!
		TBranch        *b_Ss04_kmul;   //!
		TBranch        *b_Ss04_ke;   //!
		TBranch        *b_Ss04_kpos;   //!
		TBranch        *b_Ss04_kbw;   //!
		TBranch        *b_Ss04_karea;   //!
		TBranch        *b_Ss04_keta;   //!
		TBranch        *b_Ss05_smul;   //!
		TBranch        *b_Ss05_se;   //!
		TBranch        *b_Ss05_spos;   //!
		TBranch        *b_Ss05_sbw;   //!
		TBranch        *b_Ss05_sarea;   //!
		TBranch        *b_Ss05_seta;   //!
		TBranch        *b_Ss05_kmul;   //!
		TBranch        *b_Ss05_ke;   //!
		TBranch        *b_Ss05_kpos;   //!
		TBranch        *b_Ss05_kbw;   //!
		TBranch        *b_Ss05_karea;   //!
		TBranch        *b_Ss05_keta;   //!
		TBranch        *b_Ss06_smul;   //!
		TBranch        *b_Ss06_se;   //!
		TBranch        *b_Ss06_spos;   //!
		TBranch        *b_Ss06_sbw;   //!
		TBranch        *b_Ss06_sarea;   //!
		TBranch        *b_Ss06_seta;   //!
		TBranch        *b_Ss06_kmul;   //!
		TBranch        *b_Ss06_ke;   //!
		TBranch        *b_Ss06_kpos;   //!
		TBranch        *b_Ss06_kbw;   //!
		TBranch        *b_Ss06_karea;   //!
		TBranch        *b_Ss06_keta;   //!
		TBranch        *b_chi2_f;   //!
		TBranch        *b_Z;   //!
		TBranch        *b_A;   //!
		TBranch        *b_beta;   //!
		TBranch        *b_p;   //!
		TBranch        *b_fdx;   //!
		TBranch        *b_fdy;   //!
		TBranch        *b_fdz;   //!
		TBranch        *b_fdxsst;   //!
		TBranch        *b_fdysst;   //!
		TBranch        *b_fdzsst;   //!
		TBranch        *b_flightpath;   //!
		TBranch        *b_x0;   //!
		TBranch        *b_y0;   //!
		TBranch        *b_sst1x;   //!
		TBranch        *b_sst1y;   //!
		TBranch        *b_sst2x;   //!
		TBranch        *b_sst2y;   //!
		TBranch        *b_ntfx;   //!
		TBranch        *b_ntfy;   //!
		TBranch        *b_ntft;   //!
		TBranch        *b_gfi1x;   //!
		TBranch        *b_gfi2x;   //!
		TBranch        *b_psp1x;   //!
		TBranch        *b_psp1y;   //!
		TBranch        *b_psp2x;   //!
		TBranch        *b_psp2y;   //!
		TBranch        *b_p_mul;   //!
		TBranch        *b_p_chi2;   //!
		TBranch        *b_p_x0;   //!
		TBranch        *b_p_y0;   //!
		TBranch        *b_p_z0;   //!
		TBranch        *b_p_dx;   //!
		TBranch        *b_p_dy;   //!
		TBranch        *b_p_dz;   //!
		TBranch        *b_p_beta;   //!
		TBranch        *b_p_beta_tof;   //!
		TBranch        *b_p_flightpath;   //!
		TBranch        *b_p_measured_tof;   //!
		TBranch        *b_p_ekin;   //!
		TBranch        *b_p_p;   //!
		TBranch        *b_p_pid;   //!
		TBranch        *b_p_ss1i;   //!
		TBranch        *b_p_ss2i;   //!
		TBranch        *b_p_pd1i;   //!
		TBranch        *b_p_pd2i;   //!
		TBranch        *b_p_tfi;   //!
		TBranch        *b_p_rx0;   //!
		TBranch        *b_p_ry0;   //!
		TBranch        *b_p_rss1x;   //!
		TBranch        *b_p_rss1y;   //!
		TBranch        *b_p_rss2x;   //!
		TBranch        *b_p_rss2y;   //!
		TBranch        *b_p_rpd1x;   //!
		TBranch        *b_p_rpd1y;   //!
		TBranch        *b_p_rpd2x;   //!
		TBranch        *b_p_rpd2y;   //!
		TBranch        *b_p_rtfx;   //!
		TBranch        *b_p_rtfy;   //!
		TBranch        *b_p_rtft;   //!



		track() {}
		track(TTree *tree);
		virtual void     AddChain(TTree* tree);
		virtual ~track();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		//virtual void     Analyse();
		virtual void     Analyse(Char_t *out_file);
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);
};

#endif

