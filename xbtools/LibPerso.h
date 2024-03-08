#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sstream>
#include <iostream>
//#include <fstream.h>



#include        "TROOT.h"
#include        "TAttText.h"
#include        "TAxis.h"
#include        "TCanvas.h"
#include        "TChain.h"
//#include        "TCutG.h"
#include        "TCut.h"
#include        "TF1.h"
#include        "TFile.h"
#include        "TGraph.h"
#include        "TGraphAsymmErrors.h"
#include        "TGraphErrors.h"
#include        "TH1.h"
#include        "TH2.h"
#include        "THistPainter.h"
#include        "TKey.h"
#include        "TLatex.h"
#include        "TLegend.h"
#include        "TMath.h"
#include        "TMatrixD.h"
#include        "TMinuit.h"
#include        "TMultiGraph.h"
#include        "TNtuple.h"
#include        "TPave.h"
#include        "TPaveText.h"
#include        "TPoint.h"
#include        "TRandom.h"
#include        "TRint.h"
#include        "TStyle.h"
#include        "TString.h"
#include        "TTree.h"

#include        "TH1F.h"
#ifndef __CINT__

#include        "TAttMarker.h"
#include        "Buttons.h"
#include        "TFrame.h"
#include        "TLine.h"
#include        "TPad.h"
#include        "TPostScript.h"
#include        "TStopwatch.h" 
#include        "TStyle.h"
#include        "TSystem.h"
#include        "TUnixSystem.h"
#include        "RQ_OBJECT.h"

#endif
 
//void SetStyle(void) 
//{ 
//  gStyle->SetFrameFillColor(0);
//  gStyle->SetFrameFillStyle(6);
//  gStyle->SetFrameBorderMode(0);
//  gStyle->SetFrameBorderSize(0);
//  gStyle->SetCanvasBorderMode(0);
//  gStyle->SetCanvasBorderSize(0);
//  gStyle->SetCanvasColor(0);
//  gStyle->SetPadColor(0);
//  gStyle->SetPadBorderMode(0);
//  gStyle->SetPadBorderSize(0);
//  gStyle->SetOptStat(0);
//  gStyle->SetOptFit(1110);
//  gStyle->SetOptTitle(1);
//  gStyle->SetOptDate(0);
//  gStyle->SetOptLogx(0);
//  gStyle->SetOptLogy(0);
//  gStyle->SetPalette(1,0);
//  gStyle->SetLabelSize(0.06);
//  gStyle->SetLabelOffset(0.01,"xyz");
//  gStyle->SetTitleOffset(1.1,"xyz");
//  //        TH1::AddDirectory (kFALSE); // don't keep histograms in memory
//}
//
//
//void ProcEvent(Int_t event, Int_t px, Int_t py, TObject *sel)
//{
//   //  print event type and current cursor position
//   TCanvas *c = (TCanvas *) gTQSender;
//   TPad *pad = (TPad *) c->GetSelectedPad();
//   
//    if(!pad) return;
//    gROOT->SetEditHistograms(kFALSE);
////    printf("event=%d, px=%d, py=%d\n", event, px, py);
////    Float_t x = pad->AbsPixeltoX(px);
////    Float_t y = pad->AbsPixeltoY(py);
////    x = pad->PadtoX(x);
////    y = pad->PadtoY(y);
////    printf("x=%.3g, y=%.3g\n",x,y);
//   if(event==kButton1Double){ 
//     pad->Pop();
//     //printf("%s %d\n",pad->GetName(),pad->GetNumber());
//     pad->cd();
//     TCanvas *c_blow;
//     TIter next(pad->GetListOfPrimitives());
//     if((TCanvas*)gROOT->GetListOfCanvases()->FindObject("c_blow")){
//       c_blow=(TCanvas*)gROOT->GetListOfCanvases()->FindObject("c_blow");
//       c_blow->Show();
//       c_blow->Clear();
//     }
//     else{
//       c_blow = new TCanvas("c_blow","blowup",750,0,743,525);
//       c_blow->SetFillColor(10); //white
//       c_blow->ToggleEventStatus();
//       //c_blow->SetCrosshair();
//       c_blow->Draw();
//     }
//     c_blow->cd();
//
//     gROOT->SetSelectedPad((TPad*)c_blow);
//     TPad *clone =(TPad*)pad->Clone();
//     clone->SetPad(0,0,1,1);
//     clone->Draw();
//     clone->Modified();
//     clone->Update();
//  }
//}
////___________________________________________________________________

