#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include <TF1.h>
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPad.h>

#define nCBins 5
#define nptBins 55
#define ntrkbins 4
#define nkbins 4

const Double_t pt_low = 120.;
const Double_t pt_high = 600.;

bool do_fit = false;

char saythis[500];

using namespace std;

TString cent[4] = {"0","1","2","3","4"};
TString trk[4] = {"0","1","2","3"};
TString kbin[4] = {"0","1","2","3"};
TString pt[56] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55"};

TString cent_tag[] = {"0-10%","10-30%","30-50%", "50-100%"};
TString trk_tag[] = {"p_{T}^{trk} > 0.7 GeV","p_{T}^{trk} > 2 GeV","p_{T}^{trk} > 4 GeV", "p_{T}^{trk} > 5 GeV"};

/////////// pp //////////////////

Double_t feta_reco_cent0(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_eta_q_cent0->FindBin(xx);
   Double_t q = par[0]*h_eta_q_cent0->GetBinContent(bin);
   Double_t g = (1-par[0])*h_eta_g_cent0->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_ref_cent0(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_q_0_0->FindBin(xx);
   Double_t q = par[0]*h_chg_ref_q_0_0->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g_0_0->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent0(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_0_0->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_0_0->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_0_0->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent0_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_0_1->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_0_1->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_0_1->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent0_pt4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_0_2->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_0_2->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_0_2->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent0_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_0_3->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_0_3->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_0_3->GetBinContent(bin);
   return q+g; 
}

Double_t f_udg_reco_cent0_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g_0_3->FindBin(xx);
   //Double_t q = par[0]*h_chg_reco_q_0_3->GetBinContent(bin);
   Double_t up = (par[0])*h_chg_reco_up_0_3->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down_0_3->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_g_0_3->GetBinContent(bin);
   //Double_t g = (par[2])*h_chg_reco_g_0_3->GetBinContent(bin);
   return up+dn+g; 
}
/*
Double_t f_udg_refl_reco_cent0_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g_0_1->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_uubar_0_1->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_ddbar_0_1->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_g_0_1->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_refl_reco_cent0_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g_0_3->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_uubar_0_3->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_ddbar_0_3->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_g_0_3->GetBinContent(bin);
   return up+dn+g; 
}
*/
/////////// cent 0 -10% //////////////////

Double_t feta_reco_cent1(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_eta_q_cent1->FindBin(xx);
   Double_t q = par[0]*h_eta_q_cent1->GetBinContent(bin);
   Double_t g = (1-par[0])*h_eta_g_cent1->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_ref_cent1(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_q_1_0->FindBin(xx);
   Double_t q = par[0]*h_chg_ref_q_1_0->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g_1_0->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent1(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_1_0->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_1_0->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_1_0->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent1_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_1_1->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_1_1->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_1_1->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent1_pt4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_1_2->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_1_2->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_1_2->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent1_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_1_3->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_1_3->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_1_3->GetBinContent(bin);
   return q+g; 
}

Double_t f_udg_reco_cent1_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g_1_3->FindBin(xx);
   //Double_t q = par[0]*h_chg_reco_q_0_3->GetBinContent(bin);
   Double_t up = (par[0])*h_chg_reco_up_1_3->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down_1_3->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_g_1_3->GetBinContent(bin);
   return up+dn+g; 
}

/////////// cent 10 -30% //////////////////

Double_t feta_reco_cent2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_eta_q_cent2->FindBin(xx);
   Double_t q = par[0]*h_eta_q_cent2->GetBinContent(bin);
   Double_t g = (1-par[0])*h_eta_g_cent2->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_ref_cent2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_q_2_0->FindBin(xx);
   Double_t q = par[0]*h_chg_ref_q_2_0->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g_2_0->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_2_0->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_2_0->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_2_0->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent2_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_2_1->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_2_1->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_2_1->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent2_pt4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_2_2->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_2_2->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_2_2->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent2_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_2_3->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_2_3->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_2_3->GetBinContent(bin);
   return q+g; 
}

Double_t f_udg_reco_cent2_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g_2_3->FindBin(xx);
   //Double_t q = par[0]*h_chg_reco_q_0_3->GetBinContent(bin);
   Double_t up = (par[0])*h_chg_reco_up_2_3->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down_2_3->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_g_2_3->GetBinContent(bin);
   return up+dn+g; 
}

/////////// cent 30 -50% //////////////////

Double_t feta_reco_cent3(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_eta_q_cent3->FindBin(xx);
   Double_t q = par[0]*h_eta_q_cent3->GetBinContent(bin);
   Double_t g = (1-par[0])*h_eta_g_cent3->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_ref_cent3(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_q_3_0->FindBin(xx);
   Double_t q = par[0]*h_chg_ref_q_3_0->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g_3_0->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent3(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_3_0->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_3_0->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_3_0->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent3_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_3_1->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_3_1->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_3_1->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent3_pt4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_3_2->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_3_2->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_3_2->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent3_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_3_3->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_3_3->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_3_3->GetBinContent(bin);
   return q+g; 
}

Double_t f_udg_reco_cent3_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g_3_3->FindBin(xx);
   //Double_t q = par[0]*h_chg_reco_q_0_3->GetBinContent(bin);
   Double_t up = (par[0])*h_chg_reco_up_3_3->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down_3_3->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_g_3_3->GetBinContent(bin);
   return up+dn+g; 
}

/////////// cent 50 -100% //////////////////

Double_t feta_reco_cent4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_eta_q_cent4->FindBin(xx);
   Double_t q = par[0]*h_eta_q_cent4->GetBinContent(bin);
   Double_t g = (1-par[0])*h_eta_g_cent4->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_ref_cent4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_q_4_0->FindBin(xx);
   Double_t q = par[0]*h_chg_ref_q_4_0->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g_4_0->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_4_0->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_4_0->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_4_0->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent4_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_4_1->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_4_1->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_4_1->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent4_pt4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_4_2->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_4_2->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_4_2->GetBinContent(bin);
   return q+g; 
}

Double_t ftotal_reco_cent4_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_q_4_3->FindBin(xx);
   Double_t q = par[0]*h_chg_reco_q_4_3->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_reco_g_4_3->GetBinContent(bin);
   return q+g; 
}

Double_t f_udg_reco_cent4_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g_4_3->FindBin(xx);
   //Double_t q = par[0]*h_chg_reco_q_0_3->GetBinContent(bin);
   Double_t up = (par[0])*h_chg_reco_up_4_3->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down_4_3->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_g_4_3->GetBinContent(bin);
   return up+dn+g; 
}

void jetchg_fit(){

  //TFile *closure_histos_pp_MC_kscan = TFile::Open("/home/dhanush/Documents/jet_chg/Pythia6_jetchg_Mar31.root");
  TFile *closure_histos_pp_MC = TFile::Open("/home/dhanush/Documents/jet_chg/Pythia6_jetchg_bkgsub_eta0p5_1p5_Apr16.root");
  TFile *closure_histos_PbPb_MC = TFile::Open("/data/jetchg_histos/P+H_jetchg_cymbal_bkgsub_eta0p5_1p5_Apr16.root");
  TFile *closure_histos_pp_data = TFile::Open("/home/dhanush/Documents/jet_chg/ppdata_jetchg_bkgsub_eta0p5_1p5_Apr16.root");
  TFile *closure_histos_PbPb_data = TFile::Open("/home/dhanush/Documents/jet_chg/PbPbdata_jetchg_bkgsub_eta0p5_1p5_Apr16.root");

  TH2F *h_chg_refpt[nCBins][ntrkbins];
  TH2F *h_chg_refpt_q[nCBins][ntrkbins];
  TH2F *h_chg_refpt_g[nCBins][ntrkbins];
  TH2F *h_chg_refpt_u[nCBins][ntrkbins];

  TH2F *h_trk_corrpt_MC[nCBins];
  TH2F *h_trk_corrpt_data[nCBins];
  TH2F *h_bkgtrk_corrpt_MC[nCBins];
  TH2F *h_bkgtrk_corrpt_data[nCBins];

  TH2F *h_chg_corrpt[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_q[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_g[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_u[nCBins][ntrkbins];

  TH2F *h_chg_corrpt_up[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_down[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_upb[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_downb[nCBins][ntrkbins];

  TH2F *h_chg_corrpt_bkg[nCBins][ntrkbins];

  TH2F *h_chg_corrpt_data[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_data_bkg[nCBins][ntrkbins];

  TH1D *h_chg_ref[nCBins][ntrkbins];
  TH1D *h_chg_ref_q[nCBins][ntrkbins];
  TH1D *h_chg_ref_g[nCBins][ntrkbins];
  TH1D *h_chg_ref_u[nCBins][ntrkbins];

  TH1D *h_trk_MC[nCBins];
  TH1D *h_trk_data[nCBins];
  TH1D *h_bkgtrk_MC[nCBins];
  TH1D *h_bkgtrk_data[nCBins];

  TH1D *h_eta_MC[nCBins];
  TH1D *h_eta_data[nCBins];
  TH1D *h_eta_MC_q[nCBins];
  TH1D *h_eta_MC_g[nCBins];

  TH1D *h_eta_MC_up[nCBins];
  TH1D *h_eta_MC_upbar[nCBins];
  TH1D *h_eta_MC_d[nCBins];
  TH1D *h_eta_MC_dbar[nCBins];

  TH1D *h_eta_MC_data_ratio[nCBins];
  TH1D *h_reco_MC_data_ratio[nCBins];

  TH1D *h_reco_MC[nCBins];
  TH1D *h_reco_data[nCBins];
  TH1D *h_reco_MC_q[nCBins];
  TH1D *h_reco_MC_g[nCBins];
  TH1D *h_reco_MC_up[nCBins];
  TH1D *h_reco_MC_upbar[nCBins];
  TH1D *h_reco_MC_d[nCBins];
  TH1D *h_reco_MC_dbar[nCBins];

  TH1D *h_chg_reco[nCBins][ntrkbins];
  TH1D *h_chg_reco_q[nCBins][ntrkbins];
  TH1D *h_chg_reco_g[nCBins][ntrkbins];
  TH1D *h_chg_reco_u[nCBins][ntrkbins];
  TH1D *h_chg_reco_bkg[nCBins][ntrkbins];

  TH1D *h_chg_reco_up[nCBins][ntrkbins];
  TH1D *h_chg_reco_down[nCBins][ntrkbins];
  TH1D *h_chg_reco_upb[nCBins][ntrkbins];
  TH1D *h_chg_reco_downb[nCBins][ntrkbins];

  TH1D *h_chg_reco_q_refl[nCBins][ntrkbins];
  TH1D *h_chg_reco_g_refl[nCBins][ntrkbins];
  TH1D *h_chg_reco_MC_refl[nCBins][ntrkbins];
  TH1D *h_chg_reco_data_refl[nCBins][ntrkbins];

  TH1D *h_chg_reco_uubar[nCBins][ntrkbins];
  TH1D *h_chg_reco_ddbar[nCBins][ntrkbins];

  TH1D *h_chg_data[nCBins][ntrkbins];
  TH1D *h_chg_data_bkg[nCBins][ntrkbins];

  TH1D *h_chg_ref_ratio[nCBins][ntrkbins];
  TH1D *h_chg_reco_ratio[nCBins][ntrkbins];
  TH1D *h_chg_data_ratio[nCBins][ntrkbins];

  TH1D *h_chg_qg_ratio[nCBins][ntrkbins];

  TProfile *h_chg_ref_avg[nCBins][ntrkbins];
  TProfile *h_chg_ref_q_avg[nCBins][ntrkbins];
  TProfile *h_chg_ref_g_avg[nCBins][ntrkbins];

  TProfile *h_chg_reco_avg[nCBins][ntrkbins];
  TProfile *h_chg_reco_q_avg[nCBins][ntrkbins];
  TProfile *h_chg_reco_g_avg[nCBins][ntrkbins];

  TProfile *h_chg_reco_up_avg[nCBins][ntrkbins];
  TProfile *h_chg_reco_down_avg[nCBins][ntrkbins];
  TProfile *h_chg_reco_upb_avg[nCBins][ntrkbins];
  TProfile *h_chg_reco_downb_avg[nCBins][ntrkbins];

  TProfile *h_chg_reco_bkg_avg[nCBins][ntrkbins];

  TProfile *h_chg_data_avg[nCBins][ntrkbins];
  TProfile *h_chg_data_bkg_avg[nCBins][ntrkbins];

  TF1 *f_chg_ref[nCBins][ntrkbins];  
  TF1 *f_chg_reco[nCBins][ntrkbins];
  TF1 *f_chg_data[nCBins][ntrkbins];
  TF1 *f_chg_udg_data[nCBins][ntrkbins];

  TF1 *f_eta_reco[nCBins];

  for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0){
      h_trk_corrpt_MC[ibin] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_trk_corrpt_cent"+cent[ibin])); 
      h_bkgtrk_corrpt_MC[ibin] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_bkgtrk_corrpt_cent"+cent[ibin])); 
      h_trk_corrpt_data[ibin] = (TH2F*)closure_histos_pp_data->Get((TString)("h_trk_corrpt_cent"+cent[ibin])); 
      h_bkgtrk_corrpt_data[ibin] = (TH2F*)closure_histos_pp_data->Get((TString)("h_bkgtrk_corrpt_cent"+cent[ibin])); 

      h_eta_MC[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_cent"+cent[ibin]));       
      h_eta_MC_q[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_q_cent"+cent[ibin]))->Clone((TString)("h_eta_q_cent"+cent[ibin]));
      h_eta_MC_g[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_g_cent"+cent[ibin]))->Clone((TString)("h_eta_g_cent"+cent[ibin]));
      h_eta_data[ibin] = (TH1D*)closure_histos_pp_data->Get((TString)("h_eta_full_cent"+cent[ibin]));
      h_eta_MC_up[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_up_cent"+cent[ibin]));
      h_eta_MC_upbar[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_upbar_cent"+cent[ibin]));
      h_eta_MC_d[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_d_cent"+cent[ibin]));
      h_eta_MC_dbar[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_dbar_cent"+cent[ibin]));

      h_reco_MC[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_cent"+cent[ibin]));
      h_reco_data[ibin] = (TH1D*)closure_histos_pp_data->Get((TString)("h_reco_corr_cent"+cent[ibin]));
      h_reco_MC_q[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_q_cent"+cent[ibin]));
      h_reco_MC_g[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_g_cent"+cent[ibin]));
      h_reco_MC_up[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_up_cent"+cent[ibin]));
      h_reco_MC_upbar[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_upbar_cent"+cent[ibin]));
      h_reco_MC_d[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_d_cent"+cent[ibin]));
      h_reco_MC_dbar[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_dbar_cent"+cent[ibin]));
    }
    else{
      h_trk_corrpt_MC[ibin] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_trk_corrpt_cent"+cent[ibin-1])); 
      h_bkgtrk_corrpt_MC[ibin] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_bkgtrk_corrpt_cent"+cent[ibin-1])); 
      h_trk_corrpt_data[ibin] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_trk_corrpt_cent"+cent[ibin-1])); 
      h_bkgtrk_corrpt_data[ibin] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_bkgtrk_corrpt_cent"+cent[ibin-1])); 

      h_eta_MC[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_cent"+cent[ibin-1]));       
      h_eta_MC_q[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_q_cent"+cent[ibin-1]))->Clone((TString)("h_eta_q_cent"+cent[ibin]));
      h_eta_MC_g[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_g_cent"+cent[ibin-1]))->Clone((TString)("h_eta_g_cent"+cent[ibin]));
      h_eta_data[ibin] = (TH1D*)closure_histos_PbPb_data->Get((TString)("h_eta_full_cent"+cent[ibin-1]));
      h_eta_MC_up[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_up_cent"+cent[ibin-1]));
      h_eta_MC_upbar[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_upbar_cent"+cent[ibin-1]));
      h_eta_MC_d[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_d_cent"+cent[ibin-1]));
      h_eta_MC_dbar[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_dbar_cent"+cent[ibin-1]));

      h_reco_MC[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_cent"+cent[ibin-1]));       
      h_reco_data[ibin] = (TH1D*)closure_histos_PbPb_data->Get((TString)("h_reco_corr_cent"+cent[ibin-1]));
      h_reco_MC_q[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_q_cent"+cent[ibin-1]));
      h_reco_MC_g[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_g_cent"+cent[ibin-1]));
      h_reco_MC_up[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_up_cent"+cent[ibin-1]));
      h_reco_MC_upbar[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_upbar_cent"+cent[ibin-1]));
      h_reco_MC_d[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_d_cent"+cent[ibin-1]));
      h_reco_MC_dbar[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_dbar_cent"+cent[ibin-1]));

    }
/*    
    //normalizing eta spectra
    h_eta_MC[ibin]->Scale(1./h_eta_MC[ibin]->Integral());
    h_eta_MC_q[ibin]->Scale(1./h_eta_MC_q[ibin]->Integral());
    h_eta_MC_g[ibin]->Scale(1./h_eta_MC_g[ibin]->Integral());
    h_eta_data[ibin]->Scale(1./h_eta_data[ibin]->Integral());
    h_eta_MC_up[ibin]->Scale(1./h_eta_MC_up[ibin]->Integral());
    h_eta_MC_upbar[ibin]->Scale(1./h_eta_MC_upbar[ibin]->Integral());
    h_eta_MC_d[ibin]->Scale(1./h_eta_MC_d[ibin]->Integral());
    h_eta_MC_dbar[ibin]->Scale(1./h_eta_MC_dbar[ibin]->Integral());

    //normalizing jet spectra
    h_reco_MC[ibin]->Scale(1./h_reco_MC[ibin]->Integral());
    h_reco_data[ibin]->Scale(1./h_reco_data[ibin]->Integral());
    h_reco_MC_q[ibin]->Scale(1./h_reco_MC_q[ibin]->Integral());
    h_reco_MC_g[ibin]->Scale(1./h_reco_MC_g[ibin]->Integral());
    h_reco_MC_up[ibin]->Scale(1./h_reco_MC_up[ibin]->Integral());
    h_reco_MC_upbar[ibin]->Scale(1./h_reco_MC_upbar[ibin]->Integral());
    h_reco_MC_d[ibin]->Scale(1./h_reco_MC_d[ibin]->Integral());
    h_reco_MC_dbar[ibin]->Scale(1./h_reco_MC_dbar[ibin]->Integral());
*/
    h_trk_MC[ibin] = h_trk_corrpt_MC[ibin]-> ProjectionY(Form("h_trk_MC_%d",ibin),h_trk_corrpt_MC[ibin]->GetXaxis()->FindBin(pt_low),h_trk_corrpt_MC[ibin]->GetXaxis()->FindBin(pt_high),"");    
    h_bkgtrk_MC[ibin] = h_bkgtrk_corrpt_MC[ibin]-> ProjectionY(Form("h_bkgtrk_MC_%d",ibin),h_bkgtrk_corrpt_MC[ibin]->GetXaxis()->FindBin(pt_low),h_bkgtrk_corrpt_MC[ibin]->GetXaxis()->FindBin(pt_high),"");
    h_trk_data[ibin] = h_trk_corrpt_data[ibin]-> ProjectionY(Form("h_trk_data_%d",ibin),h_trk_corrpt_data[ibin]->GetXaxis()->FindBin(pt_low),h_trk_corrpt_data[ibin]->GetXaxis()->FindBin(pt_high),"");    
    h_bkgtrk_data[ibin] = h_bkgtrk_corrpt_data[ibin]-> ProjectionY(Form("h_bkgtrk_data_%d",ibin),h_bkgtrk_corrpt_data[ibin]->GetXaxis()->FindBin(pt_low),h_bkgtrk_corrpt_data[ibin]->GetXaxis()->FindBin(pt_high),"");

    h_trk_MC[ibin]->Scale(1./h_trk_MC[ibin]->Integral());
    h_trk_data[ibin]->Scale(1./h_trk_data[ibin]->Integral());
    h_bkgtrk_MC[ibin]->Scale(1./h_bkgtrk_MC[ibin]->Integral());
    h_bkgtrk_data[ibin]->Scale(1./h_bkgtrk_data[ibin]->Integral());

  }

  for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    for(int ibin=0;ibin<nCBins;ibin++){
      if(ibin==0){
  	    h_chg_refpt[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
  	    h_chg_refpt_q[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_q_cent"+cent[ibin]+"_trk"+trk[ibin2]));
  	    h_chg_refpt_g[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_g_cent"+cent[ibin]+"_trk"+trk[ibin2]));
  	    h_chg_refpt_u[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_u_cent"+cent[ibin]+"_trk"+trk[ibin2]));

  	    h_chg_corrpt[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
  	    h_chg_corrpt_q[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_q_cent"+cent[ibin]+"_trk"+trk[ibin2]));
  	    h_chg_corrpt_g[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_g_cent"+cent[ibin]+"_trk"+trk[ibin2]));
  	    h_chg_corrpt_u[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_u_cent"+cent[ibin]+"_trk"+trk[ibin2]));

        h_chg_corrpt_bkg[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_bkg_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
        h_chg_corrpt_data_bkg[ibin][ibin2] = (TH2F*)closure_histos_pp_data->Get((TString)("h_chg_corrpt_bkg_cent"+cent[ibin]+"_trk"+trk[ibin2]));

  	    h_chg_corrpt_data[ibin][ibin2] = (TH2F*)closure_histos_pp_data->Get((TString)("h_chg_corrpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 

        h_chg_corrpt_up[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_up_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_corrpt_down[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_d_cent"+cent[ibin]+"_trk"+trk[ibin2]));      
        h_chg_corrpt_upb[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_upbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_corrpt_downb[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_dbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));
      }
      else{
  	    h_chg_refpt[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2]))->Clone((TString)("h_chg_refpt_"+cent[ibin-1]+"_"+trk[ibin2])); 
  	    h_chg_refpt_q[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_q_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
  	    h_chg_refpt_g[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_g_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
  	    h_chg_refpt_u[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_u_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

  	    h_chg_corrpt[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
  	    h_chg_corrpt_q[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_q_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
  	    h_chg_corrpt_g[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_g_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
  	    h_chg_corrpt_u[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_u_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

        h_chg_corrpt_bkg[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_bkg_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
        h_chg_corrpt_data_bkg[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_chg_corrpt_bkg_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

  	    h_chg_corrpt_data[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_chg_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 

        h_chg_corrpt_up[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_up_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_corrpt_down[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_d_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_corrpt_upb[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_upbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_corrpt_downb[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_dbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));      
      }

        h_chg_ref[ibin][ibin2] = h_chg_refpt[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_%d_%d",ibin,ibin2),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref[ibin][ibin2]->Scale(1./(h_chg_ref[ibin][ibin2]->Integral()));
        h_chg_ref_q[ibin][ibin2] = h_chg_refpt_q[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_q_%d_%d",ibin,ibin2),h_chg_refpt_q[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_q[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_q[ibin][ibin2]->Scale(1./(h_chg_ref_q[ibin][ibin2]->Integral()));
        h_chg_ref_g[ibin][ibin2] = h_chg_refpt_g[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_g_%d_%d",ibin,ibin2),h_chg_refpt_g[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_g[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_g[ibin][ibin2]->Scale(1./(h_chg_ref_g[ibin][ibin2]->Integral()));
        h_chg_ref_u[ibin][ibin2] = h_chg_refpt_u[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_u_%d_%d",ibin,ibin2),h_chg_refpt_u[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_u[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_u[ibin][ibin2]->Scale(1./(h_chg_ref_u[ibin][ibin2]->Integral()));

        h_chg_reco[ibin][ibin2] = h_chg_corrpt[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_%d_%d",ibin,ibin2),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_reco_q[ibin][ibin2] = h_chg_corrpt_q[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_q_%d_%d",ibin,ibin2),h_chg_corrpt_q[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_q[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_reco_g[ibin][ibin2] = h_chg_corrpt_g[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_g_%d_%d",ibin,ibin2),h_chg_corrpt_g[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_g[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_reco_u[ibin][ibin2] = h_chg_corrpt_u[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_u_%d_%d",ibin,ibin2),h_chg_corrpt_u[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_u[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_reco_up[ibin][ibin2] = h_chg_corrpt_up[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_up_%d_%d",ibin,ibin2),h_chg_corrpt_up[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_up[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_reco_down[ibin][ibin2] = h_chg_corrpt_down[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_down_%d_%d",ibin,ibin2),h_chg_corrpt_down[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_down[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_reco_upb[ibin][ibin2] = h_chg_corrpt_upb[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_upb_%d_%d",ibin,ibin2),h_chg_corrpt_upb[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_upb[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_reco_downb[ibin][ibin2] = h_chg_corrpt_downb[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_downb_%d_%d",ibin,ibin2),h_chg_corrpt_downb[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_downb[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_data[ibin][ibin2] = h_chg_corrpt_data[ibin][ibin2]-> ProjectionY(Form("h_chg_data_%d_%d",ibin,ibin2),h_chg_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
/*
        //normalizing chg histos
        h_chg_reco[ibin][ibin2]->Scale(1./(h_chg_reco[ibin][ibin2]->Integral()));
        h_chg_reco_q[ibin][ibin2]->Scale(1./(h_chg_reco_q[ibin][ibin2]->Integral()));
        h_chg_reco_g[ibin][ibin2]->Scale(1./(h_chg_reco_g[ibin][ibin2]->Integral()));
        h_chg_reco_u[ibin][ibin2]->Scale(1./(h_chg_reco_u[ibin][ibin2]->Integral()));                        
        h_chg_reco_up[ibin][ibin2]->Scale(1./(h_chg_reco_up[ibin][ibin2]->Integral()));
        h_chg_reco_down[ibin][ibin2]->Scale(1./(h_chg_reco_down[ibin][ibin2]->Integral()));
        h_chg_reco_upb[ibin][ibin2]->Scale(1./(h_chg_reco_upb[ibin][ibin2]->Integral()));        
        h_chg_reco_downb[ibin][ibin2]->Scale(1./(h_chg_reco_downb[ibin][ibin2]->Integral()));
        h_chg_data[ibin][ibin2]->Scale(1./(h_chg_data[ibin][ibin2]->Integral()));
*/        
        sprintf(saythis,"h_chg_reco_data_refl_%d_%d",ibin,ibin2);
        h_chg_reco_data_refl[ibin][ibin2] = new TH1D(saythis,"",250,-2.5,2.5);
        h_chg_reco_data_refl[ibin][ibin2]->Sumw2();

        sprintf(saythis,"h_chg_reco_MC_refl_%d_%d",ibin,ibin2);
        h_chg_reco_MC_refl[ibin][ibin2] = new TH1D(saythis,"",250,-2.5,2.5);
        h_chg_reco_MC_refl[ibin][ibin2]->Sumw2();

        sprintf(saythis,"h_chg_reco_q_refl_%d_%d",ibin,ibin2);
        h_chg_reco_q_refl[ibin][ibin2] = new TH1D(saythis,"",250,-2.5,2.5);
        h_chg_reco_q_refl[ibin][ibin2]->Sumw2();

        sprintf(saythis,"h_chg_reco_g_refl_%d_%d",ibin,ibin2);
        h_chg_reco_g_refl[ibin][ibin2] = new TH1D(saythis,"",250,-2.5,2.5);
        h_chg_reco_g_refl[ibin][ibin2]->Sumw2();

        for(int ibin3=h_chg_reco_q_refl[ibin][ibin2]->GetXaxis()->GetNbins(); ibin3>0; ibin3--){
          Double_t nbins = h_chg_reco_q[ibin][ibin2]->GetXaxis()->GetNbins();
          //if(ibin3>nbins/2){
            h_chg_reco_q_refl[ibin][ibin2]->SetBinContent(ibin3,(h_chg_reco_q[ibin][ibin2]->GetBinContent(ibin3))-(h_chg_reco_q[ibin][ibin2]->GetBinContent(nbins-ibin3+1)));
            h_chg_reco_q_refl[ibin][ibin2]->SetBinError(ibin3,(h_chg_reco_q[ibin][ibin2]->GetBinError(ibin3))+(h_chg_reco_q[ibin][ibin2]->GetBinError(nbins-ibin3+1)));
            h_chg_reco_g_refl[ibin][ibin2]->SetBinContent(ibin3,(h_chg_reco_g[ibin][ibin2]->GetBinContent(ibin3))-(h_chg_reco_g[ibin][ibin2]->GetBinContent(nbins-ibin3+1)));
            h_chg_reco_g_refl[ibin][ibin2]->SetBinError(ibin3,(h_chg_reco_g[ibin][ibin2]->GetBinError(ibin3))+(h_chg_reco_g[ibin][ibin2]->GetBinError(nbins-ibin3+1)));
            h_chg_reco_MC_refl[ibin][ibin2]->SetBinContent(ibin3,(h_chg_reco[ibin][ibin2]->GetBinContent(ibin3))-(h_chg_reco[ibin][ibin2]->GetBinContent(nbins-ibin3+1)));
            h_chg_reco_MC_refl[ibin][ibin2]->SetBinError(ibin3,(h_chg_reco[ibin][ibin2]->GetBinError(ibin3))+(h_chg_reco[ibin][ibin2]->GetBinError(nbins-ibin3+1)));
            h_chg_reco_data_refl[ibin][ibin2]->SetBinContent(ibin3,(h_chg_data[ibin][ibin2]->GetBinContent(ibin3))-(h_chg_data[ibin][ibin2]->GetBinContent(nbins-ibin3+1)));
            h_chg_reco_data_refl[ibin][ibin2]->SetBinError(ibin3,(h_chg_data[ibin][ibin2]->GetBinError(ibin3))+(h_chg_data[ibin][ibin2]->GetBinError(nbins-ibin3+1)));
          //}
/*
          else{
            h_chg_reco_q_refl[ibin][ibin2]->SetBinContent(ibin3,);
            h_chg_reco_q_refl[ibin][ibin2]->SetBinError(ibin3,);
            h_chg_reco_g_refl[ibin][ibin2]->SetBinContent(ibin3,(h_chg_reco_g[ibin][ibin2]->GetBinContent(ibin3))-(h_chg_reco_g[ibin][ibin2]->GetBinContent(nbins-ibin3+1)));
            h_chg_reco_g_refl[ibin][ibin2]->SetBinError(ibin3,(h_chg_reco_g[ibin][ibin2]->GetBinError(ibin3))+(h_chg_reco_g[ibin][ibin2]->GetBinError(nbins-ibin3+1)));
            h_chg_reco_MC_refl[ibin][ibin2]->SetBinContent(ibin3,(h_chg_reco[ibin][ibin2]->GetBinContent(ibin3))-(h_chg_reco[ibin][ibin2]->GetBinContent(nbins-ibin3+1)));
            h_chg_reco_MC_refl[ibin][ibin2]->SetBinError(ibin3,(h_chg_reco[ibin][ibin2]->GetBinError(ibin3))+(h_chg_reco[ibin][ibin2]->GetBinError(nbins-ibin3+1)));
            h_chg_reco_data_refl[ibin][ibin2]->SetBinContent(ibin3,(h_chg_data[ibin][ibin2]->GetBinContent(ibin3))-(h_chg_data[ibin][ibin2]->GetBinContent(nbins-ibin3+1)));
            h_chg_reco_data_refl[ibin][ibin2]->SetBinError(ibin3,(h_chg_data[ibin][ibin2]->GetBinError(ibin3))+(h_chg_data[ibin][ibin2]->GetBinError(nbins-ibin3+1)));            
          } 
*/
        }
/*
        sprintf(saythis,"h_chg_reco_uubar_%d_%d",ibin,ibin2);
        h_chg_reco_uubar[ibin][ibin2] = new TH1D(saythis,"",250,-2.5,2.5);
        h_chg_reco_uubar[ibin][ibin2]->Sumw2();

        sprintf(saythis,"h_chg_reco_ddbar_%d_%d",ibin,ibin2);
        h_chg_reco_ddbar[ibin][ibin2] = new TH1D(saythis,"",250,-2.5,2.5);
        h_chg_reco_ddbar[ibin][ibin2]->Sumw2();

        for(int ibin3=1; ibin3<(h_chg_reco_uubar[ibin][ibin2]->GetXaxis()->GetNbins()+1); ibin3++){
          h_chg_reco_uubar[ibin][ibin2]->SetBinContent(ibin3,(h_chg_reco_up[ibin][ibin2]->GetBinContent(ibin3))+(h_chg_reco_upb[ibin][ibin2]->GetBinContent(ibin3)));
          h_chg_reco_uubar[ibin][ibin2]->SetBinError(ibin3,(h_chg_reco_up[ibin][ibin2]->GetBinError(ibin3))+(h_chg_reco_upb[ibin][ibin2]->GetBinError(ibin3)));
          h_chg_reco_ddbar[ibin][ibin2]->SetBinContent(ibin3,(h_chg_reco_down[ibin][ibin2]->GetBinContent(ibin3))+(h_chg_reco_downb[ibin][ibin2]->GetBinContent(ibin3)));
          h_chg_reco_ddbar[ibin][ibin2]->SetBinError(ibin3,(h_chg_reco_down[ibin][ibin2]->GetBinError(ibin3))+(h_chg_reco_downb[ibin][ibin2]->GetBinError(ibin3)));
        } 
        for(int ibin3=1; ibin3<(h_chg_reco_uubar[ibin][ibin2]->GetXaxis()->GetNbins()+1); ibin3++){
          Double_t nbins = h_chg_reco_uubar[ibin][ibin2]->GetXaxis()->GetNbins();
          if(ibin3<=nbins/2){
            h_chg_reco_uubar[ibin][ibin2]->SetBinContent(ibin3,((h_chg_reco_uubar[ibin][ibin2]->GetBinContent(ibin3))+(h_chg_reco_uubar[ibin][ibin2]->GetBinContent(nbins-ibin3)))/2.);
            h_chg_reco_uubar[ibin][ibin2]->SetBinError(ibin3,((h_chg_reco_uubar[ibin][ibin2]->GetBinError(ibin3))+(h_chg_reco_uubar[ibin][ibin2]->GetBinError(ibin3)))/2.);
            h_chg_reco_ddbar[ibin][ibin2]->SetBinContent(ibin3,((h_chg_reco_ddbar[ibin][ibin2]->GetBinContent(ibin3))+(h_chg_reco_ddbar[ibin][ibin2]->GetBinContent(ibin3)))/2.);
            h_chg_reco_ddbar[ibin][ibin2]->SetBinError(ibin3,((h_chg_reco_ddbar[ibin][ibin2]->GetBinError(ibin3))+(h_chg_reco_ddbar[ibin][ibin2]->GetBinError(ibin3)))/2.);
          }
          else{
            h_chg_reco_uubar[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_uubar[ibin][ibin2]->GetBinContent(nbins-ibin3+1));
            h_chg_reco_uubar[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_uubar[ibin][ibin2]->GetBinError(nbins-ibin3+1));
            h_chg_reco_ddbar[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_ddbar[ibin][ibin2]->GetBinContent(nbins-ibin3+1));
            h_chg_reco_ddbar[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_ddbar[ibin][ibin2]->GetBinError(nbins-ibin3+1));
          }
        }
        h_chg_reco_uubar[ibin][ibin2]->Scale(1./(h_chg_reco_uubar[ibin][ibin2]->Integral()));
        h_chg_reco_ddbar[ibin][ibin2]->Scale(1./(h_chg_reco_ddbar[ibin][ibin2]->Integral()));
*/

        h_chg_ref_avg[ibin][ibin2] = h_chg_refpt[ibin][ibin2]-> ProfileX(Form("h_chg_ref_avg_%d_%d",ibin,ibin2),h_chg_refpt[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_refpt[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_ref_q_avg[ibin][ibin2] = h_chg_refpt_q[ibin][ibin2]-> ProfileX(Form("h_chg_ref_avg_q_%d_%d",ibin,ibin2),h_chg_refpt_q[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_refpt_q[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_ref_g_avg[ibin][ibin2] = h_chg_refpt_g[ibin][ibin2]-> ProfileX(Form("h_chg_ref_avg_g_%d_%d",ibin,ibin2),h_chg_refpt_g[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_refpt_g[ibin][ibin2]->GetYaxis()->FindBin(10),"");

        h_chg_reco_avg[ibin][ibin2] = h_chg_corrpt[ibin][ibin2]-> ProfileX(Form("h_chg_reco_avg_%d_%d",ibin,ibin2),h_chg_corrpt[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_reco_q_avg[ibin][ibin2] = h_chg_corrpt_q[ibin][ibin2]-> ProfileX(Form("h_chg_reco_q_avg_%d_%d",ibin,ibin2),h_chg_corrpt_q[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_q[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_reco_g_avg[ibin][ibin2] = h_chg_corrpt_g[ibin][ibin2]-> ProfileX(Form("h_chg_reco_g_avg_%d_%d",ibin,ibin2),h_chg_corrpt_g[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_g[ibin][ibin2]->GetYaxis()->FindBin(10),"");

        h_chg_reco_up_avg[ibin][ibin2] = h_chg_corrpt_up[ibin][ibin2]-> ProfileX(Form("h_chg_reco_up_avg_%d_%d",ibin,ibin2),h_chg_corrpt_up[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_up[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_reco_down_avg[ibin][ibin2] = h_chg_corrpt_down[ibin][ibin2]-> ProfileX(Form("h_chg_reco_down_avg_%d_%d",ibin,ibin2),h_chg_corrpt_down[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_down[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_reco_upb_avg[ibin][ibin2] = h_chg_corrpt_upb[ibin][ibin2]-> ProfileX(Form("h_chg_reco_upb_avg_%d_%d",ibin,ibin2),h_chg_corrpt_upb[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_upb[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_reco_downb_avg[ibin][ibin2] = h_chg_corrpt_downb[ibin][ibin2]-> ProfileX(Form("h_chg_reco_downb_avg_%d_%d",ibin,ibin2),h_chg_corrpt_downb[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_downb[ibin][ibin2]->GetYaxis()->FindBin(10),"");

        h_chg_data_avg[ibin][ibin2] = h_chg_corrpt_data[ibin][ibin2]-> ProfileX(Form("h_chg_data_avg_%d_%d",ibin,ibin2),h_chg_corrpt_data[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_data[ibin][ibin2]->GetYaxis()->FindBin(10),"");

        h_chg_reco_bkg_avg[ibin][ibin2] = h_chg_corrpt_bkg[ibin][ibin2]-> ProfileX(Form("h_chg_bkg_reco_avg_%d_%d",ibin,ibin2),h_chg_corrpt_bkg[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_bkg[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_data_bkg_avg[ibin][ibin2] = h_chg_corrpt_data_bkg[ibin][ibin2]-> ProfileX(Form("h_chg_bkg_data_avg_%d_%d",ibin,ibin2),h_chg_corrpt_data_bkg[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_data_bkg[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        
        h_chg_reco_bkg[ibin][ibin2] = h_chg_corrpt_bkg[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_bkg_%d_%d",ibin,ibin2),h_chg_corrpt_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_reco_bkg[ibin][ibin2]->Scale(1./(h_chg_reco_bkg[ibin][ibin2]->Integral()));
        h_chg_data_bkg[ibin][ibin2] = h_chg_corrpt_data_bkg[ibin][ibin2]-> ProjectionY(Form("h_chg_data_bkg_%d_%d",ibin,ibin2),h_chg_corrpt_data_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_data_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_data_bkg[ibin][ibin2]->Scale(1./(h_chg_data_bkg[ibin][ibin2]->Integral()));
    }
  }

  f_chg_ref[0][0] = new TF1("f_chg_ref_cent0",ftotal_ref_cent0,-10.,10.,1);
  f_chg_ref[1][0] = new TF1("f_chg_ref_cent1",ftotal_ref_cent1,-10.,10.,1);
  f_chg_ref[2][0] = new TF1("f_chg_ref_cent2",ftotal_ref_cent2,-10.,10.,1);
  f_chg_ref[3][0] = new TF1("f_chg_ref_cent3",ftotal_ref_cent3,-10.,10.,1);
  f_chg_ref[4][0] = new TF1("f_chg_ref_cent4",ftotal_ref_cent4,-10.,10.,1);

  f_chg_reco[0][0] = new TF1("f_chg_reco_cent0",ftotal_reco_cent0,-10.,10.,1);
  f_chg_reco[1][0] = new TF1("f_chg_reco_cent1",ftotal_reco_cent1,-10.,10.,1);
  f_chg_reco[2][0] = new TF1("f_chg_reco_cent2",ftotal_reco_cent2,-10.,10.,1);
  f_chg_reco[3][0] = new TF1("f_chg_reco_cent3",ftotal_reco_cent3,-10.,10.,1);
  f_chg_reco[4][0] = new TF1("f_chg_reco_cent4",ftotal_reco_cent4,-10.,10.,1);

  f_eta_reco[0] = new TF1("f_eta_reco_cent0",feta_reco_cent0,-1.5,1.5,1);
  f_eta_reco[1] = new TF1("f_eta_reco_cent1",feta_reco_cent1,-1.5,1.5,1);
  f_eta_reco[2] = new TF1("f_eta_reco_cent2",feta_reco_cent2,-1.5,1.5,1);
  f_eta_reco[3] = new TF1("f_eta_reco_cent3",feta_reco_cent3,-1.5,1.5,1);
  f_eta_reco[4] = new TF1("f_eta_reco_cent4",feta_reco_cent4,-1.5,1.5,1);
      
  f_chg_data[0][0] = new TF1("f_chg_data_cent0",ftotal_reco_cent0,-10.,10.,1);
  f_chg_data[0][1] = new TF1("f_chg_data_cent0_trkpt2",ftotal_reco_cent0_pt2,-10.,10.,1);
  f_chg_data[0][2] = new TF1("f_chg_data_cent0_trkpt4",ftotal_reco_cent0_pt4,-10.,10.,1);
  f_chg_data[0][3] = new TF1("f_chg_data_cent0_trkpt5",ftotal_reco_cent0_pt5,-10.,10.,1);
  f_chg_data[1][0] = new TF1("f_chg_data_cent0",ftotal_reco_cent1,-10.,10.,1);
  f_chg_data[1][1] = new TF1("f_chg_data_cent0_trkpt2",ftotal_reco_cent1_pt2,-10.,10.,1);
  f_chg_data[1][2] = new TF1("f_chg_data_cent0_trkpt4",ftotal_reco_cent1_pt4,-10.,10.,1);
  f_chg_data[1][3] = new TF1("f_chg_data_cent0_trkpt5",ftotal_reco_cent1_pt5,-10.,10.,1);
  f_chg_data[2][0] = new TF1("f_chg_data_cent0",ftotal_reco_cent2,-10.,10.,1);
  f_chg_data[2][1] = new TF1("f_chg_data_cent0_trkpt2",ftotal_reco_cent2_pt2,-10.,10.,1);
  f_chg_data[2][2] = new TF1("f_chg_data_cent0_trkpt4",ftotal_reco_cent2_pt4,-10.,10.,1);
  f_chg_data[2][3] = new TF1("f_chg_data_cent0_trkpt5",ftotal_reco_cent2_pt5,-10.,10.,1);
  f_chg_data[3][0] = new TF1("f_chg_data_cent0",ftotal_reco_cent3,-10.,10.,1);
  f_chg_data[3][1] = new TF1("f_chg_data_cent0_trkpt2",ftotal_reco_cent3_pt2,-10.,10.,1);
  f_chg_data[3][2] = new TF1("f_chg_data_cent0_trkpt4",ftotal_reco_cent3_pt4,-10.,10.,1);
  f_chg_data[3][3] = new TF1("f_chg_data_cent0_trkpt5",ftotal_reco_cent3_pt5,-10.,10.,1);
  f_chg_data[4][0] = new TF1("f_chg_data_cent0",ftotal_reco_cent4,-10.,10.,1);
  f_chg_data[4][1] = new TF1("f_chg_data_cent0_trkpt2",ftotal_reco_cent4_pt2,-10.,10.,1);
  f_chg_data[4][2] = new TF1("f_chg_data_cent0_trkpt4",ftotal_reco_cent4_pt4,-10.,10.,1);
  f_chg_data[4][3] = new TF1("f_chg_data_cent0_trkpt5",ftotal_reco_cent4_pt5,-10.,10.,1);

  f_chg_udg_data[0][3] = new TF1("f_chg_udg_data_cent0_trkpt5",f_udg_reco_cent0_pt5,-10.,10.,2);
  f_chg_udg_data[1][3] = new TF1("f_chg_udg_data_cent1_trkpt5",f_udg_reco_cent1_pt5,-10.,10.,2);
  f_chg_udg_data[2][3] = new TF1("f_chg_udg_data_cent2_trkpt5",f_udg_reco_cent2_pt5,-10.,10.,2);
  f_chg_udg_data[3][3] = new TF1("f_chg_udg_data_cent3_trkpt5",f_udg_reco_cent3_pt5,-10.,10.,2);
  f_chg_udg_data[4][3] = new TF1("f_chg_udg_data_cent4_trkpt5",f_udg_reco_cent4_pt5,-10.,10.,2);

  for(int ibin=0;ibin<nCBins;ibin++){
    f_chg_ref[ibin][0]->SetLineColor(kBlack);
    f_chg_reco[ibin][0]->SetLineColor(kBlack);
    f_eta_reco[ibin]->SetLineColor(kBlack);
    for(int ibin2=0;ibin2<ntrkbins;ibin2++){
        f_chg_data[ibin][ibin2]->SetLineColor(kGreen);
    }
  }

  TLine *tl1 = new TLine(-1.5,1.,1.5,1.);
  tl1->SetLineStyle(2);

  TLine *tl2 = new TLine(120.,0.,300.,0.);
  tl2->SetLineStyle(2);
/*
  TCanvas *c_chg_MC_ref = new TCanvas("c_chg_MC_ref","",600,1200);
  c_chg_MC_ref->Divide(1,2);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins;ibin++){
  	for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    c_chg_MC_ref->cd(ibin+1);   
    h_chg_ref[ibin][ibin2]->GetXaxis()->SetTitle("gen jet chg");
    h_chg_ref[ibin][ibin2]->GetYaxis()->SetTitle("# of jets");
    h_chg_ref[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_ref[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_ref[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
    h_chg_ref[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_ref[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_ref[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_ref[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
    h_chg_ref[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    //h_chg_ref[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.4);
    h_chg_ref[ibin][ibin2]->SetLineColor(kBlack); h_chg_ref[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_ref[ibin][ibin2]->SetMarkerStyle(4); 
    h_chg_ref_q[ibin][ibin2]->SetLineColor(kBlue); h_chg_ref_q[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_ref_g[ibin][ibin2]->SetLineColor(kRed); h_chg_ref_g[ibin][ibin2]->SetMarkerColor(kRed);
    if(ibin2==0 && ibin==0){
      h_chg_ref[ibin][ibin2]->Draw("e0 same");
      h_chg_ref_q[ibin][ibin2]->Draw("e0 same");
      h_chg_ref_g[ibin][ibin2]->Draw("e0 same");
    }

    if(ibin2==0){ 
    h_chg_ref[ibin][0]->Fit(f_chg_ref[ibin][0],"Q N M R","sames",-2.,2.);

    c_chg_MC_ref->cd(ibin+2);   
    h_chg_ref_ratio[ibin][0] = (TH1D*)h_chg_ref[ibin][0]->Clone("h_chg_ref_ratio");
    h_chg_ref_ratio[ibin][0]->Divide(f_chg_ref[ibin][0]);
    h_chg_ref_ratio[ibin][0]->GetXaxis()->SetTitle("gen jet charge");
    h_chg_ref_ratio[ibin][0]->GetYaxis()->SetTitle("Gen / Fit");
    h_chg_ref_ratio[ibin][0]->GetYaxis()->SetNdivisions(505);
    h_chg_ref_ratio[ibin][0]->GetYaxis()->SetTitleSize(0.05);
    h_chg_ref_ratio[ibin][0]->GetYaxis()->SetLabelSize(0.09);
    h_chg_ref_ratio[ibin][0]->GetXaxis()->CenterTitle();
    h_chg_ref_ratio[ibin][0]->GetXaxis()->SetNdivisions(505);
    h_chg_ref_ratio[ibin][0]->GetXaxis()->SetTitleSize(0.05);
    h_chg_ref_ratio[ibin][0]->GetXaxis()->SetLabelSize(0.09);
    h_chg_ref_ratio[ibin][0]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_chg_ref_ratio[ibin][0]->GetYaxis()->SetRangeUser(0.,2.);
    h_chg_ref_ratio[ibin][0]->SetLineColor(kBlack); h_chg_ref_ratio[ibin][0]->SetMarkerColor(kBlack);
    h_chg_ref_ratio[ibin][0]->SetMarkerStyle(4); h_chg_ref_ratio[ibin][0]->Draw("e0 same");
    tl1->Draw("same");
    }}
  }
*/
/*
  TCanvas *c_chg_MC_reco = new TCanvas("c_chg_MC_reco","",1500,600);
  c_chg_MC_reco->Divide(5,2);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins;ibin++){
    for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    if(ibin==0) c_chg_MC_reco->cd(1);
    else c_chg_MC_reco->cd(6-ibin);   
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetTitle("reco jet chg");
    h_chg_reco[ibin][ibin2]->GetYaxis()->SetTitle("# of jets");
    h_chg_reco[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_reco[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_reco[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
    h_chg_reco[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    //h_chg_reco[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);
    h_chg_reco[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_reco[ibin][ibin2]->SetMarkerStyle(4);
    h_chg_reco_q[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_q[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_g[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_g[ibin][ibin2]->SetMarkerColor(kRed);
    if(ibin2==3){
      h_chg_reco[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_q[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_g[ibin][ibin2]->Draw("e0 same");
    }

    if(ibin2==3){
      h_chg_reco[ibin][ibin2]->Fit(f_chg_data[ibin][ibin2],"Q N M R","sames",-2.,2.);
      //cout<<f_chg_data[ibin][ibin2]->GetChisquare()/f_chg_data[ibin][ibin2]->GetNDF()<<endl;      
    }

    if(ibin2==3){
    if(ibin==0) c_chg_MC_reco->cd(6);
    else c_chg_MC_reco->cd(11-ibin);   
    h_chg_reco_ratio[ibin][ibin2] = (TH1D*)h_chg_reco[ibin][ibin2]->Clone("h_chg_reco_ratio");
    h_chg_reco_ratio[ibin][ibin2]->Divide(f_chg_data[ibin][ibin2]);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetTitle("reco jet charge");
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetTitle("Reco / Fit");
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetLabelSize(0.09);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetLabelSize(0.09);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,2.);
    h_chg_reco_ratio[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio[ibin][ibin2]->SetMarkerColor(kBlack);
    h_chg_reco_ratio[ibin][ibin2]->SetMarkerStyle(1); h_chg_reco_ratio[ibin][ibin2]->Draw("e0 same");
    tl1->Draw("same");
    }}
  }
*/
  TCanvas *c_chg_udg_MC_reco = new TCanvas("c_chg_udg_MC_reco","",1500,600);
  c_chg_udg_MC_reco->Divide(5,2);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins;ibin++){
    for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    if(ibin==0) c_chg_udg_MC_reco->cd(1);
    else c_chg_udg_MC_reco->cd(6-ibin);   
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetTitle("reco jet chg");
    h_chg_reco[ibin][ibin2]->GetYaxis()->SetTitle("# of jets");
    h_chg_reco[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_reco[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_reco[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
    h_chg_reco[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_chg_reco[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);
    h_chg_reco[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_reco[ibin][ibin2]->SetMarkerStyle(4);
    h_chg_reco_up[ibin][ibin2]->SetMarkerStyle(22); h_chg_reco_up[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_up[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_down[ibin][ibin2]->SetMarkerStyle(23); h_chg_reco_down[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_down[ibin][ibin2]->SetMarkerColor(kBlue);
    //h_chg_reco_uubar[ibin][ibin2]->SetMarkerStyle(22); h_chg_reco_uubar[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_uubar[ibin][ibin2]->SetMarkerColor(kBlue);
    //h_chg_reco_ddbar[ibin][ibin2]->SetMarkerStyle(23); h_chg_reco_ddbar[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_ddbar[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_g[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_g[ibin][ibin2]->SetMarkerColor(kRed);
    if(ibin2==3){
      h_chg_reco[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_up[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_down[ibin][ibin2]->Draw("e0 same");
      //h_chg_reco_uubar[ibin][ibin2]->Rebin(10); h_chg_reco_uubar[ibin][ibin2]->Draw("e0 same");
      //h_chg_reco_ddbar[ibin][ibin2]->Rebin(10); h_chg_reco_ddbar[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_g[ibin][ibin2]->Draw("e0 same");
    }
    if(do_fit==true){
    if(ibin2==3){
      h_chg_reco[ibin][ibin2]->Fit(f_chg_udg_data[ibin][ibin2],"N M R","sames",-2.,2.);
      cout<<f_chg_udg_data[ibin][ibin2]->GetChisquare()/f_chg_udg_data[ibin][ibin2]->GetNDF()<<endl;      
    }

    if(ibin2==3){
    if(ibin==0) c_chg_udg_MC_reco->cd(6);
    else c_chg_udg_MC_reco->cd(11-ibin);   
    h_chg_reco_ratio[ibin][ibin2] = (TH1D*)h_chg_reco[ibin][ibin2]->Clone("h_chg_reco_ratio");
    h_chg_reco_ratio[ibin][ibin2]->Divide(f_chg_udg_data[ibin][ibin2]);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetTitle("reco jet charge");
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetTitle("Reco / Fit");
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetLabelSize(0.09);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetLabelSize(0.09);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,2.);
    h_chg_reco_ratio[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio[ibin][ibin2]->SetMarkerColor(kBlack);
    h_chg_reco_ratio[ibin][ibin2]->SetMarkerStyle(1); h_chg_reco_ratio[ibin][ibin2]->Draw("e0 same");
    tl1->Draw("same");
    }}}
  }
/*
  TH1D *h_reco_q_cent_ratio[nCBins][ntrkbins];
  TH1D *h_reco_g_cent_ratio[nCBins][ntrkbins];

  for(int ibin=0;ibin<nCBins;ibin++){
    for(int ibin2=0;ibin2<ntrkbins;ibin2++){
      h_chg_reco_q[ibin][ibin2]->Rebin(10);
      h_chg_reco_g[ibin][ibin2]->Rebin(10);
    }
  }
 
  TCanvas *c_chg_qg_ratio_reco = new TCanvas("c_chg_qg_ratio_reco","",600,300);
  c_chg_qg_ratio_reco->Divide(2,1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins-1;ibin++){
    c_chg_qg_ratio_reco->cd(1);    
    h_reco_q_cent_ratio[ibin][3] = (TH1D*)h_chg_reco_q[4][3]->Clone((TString)("h_reco_q_ratio_cent_"+cent[ibin]));
    h_reco_q_cent_ratio[ibin][3]->Divide(h_chg_reco_q[ibin][3]);
    h_reco_q_cent_ratio[ibin][3]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_reco_q_cent_ratio[ibin][3]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_reco_q_cent_ratio[ibin][3]->SetMarkerStyle(ibin);
    h_reco_q_cent_ratio[ibin][3]->Draw("e0 same");
    tl1->Draw("same");

    c_chg_qg_ratio_reco->cd(2);    
    h_reco_g_cent_ratio[ibin][3] = (TH1D*)h_chg_reco_g[4][3]->Clone((TString)("h_reco_g_ratio_cent_"+cent[ibin]));
    h_reco_g_cent_ratio[ibin][3]->Divide(h_chg_reco_g[ibin][3]);
    h_reco_g_cent_ratio[ibin][3]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_reco_g_cent_ratio[ibin][3]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_reco_g_cent_ratio[ibin][3]->SetMarkerStyle(ibin);
    h_reco_g_cent_ratio[ibin][3]->Draw("e0 same");
    tl1->Draw("same");
  }
*/
/*  
  TCanvas *c_chg_MC_data = new TCanvas("c_chg_MC_data","",1500,600);
  c_chg_MC_data->Divide(5,2);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins;ibin++){
    for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    if(ibin==0) c_chg_MC_data->cd(1);
    else c_chg_MC_data->cd(6-ibin);    	
    h_chg_data[ibin][ibin2]->GetXaxis()->SetTitle("data jet chg");
    h_chg_data[ibin][ibin2]->GetYaxis()->SetTitle("# of jets");
    h_chg_data[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_data[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_data[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
    h_chg_data[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_data[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_data[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_data[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
    h_chg_data[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_chg_data[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);
    h_chg_data[ibin][ibin2]->SetLineColor(kGreen); h_chg_data[ibin][ibin2]->SetMarkerColor(kGreen); h_chg_data[ibin][ibin2]->SetMarkerStyle(4);
    h_chg_reco_q[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_q[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_g[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_g[ibin][ibin2]->SetMarkerColor(kRed);
    
    if(ibin2==3){
      h_chg_data[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_q[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_g[ibin][ibin2]->Draw("e0 same");
    }
    if(ibin2==3){
      h_chg_data[ibin][ibin2]->Fit(f_chg_data[ibin][ibin2],"Q N M R","sames",-2.,2.);
      //cout<<f_chg_data[ibin][ibin2]->GetChisquare()/f_chg_data[ibin][ibin2]->GetNDF()<<endl;      
    }
    if(ibin2==3){
    if(ibin==0) c_chg_MC_data->cd(6);
    else c_chg_MC_data->cd(11-ibin);   
    h_chg_data_ratio[ibin][ibin2] = (TH1D*)h_chg_data[ibin][ibin2]->Clone("h_chg_data_ratio");
    h_chg_data_ratio[ibin][ibin2]->Divide(f_chg_data[ibin][ibin2]);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetTitle("data jet charge");
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetTitle("Data / Fit");
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetLabelSize(0.09);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetLabelSize(0.09);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,2.);
    h_chg_data_ratio[ibin][ibin2]->SetLineColor(kBlack); h_chg_data_ratio[ibin][ibin2]->SetMarkerColor(kBlack);
    h_chg_data_ratio[ibin][ibin2]->SetMarkerStyle(4); h_chg_data_ratio[ibin][ibin2]->Draw("e0 same");
    tl1->Draw("same");
    }}
  }
*/
  TCanvas *c_chg_udg_MC_data = new TCanvas("c_chg_udg_MC_data","",1500,600);
  c_chg_udg_MC_data->Divide(5,2);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins;ibin++){
    for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    if(ibin==0) c_chg_udg_MC_data->cd(1);
    else c_chg_udg_MC_data->cd(6-ibin);   
    h_chg_data[ibin][ibin2]->GetXaxis()->SetTitle("data jet chg");
    h_chg_data[ibin][ibin2]->GetYaxis()->SetTitle("# of jets");
    h_chg_data[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_data[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_data[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
    h_chg_data[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_data[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_data[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_data[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
    h_chg_data[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_chg_data[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);
    h_chg_data[ibin][ibin2]->SetLineColor(kBlack); h_chg_data[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_data[ibin][ibin2]->SetMarkerStyle(4);
    h_chg_reco_up[ibin][ibin2]->SetMarkerStyle(22); h_chg_reco_up[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_up[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_down[ibin][ibin2]->SetMarkerStyle(23); h_chg_reco_down[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_down[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_g[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_g[ibin][ibin2]->SetMarkerColor(kRed);
    if(ibin2==3){
      h_chg_data[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_up[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_down[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_g[ibin][ibin2]->Draw("e0 same");
    }
    if(do_fit==true){
    if(ibin2==3){
      h_chg_data[ibin][ibin2]->Fit(f_chg_udg_data[ibin][ibin2],"Q N M R","sames",-2.,2.);
      //cout<<f_chg_udg_data[ibin][ibin2]->GetChisquare()/f_chg_udg_data[ibin][ibin2]->GetNDF()<<endl;      
    }

    if(ibin2==3){
    if(ibin==0) c_chg_udg_MC_data->cd(6);
    else c_chg_udg_MC_data->cd(11-ibin);   
    h_chg_data_ratio[ibin][ibin2] = (TH1D*)h_chg_data[ibin][ibin2]->Clone("h_chg_data_ratio");
    h_chg_data_ratio[ibin][ibin2]->Divide(f_chg_udg_data[ibin][ibin2]);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetTitle("data jet charge");
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetTitle("data / Fit");
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetLabelSize(0.09);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetLabelSize(0.09);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,2.);
    h_chg_data_ratio[ibin][ibin2]->SetLineColor(kBlack); h_chg_data_ratio[ibin][ibin2]->SetMarkerColor(kBlack);
    h_chg_data_ratio[ibin][ibin2]->SetMarkerStyle(1); h_chg_data_ratio[ibin][ibin2]->Draw("e0 same");
    tl1->Draw("same");
    }}}
  }

  TCanvas *c_chg_MC_data_stk = new TCanvas("c_chg_MC_data_stk","",1800,1200);
  c_chg_MC_data_stk->Divide(3,2);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin=0;ibin<1;ibin++){
    for(int ibin2=1;ibin2<ntrkbins;ibin2++){
    c_chg_MC_data_stk->cd(ibin2);   
    h_chg_data[ibin][ibin2]->GetXaxis()->SetTitle("data jet chg");
    h_chg_data[ibin][ibin2]->GetYaxis()->SetTitle("# of jets");
    h_chg_data[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_data[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_data[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
    h_chg_data[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_data[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_data[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_data[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
    h_chg_data[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_chg_data[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);
    h_chg_data[ibin][ibin2]->SetLineColor(kGreen); h_chg_data[ibin][ibin2]->SetMarkerColor(kGreen);
    h_chg_data[ibin][ibin2]->SetMarkerStyle(4); h_chg_data[ibin][ibin2]->Draw("e0 same");
    h_chg_reco_q[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_q[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_q[ibin][ibin2]->Draw("e0 same");
    h_chg_reco_g[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_g[ibin][ibin2]->SetMarkerColor(kRed);
    h_chg_reco_g[ibin][ibin2]->Draw("e0 same");

    h_chg_data[ibin][ibin2]->Fit(f_chg_data[ibin][ibin2],"Q N M R","sames",-2.,2.);

    c_chg_MC_data_stk->cd(ibin2+3);   
    h_chg_data_ratio[ibin][ibin2] = (TH1D*)h_chg_data[ibin][ibin2]->Clone("h_chg_data_ratio");
    h_chg_data_ratio[ibin][ibin2]->Divide(f_chg_data[ibin][ibin2]);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetTitle("data jet charge");
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetTitle("Data / Fit");
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetLabelSize(0.09);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetLabelSize(0.09);
    h_chg_data_ratio[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_chg_data_ratio[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,2.);
    h_chg_data_ratio[ibin][ibin2]->SetLineColor(kBlack); h_chg_data_ratio[ibin][ibin2]->SetMarkerColor(kBlack);
    h_chg_data_ratio[ibin][ibin2]->SetMarkerStyle(4); h_chg_data_ratio[ibin][ibin2]->Draw("e0 same");
    tl1->Draw("same");

  }
  }

  TCanvas *c_chg_MC_stk = new TCanvas("c_chg_MC_stk","",1800,1200);
  c_chg_MC_stk->Divide(3,2);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin=0;ibin<1;ibin++){
    for(int ibin2=1;ibin2<ntrkbins;ibin2++){
    c_chg_MC_stk->cd(ibin2);   
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetTitle("reco jet chg");
    h_chg_reco[ibin][ibin2]->GetYaxis()->SetTitle("# of jets");
    h_chg_reco[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_reco[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_reco[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
    h_chg_reco[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
    h_chg_reco[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    //h_chg_reco[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.4);
    h_chg_reco[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco[ibin][ibin2]->SetMarkerColor(kBlack);
    h_chg_reco[ibin][ibin2]->SetMarkerStyle(4); h_chg_reco[ibin][ibin2]->Draw("e0 same");
    h_chg_reco_q[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_q[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_q[ibin][ibin2]->Draw("e0 same");
    h_chg_reco_g[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_g[ibin][ibin2]->SetMarkerColor(kRed);
    h_chg_reco_g[ibin][ibin2]->Draw("e0 same");

    h_chg_reco[ibin][ibin2]->Fit(f_chg_data[ibin][ibin2],"Q N M R","sames",-2.,2.);

    c_chg_MC_stk->cd(ibin2+3);   
    h_chg_reco_ratio[ibin][ibin2] = (TH1D*)h_chg_reco[ibin][ibin2]->Clone("h_chg_reco_ratio");
    h_chg_reco_ratio[ibin][ibin2]->Divide(f_chg_data[ibin][ibin2]);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetTitle("reco jet charge");
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetTitle("Inc / Fit");
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetLabelSize(0.09);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetLabelSize(0.09);
    h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,2.);
    h_chg_reco_ratio[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio[ibin][ibin2]->SetMarkerColor(kBlack);
    h_chg_reco_ratio[ibin][ibin2]->SetMarkerStyle(4); h_chg_reco_ratio[ibin][ibin2]->Draw("e0 same");
    tl1->Draw("same");

  }
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
  TCanvas *c_eta = new TCanvas("c_eta","",1500,700);
  c_eta->Divide(5,2,0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0) c_eta->cd(ibin+1);
    else c_eta->cd(6-ibin);   
    h_eta_MC[ibin]->GetXaxis()->SetTitle("jet #eta");
    h_eta_MC[ibin]->GetYaxis()->SetNdivisions(505);
    h_eta_MC[ibin]->GetYaxis()->SetTitleSize(0.05);
    h_eta_MC[ibin]->GetYaxis()->SetLabelSize(0.07);
    h_eta_MC[ibin]->GetXaxis()->CenterTitle();
    h_eta_MC[ibin]->GetXaxis()->SetNdivisions(505);
    h_eta_MC[ibin]->GetXaxis()->SetTitleSize(0.05);
    h_eta_MC[ibin]->GetXaxis()->SetLabelSize(0.07);
    h_eta_MC[ibin]->GetXaxis()->SetRangeUser(-1.5,1.52);
    h_eta_MC[ibin]->GetYaxis()->SetRangeUser(0.02,0.12);
    h_eta_MC[ibin]->SetLineColor(kBlack); h_eta_MC[ibin]->SetMarkerColor(kBlack);
    h_eta_MC[ibin]->SetMarkerStyle(4); h_eta_MC[ibin]->Draw("e0 same");
    h_eta_data[ibin]->SetLineColor(kGreen); h_eta_data[ibin]->SetMarkerColor(kGreen);
    h_eta_data[ibin]->SetMarkerStyle(4); h_eta_data[ibin]->Draw("e0 same");

    if(ibin>0){
      TString tmp1 = cent_tag[ibin-1];
      tx1 = new TLatex(); tx1->SetTextSize(.12);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }

    if(ibin==0) c_eta->cd(ibin+6);
    else c_eta->cd(11-ibin);   
    h_eta_MC_data_ratio[ibin] = (TH1D*)h_eta_MC[ibin]->Clone("h_eta_MC_data_ratio");
    h_eta_MC_data_ratio[ibin]->Divide(h_eta_data[ibin]);
    h_eta_MC_data_ratio[ibin]->GetYaxis()->SetNdivisions(505);
    h_eta_MC_data_ratio[ibin]->GetYaxis()->SetTitleSize(0.05);
    h_eta_MC_data_ratio[ibin]->GetYaxis()->SetLabelSize(0.07);
    h_eta_MC_data_ratio[ibin]->GetXaxis()->CenterTitle();
    h_eta_MC_data_ratio[ibin]->GetXaxis()->SetNdivisions(505);
    h_eta_MC_data_ratio[ibin]->GetXaxis()->SetTitleSize(0.05);
    h_eta_MC_data_ratio[ibin]->GetXaxis()->SetLabelSize(0.07);
    h_eta_MC_data_ratio[ibin]->GetXaxis()->SetRangeUser(-1.5,1.52);
    h_eta_MC_data_ratio[ibin]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_eta_MC_data_ratio[ibin]->SetLineColor(kBlack); h_eta_MC_data_ratio[ibin]->SetMarkerColor(kBlack);
    h_eta_MC_data_ratio[ibin]->SetMarkerStyle(20); h_eta_MC_data_ratio[ibin]->Draw("e0 same");
    tl1->Draw("same");
  }

  TCanvas *c_qg_eta = new TCanvas("c_qg_eta","",1500,350);
  c_qg_eta->Divide(5,1,0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0) c_qg_eta->cd(ibin+1);
    else c_qg_eta->cd(6-ibin);   
    h_eta_MC_g[ibin]->GetXaxis()->SetTitle("jet #eta");
    h_eta_MC_g[ibin]->GetYaxis()->SetNdivisions(505);
    h_eta_MC_g[ibin]->GetYaxis()->SetTitleSize(0.05);
    h_eta_MC_g[ibin]->GetYaxis()->SetLabelSize(0.07);
    h_eta_MC_g[ibin]->GetXaxis()->CenterTitle();
    h_eta_MC_g[ibin]->GetXaxis()->SetNdivisions(505);
    h_eta_MC_g[ibin]->GetXaxis()->SetTitleSize(0.05);
    h_eta_MC_g[ibin]->GetXaxis()->SetLabelSize(0.07);
    h_eta_MC_g[ibin]->GetXaxis()->SetRangeUser(-1.5,1.52);
    h_eta_MC_g[ibin]->GetYaxis()->SetRangeUser(0.02,0.12);
    h_eta_MC_g[ibin]->SetLineColor(kRed); h_eta_MC_g[ibin]->SetMarkerColor(kRed);
    h_eta_MC_g[ibin]->SetMarkerStyle(4); h_eta_MC_g[ibin]->Draw("e0 same");
    h_eta_MC_q[ibin]->SetLineColor(kBlue); h_eta_MC_q[ibin]->SetMarkerColor(kBlue);
    h_eta_MC_q[ibin]->SetMarkerStyle(4); h_eta_MC_q[ibin]->Draw("e0 same");
    //h_eta_data[ibin]->Draw("e0 same"); h_eta_MC[ibin]->Draw("e0 same");
    //h_eta_data[ibin]->Fit(f_eta_reco[ibin],"N M R","",-1.5,1.5);
    //cout<<f_eta_reco[ibin]->GetChisquare()/f_eta_reco[ibin]->GetNDF()<<endl;
    h_eta_MC_up[ibin]->SetLineColor(kBlue); h_eta_MC_up[ibin]->SetMarkerColor(kBlue);
    h_eta_MC_up[ibin]->SetMarkerStyle(22); h_eta_MC_up[ibin]->Draw("e0 same");
    h_eta_MC_upbar[ibin]->SetLineColor(kBlue); h_eta_MC_upbar[ibin]->SetMarkerColor(kBlue);
    h_eta_MC_upbar[ibin]->SetMarkerStyle(26); h_eta_MC_upbar[ibin]->Draw("e0 same");
    h_eta_MC_d[ibin]->SetLineColor(kBlue); h_eta_MC_d[ibin]->SetMarkerColor(kBlue);
    h_eta_MC_d[ibin]->SetMarkerStyle(23); h_eta_MC_d[ibin]->Draw("e0 same");
    h_eta_MC_dbar[ibin]->SetLineColor(kBlue); h_eta_MC_dbar[ibin]->SetMarkerColor(kBlue);
    h_eta_MC_dbar[ibin]->SetMarkerStyle(32); h_eta_MC_dbar[ibin]->Draw("e0 same");

    if(ibin>0){
      TString tmp1 = cent_tag[ibin-1];
      tx1 = new TLatex(); tx1->SetTextSize(.12);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }
  }

  TCanvas *c_trk_bkg_MC = new TCanvas("c_trk_bkg_MC","",1500,350);
  c_trk_bkg_MC->Divide(5,1,0);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(0);

  for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0) c_trk_bkg_MC->cd(ibin+1);
    else c_trk_bkg_MC->cd(6-ibin);   
    h_bkgtrk_MC[ibin]->GetXaxis()->SetTitle("# tracks");
    h_bkgtrk_MC[ibin]->GetYaxis()->SetNdivisions(505);
    h_bkgtrk_MC[ibin]->GetYaxis()->SetTitleSize(0.05);
    h_bkgtrk_MC[ibin]->GetYaxis()->SetLabelSize(0.07);
    h_bkgtrk_MC[ibin]->GetXaxis()->CenterTitle();
    h_bkgtrk_MC[ibin]->GetXaxis()->SetNdivisions(505);
    h_bkgtrk_MC[ibin]->GetXaxis()->SetTitleSize(0.05);
    h_bkgtrk_MC[ibin]->GetXaxis()->SetLabelSize(0.07);
    h_bkgtrk_MC[ibin]->GetXaxis()->SetRangeUser(0.,40.);
    //h_bkgtrk_MC[ibin]->GetYaxis()->SetRangeUser(0.,0.05e-03);
    h_bkgtrk_MC[ibin]->SetLineColor(kRed); h_bkgtrk_MC[ibin]->SetMarkerColor(kRed);
    h_bkgtrk_MC[ibin]->SetMarkerStyle(20); h_bkgtrk_MC[ibin]->Draw("e0 same");
    h_trk_MC[ibin]->SetLineColor(kBlue); h_trk_MC[ibin]->SetMarkerColor(kBlue);
    h_trk_MC[ibin]->SetMarkerStyle(20); h_trk_MC[ibin]->Draw("e0 same");

    if(ibin>0){
      TString tmp1 = cent_tag[ibin-1];
      tx1 = new TLatex(); tx1->SetTextSize(.12);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }
  }

  TCanvas *c_trk_bkg_data = new TCanvas("c_trk_bkg_data","",1500,350);
  c_trk_bkg_data->Divide(5,1,0);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(0);

  for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0) c_trk_bkg_data->cd(ibin+1);
    else c_trk_bkg_data->cd(6-ibin);   
    h_bkgtrk_data[ibin]->GetXaxis()->SetTitle("# tracks");
    h_bkgtrk_data[ibin]->GetYaxis()->SetNdivisions(505);
    h_bkgtrk_data[ibin]->GetYaxis()->SetTitleSize(0.05);
    h_bkgtrk_data[ibin]->GetYaxis()->SetLabelSize(0.07);
    h_bkgtrk_data[ibin]->GetXaxis()->CenterTitle();
    h_bkgtrk_data[ibin]->GetXaxis()->SetNdivisions(505);
    h_bkgtrk_data[ibin]->GetXaxis()->SetTitleSize(0.05);
    h_bkgtrk_data[ibin]->GetXaxis()->SetLabelSize(0.07);
    h_bkgtrk_data[ibin]->GetXaxis()->SetRangeUser(0.,40.);
    //h_bkgtrk_data[ibin]->GetYaxis()->SetRangeUser(0.,0.05e-03);
    h_bkgtrk_data[ibin]->SetLineColor(kRed); h_bkgtrk_data[ibin]->SetMarkerColor(kRed);
    h_bkgtrk_data[ibin]->SetMarkerStyle(20); h_bkgtrk_data[ibin]->Draw("e0 same");
    h_trk_data[ibin]->SetLineColor(kBlue); h_trk_data[ibin]->SetMarkerColor(kBlue);
    h_trk_data[ibin]->SetMarkerStyle(20); h_trk_data[ibin]->Draw("e0 same");

    if(ibin>0){
      TString tmp1 = cent_tag[ibin-1];
      tx1 = new TLatex(); tx1->SetTextSize(.12);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }
  }

  TCanvas *c_reco = new TCanvas("c_reco","",1500,400);
  c_reco->Divide(5,1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0) c_reco->cd(ibin+1);
    else c_reco->cd(6-ibin);   
    gPad->SetLogy();
    h_reco_MC[ibin]->GetXaxis()->SetTitle("jet pT");
    h_reco_MC[ibin]->GetYaxis()->SetNdivisions(505);
    h_reco_MC[ibin]->GetYaxis()->SetTitleSize(0.05);
    h_reco_MC[ibin]->GetYaxis()->SetLabelSize(0.07);
    h_reco_MC[ibin]->GetXaxis()->CenterTitle();
    h_reco_MC[ibin]->GetXaxis()->SetNdivisions(505);
    h_reco_MC[ibin]->GetXaxis()->SetTitleSize(0.05);
    h_reco_MC[ibin]->GetXaxis()->SetLabelSize(0.07);
    h_reco_MC[ibin]->GetXaxis()->SetRangeUser(120.,400.);
    //h_reco_MC[ibin]->GetYaxis()->SetRangeUser(0.05,0.2);
    h_reco_MC[ibin]->SetLineColor(kBlack); h_reco_MC[ibin]->SetMarkerColor(kBlack);
    h_reco_MC[ibin]->SetMarkerStyle(4); h_reco_MC[ibin]->Draw("e0 same");
    //h_reco_data[ibin]->SetLineColor(kGreen); h_reco_data[ibin]->SetMarkerColor(kGreen);
    //h_reco_data[ibin]->SetMarkerStyle(4); h_reco_data[ibin]->Draw("e0 same");
    h_reco_MC_q[ibin]->SetLineColor(kBlue); h_reco_MC_q[ibin]->SetMarkerColor(kBlue);
    h_reco_MC_q[ibin]->SetMarkerStyle(4); h_reco_MC_q[ibin]->Draw("e0 same");
    h_reco_MC_g[ibin]->SetLineColor(kRed); h_reco_MC_g[ibin]->SetMarkerColor(kRed);
    h_reco_MC_g[ibin]->SetMarkerStyle(4); h_reco_MC_g[ibin]->Draw("e0 same");

    h_reco_MC_up[ibin]->SetLineColor(kBlue); h_reco_MC_up[ibin]->SetMarkerColor(kBlue);
    h_reco_MC_up[ibin]->SetMarkerStyle(22); h_reco_MC_up[ibin]->Draw("e0 same");
    h_reco_MC_upbar[ibin]->SetLineColor(kBlue); h_reco_MC_upbar[ibin]->SetMarkerColor(kBlue);
    h_reco_MC_upbar[ibin]->SetMarkerStyle(26); h_reco_MC_upbar[ibin]->Draw("e0 same");

    h_reco_MC_d[ibin]->SetLineColor(kBlue); h_reco_MC_d[ibin]->SetMarkerColor(kBlue);
    h_reco_MC_d[ibin]->SetMarkerStyle(23); h_reco_MC_d[ibin]->Draw("e0 same");
    h_reco_MC_dbar[ibin]->SetLineColor(kBlue); h_reco_MC_dbar[ibin]->SetMarkerColor(kBlue);
    h_reco_MC_dbar[ibin]->SetMarkerStyle(32); h_reco_MC_dbar[ibin]->Draw("e0 same");

    if(ibin>0){
      TString tmp1 = cent_tag[ibin-1];
      tx1 = new TLatex(); tx1->SetTextSize(.12);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }
/*
    if(ibin==0) c_reco->cd(ibin+6);
    else c_reco->cd(11-ibin);   
    h_reco_MC_data_ratio[ibin] = (TH1D*)h_reco_MC[ibin]->Clone("h_reco_MC_data_ratio");
    h_reco_MC_data_ratio[ibin]->Divide(h_reco_data[ibin]);
    h_reco_MC_data_ratio[ibin]->GetYaxis()->SetNdivisions(505);
    h_reco_MC_data_ratio[ibin]->GetYaxis()->SetTitleSize(0.05);
    h_reco_MC_data_ratio[ibin]->GetYaxis()->SetLabelSize(0.07);
    h_reco_MC_data_ratio[ibin]->GetXaxis()->CenterTitle();
    h_reco_MC_data_ratio[ibin]->GetXaxis()->SetNdivisions(505);
    h_reco_MC_data_ratio[ibin]->GetXaxis()->SetTitleSize(0.05);
    h_reco_MC_data_ratio[ibin]->GetXaxis()->SetLabelSize(0.07);
    h_reco_MC_data_ratio[ibin]->GetXaxis()->SetRangeUser(120.,400.);
    h_reco_MC_data_ratio[ibin]->GetYaxis()->SetRangeUser(0.5,1.5);
    h_reco_MC_data_ratio[ibin]->SetLineColor(kBlack); h_reco_MC_data_ratio[ibin]->SetMarkerColor(kBlack);
    h_reco_MC_data_ratio[ibin]->SetMarkerStyle(20); h_reco_MC_data_ratio[ibin]->Draw("e0 same");
    tl2->Draw("same");

  }
*/

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

  TCanvas *c_chg_reco_bkg = new TCanvas("c_chg_reco_bkg","",1500,1200);
  c_chg_reco_bkg->Divide(5,4,0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0) c_chg_reco_bkg->cd((5*ibin2)+1);
    else c_chg_reco_bkg->cd((5*(ibin2+1))-(ibin-1));   
    h_chg_reco_q[ibin][ibin2]->GetXaxis()->SetTitle("jet chg");
    h_chg_reco_q[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_reco_q[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_reco_q[ibin][ibin2]->GetYaxis()->SetLabelSize(0.07);
    h_chg_reco_q[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_reco_q[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_reco_q[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_reco_q[ibin][ibin2]->GetXaxis()->SetLabelSize(0.07);
    h_chg_reco_q[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    //h_chg_reco_q[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);
    h_chg_reco_q[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_q[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_q[ibin][ibin2]->SetMarkerStyle(2); h_chg_reco_q[ibin][ibin2]->Draw("e0 same");
    h_chg_reco_q_refl[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_q_refl[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_q_refl[ibin][ibin2]->SetMarkerStyle(23); h_chg_reco_q_refl[ibin][ibin2]->Draw("e0 same");
    //h_chg_reco[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco[ibin][ibin2]->SetMarkerColor(kBlack);
    //h_chg_reco[ibin][ibin2]->SetMarkerStyle(20); h_chg_reco[ibin][ibin2]->Draw("e0 same");
    //h_chg_reco_up[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_up[ibin][ibin2]->SetMarkerColor(kBlue);
    //h_chg_reco_up[ibin][ibin2]->SetMarkerStyle(22); h_chg_reco_up[ibin][ibin2]->Draw("e0 same");
    //h_chg_reco_down[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_down[ibin][ibin2]->SetMarkerColor(kBlue);
    //h_chg_reco_down[ibin][ibin2]->SetMarkerStyle(23); h_chg_reco_down[ibin][ibin2]->Draw("e0 same");
    //h_chg_reco_upb[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_upb[ibin][ibin2]->SetMarkerColor(kBlue);
    //h_chg_reco_upb[ibin][ibin2]->SetMarkerStyle(26); h_chg_reco_upb[ibin][ibin2]->Draw("e0 same");
    //h_chg_reco_downb[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_downb[ibin][ibin2]->SetMarkerColor(kBlue);
    //h_chg_reco_downb[ibin][ibin2]->SetMarkerStyle(32); h_chg_reco_downb[ibin][ibin2]->Draw("e0 same");
    h_chg_reco_g_refl[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_g_refl[ibin][ibin2]->SetMarkerColor(kRed);
    h_chg_reco_g_refl[ibin][ibin2]->SetMarkerStyle(2); h_chg_reco_g_refl[ibin][ibin2]->Draw("e0 same");
    
    if(ibin>0){
      TString tmp1 = cent_tag[ibin-1];
      tx1 = new TLatex(); tx1->SetTextSize(.12);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }
    
    TString tmp = trk_tag[ibin2];
    tx = new TLatex(); tx->SetTextSize(.12);
    tx->DrawLatexNDC(0.15,0.65, tmp);
    }
  }

  TCanvas *c_chg_reco_refl = new TCanvas("c_chg_reco_refl","",1500,1200);
  c_chg_reco_refl->Divide(5,4);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0) c_chg_reco_refl->cd((5*ibin2)+1);
    else c_chg_reco_refl->cd((5*(ibin2+1))-(ibin-1));   
    h_chg_reco_q_refl[ibin][ibin2]->GetXaxis()->SetTitle("jet chg");
    h_chg_reco_q_refl[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_reco_q_refl[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_reco_q_refl[ibin][ibin2]->GetYaxis()->SetLabelSize(0.07);
    h_chg_reco_q_refl[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_reco_q_refl[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_reco_q_refl[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_reco_q_refl[ibin][ibin2]->GetXaxis()->SetLabelSize(0.07);
    h_chg_reco_q_refl[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_chg_reco_q_refl[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_q_refl[ibin][ibin2]->SetMarkerColor(kBlue);
    /*h_chg_reco_q_refl[ibin][ibin2]->Rebin(10);*/ h_chg_reco_q_refl[ibin][ibin2]->SetMarkerStyle(2); h_chg_reco_q_refl[ibin][ibin2]->Draw("e0 same");
    h_chg_reco_g_refl[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_g_refl[ibin][ibin2]->SetMarkerColor(kRed);
    /*h_chg_reco_g_refl[ibin][ibin2]->Rebin(10);*/ h_chg_reco_g_refl[ibin][ibin2]->SetMarkerStyle(2); h_chg_reco_g_refl[ibin][ibin2]->Draw("e0 same");
    h_chg_reco_MC_refl[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_MC_refl[ibin][ibin2]->SetMarkerColor(kBlack);
    /*h_chg_reco_MC_refl[ibin][ibin2]->Rebin(10);*/ h_chg_reco_MC_refl[ibin][ibin2]->SetMarkerStyle(2); h_chg_reco_MC_refl[ibin][ibin2]->Draw("e0 same");
    //h_chg_reco_data_refl[ibin][ibin2]->SetLineColor(kGreen); h_chg_reco_data_refl[ibin][ibin2]->SetMarkerColor(kGreen);
    //h_chg_reco_data_refl[ibin][ibin2]->Rebin(10); h_chg_reco_data_refl[ibin][ibin2]->SetMarkerStyle(2); h_chg_reco_data_refl[ibin][ibin2]->Draw("e0 same");
    cout<<ibin<<"  "<<ibin2<<" q fraction: "<<h_chg_reco_q[ibin][ibin2]->Integral()/h_chg_reco[ibin][ibin2]->Integral()<<endl;
    cout<<ibin<<"  "<<ibin2<<" q fraction refl: "<<h_chg_reco_MC_refl[ibin][ibin2]->Integral()/h_chg_reco[ibin][ibin2]->Integral()<<endl;

    if(ibin>0){
      TString tmp1 = cent_tag[ibin-1];
      tx1 = new TLatex(); tx1->SetTextSize(.12);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }
    
    TString tmp = trk_tag[ibin2];
    tx = new TLatex(); tx->SetTextSize(.09);
    tx->DrawLatexNDC(0.15,0.65, tmp);
    }
  }

  TCanvas *c_chg_data_bkg = new TCanvas("c_chg_data_bkg","",1500,1200);
  c_chg_data_bkg->Divide(5,4,0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0) c_chg_data_bkg->cd((5*ibin2)+1);
    else c_chg_data_bkg->cd((5*(ibin2+1))-(ibin-1));  
    h_chg_data_bkg[ibin][ibin2]->GetXaxis()->SetTitle("jet chg");
    h_chg_data_bkg[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_data_bkg[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_data_bkg[ibin][ibin2]->GetYaxis()->SetLabelSize(0.07);
    h_chg_data_bkg[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_data_bkg[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_data_bkg[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_data_bkg[ibin][ibin2]->GetXaxis()->SetLabelSize(0.07);
    h_chg_data_bkg[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
    h_chg_data_bkg[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.25);
    h_chg_data_bkg[ibin][ibin2]->SetLineColor(kRed); h_chg_data_bkg[ibin][ibin2]->SetMarkerColor(kRed);
    h_chg_data_bkg[ibin][ibin2]->SetMarkerStyle(20); h_chg_data_bkg[ibin][ibin2]->Draw("e0 same");
    h_chg_data[ibin][ibin2]->SetLineColor(kGreen); h_chg_data[ibin][ibin2]->SetMarkerColor(kGreen);
    h_chg_data[ibin][ibin2]->SetMarkerStyle(2); h_chg_data[ibin][ibin2]->Draw("e0 same");

    if(ibin>0){
      TString tmp1 = cent_tag[ibin-1];
      tx1 = new TLatex(); tx1->SetTextSize(.12);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }

    TString tmp = trk_tag[ibin2];
    tx = new TLatex(); tx->SetTextSize(.12);
    tx->DrawLatexNDC(0.15,0.65, tmp);
    }
  }

  TCanvas *c_chg_qg_avg = new TCanvas("c_chg_qg_avg","",1200,1200);
  c_chg_qg_avg->Divide(4,4,0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    h_chg_reco_avg[0][ibin2]->Rebin(20);
    h_chg_reco_avg[0][ibin2]->SetLineColor(kBlack); h_chg_reco_avg[0][ibin2]->SetMarkerColor(kBlack);
    h_chg_reco_avg[0][ibin2]->SetMarkerStyle(4);

    for(int ibin=1;ibin<nCBins;ibin++){
    c_chg_qg_avg->cd((4*(ibin2+1))-(ibin-1));   
    h_chg_reco_avg[ibin][ibin2]->Rebin(20);
    h_chg_reco_bkg_avg[ibin][ibin2]->Rebin(20);
    h_chg_reco_q_avg[ibin][ibin2]->Rebin(20);
    h_chg_reco_g_avg[ibin][ibin2]->Rebin(20);
    h_chg_reco_avg[ibin][ibin2]->GetXaxis()->SetTitle("jet pT");
    //h_chg_reco_avg[ibin][ibin2]->GetYaxis()->SetTitle("<jet chg>");
    h_chg_reco_avg[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_reco_avg[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_reco_avg[ibin][ibin2]->GetYaxis()->SetLabelSize(0.07);
    h_chg_reco_avg[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_reco_avg[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_reco_avg[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_reco_avg[ibin][ibin2]->GetXaxis()->SetLabelSize(0.07);
    h_chg_reco_avg[ibin][ibin2]->GetXaxis()->SetRangeUser(120.,300.);
    h_chg_reco_avg[ibin][ibin2]->GetYaxis()->SetRangeUser(-0.5,0.5);
    h_chg_reco_avg[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_avg[ibin][ibin2]->SetMarkerColor(kBlack);
    h_chg_reco_avg[ibin][ibin2]->SetMarkerStyle(7); h_chg_reco_avg[ibin][ibin2]->Draw("e0 same");

    h_chg_reco_q_avg[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_q_avg[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_q_avg[ibin][ibin2]->SetMarkerStyle(7); h_chg_reco_q_avg[ibin][ibin2]->Draw("e0 same");
    h_chg_reco_g_avg[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_g_avg[ibin][ibin2]->SetMarkerColor(kRed);
    h_chg_reco_g_avg[ibin][ibin2]->SetMarkerStyle(7); h_chg_reco_g_avg[ibin][ibin2]->Draw("e0 same");

    h_chg_reco_up_avg[ibin][ibin2]->Rebin(20); h_chg_reco_down_avg[ibin][ibin2]->Rebin(20);
    h_chg_reco_up_avg[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_up_avg[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_up_avg[ibin][ibin2]->SetMarkerStyle(22); h_chg_reco_up_avg[ibin][ibin2]->Draw("e0 same");
    h_chg_reco_down_avg[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_down_avg[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_down_avg[ibin][ibin2]->SetMarkerStyle(23); h_chg_reco_down_avg[ibin][ibin2]->Draw("e0 same");    

    h_chg_reco_upb_avg[ibin][ibin2]->Rebin(20); h_chg_reco_downb_avg[ibin][ibin2]->Rebin(20);
    h_chg_reco_upb_avg[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_upb_avg[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_upb_avg[ibin][ibin2]->SetMarkerStyle(26); h_chg_reco_upb_avg[ibin][ibin2]->Draw("e0 same");
    h_chg_reco_downb_avg[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_downb_avg[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_reco_downb_avg[ibin][ibin2]->SetMarkerStyle(32); h_chg_reco_downb_avg[ibin][ibin2]->Draw("e0 same");

    //h_chg_reco_avg[0][ibin2]->Draw("e0 same");
/*
    h_chg_ref_q_avg[ibin]->SetLineColor(kBlue); h_chg_ref_q_avg[ibin]->SetMarkerColor(kBlue);
    h_chg_ref_q_avg[ibin]->SetMarkerStyle(7); h_chg_ref_q_avg[ibin]->Draw("e0 same");
    h_chg_ref_g_avg[ibin]->SetLineColor(kRed); h_chg_ref_g_avg[ibin]->SetMarkerColor(kRed);
    h_chg_ref_g_avg[ibin]->SetMarkerStyle(7); h_chg_ref_g_avg[ibin]->Draw("e0 same");
*/
    //h_chg_reco_bkg_avg[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_bkg_avg[ibin][ibin2]->SetMarkerColor(kRed);
    //h_chg_reco_bkg_avg[ibin][ibin2]->SetMarkerStyle(7); h_chg_reco_bkg_avg[ibin][ibin2]->Draw("e0 same");

    TString tmp1 = cent_tag[ibin-1];
    tx1 = new TLatex(); tx1->SetTextSize(.12);
    tx1->DrawLatexNDC(0.15,0.85, tmp1);

    TString tmp = trk_tag[ibin2];
    tx = new TLatex(); tx->SetTextSize(.12);
    tx->DrawLatexNDC(0.15,0.65, tmp);

    tl2->Draw("same"); 
/*
    h_chg_ref_q_avg[ibin]->SetLineColor(kBlue); h_chg_ref_q_avg[ibin]->SetMarkerColor(kBlue);
    h_chg_ref_q_avg[ibin]->SetMarkerStyle(4); h_chg_ref_q_avg[ibin]->Draw("e0 same");
    h_chg_ref_g_avg[ibin]->SetLineColor(kRed); h_chg_ref_g_avg[ibin]->SetMarkerColor(kRed);
    h_chg_ref_g_avg[ibin]->SetMarkerStyle(4); h_chg_ref_g_avg[ibin]->Draw("e0 same");
*/
    }
  }

  TCanvas *c_chg_data_avg = new TCanvas("c_chg_data_avg","",1200,1200);
  c_chg_data_avg->Divide(4,4,0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    h_chg_data_avg[0][ibin2]->Rebin(20);
    h_chg_data_avg[0][ibin2]->SetLineColor(kBlack); h_chg_data_avg[0][ibin2]->SetMarkerColor(kBlack);
    h_chg_data_avg[0][ibin2]->SetMarkerStyle(4);

    for(int ibin=1;ibin<nCBins;ibin++){
    c_chg_data_avg->cd((4*(ibin2+1))-(ibin-1));   
    h_chg_data_avg[ibin][ibin2]->Rebin(20);
    h_chg_data_bkg_avg[ibin][ibin2]->Rebin(20);
    //if(ibin2==0)h_chg_data_avg[ibin][ibin2]->SetTitle(cent_tag[ibin-1]);
    h_chg_data_avg[ibin][ibin2]->GetXaxis()->SetTitle("jet pT");
    //h_chg_data_avg[ibin][ibin2]->GetYaxis()->SetTitle("<jet chg>");
    h_chg_data_avg[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_data_avg[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_data_avg[ibin][ibin2]->GetYaxis()->SetLabelSize(0.07);
    h_chg_data_avg[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_data_avg[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_data_avg[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_data_avg[ibin][ibin2]->GetXaxis()->SetLabelSize(0.07);
    h_chg_data_avg[ibin][ibin2]->GetXaxis()->SetRangeUser(120.,300.);
    h_chg_data_avg[ibin][ibin2]->GetYaxis()->SetRangeUser(-0.15,0.4);
    h_chg_data_avg[ibin][ibin2]->SetLineColor(kBlack); h_chg_data_avg[ibin][ibin2]->SetMarkerColor(kBlack);
    h_chg_data_avg[ibin][ibin2]->SetMarkerStyle(7); h_chg_data_avg[ibin][ibin2]->Draw("e0 same");
    h_chg_data_avg[0][ibin2]->Draw("e0 same");

    h_chg_data_bkg_avg[ibin][ibin2]->SetLineColor(kRed); h_chg_data_bkg_avg[ibin][ibin2]->SetMarkerColor(kRed);
    h_chg_data_bkg_avg[ibin][ibin2]->SetMarkerStyle(7); h_chg_data_bkg_avg[ibin][ibin2]->Draw("e0 same");

    TString tmp1 = cent_tag[ibin-1];
    tx1 = new TLatex(); tx1->SetTextSize(.12);
    tx1->DrawLatexNDC(0.15,0.85, tmp1);

    TString tmp = trk_tag[ibin2];
    tx = new TLatex(); tx->SetTextSize(.12);
    tx->DrawLatexNDC(0.15,0.65, tmp);
    }
  }

///////////////////////////////////////////////
///////////// fiiting trkpt>0.7 plots //////////
///////////////////////////////////////////////

  TString pythia = "Pythia";
  TString ppdata = "pp data 5.02 TeV";
  TString cms = "CMS";
  TString fit_ratio = "Inclusive / Fit";
  TString ncs_label = "1 / N_{jets} * N";

  char q_frac[50],q_error[50],chi2ndf[50];
  TString tmp;

  TLatex *title_tex;
  title_tex = new TLatex();
  title_tex->SetTextSize(0.55);
  title_tex->SetLineColor(kWhite);

  TLatex *ratio_tex;
  ratio_tex = new TLatex();
  ratio_tex->SetTextAngle(90);
  ratio_tex->SetTextSize(0.55);
  ratio_tex->SetLineColor(kWhite);

  TLatex *ncs_tex;
  ncs_tex = new TLatex();
  ncs_tex->SetTextAngle(90);
  ncs_tex->SetTextSize(0.65);
  ncs_tex->SetLineColor(kWhite);

  TCanvas *c_chg_MC_data_fin = new TCanvas("c_chg_MC_data_fin","",900,450);
  TPad *pad_title = new TPad("pad_title", "",0.,0.92,0.99,0.99);
  TPad *pad_label = new TPad("pad_label", "",0.,0.,0.03,0.95);
  TPad *pad1 = new TPad("pad1", "",0.02,0.35,0.99,0.98);
  TPad *pad2 = new TPad("pad2", "",0.02,0.,0.99,0.35);
  pad1->Divide(3,1);  
  pad1->Draw();
  pad2->Divide(3,1);  
  pad2->Draw();
  pad_title->Draw();
  pad_label->Draw();
  pad_title->cd();
  title_tex->DrawLatexNDC(0.3,0.5, pythia);
  title_tex->DrawLatexNDC(0.78,0.5, ppdata);
  title_tex->DrawLatexNDC(0.05,0.5, cms);  
  pad_label->cd();
  ratio_tex->DrawLatexNDC(0.7,0.1, fit_ratio);
  ncs_tex->DrawLatexNDC(0.7,0.6, ncs_label);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  for(int ibin=0;ibin<1;ibin++){
    pad1->cd(1);   
    h_chg_ref[ibin][0]->GetXaxis()->SetTitle("jet charge");
    h_chg_ref[ibin][0]->GetYaxis()->SetTitle("");
    h_chg_ref[ibin][0]->GetYaxis()->SetRangeUser(0.,1.1*h_chg_ref_g[ibin][0]->GetMaximum());
    h_chg_ref[ibin][0]->Draw("e0 same"); h_chg_ref[ibin][0]->Draw("l hist same");
    h_chg_ref_q[ibin][0]->Draw("e0 same"); h_chg_ref_q[ibin][0]->Draw("l hist same");
    h_chg_ref_g[ibin][0]->Draw("e0 same"); h_chg_ref_g[ibin][0]->Draw("l hist same");

    h_chg_ref[ibin][0]->Fit(f_chg_ref[ibin][0],"N Q M R","",-2.,2.);

    sprintf(q_frac, "quark %.3f",f_chg_ref[ibin][0]->GetParameter(0));
    sprintf(q_error, "#pm %.3f",f_chg_ref[ibin][0]->GetParError(0));
    sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_ref[ibin][0]->GetChisquare()/f_chg_ref[ibin][0]->GetNDF());
    tx1 = new TLatex(); tx1->SetTextSize(.07);
    tx1->DrawLatexNDC(0.35,0.25, q_frac);
    tx3 = new TLatex(); tx3->SetTextSize(.06);
    tx3->DrawLatexNDC(0.5,0.2, q_error);

    pad2->cd(1);  
    h_chg_ref_ratio[ibin][0]->GetXaxis()->SetTitle("");
    h_chg_ref_ratio[ibin][0]->GetYaxis()->SetTitle("");
    h_chg_ref_ratio[ibin][0]->SetMarkerStyle(7); h_chg_ref_ratio[ibin][0]->Draw("e0 same");
    tl1->Draw("same");
    tx2 = new TLatex(); tx2->SetTextSize(.10);
    tx2->DrawLatexNDC(0.3,0.2, chi2ndf);

    pad1->cd(2);   
    h_chg_reco[ibin][0]->GetXaxis()->SetTitle("jet charge");
    h_chg_reco[ibin][0]->GetYaxis()->SetTitle("");
    h_chg_reco[ibin][0]->GetYaxis()->SetRangeUser(0.,1.1*h_chg_reco_g[ibin][0]->GetMaximum());
    h_chg_reco[ibin][0]->Draw("e0 same"); h_chg_reco[ibin][0]->Draw("l hist same");
    h_chg_reco_q[ibin][0]->Draw("e0 same"); h_chg_reco_q[ibin][0]->Draw("l hist same");
    h_chg_reco_g[ibin][0]->Draw("e0 same"); h_chg_reco_g[ibin][0]->Draw("l hist same");

    h_chg_reco[ibin][0]->Fit(f_chg_reco[ibin][0],"N Q M R","",-2.,2.);

    sprintf(q_frac, "quark %.3f",f_chg_reco[ibin][0]->GetParameter(0));
    sprintf(q_error, "#pm %.3f",f_chg_reco[ibin][0]->GetParError(0));
    sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_reco[ibin][0]->GetChisquare()/f_chg_reco[ibin][0]->GetNDF());
    tx1 = new TLatex(); tx1->SetTextSize(.07);
    tx1->DrawLatexNDC(0.35,0.25, q_frac);
    tx3 = new TLatex(); tx3->SetTextSize(.06);
    tx3->DrawLatexNDC(0.5,0.2, q_error);

    pad2->cd(2);  
    h_chg_reco_ratio[ibin][0]->GetXaxis()->SetTitle("");
    h_chg_reco_ratio[ibin][0]->GetYaxis()->SetTitle("");
    h_chg_reco_ratio[ibin][0]->SetMarkerStyle(7); h_chg_reco_ratio[ibin][0]->Draw("e0 same");
    tl1->Draw("same");
    tx2 = new TLatex(); tx2->SetTextSize(.10);
    tx2->DrawLatexNDC(0.3,0.2, chi2ndf);

    pad1->cd(3);   
    h_chg_data[ibin][0]->GetXaxis()->SetTitle("jet charge");
    h_chg_data[ibin][0]->GetYaxis()->SetTitle("");
    h_chg_data[ibin][0]->GetYaxis()->SetRangeUser(0.,1.1*h_chg_reco_g[ibin][0]->GetMaximum());
    h_chg_data[ibin][0]->Draw("e0 same"); h_chg_data[ibin][0]->Draw("l hist same");
    h_chg_reco_q[ibin][0]->Draw("e0 same"); h_chg_reco_q[ibin][0]->Draw("l hist same");
    h_chg_reco_g[ibin][0]->Draw("e0 same"); h_chg_reco_g[ibin][0]->Draw("l hist same");

    h_chg_data[ibin][0]->Fit(f_chg_reco[ibin][0],"N Q M R","",-2.,2.);

    sprintf(q_frac, "quark %.3f",f_chg_reco[ibin][0]->GetParameter(0));
    sprintf(q_error, "#pm %.3f",f_chg_reco[ibin][0]->GetParError(0));
    sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_reco[ibin][0]->GetChisquare()/f_chg_reco[ibin][0]->GetNDF());
    tx1 = new TLatex(); tx1->SetTextSize(.07);
    tx1->DrawLatexNDC(0.35,0.25, q_frac);
    tx3 = new TLatex(); tx3->SetTextSize(.06);
    tx3->DrawLatexNDC(0.5,0.2, q_error);

    pad2->cd(3);  
    h_chg_data_ratio[ibin][0]->GetXaxis()->SetTitle("");
    h_chg_data_ratio[ibin][0]->GetYaxis()->SetTitle("");
    h_chg_data_ratio[ibin][0]->SetMarkerStyle(7); h_chg_data_ratio[ibin][0]->Draw("e0 same");
    tl1->Draw("same");
    tx2 = new TLatex(); tx2->SetTextSize(.10);
    tx2->DrawLatexNDC(0.3,0.2, chi2ndf);
  }  

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_chg_MC_fin = new TCanvas("c_chg_MC_fin","",900,450);
  TPad *pad_title = new TPad("pad_title", "",0.,0.92,0.99,0.99);
  TPad *pad_label = new TPad("pad_label", "",0.,0.,0.03,0.95);
  TPad *pad1 = new TPad("pad1", "",0.02,0.35,0.99,0.98);
  TPad *pad2 = new TPad("pad2", "",0.02,0.,0.99,0.35);
  pad1->Divide(3,1);  
  pad1->Draw();
  pad2->Divide(3,1);  
  pad2->Draw();
  pad_title->Draw();
  pad_label->Draw();
  pad_title->cd();
  title_tex->DrawLatexNDC(0.5,0.5, pythia);
  title_tex->DrawLatexNDC(0.05,0.5, cms);  
  pad_label->cd();
  ratio_tex->DrawLatexNDC(0.7,0.1, fit_ratio);
  ncs_tex->DrawLatexNDC(0.7,0.6, ncs_label);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  for(int ibin=1;ibin<4;ibin++){
    pad1->cd(ibin);   
    h_chg_reco[0][ibin]->GetXaxis()->SetTitle("jet charge");
    h_chg_reco[0][ibin]->GetYaxis()->SetTitle("");
    h_chg_reco[0][ibin]->GetYaxis()->SetRangeUser(0.,1.1*h_chg_reco_q[0][ibin]->GetMaximum());
    h_chg_reco[0][ibin]->Draw("e0 same"); h_chg_reco[0][ibin]->Draw("l hist same");
    h_chg_reco_q[0][ibin]->Draw("e0 same"); h_chg_reco_q[0][ibin]->Draw("l hist same");
    h_chg_reco_g[0][ibin]->Draw("e0 same"); h_chg_reco_g[0][ibin]->Draw("l hist same");

    h_chg_reco[0][ibin]->Fit(f_chg_data[0][ibin],"N Q M R","",-2.,2.);

    sprintf(q_frac, "quark %.3f",f_chg_data[0][ibin]->GetParameter(0));
    sprintf(q_error, "#pm %.3f",f_chg_data[0][ibin]->GetParError(0));
    sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_data[0][ibin]->GetChisquare()/f_chg_data[0][ibin]->GetNDF());
    tmp = trk_tag[ibin-1];
    tx = new TLatex(); tx->SetTextSize(.07);
    tx->DrawLatexNDC(0.35,0.25, tmp);

    pad2->cd(ibin);  
    h_chg_reco_ratio[0][ibin]->GetXaxis()->SetTitle("");
    h_chg_reco_ratio[0][ibin]->GetYaxis()->SetTitle("");
    h_chg_reco_ratio[0][ibin]->SetMarkerStyle(7); h_chg_reco_ratio[0][ibin]->Draw("e0 same");
    tl1->Draw("same");
    tx1 = new TLatex(); tx1->SetTextSize(.11);
    tx1->DrawLatexNDC(0.35,0.75, q_frac);
    tx3 = new TLatex(); tx3->SetTextSize(.10);
    tx3->DrawLatexNDC(0.47,0.67, q_error);
    tx2 = new TLatex(); tx2->SetTextSize(.10);
    tx2->DrawLatexNDC(0.3,0.2, chi2ndf);
  }

//////////////////////////////cent0-10%/////////////////////

  TCanvas *c_chg_MC_fin_cen = new TCanvas("c_chg_MC_fin_cen","",900,450);
  TPad *pad_title = new TPad("pad_title", "",0.,0.92,0.99,0.99);
  TPad *pad_label = new TPad("pad_label", "",0.,0.,0.03,0.95);
  TPad *pad1 = new TPad("pad1", "",0.02,0.35,0.99,0.98);
  TPad *pad2 = new TPad("pad2", "",0.02,0.,0.99,0.35);
  pad1->Divide(3,1);  
  pad1->Draw();
  pad2->Divide(3,1);  
  pad2->Draw();
  pad_title->Draw();
  pad_label->Draw();
  pad_title->cd();
  title_tex->DrawLatexNDC(0.5,0.5, pythia);
  title_tex->DrawLatexNDC(0.05,0.5, cms);  
  pad_label->cd();
  ratio_tex->DrawLatexNDC(0.7,0.1, fit_ratio);
  ncs_tex->DrawLatexNDC(0.7,0.6, ncs_label);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  for(int ibin=1;ibin<4;ibin++){
    pad1->cd(ibin);   
    h_chg_reco[0][ibin]->GetXaxis()->SetTitle("jet charge");
    h_chg_reco[0][ibin]->GetYaxis()->SetTitle("");
    h_chg_reco[0][ibin]->GetYaxis()->SetRangeUser(0.,1.1*h_chg_reco_q[0][ibin]->GetMaximum());
    h_chg_reco[0][ibin]->Draw("e0 same"); h_chg_reco[0][ibin]->Draw("l hist same");
    h_chg_reco_q[0][ibin]->Draw("e0 same"); h_chg_reco_q[0][ibin]->Draw("l hist same");
    h_chg_reco_g[0][ibin]->Draw("e0 same"); h_chg_reco_g[0][ibin]->Draw("l hist same");

    h_chg_reco[0][ibin]->Fit(f_chg_data[0][ibin],"N Q M R","",-2.,2.);

    sprintf(q_frac, "quark %.3f",f_chg_data[0][ibin]->GetParameter(0));
    sprintf(q_error, "#pm %.3f",f_chg_data[0][ibin]->GetParError(0));
    sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_data[0][ibin]->GetChisquare()/f_chg_data[0][ibin]->GetNDF());
    tmp = trk_tag[ibin-1];
    tx = new TLatex(); tx->SetTextSize(.07);
    tx->DrawLatexNDC(0.35,0.25, tmp);

    pad2->cd(ibin);  
    h_chg_reco_ratio[0][ibin]->GetXaxis()->SetTitle("");
    h_chg_reco_ratio[0][ibin]->GetYaxis()->SetTitle("");
    h_chg_reco_ratio[0][ibin]->SetMarkerStyle(7); h_chg_reco_ratio[0][ibin]->Draw("e0 same");
    tl1->Draw("same");
    tx1 = new TLatex(); tx1->SetTextSize(.11);
    tx1->DrawLatexNDC(0.35,0.75, q_frac);
    tx3 = new TLatex(); tx3->SetTextSize(.10);
    tx3->DrawLatexNDC(0.47,0.67, q_error);
    tx2 = new TLatex(); tx2->SetTextSize(.10);
    tx2->DrawLatexNDC(0.3,0.2, chi2ndf);
  }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_chg_data_fin = new TCanvas("c_chg_data_fin","",900,450);
  TPad *pad_title = new TPad("pad_title", "",0.,0.92,0.99,0.99);
  TPad *pad_label = new TPad("pad_label", "",0.,0.,0.03,0.95);
  TPad *pad1 = new TPad("pad1", "",0.02,0.35,0.99,0.98);
  TPad *pad2 = new TPad("pad2", "",0.02,0.,0.99,0.35);
  pad1->Divide(3,1);  
  pad1->Draw();
  pad2->Divide(3,1);  
  pad2->Draw();
  pad_title->Draw();
  pad_label->Draw();
  pad_title->cd();
  title_tex->DrawLatexNDC(0.5,0.5, ppdata);
  title_tex->DrawLatexNDC(0.05,0.5, cms);  
  pad_label->cd();
  ratio_tex->DrawLatexNDC(0.7,0.1, fit_ratio);
  ncs_tex->DrawLatexNDC(0.7,0.6, ncs_label);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  for(int ibin=1;ibin<4;ibin++){
    pad1->cd(ibin);   
    h_chg_data[0][ibin]->GetXaxis()->SetTitle("jet charge");
    h_chg_data[0][ibin]->GetYaxis()->SetTitle("");
    h_chg_data[0][ibin]->GetYaxis()->SetRangeUser(0.,1.1*h_chg_reco_q[0][ibin]->GetMaximum());
    h_chg_data[0][ibin]->Draw("e0 same"); h_chg_data[0][ibin]->Draw("l hist same");
    h_chg_reco_q[0][ibin]->Draw("e0 same"); h_chg_reco_q[0][ibin]->Draw("l hist same");
    h_chg_reco_g[0][ibin]->Draw("e0 same"); h_chg_reco_g[0][ibin]->Draw("l hist same");

    h_chg_data[0][ibin]->Fit(f_chg_data[0][ibin],"N Q M R","",-2.,2.);

    sprintf(q_frac, "quark %.3f",f_chg_data[0][ibin]->GetParameter(0));
    sprintf(q_error, "#pm %.3f",f_chg_data[0][ibin]->GetParError(0));
    sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_data[0][ibin]->GetChisquare()/f_chg_data[0][ibin]->GetNDF());
    tmp = trk_tag[ibin-1];
    tx = new TLatex(); tx->SetTextSize(.07);
    tx->DrawLatexNDC(0.35,0.25, tmp);

    pad2->cd(ibin);  
    h_chg_data_ratio[0][ibin]->GetXaxis()->SetTitle("");
    h_chg_data_ratio[0][ibin]->GetYaxis()->SetTitle("");
    h_chg_data_ratio[0][ibin]->SetMarkerStyle(7); h_chg_data_ratio[0][ibin]->Draw("e0 same");
    tl1->Draw("same");
    tx1 = new TLatex(); tx1->SetTextSize(.11);
    tx1->DrawLatexNDC(0.35,0.75, q_frac);
    tx3 = new TLatex(); tx3->SetTextSize(.10);
    tx3->DrawLatexNDC(0.47,0.67, q_error);
    tx2 = new TLatex(); tx2->SetTextSize(.10);
    tx2->DrawLatexNDC(0.3,0.2, chi2ndf);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_chg_qg_fin = new TCanvas("c_chg_qg_fin","",300,450);
  TPad *pad_title = new TPad("pad_title", "",0.,0.92,0.99,0.99);
  TPad *pad1 = new TPad("pad1", "",0.02,0.35,0.99,0.98);
  TPad *pad2 = new TPad("pad2", "",0.02,0.,0.99,0.35);
  pad1->Divide(1,1);
  pad1->Draw();
  pad2->Divide(1,1);
  pad2->Draw();
  pad_title->Draw();
  pad_title->cd();
  title_tex->DrawLatexNDC(0.5,0.5, pythia);
  title_tex->DrawLatexNDC(0.05,0.5, cms);  

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

    pad1->cd(1);   
    h_chg_ref[0][0]->GetXaxis()->SetTitle("jet charge");
    h_chg_ref[0][0]->GetYaxis()->SetTitle("");
    h_chg_ref[0][0]->GetYaxis()->SetRangeUser(0.,1.1*h_chg_ref_g[0][0]->GetMaximum());
    h_chg_ref[0][0]->Draw("e0 same");
    h_chg_ref_q[0][0]->Draw("e0 same");
    h_chg_ref_g[0][0]->Draw("e0 same");

    pad2->cd(1);  
    h_chg_qg_ratio[0][0] = (TH1D*)h_chg_ref_q[0][0]->Clone("h_chg_ref_qg_ratio");
    h_chg_qg_ratio[0][0]->Divide(h_chg_ref_g[0][0]);
    h_chg_qg_ratio[0][0]->GetXaxis()->SetTitle("");
    h_chg_qg_ratio[0][0]->GetYaxis()->SetTitle("");
    h_chg_qg_ratio[0][0]->SetMarkerStyle(7); h_chg_qg_ratio[0][0]->Draw("e0 same");
    tl1->Draw("same");

}