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
#include <TLine.h>
#include "fitting_templates.h"

#define nCBins 5
#define nptBins 3
#define ntrkbins 4
#define nkbins 4
#define njtptbins 10

Double_t k_scale=1.;

const Double_t pt_low = 120.;
const Double_t pt_high = 600.;

const Double_t ubar_scale = 1.167;
const Double_t dbar_scale = 0.903;

bool do_fit = true;

char saythis[500];

using namespace std;

TString cent[5] = {"0","1","2","3","4"};
TString trk[4] = {"0","1","2","3"};
TString kbin[4] = {"0","1","2","3"};
TString pt[56] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55"};

TString cent_tag[] = {"0-10%","10-30%","30-50%","50-100%"};
TString trk_tag[] = {"p_{T}^{trk} > 0.7 GeV","p_{T}^{trk} > 2 GeV","p_{T}^{trk} > 4 GeV", "p_{T}^{trk} > 5 GeV"};

Double_t jtpt_bound[njtptbins+1] = {120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.};
Double_t x_mean[ntrkbins-1] = {2.,4.,5.};
Double_t x_mean_err[ntrkbins-1] = {0.1,0.1,0.1};
Double_t gluon[nCBins],quark[nCBins],up[nCBins],upbar[nCBins],down[nCBins],downbar[nCBins],gluon_ubar_dbar[nCBins], charm[nCBins], strange[nCBins], bottom[nCBins];
Double_t gluon_gen[nCBins],quark_gen[nCBins], down_gen[nCBins], up_gen[nCBins], gluon_ubar_dbar_gen[nCBins];
Double_t ubdbcsb[nCBins][nptBins];

  TFile *closure_histos_pp_MC = TFile::Open("Pythia6_jetchg_bkgsub_eta0p5_1p5_May23.root");
  TFile *closure_histos_PbPb_MC = TFile::Open("P+H_jetchg_cymbal_bkgsub_eta0p5_1p5_May23.root");
  TFile *closure_histos_pp_data = TFile::Open("ppdata_jetchg_bkgsub_eta0p5_1p5_Apr26.root");
  TFile *closure_histos_PbPb_data = TFile::Open("PbPbdata_jetchg_bkgsub_eta0p5_1p5_Apr30.root");

  TH1D *h_dummy = new TH1D("h_dummy","",5,1,6);
  TH1D *h_dummy2 = new TH1D("h_dummy2","",5,1,6);
  TH1D *h_dummy3 = new TH1D("h_dummy3","",10,120.,220.);
  TH1D *h_dummy4 = new TH1D("h_dummy4","",10,120.,220.);

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

  TH2F *h_chg_refpt_up[nCBins][ntrkbins];
  TH2F *h_chg_refpt_down[nCBins][ntrkbins];
  TH2F *h_chg_refpt_upb[nCBins][ntrkbins];
  TH2F *h_chg_refpt_downb[nCBins][ntrkbins];

  TH2F *h_chg_corrpt_c[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_s[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_b[nCBins][ntrkbins];

  TH2F *h_chg_corrpt_up[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_down[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_upb[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_downb[nCBins][ntrkbins];

  TH2F *h_chg_corrpt_bkg[nCBins][ntrkbins];

  TH2F *h_chg_corrpt_data[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_data_bkg[nCBins][ntrkbins];

  TH1D *h_chg_ref[nCBins][ntrkbins];
  //TH1D *h_chg_ref_q[nCBins][ntrkbins];
  //TH1D *h_chg_ref_g[nCBins][ntrkbins];
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

  TH1D *h_reco_up_inc_ratio[nCBins];
  TH1D *h_reco_upbar_inc_ratio[nCBins];  
  TH1D *h_reco_d_inc_ratio[nCBins];
  TH1D *h_reco_dbar_inc_ratio[nCBins];
  TH1D *h_reco_gluon_inc_ratio[nCBins];

  TH1D *h_reco_c_inc_ratio[nCBins];
  TH1D *h_reco_s_inc_ratio[nCBins];  
  TH1D *h_reco_b_inc_ratio[nCBins];

  TH1D *h_reco_gluon_upbar_ratio[nCBins];
  TH1D *h_reco_gluon_dbar_ratio[nCBins];
  TH1D *h_reco_gluon_c_ratio[nCBins];
  TH1D *h_reco_gluon_s_ratio[nCBins];
  TH1D *h_reco_gluon_b_ratio[nCBins];

  TH1D *h_reco_MC[nCBins];
  TH1D *h_reco_data[nCBins];
  TH1D *h_reco_MC_q[nCBins];
  TH1D *h_reco_MC_g[nCBins];
  TH1D *h_reco_MC_up[nCBins];
  TH1D *h_reco_MC_upbar[nCBins];
  TH1D *h_reco_MC_d[nCBins];
  TH1D *h_reco_MC_dbar[nCBins];

  TH1D *h_reco_MC_c[nCBins];
  TH1D *h_reco_MC_s[nCBins];
  TH1D *h_reco_MC_b[nCBins];

  TH1D *h_ref_MC[nCBins];
  TH1D *h_ref_MC_q[nCBins];
  TH1D *h_ref_MC_g[nCBins];
  TH1D *h_ref_MC_up[nCBins];
  TH1D *h_ref_MC_upbar[nCBins];
  TH1D *h_ref_MC_d[nCBins];
  TH1D *h_ref_MC_dbar[nCBins];

  TH1D *h_chg_reco[nCBins][ntrkbins];
  TH1D *h_chg_reco_q[nCBins][ntrkbins];
  //TH1D *h_chg_reco_g[nCBins][ntrkbins];
  TH1D *h_chg_reco_u[nCBins][ntrkbins];
  TH1D *h_chg_reco_bkg[nCBins][ntrkbins];

  //TH1D *h_chg_reco_up[nCBins][ntrkbins];
  //TH1D *h_chg_reco_down[nCBins][ntrkbins];
  TH1D *h_chg_reco_upb[nCBins][ntrkbins];
  TH1D *h_chg_reco_downb[nCBins][ntrkbins];

  //TH1D *h_chg_ref_up[nCBins][ntrkbins];
  //TH1D *h_chg_ref_down[nCBins][ntrkbins];
  //TH1D *h_chg_ref_upb[nCBins][ntrkbins];
  //TH1D *h_chg_ref_downb[nCBins][ntrkbins];

  //TH1D *h_chg_reco_gubdb[nCBins][ntrkbins];

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
/*
  TProfile *h_chg_reco_avg[nCBins][ntrkbins];
  TProfile *h_chg_reco_q_avg[nCBins][ntrkbins];
  TProfile *h_chg_reco_g_avg[nCBins][ntrkbins];

  TProfile *h_chg_reco_up_avg[nCBins][ntrkbins];
  TProfile *h_chg_reco_down_avg[nCBins][ntrkbins];
  TProfile *h_chg_reco_upb_avg[nCBins][ntrkbins];
  TProfile *h_chg_reco_downb_avg[nCBins][ntrkbins];
*/
  TProfile *h_chg_reco_bkg_avg[nCBins][ntrkbins];

  TProfile *h_chg_data_avg[nCBins][ntrkbins];
  TProfile *h_chg_data_bkg_avg[nCBins][ntrkbins];

  TGraphErrors *up_reco[nCBins];
  TGraphErrors *down_reco[nCBins];
  TGraphErrors *gluon_reco[nCBins];
  TGraphErrors *quark_reco[nCBins];

  TGraphErrors *up_ref[nCBins];
  TGraphErrors *down_ref[nCBins];
  TGraphErrors *gluon_ref[nCBins];

  TGraphErrors *up_data[nCBins];
  TGraphErrors *down_data[nCBins];
  TGraphErrors *gluon_data[nCBins];
  TGraphErrors *quark_data[nCBins];

  TF1 *f_chg_ref[nCBins][ntrkbins];  
  TF1 *f_chg_reco[nCBins][ntrkbins];
  TF1 *f_chg_udg_data[nCBins][ntrkbins];
  TF1 *f_chg_udg_ref[nCBins][ntrkbins];

  TF1 *f_eta_reco[nCBins];

  TH1D *h_chg_reco_avg[nCBins][njtptbins][ntrkbins];
  TH1D *h_chg_reco_q_avg[nCBins][njtptbins][ntrkbins];
  TH1D *h_chg_reco_g_avg[nCBins][njtptbins][ntrkbins];
  TH1D *h_chg_reco_up_avg[nCBins][njtptbins][ntrkbins];
  TH1D *h_chg_reco_down_avg[nCBins][njtptbins][ntrkbins];

  TF1 *f_chg_reco_avg[nCBins][njtptbins][ntrkbins];  
  TF1 *f_chg_reco_q_avg[nCBins][njtptbins][ntrkbins];
  TF1 *f_chg_reco_g_avg[nCBins][njtptbins][ntrkbins];
  TF1 *f_chg_reco_up_avg[nCBins][njtptbins][ntrkbins];
  TF1 *f_chg_reco_down_avg[nCBins][njtptbins][ntrkbins];

void get_histos(){

  for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0){
      h_eta_MC[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_cent"+cent[ibin]));       
      h_eta_MC_q[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_q_cent"+cent[ibin]));
      h_eta_MC_g[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_g_cent"+cent[ibin]));
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

      h_reco_MC_c[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_c_cent"+cent[ibin]));
      h_reco_MC_s[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_s_cent"+cent[ibin]));      
      h_reco_MC_b[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_b_cent"+cent[ibin]));

      h_ref_MC[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_gen_full_cent"+cent[ibin]));
      h_ref_MC_q[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_gen_full_q_cent"+cent[ibin]));
      h_ref_MC_g[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_gen_full_g_cent"+cent[ibin]));
      h_ref_MC_up[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_gen_full_up_cent"+cent[ibin]));
      h_ref_MC_upbar[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_gen_full_upbar_cent"+cent[ibin]));
      h_ref_MC_d[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_gen_full_d_cent"+cent[ibin]));
      h_ref_MC_dbar[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_gen_full_dbar_cent"+cent[ibin]));
    }
    else{
      h_eta_MC[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_cent"+cent[ibin-1]));       
      h_eta_MC_q[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_q_cent"+cent[ibin-1]));
      h_eta_MC_g[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_g_cent"+cent[ibin-1]));
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

      h_reco_MC_c[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_c_cent"+cent[ibin-1]));
      h_reco_MC_s[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_s_cent"+cent[ibin-1]));
      h_reco_MC_b[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_b_cent"+cent[ibin-1]));

      h_ref_MC[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_gen_full_cent"+cent[ibin-1]));
      h_ref_MC_q[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_gen_full_q_cent"+cent[ibin-1]));
      h_ref_MC_g[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_gen_full_g_cent"+cent[ibin-1]));
      h_ref_MC_up[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_gen_full_up_cent"+cent[ibin-1]));
      h_ref_MC_upbar[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_gen_full_upbar_cent"+cent[ibin-1]));
      h_ref_MC_d[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_gen_full_d_cent"+cent[ibin-1]));
      h_ref_MC_dbar[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_gen_full_dbar_cent"+cent[ibin-1]));

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
  }

  for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    for(int ibin=0;ibin<nCBins;ibin++){
      if(ibin==0){
        h_chg_refpt[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
        h_chg_refpt_q[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_q_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_refpt_g[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_g_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_refpt_u[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_u_cent"+cent[ibin]+"_trk"+trk[ibin2]));

        h_chg_refpt_up[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_up_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_refpt_down[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_d_cent"+cent[ibin]+"_trk"+trk[ibin2]));      
        h_chg_refpt_upb[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_upbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_refpt_downb[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_dbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));

        h_chg_corrpt[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
        h_chg_corrpt_q[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_q_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_corrpt_g[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_g_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_corrpt_u[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_u_cent"+cent[ibin]+"_trk"+trk[ibin2]));

        h_chg_corrpt_up[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_up_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_corrpt_down[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_d_cent"+cent[ibin]+"_trk"+trk[ibin2]));      
        h_chg_corrpt_upb[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_upbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_corrpt_downb[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_dbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));

        h_chg_corrpt_c[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_c_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_corrpt_s[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_s_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_corrpt_b[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_b_cent"+cent[ibin]+"_trk"+trk[ibin2]));

        h_chg_corrpt_data[ibin][ibin2] = (TH2F*)closure_histos_pp_data->Get((TString)("h_chg_corrpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 

        h_chg_corrpt_bkg[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_bkg_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
        h_chg_corrpt_data_bkg[ibin][ibin2] = (TH2F*)closure_histos_pp_data->Get((TString)("h_chg_corrpt_bkg_cent"+cent[ibin]+"_trk"+trk[ibin2]));
      }
      else{
        h_chg_refpt[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2]))->Clone((TString)("h_chg_refpt_"+cent[ibin-1]+"_"+trk[ibin2])); 
        h_chg_refpt_q[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_q_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_refpt_g[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_g_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_refpt_u[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_u_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

        h_chg_refpt_up[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_up_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_refpt_down[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_d_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));      
        h_chg_refpt_upb[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_upbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_refpt_downb[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_dbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

        h_chg_corrpt[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
        h_chg_corrpt_q[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_q_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_corrpt_g[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_g_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_corrpt_u[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_u_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

        h_chg_corrpt_up[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_up_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_corrpt_down[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_d_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_corrpt_upb[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_upbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_corrpt_downb[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_dbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

        h_chg_corrpt_c[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_c_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_corrpt_s[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_s_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_chg_corrpt_b[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_b_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

        h_chg_corrpt_data[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_chg_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 

        h_chg_corrpt_bkg[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_bkg_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
        h_chg_corrpt_data_bkg[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_chg_corrpt_bkg_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
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
        h_chg_reco_c[ibin][ibin2] = h_chg_corrpt_c[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_c_%d_%d",ibin,ibin2),h_chg_corrpt_c[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_c[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_reco_s[ibin][ibin2] = h_chg_corrpt_s[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_s_%d_%d",ibin,ibin2),h_chg_corrpt_s[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_s[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_reco_b[ibin][ibin2] = h_chg_corrpt_b[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_b_%d_%d",ibin,ibin2),h_chg_corrpt_b[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_b[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
/*
        if(ibin!=0){
          h_chg_reco_upb[ibin][ibin2]->Scale(ubar_scale);
          h_chg_reco_downb[ibin][ibin2]->Scale(dbar_scale);
        }
*/
        //h_chg_reco_upb[ibin][ibin2]->Scale(k_scale);
        //h_chg_reco_downb[ibin][ibin2]->Scale(k_scale);
        //h_chg_reco_c[ibin][ibin2]->Scale(k_scale);
        //h_chg_reco_s[ibin][ibin2]->Scale(k_scale);
        //h_chg_reco_b[ibin][ibin2]->Scale(k_scale);

        int low_bin = h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(pt_low); int high_bin = h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(pt_high);

        ubdbcsb[ibin][ibin2] = h_chg_reco_upb[ibin][ibin2]->Integral()+h_chg_reco_downb[ibin][ibin2]->Integral()+h_chg_reco_c[ibin][ibin2]->Integral()+h_chg_reco_s[ibin][ibin2]->Integral()+h_chg_reco_b[ibin][ibin2]->Integral();
        ubdbcsb[ibin][ibin2] = ubdbcsb[ibin][ibin2]/h_chg_reco[ibin][ibin2]->Integral();
        //cout<<ubdbcsb[ibin][ibin2]<<endl;         

        h_chg_reco_gubdb[ibin][ibin2] = (TH1D*)h_chg_reco_g[ibin][ibin2]->Clone(Form("h_chg_reco_gubdb_%d_%d",ibin,ibin2));
        h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_upb[ibin][ibin2]);
        h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_downb[ibin][ibin2]);
        h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_c[ibin][ibin2]);
        h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_s[ibin][ibin2]);
        h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_b[ibin][ibin2]);

        //normalizing chg histos
        h_chg_reco[ibin][ibin2]->Scale(1./(h_chg_reco[ibin][ibin2]->Integral()));
        h_chg_reco_q[ibin][ibin2]->Scale(1./(h_chg_reco_q[ibin][ibin2]->Integral()));
        h_chg_reco_g[ibin][ibin2]->Scale(1./(h_chg_reco_g[ibin][ibin2]->Integral()));
        h_chg_reco_u[ibin][ibin2]->Scale(1./(h_chg_reco_u[ibin][ibin2]->Integral()));                        
        h_chg_reco_up[ibin][ibin2]->Scale(1./(h_chg_reco_up[ibin][ibin2]->Integral()));
        h_chg_reco_down[ibin][ibin2]->Scale(1./(h_chg_reco_down[ibin][ibin2]->Integral()));
        h_chg_reco_upb[ibin][ibin2]->Scale(1./(h_chg_reco_upb[ibin][ibin2]->Integral()));        
        h_chg_reco_downb[ibin][ibin2]->Scale(1./(h_chg_reco_downb[ibin][ibin2]->Integral()));
        h_chg_reco_gubdb[ibin][ibin2]->Scale(1./(h_chg_reco_gubdb[ibin][ibin2]->Integral()));
        h_chg_data[ibin][ibin2]->Scale(1./(h_chg_data[ibin][ibin2]->Integral()));

        for(int ibin3=0; ibin3<njtptbins; ibin3++){
          h_chg_reco_avg[ibin][ibin3][ibin2] = h_chg_corrpt[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_q_avg[ibin][ibin3][ibin2] = h_chg_corrpt_q[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_q_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_g_avg[ibin][ibin3][ibin2] = h_chg_corrpt_g[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_g_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_up_avg[ibin][ibin3][ibin2] = h_chg_corrpt_up[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_up_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_down_avg[ibin][ibin3][ibin2] = h_chg_corrpt_down[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_down_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
        
          f_chg_reco_avg[ibin][ibin3][ibin2] = new TF1(Form("f_chg_%d_%d_%d",ibin,ibin2,ibin3),"gaus",-2.,2.);
          f_chg_reco_q_avg[ibin][ibin3][ibin2] = new TF1(Form("f_chg_q_%d_%d_%d",ibin,ibin2,ibin3),"gaus",-2.,2.);
          f_chg_reco_g_avg[ibin][ibin3][ibin2] = new TF1(Form("f_chg_g_%d_%d_%d",ibin,ibin2,ibin3),"gaus",-2.,2.);
          f_chg_reco_up_avg[ibin][ibin3][ibin2] = new TF1(Form("f_chg_up_%d_%d_%d",ibin,ibin2,ibin3),"gaus",-2.,2.);
          f_chg_reco_down_avg[ibin][ibin3][ibin2] = new TF1(Form("f_chg_down_%d_%d_%d",ibin,ibin2,ibin3),"gaus",-2.,2.);
        }


        h_chg_ref[ibin][ibin2] = h_chg_refpt[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_%d_%d",ibin,ibin2),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_g[ibin][ibin2] = h_chg_refpt_g[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_g_%d_%d",ibin,ibin2),h_chg_refpt_g[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_g[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_up[ibin][ibin2] = h_chg_refpt_up[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_up_%d_%d",ibin,ibin2),h_chg_refpt_up[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_up[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_down[ibin][ibin2] = h_chg_refpt_down[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_down_%d_%d",ibin,ibin2),h_chg_refpt_down[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_down[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_upb[ibin][ibin2] = h_chg_refpt_upb[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_upb_%d_%d",ibin,ibin2),h_chg_refpt_upb[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_upb[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_downb[ibin][ibin2] = h_chg_refpt_downb[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_downb_%d_%d",ibin,ibin2),h_chg_refpt_downb[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_downb[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");

        h_chg_ref_gubdb[ibin][ibin2] = (TH1D*)h_chg_ref_g[ibin][ibin2]->Clone(Form("h_chg_ref_gubdb_%d_%d",ibin,ibin2));
        h_chg_ref_gubdb[ibin][ibin2]->Add(h_chg_ref_upb[ibin][ibin2]);
        h_chg_ref_gubdb[ibin][ibin2]->Add(h_chg_ref_downb[ibin][ibin2]);

        h_chg_ref[ibin][ibin2]->Scale(1./(h_chg_ref[ibin][ibin2]->Integral()));
        h_chg_ref_g[ibin][ibin2]->Scale(1./(h_chg_ref_g[ibin][ibin2]->Integral()));
        h_chg_ref_up[ibin][ibin2]->Scale(1./(h_chg_ref_up[ibin][ibin2]->Integral()));
        h_chg_ref_down[ibin][ibin2]->Scale(1./(h_chg_ref_down[ibin][ibin2]->Integral()));
        h_chg_ref_gubdb[ibin][ibin2]->Scale(1./(h_chg_ref_gubdb[ibin][ibin2]->Integral()));

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
/*
        h_chg_reco_avg[ibin][ibin2] = h_chg_corrpt[ibin][ibin2]-> ProfileX(Form("h_chg_reco_avg_%d_%d",ibin,ibin2),h_chg_corrpt[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_reco_q_avg[ibin][ibin2] = h_chg_corrpt_q[ibin][ibin2]-> ProfileX(Form("h_chg_reco_q_avg_%d_%d",ibin,ibin2),h_chg_corrpt_q[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_q[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_reco_g_avg[ibin][ibin2] = h_chg_corrpt_g[ibin][ibin2]-> ProfileX(Form("h_chg_reco_g_avg_%d_%d",ibin,ibin2),h_chg_corrpt_g[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_g[ibin][ibin2]->GetYaxis()->FindBin(10),"");

        h_chg_reco_up_avg[ibin][ibin2] = h_chg_corrpt_up[ibin][ibin2]-> ProfileX(Form("h_chg_reco_up_avg_%d_%d",ibin,ibin2),h_chg_corrpt_up[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_up[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_reco_down_avg[ibin][ibin2] = h_chg_corrpt_down[ibin][ibin2]-> ProfileX(Form("h_chg_reco_down_avg_%d_%d",ibin,ibin2),h_chg_corrpt_down[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_down[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_reco_upb_avg[ibin][ibin2] = h_chg_corrpt_upb[ibin][ibin2]-> ProfileX(Form("h_chg_reco_upb_avg_%d_%d",ibin,ibin2),h_chg_corrpt_upb[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_upb[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_reco_downb_avg[ibin][ibin2] = h_chg_corrpt_downb[ibin][ibin2]-> ProfileX(Form("h_chg_reco_downb_avg_%d_%d",ibin,ibin2),h_chg_corrpt_downb[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_downb[ibin][ibin2]->GetYaxis()->FindBin(10),"");
*/
        h_chg_data_avg[ibin][ibin2] = h_chg_corrpt_data[ibin][ibin2]-> ProfileX(Form("h_chg_data_avg_%d_%d",ibin,ibin2),h_chg_corrpt_data[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_data[ibin][ibin2]->GetYaxis()->FindBin(10),"");

        h_chg_reco_bkg_avg[ibin][ibin2] = h_chg_corrpt_bkg[ibin][ibin2]-> ProfileX(Form("h_chg_bkg_reco_avg_%d_%d",ibin,ibin2),h_chg_corrpt_bkg[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_bkg[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        h_chg_data_bkg_avg[ibin][ibin2] = h_chg_corrpt_data_bkg[ibin][ibin2]-> ProfileX(Form("h_chg_bkg_data_avg_%d_%d",ibin,ibin2),h_chg_corrpt_data_bkg[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_data_bkg[ibin][ibin2]->GetYaxis()->FindBin(10),"");
        
        h_chg_reco_bkg[ibin][ibin2] = h_chg_corrpt_bkg[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_bkg_%d_%d",ibin,ibin2),h_chg_corrpt_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_reco_bkg[ibin][ibin2]->Scale(1./(h_chg_reco_bkg[ibin][ibin2]->Integral()));
        h_chg_data_bkg[ibin][ibin2] = h_chg_corrpt_data_bkg[ibin][ibin2]-> ProjectionY(Form("h_chg_data_bkg_%d_%d",ibin,ibin2),h_chg_corrpt_data_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_data_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_data_bkg[ibin][ibin2]->Scale(1./(h_chg_data_bkg[ibin][ibin2]->Integral()));
    }
  }
}

void plot_charge(){

  get_histos();

  TLine *tl_q[nCBins]; TLine *tl_g[nCBins]; TLine *tl_u[nCBins]; TLine *tl_d[nCBins];
  TLine *tl_g_data; TLine *tl_u_data; TLine *tl_d_data;
  TLine *tl_ubar[nCBins]; TLine *tl_dbar[nCBins]; TLine *tl_c[nCBins]; TLine *tl_s[nCBins]; TLine *tl_b[nCBins];
  TLine *tl_g_ref[nCBins]; TLine *tl_q_ref[nCBins]; TLine *tl_u_ref[nCBins]; TLine *tl_d_ref[nCBins];

  ////////// direct MC fraction from spectra
  for(int ibin=0;ibin<nCBins;ibin++){
      int low_bin = h_reco_MC[ibin]->FindBin(pt_low); int high_bin = h_reco_MC[ibin]->FindBin(pt_high); 

      gluon_gen[ibin] = h_ref_MC_g[ibin]->Integral(low_bin,high_bin)/h_ref_MC[ibin]->Integral(low_bin,high_bin);
      quark_gen[ibin] = h_ref_MC_q[ibin]->Integral(low_bin,high_bin)/h_ref_MC[ibin]->Integral(low_bin,high_bin);
      up_gen[ibin] = h_ref_MC_up[ibin]->Integral(low_bin,high_bin)/h_ref_MC[ibin]->Integral(low_bin,high_bin);
      down_gen[ibin] = h_ref_MC_d[ibin]->Integral(low_bin,high_bin)/h_ref_MC[ibin]->Integral(low_bin,high_bin);
      gluon_ubar_dbar_gen[ibin] = (h_ref_MC_g[ibin]->Integral(low_bin,high_bin)+h_ref_MC_upbar[ibin]->Integral(low_bin,high_bin)+h_ref_MC_dbar[ibin]->Integral(low_bin,high_bin))/h_ref_MC[ibin]->Integral(low_bin,high_bin);

      gluon[ibin] = h_reco_MC_g[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
      quark[ibin] = h_reco_MC_q[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
      up[ibin] = h_reco_MC_up[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
      down[ibin] = h_reco_MC_d[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
      upbar[ibin] = h_reco_MC_upbar[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
      downbar[ibin] = h_reco_MC_dbar[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
      charm[ibin] = h_reco_MC_c[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
      strange[ibin] = h_reco_MC_s[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
      bottom[ibin] = h_reco_MC_b[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
      gluon_ubar_dbar[ibin] = (h_reco_MC_g[ibin]->Integral(low_bin,high_bin)+h_reco_MC_upbar[ibin]->Integral(low_bin,high_bin)+h_reco_MC_dbar[ibin]->Integral(low_bin,high_bin)+h_reco_MC_c[ibin]->Integral(low_bin,high_bin)+h_reco_MC_s[ibin]->Integral(low_bin,high_bin)+h_reco_MC_b[ibin]->Integral(low_bin,high_bin))/h_reco_MC[ibin]->Integral(low_bin,high_bin);

      cout<<"up: "<<up[ibin]<<"  down:  "<<down[ibin]<<"  g_ub_db:  "<<gluon_ubar_dbar[ibin]<<endl;
      cout<<"charm: "<<charm[ibin]<<"  strange:  "<<strange[ibin]<<"  bottom:  "<<bottom[ibin]<<endl;

      cout<<"total: "<<gluon_ubar_dbar[ibin]+up[ibin]+down[ibin]<<endl;
  
      tl_u[ibin] = new TLine(1.,up[ibin],6.,up[ibin]); tl_u[ibin]->SetLineStyle(2); tl_u[ibin]->SetLineColor(kBlue);
      tl_d[ibin] = new TLine(1.,down[ibin],6.,down[ibin]); tl_d[ibin]->SetLineStyle(9); tl_d[ibin]->SetLineColor(kBlue);
      tl_q[ibin] = new TLine(1.,quark[ibin],6.,quark[ibin]); tl_q[ibin]->SetLineStyle(1); tl_q[ibin]->SetLineColor(kBlue);
      tl_g[ibin] = new TLine(1.,gluon[ibin],6.,gluon[ibin]); tl_g[ibin]->SetLineStyle(1); tl_g[ibin]->SetLineColor(kRed);

      tl_ubar[ibin] = new TLine(1.,upbar[ibin],6.,upbar[ibin]); tl_ubar[ibin]->SetLineStyle(2); tl_ubar[ibin]->SetLineColor(kBlue);
      tl_dbar[ibin] = new TLine(1.,downbar[ibin],6.,downbar[ibin]); tl_dbar[ibin]->SetLineStyle(2); tl_dbar[ibin]->SetLineColor(kBlue);
      tl_c[ibin] = new TLine(1.,charm[ibin],6.,charm[ibin]); tl_c[ibin]->SetLineStyle(2); tl_c[ibin]->SetLineColor(kBlue);
      tl_s[ibin] = new TLine(1.,strange[ibin],6.,strange[ibin]); tl_s[ibin]->SetLineStyle(2); tl_s[ibin]->SetLineColor(kBlue);
      tl_b[ibin] = new TLine(1.,bottom[ibin],6.,bottom[ibin]); tl_b[ibin]->SetLineStyle(2); tl_b[ibin]->SetLineColor(kBlue);

      tl_q_ref[ibin] = new TLine(1.,quark_gen[ibin],6.,quark_gen[ibin]); tl_q_ref[ibin]->SetLineStyle(2); tl_q_ref[ibin]->SetLineColor(kBlue);
      tl_g_ref[ibin] = new TLine(1.,gluon_ubar_dbar_gen[ibin],6.,gluon_ubar_dbar_gen[ibin]); tl_g_ref[ibin]->SetLineStyle(2); tl_g_ref[ibin]->SetLineColor(kRed);
      tl_u_ref[ibin] = new TLine(1.,up_gen[ibin],6.,up_gen[ibin]); tl_u_ref[ibin]->SetLineStyle(2); tl_u_ref[ibin]->SetLineColor(kBlue);
      tl_d_ref[ibin] = new TLine(1.,down_gen[ibin],6.,down_gen[ibin]); tl_d_ref[ibin]->SetLineStyle(9); tl_d_ref[ibin]->SetLineColor(kBlue);
  }

  //tl_g_data = new TLine(1.,gluon[0],6.,gluon[0]); tl_g_data->SetLineStyle(2); tl_g_data->SetLineColor(kRed);
  tl_u_data = new TLine(1.,0.134,6.,0.134); tl_u_data->SetLineStyle(2); tl_u_data->SetLineColor(kBlue);
  tl_d_data = new TLine(1.,0.154,6.,0.154); tl_d_data->SetLineStyle(9); tl_d_data->SetLineColor(kBlue);

  f_chg_udg_data[0][1] = new TF1("f_chg_udg_data_cent0_trkpt2",f_udg_reco_cent0_pt2,-2.,2.,2);
  f_chg_udg_data[1][1] = new TF1("f_chg_udg_data_cent1_trkpt2",f_udg_reco_cent1_pt2,-2.,2.,2);
  f_chg_udg_data[2][1] = new TF1("f_chg_udg_data_cent2_trkpt2",f_udg_reco_cent2_pt2,-2.,2.,2);
  f_chg_udg_data[3][1] = new TF1("f_chg_udg_data_cent3_trkpt2",f_udg_reco_cent3_pt2,-2.,2.,2);
  f_chg_udg_data[4][1] = new TF1("f_chg_udg_data_cent4_trkpt2",f_udg_reco_cent4_pt2,-2.,2.,2);

  f_chg_udg_data[0][2] = new TF1("f_chg_udg_data_cent0_trkpt4",f_udg_reco_cent0_pt4,-2.,2.,2);
  f_chg_udg_data[1][2] = new TF1("f_chg_udg_data_cent1_trkpt4",f_udg_reco_cent1_pt4,-2.,2.,2);
  f_chg_udg_data[2][2] = new TF1("f_chg_udg_data_cent2_trkpt4",f_udg_reco_cent2_pt4,-2.,2.,2);
  f_chg_udg_data[3][2] = new TF1("f_chg_udg_data_cent3_trkpt4",f_udg_reco_cent3_pt4,-2.,2.,2);
  f_chg_udg_data[4][2] = new TF1("f_chg_udg_data_cent4_trkpt4",f_udg_reco_cent4_pt4,-2.,2.,2);

  f_chg_udg_data[0][3] = new TF1("f_chg_udg_data_cent0_trkpt5",f_udg_reco_cent0_pt5,-2.,2.,2);
  f_chg_udg_data[1][3] = new TF1("f_chg_udg_data_cent1_trkpt5",f_udg_reco_cent1_pt5,-2.,2.,2);
  f_chg_udg_data[2][3] = new TF1("f_chg_udg_data_cent2_trkpt5",f_udg_reco_cent2_pt5,-2.,2.,2);
  f_chg_udg_data[3][3] = new TF1("f_chg_udg_data_cent3_trkpt5",f_udg_reco_cent3_pt5,-2.,2.,2);
  f_chg_udg_data[4][3] = new TF1("f_chg_udg_data_cent4_trkpt5",f_udg_reco_cent4_pt5,-2.,2.,2);

  f_chg_udg_ref[0][1] = new TF1("f_chg_udg_ref_cent0_trkpt2",f_udg_ref_cent0_pt2,-2.,2.,2);
  f_chg_udg_ref[1][1] = new TF1("f_chg_udg_ref_cent1_trkpt2",f_udg_ref_cent1_pt2,-2.,2.,2);
  f_chg_udg_ref[2][1] = new TF1("f_chg_udg_ref_cent2_trkpt2",f_udg_ref_cent2_pt2,-2.,2.,2);
  f_chg_udg_ref[3][1] = new TF1("f_chg_udg_ref_cent3_trkpt2",f_udg_ref_cent3_pt2,-2.,2.,2);
  f_chg_udg_ref[4][1] = new TF1("f_chg_udg_ref_cent4_trkpt2",f_udg_ref_cent4_pt2,-2.,2.,2);

  f_chg_udg_ref[0][2] = new TF1("f_chg_udg_ref_cent0_trkpt4",f_udg_ref_cent0_pt4,-2.,2.,2);
  f_chg_udg_ref[1][2] = new TF1("f_chg_udg_ref_cent1_trkpt4",f_udg_ref_cent1_pt4,-2.,2.,2);
  f_chg_udg_ref[2][2] = new TF1("f_chg_udg_ref_cent2_trkpt4",f_udg_ref_cent2_pt4,-2.,2.,2);
  f_chg_udg_ref[3][2] = new TF1("f_chg_udg_ref_cent3_trkpt4",f_udg_ref_cent3_pt4,-2.,2.,2);
  f_chg_udg_ref[4][2] = new TF1("f_chg_udg_ref_cent4_trkpt4",f_udg_ref_cent4_pt4,-2.,2.,2);

  f_chg_udg_ref[0][3] = new TF1("f_chg_udg_ref_cent0_trkpt5",f_udg_ref_cent0_pt5,-2.,2.,2);
  f_chg_udg_ref[1][3] = new TF1("f_chg_udg_ref_cent1_trkpt5",f_udg_ref_cent1_pt5,-2.,2.,2);
  f_chg_udg_ref[2][3] = new TF1("f_chg_udg_ref_cent2_trkpt5",f_udg_ref_cent2_pt5,-2.,2.,2);
  f_chg_udg_ref[3][3] = new TF1("f_chg_udg_ref_cent3_trkpt5",f_udg_ref_cent3_pt5,-2.,2.,2);
  f_chg_udg_ref[4][3] = new TF1("f_chg_udg_ref_cent4_trkpt5",f_udg_ref_cent4_pt5,-2.,2.,2);

  TLine *tl1 = new TLine(-1.5,1.,1.5,1.);
  tl1->SetLineStyle(2);
  TLine *tl2 = new TLine(120.,1.,400.,1.);
  tl2->SetLineStyle(2);
  TLatex *tx, *tx1, *tx2, *tx3, *tx4, *tx5;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_qg_eta = new TCanvas("c_qg_eta","",1500,300);
  c_qg_eta->Divide(5,1,0);
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1);

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
    //h_eta_MC_g[ibin]->GetYaxis()->SetRangeUser(0.02,0.12);
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

  TCanvas *c_reco = new TCanvas("c_reco","",1500,600);
  c_reco->Divide(5,2,0);
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1);

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

    h_reco_MC_c[ibin]->SetLineColor(kBlue); h_reco_MC_c[ibin]->SetMarkerColor(kBlue);
    h_reco_MC_c[ibin]->SetMarkerStyle(28); h_reco_MC_c[ibin]->Draw("e0 same");
    h_reco_MC_s[ibin]->SetLineColor(kBlue); h_reco_MC_s[ibin]->SetMarkerColor(kBlue);
    h_reco_MC_s[ibin]->SetMarkerStyle(29); h_reco_MC_s[ibin]->Draw("e0 same");
    h_reco_MC_b[ibin]->SetLineColor(kBlue); h_reco_MC_b[ibin]->SetMarkerColor(kBlue);
    h_reco_MC_b[ibin]->SetMarkerStyle(30); h_reco_MC_b[ibin]->Draw("e0 same");
    
    if(ibin>0){
      TString tmp1 = cent_tag[ibin-1];
      tx1 = new TLatex(); tx1->SetTextSize(.12);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }

    if(ibin==0) c_reco->cd(ibin+6);
    else c_reco->cd(11-ibin);   
    h_reco_gluon_inc_ratio[ibin] = (TH1D*)h_reco_MC_g[ibin]->Clone("h_reco_gluon_inc_ratio");
    h_reco_gluon_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
    h_reco_up_inc_ratio[ibin] = (TH1D*)h_reco_MC_up[ibin]->Clone("h_reco_up_inc_ratio");
    h_reco_up_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
    h_reco_upbar_inc_ratio[ibin] = (TH1D*)h_reco_MC_upbar[ibin]->Clone("h_reco_upbar_inc_ratio");
    h_reco_upbar_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
    h_reco_d_inc_ratio[ibin] = (TH1D*)h_reco_MC_d[ibin]->Clone("h_reco_d_inc_ratio");
    h_reco_d_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
    h_reco_dbar_inc_ratio[ibin] = (TH1D*)h_reco_MC_dbar[ibin]->Clone("h_reco_dbar_inc_ratio");
    h_reco_dbar_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
    h_reco_c_inc_ratio[ibin] = (TH1D*)h_reco_MC_c[ibin]->Clone("h_reco_c_inc_ratio");
    h_reco_c_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
    h_reco_s_inc_ratio[ibin] = (TH1D*)h_reco_MC_s[ibin]->Clone("h_reco_s_inc_ratio");
    h_reco_s_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
    h_reco_b_inc_ratio[ibin] = (TH1D*)h_reco_MC_b[ibin]->Clone("h_reco_b_inc_ratio");
    h_reco_b_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);

    h_reco_gluon_inc_ratio[ibin]->GetYaxis()->SetNdivisions(505);
    h_reco_gluon_inc_ratio[ibin]->GetYaxis()->SetTitleSize(0.05);
    h_reco_gluon_inc_ratio[ibin]->GetYaxis()->SetLabelSize(0.07);
    h_reco_gluon_inc_ratio[ibin]->GetXaxis()->CenterTitle();
    h_reco_gluon_inc_ratio[ibin]->GetXaxis()->SetNdivisions(505);
    h_reco_gluon_inc_ratio[ibin]->GetXaxis()->SetTitleSize(0.05);
    h_reco_gluon_inc_ratio[ibin]->GetXaxis()->SetLabelSize(0.07);
    h_reco_gluon_inc_ratio[ibin]->GetXaxis()->SetRangeUser(120.,400.);
    h_reco_gluon_inc_ratio[ibin]->GetYaxis()->SetRangeUser(0.4,1.6);
    h_reco_gluon_inc_ratio[ibin]->SetLineColor(kRed); h_reco_gluon_inc_ratio[ibin]->SetMarkerColor(kRed);
    h_reco_gluon_inc_ratio[ibin]->SetMarkerStyle(20); h_reco_gluon_inc_ratio[ibin]->Draw("e0 same");
    h_reco_up_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_up_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_up_inc_ratio[ibin]->SetMarkerStyle(22); h_reco_up_inc_ratio[ibin]->Draw("e0 same");
    h_reco_upbar_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_upbar_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_upbar_inc_ratio[ibin]->SetMarkerStyle(26); h_reco_upbar_inc_ratio[ibin]->Draw("e0 same");
    h_reco_d_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_d_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_d_inc_ratio[ibin]->SetMarkerStyle(23); h_reco_d_inc_ratio[ibin]->Draw("e0 same");
    h_reco_dbar_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_dbar_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_dbar_inc_ratio[ibin]->SetMarkerStyle(32); h_reco_dbar_inc_ratio[ibin]->Draw("e0 same");
    h_reco_c_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_c_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_c_inc_ratio[ibin]->SetMarkerStyle(28); h_reco_c_inc_ratio[ibin]->Draw("e0 same");
    h_reco_s_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_s_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_s_inc_ratio[ibin]->SetMarkerStyle(29); h_reco_s_inc_ratio[ibin]->Draw("e0 same");
    h_reco_b_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_b_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_b_inc_ratio[ibin]->SetMarkerStyle(30); h_reco_b_inc_ratio[ibin]->Draw("e0 same");
    tl2->Draw("same");

  }

  TCanvas *c_reco_gb = new TCanvas("c_reco_gb","",1500,600);
  c_reco_gb->Divide(5,2,0);
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0) c_reco_gb->cd(ibin+1);
    else c_reco_gb->cd(6-ibin);   
    gPad->SetLogy();
    h_reco_MC_g[ibin]->GetXaxis()->SetTitle("jet pT");
    h_reco_MC_g[ibin]->GetYaxis()->SetNdivisions(505);
    h_reco_MC_g[ibin]->GetYaxis()->SetTitleSize(0.05);
    h_reco_MC_g[ibin]->GetYaxis()->SetLabelSize(0.07);
    h_reco_MC_g[ibin]->GetXaxis()->CenterTitle();
    h_reco_MC_g[ibin]->GetXaxis()->SetNdivisions(505);
    h_reco_MC_g[ibin]->GetXaxis()->SetTitleSize(0.05);
    h_reco_MC_g[ibin]->GetXaxis()->SetLabelSize(0.07);
    h_reco_MC_g[ibin]->GetXaxis()->SetRangeUser(120.,400.);
    h_reco_MC_g[ibin]->SetLineColor(kRed); h_reco_MC_g[ibin]->SetMarkerColor(kRed);
    h_reco_MC_g[ibin]->SetMarkerStyle(4); h_reco_MC_g[ibin]->Draw("e0 same");
    h_reco_MC_upbar[ibin]->SetLineColor(kBlue); h_reco_MC_upbar[ibin]->SetMarkerColor(kBlue);
    h_reco_MC_upbar[ibin]->SetMarkerStyle(26); h_reco_MC_upbar[ibin]->Draw("e0 same");
    h_reco_MC_dbar[ibin]->SetLineColor(kBlue); h_reco_MC_dbar[ibin]->SetMarkerColor(kBlue);
    h_reco_MC_dbar[ibin]->SetMarkerStyle(32); h_reco_MC_dbar[ibin]->Draw("e0 same");
    h_reco_MC_c[ibin]->Draw("e0 same");
    h_reco_MC_s[ibin]->Draw("e0 same");
    h_reco_MC_b[ibin]->Draw("e0 same");

    if(ibin>0){
      TString tmp1 = cent_tag[ibin-1];
      tx1 = new TLatex(); tx1->SetTextSize(.12);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }

    if(ibin==0) c_reco_gb->cd(ibin+6);
    else c_reco_gb->cd(11-ibin);   
    h_reco_gluon_upbar_ratio[ibin] = (TH1D*)h_reco_MC_upbar[ibin]->Clone("h_reco_gluon_upbar_ratio");
    h_reco_gluon_upbar_ratio[ibin]->Divide(h_reco_MC_g[ibin]);
    h_reco_gluon_dbar_ratio[ibin] = (TH1D*)h_reco_MC_dbar[ibin]->Clone("h_reco_gluon_dbar_ratio");
    h_reco_gluon_dbar_ratio[ibin]->Divide(h_reco_MC_g[ibin]);
    h_reco_gluon_c_ratio[ibin] = (TH1D*)h_reco_MC_c[ibin]->Clone("h_reco_gluon_c_ratio");
    h_reco_gluon_c_ratio[ibin]->Divide(h_reco_MC_g[ibin]);
    h_reco_gluon_s_ratio[ibin] = (TH1D*)h_reco_MC_s[ibin]->Clone("h_reco_gluon_s_ratio");
    h_reco_gluon_s_ratio[ibin]->Divide(h_reco_MC_g[ibin]);
    h_reco_gluon_b_ratio[ibin] = (TH1D*)h_reco_MC_b[ibin]->Clone("h_reco_gluon_b_ratio");
    h_reco_gluon_b_ratio[ibin]->Divide(h_reco_MC_g[ibin]);
    h_reco_gluon_dbar_ratio[ibin]->GetYaxis()->SetNdivisions(505);
    h_reco_gluon_dbar_ratio[ibin]->GetYaxis()->SetTitleSize(0.05);
    h_reco_gluon_dbar_ratio[ibin]->GetYaxis()->SetLabelSize(0.07);
    h_reco_gluon_dbar_ratio[ibin]->GetXaxis()->CenterTitle();
    h_reco_gluon_dbar_ratio[ibin]->GetXaxis()->SetNdivisions(505);
    h_reco_gluon_dbar_ratio[ibin]->GetXaxis()->SetTitleSize(0.05);
    h_reco_gluon_dbar_ratio[ibin]->GetXaxis()->SetLabelSize(0.07);
    h_reco_gluon_dbar_ratio[ibin]->GetXaxis()->SetRangeUser(120.,400.);
    h_reco_gluon_dbar_ratio[ibin]->GetYaxis()->SetRangeUser(0.,0.15);
    h_reco_gluon_dbar_ratio[ibin]->SetLineColor(kBlue); h_reco_gluon_dbar_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_gluon_dbar_ratio[ibin]->SetMarkerStyle(23); h_reco_gluon_dbar_ratio[ibin]->Draw("e0 same");
    h_reco_gluon_upbar_ratio[ibin]->SetLineColor(kBlue); h_reco_gluon_upbar_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_gluon_upbar_ratio[ibin]->SetMarkerStyle(22); h_reco_gluon_upbar_ratio[ibin]->Draw("e0 same");
    h_reco_gluon_c_ratio[ibin]->SetLineColor(kBlue); h_reco_gluon_c_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_gluon_c_ratio[ibin]->SetMarkerStyle(28); h_reco_gluon_c_ratio[ibin]->Draw("e0 same");
    h_reco_gluon_s_ratio[ibin]->SetLineColor(kBlue); h_reco_gluon_s_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_gluon_s_ratio[ibin]->SetMarkerStyle(29); h_reco_gluon_s_ratio[ibin]->Draw("e0 same");
    h_reco_gluon_b_ratio[ibin]->SetLineColor(kBlue); h_reco_gluon_b_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_gluon_b_ratio[ibin]->SetMarkerStyle(30); h_reco_gluon_b_ratio[ibin]->Draw("e0 same");

    tl2->Draw("same");

  }

  char g_frac[50],g_error[50],chi2ndf[50], up_frac[50], down_frac[50];
  TString tmp;

  TCanvas *c_chg_udg_MC_gen = new TCanvas("c_chg_udg_MC_gen","",1500,600);
  c_chg_udg_MC_gen->Divide(5,2,0);
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins;ibin++){
    for(int ibin2=1;ibin2<ntrkbins;ibin2++){
    if(ibin==0) c_chg_udg_MC_gen->cd(1);
    else c_chg_udg_MC_gen->cd(6-ibin);   
    h_chg_ref[ibin][ibin2]->GetXaxis()->SetTitle("ref jet chg");
    h_chg_ref[ibin][ibin2]->GetYaxis()->SetTitle("# of jets");
    h_chg_ref[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
    h_chg_ref[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
    h_chg_ref[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
    h_chg_ref[ibin][ibin2]->GetXaxis()->CenterTitle();
    h_chg_ref[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
    h_chg_ref[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
    h_chg_ref[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
    h_chg_ref[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
    //h_chg_ref[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);
    h_chg_ref[ibin][ibin2]->SetLineColor(kBlack); h_chg_ref[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_ref[ibin][ibin2]->SetMarkerStyle(2);
    //h_chg_ref_q[ibin][ibin2]->SetMarkerStyle(22); h_chg_ref_q[ibin][ibin2]->SetLineColor(kBlue); h_chg_ref_q[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_ref_up[ibin][ibin2]->SetMarkerStyle(22); h_chg_ref_up[ibin][ibin2]->SetLineColor(kBlue); h_chg_ref_up[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_ref_down[ibin][ibin2]->SetMarkerStyle(23); h_chg_ref_down[ibin][ibin2]->SetLineColor(kBlue); h_chg_ref_down[ibin][ibin2]->SetMarkerColor(kBlue);
    //h_chg_ref_uubar[ibin][ibin2]->SetMarkerStyle(22); h_chg_ref_uubar[ibin][ibin2]->SetLineColor(kBlue); h_chg_ref_uubar[ibin][ibin2]->SetMarkerColor(kBlue);
    //h_chg_ref_ddbar[ibin][ibin2]->SetMarkerStyle(23); h_chg_ref_ddbar[ibin][ibin2]->SetLineColor(kBlue); h_chg_ref_ddbar[ibin][ibin2]->SetMarkerColor(kBlue);
    h_chg_ref_gubdb[ibin][ibin2]->SetLineColor(kRed); h_chg_ref_gubdb[ibin][ibin2]->SetMarkerColor(kRed); h_chg_ref[ibin][ibin2]->SetMarkerStyle(2);
    //h_chg_ref_g[ibin][ibin2]->SetLineColor(kRed); h_chg_ref_g[ibin][ibin2]->SetMarkerColor(kRed); h_chg_ref_g[ibin][ibin2]->SetMarkerStyle(2);
    //if(ibin2==3){
      h_chg_ref[ibin][ibin2]->Draw("e0 same");
      h_chg_ref_up[ibin][ibin2]->Draw("e0 same");
      h_chg_ref_down[ibin][ibin2]->Draw("e0 same");
      //h_chg_ref_uubar[ibin][ibin2]->Rebin(10); h_chg_ref_uubar[ibin][ibin2]->Draw("e0 same");
      //h_chg_ref_ddbar[ibin][ibin2]->Rebin(10); h_chg_ref_ddbar[ibin][ibin2]->Draw("e0 same");
      h_chg_ref_gubdb[ibin][ibin2]->Draw("e0 same");

      if(ibin>0){
        TString tmp1 = cent_tag[ibin-1];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.15,0.85, tmp1);
      }
    
      TString tmp = trk_tag[ibin2];
      tx = new TLatex(); tx->SetTextSize(.1);
      tx->DrawLatexNDC(0.15,0.65, tmp);
    
    if(do_fit==true){
    
      //f_chg_udg_ref[ibin][ibin2]->SetParLimits(0,0.,1.);
      //f_chg_udg_ref[ibin][ibin2]->SetParLimits(1,0.,1.);
      h_chg_ref[ibin][ibin2]->Fit(f_chg_udg_ref[ibin][ibin2],"Q N M R L","sames",-2.,2.);
    
      if(ibin==0) c_chg_udg_MC_gen->cd(6);
      else c_chg_udg_MC_gen->cd(11-ibin);   
      h_chg_ref_ratio[ibin][ibin2] = (TH1D*)h_chg_ref[ibin][ibin2]->Clone("h_chg_ref_ratio");
      h_chg_ref_ratio[ibin][ibin2]->Divide(f_chg_udg_ref[ibin][ibin2]);
      h_chg_ref_ratio[ibin][ibin2]->Rebin(5); h_chg_ref_ratio[ibin][ibin2]->Scale(1./5.);
      h_chg_ref_ratio[ibin][ibin2]->GetXaxis()->SetTitle("ref jet charge");
      h_chg_ref_ratio[ibin][ibin2]->GetYaxis()->SetTitle("ref / Fit");
      h_chg_ref_ratio[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
      h_chg_ref_ratio[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
      h_chg_ref_ratio[ibin][ibin2]->GetYaxis()->SetLabelSize(0.09);
      h_chg_ref_ratio[ibin][ibin2]->GetXaxis()->CenterTitle();
      h_chg_ref_ratio[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
      h_chg_ref_ratio[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
      h_chg_ref_ratio[ibin][ibin2]->GetXaxis()->SetLabelSize(0.09);
      h_chg_ref_ratio[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
      h_chg_ref_ratio[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,2.);
      h_chg_ref_ratio[ibin][ibin2]->SetLineColor(kBlack); h_chg_ref_ratio[ibin][ibin2]->SetMarkerColor(kBlack);
      h_chg_ref_ratio[ibin][ibin2]->SetMarkerStyle(1); h_chg_ref_ratio[ibin][ibin2]->Draw("e0 same");
      tl1->Draw("same");
      sprintf(g_frac, "g (MC/Fit) %.3f",(1.-(f_chg_udg_ref[ibin][ibin2]->GetParameter(0)+f_chg_udg_ref[ibin][ibin2]->GetParameter(1))));
      sprintf(g_error, "#pm %.3f",f_chg_udg_ref[ibin][ibin2]->GetParError(0));
      sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_udg_ref[ibin][ibin2]->GetChisquare()/f_chg_udg_ref[ibin][ibin2]->GetNDF());
      sprintf(up_frac, "up (MC/Fit) %.3f",(f_chg_udg_ref[ibin][ibin2]->GetParameter(0)));
      sprintf(down_frac, "down (MC/Fit) %.3f",(f_chg_udg_ref[ibin][ibin2]->GetParameter(1)));
      tx1 = new TLatex(); tx1->SetTextSize(.07);
      tx1->DrawLatexNDC(0.15,0.75, g_frac);
      tx3 = new TLatex(); tx3->SetTextSize(.06);
      tx3->DrawLatexNDC(0.65,0.7, g_error);
      tx4 = new TLatex(); tx4->SetTextSize(.07);
      tx4->DrawLatexNDC(0.3,0.4, up_frac);
      //tx5 = new TLatex(); tx5->SetTextSize(.07);
      //tx5->DrawLatexNDC(0.3,0.3, down_frac);
      tx2 = new TLatex(); tx2->SetTextSize(.07);
      tx2->DrawLatexNDC(0.3,0.2, chi2ndf);
    }}
  }

  h_dummy->GetXaxis()->SetTitle("p_{T} cut");
  h_dummy->GetYaxis()->SetTitle("1/N_{jets}");
  h_dummy->GetYaxis()->SetNdivisions(505);
  h_dummy->GetYaxis()->SetTitleSize(0.05);
  h_dummy->GetYaxis()->SetTitleOffset(0.75);
  h_dummy->GetYaxis()->SetLabelSize(0.07);
  h_dummy->GetXaxis()->CenterTitle();
  h_dummy->GetXaxis()->SetNdivisions(505);
  h_dummy->GetXaxis()->SetTitleSize(0.05);
  h_dummy->GetXaxis()->SetLabelSize(0.07);
  h_dummy->GetXaxis()->SetRangeUser(1.,6.);
  h_dummy->GetYaxis()->SetRangeUser(0.15,0.85);

  h_dummy2->GetXaxis()->SetTitle("p_{T} cut");
  h_dummy2->GetYaxis()->SetTitle("1/N_{jets}");
  h_dummy2->GetYaxis()->SetNdivisions(505);
  h_dummy2->GetYaxis()->SetTitleSize(0.05);
  h_dummy2->GetYaxis()->SetTitleOffset(0.75);
  h_dummy2->GetYaxis()->SetLabelSize(0.07);
  h_dummy2->GetXaxis()->CenterTitle();
  h_dummy2->GetXaxis()->SetNdivisions(505);
  h_dummy2->GetXaxis()->SetTitleSize(0.05);
  h_dummy2->GetXaxis()->SetLabelSize(0.07);
  h_dummy2->GetXaxis()->SetRangeUser(1.,6.);
  h_dummy2->GetYaxis()->SetRangeUser(0.,0.4);

  h_dummy3->GetXaxis()->SetTitle("jet p_{T}");
  h_dummy3->GetYaxis()->SetTitle("mean");
  h_dummy3->GetYaxis()->SetNdivisions(505);
  h_dummy3->GetYaxis()->SetTitleSize(0.05);
  h_dummy3->GetYaxis()->SetTitleOffset(0.75);
  h_dummy3->GetYaxis()->SetLabelSize(0.07);
  h_dummy3->GetXaxis()->CenterTitle();
  h_dummy3->GetXaxis()->SetNdivisions(505);
  h_dummy3->GetXaxis()->SetTitleSize(0.05);
  h_dummy3->GetXaxis()->SetLabelSize(0.07);
  h_dummy3->GetXaxis()->SetRangeUser(120.,200.);
  h_dummy3->GetYaxis()->SetRangeUser(-0.3,0.4);

  h_dummy4->GetXaxis()->SetTitle("jet p_{T}");
  h_dummy4->GetYaxis()->SetTitle("sigma");
  h_dummy4->GetYaxis()->SetNdivisions(505);
  h_dummy4->GetYaxis()->SetTitleSize(0.05);
  h_dummy4->GetYaxis()->SetTitleOffset(0.75);
  h_dummy4->GetYaxis()->SetLabelSize(0.07);
  h_dummy4->GetXaxis()->CenterTitle();
  h_dummy4->GetXaxis()->SetNdivisions(505);
  h_dummy4->GetXaxis()->SetTitleSize(0.05);
  h_dummy4->GetXaxis()->SetLabelSize(0.07);
  h_dummy4->GetXaxis()->SetRangeUser(120.,200.);
  h_dummy4->GetYaxis()->SetRangeUser(0.35,0.65);

  Double_t up_fit[ntrkbins-1]; Double_t up_fiterr[ntrkbins-1];
  Double_t down_fit[ntrkbins-1]; Double_t down_fiterr[ntrkbins-1];
  Double_t gluon_fit[ntrkbins-1]; Double_t gluon_fiterr[ntrkbins-1];
  Double_t quark_fit[ntrkbins-1]; Double_t quark_fiterr[ntrkbins-1];

  for(int ibin=0; ibin<nCBins; ibin++){
    for(int ibin2=0;ibin2<ntrkbins-1;ibin2++){
      up_fit[ibin2] = f_chg_udg_ref[ibin][ibin2+1]->GetParameter(0); up_fiterr[ibin2] = 0;// f_chg_udg_ref[ibin][ibin2+1]->GetParError(0);
      down_fit[ibin2] = f_chg_udg_ref[ibin][ibin2+1]->GetParameter(1); down_fiterr[ibin2] = 0;// f_chg_udg_ref[ibin][ibin2+1]->GetParError(0);
      gluon_fit[ibin2] = 1.-(f_chg_udg_ref[ibin][ibin2+1]->GetParameter(0)+f_chg_udg_ref[ibin][ibin2+1]->GetParameter(1)); gluon_fiterr[ibin2] = 0;// f_chg_udg_ref[ibin][ibin2+1]->GetParError(0);
    }
    up_ref[ibin] = new TGraphErrors(ntrkbins-1,x_mean,up_fit,x_mean_err,up_fiterr);
    up_ref[ibin]->SetName((TString)("up_ref_cent"+cent[ibin]));
    down_ref[ibin] = new TGraphErrors(ntrkbins-1,x_mean,down_fit,x_mean_err,down_fiterr);
    down_ref[ibin]->SetName((TString)("down_ref_cent"+cent[ibin]));
    gluon_ref[ibin] = new TGraphErrors(ntrkbins-1,x_mean,gluon_fit,x_mean_err,gluon_fiterr);
    gluon_ref[ibin]->SetName((TString)("gluon_ref_cent"+cent[ibin]));
  }

  TCanvas *c_ref_graph = new TCanvas("c_ref_graph","c_ref_graph",1500,350);
  c_ref_graph->Divide(5,1,0);
  
  for(int ibin=0; ibin<nCBins; ibin++){
    if(ibin==0) c_ref_graph->cd(1);
    else c_ref_graph->cd(6-ibin);
    h_dummy->Draw("same");
    up_ref[ibin]->SetMarkerColor(kBlue); up_ref[ibin]->SetMarkerStyle(20); up_ref[ibin]->SetMarkerSize(1.2);
    up_ref[ibin]->Draw("PL");
    down_ref[ibin]->SetMarkerColor(kBlue); down_ref[ibin]->SetMarkerStyle(20); down_ref[ibin]->SetMarkerSize(1.2);
    down_ref[ibin]->Draw("PL");
    gluon_ref[ibin]->SetMarkerColor(kRed); gluon_ref[ibin]->SetMarkerStyle(20); gluon_ref[ibin]->SetMarkerSize(1.2);
    gluon_ref[ibin]->Draw("PL");
  
    //tl_q_ref[ibin]->Draw("same"); tl_g_ref[ibin]->Draw("same");
    tl_g_ref[ibin]->Draw("same"); tl_u_ref[ibin]->Draw("same"); tl_d_ref[ibin]->Draw("same");

    if(ibin>0){
        TString tmp1 = cent_tag[ibin-1];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.55,0.85, tmp1);
    }
  }

  TCanvas *c_chg_udg_MC_reco[ntrkbins];
  for(int ibin2=1;ibin2<ntrkbins;ibin2++){
    c_chg_udg_MC_reco[ibin2] = new TCanvas((TString)("c_chg_udg_MC_reco"+trk[ibin2]),(TString)("c_chg_udg_MC_reco"+trk[ibin2]),1500,600);
    c_chg_udg_MC_reco[ibin2]->Divide(5,2,0);
    gStyle->SetOptStat(0);
    for(int ibin=0;ibin<nCBins;ibin++){
      if(ibin==0) c_chg_udg_MC_reco[ibin2]->cd(1);
      else c_chg_udg_MC_reco[ibin2]->cd(6-ibin);   
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
      h_chg_reco[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_reco[ibin][ibin2]->SetMarkerStyle(2);
      h_chg_reco_up[ibin][ibin2]->SetMarkerStyle(22); h_chg_reco_up[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_up[ibin][ibin2]->SetMarkerColor(kBlue);
      h_chg_reco_down[ibin][ibin2]->SetMarkerStyle(23); h_chg_reco_down[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_down[ibin][ibin2]->SetMarkerColor(kBlue);
      h_chg_reco_gubdb[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_gubdb[ibin][ibin2]->SetMarkerColor(kRed); h_chg_reco_gubdb[ibin][ibin2]->SetMarkerStyle(2);
    
      h_chg_reco[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_up[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_down[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_gubdb[ibin][ibin2]->Draw("e0 same");

      if(ibin>0){
        TString tmp1 = cent_tag[ibin-1];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.15,0.85, tmp1);
      }
    
      TString tmp = trk_tag[ibin2];
      tx = new TLatex(); tx->SetTextSize(.1);
      tx->DrawLatexNDC(0.15,0.65, tmp);
    
      if(do_fit==true){
      h_chg_reco[ibin][ibin2]->Fit(f_chg_udg_data[ibin][ibin2],"N M R","sames",-2.,2.);
    
      if(ibin==0) c_chg_udg_MC_reco[ibin2]->cd(6);
      else c_chg_udg_MC_reco[ibin2]->cd(11-ibin);   
      h_chg_reco_ratio[ibin][ibin2] = (TH1D*)h_chg_reco[ibin][ibin2]->Clone("h_chg_reco_ratio");
      h_chg_reco_ratio[ibin][ibin2]->Divide(f_chg_udg_data[ibin][ibin2]);
      h_chg_reco_ratio[ibin][ibin2]->Rebin(5); h_chg_reco_ratio[ibin][ibin2]->Scale(1./5.);
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
      sprintf(g_frac, "g+ubar+dbar (MC/Fit) %.3f",(1.-(f_chg_udg_data[ibin][ibin2]->GetParameter(0)+f_chg_udg_data[ibin][ibin2]->GetParameter(1)))/gluon_ubar_dbar[ibin]);
      sprintf(g_error, "#pm %.3f",f_chg_udg_data[ibin][ibin2]->GetParError(0));
      sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_udg_data[ibin][ibin2]->GetChisquare()/f_chg_udg_data[ibin][ibin2]->GetNDF());
      sprintf(up_frac, "up (MC/Fit) %.3f",(f_chg_udg_data[ibin][ibin2]->GetParameter(0)/up[ibin]));
      sprintf(down_frac, "down (MC/Fit) %.3f",(f_chg_udg_data[ibin][ibin2]->GetParameter(1)/down[ibin]));
      tx1 = new TLatex(); tx1->SetTextSize(.07);
      tx1->DrawLatexNDC(0.15,0.75, g_frac);
      tx3 = new TLatex(); tx3->SetTextSize(.06);
      tx3->DrawLatexNDC(0.65,0.7, g_error);
      tx4 = new TLatex(); tx4->SetTextSize(.07);
      tx4->DrawLatexNDC(0.3,0.4, up_frac);
      tx5 = new TLatex(); tx5->SetTextSize(.07);
      tx5->DrawLatexNDC(0.3,0.3, down_frac);
      tx2 = new TLatex(); tx2->SetTextSize(.07);
      tx2->DrawLatexNDC(0.3,0.2, chi2ndf);
    }}
  }

  for(int ibin=0; ibin<nCBins; ibin++){
    for(int ibin2=0;ibin2<ntrkbins-1;ibin2++){
      up_fit[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParameter(0); up_fiterr[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParError(0);
      down_fit[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParameter(1); down_fiterr[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParError(1);
      gluon_fit[ibin2] = 1.-(f_chg_udg_data[ibin][ibin2+1]->GetParameter(0)+f_chg_udg_data[ibin][ibin2+1]->GetParameter(1)+ubdbcsb[ibin][ibin2+1]); gluon_fiterr[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParError(1);
      quark_fit[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParameter(0)+f_chg_udg_data[ibin][ibin2+1]->GetParameter(1)+ubdbcsb[ibin][ibin2+1]; quark_fiterr[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParError(1);

    }
    up_reco[ibin] = new TGraphErrors(ntrkbins-1,x_mean,up_fit,x_mean_err,up_fiterr);
    up_reco[ibin]->SetName((TString)("up_reco_cent"+cent[ibin]));
    down_reco[ibin] = new TGraphErrors(ntrkbins-1,x_mean,down_fit,x_mean_err,down_fiterr);
    down_reco[ibin]->SetName((TString)("down_reco_cent"+cent[ibin]));
    gluon_reco[ibin] = new TGraphErrors(ntrkbins-1,x_mean,gluon_fit,x_mean_err,gluon_fiterr);
    gluon_reco[ibin]->SetName((TString)("gluon_reco_cent"+cent[ibin]));
    quark_reco[ibin] = new TGraphErrors(ntrkbins-1,x_mean,quark_fit,x_mean_err,quark_fiterr);
    quark_reco[ibin]->SetName((TString)("quark_reco_cent"+cent[ibin]));
  }


  TCanvas *c_reco_graph = new TCanvas("c_reco_graph","c_reco_graph",1500,350);
  c_reco_graph->Divide(5,1,0);
  
  for(int ibin=0; ibin<nCBins; ibin++){
    if(ibin==0) c_reco_graph->cd(1);
    else c_reco_graph->cd(6-ibin);
    h_dummy->Draw("same");
    up_reco[ibin]->SetMarkerColor(kBlue); up_reco[ibin]->SetMarkerStyle(26); up_reco[ibin]->SetMarkerSize(1.3);
    up_reco[ibin]->Draw("PL");
    down_reco[ibin]->SetMarkerColor(kBlue); down_reco[ibin]->SetMarkerStyle(32); down_reco[ibin]->SetMarkerSize(1.3);
    down_reco[ibin]->Draw("PL");
    gluon_reco[ibin]->SetMarkerColor(kRed); gluon_reco[ibin]->SetMarkerStyle(20); gluon_reco[ibin]->SetMarkerSize(1.2);
    gluon_reco[ibin]->Draw("PL");
    quark_reco[ibin]->SetMarkerColor(kBlue); quark_reco[ibin]->SetMarkerStyle(20); quark_reco[ibin]->SetMarkerSize(1.2);
    quark_reco[ibin]->Draw("PL");
  
    tl_u[ibin]->Draw("same"); tl_d[ibin]->Draw("same"); tl_g[ibin]->Draw("same"); tl_q[ibin]->Draw("same");

    if(ibin>0){
        TString tmp1 = cent_tag[ibin-1];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.55,0.85, tmp1);
    }
  }

  TCanvas *c_chg_udg_MC_data[ntrkbins];

  for(int ibin2=1;ibin2<ntrkbins;ibin2++){
    c_chg_udg_MC_data[ibin2] = new TCanvas((TString)("c_chg_udg_MC_data"+trk[ibin2]),(TString)("c_chg_udg_MC_data"+trk[ibin2]),1500,600);
    c_chg_udg_MC_data[ibin2]->Divide(5,2,0);
    //gStyle->SetOptStat(0);
    //gStyle->SetOptFit(1);
    for(int ibin=0;ibin<nCBins;ibin++){
  
      if(ibin==0) c_chg_udg_MC_data[ibin2]->cd(1);
      else c_chg_udg_MC_data[ibin2]->cd(6-ibin);   
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
      //h_chg_data[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);
      h_chg_data[ibin][ibin2]->SetLineColor(kGreen); h_chg_data[ibin][ibin2]->SetMarkerColor(kGreen); h_chg_data[ibin][ibin2]->SetMarkerStyle(2);
    
      h_chg_data[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_up[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_down[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_gubdb[ibin][ibin2]->Draw("e0 same");
    
      if(ibin>0){
        TString tmp1 = cent_tag[ibin-1];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.15,0.85, tmp1);
      }
    
      TString tmp = trk_tag[ibin2];
      tx = new TLatex(); tx->SetTextSize(.1);
      tx->DrawLatexNDC(0.15,0.65, tmp);

      if(do_fit==true){  
        h_chg_data[ibin][ibin2]->Fit(f_chg_udg_data[ibin][ibin2],"Q N M R","sames",-2.,2.);
    
        if(ibin==0) c_chg_udg_MC_data[ibin2]->cd(6);
        else c_chg_udg_MC_data[ibin2]->cd(11-ibin);   
      h_chg_reco_ratio[ibin][ibin2] = (TH1D*)h_chg_data[ibin][ibin2]->Clone("h_chg_data_ratio");
      h_chg_reco_ratio[ibin][ibin2]->Divide(f_chg_udg_data[ibin][ibin2]);
      h_chg_reco_ratio[ibin][ibin2]->Rebin(5); h_chg_reco_ratio[ibin][ibin2]->Scale(1./5.);
      h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetTitle("data jet charge");
      h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetTitle("Data / Fit");
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
      sprintf(g_frac, "g+ubar+dbar %.3f",(1.-(f_chg_udg_data[ibin][ibin2]->GetParameter(0)+f_chg_udg_data[ibin][ibin2]->GetParameter(1))));
      sprintf(g_error, "#pm %.3f",f_chg_udg_data[ibin][ibin2]->GetParError(0));
      sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_udg_data[ibin][ibin2]->GetChisquare()/f_chg_udg_data[ibin][ibin2]->GetNDF());
      sprintf(up_frac, "up %.3f",f_chg_udg_data[ibin][ibin2]->GetParameter(0));
      sprintf(down_frac, "down %.3f",f_chg_udg_data[ibin][ibin2]->GetParameter(1));
      tx1 = new TLatex(); tx1->SetTextSize(.08);
      tx1->DrawLatexNDC(0.25,0.85, g_frac);
      tx3 = new TLatex(); tx3->SetTextSize(.06);
      tx3->DrawLatexNDC(0.65,0.8, g_error);
      tx4 = new TLatex(); tx4->SetTextSize(.07);
      tx4->DrawLatexNDC(0.3,0.4, up_frac);
      tx5 = new TLatex(); tx5->SetTextSize(.07);
      tx5->DrawLatexNDC(0.3,0.3, down_frac);      
      tx2 = new TLatex(); tx2->SetTextSize(.07);
      tx2->DrawLatexNDC(0.3,0.2, chi2ndf);
    }}
  }

  for(int ibin=0; ibin<nCBins; ibin++){
    for(int ibin2=0;ibin2<ntrkbins-1;ibin2++){
      up_fit[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParameter(0); up_fiterr[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParError(0);
      down_fit[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParameter(1); down_fiterr[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParError(1);
      gluon_fit[ibin2] = 1.-(f_chg_udg_data[ibin][ibin2+1]->GetParameter(0)+f_chg_udg_data[ibin][ibin2+1]->GetParameter(1)+ubdbcsb[ibin][ibin2+1]); gluon_fiterr[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParError(1);
      quark_fit[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParameter(0)+f_chg_udg_data[ibin][ibin2+1]->GetParameter(1)+ubdbcsb[ibin][ibin2+1]; quark_fiterr[ibin2] = f_chg_udg_data[ibin][ibin2+1]->GetParError(1);

      cout<<"quark["<<ibin<<"]["<<ibin2<<"] = "<<quark_fit[ibin2]<<";"<<endl;
      cout<<"gluon["<<ibin<<"]["<<ibin2<<"] = "<<gluon_fit[ibin2]<<";"<<endl;
      //cout<<"cent: "<<ibin<<".  trk: "<<ibin2<<".   gluon: "<<gluon_fit[ibin2]<<".  quark:"<<quark_fit[ibin2]<<endl;
    }
    up_data[ibin] = new TGraphErrors(ntrkbins-1,x_mean,up_fit,x_mean_err,up_fiterr);
    up_data[ibin]->SetName((TString)("up_data_cent"+cent[ibin]));
    down_data[ibin] = new TGraphErrors(ntrkbins-1,x_mean,down_fit,x_mean_err,down_fiterr);
    down_data[ibin]->SetName((TString)("down_data_cent"+cent[ibin]));
    gluon_data[ibin] = new TGraphErrors(ntrkbins-1,x_mean,gluon_fit,x_mean_err,gluon_fiterr);
    gluon_data[ibin]->SetName((TString)("gluon_data_cent"+cent[ibin]));
    quark_data[ibin] = new TGraphErrors(ntrkbins-1,x_mean,quark_fit,x_mean_err,quark_fiterr);
    quark_data[ibin]->SetName((TString)("quark_data_cent"+cent[ibin]));
  }

  TCanvas *c_data_graph = new TCanvas("c_data_graph","c_data_graph",1500,350);
  c_data_graph->Divide(5,1,0);
  
  for(int ibin=0; ibin<nCBins; ibin++){
    if(ibin==0) c_data_graph->cd(1);
    else c_data_graph->cd(6-ibin);
    h_dummy->Draw("same");
    up_data[ibin]->SetTitle(" ; pT cut; fraction; ");
    up_data[ibin]->GetXaxis()->SetTitle(" ; pT cut; fraction; ");
    up_data[ibin]->SetMarkerColor(kBlue); up_data[ibin]->SetMarkerStyle(26); up_data[ibin]->SetMarkerSize(1.3);
    up_data[ibin]->Draw("PL");
    down_data[ibin]->SetMarkerColor(kBlue); down_data[ibin]->SetMarkerStyle(32); down_data[ibin]->SetMarkerSize(1.3);
    down_data[ibin]->Draw("PL");
    gluon_data[ibin]->SetMarkerColor(kRed); gluon_data[ibin]->SetMarkerStyle(20); gluon_data[ibin]->SetMarkerSize(1.2);
    gluon_data[ibin]->Draw("PL");
    quark_data[ibin]->SetMarkerColor(kBlue); quark_data[ibin]->SetMarkerStyle(20); quark_data[ibin]->SetMarkerSize(1.2);
    quark_data[ibin]->Draw("PL");

    if(ibin==0) {tl_u[ibin]->Draw("same"); tl_d[ibin]->Draw("same"); tl_g[ibin]->Draw("same"); tl_q[ibin]->Draw("same");}
    else {tl_u_data->Draw("same"); tl_d_data->Draw("same"); tl_g[ibin]->Draw("same"); tl_q[ibin]->Draw("same");}
    
    if(ibin>0){
        TString tmp1 = cent_tag[ibin-1];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }
  }

  TCanvas *c_data_qg_graph = new TCanvas("c_data_qg_graph","c_data_qg_graph",1500,350);
  c_data_qg_graph->Divide(5,1,0);
  
  for(int ibin=0; ibin<nCBins; ibin++){
    if(ibin==0) c_data_qg_graph->cd(1);
    else c_data_qg_graph->cd(6-ibin);
    h_dummy->Draw("same");
    gluon_data[ibin]->SetMarkerColor(kRed); gluon_data[ibin]->SetMarkerStyle(20); gluon_data[ibin]->SetMarkerSize(1.2);
    gluon_data[ibin]->Draw("PL");
    quark_data[ibin]->SetMarkerColor(kBlue); quark_data[ibin]->SetMarkerStyle(20); quark_data[ibin]->SetMarkerSize(1.2);
    quark_data[ibin]->Draw("PL");

    tl_g[ibin]->Draw("same"); tl_q[ibin]->Draw("same");
    
    if(ibin>0){
        TString tmp1 = cent_tag[ibin-1];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }
  }

  TCanvas *c_data_q_graph = new TCanvas("c_data_q_graph","c_data_q_graph",1500,350);
  c_data_q_graph->Divide(5,1,0);
  
  for(int ibin=0; ibin<nCBins; ibin++){
    if(ibin==0) c_data_q_graph->cd(1);
    else c_data_q_graph->cd(6-ibin);
    h_dummy2->Draw("same");
    up_data[ibin]->SetMarkerColor(kBlue); up_data[ibin]->SetMarkerStyle(22); up_data[ibin]->SetMarkerSize(1.3);
    up_data[ibin]->Draw("PL");
    down_data[ibin]->SetMarkerColor(kBlue); down_data[ibin]->SetMarkerStyle(23); down_data[ibin]->SetMarkerSize(1.3);
    down_data[ibin]->Draw("PL");

    if(ibin==0) {tl_u[ibin]->Draw("same"); tl_d[ibin]->Draw("same"); tl_ubar[ibin]->Draw("same"); tl_dbar[ibin]->Draw("same"); tl_c[ibin]->Draw("same"); tl_s[ibin]->Draw("same");tl_b[ibin]->Draw("same");}
    else {tl_u_data->Draw("same"); tl_d_data->Draw("same");}
    
    if(ibin>0){
        TString tmp1 = cent_tag[ibin-1];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }
  }






  TCanvas *c_chg_avg_MC[ntrkbins];

  for(int ibin2=1;ibin2<ntrkbins;ibin2++){
    c_chg_avg_MC[ibin2] = new TCanvas((TString)("c_chg_avg_MC"+trk[ibin2]),(TString)("c_chg_avg_MC"+trk[ibin2]),1500,600);
    c_chg_avg_MC[ibin2]->Divide(5,2,0);
    //gStyle->SetOptStat(0);
    //gStyle->SetOptFit(1);
    for(int ibin=0;ibin<1;ibin++){
     for(int ibin3=0;ibin3<njtptbins;ibin3++){  
      c_chg_avg_MC[ibin2]->cd(ibin3+1);   
      h_chg_reco_avg[ibin][ibin3][ibin2]->GetXaxis()->SetTitle("jet chg");
      h_chg_reco_avg[ibin][ibin3][ibin2]->GetYaxis()->SetTitle("# of jets");
      h_chg_reco_avg[ibin][ibin3][ibin2]->GetYaxis()->SetNdivisions(505);
      h_chg_reco_avg[ibin][ibin3][ibin2]->GetYaxis()->SetTitleSize(0.05);
      h_chg_reco_avg[ibin][ibin3][ibin2]->GetYaxis()->SetLabelSize(0.05);
      h_chg_reco_avg[ibin][ibin3][ibin2]->GetXaxis()->CenterTitle();
      h_chg_reco_avg[ibin][ibin3][ibin2]->GetXaxis()->SetNdivisions(505);
      h_chg_reco_avg[ibin][ibin3][ibin2]->GetXaxis()->SetTitleSize(0.05);
      h_chg_reco_avg[ibin][ibin3][ibin2]->GetXaxis()->SetLabelSize(0.05);
      h_chg_reco_avg[ibin][ibin3][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
      h_chg_reco_avg[ibin][ibin3][ibin2]->SetLineColor(kBlack); h_chg_reco_avg[ibin][ibin3][ibin2]->SetMarkerColor(kBlack); h_chg_reco_avg[ibin][ibin3][ibin2]->SetMarkerStyle(2);
      h_chg_reco_q_avg[ibin][ibin3][ibin2]->SetLineColor(kBlue); h_chg_reco_q_avg[ibin][ibin3][ibin2]->SetMarkerColor(kBlue); h_chg_reco_q_avg[ibin][ibin3][ibin2]->SetMarkerStyle(2);    
      h_chg_reco_g_avg[ibin][ibin3][ibin2]->SetLineColor(kRed); h_chg_reco_g_avg[ibin][ibin3][ibin2]->SetMarkerColor(kRed); h_chg_reco_g_avg[ibin][ibin3][ibin2]->SetMarkerStyle(2);
      h_chg_reco_up_avg[ibin][ibin3][ibin2]->SetLineColor(kBlue); h_chg_reco_up_avg[ibin][ibin3][ibin2]->SetMarkerColor(kBlue); h_chg_reco_up_avg[ibin][ibin3][ibin2]->SetMarkerStyle(22);
      h_chg_reco_down_avg[ibin][ibin3][ibin2]->SetLineColor(kBlue); h_chg_reco_down_avg[ibin][ibin3][ibin2]->SetMarkerColor(kBlue); h_chg_reco_down_avg[ibin][ibin3][ibin2]->SetMarkerStyle(23);

      h_chg_reco_avg[ibin][ibin3][ibin2]->Draw("e0 same");
      h_chg_reco_q_avg[ibin][ibin3][ibin2]->Draw("e0 same");
      h_chg_reco_g_avg[ibin][ibin3][ibin2]->Draw("e0 same");
      h_chg_reco_up_avg[ibin][ibin3][ibin2]->Draw("e0 same");
      h_chg_reco_down_avg[ibin][ibin3][ibin2]->Draw("e0 same");
      
      h_chg_reco_avg[ibin][ibin3][ibin2]->Fit(f_chg_reco_avg[ibin][ibin3][ibin2],"Q M R","sames",-2.,2.);
      h_chg_reco_q_avg[ibin][ibin3][ibin2]->Fit(f_chg_reco_q_avg[ibin][ibin3][ibin2],"Q M R","sames",-2.,2.);
      h_chg_reco_g_avg[ibin][ibin3][ibin2]->Fit(f_chg_reco_g_avg[ibin][ibin3][ibin2],"Q M R","sames",-2.,2.);
      h_chg_reco_up_avg[ibin][ibin3][ibin2]->Fit(f_chg_reco_up_avg[ibin][ibin3][ibin2],"Q M R","sames",-2.,2.);
      h_chg_reco_down_avg[ibin][ibin3][ibin2]->Fit(f_chg_reco_down_avg[ibin][ibin3][ibin2],"Q M R","sames",-2.,2.);
     }
   }
 }

  Double_t inclusive_mean[njtptbins]; Double_t quark_mean[njtptbins]; Double_t gluon_mean[njtptbins]; Double_t up_mean[njtptbins]; Double_t down_mean[njtptbins];
  Double_t inclusive_res[njtptbins]; Double_t quark_res[njtptbins]; Double_t gluon_res[njtptbins]; Double_t up_res[njtptbins]; Double_t down_res[njtptbins];

  TGraph *inclusive_avg[nCBins][ntrkbins]; TGraph *quark_avg[nCBins][ntrkbins]; TGraph *gluon_avg[nCBins][ntrkbins]; TGraph *up_avg[nCBins][ntrkbins]; TGraph *down_avg[nCBins][ntrkbins];
  TGraph *inclusive_sig[nCBins][ntrkbins]; TGraph *quark_sig[nCBins][ntrkbins]; TGraph *gluon_sig[nCBins][ntrkbins]; TGraph *up_sig[nCBins][ntrkbins]; TGraph *down_sig[nCBins][ntrkbins];

  for(int ibin=0; ibin<1; ibin++){
    for(int ibin2=0;ibin2<ntrkbins-1;ibin2++){
      for(int ibin3=0;ibin3<njtptbins;ibin3++){
        inclusive_mean[ibin3] = f_chg_reco_avg[ibin][ibin3][ibin2+1]->GetParameter(1);
        quark_mean[ibin3] = f_chg_reco_q_avg[ibin][ibin3][ibin2+1]->GetParameter(1);
        gluon_mean[ibin3] = f_chg_reco_g_avg[ibin][ibin3][ibin2+1]->GetParameter(1);
        up_mean[ibin3] = f_chg_reco_up_avg[ibin][ibin3][ibin2+1]->GetParameter(1);
        down_mean[ibin3] = f_chg_reco_down_avg[ibin][ibin3][ibin2+1]->GetParameter(1);

        inclusive_res[ibin3] = f_chg_reco_avg[ibin][ibin3][ibin2+1]->GetParameter(2);
        quark_res[ibin3] = f_chg_reco_q_avg[ibin][ibin3][ibin2+1]->GetParameter(2);
        gluon_res[ibin3] = f_chg_reco_g_avg[ibin][ibin3][ibin2+1]->GetParameter(2);
        up_res[ibin3] = f_chg_reco_up_avg[ibin][ibin3][ibin2+1]->GetParameter(2);
        down_res[ibin3] = f_chg_reco_down_avg[ibin][ibin3][ibin2+1]->GetParameter(2);
      }

    inclusive_avg[ibin][ibin2] = new TGraph(njtptbins-1,jtpt_bound,inclusive_mean);
    inclusive_avg[ibin][ibin2]->SetName((TString)("inclusive_mean_"+cent[ibin]+"_"+trk[ibin2]));
    quark_avg[ibin][ibin2] = new TGraph(njtptbins-1,jtpt_bound,quark_mean);
    quark_avg[ibin][ibin2]->SetName((TString)("quark_mean_"+cent[ibin]+"_"+trk[ibin2]));
    gluon_avg[ibin][ibin2] = new TGraph(njtptbins-1,jtpt_bound,gluon_mean);
    gluon_avg[ibin][ibin2]->SetName((TString)("gluon_mean_"+cent[ibin]+"_"+trk[ibin2]));
    up_avg[ibin][ibin2] = new TGraph(njtptbins-1,jtpt_bound,up_mean);
    up_avg[ibin][ibin2]->SetName((TString)("up_mean_"+cent[ibin]+"_"+trk[ibin2]));
    down_avg[ibin][ibin2] = new TGraph(njtptbins-1,jtpt_bound,down_mean);
    down_avg[ibin][ibin2]->SetName((TString)("down_mean_"+cent[ibin]+"_"+trk[ibin2]));

    inclusive_sig[ibin][ibin2] = new TGraph(njtptbins-1,jtpt_bound,inclusive_res);
    inclusive_sig[ibin][ibin2]->SetName((TString)("inclusive_res_"+cent[ibin]+"_"+trk[ibin2]));
    quark_sig[ibin][ibin2] = new TGraph(njtptbins-1,jtpt_bound,quark_res);
    quark_sig[ibin][ibin2]->SetName((TString)("quark_res_"+cent[ibin]+"_"+trk[ibin2]));
    gluon_sig[ibin][ibin2] = new TGraph(njtptbins-1,jtpt_bound,gluon_res);
    gluon_sig[ibin][ibin2]->SetName((TString)("gluon_res_"+cent[ibin]+"_"+trk[ibin2]));
    up_sig[ibin][ibin2] = new TGraph(njtptbins-1,jtpt_bound,up_res);
    up_sig[ibin][ibin2]->SetName((TString)("up_res_"+cent[ibin]+"_"+trk[ibin2]));
    down_sig[ibin][ibin2] = new TGraph(njtptbins-1,jtpt_bound,down_res);
    down_sig[ibin][ibin2]->SetName((TString)("down_res_"+cent[ibin]+"_"+trk[ibin2]));
    }
  }

  TCanvas *c_MC_avg_graph = new TCanvas("c_MC_avg_graph","c_MC_avg_graph",1500,500);
  c_MC_avg_graph->Divide(3,1,0);
  
  for(int ibin=0; ibin<1; ibin++){
   for(int ibin2=0; ibin2<ntrkbins-1; ibin2++){
    c_MC_avg_graph->cd(ibin2+1);
    h_dummy3->Draw("same");
    inclusive_avg[ibin][ibin2]->SetMarkerColor(kBlack); inclusive_avg[ibin][ibin2]->SetMarkerStyle(2); inclusive_avg[ibin][ibin2]->SetMarkerSize(1.3);
    inclusive_avg[ibin][ibin2]->Draw("PL");
    quark_avg[ibin][ibin2]->SetMarkerColor(kBlue); quark_avg[ibin][ibin2]->SetMarkerStyle(2); quark_avg[ibin][ibin2]->SetMarkerSize(1.3);
    quark_avg[ibin][ibin2]->Draw("PL");
    gluon_avg[ibin][ibin2]->SetMarkerColor(kRed); gluon_avg[ibin][ibin2]->SetMarkerStyle(2); gluon_avg[ibin][ibin2]->SetMarkerSize(1.3);
    gluon_avg[ibin][ibin2]->Draw("PL");
    up_avg[ibin][ibin2]->SetMarkerColor(kBlue); up_avg[ibin][ibin2]->SetMarkerStyle(22); up_avg[ibin][ibin2]->SetMarkerSize(1.3);
    up_avg[ibin][ibin2]->Draw("PL");
    down_avg[ibin][ibin2]->SetMarkerColor(kBlue); down_avg[ibin][ibin2]->SetMarkerStyle(23); down_avg[ibin][ibin2]->SetMarkerSize(1.3);
    down_avg[ibin][ibin2]->Draw("PL");
      
      TString tmp = trk_tag[ibin2+1];
      tx = new TLatex(); tx->SetTextSize(.1);
      tx->DrawLatexNDC(0.15,0.85, tmp);
   }
  }

  TCanvas *c_MC_res_graph = new TCanvas("c_MC_res_graph","c_MC_res_graph",1500,500);
  c_MC_res_graph->Divide(3,1,0);
  
  for(int ibin=0; ibin<1; ibin++){
   for(int ibin2=0; ibin2<ntrkbins-1; ibin2++){
    c_MC_res_graph->cd(ibin2+1);
    h_dummy4->Draw("same");
    inclusive_sig[ibin][ibin2]->SetMarkerColor(kBlack); inclusive_sig[ibin][ibin2]->SetMarkerStyle(2); inclusive_sig[ibin][ibin2]->SetMarkerSize(1.3);
    inclusive_sig[ibin][ibin2]->Draw("PL");
    quark_sig[ibin][ibin2]->SetMarkerColor(kBlue); quark_sig[ibin][ibin2]->SetMarkerStyle(2); quark_sig[ibin][ibin2]->SetMarkerSize(1.3);
    quark_sig[ibin][ibin2]->Draw("PL");
    gluon_sig[ibin][ibin2]->SetMarkerColor(kRed); gluon_sig[ibin][ibin2]->SetMarkerStyle(2); gluon_sig[ibin][ibin2]->SetMarkerSize(1.3);
    gluon_sig[ibin][ibin2]->Draw("PL");
    up_sig[ibin][ibin2]->SetMarkerColor(kBlue); up_sig[ibin][ibin2]->SetMarkerStyle(22); up_sig[ibin][ibin2]->SetMarkerSize(1.3);
    up_sig[ibin][ibin2]->Draw("PL");
    down_sig[ibin][ibin2]->SetMarkerColor(kBlue); down_sig[ibin][ibin2]->SetMarkerStyle(23); down_sig[ibin][ibin2]->SetMarkerSize(1.3);
    down_sig[ibin][ibin2]->Draw("PL");
      
      TString tmp = trk_tag[ibin2+1];
      tx = new TLatex(); tx->SetTextSize(.1);
      tx->DrawLatexNDC(0.15,0.85, tmp);
   }
  }

}