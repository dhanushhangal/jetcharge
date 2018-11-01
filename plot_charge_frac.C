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
//#include "assert.h"
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
#define njtptbins 5
#define trkeffbins 11

Double_t k_scale=1.;

const Double_t pt_low = 120.;
const Double_t pt_high = 600.;

const Double_t pt_low_up = 126.;
const Double_t pt_low_dn = 114.;

bool do_fit = true;
bool do_systematics = true;
bool do_ptcut_fitting = true;
bool do_kappa_fitting = false;
bool draw_tracking = false;
bool do_leading_track = false;
bool do_mean_sigma = false;
bool do_gen = false;

char saythis[500];

using namespace std;

TString cent[6] = {"0","1","2","3","4","5"};
TString trk[4] = {"0","1","2","3"};
TString kbin[4] = {"0","1","2","3"};
TString pt[56] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55"};

TString trkeff[11] = {"0","1","2","3","4","5","6","7","8","9","10"};
TString trkeffvary[trkeffbins] = {"trk eff -5%","trk eff -4%","trk eff -3%","trk eff -2%","trk eff -1%","trk eff nom.","trk eff +1%","trk eff +2%","trk eff +3%","trk eff +4%","trk eff +5%"};

TString trk_pt[4] = {"","2","4","5"};
TString kappa_str[4] = {"1","1p5","2","2p5"};

TString cent_tag_pp[] = {"pp ref","0-10% PbPb","10-30% PbPb","30-50% PbPb","50-100% PbPb"};
TString cent_tag_pp_MC[] = {"PYTHIA ref","0-10% P+H","10-30% P+H","30-50% P+H","50-100% P+H"};
TString x_label[] = {"x1","x2","x3"};
TString x_label_kappa[] = {"pT > 2 (#kappa=0.5); #kappa=0.5 (pT > 2)","pT > 4 (#kappa=0.5); #kappa=1 (pT > 2)","pT > 5 (#kappa=0.5); #kappa=2 (pT > 2)"};
TString cent_tag[] = {"0-10%","10-30%","30-50%","50-100%","70-100%"};
TString trk_tag[] = {"p_{T}^{trk} > 2 GeV, #kappa = 0.5","p_{T}^{trk} > 3 GeV","p_{T}^{trk} > 4 GeV, #kappa = 0.5", "p_{T}^{trk} > 5 GeV, #kappa = 0.5"};
//TString kappa_tag[] = {"#kappa = 0.5 , p_{T}^{trk} > 2 GeV","#kappa = 1 , p_{T}^{trk} > 2 GeV","#kappa = 1 , p_{T}^{trk} > 2 GeV", "#kappa = 2 , p_{T}^{trk} > 2 GeV"};
TString kappa_tag[] = {"#kappa = 0.3 , p_{T}^{trk} > 2 GeV","#kappa = 0.5 , p_{T}^{trk} > 2 GeV","#kappa = 0.5 , p_{T}^{trk} > 2 GeV", "#kappa = 0.7 , p_{T}^{trk} > 2 GeV"};

//Double_t jtpt_bound[njtptbins+1] = {120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.};
Double_t jtpt_bound[njtptbins+1] = {120.,130.,140.,160.,190.,500.};
Double_t jtpt_bound_err[njtptbins+1] = {0.,0.,0.,0.,0.,0.};
//Double_t jtpt_bound[njtptbins+1] = {120.,300.};
Double_t x_mean[ntrkbins-1] = {2.,4.,5.};
Double_t x_mean_err[ntrkbins-1] = {0.,0.,0.};
Double_t kappa[ntrkbins-1] = {0.3,0.5,0.7};
Double_t kappa_err[ntrkbins-1] = {0.,0.,0.};
Double_t cent_mean[nCBins] = {1.,2.,3.,4.,5.};
Double_t cent_mean_err[nCBins] = {0.1,0.1,0.1,0.1,0.1};

Double_t trkpt_bounds[ntrkbins] = {2.,4.,5.,6.};
Double_t kappa_bounds[ntrkbins] = {0.2,0.4,0.6,0.8};

Double_t diff_0trk_data_MC[nCBins][ntrkbins];
Double_t gluon[nCBins],quark[nCBins],up[nCBins],upbar[nCBins],down[nCBins],downbar[nCBins],gluon_ubar_dbar[nCBins], charm[nCBins], strange[nCBins], bottom[nCBins];
Double_t gluon_gen[nCBins],quark_gen[nCBins], down_gen[nCBins], up_gen[nCBins], gluon_ubar_dbar_gen[nCBins];
Double_t ubdbcsb[nCBins],gluon_ubdbcsb[nCBins], u_ud[nCBins], d_ud[nCBins];
//Double_t gluon_corr_factor[nCBins], quark_corr_factor[nCBins];
Double_t gluon_corr_factor[nCBins][ntrkbins], quark_corr_factor[nCBins][ntrkbins];
Double_t gluon_gen_int[nCBins][ntrkbins], quark_gen_int[nCBins][ntrkbins];
Double_t gluon_reco_int[nCBins][ntrkbins], quark_reco_int[nCBins][ntrkbins];
Double_t gluon_no0trk_int[nCBins][ntrkbins], quark_no0trk_int[nCBins][ntrkbins], quark_0trk_int[nCBins][ntrkbins], gluon_0trk_int[nCBins][ntrkbins]; Double_t gluon_no0trk_gen_int[nCBins][ntrkbins], quark_no0trk_gen_int[nCBins][ntrkbins];
Double_t data_MC_0trk_factor[nCBins][ntrkbins], data_int[nCBins], quark_data_0trk_int[nCBins][ntrkbins], gluon_data_0trk_int[nCBins][ntrkbins], data_0trk_int[nCBins][ntrkbins], data_no0trk_int[nCBins][ntrkbins];

  //TFile *closure_histos_pp_MC = TFile::Open("Pythia6_jetchg_bkgsub_kapppa_no0trks_eta0p5_1p5_Aug1.root");
  //TFile *closure_histos_pp_MC = TFile::Open("Pythia6_jetchg_kappa_no0trksanypT_Aug30.root");
  TFile *closure_histos_pp_MC = TFile::Open("Pythia6_jetchg_no0trkjets_20181022.root");
  //TFile *closure_histos_pp_MC = TFile::Open("Pythia6_photonjet_jetchg_bkgsub_kapppa_no0trks_eta0p5_1p5_Aug26.root");
  //TFile *closure_histos_PbPb_MC = TFile::Open("P+H_jetchg_cymbal_bkgsub_kapppa_recogen_no0trks_eta0p5_1p5_xiaotrkcorr_Aug3.root");
  //TFile *closure_histos_PbPb_MC = TFile::Open("P+H_jetchg_trkeffconst_ptcut_no0trkjets_20180913.root");
  TFile *closure_histos_PbPb_MC = TFile::Open("P+H_jetchg_trkeffconst_no0trkjets_20181021.root");
  //TFile *closure_histos_PbPb_MC = TFile::Open("P+H_jetchg_trkeffconst_ptcut_alltrkjets_20180913.root");
  //TFile *closure_histos_pp_data = TFile::Open("ppdata_jetchg_bkgsub_kapppa_no0trks_eta0p5_1p5_Aug1.root");
  //TFile *closure_histos_pp_data = TFile::Open("ppdata_jetchg_bkgsub_trkeffvary_kapppa_no0trks_eta0p5_1p5_Aug16.root");
  //TFile *closure_histos_pp_data = TFile::Open("ppdata_jetchg_kappa_no0trksanypT_Aug30.root");
  TFile *closure_histos_pp_data = TFile::Open("ppdata_jetchg_trkeffvary_no0trkjets_20181022.root");
  //TFile *closure_histos_pp_data = TFile::Open("ppdata_jetchg_postrkup_bkgsub_kapppa_no0trks_eta0p5_1p5_Aug13.root");
  //TFile *closure_histos_PbPb_data = TFile::Open("PbPbdata_jetchg_bkgsub_kapppa_no0trks_eta0p5_1p5_xiaotrkcorr_Aug3.root");
  //TFile *closure_histos_PbPb_data = TFile::Open("PbPbdata_jetchg_bkgsub_trkeffvary_peripheral5_kapppa_no0trks_eta0p5_1p5_xiaotrkcorr_Aug16.root");
  //TFile *closure_histos_PbPb_data = TFile::Open("PbPbdata_jetchg_trkeffvary_ptcut_no0trkjets_20180913.root");
  TFile *closure_histos_PbPb_data = TFile::Open("PbPbdata_jetchg_trkeffvary_no0trkjets_20181022.root");
  //TFile *closure_histos_PbPb_data = TFile::Open("PbPbdata_jetchg_trkeffvary_ptcut_no0trkjets_20180917.root");
  //TFile *closure_histos_PbPb_data = TFile::Open("PbPbdata_jetchg_trkeffconst_ptcut_alltrkjets_20180913.root");

/*
  //kappa
  TFile *closure_histos_pp_MC = TFile::Open("Pythia6_jetchg_kappa_no0trkjets_20181019.root");
  TFile *closure_histos_PbPb_MC = TFile::Open("P+H_jetchg_trkeffconst_kappa0p5_no0trkjets_20180917.root");
  TFile *closure_histos_pp_data = TFile::Open("ppdata_jetchg_trkeffvary_kappa_no0trkjets_20181019.root");
  TFile *closure_histos_PbPb_data = TFile::Open("PbPbdata_jetchg_trkeffvary_kappa0p5_no0trkjets_20181018.root");
*/
  TFile *closure_histos_pp_data_trkuncpos = TFile::Open("ppdata_jetchg_postrkup_bkgsub_kapppa_no0trks_eta0p5_1p5_Aug13.root");
  TFile *closure_histos_PbPb_data_trkuncpos = TFile::Open("PbPbdata_jetchg_trkeffposupvary_ptcut_no0trkjets_20180913.root");

  TFile *closure_histos_PbPb_data_kappa_trkuncpos = TFile::Open("PbPbdata_jetchg_trkeffposupvary_kappa0p5_no0trkjets_20180917.root");

  TFile *closure_histos_PbPb_data_jer = TFile::Open("PbPbdata_jetchg_JER_trkeffvary_ptcut_no0trkjets_20180913.root");
  TFile *closure_histos_PbPb_data_kappa_jer = TFile::Open("PbPbdata_jetchg_JER_trkeffvary_kappa0p5_no0trkjets_20180917.root");

  TFile *closure_histos_PbPb_MC_gensube0 = TFile::Open("P+H_jetchg_cymbal_bkgsub_kapppa_no0trks_eta0p5_1p5_Jul21.root");

  TFile *closure_histos_pp_pythia_gen = TFile::Open("Pythia6_gen_jetchg_bkgsub_kapppa_no0trks_eta0p5_1p5_Aug3.root");
  TFile *closure_histos_pp_herwig_gen = TFile::Open("HERWIG_gen_jetchg_bkgsub_kapppa_no0trks_eta0p5_1p5_Aug2.root");

  TFile *f_trk_pp_MC = TFile::Open("skims/preapp/Pythia6_chg_Calojets_Jul27.root");
  //TFile *f_trk_PbPb_MC = TFile::Open("skims/preapp/Pythia6Hydjet_PbPbMC_xiaotrkcorr_jetchg_Aug3.root");
  TFile *f_trk_pp_data = TFile::Open("skims/preapp/ppdata_chg_Calojets_Jul31.root");
  TFile *f_trk_PbPb_data = TFile::Open("skims/preapp/PbPbdata_xiaotrkcorr_jetchg_Aug3.root");

  TFile *f_trk_PbPb_data_MB = TFile::Open("skims/preapp/PbPbdata_chg_MB_notrkcorr.root");
  TFile *f_trk_PbPb_MC_MB = TFile::Open("skims/preapp/HydJet_MB_Aug3.root");

  TH1D *h_dummy = new TH1D("h_dummy","",10,1,6);
  TH1D *h_dummy2 = new TH1D("h_dummy2","",5,1,6);

  TH1D *h_dummy3 = new TH1D("h_dummy3","",njtptbins,jtpt_bound);
  TH1D *h_dummy4 = new TH1D("h_dummy4","",njtptbins,jtpt_bound);
  
  TH1D *h_dummy5 = new TH1D("h_dummy5","",5,0.,2.5);

  TH1D *h_dummy7 = new TH1D("h_dummy7","",5,1,6);
  TH1D *h_dummy8 = new TH1D("h_dummy8","",5,1,6);
  TH1D *h_dummy9 = new TH1D("h_dummy9","",5,1,6);
  TH1D *h_dummy10 = new TH1D("h_dummy10","",5,1,6);

  TH1D *h_0trks_data[nCBins];
  TH1D *h_0trks_MC[nCBins];  
  TH1D *h_corr_factor[nCBins];  
  TH1D *h_syst[nCBins];  
  TH1D *h_syst_trk[nCBins];  
  TH1D *h_syst_posneg[nCBins];  
  TH1D *h_syst_JER[nCBins];  
  TH1D *h_syst_0trk[nCBins];  

  TH2F *h_chg_refpt[nCBins][ntrkbins];
  TH2F *h_chg_refpt_q[nCBins][ntrkbins];
  TH2F *h_chg_refpt_g[nCBins][ntrkbins];
  TH2F *h_chg_refpt_u[nCBins][ntrkbins];
  TH2F *h_chg_refpt_sube0[nCBins][ntrkbins];

  TH2F *h_chg_refpt_pythia[nCBins][ntrkbins];
  TH2F *h_chg_refpt_herwig[nCBins][ntrkbins];

  TH2F *h_trk_corrpt_MC[nCBins][ntrkbins];
  TH2F *h_trk_corrpt_data[nCBins][ntrkbins];
  TH2F *h_bkgtrk_corrpt_MC[nCBins][ntrkbins];
  TH2F *h_bkgtrk_corrpt_data[nCBins][ntrkbins];

  TH2F *h_bkggen_corrpt_MC[nCBins];
  TH2F *h_subegen_corrpt_MC[nCBins];

  TH2F *h_chg_corrpt[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_q[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_g[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_u[nCBins][ntrkbins];

  TH2F *h_chg_corrpt_kappa[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_q_kappa[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_g_kappa[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_u_kappa[nCBins][ntrkbins];

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

  TH2F *h_chg_corrpt_c_kappa[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_s_kappa[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_b_kappa[nCBins][ntrkbins];

  TH2F *h_chg_corrpt_up_kappa[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_down_kappa[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_upb_kappa[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_downb_kappa[nCBins][ntrkbins];

  TH2F *h_chg_corrpt_bkg[nCBins][ntrkbins];
  TH2F *h_chg_refpt_bkg[nCBins][ntrkbins];
  TH2F *h_chg_refpt_bkg_sube0[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_bkg_kappa[nCBins][ntrkbins];

  TH2F *h_chg_corrpt_highest[nCBins];
  TH2F *h_chg_corrpt_highest_2[nCBins];

  TH2F *h_chg_corrpt_data[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_data_bkg[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_spillover[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_spilloverdn[nCBins][ntrkbins];

  TH2F *h_chg_corrpt_data_trkeff[nCBins][ntrkbins][trkeffbins];
  TH1D *h_chg_data_trkeff[nCBins][ntrkbins][trkeffbins];
  TH2F *h_chg_corrpt_data_kappa_trkeff[nCBins][ntrkbins][trkeffbins];
  TH1D *h_chg_data_kappa_trkeff[nCBins][ntrkbins][trkeffbins];

  TH2F *h_chg_corrpt_data_kappa[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_data_bkg_kappa[nCBins][ntrkbins];

  TH2F *h_chg_corrpt_data_trkpos[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_data_jer[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_data_trkpos_kappa[nCBins][ntrkbins];
  TH2F *h_chg_corrpt_data_jer_kappa[nCBins][ntrkbins];

  TH1D *h_chg_ref[nCBins][ntrkbins];
  TH1D *h_chg_ref_sube0[nCBins][ntrkbins];
  //TH1D *h_chg_ref_q[nCBins][ntrkbins];
  //TH1D *h_chg_ref_g[nCBins][ntrkbins];
  TH1D *h_chg_ref_u[nCBins][ntrkbins];

  TH1D *h_chg_ref_pythia[nCBins][ntrkbins];
  TH1D *h_chg_ref_herwig[nCBins][ntrkbins];
  TH1D *h_chg_ref_pythia_herwig[nCBins][ntrkbins];

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

  TH1D *h_reco_corr_no0trks[nCBins][ntrkbins];
  TH1D *h_reco_corr_no0trks_q[nCBins][ntrkbins];
  TH1D *h_reco_corr_no0trks_g[nCBins][ntrkbins];
  TH1D *h_gen_full_no0trks[nCBins][ntrkbins];
  TH1D *h_gen_full_no0trks_q[nCBins][ntrkbins];
  TH1D *h_gen_full_no0trks_g[nCBins][ntrkbins];

  TH1D *h_reco_corr_0trks[nCBins][ntrkbins];
  TH1D *h_reco_corr_0trks_q[nCBins][ntrkbins];
  TH1D *h_reco_corr_0trks_g[nCBins][ntrkbins];

  TH1D *h_gen_full_0trks[nCBins][ntrkbins];
  TH1D *h_gen_full_0trks_q[nCBins][ntrkbins];
  TH1D *h_gen_full_0trks_g[nCBins][ntrkbins];

  TH1D *h_reco_data_no0trks[nCBins][ntrkbins];
  TH1D *h_reco_data_0trks[nCBins][ntrkbins];

  TH1D *h_gen_pt[nCBins];
  TH1D *h_gen_pt_q[nCBins];
  TH1D *h_gen_pt_g[nCBins];

  TH1D *h_reco_MC_c[nCBins];
  TH1D *h_reco_MC_s[nCBins];
  TH1D *h_reco_MC_b[nCBins];
  TH1D *h_reco_MC_ubdbcsb[nCBins];

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

  TH1D *h_chg_ref_bkg[nCBins][ntrkbins];
  TH1D *h_chg_ref_bkg_sube0[nCBins][ntrkbins];
  TH1D *h_chg_spill[nCBins][ntrkbins];

  TH1D *h_chg_highest[nCBins];
  TH1D *h_chg_highest_2[nCBins];

  //TH1D *h_chg_reco_up[nCBins][ntrkbins];
  //TH1D *h_chg_reco_down[nCBins][ntrkbins];
  TH1D *h_chg_reco_upb[nCBins][ntrkbins];
  TH1D *h_chg_reco_downb[nCBins][ntrkbins];

  //TH1D *h_chg_ref_up[nCBins][ntrkbins];
  //TH1D *h_chg_ref_down[nCBins][ntrkbins];
  //TH1D *h_chg_ref_upb[nCBins][ntrkbins];
  //TH1D *h_chg_ref_downb[nCBins][ntrkbins];

  //TH1D *h_chg_reco_gubdb[nCBins][ntrkbins];
  TH1D *h_chg_reco_others[nCBins][ntrkbins];

  TH1D *h_chg_reco_kappa[nCBins][ntrkbins];
  TH1D *h_chg_reco_q_kappa[nCBins][ntrkbins];
  //TH1D *h_chg_reco_g_kappa[nCBins][ntrkbins];
  TH1D *h_chg_reco_u_kappa[nCBins][ntrkbins];
  TH1D *h_chg_reco_bkg_kappa[nCBins][ntrkbins];

  //TH1D *h_chg_reco_up_kappa[nCBins][ntrkbins];
  //TH1D *h_chg_reco_down_kappa[nCBins][ntrkbins];
  TH1D *h_chg_reco_upb_kappa[nCBins][ntrkbins];
  TH1D *h_chg_reco_downb_kappa[nCBins][ntrkbins];
  TH1D *h_chg_reco_others_kappa[nCBins][ntrkbins];

  TH1D *h_chg_reco_q_refl[nCBins][ntrkbins];
  TH1D *h_chg_reco_g_refl[nCBins][ntrkbins];
  TH1D *h_chg_reco_MC_refl[nCBins][ntrkbins];
  TH1D *h_chg_reco_data_refl[nCBins][ntrkbins];

  TH1D *h_chg_reco_uubar[nCBins][ntrkbins];
  TH1D *h_chg_reco_ddbar[nCBins][ntrkbins];

  TH1D *h_chg_data[nCBins][ntrkbins];
  TH1D *h_chg_data_bkg[nCBins][ntrkbins];
  TH1D *h_chg_data_spillover[nCBins][ntrkbins];
  TH1D *h_chg_data_spilloverdn[nCBins][ntrkbins];

  TH1D *h_chg_data_sys[nCBins][ntrkbins];

  TH1D *h_chg_data_trkpos[nCBins][ntrkbins];
  TH1D *h_chg_data_jer[nCBins][ntrkbins];
  TH1D *h_chg_data_trkpos_kappa[nCBins][ntrkbins];
  TH1D *h_chg_data_jer_kappa[nCBins][ntrkbins];

  TH1D *h_chg_data_jeup[nCBins][ntrkbins];
  TH1D *h_chg_data_jedn[nCBins][ntrkbins];

  TH1D *h_chg_data_kappa[nCBins][ntrkbins];
  TH1D *h_chg_data_bkg_kappa[nCBins][ntrkbins];

  TH1D *h_chg_reco_cop[nCBins][ntrkbins];
  TH1D *h_chg_reco_up_cop[nCBins][ntrkbins];
  TH1D *h_chg_reco_down_cop[nCBins][ntrkbins];
  TH1D *h_chg_reco_gubdb_cop[nCBins][ntrkbins];
  TH1D *h_chg_data_cop[nCBins][ntrkbins];
  TH1D *h_chg_reco_g_cop[nCBins][ntrkbins];
  TH1D *h_chg_reco_others_cop[nCBins][ntrkbins];

  TH1D *h_chg_reco_data_ratio[nCBins][ntrkbins];

  TH1D *h_chg_ref_ratio[nCBins][ntrkbins];
  TH1D *h_chg_reco_ratio[nCBins][ntrkbins];
  TH1D *h_chg_reco_ratio_spillover[nCBins][ntrkbins];
  TH1D *h_chg_reco_ratio_spilloverdn[nCBins][ntrkbins];
  TH1D *h_chg_reco_ratio_2param[nCBins][ntrkbins];
  TH1D *h_chg_data_ratio[nCBins][ntrkbins];
  TH1D *h_chg_reco_ratio_trkeff[nCBins][ntrkbins][trkeffbins];
  TH1D *h_chg_reco_ratio_kappa_trkeff[nCBins][ntrkbins][trkeffbins];

  TH1D *h_chg_reco_ratio_jer[nCBins][ntrkbins];

  TH1D *h_chg_reco_cop_kappa[nCBins][ntrkbins];
  TH1D *h_chg_reco_up_cop_kappa[nCBins][ntrkbins];
  TH1D *h_chg_reco_down_cop_kappa[nCBins][ntrkbins];
  TH1D *h_chg_reco_gubdb_cop_kappa[nCBins][ntrkbins];
  TH1D *h_chg_data_cop_kappa[nCBins][ntrkbins];
  TH1D *h_chg_reco_g_cop_kappa[nCBins][ntrkbins];
  TH1D *h_chg_reco_others_cop_kappa[nCBins][ntrkbins];
  TH1D *h_chg_reco_ratio_kappa[nCBins][ntrkbins];
  TH1D *h_chg_data_ratio_kappa[nCBins][ntrkbins];

  TH1D *h_pt_reco_q[nCBins][ntrkbins];
  TH1D *h_pt_reco_g[nCBins][ntrkbins];
  TH1D *h_pt_ref_q[nCBins][ntrkbins];
  TH1D *h_pt_ref_g[nCBins][ntrkbins];

  TH1D *h_chg_qg_ratio[nCBins][ntrkbins];

  TH1D *h_trk_MC[nCBins][ntrkbins];
  TH1D *h_trk_data[nCBins][ntrkbins];
  TH1D *h_bkgtrk_MC[nCBins][ntrkbins];
  TH1D *h_bkgtrk_data[nCBins][ntrkbins];
  TH1D *h_sube0trk_MC[nCBins][ntrkbins];

  TH1D *h_bkggen_MC[nCBins];
  TH1D *h_subegen_MC[nCBins];

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

  TGraphErrors *q_reco_0trk[nCBins];

  TGraphErrors *up_reco[nCBins];
  TGraphErrors *down_reco[nCBins];
  TGraphErrors *gluon_reco[nCBins];
  TGraphErrors *quark_reco[nCBins];

  TGraphErrors *up_reco_kappa[nCBins];
  TGraphErrors *down_reco_kappa[nCBins];
  TGraphErrors *gluon_reco_kappa[nCBins];
  TGraphErrors *quark_reco_kappa[nCBins];

  TGraphErrors *up_ref[nCBins];
  TGraphErrors *down_ref[nCBins];
  TGraphErrors *gluon_ref[nCBins];

  TGraphErrors *up_data[nCBins];
  TGraphErrors *down_data[nCBins];
  TGraphErrors *gluon_data[nCBins];
  TGraphErrors *quark_data[nCBins];

  TGraphAsymmErrors *gluon_data_asymm_sys[nCBins];
  TGraphErrors *gluon_data_sys[nCBins];

  TGraphErrors *up_data_kappa[nCBins];
  TGraphErrors *down_data_kappa[nCBins];
  TGraphErrors *gluon_data_kappa[nCBins];
  TGraphErrors *quark_data_kappa[nCBins];

  TGraphErrors *gluon_ratio_data[nCBins];
  TGraphErrors *quark_ratio_data[nCBins];

  TGraphErrors *gluon_ratio_reco[nCBins];
  TGraphErrors *quark_ratio_reco[nCBins];

  TGraphErrors *gluon_sys[nCBins];
  TGraphErrors *quark_sys[nCBins];
  TGraphErrors *gluon_ncs[nCBins];
  TGraphErrors *quark_ncs[nCBins];

  TGraphErrors *gluon_cent[ntrkbins-1];
  TGraphErrors *quark_cent[ntrkbins-1];
  TGraphErrors *gluon_centsys_gr[ntrkbins-1];
  TGraphErrors *quark_centsys_gr[ntrkbins-1];

  TGraph *gluon_ncs_data[nCBins];
  TGraph *quark_ncs_data[nCBins];

  TGraph *gluon_trkup_gr[nCBins];
  TGraph *quark_trkup_gr[nCBins];
  TGraph *gluon_trkdn_gr[nCBins];
  TGraph *quark_trkdn_gr[nCBins];
  TGraph *gluon_trkpos_gr[nCBins];
  TGraph *quark_trkpos_gr[nCBins];

  TGraph *gluon_spillup_gr[nCBins];
  TGraph *gluon_spilldn_gr[nCBins];

  TF1 *f_chg_ref[nCBins][ntrkbins];  
  TF1 *f_chg_reco[nCBins][ntrkbins];
  TF1 *f_chg_qg_data[nCBins][ntrkbins];
  TF1 *f_chg_udg_data[nCBins][ntrkbins];
  TF1 *f_chg_udg_ref[nCBins][ntrkbins];
  TF1 *f_chg_udg_kappa[nCBins][ntrkbins];
  TF1 *f_chg_qg_kappa[nCBins][ntrkbins];

  TF1 *f_eta_reco[nCBins];

  //TH1D *h_chg_data_avg[nCBins][njtptbins][ntrkbins];
  TH1D *h_chg_reco_avg[nCBins][njtptbins][ntrkbins];
  TH1D *h_chg_reco_q_avg[nCBins][njtptbins][ntrkbins];
  TH1D *h_chg_reco_g_avg[nCBins][njtptbins][ntrkbins];
  TH1D *h_chg_reco_up_avg[nCBins][njtptbins][ntrkbins];
  TH1D *h_chg_reco_down_avg[nCBins][njtptbins][ntrkbins];

  TF1 *f_chg_data_avg[nCBins][njtptbins][ntrkbins];  
  TF1 *f_chg_reco_avg[nCBins][njtptbins][ntrkbins];  
  TF1 *f_chg_reco_q_avg[nCBins][njtptbins][ntrkbins];
  TF1 *f_chg_reco_g_avg[nCBins][njtptbins][ntrkbins];
  TF1 *f_chg_reco_up_avg[nCBins][njtptbins][ntrkbins];
  TF1 *f_chg_reco_down_avg[nCBins][njtptbins][ntrkbins];

  THStack *h_chg_data_stk[nCBins][ntrkbins];
  THStack *h_chg_reco_stk[nCBins][ntrkbins];

  THStack *h_chg_data_stk_kappa[nCBins][ntrkbins];
/////trk histos
  TH1D *h_genpt_MC[nCBins];
  TH1D *h_bkggenpt_MC[nCBins];
  TH1D *h_genpt_MC_full[nCBins];
  TH1D *h_genpt_MC_pos[nCBins];
  TH1D *h_genpt_MC_neg[nCBins];

  TH1D *h_trkpt_MC_full[nCBins];
  TH1D *h_trkpt_MC_pos[nCBins];
  TH1D *h_trkpt_MC_neg[nCBins];
  TH1D *h_trkpt_data_pos[nCBins];
  TH1D *h_trkpt_data_neg[nCBins];

  TH1D *h_trkpt_MC[nCBins][ntrkbins];
  TH1D *h_trkpt_data[nCBins][ntrkbins];
  TH1D *h_trkpt_MC_sub[nCBins][ntrkbins];
  TH1D *h_bkgtrkpt_MC[nCBins][ntrkbins];
  TH1D *h_bkgtrkpt_data[nCBins][ntrkbins]; 
  TH1D *h_trkpt_data_sub[nCBins][ntrkbins];
  
  TH1D *h_trkk_MC[nCBins][ntrkbins];
  TH1D *h_trkk_data[nCBins][ntrkbins];
  TH1D *h_trkk_MC_sub[nCBins][ntrkbins];
  TH1D *h_bkgtrkk_MC[nCBins][ntrkbins];
  TH1D *h_bkgtrkk_data[nCBins][ntrkbins];
  TH1D *h_trkk_data_sub[nCBins][ntrkbins];

  TH1D *h_trkpt_MC_full_ratio[nCBins];
  TH1D *h_trkpt_MC_pos_ratio[nCBins];
  TH1D *h_trkpt_MC_neg_ratio[nCBins];

  TH1D *h_ntrk_MC_MB;
  TH1D *h_ntrk_data_MB;

void get_histos(bool do_data = false){

  h_ntrk_MC_MB = (TH1D*)f_trk_PbPb_MC_MB->Get((TString)("h_ntrk"))->Clone("h_ntrk_MC_MB");
  h_ntrk_data_MB = (TH1D*)f_trk_PbPb_data_MB->Get((TString)("h_ntrk"))->Clone("h_ntrk_data_MB");

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

      h_genpt_MC_full[ibin] = (TH1D*)f_trk_pp_MC->Get((TString)("h_trk_pt_gen_full"));
      h_genpt_MC_pos[ibin] = (TH1D*)f_trk_pp_MC->Get((TString)("h_trk_pt_gen_pos"));
      h_genpt_MC_neg[ibin] = (TH1D*)f_trk_pp_MC->Get((TString)("h_trk_pt_gen_neg"));

      h_trkpt_MC_full[ibin] = (TH1D*)f_trk_pp_MC->Get((TString)("h_trk_pt_reco_full"));
      h_trkpt_MC_pos[ibin] = (TH1D*)f_trk_pp_MC->Get((TString)("h_trk_pt_reco_pos"));
      h_trkpt_MC_neg[ibin] = (TH1D*)f_trk_pp_MC->Get((TString)("h_trk_pt_reco_neg"));

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
/*
      h_genpt_MC_full[ibin] = (TH1D*)f_trk_PbPb_MC->Get((TString)("h_trk_pt_gen_full_"+cent[ibin-1]));
      h_genpt_MC_pos[ibin] = (TH1D*)f_trk_PbPb_MC->Get((TString)("h_trk_pt_gen_pos_"+cent[ibin-1]));
      h_genpt_MC_neg[ibin] = (TH1D*)f_trk_PbPb_MC->Get((TString)("h_trk_pt_gen_neg_"+cent[ibin-1]));

      h_trkpt_MC_full[ibin] = (TH1D*)f_trk_PbPb_MC->Get((TString)("h_trk_pt_reco_full"+cent[ibin-1]));
      h_trkpt_MC_pos[ibin] = (TH1D*)f_trk_PbPb_MC->Get((TString)("h_trk_pt_reco_pos"+cent[ibin-1]));
      h_trkpt_MC_neg[ibin] = (TH1D*)f_trk_PbPb_MC->Get((TString)("h_trk_pt_reco_neg"+cent[ibin-1]));
*/
      h_bkggen_corrpt_MC[ibin] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_bkgtrk_refpt_cent"+cent[ibin-1]));
      h_subegen_corrpt_MC[ibin] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_sube0trk_refpt_cent"+cent[ibin-1])); 
      h_bkggen_MC[ibin] = h_bkggen_corrpt_MC[ibin]-> ProjectionY(Form("h_bkggen_MC_%d",ibin),h_bkggen_corrpt_MC[ibin]->GetXaxis()->FindBin(pt_low),h_bkggen_corrpt_MC[ibin]->GetXaxis()->FindBin(pt_high),"");
      h_subegen_MC[ibin] = h_subegen_corrpt_MC[ibin]-> ProjectionY(Form("h_subegen_MC_%d",ibin),h_subegen_corrpt_MC[ibin]->GetXaxis()->FindBin(pt_low),h_subegen_corrpt_MC[ibin]->GetXaxis()->FindBin(pt_high),"");
 
      if(do_leading_track){
        h_chg_corrpt_highest[ibin] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_highest_cent"+cent[ibin-1])); 
        h_chg_corrpt_highest_2[ibin] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_highest_2_cent"+cent[ibin-1])); 
  
        h_chg_highest[ibin] = h_chg_corrpt_highest[ibin]-> ProjectionY(Form("h_chg_highest_%d",ibin),h_chg_corrpt_highest[ibin]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_highest[ibin]->GetXaxis()->FindBin(pt_high),"");
        h_chg_highest_2[ibin] = h_chg_corrpt_highest_2[ibin]-> ProjectionY(Form("h_chg_highest_2_%d",ibin),h_chg_corrpt_highest[ibin]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_highest[ibin]->GetXaxis()->FindBin(pt_high),"");
      }
    }

    //normalizing eta spectra
    h_eta_MC[ibin]->Scale(1./h_eta_MC[ibin]->Integral());
    h_eta_MC_q[ibin]->Scale(1./h_eta_MC_q[ibin]->Integral());
    h_eta_MC_g[ibin]->Scale(1./h_eta_MC_g[ibin]->Integral());
    h_eta_data[ibin]->Scale(1./h_eta_data[ibin]->Integral());
    h_eta_MC_up[ibin]->Scale(1./h_eta_MC_up[ibin]->Integral());
    h_eta_MC_upbar[ibin]->Scale(1./h_eta_MC_upbar[ibin]->Integral());
    h_eta_MC_d[ibin]->Scale(1./h_eta_MC_d[ibin]->Integral());
    h_eta_MC_dbar[ibin]->Scale(1./h_eta_MC_dbar[ibin]->Integral());

    h_reco_MC_ubdbcsb[ibin] = (TH1D*)h_reco_MC_upbar[ibin]->Clone((TString)("h_reco_MC_ubdbcsb"+cent[ibin]));
    h_reco_MC_ubdbcsb[ibin]->Add(h_reco_MC_dbar[ibin]);
    h_reco_MC_ubdbcsb[ibin]->Add(h_reco_MC_c[ibin]);
    h_reco_MC_ubdbcsb[ibin]->Add(h_reco_MC_s[ibin]);
    h_reco_MC_ubdbcsb[ibin]->Add(h_reco_MC_b[ibin]);
/*
    //normalizing jet spectra
    h_reco_MC[ibin]->Scale(1./h_reco_MC[ibin]->Integral());
    h_reco_data[ibin]->Scale(1./h_reco_data[ibin]->Integral());
    h_reco_MC_q[ibin]->Scale(1./h_reco_MC_q[ibin]->Integral());
    h_reco_MC_g[ibin]->Scale(1./h_reco_MC_g[ibin]->Integral());
    h_reco_MC_up[ibin]->Scale(1./h_reco_MC_up[ibin]->Integral());
    h_reco_MC_upbar[ibin]->Scale(1./h_reco_MC_upbar[ibin]->Integral());
    h_reco_MC_d[ibin]->Scale(1./h_reco_MC_d[ibin]->Integral());
    h_reco_MC_dbar[ibin]->Scale(1./h_reco_MC_dbar[ibin]->Integral());
    h_reco_MC_c[ibin]->Scale(1./h_reco_MC_c[ibin]->Integral());
    h_reco_MC_s[ibin]->Scale(1./h_reco_MC_s[ibin]->Integral());
    h_reco_MC_b[ibin]->Scale(1./h_reco_MC_b[ibin]->Integral());
*/
  }

  for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    for(int ibin=0;ibin<nCBins;ibin++){
      if(ibin==0){
        h_reco_corr_no0trks[ibin][ibin2] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_no0trks_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_reco_corr_no0trks_q[ibin][ibin2] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_no0trks_q_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_reco_corr_no0trks_g[ibin][ibin2] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_no0trks_g_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_gen_full_no0trks[ibin][ibin2] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_gen_full_no0trks_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_gen_full_no0trks_q[ibin][ibin2] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_gen_full_no0trks_q_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_gen_full_no0trks_g[ibin][ibin2] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_gen_full_no0trks_g_cent"+cent[ibin]+"_trk"+trk[ibin2]));

        h_reco_data_no0trks[ibin][ibin2] = (TH1D*)closure_histos_pp_data->Get((TString)("h_reco_corr_no0trks_cent"+cent[ibin]+"_trk"+trk[ibin2]));

        h_chg_refpt[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
        h_chg_refpt_q[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_q_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_refpt_g[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_g_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_refpt_u[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_u_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_refpt_bkg[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_bkg_cent"+cent[ibin]+"_trk"+trk[ibin2])); 

        h_chg_refpt_up[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_up_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_refpt_down[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_d_cent"+cent[ibin]+"_trk"+trk[ibin2]));      
        h_chg_refpt_upb[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_upbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        h_chg_refpt_downb[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_refpt_dbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));

        h_chg_refpt_pythia[ibin][ibin2] = (TH2F*)closure_histos_pp_pythia_gen->Get((TString)("h_chg_refpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
        h_chg_refpt_herwig[ibin][ibin2] = (TH2F*)closure_histos_pp_herwig_gen->Get((TString)("h_chg_refpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 

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
        h_chg_corrpt_data_trkpos[ibin][ibin2] = (TH2F*)closure_histos_pp_data_trkuncpos->Get((TString)("h_chg_corrpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
         
        for(int ibin3=0; ibin3<trkeffbins;ibin3++){
            h_chg_corrpt_data_trkeff[ibin][ibin2][ibin3] = (TH2F*)closure_histos_pp_data->Get((TString)("h_chg_corrpt_trkeff_cent"+cent[ibin]+"_trk"+trk[ibin2]+"_eff"+trkeff[ibin3]));           
            if(do_kappa_fitting) h_chg_corrpt_data_kappa_trkeff[ibin][ibin2][ibin3] = (TH2F*)closure_histos_pp_data->Get((TString)("h_chg_corrpt_kappa_trkeff_cent"+cent[ibin]+"_trk"+trk[ibin2]+"_eff"+trkeff[ibin3]));           
        }

        h_chg_corrpt_bkg[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_bkg_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
        h_chg_corrpt_data_bkg[ibin][ibin2] = (TH2F*)closure_histos_pp_data->Get((TString)("h_chg_corrpt_bkg_cent"+cent[ibin]+"_trk"+trk[ibin2]));

        if(draw_tracking){
            h_trk_corrpt_data[ibin][ibin2] = (TH2F*)closure_histos_pp_data->Get((TString)("h_trk_corrpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
            h_trk_corrpt_MC[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_trk_corrpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
            h_bkgtrk_corrpt_data[ibin][ibin2] = (TH2F*)closure_histos_pp_data->Get((TString)("h_bkgtrk_corrpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
            h_bkgtrk_corrpt_MC[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_bkgtrk_corrpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
        }

        h_chg_ref_pythia[ibin][ibin2] = h_chg_refpt_pythia[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_pythia_%d_%d",ibin,ibin2),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_herwig[ibin][ibin2] = h_chg_refpt_herwig[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_herwig_%d_%d",ibin,ibin2),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_pythia[ibin][ibin2]->Scale(1./(h_chg_ref_pythia[ibin][ibin2]->Integral()));
        h_chg_ref_herwig[ibin][ibin2]->Scale(1./(h_chg_ref_herwig[ibin][ibin2]->Integral()));

        h_trkpt_MC[ibin][ibin2] = (TH1D*)f_trk_pp_MC->Get((TString)("h_trk_pt"+trk_pt[ibin2]+"_reco"));
        h_bkgtrkpt_MC[ibin][ibin2] = (TH1D*)f_trk_pp_MC->Get((TString)("h_trk_pt"+trk_pt[ibin2]+"_reco_refl"));
        h_trkk_MC[ibin][ibin2] = (TH1D*)f_trk_pp_MC->Get((TString)("h_trk_ptk"+kappa_str[ibin2]+"_reco"));
        h_bkgtrkk_MC[ibin][ibin2] = (TH1D*)f_trk_pp_MC->Get((TString)("h_trk_ptk"+kappa_str[ibin2]+"_reco_refl"));

        h_trkpt_data[ibin][ibin2] = (TH1D*)f_trk_pp_data->Get((TString)("h_trk_pt"+trk_pt[ibin2]+"_reco"));
        h_bkgtrkpt_data[ibin][ibin2] = (TH1D*)f_trk_pp_data->Get((TString)("h_trk_pt"+trk_pt[ibin2]+"_reco_refl"));
        h_trkk_data[ibin][ibin2] = (TH1D*)f_trk_pp_data->Get((TString)("h_trk_ptk"+kappa_str[ibin2]+"_reco"));
        h_bkgtrkk_data[ibin][ibin2] = (TH1D*)f_trk_pp_data->Get((TString)("h_trk_ptk"+kappa_str[ibin2]+"_reco_refl"));

        if(do_kappa_fitting){
            h_chg_corrpt_kappa[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_kappa_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
            h_chg_corrpt_q_kappa[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_kappa_q_cent"+cent[ibin]+"_trk"+trk[ibin2]));
            h_chg_corrpt_g_kappa[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_kappa_g_cent"+cent[ibin]+"_trk"+trk[ibin2]));
            
            h_chg_corrpt_up_kappa[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_kappa_up_cent"+cent[ibin]+"_trk"+trk[ibin2]));
            h_chg_corrpt_down_kappa[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_kappa_d_cent"+cent[ibin]+"_trk"+trk[ibin2]));      
            h_chg_corrpt_upb_kappa[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_kappa_upbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));
            h_chg_corrpt_downb_kappa[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_kappa_dbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));

            h_chg_corrpt_c_kappa[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_kappa_c_cent"+cent[ibin]+"_trk"+trk[ibin2]));
            h_chg_corrpt_s_kappa[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_kappa_s_cent"+cent[ibin]+"_trk"+trk[ibin2]));
            h_chg_corrpt_b_kappa[ibin][ibin2] = (TH2F*)closure_histos_pp_MC->Get((TString)("h_chg_corrpt_kappa_b_cent"+cent[ibin]+"_trk"+trk[ibin2]));

            h_chg_corrpt_data_kappa[ibin][ibin2] = (TH2F*)closure_histos_pp_data->Get((TString)("h_chg_corrpt_kappa_cent"+cent[ibin]+"_trk"+trk[ibin2]));
        }

      }
      else{
        h_reco_corr_no0trks[ibin][ibin2] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_no0trks_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_reco_corr_no0trks_q[ibin][ibin2] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_no0trks_q_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_reco_corr_no0trks_g[ibin][ibin2] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_no0trks_g_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_gen_full_no0trks[ibin][ibin2] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_gen_full_no0trks_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_gen_full_no0trks_q[ibin][ibin2] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_gen_full_no0trks_q_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        h_gen_full_no0trks_g[ibin][ibin2] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_gen_full_no0trks_g_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

        h_reco_data_no0trks[ibin][ibin2] = (TH1D*)closure_histos_PbPb_data->Get((TString)("h_reco_corr_no0trks_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

        h_chg_refpt[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2]))->Clone((TString)("h_chg_refpt_"+cent[ibin-1]+"_"+trk[ibin2])); 
        h_chg_refpt_bkg[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_refpt_bkg_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
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
        h_chg_corrpt_data_trkpos[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data_trkuncpos->Get((TString)("h_chg_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 

        h_chg_corrpt_data_jer[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data_jer->Get((TString)("h_chg_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));  
 
        h_chg_corrpt_bkg[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_bkg_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
        h_chg_corrpt_data_bkg[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_chg_corrpt_bkg_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

        h_chg_data_jer[ibin][ibin2] = h_chg_corrpt_data_jer[ibin][ibin2]-> ProjectionY(Form("h_chg_data_jer_%d_%d",ibin,ibin2),h_chg_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_data_jer[ibin][ibin2]->Scale(1./h_chg_data_jer[ibin][ibin2]->Integral());
 
        for(int ibin3=0; ibin3<trkeffbins;ibin3++){
            h_chg_corrpt_data_trkeff[ibin][ibin2][ibin3] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_chg_corrpt_trkeff_cent"+cent[ibin-1]+"_trk"+trk[ibin2]+"_eff"+trkeff[ibin3]));           
            if(do_kappa_fitting) h_chg_corrpt_data_kappa_trkeff[ibin][ibin2][ibin3] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_chg_corrpt_kappa_trkeff_cent"+cent[ibin-1]+"_trk"+trk[ibin2]+"_eff"+trkeff[ibin3]));           
        }

        if(draw_tracking){
            h_trk_corrpt_data[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_trk_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
            h_trk_corrpt_MC[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_trk_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
            h_bkgtrk_corrpt_data[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_bkgtrk_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
            h_bkgtrk_corrpt_MC[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_bkgtrk_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
        }
/*
        h_trkpt_MC[ibin][ibin2] = (TH1D*)f_trk_PbPb_MC->Get((TString)("h_trk_pt"+trk_pt[ibin2]+"_reco_"+cent[ibin-1]));
        if(ibin2==0) h_bkgtrkpt_MC[ibin][ibin2] = (TH1D*)f_trk_PbPb_MC->Get((TString)("h_trk_pt"+trk_pt[ibin2]+"_reco_refl_"+cent[ibin-1]));
        else h_bkgtrkpt_MC[ibin][ibin2] = (TH1D*)f_trk_PbPb_MC->Get((TString)("h_trk_pt"+trk_pt[ibin2]+"_reco_refl"+cent[ibin-1]));
        h_trkk_MC[ibin][ibin2] = (TH1D*)f_trk_PbPb_MC->Get((TString)("h_trk_ptk"+kappa_str[ibin2]+"_reco_"+cent[ibin-1]));
        h_bkgtrkk_MC[ibin][ibin2] = (TH1D*)f_trk_PbPb_MC->Get((TString)("h_trk_ptk"+kappa_str[ibin2]+"_reco_refl_"+cent[ibin-1]));
*/
        h_trkpt_data[ibin][ibin2] = (TH1D*)f_trk_PbPb_data->Get((TString)("h_trk_pt"+trk_pt[ibin2]+"_reco_"+cent[ibin-1]));
        if(ibin2==0) h_bkgtrkpt_data[ibin][ibin2] = (TH1D*)f_trk_PbPb_data->Get((TString)("h_trk_pt"+trk_pt[ibin2]+"_reco_refl_"+cent[ibin-1]));
        else h_bkgtrkpt_data[ibin][ibin2] = (TH1D*)f_trk_PbPb_data->Get((TString)("h_trk_pt"+trk_pt[ibin2]+"_reco_refl"+cent[ibin-1]));
        h_trkk_data[ibin][ibin2] = (TH1D*)f_trk_PbPb_data->Get((TString)("h_trk_ptk"+kappa_str[ibin2]+"_reco_"+cent[ibin-1]));
        h_bkgtrkk_data[ibin][ibin2] = (TH1D*)f_trk_PbPb_data->Get((TString)("h_trk_ptk"+kappa_str[ibin2]+"_reco_refl_"+cent[ibin-1]));

        if(do_kappa_fitting){
            h_chg_corrpt_kappa[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
            h_chg_corrpt_q_kappa[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_q_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
            h_chg_corrpt_g_kappa[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_g_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
            
            h_chg_corrpt_up_kappa[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_up_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
            h_chg_corrpt_down_kappa[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_d_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));      
            h_chg_corrpt_upb_kappa[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_upbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
            h_chg_corrpt_downb_kappa[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_dbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

            h_chg_corrpt_c_kappa[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_c_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
            h_chg_corrpt_s_kappa[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_s_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
            h_chg_corrpt_b_kappa[ibin][ibin2] = (TH2F*)closure_histos_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_b_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

            h_chg_corrpt_data_kappa[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_chg_corrpt_kappa_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

            h_chg_corrpt_data_jer_kappa[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data_kappa_jer->Get((TString)("h_chg_corrpt_kappa_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));  
            h_chg_corrpt_data_trkpos_kappa[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data_kappa_trkuncpos->Get((TString)("h_chg_corrpt_kappa_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));  

            h_chg_data_jer_kappa[ibin][ibin2] = h_chg_corrpt_data_jer_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_data_kappa_jer_%d_%d",ibin,ibin2),h_chg_corrpt_data_jer_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_data_jer_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_data_jer_kappa[ibin][ibin2]->Scale(1./h_chg_data_jer_kappa[ibin][ibin2]->Integral());
            h_chg_data_trkpos_kappa[ibin][ibin2] = h_chg_corrpt_data_trkpos_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_data_trkpos_kappa_%d_%d",ibin,ibin2),h_chg_corrpt_data_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_data_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_data_trkpos_kappa[ibin][ibin2]->Scale(1./h_chg_data_trkpos_kappa[ibin][ibin2]->Integral());
        }
/*
        h_chg_corrpt_spillover[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data_spillover->Get((TString)("h_chg_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
        h_chg_data_spillover[ibin][ibin2] = h_chg_corrpt_spillover[ibin][ibin2]-> ProjectionY(Form("h_chg_data_spillover_%d_%d",ibin,ibin2),h_chg_corrpt_spillover[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_spillover[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_data_spillover[ibin][ibin2]->Scale(1./h_chg_data_spillover[ibin][ibin2]->Integral());
        h_chg_data_spillover[ibin][ibin2]->Rebin(10);

        h_chg_corrpt_spilloverdn[ibin][ibin2] = (TH2F*)closure_histos_PbPb_data_spillover_dn->Get((TString)("h_chg_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
        h_chg_data_spilloverdn[ibin][ibin2] = h_chg_corrpt_spilloverdn[ibin][ibin2]-> ProjectionY(Form("h_chg_data_spilloverdn_%d_%d",ibin,ibin2),h_chg_corrpt_spilloverdn[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_spilloverdn[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_data_spilloverdn[ibin][ibin2]->Scale(1./h_chg_data_spilloverdn[ibin][ibin2]->Integral());
        h_chg_data_spilloverdn[ibin][ibin2]->Rebin(10);     
*/
      }

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
      
      h_pt_reco_q[ibin][ibin2] = h_chg_corrpt_q[ibin][ibin2]-> ProjectionX(Form("h_pt_reco_q_%d_%d",ibin,ibin2),h_chg_corrpt_q[ibin][ibin2]->GetYaxis()->FindBin(-10.),h_chg_corrpt_q[ibin][ibin2]->GetYaxis()->FindBin(10.),"");
      h_pt_reco_g[ibin][ibin2] = h_chg_corrpt_g[ibin][ibin2]-> ProjectionX(Form("h_pt_reco_g_%d_%d",ibin,ibin2),h_chg_corrpt_g[ibin][ibin2]->GetYaxis()->FindBin(-10.),h_chg_corrpt_g[ibin][ibin2]->GetYaxis()->FindBin(10.),"");

      quark_reco_int[ibin][ibin2] = h_pt_reco_q[ibin][ibin2]->Integral(h_pt_reco_q[ibin][ibin2]->FindBin(120.),h_pt_reco_q[ibin][ibin2]->FindBin(600.));
      gluon_reco_int[ibin][ibin2] = h_pt_reco_g[ibin][ibin2]->Integral(h_pt_reco_g[ibin][ibin2]->FindBin(120.),h_pt_reco_g[ibin][ibin2]->FindBin(600.)); 

      h_chg_data_trkpos[ibin][ibin2] = h_chg_corrpt_data_trkpos[ibin][ibin2]-> ProjectionY(Form("h_chg_data_trkpos_%d_%d",ibin,ibin2),h_chg_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");

      for(int ibin3=0; ibin3<trkeffbins;ibin3++){        
          h_chg_data_trkeff[ibin][ibin2][ibin3] = h_chg_corrpt_data_trkeff[ibin][ibin2][ibin3]-> ProjectionY(Form("h_chg_data_trkeff_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt_data_trkeff[ibin][ibin2][ibin3]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_data_trkeff[ibin][ibin2][ibin3]->GetXaxis()->FindBin(pt_high),"");
          h_chg_data_trkeff[ibin][ibin2][ibin3]->Scale(1./(h_chg_data_trkeff[ibin][ibin2][ibin3]->Integral()));
          if(do_kappa_fitting){
            h_chg_data_kappa_trkeff[ibin][ibin2][ibin3] = h_chg_corrpt_data_kappa_trkeff[ibin][ibin2][ibin3]-> ProjectionY(Form("h_chg_data_kappa_trkeff_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt_data_trkeff[ibin][ibin2][ibin3]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_data_trkeff[ibin][ibin2][ibin3]->GetXaxis()->FindBin(pt_high),"");
            h_chg_data_kappa_trkeff[ibin][ibin2][ibin3]->Scale(1./(h_chg_data_kappa_trkeff[ibin][ibin2][ibin3]->Integral()));
          }
      }

      h_chg_data_jeup[ibin][ibin2] = h_chg_corrpt_data[ibin][ibin2]-> ProjectionY(Form("h_chg_data_jeup_%d_%d",ibin,ibin2),h_chg_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_low_up),h_chg_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
      h_chg_data_jedn[ibin][ibin2] = h_chg_corrpt_data[ibin][ibin2]-> ProjectionY(Form("h_chg_data_jedn_%d_%d",ibin,ibin2),h_chg_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_low_dn),h_chg_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");

      h_chg_reco_gubdb[ibin][ibin2] = (TH1D*)h_chg_reco_g[ibin][ibin2]->Clone(Form("h_chg_reco_gubdb_%d_%d",ibin,ibin2));
      h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_upb[ibin][ibin2]);
      h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_downb[ibin][ibin2]);
      h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_c[ibin][ibin2]);
      h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_s[ibin][ibin2]);
      h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_b[ibin][ibin2]);

      h_chg_reco_others[ibin][ibin2] = (TH1D*)h_chg_reco_upb[ibin][ibin2]->Clone(Form("h_chg_reco_others_%d_%d",ibin,ibin2));
      h_chg_reco_others[ibin][ibin2]->Add(h_chg_reco_downb[ibin][ibin2]);
      h_chg_reco_others[ibin][ibin2]->Add(h_chg_reco_c[ibin][ibin2]);
      h_chg_reco_others[ibin][ibin2]->Add(h_chg_reco_s[ibin][ibin2]);
      h_chg_reco_others[ibin][ibin2]->Add(h_chg_reco_b[ibin][ibin2]);

      //normalizing chg histos
      h_chg_reco_q[ibin][ibin2]->Scale(1./(h_chg_reco_q[ibin][ibin2]->Integral()));
      h_chg_reco_g[ibin][ibin2]->Scale(1./(h_chg_reco_g[ibin][ibin2]->Integral()));
      h_chg_reco_u[ibin][ibin2]->Scale(1./(h_chg_reco_u[ibin][ibin2]->Integral()));                        
      h_chg_reco_upb[ibin][ibin2]->Scale(1./(h_chg_reco_upb[ibin][ibin2]->Integral()));        
      h_chg_reco_downb[ibin][ibin2]->Scale(1./(h_chg_reco_downb[ibin][ibin2]->Integral()));
      h_chg_reco_gubdb[ibin][ibin2]->Scale(1./(h_chg_reco_gubdb[ibin][ibin2]->Integral()));
      h_chg_reco_others[ibin][ibin2]->Scale(1./(h_chg_reco_others[ibin][ibin2]->Integral()));
      h_chg_data[ibin][ibin2]->Scale(1./(h_chg_data[ibin][ibin2]->Integral()));
      h_chg_reco[ibin][ibin2]->Scale(1./(h_chg_reco[ibin][ibin2]->Integral()));

      h_chg_data_trkpos[ibin][ibin2]->Scale(1./(h_chg_data_trkpos[ibin][ibin2]->Integral()));
      h_chg_data_jeup[ibin][ibin2]->Scale(1./(h_chg_data_jeup[ibin][ibin2]->Integral()));
      h_chg_data_jedn[ibin][ibin2]->Scale(1./(h_chg_data_jedn[ibin][ibin2]->Integral()));

      //for(int ibin4 = 1; ibin4 < h_chg_data[ibin][ibin2]->GetNbinsX()+1; ibin4++){
      //  h_chg_data_sys[ibin][ibin2]->SetBinContent(ibin4,h_chg_data[ibin][ibin2]->GetBinContent(ibin4));
      //}

      if(do_data){
        //scaling up and down to match data
        h_chg_reco_up[ibin][ibin2]->Scale(1./(h_chg_reco_up[ibin][ibin2]->Integral()));
        h_chg_reco_down[ibin][ibin2]->Scale(1./(h_chg_reco_down[ibin][ibin2]->Integral()));

        h_chg_reco_ud_scaled[ibin][ibin2] = (TH1D*)h_chg_reco_down[ibin][ibin2]->Clone((TString)("h_chg_reco_ud_scaled_"+cent[ibin]+"_"+trk[ibin2]));
        if(ibin==0) h_chg_reco_up[ibin][ibin2]->Scale(1.896);
        else h_chg_reco_up[ibin][ibin2]->Scale(0.89);
        h_chg_reco_ud_scaled[ibin][ibin2]->Add(h_chg_reco_up[ibin][ibin2]); 
        h_chg_reco_ud_scaled[ibin][ibin2]->Scale(1./h_chg_reco_ud_scaled[ibin][ibin2]->Integral());

        if(ibin==0) h_chg_reco_up[ibin][ibin2]->Scale(1./1.896);
        else h_chg_reco_up[ibin][ibin2]->Scale(1./0.89);
      }
      else{
        //MC 
        h_chg_reco_up[ibin][ibin2]->Scale(1./(h_chg_reco[ibin][ibin2]->Integral()));
        h_chg_reco_down[ibin][ibin2]->Scale(1./(h_chg_reco[ibin][ibin2]->Integral()));

        h_chg_reco_ud_scaled[ibin][ibin2] = (TH1D*)h_chg_reco_down[ibin][ibin2]->Clone((TString)("h_chg_reco_ud_scaled_"+cent[ibin]+"_"+trk[ibin2]));
        h_chg_reco_ud_scaled[ibin][ibin2]->Add(h_chg_reco_up[ibin][ibin2]); 
        h_chg_reco_ud_scaled[ibin][ibin2]->Scale(1./h_chg_reco_ud_scaled[ibin][ibin2]->Integral());

        h_chg_reco_up[ibin][ibin2]->Scale(1./(h_chg_reco_up[ibin][ibin2]->Integral()));
        h_chg_reco_down[ibin][ibin2]->Scale(1./(h_chg_reco_down[ibin][ibin2]->Integral()));
      }

      //cout<<"up mean: "<<h_chg_reco_up[ibin][ibin2]->GetMean()<<endl;
      //cout<<"down mean: "<<h_chg_reco_down[ibin][ibin2]->GetMean()<<endl;
      //cout<<"gluon mean: "<<h_chg_reco_g[ibin][ibin2]->GetMean()<<endl;

      ///kappa projections
      if(do_kappa_fitting){
          h_chg_reco_kappa[ibin][ibin2] = h_chg_corrpt_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_kappa_%d_%d",ibin,ibin2),h_chg_corrpt_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
          h_chg_reco_q_kappa[ibin][ibin2] = h_chg_corrpt_q_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_kappa_q_%d_%d",ibin,ibin2),h_chg_corrpt_q_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_q_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
          h_chg_reco_g_kappa[ibin][ibin2] = h_chg_corrpt_g_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_kappa_g_%d_%d",ibin,ibin2),h_chg_corrpt_g_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_g_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
          h_chg_reco_up_kappa[ibin][ibin2] = h_chg_corrpt_up_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_kappa_up_%d_%d",ibin,ibin2),h_chg_corrpt_up_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_up_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
          h_chg_reco_down_kappa[ibin][ibin2] = h_chg_corrpt_down_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_kappa_down_%d_%d",ibin,ibin2),h_chg_corrpt_down_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_down_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
          h_chg_reco_upb_kappa[ibin][ibin2] = h_chg_corrpt_upb_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_kappa_upb_%d_%d",ibin,ibin2),h_chg_corrpt_upb_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_upb_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
          h_chg_reco_downb_kappa[ibin][ibin2] = h_chg_corrpt_downb_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_kappa_downb_%d_%d",ibin,ibin2),h_chg_corrpt_downb_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_downb_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
          h_chg_data_kappa[ibin][ibin2] = h_chg_corrpt_data_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_data_kappa_%d_%d",ibin,ibin2),h_chg_corrpt_data_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_data_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
          h_chg_reco_c_kappa[ibin][ibin2] = h_chg_corrpt_c_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_kappa_c_%d_%d",ibin,ibin2),h_chg_corrpt_c_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_c_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
          h_chg_reco_s_kappa[ibin][ibin2] = h_chg_corrpt_s_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_kappa_s_%d_%d",ibin,ibin2),h_chg_corrpt_s_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_s_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
          h_chg_reco_b_kappa[ibin][ibin2] = h_chg_corrpt_b_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_kappa_b_%d_%d",ibin,ibin2),h_chg_corrpt_b_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_b_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");

          h_chg_reco_gubdb_kappa[ibin][ibin2] = (TH1D*)h_chg_reco_g_kappa[ibin][ibin2]->Clone(Form("h_chg_reco_kappa_gubdb_%d_%d",ibin,ibin2));
          h_chg_reco_gubdb_kappa[ibin][ibin2]->Add(h_chg_reco_upb_kappa[ibin][ibin2]);
          h_chg_reco_gubdb_kappa[ibin][ibin2]->Add(h_chg_reco_downb_kappa[ibin][ibin2]);
          h_chg_reco_gubdb_kappa[ibin][ibin2]->Add(h_chg_reco_c_kappa[ibin][ibin2]);
          h_chg_reco_gubdb_kappa[ibin][ibin2]->Add(h_chg_reco_s_kappa[ibin][ibin2]);
          h_chg_reco_gubdb_kappa[ibin][ibin2]->Add(h_chg_reco_b_kappa[ibin][ibin2]);

          h_chg_reco_others_kappa[ibin][ibin2] = (TH1D*)h_chg_reco_upb_kappa[ibin][ibin2]->Clone(Form("h_chg_reco_kappa_others_%d_%d",ibin,ibin2));
          h_chg_reco_others_kappa[ibin][ibin2]->Add(h_chg_reco_downb_kappa[ibin][ibin2]);
          h_chg_reco_others_kappa[ibin][ibin2]->Add(h_chg_reco_c_kappa[ibin][ibin2]);
          h_chg_reco_others_kappa[ibin][ibin2]->Add(h_chg_reco_s_kappa[ibin][ibin2]);
          h_chg_reco_others_kappa[ibin][ibin2]->Add(h_chg_reco_b_kappa[ibin][ibin2]);

          //normalizing chg histos
          h_chg_reco_q_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_q_kappa[ibin][ibin2]->Integral()));
          h_chg_reco_g_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_g_kappa[ibin][ibin2]->Integral()));
          h_chg_reco_upb_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_upb_kappa[ibin][ibin2]->Integral()));        
          h_chg_reco_downb_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_downb_kappa[ibin][ibin2]->Integral()));
          h_chg_reco_gubdb_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_gubdb_kappa[ibin][ibin2]->Integral()));
          h_chg_reco_others_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_others_kappa[ibin][ibin2]->Integral()));
          h_chg_data_kappa[ibin][ibin2]->Scale(1./(h_chg_data_kappa[ibin][ibin2]->Integral()));
          h_chg_reco_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_kappa[ibin][ibin2]->Integral()));

          if(do_data){
            //scaling up and down to match data
            h_chg_reco_up_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_up_kappa[ibin][ibin2]->Integral()));
            h_chg_reco_down_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_down_kappa[ibin][ibin2]->Integral()));

            h_chg_reco_kappa_ud_scaled[ibin][ibin2] = (TH1D*)h_chg_reco_down_kappa[ibin][ibin2]->Clone((TString)("h_chg_reco_kappa_ud_scaled_"+cent[ibin]+"_"+trk[ibin2]));
            if(ibin==0) h_chg_reco_up_kappa[ibin][ibin2]->Scale(1.896);
            else h_chg_reco_up_kappa[ibin][ibin2]->Scale(0.89);
            h_chg_reco_kappa_ud_scaled[ibin][ibin2]->Add(h_chg_reco_up_kappa[ibin][ibin2]); 
            h_chg_reco_kappa_ud_scaled[ibin][ibin2]->Scale(1./h_chg_reco_kappa_ud_scaled[ibin][ibin2]->Integral());

            if(ibin==0) h_chg_reco_up_kappa[ibin][ibin2]->Scale(1./1.896);
            else h_chg_reco_up_kappa[ibin][ibin2]->Scale(1./0.89);
          }
          else{
            //MC
            h_chg_reco_up_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_kappa[ibin][ibin2]->Integral()));
            h_chg_reco_down_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_kappa[ibin][ibin2]->Integral()));

            h_chg_reco_kappa_ud_scaled[ibin][ibin2] = (TH1D*)h_chg_reco_down_kappa[ibin][ibin2]->Clone((TString)("h_chg_reco_kappa_ud_scaled_"+cent[ibin]+"_"+trk[ibin2]));
            h_chg_reco_kappa_ud_scaled[ibin][ibin2]->Add(h_chg_reco_up_kappa[ibin][ibin2]); 
            h_chg_reco_kappa_ud_scaled[ibin][ibin2]->Scale(1./h_chg_reco_kappa_ud_scaled[ibin][ibin2]->Integral());

            h_chg_reco_up_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_up_kappa[ibin][ibin2]->Integral()));
            h_chg_reco_down_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_down_kappa[ibin][ibin2]->Integral()));
          }
      }

      ///////////////////width of jet chg distributions

      for(int ibin3=0; ibin3<njtptbins; ibin3++){
        
        if(do_ptcut_fitting){
          //h_chg_data_avg[ibin][ibin3][ibin2] = h_chg_corrpt_data[ibin][ibin2]-> ProjectionY(Form("h_chg_data_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_avg[ibin][ibin3][ibin2] = h_chg_corrpt[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_q_avg[ibin][ibin3][ibin2] = h_chg_corrpt_q[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_q_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_g_avg[ibin][ibin3][ibin2] = h_chg_corrpt_g[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_g_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_up_avg[ibin][ibin3][ibin2] = h_chg_corrpt_up[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_up_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_down_avg[ibin][ibin3][ibin2] = h_chg_corrpt_down[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_down_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
        }
        else if(do_kappa_fitting){
          //h_chg_data_avg[ibin][ibin3][ibin2] = h_chg_corrpt_data_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_data_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_avg[ibin][ibin3][ibin2] = h_chg_corrpt_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_q_avg[ibin][ibin3][ibin2] = h_chg_corrpt_q_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_q_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_g_avg[ibin][ibin3][ibin2] = h_chg_corrpt_g_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_g_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_up_avg[ibin][ibin3][ibin2] = h_chg_corrpt_up_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_up_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");
          h_chg_reco_down_avg[ibin][ibin3][ibin2] = h_chg_corrpt_down_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_down_avg_%d_%d_%d",ibin,ibin2,ibin3),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3]),h_chg_corrpt[ibin][ibin2]->GetXaxis()->FindBin(jtpt_bound[ibin3+1]),"");          
        }
        
        //h_chg_data_avg[ibin][ibin3][ibin2]->Scale(1./h_chg_data_avg[ibin][ibin3][ibin2]->Integral()); 
        h_chg_reco_avg[ibin][ibin3][ibin2]->Scale(1./h_chg_reco_avg[ibin][ibin3][ibin2]->Integral()); 
        h_chg_reco_q_avg[ibin][ibin3][ibin2]->Scale(1./h_chg_reco_q_avg[ibin][ibin3][ibin2]->Integral()); 
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->Scale(1./h_chg_reco_g_avg[ibin][ibin3][ibin2]->Integral()); 
        h_chg_reco_up_avg[ibin][ibin3][ibin2]->Scale(1./h_chg_reco_up_avg[ibin][ibin3][ibin2]->Integral()); 
        h_chg_reco_down_avg[ibin][ibin3][ibin2]->Scale(1./h_chg_reco_down_avg[ibin][ibin3][ibin2]->Integral()); 

        f_chg_data_avg[ibin][ibin3][ibin2] = new TF1(Form("f_chg_data_%d_%d_%d",ibin,ibin2,ibin3),"gaus",-2.,2.);
        f_chg_reco_avg[ibin][ibin3][ibin2] = new TF1(Form("f_chg_%d_%d_%d",ibin,ibin2,ibin3),"gaus",-2.,2.);
        f_chg_reco_q_avg[ibin][ibin3][ibin2] = new TF1(Form("f_chg_q_%d_%d_%d",ibin,ibin2,ibin3),"gaus",-2.,2.);
        f_chg_reco_g_avg[ibin][ibin3][ibin2] = new TF1(Form("f_chg_g_%d_%d_%d",ibin,ibin2,ibin3),"gaus",-2.,2.);
        f_chg_reco_up_avg[ibin][ibin3][ibin2] = new TF1(Form("f_chg_up_%d_%d_%d",ibin,ibin2,ibin3),"gaus",-2.,2.);
        f_chg_reco_up_avg[ibin][ibin3][ibin2]->SetLineColor(kBlue);
        f_chg_reco_down_avg[ibin][ibin3][ibin2] = new TF1(Form("f_chg_down_%d_%d_%d",ibin,ibin2,ibin3),"gaus",-2.,2.);
        f_chg_reco_down_avg[ibin][ibin3][ibin2]->SetLineColor(kGreen-2);
      }

      if(do_gen){
        //gen jet charge projections
        h_chg_ref[ibin][ibin2] = h_chg_refpt[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_%d_%d",ibin,ibin2),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_g[ibin][ibin2] = h_chg_refpt_g[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_g_%d_%d",ibin,ibin2),h_chg_refpt_g[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_g[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_q[ibin][ibin2] = h_chg_refpt_q[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_q_%d_%d",ibin,ibin2),h_chg_refpt_q[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_q[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_up[ibin][ibin2] = h_chg_refpt_up[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_up_%d_%d",ibin,ibin2),h_chg_refpt_up[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_up[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_down[ibin][ibin2] = h_chg_refpt_down[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_down_%d_%d",ibin,ibin2),h_chg_refpt_down[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_down[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_upb[ibin][ibin2] = h_chg_refpt_upb[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_upb_%d_%d",ibin,ibin2),h_chg_refpt_upb[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_upb[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
        h_chg_ref_downb[ibin][ibin2] = h_chg_refpt_downb[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_downb_%d_%d",ibin,ibin2),h_chg_refpt_downb[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_downb[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");

        h_pt_ref_q[ibin][ibin2] = h_chg_refpt_q[ibin][ibin2]-> ProjectionX(Form("h_pt_ref_q_%d_%d",ibin,ibin2),h_chg_corrpt_q[ibin][ibin2]->GetYaxis()->FindBin(-10.),h_chg_corrpt_q[ibin][ibin2]->GetYaxis()->FindBin(10.),"");
        h_pt_ref_g[ibin][ibin2] = h_chg_refpt_g[ibin][ibin2]-> ProjectionX(Form("h_pt_ref_g_%d_%d",ibin,ibin2),h_chg_corrpt_g[ibin][ibin2]->GetYaxis()->FindBin(-10.),h_chg_corrpt_g[ibin][ibin2]->GetYaxis()->FindBin(10.),"");

        quark_gen_int[ibin][ibin2] = h_pt_ref_q[ibin][ibin2]->Integral(h_pt_ref_q[ibin][ibin2]->FindBin(120.),h_pt_ref_q[ibin][ibin2]->FindBin(600.));
        gluon_gen_int[ibin][ibin2] = h_pt_ref_g[ibin][ibin2]->Integral(h_pt_ref_g[ibin][ibin2]->FindBin(120.),h_pt_ref_g[ibin][ibin2]->FindBin(600.)); 

        h_chg_ref_bkg[ibin][ibin2] = h_chg_refpt_bkg[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_bkg_%d_%d",ibin,ibin2),h_chg_refpt_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");

        h_chg_ref_gubdb[ibin][ibin2] = (TH1D*)h_chg_ref_g[ibin][ibin2]->Clone(Form("h_chg_ref_gubdb_%d_%d",ibin,ibin2));
        h_chg_ref_gubdb[ibin][ibin2]->Add(h_chg_ref_upb[ibin][ibin2]);
        h_chg_ref_gubdb[ibin][ibin2]->Add(h_chg_ref_downb[ibin][ibin2]);

        h_chg_ref[ibin][ibin2]->Scale(1./(h_chg_ref[ibin][ibin2]->Integral()));
        h_chg_ref_g[ibin][ibin2]->Scale(1./(h_chg_ref_g[ibin][ibin2]->Integral()));
        h_chg_ref_up[ibin][ibin2]->Scale(1./(h_chg_ref_up[ibin][ibin2]->Integral()));
        h_chg_ref_down[ibin][ibin2]->Scale(1./(h_chg_ref_down[ibin][ibin2]->Integral()));
        h_chg_ref_gubdb[ibin][ibin2]->Scale(1./(h_chg_ref_gubdb[ibin][ibin2]->Integral()));
        h_chg_ref_bkg[ibin][ibin2]->Scale(1./(h_chg_ref_bkg[ibin][ibin2]->Integral()));
      }

      if(draw_tracking){
          h_trk_data[ibin][ibin2] = h_trk_corrpt_data[ibin][ibin2]-> ProjectionY(Form("h_trk_data_%d_%d",ibin,ibin2),h_trk_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_trk_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
          h_trk_MC[ibin][ibin2] = h_trk_corrpt_MC[ibin][ibin2]-> ProjectionY(Form("h_trk_MC_%d_%d",ibin,ibin2),h_trk_corrpt_MC[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_trk_corrpt_MC[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
          h_bkgtrk_data[ibin][ibin2] = h_bkgtrk_corrpt_data[ibin][ibin2]-> ProjectionY(Form("h_bkgtrk_data_%d_%d",ibin,ibin2),h_trk_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_trk_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
          h_bkgtrk_MC[ibin][ibin2] = h_bkgtrk_corrpt_MC[ibin][ibin2]-> ProjectionY(Form("h_bkgtrk_MC_%d_%d",ibin,ibin2),h_trk_corrpt_MC[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_trk_corrpt_MC[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
      }

      h_chg_data_avg[ibin][ibin2] = h_chg_corrpt_data[ibin][ibin2]-> ProfileX(Form("h_chg_data_avg_%d_%d",ibin,ibin2),h_chg_corrpt_data[ibin][ibin2]->GetYaxis()->FindBin(-2.),h_chg_corrpt_data[ibin][ibin2]->GetYaxis()->FindBin(2.),"");

      h_chg_reco_bkg_avg[ibin][ibin2] = h_chg_corrpt_bkg[ibin][ibin2]-> ProfileX(Form("h_chg_bkg_reco_avg_%d_%d",ibin,ibin2),h_chg_corrpt_bkg[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_bkg[ibin][ibin2]->GetYaxis()->FindBin(10),"");
      h_chg_data_bkg_avg[ibin][ibin2] = h_chg_corrpt_data_bkg[ibin][ibin2]-> ProfileX(Form("h_chg_bkg_data_avg_%d_%d",ibin,ibin2),h_chg_corrpt_data_bkg[ibin][ibin2]->GetYaxis()->FindBin(-10),h_chg_corrpt_data_bkg[ibin][ibin2]->GetYaxis()->FindBin(10),"");
      
      h_chg_reco_bkg[ibin][ibin2] = h_chg_corrpt_bkg[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_bkg_%d_%d",ibin,ibin2),h_chg_corrpt_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
      h_chg_reco_bkg[ibin][ibin2]->Scale(1./(h_chg_reco_bkg[ibin][ibin2]->Integral()));
      h_chg_data_bkg[ibin][ibin2] = h_chg_corrpt_data_bkg[ibin][ibin2]-> ProjectionY(Form("h_chg_data_bkg_%d_%d",ibin,ibin2),h_chg_corrpt_data_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_data_bkg[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
      h_chg_data_bkg[ibin][ibin2]->Scale(1./(h_chg_data_bkg[ibin][ibin2]->Integral()));
  
      h_reco_corr_0trks[ibin][ibin2] = (TH1D*)h_reco_MC[ibin]->Clone(Form("h_reco_corr_0trk_%d_%d",ibin,ibin2));  
      h_reco_corr_0trks[ibin][ibin2]->Add(h_reco_corr_no0trks[ibin][ibin2],-1);
      h_reco_corr_0trks_q[ibin][ibin2] = (TH1D*)h_reco_MC_q[ibin]->Clone(Form("h_reco_corr_0trk_q_%d_%d",ibin,ibin2));  
      h_reco_corr_0trks_q[ibin][ibin2]->Add(h_reco_corr_no0trks_q[ibin][ibin2],-1);
      h_reco_corr_0trks_g[ibin][ibin2] = (TH1D*)h_reco_MC_g[ibin]->Clone(Form("h_reco_corr_0trk_g_%d_%d",ibin,ibin2));  
      h_reco_corr_0trks_g[ibin][ibin2]->Add(h_reco_corr_no0trks_g[ibin][ibin2],-1);

      h_gen_full_0trks[ibin][ibin2] = (TH1D*)h_ref_MC[ibin]->Clone(Form("h_gen_full_0trk_%d_%d",ibin,ibin2));  
      h_gen_full_0trks[ibin][ibin2]->Add(h_gen_full_no0trks[ibin][ibin2],-1);
      h_gen_full_0trks_q[ibin][ibin2] = (TH1D*)h_ref_MC_q[ibin]->Clone(Form("h_gen_full_0trk_q_%d_%d",ibin,ibin2));  
      h_gen_full_0trks_q[ibin][ibin2]->Add(h_gen_full_no0trks_q[ibin][ibin2],-1);
      h_gen_full_0trks_g[ibin][ibin2] = (TH1D*)h_ref_MC_g[ibin]->Clone(Form("h_gen_full_0trk_g_%d_%d",ibin,ibin2));  
      h_gen_full_0trks_g[ibin][ibin2]->Add(h_gen_full_no0trks_g[ibin][ibin2],-1);

      h_reco_data_0trks[ibin][ibin2] = (TH1D*)h_reco_data[ibin]->Clone(Form("h_reco_data_0trk_%d_%d",ibin,ibin2));  
      h_reco_data_0trks[ibin][ibin2]->Add(h_reco_data_no0trks[ibin][ibin2],-1);
    }
  }

  for(int ibin=0; ibin<nCBins; ibin++){

    sprintf(saythis,"h_0trks_data_cent%d",ibin);
    h_0trks_data[ibin] = new TH1D(saythis,"",4,2.,6.);
    h_0trks_data[ibin]->Sumw2(); 

    sprintf(saythis,"h_0trks_MC_cent%d",ibin);
    h_0trks_MC[ibin] = new TH1D(saythis,"",4,2.,6.);
    h_0trks_MC[ibin]->Sumw2(); 

    sprintf(saythis,"h_corr_factor_cent%d",ibin);
    h_corr_factor[ibin] = new TH1D(saythis,"",4,2.,6.);
    h_corr_factor[ibin]->Sumw2();

    sprintf(saythis,"h_syst_cent%d",ibin);
    h_syst[ibin] = new TH1D(saythis,"",3,kappa_bounds);
    h_syst[ibin]->Sumw2();

    sprintf(saythis,"h_syst_JER_cent%d",ibin);
    h_syst_JER[ibin] = new TH1D(saythis,"",3,kappa_bounds);
    h_syst_JER[ibin]->Sumw2();

    sprintf(saythis,"h_syst_trk_cent%d",ibin);
    h_syst_trk[ibin] = new TH1D(saythis,"",3,kappa_bounds);
    h_syst_trk[ibin]->Sumw2();

    sprintf(saythis,"h_syst_posneg_cent%d",ibin);
    h_syst_posneg[ibin] = new TH1D(saythis,"",3,kappa_bounds);
    h_syst_posneg[ibin]->Sumw2();

    sprintf(saythis,"h_syst_0trk_cent%d",ibin);
    h_syst_0trk[ibin] = new TH1D(saythis,"",3,kappa_bounds);
    h_syst_0trk[ibin]->Sumw2();

    for(int ibin2=0; ibin2<ntrkbins; ibin2++){
      h_0trks_data[ibin]->SetBinContent(ibin2+1,h_reco_data_0trks[ibin][ibin2]->Integral(h_reco_data_0trks[ibin][ibin2]->FindBin(pt_low),h_reco_data_0trks[ibin][ibin2]->FindBin(pt_high))/h_reco_data[ibin]->Integral(h_reco_data_0trks[ibin][ibin2]->FindBin(pt_low),h_reco_data_0trks[ibin][ibin2]->FindBin(pt_high)));  
      h_0trks_data[ibin]->SetBinError(ibin2+1,0.001);  
      h_0trks_MC[ibin]->SetBinContent(ibin2+1,h_reco_corr_0trks[ibin][ibin2]->Integral(h_reco_data_0trks[ibin][ibin2]->FindBin(pt_low),h_reco_data_0trks[ibin][ibin2]->FindBin(pt_high))/h_reco_MC[ibin]->Integral(h_reco_data_0trks[ibin][ibin2]->FindBin(pt_low),h_reco_data_0trks[ibin][ibin2]->FindBin(pt_high)));  
      h_0trks_MC[ibin]->SetBinError(ibin2+1,0.001);  
      diff_0trk_data_MC[ibin][ibin2] = h_0trks_data[ibin]->GetBinContent(ibin2+1) - h_0trks_MC[ibin]->GetBinContent(ibin2+1);
      //cout<<ibin<<". "<<ibin2<<" "<<diff_0trk_data_MC[ibin][ibin2]<<endl;
    }
  }

  //replicas for plotting 
  for (int ibin=0;ibin<nCBins;ibin++){
    for (int ibin2=0;ibin2<ntrkbins;ibin2++){
      sprintf(saythis,"h_chg_data_cop_cent%d_trk%d",ibin,ibin2);
      h_chg_data_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_data_cop[ibin][ibin2]->Sumw2();
      sprintf(saythis,"h_chg_reco_cop_cent%d_trk%d",ibin,ibin2);
      h_chg_reco_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_reco_cop[ibin][ibin2]->Sumw2();
      sprintf(saythis,"h_chg_reco_up_cop_cent%d_trk%d",ibin,ibin2);
      h_chg_reco_up_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_reco_up_cop[ibin][ibin2]->Sumw2();
      sprintf(saythis,"h_chg_reco_down_cop_cent%d_trk%d",ibin,ibin2);
      h_chg_reco_down_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_reco_down_cop[ibin][ibin2]->Sumw2();
      sprintf(saythis,"h_chg_reco_gubdb_cop_cent%d_trk%d",ibin,ibin2);
      h_chg_reco_gubdb_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_reco_gubdb_cop[ibin][ibin2]->Sumw2();
      sprintf(saythis,"h_chg_reco_g_cop_cent%d_trk%d",ibin,ibin2);
      h_chg_reco_g_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_reco_g_cop[ibin][ibin2]->Sumw2();
      sprintf(saythis,"h_chg_reco_others_cop_cent%d_trk%d",ibin,ibin2);
      h_chg_reco_others_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_reco_others_cop[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_data_cop_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_data_cop_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_data_cop_kappa[ibin][ibin2]->Sumw2();
      sprintf(saythis,"h_chg_reco_cop_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_reco_cop_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_reco_cop_kappa[ibin][ibin2]->Sumw2();
      sprintf(saythis,"h_chg_reco_up_cop_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_reco_up_cop_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_reco_up_cop_kappa[ibin][ibin2]->Sumw2();
      sprintf(saythis,"h_chg_reco_down_cop_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_reco_down_cop_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_reco_down_cop_kappa[ibin][ibin2]->Sumw2();
      sprintf(saythis,"h_chg_reco_gubdb_cop_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_reco_gubdb_cop_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_reco_gubdb_cop_kappa[ibin][ibin2]->Sumw2();
      sprintf(saythis,"h_chg_reco_g_cop_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_reco_g_cop_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_reco_g_cop_kappa[ibin][ibin2]->Sumw2();
      sprintf(saythis,"h_chg_reco_others_cop_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_reco_others_cop_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_reco_others_cop_kappa[ibin][ibin2]->Sumw2();
    }
  }

  for(int ibin2=0;ibin2<ntrkbins;ibin2++){
    for(int ibin=0;ibin<nCBins;ibin++){
      for(int ibin3=1;ibin3<(h_chg_data[ibin][ibin2]->GetXaxis()->GetNbins())+1;ibin3++){
        if(do_ptcut_fitting){
            h_chg_data_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_data[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_data_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_data[ibin][ibin2]->GetBinError(ibin3));
            h_chg_reco_up_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_up[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_reco_up_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_up[ibin][ibin2]->GetBinError(ibin3));
            h_chg_reco_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_reco_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_reco[ibin][ibin2]->GetBinError(ibin3));
            h_chg_reco_down_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_down[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_reco_down_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_down[ibin][ibin2]->GetBinError(ibin3));
            h_chg_reco_g_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_g[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_reco_g_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_g[ibin][ibin2]->GetBinError(ibin3));
            h_chg_reco_gubdb_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_gubdb[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_reco_gubdb_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_gubdb[ibin][ibin2]->GetBinError(ibin3));
            h_chg_reco_others_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_others[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_reco_others_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_others[ibin][ibin2]->GetBinError(ibin3));
        }

        if(do_kappa_fitting){
            h_chg_data_cop_kappa[ibin][ibin2]->SetBinContent(ibin3,h_chg_data_kappa[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_data_cop_kappa[ibin][ibin2]->SetBinError(ibin3,h_chg_data_kappa[ibin][ibin2]->GetBinError(ibin3));
            h_chg_reco_up_cop_kappa[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_up_kappa[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_reco_up_cop_kappa[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_up_kappa[ibin][ibin2]->GetBinError(ibin3));
            h_chg_reco_cop_kappa[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_kappa[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_reco_cop_kappa[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_kappa[ibin][ibin2]->GetBinError(ibin3));
            h_chg_reco_down_cop_kappa[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_down_kappa[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_reco_down_cop_kappa[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_down_kappa[ibin][ibin2]->GetBinError(ibin3));
            h_chg_reco_g_cop_kappa[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_g_kappa[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_reco_g_cop_kappa[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_g_kappa[ibin][ibin2]->GetBinError(ibin3));
            h_chg_reco_gubdb_cop_kappa[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_gubdb_kappa[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_reco_gubdb_cop_kappa[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_gubdb_kappa[ibin][ibin2]->GetBinError(ibin3));
            h_chg_reco_others_cop_kappa[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_others_kappa[ibin][ibin2]->GetBinContent(ibin3));
            h_chg_reco_others_cop_kappa[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_others_kappa[ibin][ibin2]->GetBinError(ibin3));      
        }
      }
      if(do_ptcut_fitting){
          h_chg_data_cop[ibin][ibin2]->SetLineColor(kBlack); h_chg_data_cop[ibin][ibin2]->SetFillColor(kGray); h_chg_data_cop[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_data_cop[ibin][ibin2]->SetMarkerStyle(20);
          h_chg_reco_cop[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_cop[ibin][ibin2]->SetFillColor(kGray); h_chg_reco_cop[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_reco_cop[ibin][ibin2]->SetMarkerStyle(20);
          h_chg_reco_up_cop[ibin][ibin2]->SetLineColor(kRed+2); h_chg_reco_up_cop[ibin][ibin2]->SetFillColor(kRed+2); h_chg_reco_up_cop[ibin][ibin2]->SetMarkerColor(kRed+2); h_chg_reco_up_cop[ibin][ibin2]->SetMarkerStyle(0);    
          h_chg_reco_down_cop[ibin][ibin2]->SetLineColor(kYellow-7); h_chg_reco_down_cop[ibin][ibin2]->SetFillColorAlpha(kYellow-7,0.7); h_chg_reco_down_cop[ibin][ibin2]->SetMarkerColor(kYellow-7); h_chg_reco_down_cop[ibin][ibin2]->SetMarkerStyle(0); 
          h_chg_reco_g_cop[ibin][ibin2]->SetLineColor(kAzure+1); h_chg_reco_g_cop[ibin][ibin2]->SetFillColorAlpha(kAzure+1,0.7); h_chg_reco_g_cop[ibin][ibin2]->SetMarkerColor(kAzure+1); h_chg_reco_g_cop[ibin][ibin2]->SetMarkerStyle(0);      
          h_chg_reco_others_cop[ibin][ibin2]->SetLineColor(kGreen-2); h_chg_reco_others_cop[ibin][ibin2]->SetFillColorAlpha(kGreen-2,0.7); h_chg_reco_others_cop[ibin][ibin2]->SetMarkerColor(kGreen-2); h_chg_reco_others_cop[ibin][ibin2]->SetMarkerStyle(0);

          h_chg_reco_gubdb_cop[ibin][ibin2]->Rebin(20);
          h_chg_reco_g_cop[ibin][ibin2]->Rebin(20);
          h_chg_reco_others_cop[ibin][ibin2]->Rebin(20);
          h_chg_reco_down_cop[ibin][ibin2]->Rebin(20);
          h_chg_reco_cop[ibin][ibin2]->Rebin(20);
          h_chg_data_cop[ibin][ibin2]->Rebin(20);
          h_chg_reco_up_cop[ibin][ibin2]->Rebin(20);
      }
      if(do_kappa_fitting){
          h_chg_data_cop_kappa[ibin][ibin2]->SetLineColor(kBlack); h_chg_data_cop_kappa[ibin][ibin2]->SetFillColor(kGray); h_chg_data_cop_kappa[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_data_cop_kappa[ibin][ibin2]->SetMarkerStyle(2);
          h_chg_reco_up_cop_kappa[ibin][ibin2]->SetLineColor(kRed+2); h_chg_reco_up_cop_kappa[ibin][ibin2]->SetFillColor(kRed+2); h_chg_reco_up_cop_kappa[ibin][ibin2]->SetMarkerColor(kRed+2); h_chg_reco_up_cop_kappa[ibin][ibin2]->SetMarkerStyle(0);    
          h_chg_reco_down_cop_kappa[ibin][ibin2]->SetLineColor(kYellow-7); h_chg_reco_down_cop_kappa[ibin][ibin2]->SetFillColorAlpha(kYellow-7,0.7); h_chg_reco_down_cop_kappa[ibin][ibin2]->SetMarkerColor(kYellow-7); h_chg_reco_down_cop_kappa[ibin][ibin2]->SetMarkerStyle(0); 
          h_chg_reco_g_cop_kappa[ibin][ibin2]->SetLineColor(kAzure+1); h_chg_reco_g_cop_kappa[ibin][ibin2]->SetFillColorAlpha(kAzure+1,0.7); h_chg_reco_g_cop_kappa[ibin][ibin2]->SetMarkerColor(kAzure+1); h_chg_reco_g_cop_kappa[ibin][ibin2]->SetMarkerStyle(0);      
          h_chg_reco_others_cop_kappa[ibin][ibin2]->SetLineColor(kGreen-2); h_chg_reco_others_cop_kappa[ibin][ibin2]->SetFillColorAlpha(kGreen-2,0.7); h_chg_reco_others_cop_kappa[ibin][ibin2]->SetMarkerColor(kGreen-2); h_chg_reco_others_cop_kappa[ibin][ibin2]->SetMarkerStyle(0);
      }
    }
  }
}

TPad *pad_titlex;
TPad *pad_labelx;
TPad *padx;
TPad *pady;
TPad *pad_titlec;
TPad *pad_labelc;
TPad *pad1c;
TPad *pad2c;
TPad *pad_titleg;
TPad *pad_labelg;
TPad *pada;

TLine *tl;

void drawline(double x1, double y1, double x2, double y2){
    tl = new TLine(x1,y1,x2,y2);  
    tl->SetLineStyle(2);
    tl->Draw("same");
}

void drawPad(){

    pad_titlex = new TPad("pad_titlex", "",0.,0.94,1.,1.);
    pad_labelx = new TPad("pad_labelx", "",0.,0.,0.02,0.95);
    padx = new TPad("padx", "",0.02,0.4,1.,0.98);
    padx->Divide(5,1,0);  
    padx->Draw();
    pady = new TPad("pady", "",0.02,0.,1.,0.4);
    pady->Divide(5,1,0);  
    pady->Draw();
    pad_titlex->Draw();
    pad_labelx->Draw();

}

void drawgraphPad(){

  pad_titleg = new TPad("pad_titleg", "",0.,0.93,1.,1.);
  pad_labelg = new TPad("pad_labelg", "",0.,0.,0.02,0.95);
  pada = new TPad("pada", "",0.02,0.,1.,0.98);
  pada->Divide(5,1,0);  
  pada->Draw();
  pad_titleg->Draw();
  pad_labelg->Draw();

}

void drawchgPad(){
  pad_titlec = new TPad("pad_titlec", "",0.,0.94,1.,1.);
  pad_labelc = new TPad("pad_labelc", "",0.,0.,0.02,0.95);
  pad1c = new TPad("pad1c", "",0.02,0.3,1.,0.98);
  pad1c->Divide(5,1,0);  
  pad1c->Draw();
  pad2c = new TPad("pad2c", "",0.02,0.,1.,0.3);
  pad2c->Divide(5,1,0);  
  pad2c->Draw();
  pad_titlec->Draw();
  pad_labelc->Draw();
}

void h_cosm(TH1D* h){
  h->GetYaxis()->SetNdivisions(505);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.7);
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetNdivisions(505);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->CenterTitle();
}

void plot_charge_frac(bool do_data = false){

  get_histos(do_data);

  TLine *tl_q[nCBins]; TLine *tl_g[nCBins]; TLine *tl_u[nCBins]; TLine *tl_d[nCBins];
  TLine *tl_q_kappa[nCBins]; TLine *tl_g_kappa[nCBins];
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

      data_int[ibin] = h_reco_data[ibin]->Integral(low_bin,high_bin);

      ubdbcsb[ibin] = upbar[ibin]+downbar[ibin]+charm[ibin]+strange[ibin]+bottom[ibin];
      gluon_ubdbcsb[ibin] = gluon[ibin]/(gluon[ibin]+ubdbcsb[ibin]);
      u_ud[ibin] = up[ibin]/(up[ibin]+down[ibin]);
      d_ud[ibin] = down[ibin]/(up[ibin]+down[ibin]);

      cout<<ibin<<endl;

      cout<<"up: "<<up[ibin]<<"  down:  "<<down[ibin]<<"  g_ub_db:  "<<gluon_ubar_dbar[ibin]<<endl;
      cout<<"charm: "<<charm[ibin]<<"  strange:  "<<strange[ibin]<<"  bottom:  "<<bottom[ibin]<<endl;

      cout<<"total: "<<gluon_ubar_dbar[ibin]+up[ibin]+down[ibin]<<endl;

      cout<<"gluongen / gluonreco"<<gluon_gen[ibin]/gluon[ibin]<<".  quarkgen / quarkreco"<<quark_gen[ibin]/quark[ibin]<<endl;

      for(int ibin2=0; ibin2<ntrkbins; ibin2++){

        gluon_no0trk_gen_int[ibin][ibin2] = h_gen_full_no0trks_g[ibin][ibin2]->Integral(low_bin,high_bin)/h_gen_full_no0trks[ibin][ibin2]->Integral(low_bin,high_bin);
        quark_no0trk_gen_int[ibin][ibin2] = h_gen_full_no0trks_q[ibin][ibin2]->Integral(low_bin,high_bin)/h_gen_full_no0trks[ibin][ibin2]->Integral(low_bin,high_bin);

        gluon_no0trk_int[ibin][ibin2] = h_reco_corr_no0trks_g[ibin][ibin2]->Integral(low_bin,high_bin)/h_reco_corr_no0trks[ibin][ibin2]->Integral(low_bin,high_bin);
        quark_no0trk_int[ibin][ibin2] = h_reco_corr_no0trks_q[ibin][ibin2]->Integral(low_bin,high_bin)/h_reco_corr_no0trks[ibin][ibin2]->Integral(low_bin,high_bin);

        //quark_0trk_int[ibin][ibin2] = h_reco_corr_0trks_q[ibin][ibin2]->Integral(low_bin,high_bin)/h_reco_corr_0trks[ibin][ibin2]->Integral(low_bin,high_bin);
        quark_0trk_int[ibin][ibin2] = h_gen_full_0trks_q[ibin][ibin2]->Integral(low_bin,high_bin)/h_gen_full_0trks[ibin][ibin2]->Integral(low_bin,high_bin);
        gluon_0trk_int[ibin][ibin2] = h_gen_full_0trks_g[ibin][ibin2]->Integral(low_bin,high_bin)/h_gen_full_0trks[ibin][ibin2]->Integral(low_bin,high_bin);

        data_0trk_int[ibin][ibin2] = h_reco_data_0trks[ibin][ibin2]->Integral(low_bin,high_bin);
        data_no0trk_int[ibin][ibin2] = h_reco_data_no0trks[ibin][ibin2]->Integral(low_bin,high_bin);
        quark_data_0trk_int[ibin][ibin2] = quark_0trk_int[ibin][ibin2]*data_0trk_int[ibin][ibin2];
        gluon_data_0trk_int[ibin][ibin2] = (1-quark_0trk_int[ibin][ibin2])*data_0trk_int[ibin][ibin2];

        //gluon_corr_factor[ibin][ibin2] = gluon_no0trk_gen_int[ibin][ibin2]/gluon_no0trk_int[ibin][ibin2];
        //quark_corr_factor[ibin][ibin2] = quark_no0trk_gen_int[ibin][ibin2]/quark_no0trk_int[ibin][ibin2];

        gluon_corr_factor[ibin][ibin2] = gluon_gen[ibin]/gluon[ibin];
        quark_corr_factor[ibin][ibin2] = quark_gen[ibin]/quark[ibin];
        
        h_corr_factor[ibin]->SetBinContent(ibin2+1,quark_corr_factor[ibin][ibin2]);  
        h_corr_factor[ibin]->SetBinError(ibin2+1,0.001);  
      }

      tl_u[ibin] = new TLine(1.,up[ibin],6.,up[ibin]); tl_u[ibin]->SetLineStyle(2); tl_u[ibin]->SetLineColor(kBlue);
      tl_d[ibin] = new TLine(1.,down[ibin],6.,down[ibin]); tl_d[ibin]->SetLineStyle(9); tl_d[ibin]->SetLineColor(kBlue);
      //tl_q[ibin] = new TLine(1.,quark[ibin]*quark_corr_factor[ibin][0],6.,quark[ibin]*quark_corr_factor[ibin][0]); tl_q[ibin]->SetLineStyle(1); tl_q[ibin]->SetLineColor(kBlue);
      //tl_g[ibin] = new TLine(1.,gluon[ibin]*gluon_corr_factor[ibin][0],6.,gluon[ibin]*gluon_corr_factor[ibin][0]); tl_g[ibin]->SetLineStyle(1); tl_g[ibin]->SetLineColor(kRed);

      tl_q[ibin] = new TLine(1.,quark_gen[ibin],6.,quark_gen[ibin]); tl_q[ibin]->SetLineStyle(1); tl_q[ibin]->SetLineColor(kBlue);
      tl_g[ibin] = new TLine(1.,gluon_gen[ibin],6.,gluon_gen[ibin]); tl_g[ibin]->SetLineStyle(2); tl_g[ibin]->SetLineColor(kRed);

      tl_q_kappa[ibin] = new TLine(0.2,quark_gen[ibin],0.8,quark_gen[ibin]); tl_q_kappa[ibin]->SetLineStyle(2); tl_q_kappa[ibin]->SetLineColor(kBlue);
      tl_g_kappa[ibin] = new TLine(0.2,gluon_gen[ibin],0.8,gluon_gen[ibin]); tl_g_kappa[ibin]->SetLineStyle(2); tl_g_kappa[ibin]->SetLineColor(kRed);

      tl_ubar[ibin] = new TLine(1.,upbar[ibin],6.,upbar[ibin]); tl_ubar[ibin]->SetLineStyle(3); tl_ubar[ibin]->SetLineColor(kBlue);
      tl_dbar[ibin] = new TLine(1.,downbar[ibin],6.,downbar[ibin]); tl_dbar[ibin]->SetLineStyle(4); tl_dbar[ibin]->SetLineColor(kBlue);
      tl_c[ibin] = new TLine(1.,charm[ibin],6.,charm[ibin]); tl_c[ibin]->SetLineStyle(5); tl_c[ibin]->SetLineColor(kBlue);
      tl_s[ibin] = new TLine(1.,strange[ibin],6.,strange[ibin]); tl_s[ibin]->SetLineStyle(6); tl_s[ibin]->SetLineColor(kBlue);
      tl_b[ibin] = new TLine(1.,bottom[ibin],6.,bottom[ibin]); tl_b[ibin]->SetLineStyle(7); tl_b[ibin]->SetLineColor(kBlue);

      tl_q_ref[ibin] = new TLine(1.,quark_gen[ibin],6.,quark_gen[ibin]); tl_q_ref[ibin]->SetLineStyle(2); tl_q_ref[ibin]->SetLineColor(kBlue);
      tl_g_ref[ibin] = new TLine(1.,gluon_ubar_dbar_gen[ibin],6.,gluon_ubar_dbar_gen[ibin]); tl_g_ref[ibin]->SetLineStyle(2); tl_g_ref[ibin]->SetLineColor(kRed);
      tl_u_ref[ibin] = new TLine(1.,up_gen[ibin],6.,up_gen[ibin]); tl_u_ref[ibin]->SetLineStyle(2); tl_u_ref[ibin]->SetLineColor(kBlue);
      tl_d_ref[ibin] = new TLine(1.,down_gen[ibin],6.,down_gen[ibin]); tl_d_ref[ibin]->SetLineStyle(9); tl_d_ref[ibin]->SetLineColor(kBlue);
  }  

  tl_u_data = new TLine(1.,0.134,6.,0.134); tl_u_data->SetLineStyle(2); tl_u_data->SetLineColor(kBlue);
  tl_d_data = new TLine(1.,0.154,6.,0.154); tl_d_data->SetLineStyle(9); tl_d_data->SetLineColor(kBlue);

////systematics

  Double_t trkeff_up_sys, trkeff_dn_sys, trkeff_pos_sys, jes_up_sys, jes_dn_sys;
  Double_t trkeff_up_kappa_sys, trkeff_dn_kappa_sys, trkeff_pos_kappa_sys;
  Double_t trkeff_sys, jes_res_sys, total_err;
  Double_t trkeff_kappa_sys, jes_res_kappa_sys, total_kappa_err;
  
  for(int ibin = 0; ibin < nCBins; ibin++){
    for(int ibin2 = 0; ibin2 < ntrkbins; ibin2++){
      for(int ibin3 = 1; ibin3 < h_chg_data[0][0]->GetNbinsX()+1; ibin3++){
        if(do_ptcut_fitting && do_systematics){
          //pt cut
          trkeff_up_sys = fabs(h_chg_data[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_trkeff[ibin][ibin2][7]->GetBinContent(ibin3));
          trkeff_dn_sys = fabs(h_chg_data[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_trkeff[ibin][ibin2][3]->GetBinContent(ibin3));
          trkeff_sys = 0.;
          trkeff_sys = (trkeff_up_sys+trkeff_dn_sys)/2.;
          trkeff_pos_sys = 0.;   
          //trkeff_pos_sys = fabs(h_chg_data[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_trkpos[ibin][ibin2]->GetBinContent(ibin3));       
          jes_res_sys = 0.;
          //jes_up_sys = fabs(h_chg_data[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_jeup[ibin][ibin2]->GetBinContent(ibin3));
          //jes_dn_sys = fabs(h_chg_data[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_jedn[ibin][ibin2]->GetBinContent(ibin3));
          //jes_res_sys = (jes_up_sys+jes_dn_sys)/2.;
          //if(ibin>0) jes_res_sys = fabs(h_chg_data[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_jer[ibin][ibin2]->GetBinContent(ibin3));
          //if(ibin3 == 251) cout<<"trkeff_sys: "<<trkeff_sys<<"  trkeff_pos_sys: "<<trkeff_pos_sys<<"  jes_res_sys: "<<jes_res_sys<<endl;
          total_err = sqrt(pow(trkeff_pos_sys,2)+pow(trkeff_sys,2)+pow(jes_res_sys,2)+pow(h_chg_data[ibin][ibin2]->GetBinError(ibin3),2));
          h_chg_data[ibin][ibin2]->SetBinError(ibin3, total_err);
          for(int ibin4=0; ibin4<trkeffbins;ibin4++){
              h_chg_data_trkeff[ibin][ibin2][ibin4]->SetBinError(ibin3, total_err);
          }
        }

        if(do_kappa_fitting && do_systematics){  
          //kappa
          trkeff_up_kappa_sys = fabs(h_chg_data_kappa[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_kappa_trkeff[ibin][ibin2][7]->GetBinContent(ibin3));
          trkeff_dn_kappa_sys = fabs(h_chg_data_kappa[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_kappa_trkeff[ibin][ibin2][3]->GetBinContent(ibin3));
          trkeff_kappa_sys = 0.;
          trkeff_kappa_sys = (trkeff_up_kappa_sys+trkeff_dn_kappa_sys)/2.;
          trkeff_pos_kappa_sys = 0.; 
          //if(ibin>0) trkeff_pos_sys = fabs(h_chg_data_kappa[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_trkpos_kappa[ibin][ibin2]->GetBinContent(ibin3));         
          //jes_up_sys = fabs(h_chg_data[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_jeup[ibin][ibin2]->GetBinContent(ibin3));
          //jes_dn_sys = fabs(h_chg_data[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_jedn[ibin][ibin2]->GetBinContent(ibin3));
          jes_res_kappa_sys = 0.;
          //jes_res_sys = (jes_up_sys+jes_dn_sys)/2.;
          //if(ibin>0) jes_res_sys = fabs(h_chg_data_kappa[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_jer_kappa[ibin][ibin2]->GetBinContent(ibin3));
          total_kappa_err = sqrt(pow(trkeff_pos_kappa_sys,2)+pow(trkeff_kappa_sys,2)+pow(jes_res_kappa_sys,2)+pow(h_chg_data_kappa[ibin][ibin2]->GetBinError(ibin3),2));
          h_chg_data_kappa[ibin][ibin2]->SetBinError(ibin3, total_kappa_err);
          for(int ibin4=0; ibin4<trkeffbins;ibin4++){
              h_chg_data_kappa_trkeff[ibin][ibin2][ibin4]->SetBinError(ibin3, total_kappa_err);
          }
        }
      }
    }
  }
  
  ///defining fitting templates

  f_chg_qg_data[0][0] = new TF1("f_chg_qg_data_cent0_trkpt2",f_qg_reco_cent0_pt2,-2.,2.,1);
  f_chg_qg_data[1][0] = new TF1("f_chg_qg_data_cent1_trkpt2",f_qg_reco_cent1_pt2,-2.,2.,1);
  f_chg_qg_data[2][0] = new TF1("f_chg_qg_data_cent2_trkpt2",f_qg_reco_cent2_pt2,-2.,2.,1);
  f_chg_qg_data[3][0] = new TF1("f_chg_qg_data_cent3_trkpt2",f_qg_reco_cent3_pt2,-2.,2.,1);
  f_chg_qg_data[4][0] = new TF1("f_chg_qg_data_cent4_trkpt2",f_qg_reco_cent4_pt2,-2.,2.,1);

  f_chg_qg_data[0][1] = new TF1("f_chg_qg_data_cent0_trkpt3",f_qg_reco_cent0_pt3,-2.,2.,1);
  f_chg_qg_data[1][1] = new TF1("f_chg_qg_data_cent1_trkpt3",f_qg_reco_cent1_pt3,-2.,2.,1);
  f_chg_qg_data[2][1] = new TF1("f_chg_qg_data_cent2_trkpt3",f_qg_reco_cent2_pt3,-2.,2.,1);
  f_chg_qg_data[3][1] = new TF1("f_chg_qg_data_cent3_trkpt3",f_qg_reco_cent3_pt3,-2.,2.,1);
  f_chg_qg_data[4][1] = new TF1("f_chg_qg_data_cent4_trkpt3",f_qg_reco_cent4_pt3,-2.,2.,1);

  f_chg_qg_data[0][2] = new TF1("f_chg_qg_data_cent0_trkpt4",f_qg_reco_cent0_pt4,-2.,2.,1);
  f_chg_qg_data[1][2] = new TF1("f_chg_qg_data_cent1_trkpt4",f_qg_reco_cent1_pt4,-2.,2.,1);
  f_chg_qg_data[2][2] = new TF1("f_chg_qg_data_cent2_trkpt4",f_qg_reco_cent2_pt4,-2.,2.,1);
  f_chg_qg_data[3][2] = new TF1("f_chg_qg_data_cent3_trkpt4",f_qg_reco_cent3_pt4,-2.,2.,1);
  f_chg_qg_data[4][2] = new TF1("f_chg_qg_data_cent4_trkpt4",f_qg_reco_cent4_pt4,-2.,2.,1);

  f_chg_qg_data[0][3] = new TF1("f_chg_qg_data_cent0_trkpt5",f_qg_reco_cent0_pt5,-2.,2.,1);
  f_chg_qg_data[1][3] = new TF1("f_chg_qg_data_cent1_trkpt5",f_qg_reco_cent1_pt5,-2.,2.,1);
  f_chg_qg_data[2][3] = new TF1("f_chg_qg_data_cent2_trkpt5",f_qg_reco_cent2_pt5,-2.,2.,1);
  f_chg_qg_data[3][3] = new TF1("f_chg_qg_data_cent3_trkpt5",f_qg_reco_cent3_pt5,-2.,2.,1);
  f_chg_qg_data[4][3] = new TF1("f_chg_qg_data_cent4_trkpt5",f_qg_reco_cent4_pt5,-2.,2.,1);

///kappa templates
  f_chg_qg_kappa[0][0] = new TF1("f_chg_qg_kappa_cent0_trkpt2",f_qg_reco_kappa_cent0_pt2,-2.,2.,1);
  f_chg_qg_kappa[1][0] = new TF1("f_chg_qg_kappa_cent1_trkpt2",f_qg_reco_kappa_cent1_pt2,-2.,2.,1);
  f_chg_qg_kappa[2][0] = new TF1("f_chg_qg_kappa_cent2_trkpt2",f_qg_reco_kappa_cent2_pt2,-2.,2.,1);
  f_chg_qg_kappa[3][0] = new TF1("f_chg_qg_kappa_cent3_trkpt2",f_qg_reco_kappa_cent3_pt2,-2.,2.,1);
  f_chg_qg_kappa[4][0] = new TF1("f_chg_qg_kappa_cent4_trkpt2",f_qg_reco_kappa_cent4_pt2,-2.,2.,1);

  f_chg_qg_kappa[0][1] = new TF1("f_chg_qg_kappa_cent0_trkpt3",f_qg_reco_kappa_cent0_pt3,-2.,2.,1);
  f_chg_qg_kappa[1][1] = new TF1("f_chg_qg_kappa_cent1_trkpt3",f_qg_reco_kappa_cent1_pt3,-2.,2.,1);
  f_chg_qg_kappa[2][1] = new TF1("f_chg_qg_kappa_cent2_trkpt3",f_qg_reco_kappa_cent2_pt3,-2.,2.,1);
  f_chg_qg_kappa[3][1] = new TF1("f_chg_qg_kappa_cent3_trkpt3",f_qg_reco_kappa_cent3_pt3,-2.,2.,1);
  f_chg_qg_kappa[4][1] = new TF1("f_chg_qg_kappa_cent4_trkpt3",f_qg_reco_kappa_cent4_pt3,-2.,2.,1);

  f_chg_qg_kappa[0][2] = new TF1("f_chg_qg_kappa_cent0_trkpt4",f_qg_reco_kappa_cent0_pt4,-2.,2.,1);
  f_chg_qg_kappa[1][2] = new TF1("f_chg_qg_kappa_cent1_trkpt4",f_qg_reco_kappa_cent1_pt4,-2.,2.,1);
  f_chg_qg_kappa[2][2] = new TF1("f_chg_qg_kappa_cent2_trkpt4",f_qg_reco_kappa_cent2_pt4,-2.,2.,1);
  f_chg_qg_kappa[3][2] = new TF1("f_chg_qg_kappa_cent3_trkpt4",f_qg_reco_kappa_cent3_pt4,-2.,2.,1);
  f_chg_qg_kappa[4][2] = new TF1("f_chg_qg_kappa_cent4_trkpt4",f_qg_reco_kappa_cent4_pt4,-2.,2.,1);

  f_chg_qg_kappa[0][3] = new TF1("f_chg_qg_kappa_cent0_trkpt5",f_qg_reco_kappa_cent0_pt5,-2.,2.,1);
  f_chg_qg_kappa[1][3] = new TF1("f_chg_qg_kappa_cent1_trkpt5",f_qg_reco_kappa_cent1_pt5,-2.,2.,1);
  f_chg_qg_kappa[2][3] = new TF1("f_chg_qg_kappa_cent2_trkpt5",f_qg_reco_kappa_cent2_pt5,-2.,2.,1);
  f_chg_qg_kappa[3][3] = new TF1("f_chg_qg_kappa_cent3_trkpt5",f_qg_reco_kappa_cent3_pt5,-2.,2.,1);
  f_chg_qg_kappa[4][3] = new TF1("f_chg_qg_kappa_cent4_trkpt5",f_qg_reco_kappa_cent4_pt5,-2.,2.,1);

  TLine *tl1 = new TLine(-2.,1.,2.,1.);
  tl1->SetLineStyle(2);
  TLine *tl2 = new TLine(120.,1.,400.,1.);
  tl2->SetLineStyle(2);
  TLine *tl3 = new TLine(-1.5,0.,1.5,0.);
  tl3->SetLineStyle(2);

  TLatex *tx, *tx1, *tx2, *tx3, *tx4, *tx5;

  char g_frac[50],g_error[50],chi2ndf[50], up_frac[50], down_frac[50];

  TString tmp;
  TString pythia = "MC";
  TString ppdata = "Data";
  TString data = "pp 27.4 pb^{-1} (5.02 TeV)  PbPb 404 #mub^{-1} (5.02 TeV)";
  TString cms = "CMS";
  TString gen_tex = "MC Gen";
  TString chg_label = "1/N_{jets}";
  TString ratio = "Data / MC";
  TString fraction = "gluon like jet fraction";
  TString y_label =  "1/N_{jets} (dN/dQ) [1/e]";
  TString datafit =  "Fit / Data";
  TString diff =  "diff";

  TLatex *title_tex;
  title_tex = new TLatex();
  title_tex->SetTextSize(0.75);
  title_tex->SetLineColor(kWhite);

  TLatex *xlabel_tex;
  xlabel_tex = new TLatex();
  xlabel_tex->SetTextSize(0.75);
  xlabel_tex->SetLineColor(kWhite);

  TLatex *ylabel_tex;
  ylabel_tex = new TLatex();
  ylabel_tex->SetTextAngle(90);
  ylabel_tex->SetTextSize(0.65);
  ylabel_tex->SetLineColor(kWhite);
  
  TLatex *chg_tex;
  chg_tex = new TLatex();
  chg_tex->SetTextAngle(90);
  chg_tex->SetTextSize(0.65);
  chg_tex->SetLineColor(kWhite);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_qg_eta = new TCanvas("c_qg_eta","",1500,300);
  c_qg_eta->Divide(5,1,0);
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0) c_qg_eta->cd(ibin+1);
    else c_qg_eta->cd(6-ibin);   
    h_eta_MC_g[ibin]->GetXaxis()->SetTitle("jet #eta");
    h_cosm(h_eta_MC_g[ibin]);
    h_eta_MC_g[ibin]->GetXaxis()->SetRangeUser(-1.5,1.52);
    h_eta_MC_g[ibin]->GetYaxis()->SetRangeUser(0.05,0.15);
    h_eta_MC_g[ibin]->SetLineColor(kRed); h_eta_MC_g[ibin]->SetMarkerColor(kRed);
    h_eta_MC_g[ibin]->SetMarkerStyle(4); h_eta_MC_g[ibin]->Draw("e0 same");
    h_eta_MC_q[ibin]->SetLineColor(kBlue); h_eta_MC_q[ibin]->SetMarkerColor(kBlue);
    h_eta_MC_q[ibin]->SetMarkerStyle(4); h_eta_MC_q[ibin]->Draw("e0 same");
    h_eta_MC[ibin]->SetLineColor(kBlack); h_eta_MC[ibin]->SetMarkerColor(kBlack);
    h_eta_MC[ibin]->SetMarkerStyle(20); h_eta_MC[ibin]->Draw("e0 same");
    h_eta_MC_up[ibin]->SetLineColor(kBlue); h_eta_MC_up[ibin]->SetMarkerColor(kBlue);
    h_eta_MC_up[ibin]->SetMarkerStyle(22); h_eta_MC_up[ibin]->Draw("e0 same");
    h_eta_MC_d[ibin]->SetLineColor(kBlue); h_eta_MC_d[ibin]->SetMarkerColor(kBlue);
    h_eta_MC_d[ibin]->SetMarkerStyle(23); h_eta_MC_d[ibin]->Draw("e0 same");

    if(ibin==0){
      TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry(h_eta_MC[ibin], "inc", "lepf");
      legend ->AddEntry(h_eta_MC_q[ibin], "quark", "lepf");
      legend ->AddEntry(h_eta_MC_g[ibin], "gluon", "lepf");
      legend ->AddEntry(h_eta_MC_up[ibin], "up", "lepf");
      legend ->AddEntry(h_eta_MC_d[ibin], "down", "lepf");
      legend ->Draw("same");
    } 

    TString tmp1 = cent_tag_pp_MC[ibin];
    tx1 = new TLatex(); tx1->SetTextSize(.12);
    tx1->DrawLatexNDC(0.35,0.85, tmp1);
  }

  TCanvas *c_reco = new TCanvas("c_reco","c_reco",1500,600);
  c_reco->Divide(5,2,0);
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0) c_reco->cd(ibin+1);
    else c_reco->cd(6-ibin);   
    gPad->SetLogy();
/*
    h_reco_MC[ibin]->Rebin(2);
    h_reco_MC_up[ibin]->Rebin(2);
    h_reco_MC_d[ibin]->Rebin(2);
    h_reco_MC_g[ibin]->Rebin(5);
    h_reco_MC_ubdbcsb[ibin]->Rebin(2);
*/
    h_reco_MC[ibin]->GetXaxis()->SetTitle("jet pT");
    h_cosm(h_reco_MC[ibin]);
    h_reco_MC[ibin]->GetXaxis()->SetRangeUser(120.,400.);
    h_reco_MC[ibin]->SetLineColor(kGreen-2); h_reco_MC[ibin]->SetMarkerColor(kBlack); h_reco_MC[ibin]->SetFillColor(kGreen-2);
    h_reco_MC[ibin]->SetMarkerStyle(4); h_reco_MC[ibin]->Draw("hist same");
    //h_reco_data[ibin]->SetLineColor(kGreen); h_reco_data[ibin]->SetMarkerColor(kGreen);
    //h_reco_data[ibin]->SetMarkerStyle(4); h_reco_data[ibin]->Draw("e0 same");
    h_reco_MC_q[ibin]->SetLineColor(kBlue); h_reco_MC_q[ibin]->SetMarkerColor(kBlue);
    h_reco_MC_q[ibin]->SetMarkerStyle(4); //h_reco_MC_q[ibin]->Draw("e0 same");
    h_reco_MC_g[ibin]->SetLineColor(kRed); h_reco_MC_g[ibin]->SetMarkerColor(kRed); h_reco_MC_g[ibin]->SetFillColor(kRed);
    h_reco_MC_g[ibin]->SetMarkerStyle(4); h_reco_MC_g[ibin]->Draw("hist same");

    h_reco_MC_up[ibin]->SetLineColor(kCyan); h_reco_MC_up[ibin]->SetMarkerColor(kCyan); h_reco_MC_up[ibin]->SetFillColor(kCyan);
    h_reco_MC_up[ibin]->SetMarkerStyle(22); h_reco_MC_up[ibin]->Draw("hist same");
    h_reco_MC_upbar[ibin]->SetLineColor(kBlue); h_reco_MC_upbar[ibin]->SetMarkerColor(kBlue);
    //h_reco_MC_upbar[ibin]->SetMarkerStyle(26); h_reco_MC_upbar[ibin]->Draw("e0 same");

    h_reco_MC_d[ibin]->SetLineColor(kMagenta); h_reco_MC_d[ibin]->SetMarkerColor(kMagenta); h_reco_MC_d[ibin]->SetFillColor(kMagenta);
    h_reco_MC_d[ibin]->SetMarkerStyle(23); h_reco_MC_d[ibin]->Draw("hist same");
    h_reco_MC_dbar[ibin]->SetLineColor(kBlue); h_reco_MC_dbar[ibin]->SetMarkerColor(kBlue);
    //h_reco_MC_dbar[ibin]->SetMarkerStyle(32); h_reco_MC_dbar[ibin]->Draw("e0 same");
    h_reco_MC_c[ibin]->SetLineColor(kOrange+1); h_reco_MC_c[ibin]->SetMarkerColor(kOrange+1); h_reco_MC_c[ibin]->SetFillColor(kOrange+1);
    //h_reco_MC_c[ibin]->SetMarkerStyle(28); h_reco_MC_c[ibin]->Draw("e0 same");
    h_reco_MC_s[ibin]->SetLineColor(kAzure+1); h_reco_MC_s[ibin]->SetMarkerColor(kAzure+1); h_reco_MC_s[ibin]->SetFillColor(kAzure+1);
    //h_reco_MC_s[ibin]->SetMarkerStyle(29); h_reco_MC_s[ibin]->Draw("e0 same");
    h_reco_MC_b[ibin]->SetLineColor(kGreen-6); h_reco_MC_b[ibin]->SetMarkerColor(kGreen-6); h_reco_MC_b[ibin]->SetFillColor(kGreen-6);
    //h_reco_MC_b[ibin]->SetMarkerStyle(30); h_reco_MC_b[ibin]->Draw("e0 same");

    h_reco_MC_ubdbcsb[ibin]->SetLineColor(kOrange+1); h_reco_MC_ubdbcsb[ibin]->SetMarkerColor(kOrange+1); h_reco_MC_ubdbcsb[ibin]->SetFillColor(kOrange+1);
    h_reco_MC_ubdbcsb[ibin]->SetMarkerStyle(28); h_reco_MC_ubdbcsb[ibin]->Draw("hist same");

    if(ibin==0){
      TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
      legend ->SetLineColor(kWhite);
      //legend ->AddEntry(h_reco_data[ibin], "data inclusive", "lepf");
      legend ->AddEntry(h_reco_MC[ibin], "inclusive", "f");
      //legend ->AddEntry(h_reco_MC_q[ibin], "quark", "lepf");
      legend ->AddEntry(h_reco_MC_g[ibin], "gluon", "f");
      legend ->AddEntry(h_reco_MC_up[ibin], "up quark", "f");
      legend ->AddEntry(h_reco_MC_d[ibin], "down quark", "f");
      legend ->AddEntry(h_reco_MC_ubdbcsb[ibin], "other quarks", "f");
      //legend ->AddEntry(h_reco_MC_upbar[ibin], "upbar", "lepf");
      //legend ->AddEntry(h_reco_MC_dbar[ibin], "downbar", "lepf");
      //legend ->AddEntry(h_reco_MC_c[ibin], "c+cbar", "lepf");
      //legend ->AddEntry(h_reco_MC_s[ibin], "s+sbar", "lepf");
      //legend ->AddEntry(h_reco_MC_b[ibin], "b+bbar", "lepf");
      legend ->Draw("same");
    } 
    
    TString tmp1 = cent_tag_pp_MC[ibin];
    tx1 = new TLatex(); tx1->SetTextSize(.1);
    tx1->DrawLatexNDC(0.15,0.85, tmp1);
      
    if(ibin==0) c_reco->cd(ibin+6);
    else c_reco->cd(11-ibin);   
    h_reco_gluon_inc_ratio[ibin] = (TH1D*)h_reco_MC_g[ibin]->Clone("h_reco_gluon_inc_ratio");
    h_reco_gluon_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
    h_reco_up_inc_ratio[ibin] = (TH1D*)h_reco_MC_up[ibin]->Clone("h_reco_up_inc_ratio");
    h_reco_up_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
    //h_reco_upbar_inc_ratio[ibin] = (TH1D*)h_reco_MC_upbar[ibin]->Clone("h_reco_upbar_inc_ratio");
    //h_reco_upbar_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
    h_reco_d_inc_ratio[ibin] = (TH1D*)h_reco_MC_d[ibin]->Clone("h_reco_d_inc_ratio");
    h_reco_d_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
    h_reco_dbar_inc_ratio[ibin] = (TH1D*)h_reco_MC_ubdbcsb[ibin]->Clone("h_reco_ubdbcsb_inc_ratio");
    h_reco_dbar_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
/*
    h_reco_c_inc_ratio[ibin] = (TH1D*)h_reco_MC_c[ibin]->Clone("h_reco_c_inc_ratio");
    h_reco_c_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
    h_reco_s_inc_ratio[ibin] = (TH1D*)h_reco_MC_s[ibin]->Clone("h_reco_s_inc_ratio");
    h_reco_s_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
    h_reco_b_inc_ratio[ibin] = (TH1D*)h_reco_MC_b[ibin]->Clone("h_reco_b_inc_ratio");
    h_reco_b_inc_ratio[ibin]->Divide(h_reco_MC[ibin]);
*/
    h_reco_gluon_inc_ratio[ibin]->GetYaxis()->SetNdivisions(505);
    h_reco_gluon_inc_ratio[ibin]->GetYaxis()->SetTitleSize(0.05);
    h_cosm(h_reco_gluon_inc_ratio[ibin]);
    h_reco_gluon_inc_ratio[ibin]->GetXaxis()->SetRangeUser(120.,400.);
    h_reco_gluon_inc_ratio[ibin]->GetYaxis()->SetRangeUser(0.,0.8);
    h_reco_gluon_inc_ratio[ibin]->GetYaxis()->SetTitle("Ratio");
    h_reco_gluon_inc_ratio[ibin]->GetXaxis()->SetTitle("jet p_{T}");
    h_reco_gluon_inc_ratio[ibin]->SetLineColor(kRed); h_reco_gluon_inc_ratio[ibin]->SetMarkerColor(kRed);
    h_reco_gluon_inc_ratio[ibin]->SetMarkerStyle(20); h_reco_gluon_inc_ratio[ibin]->Draw("e0 same");
    h_reco_up_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_up_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_up_inc_ratio[ibin]->SetMarkerStyle(22); h_reco_up_inc_ratio[ibin]->Draw("e0 same");
    h_reco_d_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_d_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_d_inc_ratio[ibin]->SetMarkerStyle(23); h_reco_d_inc_ratio[ibin]->Draw("e0 same");
    h_reco_dbar_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_dbar_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_dbar_inc_ratio[ibin]->SetMarkerStyle(20); h_reco_dbar_inc_ratio[ibin]->Draw("e0 same");
/*
    h_reco_upbar_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_upbar_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_upbar_inc_ratio[ibin]->SetMarkerStyle(26); h_reco_upbar_inc_ratio[ibin]->Draw("e0 same");
    h_reco_c_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_c_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_c_inc_ratio[ibin]->SetMarkerStyle(28); h_reco_c_inc_ratio[ibin]->Draw("e0 same");
    h_reco_s_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_s_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_s_inc_ratio[ibin]->SetMarkerStyle(29); h_reco_s_inc_ratio[ibin]->Draw("e0 same");
    h_reco_b_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_b_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_b_inc_ratio[ibin]->SetMarkerStyle(30); h_reco_b_inc_ratio[ibin]->Draw("e0 same");
*/
    drawline(120.,1.,400.,1.);

    if(ibin==0){
      TLegend *legendp = new TLegend(0.25,0.75,0.75,0.9);
      legendp ->SetLineColor(kWhite);
      legendp ->AddEntry(h_reco_gluon_inc_ratio[ibin], "gluon / inclusive", "lep");
      legendp ->AddEntry(h_reco_up_inc_ratio[ibin], "up quark / inclusive", "lep");
      legendp ->AddEntry(h_reco_d_inc_ratio[ibin], "down quark / inclusive", "lep");
      legendp ->AddEntry(h_reco_dbar_inc_ratio[ibin], "other quarks / inclusive", "lep");
      legendp ->Draw("same");
    } 

  }

  TCanvas *c_reco_gb = new TCanvas("c_reco_gb","c_reco_gb",1500,600);
  c_reco_gb->Divide(5,2,0);
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1);

  for(int ibin=0;ibin<nCBins;ibin++){
    if(ibin==0) c_reco_gb->cd(ibin+1);
    else c_reco_gb->cd(6-ibin);   
    gPad->SetLogy();
/*
    h_reco_MC_upbar[ibin]->Rebin(2);   
    h_reco_MC_dbar[ibin]->Rebin(2);    
    h_reco_MC_c[ibin]->Rebin(2);    
    h_reco_MC_s[ibin]->Rebin(2);    
    h_reco_MC_b[ibin]->Rebin(2);    
    h_reco_MC_g[ibin]->Rebin(2);    
*/
    h_reco_MC_g[ibin]->GetXaxis()->SetTitle("jet pT");
    h_cosm(h_reco_MC_g[ibin]);
    h_reco_MC_g[ibin]->GetXaxis()->SetRangeUser(120.,400.);
    h_reco_MC_g[ibin]->SetLineColor(kRed-2); h_reco_MC_g[ibin]->SetMarkerColor(kRed-2); h_reco_MC_g[ibin]->SetFillColor(kRed-2);
    h_reco_MC_g[ibin]->SetMarkerStyle(4); h_reco_MC_g[ibin]->Draw("hist same");
    h_reco_MC_s[ibin]->Draw("hist same");
    h_reco_MC_dbar[ibin]->SetLineColor(kMagenta); h_reco_MC_dbar[ibin]->SetMarkerColor(kMagenta); h_reco_MC_dbar[ibin]->SetMarkerSize(1.5); h_reco_MC_dbar[ibin]->SetFillColor(kMagenta);
    h_reco_MC_dbar[ibin]->SetMarkerStyle(23); h_reco_MC_dbar[ibin]->Draw("hist same");
    h_reco_MC_upbar[ibin]->SetLineColor(kCyan); h_reco_MC_upbar[ibin]->SetMarkerColor(kCyan); h_reco_MC_upbar[ibin]->SetMarkerSize(1.5); h_reco_MC_upbar[ibin]->SetFillColor(kCyan);
    h_reco_MC_upbar[ibin]->SetMarkerStyle(22); h_reco_MC_upbar[ibin]->Draw("hist same");
    h_reco_MC_c[ibin]->Draw("hist same");
    h_reco_MC_b[ibin]->Draw("hist same");

    if(ibin==0){
      TLegend *legend = new TLegend(0.55,0.75,0.75,0.9);
      legend ->SetLineColor(kWhite);
      legend ->AddEntry(h_reco_MC_g[ibin], "gluon", "f");
      legend ->AddEntry(h_reco_MC_upbar[ibin], "upbar", "f");
      legend ->AddEntry(h_reco_MC_dbar[ibin], "downbar", "f");
      legend ->AddEntry(h_reco_MC_c[ibin], "c+cbar", "f");
      legend ->AddEntry(h_reco_MC_s[ibin], "s+sbar", "f");
      legend ->AddEntry(h_reco_MC_b[ibin], "b+bbar", "f");
      legend ->Draw("same");
    } 
    
    TString tmp1 = cent_tag_pp_MC[ibin];
    tx1 = new TLatex(); tx1->SetTextSize(.1);
    tx1->DrawLatexNDC(0.15,0.85, tmp1);   

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
    h_reco_gluon_dbar_ratio[ibin]->GetXaxis()->SetTitle("jet pT");
    h_reco_gluon_dbar_ratio[ibin]->GetXaxis()->CenterTitle();
    h_reco_gluon_dbar_ratio[ibin]->GetYaxis()->SetNdivisions(505);
    h_reco_gluon_dbar_ratio[ibin]->GetYaxis()->SetTitleSize(0.05);
    h_reco_gluon_dbar_ratio[ibin]->GetYaxis()->SetLabelSize(0.07);
    h_reco_gluon_dbar_ratio[ibin]->GetXaxis()->CenterTitle();
    h_reco_gluon_dbar_ratio[ibin]->GetXaxis()->SetNdivisions(505);
    h_reco_gluon_dbar_ratio[ibin]->GetXaxis()->SetTitleSize(0.05);
    h_reco_gluon_dbar_ratio[ibin]->GetXaxis()->SetLabelSize(0.07);
    h_reco_gluon_dbar_ratio[ibin]->GetXaxis()->SetRangeUser(120.,400.);
    h_reco_gluon_dbar_ratio[ibin]->GetYaxis()->SetRangeUser(0.,0.15);
    //h_reco_gluon_dbar_ratio[ibin]->SetLineColor(kBlue); h_reco_gluon_dbar_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_gluon_dbar_ratio[ibin]->SetMarkerStyle(23);  h_reco_gluon_dbar_ratio[ibin]->Draw("e0 same");
    //h_reco_gluon_upbar_ratio[ibin]->SetLineColor(kBlue); h_reco_gluon_upbar_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_gluon_upbar_ratio[ibin]->SetMarkerStyle(22); h_reco_gluon_upbar_ratio[ibin]->Draw("e0 same");
    //h_reco_gluon_c_ratio[ibin]->SetLineColor(kBlue); h_reco_gluon_c_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_gluon_c_ratio[ibin]->SetMarkerStyle(28); h_reco_gluon_c_ratio[ibin]->SetMarkerSize(1.5); h_reco_gluon_c_ratio[ibin]->Draw("e0 same");
    //h_reco_gluon_s_ratio[ibin]->SetLineColor(kBlue); h_reco_gluon_s_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_gluon_s_ratio[ibin]->SetMarkerStyle(29); h_reco_gluon_s_ratio[ibin]->SetMarkerSize(1.5); h_reco_gluon_s_ratio[ibin]->Draw("e0 same");
    //h_reco_gluon_b_ratio[ibin]->SetLineColor(kBlue); h_reco_gluon_b_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_gluon_b_ratio[ibin]->SetMarkerStyle(30); h_reco_gluon_b_ratio[ibin]->SetMarkerSize(1.5); h_reco_gluon_b_ratio[ibin]->Draw("e0 same");

    drawline(120.,1.,400.,1.);

    if(ibin==0){
      TLegend *legendp = new TLegend(0.55,0.75,0.75,0.9);
      legendp ->SetLineColor(kWhite);
      legendp ->AddEntry(h_reco_gluon_upbar_ratio[ibin], "upbar / gluon", "lep");
      legendp ->AddEntry(h_reco_gluon_dbar_ratio[ibin], "downbar / gluon", "lep");
      legendp ->AddEntry(h_reco_gluon_c_ratio[ibin], "c+cbar / gluon", "lep");
      legendp ->AddEntry(h_reco_gluon_s_ratio[ibin], "s+sbar / gluon", "lep");
      legendp ->AddEntry(h_reco_gluon_b_ratio[ibin], "b+bbar / gluon", "lep");
      legendp ->Draw("same");
    } 

  }


///////////////////////////////////// fitting /////////////////////////////////////////

  h_dummy->GetXaxis()->SetTitle("p_{T} cut");
  //h_dummy->GetXaxis()->SetTitle("");
  h_dummy->GetYaxis()->SetTitle("");
  h_dummy->GetYaxis()->SetNdivisions(505);
  h_dummy->GetYaxis()->SetTitleSize(0.05);
  h_dummy->GetYaxis()->SetTitleOffset(0.75);
  h_dummy->GetYaxis()->SetLabelSize(0.07);
  h_dummy->GetXaxis()->CenterTitle();
  h_dummy->GetXaxis()->SetNdivisions(505);
  h_dummy->GetXaxis()->SetTitleSize(0.05);
  h_dummy->GetXaxis()->SetLabelSize(0.07);
  if(do_ptcut_fitting) h_dummy->GetXaxis()->SetRangeUser(1.,6.);
  if(do_kappa_fitting) h_dummy->GetXaxis()->SetRangeUser(0.,1.);
  h_dummy->GetYaxis()->SetRangeUser(0.1,1.);

  h_dummy2->GetXaxis()->SetTitle("p_{T} cut");
  h_dummy2->GetYaxis()->SetTitle("");
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
  h_dummy3->GetXaxis()->SetRangeUser(120.,190.);
  h_dummy3->GetYaxis()->SetRangeUser(-0.05,0.1);

  h_dummy4->GetXaxis()->SetTitle("jet p_{T}");
  h_dummy4->GetYaxis()->SetTitle("#sigma");
  h_dummy4->GetYaxis()->SetNdivisions(505);
  h_dummy4->GetYaxis()->SetTitleSize(0.05);
  h_dummy4->GetYaxis()->SetTitleOffset(0.75);
  h_dummy4->GetYaxis()->SetLabelSize(0.07);
  h_dummy4->GetXaxis()->CenterTitle();
  h_dummy4->GetXaxis()->SetNdivisions(505);
  h_dummy4->GetXaxis()->SetTitleSize(0.05);
  h_dummy4->GetXaxis()->SetLabelSize(0.07);
  h_dummy4->GetXaxis()->SetRangeUser(120.,200.);
  h_dummy4->GetYaxis()->SetRangeUser(0.,0.6);

  h_dummy5->GetXaxis()->SetTitle("#kappa");
  h_dummy5->GetYaxis()->SetTitle("");
  h_dummy5->GetYaxis()->SetNdivisions(505);
  h_dummy5->GetYaxis()->SetTitleSize(0.05);
  h_dummy5->GetYaxis()->SetTitleOffset(0.75);
  h_dummy5->GetYaxis()->SetLabelSize(0.07);
  h_dummy5->GetXaxis()->CenterTitle();
  h_dummy5->GetXaxis()->SetNdivisions(505);
  h_dummy5->GetXaxis()->SetTitleSize(0.05);
  h_dummy5->GetXaxis()->SetLabelSize(0.07);
  h_dummy5->GetXaxis()->SetRangeUser(2.,2.5);
  h_dummy5->GetYaxis()->SetRangeUser(0.15,0.85);

  h_dummy7->GetXaxis()->SetTitle("Track p_{T} cut (GeV)");
  h_dummy7->GetYaxis()->SetTitle("Quark jet fraction (PbPb/pp)");
  h_dummy7->GetYaxis()->SetNdivisions(505);
  h_dummy7->GetYaxis()->SetTitleSize(0.08);
  h_dummy7->GetYaxis()->SetTitleOffset(0.75);
  h_dummy7->GetXaxis()->SetTitleOffset(0.75);
  h_dummy7->GetYaxis()->SetLabelSize(0.07);
  h_dummy7->GetXaxis()->CenterTitle();
  h_dummy7->GetYaxis()->CenterTitle();
  h_dummy7->GetXaxis()->SetNdivisions(505);
  h_dummy7->GetXaxis()->SetTitleSize(0.06);
  h_dummy7->GetXaxis()->SetLabelSize(0.05);
  h_dummy7->GetXaxis()->SetRangeUser(1.8,5.2);
  h_dummy7->GetXaxis()->SetRangeUser(0.,1.);
  //h_dummy7->GetYaxis()->SetRangeUser(0.3,1.9);

  h_dummy8->GetXaxis()->SetTitle("Track p_{T} cut (GeV)");
  h_dummy8->GetYaxis()->SetTitle("Gluon jet fraction (PbPb/pp)");
  h_dummy8->GetYaxis()->SetNdivisions(505);
  h_dummy8->GetYaxis()->SetTitleSize(0.08);
  h_dummy8->GetYaxis()->SetTitleOffset(0.75);
  h_dummy8->GetXaxis()->SetTitleOffset(0.75);
  h_dummy8->GetYaxis()->SetLabelSize(0.07);
  h_dummy8->GetXaxis()->CenterTitle();
  h_dummy8->GetYaxis()->CenterTitle();
  h_dummy8->GetXaxis()->SetNdivisions(505);
  h_dummy8->GetXaxis()->SetTitleSize(0.06);
  h_dummy8->GetXaxis()->SetLabelSize(0.05);
  h_dummy8->GetXaxis()->SetRangeUser(1.8,5.2);
  h_dummy8->GetXaxis()->SetRangeUser(0.,1.);
  h_dummy8->GetYaxis()->SetRangeUser(0.6,1.3);

  h_dummy9->GetXaxis()->SetTitle("Track p_{T} cut (GeV)");
  h_dummy9->GetYaxis()->SetTitle("Quark jet fraction");
  h_dummy9->GetYaxis()->SetNdivisions(505);
  h_dummy9->GetYaxis()->SetTitleSize(0.08);
  h_dummy9->GetYaxis()->SetTitleOffset(0.75);
  h_dummy9->GetXaxis()->SetTitleOffset(0.75);
  h_dummy9->GetYaxis()->SetLabelSize(0.06);
  h_dummy9->GetXaxis()->CenterTitle();
  h_dummy9->GetYaxis()->CenterTitle();
  h_dummy9->GetXaxis()->SetNdivisions(505);
  h_dummy9->GetXaxis()->SetTitleSize(0.06);
  h_dummy9->GetXaxis()->SetLabelSize(0.05);
  h_dummy9->GetXaxis()->SetRangeUser(1.8,5.2);
  //h_dummy9->GetXaxis()->SetRangeUser(0.3,2.);
  h_dummy9->GetYaxis()->SetRangeUser(0.3,0.85);

  h_dummy10->GetXaxis()->SetTitle("Track p_{T} cut (GeV)");
  h_dummy10->GetYaxis()->SetTitle("Gluon jet fraction");
  h_dummy10->GetYaxis()->SetNdivisions(505);
  h_dummy10->GetYaxis()->SetTitleSize(0.08);
  h_dummy10->GetYaxis()->SetTitleOffset(0.75);
  h_dummy10->GetXaxis()->SetTitleOffset(0.75);
  h_dummy10->GetYaxis()->SetLabelSize(0.06);
  h_dummy10->GetXaxis()->CenterTitle();
  h_dummy10->GetYaxis()->CenterTitle();
  h_dummy10->GetXaxis()->SetNdivisions(505);
  h_dummy10->GetXaxis()->SetTitleSize(0.06);
  h_dummy10->GetXaxis()->SetLabelSize(0.05);
  h_dummy10->GetXaxis()->SetRangeUser(1.8,5.2);
  //h_dummy10->GetXaxis()->SetRangeUser(0.3,2.);
  h_dummy10->GetYaxis()->SetRangeUser(0.3,0.85);
/*
  h_dummy->GetXaxis()->SetBinLabel(2,x_label[0]);
  h_dummy->GetXaxis()->SetBinLabel(5,x_label[1]);
  h_dummy->GetXaxis()->SetBinLabel(8,x_label[2]);
*/
/*
  for(int ibin=0; ibin<nCBins; ibin++) {
    for(int ibin2=0; ibin2<ntrkbins; ibin2++) {
      if(ibin==0) h_reco_gen_q[ibin2]->GetXaxis()->SetBinLabel(ibin+2,cent_tag_pp[ibin]);
      else h_reco_gen_q[ibin2]->GetXaxis()->SetBinLabel(7-ibin,cent_tag_pp[ibin]);
      if(ibin==0) h_reco_gen_g[ibin2]->GetXaxis()->SetBinLabel(ibin+2,cent_tag_pp[ibin]);
      else h_reco_gen_g[ibin2]->GetXaxis()->SetBinLabel(7-ibin,cent_tag_pp[ibin]);
    }
  }
*/
  Double_t up_fit[nCBins][ntrkbins-1]; Double_t up_fiterr[nCBins][ntrkbins-1];
  Double_t down_fit[nCBins][ntrkbins-1]; Double_t down_fiterr[nCBins][ntrkbins-1];
  Double_t gluon_fit[nCBins][ntrkbins-1]; Double_t gluon_fiterr[nCBins][ntrkbins-1];
  Double_t quark_fit[nCBins][ntrkbins-1]; Double_t quark_fiterr[nCBins][ntrkbins-1];

  Double_t quark_fit_high[nCBins][ntrkbins]; Double_t quark_fit_high_diff[nCBins][ntrkbins-1];
  Double_t quark_fit_low[nCBins][ntrkbins]; Double_t quark_fit_low_diff[nCBins][ntrkbins-1];
  Double_t gluon_fit_high[nCBins][ntrkbins]; Double_t gluon_fit_high_diff[nCBins][ntrkbins-1];
  Double_t gluon_fit_low[nCBins][ntrkbins]; Double_t gluon_fit_low_diff[nCBins][ntrkbins-1];

  Double_t quark_fit_kappa_high[nCBins][ntrkbins]; Double_t quark_fit_kappa_high_diff[nCBins][ntrkbins-1];
  Double_t quark_fit_kappa_low[nCBins][ntrkbins]; Double_t quark_fit_kappa_low_diff[nCBins][ntrkbins-1];
  Double_t gluon_fit_kappa_high[nCBins][ntrkbins]; Double_t gluon_fit_kappa_high_diff[nCBins][ntrkbins-1];
  Double_t gluon_fit_kappa_low[nCBins][ntrkbins]; Double_t gluon_fit_kappa_low_diff[nCBins][ntrkbins-1];

  Double_t up_fit_kappa[nCBins][ntrkbins-1]; Double_t up_fiterr_kappa[nCBins][ntrkbins-1];
  Double_t down_fit_kappa[nCBins][ntrkbins-1]; Double_t down_fiterr_kappa[nCBins][ntrkbins-1];
  Double_t gluon_fit_kappa[nCBins][ntrkbins-1]; Double_t gluon_fiterr_kappa[nCBins][ntrkbins-1];
  Double_t quark_fit_kappa[nCBins][ntrkbins-1]; Double_t quark_fiterr_kappa[nCBins][ntrkbins-1];

  Double_t gluon_fitsys[nCBins][ntrkbins-1]; Double_t gluon_fitsyserr[nCBins][ntrkbins-1];
  Double_t quark_fitsys[nCBins][ntrkbins-1]; Double_t quark_fitsyserr[nCBins][ntrkbins-1];
  Double_t gluon_reco_fitsys[nCBins][ntrkbins-1]; Double_t gluon_reco_fitsyserr[nCBins][ntrkbins-1];
  Double_t quark_reco_fitsys[nCBins][ntrkbins-1]; Double_t quark_reco_fitsyserr[nCBins][ntrkbins-1];
  
  Double_t gluon_ratio[nCBins][ntrkbins-1]; Double_t quark_ratio[nCBins][ntrkbins-1];
  Double_t gluon_recofit_ratio[nCBins][ntrkbins-1]; Double_t quark_recofit_ratio[nCBins][ntrkbins-1];


  if(!do_data){
  ///two parameter fitting
    if(do_ptcut_fitting){
      TCanvas *c_chg_qg_MC_reco[ntrkbins];
      for(int ibin2=0;ibin2<ntrkbins;ibin2++){
        c_chg_qg_MC_reco[ibin2] = new TCanvas((TString)("c_chg_qg_MC_reco"+trk[ibin2]),(TString)("c_chg_qg_MC_reco"+trk[ibin2]),1800,600);
        c_chg_qg_MC_reco[ibin2]->Divide(5,2,0);
        gStyle->SetOptStat(0);
        for(int ibin=0;ibin<nCBins;ibin++){
          if(ibin==0) c_chg_qg_MC_reco[ibin2]->cd(1);
          else c_chg_qg_MC_reco[ibin2]->cd(6-ibin);   
          h_chg_reco[ibin][ibin2]->GetXaxis()->SetTitle("reco jet chg");
          h_chg_reco[ibin][ibin2]->GetYaxis()->SetTitle("# of jets");
          h_cosm(h_chg_reco[ibin][ibin2]);
          h_chg_reco[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
          //h_chg_reco[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);

          h_chg_reco[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_reco[ibin][ibin2]->SetMarkerStyle(2);
          h_chg_reco_gubdb[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_gubdb[ibin][ibin2]->SetMarkerColor(kRed); h_chg_reco_gubdb[ibin][ibin2]->SetMarkerStyle(2);
          h_chg_reco_ud_scaled[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_ud_scaled[ibin][ibin2]->SetMarkerColor(kBlue); h_chg_reco_ud_scaled[ibin][ibin2]->SetMarkerStyle(2);

          h_chg_reco[ibin][ibin2]->Draw("e0 same");
          h_chg_reco_ud_scaled[ibin][ibin2]->Draw("e0 same");
          h_chg_reco_gubdb[ibin][ibin2]->Draw("e0 same");
       
          TString tmp1 = cent_tag_pp_MC[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.85, tmp1);
        
          TString tmp = trk_tag[ibin2];
          tx = new TLatex(); tx->SetTextSize(.08);
          tx->DrawLatexNDC(0.15,0.65, tmp);
        
          if(do_fit==true){
            h_chg_reco[ibin][ibin2]->Fit(f_chg_qg_data[ibin][ibin2],"Q N M R","sames",-2.,2.);
          
            if(ibin==0) c_chg_qg_MC_reco[ibin2]->cd(6);
            else c_chg_qg_MC_reco[ibin2]->cd(11-ibin);   
            h_chg_reco_ratio_2param[ibin][ibin2] = (TH1D*)h_chg_reco[ibin][ibin2]->Clone("h_chg_reco_ratio");
            h_chg_reco_ratio_2param[ibin][ibin2]->Divide(f_chg_qg_data[ibin][ibin2]);
            h_chg_reco_ratio_2param[ibin][ibin2]->Rebin(10); h_chg_reco_ratio_2param[ibin][ibin2]->Scale(1./10.);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetTitle("reco jet charge");
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetTitle("Reco / Fit");
            h_cosm(h_chg_reco_ratio_2param[ibin][ibin2]);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,2.);
            h_chg_reco_ratio_2param[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerColor(kBlack);
            h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerStyle(1); h_chg_reco_ratio_2param[ibin][ibin2]->Draw("e0 same");
            drawline(-2.,1.,2.,1.);
            sprintf(g_frac, "g+others (MC/Fit) %.3f",(1.-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)))/gluon_ubar_dbar[ibin]);
            sprintf(g_error, "#pm %.3f",f_chg_qg_data[ibin][ibin2]->GetParError(0));
            sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_qg_data[ibin][ibin2]->GetChisquare()/f_chg_qg_data[ibin][ibin2]->GetNDF());
            sprintf(up_frac, "up+down (MC/Fit) %.3f",(f_chg_qg_data[ibin][ibin2]->GetParameter(0)/(up[ibin]+down[ibin])));
            tx1 = new TLatex(); tx1->SetTextSize(.07);
            tx1->DrawLatexNDC(0.15,0.75, g_frac);
            tx3 = new TLatex(); tx3->SetTextSize(.06);
            tx3->DrawLatexNDC(0.65,0.7, g_error);
            tx4 = new TLatex(); tx4->SetTextSize(.07);
            tx4->DrawLatexNDC(0.3,0.4, up_frac);
            tx2 = new TLatex(); tx2->SetTextSize(.07);
            tx2->DrawLatexNDC(0.3,0.2, chi2ndf);
          }
        }
      }
    }

  //////2 param kappa MC fitting

    if(do_kappa_fitting){
      TCanvas *c_chg_udg_MC_reco_kappa[ntrkbins];
      for(int ibin2=0;ibin2<ntrkbins;ibin2++){
        c_chg_udg_MC_reco_kappa[ibin2] = new TCanvas((TString)("c_chg_udg_MC_reco_kappa"+trk[ibin2]),(TString)("c_chg_udg_MC_reco_kappa"+trk[ibin2]),1500,600);
        c_chg_udg_MC_reco_kappa[ibin2]->Divide(5,2,0);
        gStyle->SetOptStat(0);
        for(int ibin=0;ibin<nCBins;ibin++){
          if(ibin==0) c_chg_udg_MC_reco_kappa[ibin2]->cd(1);
          else c_chg_udg_MC_reco_kappa[ibin2]->cd(6-ibin);   
          h_chg_reco_kappa[ibin][ibin2]->GetXaxis()->SetTitle("reco jet chg");
          h_chg_reco_kappa[ibin][ibin2]->GetYaxis()->SetTitle("1/N_{jets} (dN/dQ) [1/e]");
          h_cosm(h_chg_reco_kappa[ibin][ibin2]);
          h_chg_reco_kappa[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
          //h_chg_reco_kappa[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);

          h_chg_reco_kappa[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_kappa[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_reco_kappa[ibin][ibin2]->SetMarkerStyle(2);
          h_chg_reco_gubdb_kappa[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_gubdb_kappa[ibin][ibin2]->SetMarkerColor(kRed); h_chg_reco_gubdb_kappa[ibin][ibin2]->SetMarkerStyle(2);
          h_chg_reco_kappa_ud_scaled[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_kappa_ud_scaled[ibin][ibin2]->SetMarkerColor(kBlue); h_chg_reco_kappa_ud_scaled[ibin][ibin2]->SetMarkerStyle(2);

          h_chg_reco_kappa[ibin][ibin2]->Draw("e0 same");
          h_chg_reco_gubdb_kappa[ibin][ibin2]->Draw("e0 same");
          h_chg_reco_kappa_ud_scaled[ibin][ibin2]->Draw("e0 same");

          TString tmp1 = cent_tag_pp_MC[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.85, tmp1);
        
          TString tmp = kappa_tag[ibin2];
          tx = new TLatex(); tx->SetTextSize(.1);
          tx->DrawLatexNDC(0.15,0.65, tmp);
        
          if(do_fit==true){
            f_chg_qg_kappa[ibin][ibin2]->SetLineColor(kBlack);
            h_chg_reco_kappa[ibin][ibin2]->Fit(f_chg_qg_kappa[ibin][ibin2],"Q N M R","sames",-2.,2.);
          
            if(ibin==0) c_chg_udg_MC_reco_kappa[ibin2]->cd(6);
            else c_chg_udg_MC_reco_kappa[ibin2]->cd(11-ibin);   
            h_chg_reco_ratio_kappa[ibin][ibin2] = (TH1D*)h_chg_reco_kappa[ibin][ibin2]->Clone("h_chg_reco_ratio");
            h_chg_reco_ratio_kappa[ibin][ibin2]->Divide(f_chg_qg_kappa[ibin][ibin2]);
            h_chg_reco_ratio_kappa[ibin][ibin2]->Rebin(20); h_chg_reco_ratio_kappa[ibin][ibin2]->Scale(1./20.);
            h_chg_reco_ratio_kappa[ibin][ibin2]->GetXaxis()->SetTitle("reco jet charge");
            h_chg_reco_ratio_kappa[ibin][ibin2]->GetYaxis()->SetTitle("Reco / Fit");
            h_cosm(h_chg_reco_ratio_kappa[ibin][ibin2]); 
            h_chg_reco_ratio_kappa[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
            h_chg_reco_ratio_kappa[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,2.);
            h_chg_reco_ratio_kappa[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio_kappa[ibin][ibin2]->SetMarkerColor(kBlack);
            h_chg_reco_ratio_kappa[ibin][ibin2]->SetMarkerStyle(1); h_chg_reco_ratio_kappa[ibin][ibin2]->Draw("e0 same");
            drawline(-2.,1.,2.,1.);
            sprintf(g_frac, "g+others (MC/Fit) %.3f",(1.-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)))/gluon_ubar_dbar[ibin]);
            sprintf(g_error, "#pm %.3f",f_chg_qg_kappa[ibin][ibin2]->GetParError(0));
            sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_qg_kappa[ibin][ibin2]->GetChisquare()/f_chg_qg_kappa[ibin][ibin2]->GetNDF());
            sprintf(up_frac, "up+down (MC/Fit) %.3f",(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)/(up[ibin]+down[ibin])));
            tx1 = new TLatex(); tx1->SetTextSize(.07);
            tx1->DrawLatexNDC(0.15,0.75, g_frac);
            tx3 = new TLatex(); tx3->SetTextSize(.06);
            tx3->DrawLatexNDC(0.65,0.7, g_error);
            tx4 = new TLatex(); tx4->SetTextSize(.07);
            tx4->DrawLatexNDC(0.3,0.4, up_frac);
            tx2 = new TLatex(); tx2->SetTextSize(.07);
            tx2->DrawLatexNDC(0.3,0.2, chi2ndf);
          }
        }
      }
    }      

  ///////////gen//////////

    if(do_gen){
      TCanvas *c_chg_udg_MC_ref_rebin[ntrkbins];

      for(int ibin2=0;ibin2<1;ibin2++){

        c_chg_udg_MC_ref_rebin[ibin2] = new TCanvas((TString)("c_chg_udg_MC_ref_rebin"+trk[ibin2]),(TString)("c_chg_udg_MC_ref_rebin"+trk[ibin2]),1500,350);

        drawchgPad();
        pad_titlec->cd();
        title_tex->DrawLatexNDC(0.05,0.1, cms);
        title_tex->DrawLatexNDC(0.4,0.3, gen_tex);  
        title_tex->DrawLatexNDC(0.8,0.3, trk_tag[ibin2]);  
        pad_labelc->cd();
        ylabel_tex->DrawLatexNDC(0.7,0.2, y_label);

        for(int ibin=0;ibin<nCBins;ibin++){
      
          if(ibin==0) pad1c->cd(1);
          else pad1c->cd(6-ibin);   
          h_chg_ref[ibin][ibin2]->Rebin(20);
          //h_chg_reco_cop[ibin][ibin2]->Rebin(10);
          //gPad->SetLogy(); 
          h_chg_ref[ibin][ibin2]->GetXaxis()->SetTitle("jet charge");
          h_cosm(h_chg_ref[ibin][ibin2]);
          h_chg_ref[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
          h_chg_ref[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.25);

          h_chg_ref[ibin][ibin2]->SetLineColor(kGray); h_chg_ref[ibin][ibin2]->SetFillColor(kGray); h_chg_ref[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_ref[ibin][ibin2]->SetMarkerStyle(4);
          h_chg_ref[ibin][ibin2]->Draw("hist same");

          h_chg_reco_cop[ibin][ibin2]->Draw("e0 same");

          TString tmp1 = cent_tag_pp_MC[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.85, tmp1);
          
          if(ibin==0){
            TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
            legendx ->SetLineColor(kWhite);
            legendx ->AddEntry(h_chg_ref[0][0], "gen", "f");
            legendx ->AddEntry(h_chg_reco_cop[0][0], "reco", "lep");
            legendx ->Draw("same");
          }
        }
      }
    }  

    if(do_leading_track){
    ///////highest trk chg plot

      TCanvas *c_chg_highest_rebin = new TCanvas("c_chg_highest_rebin","c_chg_highest_rebin",1200,350);

      TPad *pad_titles = new TPad("pad_title", "",0.,0.92,1.,1.);
      TPad *pad_labels = new TPad("pad_label", "",0.,0.,0.02,0.95);
      TPad *pad1s = new TPad("pad1", "",0.02,0.02,1.,0.98);
      pad1s->Divide(4,1,0);  
      pad1s->Draw();
      pad_titles->Draw();
      pad_labels->Draw();
      pad_titles->cd();
      //title_tex->DrawLatexNDC(0.4,0.1, data);
      title_tex->DrawLatexNDC(0.05,0.1, cms);
      //title_tex->DrawLatexNDC(0.4,0.3, MC);  
      title_tex->DrawLatexNDC(0.8,0.3, trk_tag[0]);  
      pad_labels->cd();
      ylabel_tex->DrawLatexNDC(0.7,0.2, y_label);

      for(int ibin=1;ibin<nCBins;ibin++){

        if(ibin==0) pad1s->cd(1);
        else pad1s->cd(5-ibin);   

        h_chg_highest[ibin]->Rebin(2);
        h_chg_highest_2[ibin]->Rebin(2);

        h_chg_highest[ibin]->Scale(1./h_chg_highest[ibin]->Integral());
        h_chg_highest_2[ibin]->Scale(1./h_chg_highest_2[ibin]->Integral());
        
        //gPad->SetLogy(); 
        h_chg_highest[ibin]->GetXaxis()->SetTitle("jet charge");
        h_cosm(h_chg_highest[ibin]); 
        h_chg_highest[ibin]->GetXaxis()->SetRangeUser(-1.5,1.5);
        h_chg_highest[ibin]->GetYaxis()->SetRangeUser(0.,0.14);

        h_chg_highest[ibin]->SetLineColor(kRed); h_chg_highest[ibin]->SetMarkerColor(kRed); h_chg_highest[ibin]->SetMarkerStyle(20);
        h_chg_highest[ibin]->Draw("e0 same");
        h_chg_highest_2[ibin]->SetLineColor(kBlue); h_chg_highest_2[ibin]->SetMarkerColor(kBlue); h_chg_highest_2[ibin]->SetMarkerStyle(20);
        h_chg_highest_2[ibin]->Draw("e0 same");

        TString tmp1 = cent_tag_pp_MC[ibin];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.15,0.85, tmp1);
        
        if(ibin==4){
          TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
          legendx ->SetLineColor(kWhite);
          legendx ->AddEntry(h_chg_highest[4], "n = 1", "lepf");
          legendx ->AddEntry(h_chg_highest_2[4], "n = 2", "lepf");

          legendx ->Draw("same");
        }
      }
    }

///////////final reco plot//////////

    if(do_ptcut_fitting){
      TCanvas *c_chg_udg_MC_reco_rebin[ntrkbins];

      for(int ibin2=0;ibin2<1;ibin2++){
        c_chg_udg_MC_reco_rebin[ibin2] = new TCanvas((TString)("c_chg_udg_MC_reco_rebin"+trk[ibin2]),(TString)("c_chg_udg_MC_reco_rebin"+trk[ibin2]),1500,450);

        drawchgPad();
        pad_titlec->cd();
        //title_tex->DrawLatexNDC(0.4,0.1, data);
        title_tex->DrawLatexNDC(0.05,0.1, cms);
        title_tex->DrawLatexNDC(0.8,0.3, trk_tag[ibin2]);  
        pad_labelc->cd();
        ylabel_tex->DrawLatexNDC(0.7,0.5, y_label);
        ylabel_tex->DrawLatexNDC(0.7,0.06, datafit);

        for(int ibin=0;ibin<nCBins;ibin++){
      
          if(ibin==0) pad1c->cd(1);
          else pad1c->cd(6-ibin);   
          
          //gPad->SetLogy(); 
          h_chg_reco_cop[ibin][ibin2]->GetXaxis()->SetTitle("jet charge");
          h_cosm(h_chg_reco_cop[ibin][ibin2]);
          h_chg_reco_cop[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
          h_chg_reco_cop[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.22);

          h_chg_reco_cop[ibin][ibin2]->Draw("e0 same");

          //scaling according to MC fractions
          h_chg_reco_up_cop[ibin][ibin2]->Scale(f_chg_qg_data[ibin][ibin2]->GetParameter(0)*u_ud[ibin]);
          h_chg_reco_down_cop[ibin][ibin2]->Scale(f_chg_qg_data[ibin][ibin2]->GetParameter(0)*d_ud[ibin]); 
          h_chg_reco_g_cop[ibin][ibin2]->Scale((1.-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)))*gluon_ubdbcsb[ibin]);
          h_chg_reco_others_cop[ibin][ibin2]->Scale((1.-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)))*(1-gluon_ubdbcsb[ibin]));

          h_chg_reco_stk[ibin][ibin2] = new THStack((TString)("h_chg_reco_stk"+cent[ibin]+"_"+trk[ibin2]),"");
          h_chg_reco_stk[ibin][ibin2]->Add(h_chg_reco_others_cop[ibin][ibin2]);
          h_chg_reco_stk[ibin][ibin2]->Add(h_chg_reco_down_cop[ibin][ibin2]);
          h_chg_reco_stk[ibin][ibin2]->Add(h_chg_reco_up_cop[ibin][ibin2]);
          h_chg_reco_stk[ibin][ibin2]->Add(h_chg_reco_g_cop[ibin][ibin2]);

          h_chg_reco_stk[ibin][ibin2]->Draw("hist same");      
          h_chg_reco_cop[ibin][ibin2]->Draw("e0 same");

          TString tmp1 = cent_tag_pp_MC[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.85, tmp1);
          
          if(ibin==0){
            TLegend *legendx = new TLegend(0.15,0.7,0.5,0.8);
            legendx ->SetLineColor(kWhite);
            legendx ->AddEntry(h_chg_reco_cop[0][ibin2], "MC (reco)", "lep");
            //legendx ->AddEntry(h_chg_data_jeup[0][ibin2], "Data (jet p_{T} > 126 GeV)", "lep");
            //legendx ->AddEntry(h_chg_data_jedn[0][ibin2], "Data (jet p_{T} > 114 GeV)", "lep");
            legendx ->Draw("same");

            TLegend *legendy = new TLegend(0.62,0.62,0.98,0.98);
            legendy ->SetLineColor(kWhite);
            legendy ->AddEntry((TObject*)0, "Fitting results", "");
            legendy ->AddEntry(h_chg_reco_g_cop[0][ibin2], "Gluon", "lepf");
            legendy ->AddEntry(h_chg_reco_up_cop[0][ibin2], "Up quark", "lepf");
            legendy ->AddEntry(h_chg_reco_down_cop[0][ibin2], "Down quark", "lepf");
            legendy ->AddEntry(h_chg_reco_others_cop[0][ibin2], "Others", "lepf");
            legendy ->Draw("same");
          }

          if(do_fit==true){  
        
            if(ibin==0) pad2c->cd(1);
            else pad2c->cd(6-ibin);   

            h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetTitle("");
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetTitle("");
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetLabelSize(0.13);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->CenterTitle();
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->CenterTitle();
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetTitleOffset(0.7);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetLabelSize(0.06);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetRangeUser(0.5,1.5);
            h_chg_reco_ratio_2param[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerColor(kBlack);
            h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerStyle(20); h_chg_reco_ratio_2param[ibin][ibin2]->Draw("e0 same");
            tl1->Draw("same");
          }
        }
      }
    }

    for(int ibin=0; ibin<nCBins; ibin++){
      for(int ibin2=0;ibin2<ntrkbins-1;ibin2++){
        if(ibin2==0){

          if(do_ptcut_fitting){
              up_fit[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParameter(0)*u_ud[ibin]; up_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
              down_fit[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParameter(0)*d_ud[ibin]; down_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
              quark_fit[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParameter(0)+ubdbcsb[ibin]; quark_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
              quark_fit[ibin][ibin2] = quark_fit[ibin][ibin2]*quark_corr_factor[ibin][ibin2];
              gluon_fit[ibin][ibin2] = 1.-(quark_fit[ibin][ibin2]); gluon_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(1);
          }

          if(do_kappa_fitting){
              up_fit_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)*u_ud[ibin]; up_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParError(0);
              down_fit_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)*d_ud[ibin]; down_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParError(0);
              quark_fit_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)+ubdbcsb[ibin]; quark_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParError(0);
              quark_fit_kappa[ibin][ibin2] = quark_fit_kappa[ibin][ibin2]*quark_corr_factor[ibin][ibin2];
              gluon_fit_kappa[ibin][ibin2] = 1.-(quark_fit_kappa[ibin][ibin2]); gluon_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParError(0);
          }
        }
        else{

          if(do_ptcut_fitting){
              up_fit[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParameter(0)*u_ud[ibin]; up_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParError(0);
              down_fit[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParameter(0)*d_ud[ibin]; down_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParError(0);
              quark_fit[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParameter(0)+ubdbcsb[ibin]; quark_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParError(0);       
              quark_fit[ibin][ibin2] = quark_fit[ibin][ibin2]*quark_corr_factor[ibin][ibin2];
              gluon_fit[ibin][ibin2] = 1.-(quark_fit[ibin][ibin2]); gluon_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParError(0);
          } 

          if(do_kappa_fitting){
              up_fit_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParameter(0)*u_ud[ibin]; up_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParError(0);
              down_fit_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParameter(0)*d_ud[ibin]; down_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParError(0);
              quark_fit_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParameter(0)+ubdbcsb[ibin]; quark_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParError(0);       
              quark_fit_kappa[ibin][ibin2] = quark_fit_kappa[ibin][ibin2]*quark_corr_factor[ibin][ibin2];
              gluon_fit_kappa[ibin][ibin2] = 1.-(quark_fit_kappa[ibin][ibin2]); gluon_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParError(0);
          }
        }

        gluon_reco_fitsys[ibin][ibin2] = gluon_fit[ibin][ibin2]; 
        quark_reco_fitsys[ibin][ibin2] = quark_fit[ibin][ibin2];

        quark_recofit_ratio[ibin][ibin2] = quark_reco_fitsys[ibin][ibin2]/quark_reco_fitsys[0][ibin2];
        gluon_recofit_ratio[ibin][ibin2] = gluon_reco_fitsys[ibin][ibin2]/gluon_reco_fitsys[0][ibin2];

      }

      //MC graphs
      if(do_ptcut_fitting){
          up_reco[ibin] = new TGraphErrors(ntrkbins-1,x_mean,up_fit[ibin],x_mean_err,up_fiterr[ibin]);
          up_reco[ibin]->SetName((TString)("up_reco_cent"+cent[ibin]));
          down_reco[ibin] = new TGraphErrors(ntrkbins-1,x_mean,down_fit[ibin],x_mean_err,down_fiterr[ibin]);
          down_reco[ibin]->SetName((TString)("down_reco_cent"+cent[ibin]));
          gluon_reco[ibin] = new TGraphErrors(ntrkbins-1,x_mean,gluon_fit[ibin],x_mean_err,gluon_fiterr[ibin]);
          gluon_reco[ibin]->SetName((TString)("gluon_reco_cent"+cent[ibin]));
          quark_reco[ibin] = new TGraphErrors(ntrkbins-1,x_mean,quark_fit[ibin],x_mean_err,quark_fiterr[ibin]);
          quark_reco[ibin]->SetName((TString)("quark_reco_cent"+cent[ibin]));
      }

      if(do_kappa_fitting){
          up_reco_kappa[ibin] = new TGraphErrors(ntrkbins-1,x_mean,up_fit_kappa[ibin],x_mean_err,up_fiterr_kappa[ibin]);
          up_reco_kappa[ibin]->SetName((TString)("up_reco_kappa_cent"+cent[ibin]));
          down_reco_kappa[ibin] = new TGraphErrors(ntrkbins-1,x_mean,down_fit_kappa[ibin],x_mean_err,down_fiterr_kappa[ibin]);
          down_reco_kappa[ibin]->SetName((TString)("down_reco_kappa_cent"+cent[ibin]));
          gluon_reco_kappa[ibin] = new TGraphErrors(ntrkbins-1,x_mean,gluon_fit_kappa[ibin],x_mean_err,gluon_fiterr_kappa[ibin]);
          gluon_reco_kappa[ibin]->SetName((TString)("gluon_reco_kappa_cent"+cent[ibin]));
          quark_reco_kappa[ibin] = new TGraphErrors(ntrkbins-1,x_mean,quark_fit_kappa[ibin],x_mean_err,quark_fiterr_kappa[ibin]);
          quark_reco_kappa[ibin]->SetName((TString)("quark_reco_kappa_cent"+cent[ibin]));
      }

      q_reco_0trk[ibin] = new TGraphErrors(ntrkbins-1,x_mean,quark_0trk_int[ibin],x_mean_err,x_mean_err);
      q_reco_0trk[ibin]->SetName((TString)("q_reco_0trk_cent"+cent[ibin]));
    }

    TCanvas *c_reco_graph = new TCanvas("c_reco_graph","c_reco_graph",1500,350);
    drawgraphPad();

    pad_titleg->cd();
    title_tex->DrawLatexNDC(0.5,0.1, pythia);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    pad_labelg->cd();
    chg_tex->DrawLatexNDC(0.7,0.4, fraction);

    for(int ibin=0; ibin<nCBins; ibin++){
      if(ibin==0) pada->cd(1);
      else pada->cd(6-ibin);
      h_dummy->Draw("same");

      if(do_ptcut_fitting){
          gluon_reco[ibin]->SetMarkerColor(kRed); gluon_reco[ibin]->SetMarkerStyle(5); gluon_reco[ibin]->SetMarkerSize(1.2);
          gluon_reco[ibin]->Draw("PL");
          quark_reco[ibin]->SetMarkerColor(kBlue); quark_reco[ibin]->SetMarkerStyle(5); quark_reco[ibin]->SetMarkerSize(1.2);
          quark_reco[ibin]->Draw("PL");  
      }

      if(do_kappa_fitting){
          gluon_reco_kappa[ibin]->SetMarkerColor(kRed); gluon_reco_kappa[ibin]->SetMarkerStyle(4); gluon_reco_kappa[ibin]->SetMarkerSize(1.2);
          gluon_reco_kappa[ibin]->Draw("PL");
          quark_reco_kappa[ibin]->SetMarkerColor(kBlue); quark_reco_kappa[ibin]->SetMarkerStyle(4); quark_reco_kappa[ibin]->SetMarkerSize(1.2);
          quark_reco_kappa[ibin]->Draw("PL");
      }  
      tl_g[ibin]->Draw("same"); tl_q[ibin]->Draw("same");

          TString tmp1 = cent_tag_pp[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.25,0.85, tmp1);

      if(ibin==0){
        TLegend *legend = new TLegend(0.35,0.15,0.75,0.4);
        legend ->SetLineColor(kWhite);
        if(do_ptcut_fitting){
          legend ->AddEntry(quark_reco[ibin], "quark (pT cut)", "lepf");
          legend ->AddEntry(gluon_reco[ibin], "gluon (pT cut)", "lepf");
        }
        if(do_kappa_fitting){
          legend ->AddEntry(quark_reco_kappa[ibin], "quark (kappa)", "lepf");
          legend ->AddEntry(gluon_reco_kappa[ibin], "gluon (kappa)", "lepf");
        }
        legend ->Draw("same");
      }    
      if(ibin==4){
        if(do_kappa_fitting){
          TLegend *legenda = new TLegend(0.35,0.15,0.75,0.4);
          legenda ->SetLineColor(kWhite);
          legenda ->AddEntry((TObject*)0, "x1 : pT > 2 (#kappa=0.5); #kappa=0.3 (pT > 2)", "");
          legenda ->AddEntry((TObject*)0, "x2 : pT > 4 (#kappa=0.5); #kappa=0.5 (pT > 2)", "");
          legenda ->AddEntry((TObject*)0, "x3 : pT > 5 (#kappa=0.5); #kappa=0.7 (pT > 2)", "");
          legenda ->Draw("same");
        }
      }
    }

    TCanvas *c_reco_q_graph = new TCanvas("c_reco_q_graph","c_reco_q_graph",1800,350);
    drawgraphPad();

    pad_titleg->cd();
    title_tex->DrawLatexNDC(0.5,0.1, pythia);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    pad_labelg->cd();
    chg_tex->DrawLatexNDC(0.7,0.4, fraction);

    for(int ibin=0; ibin<nCBins; ibin++){
      if(ibin==0) pada->cd(1);
      else pada->cd(6-ibin);
      h_dummy2->Draw("same");

      if(do_ptcut_fitting){
        up_reco[ibin]->SetMarkerColor(kAzure+2); up_reco[ibin]->SetMarkerStyle(22); up_reco[ibin]->SetMarkerSize(1.5);
        up_reco[ibin]->Draw("PL");
        down_reco[ibin]->SetMarkerColor(kGreen-2); down_reco[ibin]->SetMarkerStyle(23); down_reco[ibin]->SetMarkerSize(1.5);
        down_reco[ibin]->Draw("PL");
      }

      if(do_kappa_fitting){
        up_reco_kappa[ibin]->SetMarkerColor(kAzure+2); up_reco_kappa[ibin]->SetMarkerStyle(22); up_reco_kappa[ibin]->SetMarkerSize(1.5);
        up_reco_kappa[ibin]->Draw("PL");
        down_reco_kappa[ibin]->SetMarkerColor(kGreen-2); down_reco_kappa[ibin]->SetMarkerStyle(23); down_reco_kappa[ibin]->SetMarkerSize(1.5);
        down_reco_kappa[ibin]->Draw("PL");
      }

      tl_u[ibin]->Draw("same"); tl_d[ibin]->Draw("same"); //tl_ubar[ibin]->Draw("same"); tl_dbar[ibin]->Draw("same"); tl_c[ibin]->Draw("same"); tl_s[ibin]->Draw("same");tl_b[ibin]->Draw("same");

      TString tmp1 = cent_tag_pp_MC[ibin];
      tx1 = new TLatex(); tx1->SetTextSize(.1);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);

      if(ibin==0){
        TLegend *legend = new TLegend(0.15,0.6,0.4,0.9);
        legend ->SetLineColor(kWhite);
        if(do_ptcut_fitting){
          legend ->AddEntry(up_reco[ibin], "up", "lepf");
          legend ->AddEntry(down_reco[ibin], "down", "lepf");
        }
        if(do_kappa_fitting){
          legend ->AddEntry(up_reco_kappa[ibin], "up", "lepf");
          legend ->AddEntry(down_reco_kappa[ibin], "down", "lepf");
        }
        legend ->Draw("same");
      }    
    } 
  }

////data fitting  

  if(do_data){
    if(do_ptcut_fitting){
      //////2 parameter fitting data
      TCanvas *c_chg_qg_MC_data[ntrkbins];

      for(int ibin2=0;ibin2<ntrkbins;ibin2++){
        c_chg_qg_MC_data[ibin2] = new TCanvas((TString)("c_chg_qg_MC_data"+trk[ibin2]),(TString)("c_chg_qg_MC_data"+trk[ibin2]),1800,600);
        c_chg_qg_MC_data[ibin2]->Divide(5,2,0);
        //gStyle->SetOptStat(0);
        gStyle->SetOptFit(1);
        for(int ibin=0;ibin<nCBins;ibin++){
      
          if(ibin==0) c_chg_qg_MC_data[ibin2]->cd(1);
          else c_chg_qg_MC_data[ibin2]->cd(6-ibin);   
          h_chg_data[ibin][ibin2]->GetXaxis()->SetTitle("data jet chg");
          h_chg_data[ibin][ibin2]->GetYaxis()->SetTitle("# of jets");
          h_cosm(h_chg_data[ibin][ibin2]);
          h_chg_data[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
          //h_chg_data[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);
          h_chg_data[ibin][ibin2]->SetLineColor(kGreen); h_chg_data[ibin][ibin2]->SetMarkerColor(kGreen); h_chg_data[ibin][ibin2]->SetMarkerStyle(2);
        
          h_chg_data[ibin][ibin2]->Draw("e0 same");
          h_chg_reco_ud_scaled[ibin][ibin2]->Draw("e0 same");
          h_chg_reco_gubdb[ibin][ibin2]->Draw("e0 same");
          
          TString tmp1 = cent_tag_pp[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.85, tmp1);
        
          TString tmp = trk_tag[ibin2];
          tx = new TLatex(); tx->SetTextSize(.1);
          tx->DrawLatexNDC(0.15,0.65, tmp);

          if(do_fit==true){  
            h_chg_data[ibin][ibin2]->Fit(f_chg_qg_data[ibin][ibin2],"N M R","sames",-2.,2.);
        
            if(ibin==0) c_chg_qg_MC_data[ibin2]->cd(6);
            else c_chg_qg_MC_data[ibin2]->cd(11-ibin);   
            h_chg_reco_ratio[ibin][ibin2] = (TH1D*)h_chg_data[ibin][ibin2]->Clone("h_chg_data_ratio");
            h_chg_reco_ratio[ibin][ibin2]->Divide(f_chg_qg_data[ibin][ibin2]);
            /*
            h_chg_reco_ratio[ibin][ibin2] = (TH1D*)h_chg_reco_up[ibin][ibin2]->Clone("h_chg_data_ratio");
            h_chg_reco_ratio[ibin][ibin2]->Add(h_chg_reco_down[ibin][ibin2]);
            h_chg_reco_ratio[ibin][ibin2]->Add(h_chg_reco_g[ibin][ibin2]);
            h_chg_reco_ratio[ibin][ibin2]->Add(h_chg_reco_others[ibin][ibin2]);
            h_chg_reco_ratio[ibin][ibin2]->Divide(h_chg_data[ibin][ibin2]);
            */
            //h_chg_reco_ratio[ibin][ibin2] = (TH1D*)h_chg_data[ibin][ibin2]->Clone("h_chg_data_ratio");
            //h_chg_reco_ratio[ibin][ibin2]->Divide(f_chg_qg_data[ibin][ibin2]);
            h_chg_reco_ratio[ibin][ibin2]->Rebin(10); h_chg_reco_ratio[ibin][ibin2]->Scale(1./10.);
            h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetTitle("data jet charge");
            h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetTitle("Data / Fit");
            h_cosm(h_chg_reco_ratio[ibin][ibin2]);  
            h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.8,1.8);
            h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetRangeUser(0.8,1.2);
            h_chg_reco_ratio[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_ratio[ibin][ibin2]->SetMarkerColor(kRed);
            h_chg_reco_ratio[ibin][ibin2]->SetMarkerStyle(21); h_chg_reco_ratio[ibin][ibin2]->Draw("e0 same");
            tl1->Draw("same");
            sprintf(g_frac, "g+others %.3f",(1.-f_chg_qg_data[ibin][ibin2]->GetParameter(0)));
            sprintf(g_error, "#pm %.3f",f_chg_qg_data[ibin][ibin2]->GetParError(0));
            sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_qg_data[ibin][ibin2]->GetChisquare()/f_chg_qg_data[ibin][ibin2]->GetNDF());
            sprintf(up_frac, "up+down %.3f",f_chg_qg_data[ibin][ibin2]->GetParameter(0));
            tx1 = new TLatex(); tx1->SetTextSize(.08);
            tx1->DrawLatexNDC(0.25,0.85, g_frac);
            tx3 = new TLatex(); tx3->SetTextSize(.06);
            tx3->DrawLatexNDC(0.65,0.8, g_error);
            tx4 = new TLatex(); tx4->SetTextSize(.07);
            tx4->DrawLatexNDC(0.3,0.4, up_frac);
            tx2 = new TLatex(); tx2->SetTextSize(.07);
            tx2->DrawLatexNDC(0.3,0.2, chi2ndf);
          }
        }
      }

      ////pull distribution

      if(do_ptcut_fitting){ 
        Double_t pull;
        TH1D *h_pull[nCBins][ntrkbins];
        TF1 *f_pull[nCBins][ntrkbins];

        TCanvas *c_pull[ntrkbins];

        for(int ibin2=0;ibin2<1;ibin2++){
          c_pull[ibin2] = new TCanvas((TString)("c_pull"+trk[ibin2]),(TString)("c_pull"+trk[ibin2]),1500,300);
          c_pull[ibin2]->Divide(5,1,0);
          gStyle->SetOptFit(1);

          for(int ibin=0;ibin<nCBins;ibin++){  

            sprintf(saythis,"h_pull_cent%d_trk%d",ibin,ibin2);
            h_pull[ibin][ibin2] = new TH1D(saythis,saythis,100,-0.1,0.1);
            h_pull[ibin][ibin2]->Sumw2(); 

            sprintf(saythis,"f_pull_cent%d_trk%d",ibin,ibin2);
            f_pull[ibin][ibin2] = new TF1(saythis,"gaus",-0.05,0.05);

            pull = 0.;
            for(int ibin3=0;ibin3<h_chg_data[ibin][ibin2]->GetNbinsX();ibin3++){
              pull = (h_chg_data[ibin][ibin2]->GetBinContent(ibin3+1) - f_chg_qg_data[ibin][ibin2]->Eval(h_chg_data[ibin][ibin2]->GetBinCenter(ibin3+1)));
              h_pull[ibin][ibin2]->Fill(pull/f_chg_qg_data[ibin][ibin2]->GetParError(0));
            }

            if(ibin==0) c_pull[ibin2]->cd(1);
            else c_pull[ibin2]->cd(6-ibin);  
            h_cosm(h_pull[ibin][ibin2]);
            h_pull[ibin][ibin2]->Draw("e0 same");
            h_pull[ibin][ibin2]->Fit(f_pull[ibin][ibin2],"Q M R","sames",-0.05,0.05);

          }
        }
      }
     
      //////vary trk eff fitting
      
      TH1D* h_chi2ndf[nCBins][ntrkbins];
      TH1D* h_least_chi2ndf[ntrkbins];

      for(int ibin2=0; ibin2<ntrkbins; ibin2++){
          h_chi2ndf[0][ibin2] = new TH1D(Form("h_chi2ndf_%d_%d",0,ibin2),"",11,-6.,4.);
          h_chi2ndf[1][ibin2] = new TH1D(Form("h_chi2ndf_%d_%d",1,ibin2),"",11,-6.,4.);
          h_chi2ndf[2][ibin2] = new TH1D(Form("h_chi2ndf_%d_%d",2,ibin2),"",11,-5.,5.);
          h_chi2ndf[3][ibin2] = new TH1D(Form("h_chi2ndf_%d_%d",3,ibin2),"",11,-5.,5.);
          h_chi2ndf[4][ibin2] = new TH1D(Form("h_chi2ndf_%d_%d",4,ibin2),"",11,-5.,5.);

          h_least_chi2ndf[ibin2] = new TH1D(Form("h_least_chi2ndf_%d",ibin2),"",5,0.,5.);
      } 

      Double_t chi2NDF[nCBins][ntrkbins][trkeffbins+1];

      //////two parameter fitting

      TCanvas *c_chg_qg_MC_data_trkeff[ntrkbins];

      for(int ibin2=0;ibin2<ntrkbins;ibin2++){
        c_chg_qg_MC_data_trkeff[ibin2] = new TCanvas((TString)("c_chg_qg_MC_data_trkeff"+trk[ibin2]),(TString)("c_chg_qg_MC_data_trkeff"+trk[ibin2]),1500,300);
        c_chg_qg_MC_data_trkeff[ibin2]->Divide(5,1,0);
        for(int ibin=0;ibin<nCBins;ibin++){
      
          if(ibin==0) c_chg_qg_MC_data_trkeff[ibin2]->cd(1);
          else c_chg_qg_MC_data_trkeff[ibin2]->cd(6-ibin);   

          if(do_fit==true){  
            for(int ibin3=0; ibin3<trkeffbins;ibin3++){
               h_chg_data_trkeff[ibin][ibin2][ibin3]->Fit(f_chg_qg_data[ibin][ibin2],"Q N M R","sames",-2.,2.);

               if(ibin3==7){
                   quark_fit_high[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParameter(0)+(1-gluon_ubdbcsb[ibin])*(1-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)));
                   quark_fit_high[ibin][ibin2] = quark_fit_high[ibin][ibin2]*quark_corr_factor[ibin][ibin2];
                   quark_fit_high[ibin][ibin2] = (quark_fit_high[ibin][ibin2]*data_no0trk_int[ibin][ibin2] + quark_data_0trk_int[ibin][ibin2])/data_int[ibin];
                   gluon_fit_high[ibin][ibin2] = 1.-(quark_fit_high[ibin][ibin2]);
               }

               if(ibin3==3){
                   quark_fit_low[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParameter(0)+(1-gluon_ubdbcsb[ibin])*(1-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)));
                   quark_fit_low[ibin][ibin2] = quark_fit_low[ibin][ibin2]*quark_corr_factor[ibin][ibin2];
                   quark_fit_low[ibin][ibin2] = (quark_fit_low[ibin][ibin2]*data_no0trk_int[ibin][ibin2] + quark_data_0trk_int[ibin][ibin2])/data_int[ibin];
                   gluon_fit_low[ibin][ibin2] = 1.-(quark_fit_low[ibin][ibin2]);
               }

               h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3] = (TH1D*)h_chg_data_trkeff[ibin][ibin2][ibin3]->Clone("h_chg_data_ratio_trkeff");
               h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->Divide(f_chg_qg_data[ibin][ibin2]);
               h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->Rebin(10); h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->Scale(1./10.);
            
               chi2NDF[ibin][ibin2][ibin3] = f_chg_qg_data[ibin][ibin2]->GetChisquare()/f_chg_qg_data[ibin][ibin2]->GetNDF();
               h_chi2ndf[ibin][ibin2]->SetBinContent(ibin3+1,chi2NDF[ibin][ibin2][ibin3]);
            }

            h_chg_data[ibin][ibin2]->Fit(f_chg_qg_data[ibin][ibin2],"Q N M R","sames",-2.,2.);

            if(ibin==0) c_chg_qg_MC_data_trkeff[ibin2]->cd(1);
            else c_chg_qg_MC_data_trkeff[ibin2]->cd(6-ibin);   

            for(int ibin3=0; ibin3<trkeffbins;ibin3++){
                if(ibin3<5){
                  h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->SetLineColor(kCyan-5+ibin3); h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->SetMarkerColor(kCyan-5+ibin3);
                }
                else{
                  h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->SetLineColor(kGreen-5+ibin3); h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->SetMarkerColor(kGreen-5+ibin3);
                }
                h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->SetMarkerStyle(20); h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->Draw("HIST L same");        
            }
            h_chg_reco_ratio[ibin][ibin2]->Draw("HIST L same");
           
            if(ibin==0){
                TLegend *legendz = new TLegend(0.65,0.55,0.98,0.98);
                legendz ->SetLineColor(kWhite);
                legendz ->AddEntry(h_chg_reco_ratio[0][0], "Data (nominal)", "l");
                legendz ->Draw("same");
            }
            if(ibin==4){
                TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
                legendx ->SetLineColor(kWhite);
                for(int ibin3=0; ibin3<5;ibin3++){
                  legendx ->AddEntry(h_chg_reco_ratio_trkeff[0][0][ibin3], trkeffvary[ibin3], "l");
                }
                legendx ->Draw("same");
            }
            if(ibin==3){
                TLegend *legendy = new TLegend(0.65,0.55,0.98,0.98);
                legendy ->SetLineColor(kWhite);
                for(int ibin3=5; ibin3<10;ibin3++){
                  legendy ->AddEntry(h_chg_reco_ratio_trkeff[0][0][ibin3], trkeffvary[ibin3], "l");
                }
                legendy ->Draw("same");
            }

            TString tmp1 = cent_tag_pp[ibin];
            tx1 = new TLatex(); tx1->SetTextSize(.1);
            tx1->DrawLatexNDC(0.35,0.85, tmp1);

            tl1->Draw("same");
          }
        }
      }

      ///chi2/ndf plot for varying trk eff
      TCanvas *c_chg_qg_MC_data_chi2ndf;
      c_chg_qg_MC_data_chi2ndf = new TCanvas("c_chg_qg_MC_data_chi2ndf","c_chg_qg_MC_data_chi2ndf",1800,300);
      c_chg_qg_MC_data_chi2ndf->Divide(5,1);

      for(int ibin=0;ibin<nCBins;ibin++){        
        for(int ibin2=0;ibin2<4;ibin2++){
          if(ibin2==1) continue;
          //gStyle->SetOptStat(0);
          //gStyle->SetOptFit(1);
          h_least_chi2ndf[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp[ibin]);

          h_least_chi2ndf[ibin2]->SetBinContent(ibin+1,h_chi2ndf[ibin][ibin2]->GetBinContent(6));
          h_least_chi2ndf[ibin2]->SetBinError(ibin+1,0.);

          h_chi2ndf[ibin][ibin2]->GetXaxis()->SetTitle("trk. eff. vary (%)");
          h_chi2ndf[ibin][ibin2]->GetYaxis()->SetTitle("#chi^{2}/NDF");
          h_cosm(h_chi2ndf[ibin][ibin2]);
          h_chi2ndf[ibin][ibin2]->SetLineColor(ibin2+2); h_chi2ndf[ibin][ibin2]->SetMarkerColor(ibin2+2);
          h_chi2ndf[ibin][ibin2]->SetMarkerStyle(34); h_chi2ndf[ibin][ibin2]->SetMarkerSize(1.5);

          if(ibin==0) c_chg_qg_MC_data_chi2ndf->cd(1);
          else c_chg_qg_MC_data_chi2ndf->cd(6-ibin);   
          h_chi2ndf[ibin][ibin2]->Draw("p same");

          if(ibin==0){
              TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
              legendx ->SetLineColor(kWhite);
              legendx ->AddEntry(h_chi2ndf[ibin][ibin2], trk_tag[ibin2], "lepf");
              legendx ->Draw("same");
          }
        }
          TString tmp1 = cent_tag_pp[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.35,0.85, tmp1);
      }

      TCanvas *c_chg_qg_least_chi2ndf;
      c_chg_qg_least_chi2ndf = new TCanvas("c_chg_qg_least_chi2ndf","c_chg_qg_least_chi2ndf",500,500);
      c_chg_qg_least_chi2ndf->cd(1);

      TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
      legendx ->SetLineColor(kWhite);

      for(int ibin2=0;ibin2<4;ibin2++){
        if(ibin2==1) continue;

        h_least_chi2ndf[ibin2]->GetXaxis()->SetTitle("");
        h_least_chi2ndf[ibin2]->GetYaxis()->SetTitle("#chi^{2}/NDF");
        h_cosm(h_least_chi2ndf[ibin2]);
        h_least_chi2ndf[ibin2]->SetLineColor(ibin2+1); h_least_chi2ndf[ibin2]->SetMarkerColor(ibin2+1);
        h_least_chi2ndf[ibin2]->SetMarkerStyle(34); h_least_chi2ndf[ibin2]->SetMarkerSize(1.5);
        h_least_chi2ndf[ibin2]->Draw("p same");
        legendx ->AddEntry(h_least_chi2ndf[ibin2], trk_tag[ibin2], "lepf");
      }
      legendx ->Draw("same");
    }

    ////2 param kappa data fitting

    if(do_kappa_fitting){
      TCanvas *c_chg_udg_MC_data_kappa[ntrkbins];
      for(int ibin2=0;ibin2<ntrkbins;ibin2++){
        c_chg_udg_MC_data_kappa[ibin2] = new TCanvas((TString)("c_chg_udg_MC_data_kappa"+trk[ibin2]),(TString)("c_chg_udg_MC_data_kappa"+trk[ibin2]),1500,600);
        c_chg_udg_MC_data_kappa[ibin2]->Divide(5,2,0);
        gStyle->SetOptStat(0);
        for(int ibin=0;ibin<nCBins;ibin++){
          if(ibin==0) c_chg_udg_MC_data_kappa[ibin2]->cd(1);
          else c_chg_udg_MC_data_kappa[ibin2]->cd(6-ibin);   
          h_chg_data_kappa[ibin][ibin2]->GetXaxis()->SetTitle("reco jet chg");
          h_chg_data_kappa[ibin][ibin2]->GetYaxis()->SetTitle("1/N_{jets} (dN/dQ) [1/e]");
          h_cosm(h_chg_data_kappa[ibin][ibin2]); 
          h_chg_data_kappa[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
          h_chg_data_kappa[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,1.2*h_chg_reco_gubdb_kappa[ibin][ibin2]->GetMaximum());

          h_chg_data_kappa[ibin][ibin2]->SetLineColor(kBlack); h_chg_data_kappa[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_data_kappa[ibin][ibin2]->SetMarkerStyle(2);
          h_chg_reco_up_kappa[ibin][ibin2]->SetMarkerStyle(22); h_chg_reco_up_kappa[ibin][ibin2]->SetLineColor(kAzure+1); h_chg_reco_up_kappa[ibin][ibin2]->SetMarkerColor(kAzure+1);
          h_chg_reco_down_kappa[ibin][ibin2]->SetMarkerStyle(23); h_chg_reco_down_kappa[ibin][ibin2]->SetLineColor(kGreen-2); h_chg_reco_down_kappa[ibin][ibin2]->SetMarkerColor(kGreen-2);
          h_chg_reco_gubdb_kappa[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_gubdb_kappa[ibin][ibin2]->SetMarkerColor(kRed); h_chg_reco_gubdb_kappa[ibin][ibin2]->SetMarkerStyle(2);

          h_chg_data_kappa[ibin][ibin2]->Draw("e0 same");
          h_chg_reco_up_kappa[ibin][ibin2]->Draw("e0 same");
          h_chg_reco_down_kappa[ibin][ibin2]->Draw("e0 same");
          h_chg_reco_gubdb_kappa[ibin][ibin2]->Draw("e0 same");
       
          TString tmp1 = cent_tag_pp[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.85, tmp1);
        
          TString tmp = kappa_tag[ibin2];
          tx = new TLatex(); tx->SetTextSize(.1);
          tx->DrawLatexNDC(0.15,0.65, tmp);
        
          if(do_fit==true){
            h_chg_data_kappa[ibin][ibin2]->Fit(f_chg_qg_kappa[ibin][ibin2],"Q M R N","sames",-2.,2.);
          
            if(ibin==0) c_chg_udg_MC_data_kappa[ibin2]->cd(6);
            else c_chg_udg_MC_data_kappa[ibin2]->cd(11-ibin);   
            h_chg_data_ratio_kappa[ibin][ibin2] = (TH1D*)h_chg_data_kappa[ibin][ibin2]->Clone("h_chg_data_kappa_ratio");
            h_chg_data_ratio_kappa[ibin][ibin2]->Divide(f_chg_qg_kappa[ibin][ibin2]);
            //h_chg_data_ratio_kappa[ibin][ibin2]->Rebin(20); h_chg_data_ratio_kappa[ibin][ibin2]->Scale(1./20.);
            h_chg_data_ratio_kappa[ibin][ibin2]->GetXaxis()->SetTitle("reco jet charge");
            h_chg_data_ratio_kappa[ibin][ibin2]->GetYaxis()->SetTitle("Data / Fit");
            h_cosm(h_chg_data_ratio_kappa[ibin][ibin2]);
            h_chg_data_ratio_kappa[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
            h_chg_data_ratio_kappa[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,2.);
            h_chg_data_ratio_kappa[ibin][ibin2]->SetLineColor(kBlack); h_chg_data_ratio_kappa[ibin][ibin2]->SetMarkerColor(kBlack);
            h_chg_data_ratio_kappa[ibin][ibin2]->SetMarkerStyle(1); h_chg_data_ratio_kappa[ibin][ibin2]->Draw("e0 same");
            tl1->Draw("same");
            sprintf(g_frac, "g+others %.3f",(1.-f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)));
            sprintf(g_error, "#pm %.3f",f_chg_qg_kappa[ibin][ibin2]->GetParError(0));
            sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_qg_kappa[ibin][ibin2]->GetChisquare()/f_chg_qg_kappa[ibin][ibin2]->GetNDF());
            sprintf(up_frac, "up+down %.3f",f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));
            tx1 = new TLatex(); tx1->SetTextSize(.07);
            tx1->DrawLatexNDC(0.15,0.75, g_frac);
            tx3 = new TLatex(); tx3->SetTextSize(.06);
            tx3->DrawLatexNDC(0.65,0.7, g_error);
            tx4 = new TLatex(); tx4->SetTextSize(.07);
            tx4->DrawLatexNDC(0.3,0.4, up_frac);
            tx2 = new TLatex(); tx2->SetTextSize(.07);
            tx2->DrawLatexNDC(0.3,0.2, chi2ndf);
            if(ibin==0){
                TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
                legendx ->SetLineColor(kWhite);
                legendx ->AddEntry(h_chg_data_kappa[0][0], "Data", "lepf");
                legendx ->AddEntry(h_chg_reco_gubdb_kappa[0][0], "Gluon + Others", "lepf");
                legendx ->AddEntry(h_chg_reco_up_kappa[0][0], "Up", "lepf");
                legendx ->AddEntry(h_chg_reco_down_kappa[0][0], "Down", "lepf");
                legendx ->Draw("same");
            }
          }
        }
      }

      //////vary trk eff fitting
      
      TH1D* h_chi2ndf_kappa[nCBins][ntrkbins];
      TH1D* h_least_chi2ndf_kappa[ntrkbins];      

      for(int ibin2=0; ibin2<ntrkbins; ibin2++){
          h_chi2ndf_kappa[0][ibin2] = new TH1D(Form("h_chi2ndf_kappa_%d_%d",0,ibin2),"",11,-7.,3.);
          h_chi2ndf_kappa[1][ibin2] = new TH1D(Form("h_chi2ndf_kappa_%d_%d",1,ibin2),"",11,-7.,3.);
          h_chi2ndf_kappa[2][ibin2] = new TH1D(Form("h_chi2ndf_kappa_%d_%d",2,ibin2),"",11,-7.,3.);
          h_chi2ndf_kappa[3][ibin2] = new TH1D(Form("h_chi2ndf_kappa_%d_%d",3,ibin2),"",11,-8.,2.);
          h_chi2ndf_kappa[4][ibin2] = new TH1D(Form("h_chi2ndf_kappa_%d_%d",4,ibin2),"",11,-5.,5.);

          h_least_chi2ndf_kappa[ibin2] = new TH1D(Form("h_least_chi2ndf_kappa_%d",ibin2),"",5,0.,5.);
      } 

      Double_t chi2NDF_kappa[nCBins][ntrkbins][trkeffbins+1];

      //////kappa two parameter fitting trkeff
      TCanvas *c_chg_qg_MC_data_kappa_trkeff[ntrkbins];

      for(int ibin2=0;ibin2<ntrkbins;ibin2++){
        c_chg_qg_MC_data_kappa_trkeff[ibin2] = new TCanvas((TString)("c_chg_qg_MC_data_kappa_trkeff"+trk[ibin2]),(TString)("c_chg_qg_MC_data_kappa_trkeff"+trk[ibin2]),1500,300);
        c_chg_qg_MC_data_kappa_trkeff[ibin2]->Divide(5,1,0);
        //gStyle->SetOptStat(0);
        //gStyle->SetOptFit(1);
        for(int ibin=0;ibin<nCBins;ibin++){

          if(ibin==0) c_chg_qg_MC_data_kappa_trkeff[ibin2]->cd(1);
          else c_chg_qg_MC_data_kappa_trkeff[ibin2]->cd(6-ibin);   

          for(int ibin3=0; ibin3<trkeffbins;ibin3++){
             h_chg_data_kappa_trkeff[ibin][ibin2][ibin3]->Fit(f_chg_qg_kappa[ibin][ibin2],"Q N M R","sames",-2.,2.);

             if(ibin3==7){
                 quark_fit_kappa_high[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)+(1-gluon_ubdbcsb[ibin])*(1-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)));
                 quark_fit_kappa_high[ibin][ibin2] = quark_fit_kappa_high[ibin][ibin2]*quark_corr_factor[ibin][ibin2];
                 quark_fit_kappa_high[ibin][ibin2] = (quark_fit_kappa_high[ibin][ibin2]*data_no0trk_int[ibin][0] + quark_data_0trk_int[ibin][0])/data_int[ibin];
                 gluon_fit_kappa_high[ibin][ibin2] = 1.-(quark_fit_kappa_high[ibin][ibin2]);
             }

             if(ibin3==3){
                 quark_fit_kappa_low[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)+(1-gluon_ubdbcsb[ibin])*(1-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)));
                 quark_fit_kappa_low[ibin][ibin2] = quark_fit_kappa_low[ibin][ibin2]*quark_corr_factor[ibin][ibin2];
                 quark_fit_kappa_low[ibin][ibin2] = (quark_fit_kappa_low[ibin][ibin2]*data_no0trk_int[ibin][0] + quark_data_0trk_int[ibin][0])/data_int[ibin];
                 gluon_fit_kappa_low[ibin][ibin2] = 1.-(quark_fit_kappa_low[ibin][ibin2]);
             }

             h_chg_reco_ratio_kappa_trkeff[ibin][ibin2][ibin3] = (TH1D*)h_chg_data_kappa_trkeff[ibin][ibin2][ibin3]->Clone("h_chg_data_kappa_ratio_trkeff");
             h_chg_reco_ratio_kappa_trkeff[ibin][ibin2][ibin3]->Divide(f_chg_qg_kappa[ibin][ibin2]);
             //h_chg_reco_ratio_kappa_trkeff[ibin][ibin2][ibin3]->Rebin(20); h_chg_reco_ratio_kappa_trkeff[ibin][ibin2][ibin3]->Scale(1./20.);
          
             chi2NDF_kappa[ibin][ibin2][ibin3] = f_chg_qg_kappa[ibin][ibin2]->GetChisquare()/f_chg_qg_kappa[ibin][ibin2]->GetNDF();
             h_chi2ndf_kappa[ibin][ibin2]->SetBinContent(ibin3+1,chi2NDF_kappa[ibin][ibin2][ibin3]);
          }

          h_chg_data_kappa[ibin][ibin2]->Fit(f_chg_qg_kappa[ibin][ibin2],"Q N M R","sames",-2.,2.);

          if(ibin==0) c_chg_qg_MC_data_kappa_trkeff[ibin2]->cd(1);
          else c_chg_qg_MC_data_kappa_trkeff[ibin2]->cd(6-ibin);   
          h_chg_data_ratio_kappa[ibin][ibin2]->Draw("HIST L same");

          for(int ibin3=0; ibin3<trkeffbins;ibin3++){
              if(ibin3<5){
                h_chg_reco_ratio_kappa_trkeff[ibin][ibin2][ibin3]->SetLineColor(kCyan-5+ibin3); h_chg_reco_ratio_kappa_trkeff[ibin][ibin2][ibin3]->SetMarkerColor(kCyan-5+ibin3);
              }
              else{
                h_chg_reco_ratio_kappa_trkeff[ibin][ibin2][ibin3]->SetLineColor(kGreen-5+ibin3); h_chg_reco_ratio_kappa_trkeff[ibin][ibin2][ibin3]->SetMarkerColor(kGreen-5+ibin3);
              }
              h_chg_reco_ratio_kappa_trkeff[ibin][ibin2][ibin3]->Rebin(10);
              h_chg_reco_ratio_kappa_trkeff[ibin][ibin2][ibin3]->Scale(1./10.);
              h_chg_reco_ratio_kappa_trkeff[ibin][ibin2][ibin3]->SetMarkerStyle(20); h_chg_reco_ratio_kappa_trkeff[ibin][ibin2][ibin3]->Draw("HIST L same");        
          }
          h_chg_data_ratio_kappa[ibin][ibin2]->Rebin(10);
          h_chg_data_ratio_kappa[ibin][ibin2]->Scale(1./10.);
          h_chg_data_ratio_kappa[ibin][ibin2]->Draw("HIST L same");
         
          if(ibin==0){
              TLegend *legendz = new TLegend(0.65,0.55,0.98,0.98);
              legendz ->SetLineColor(kWhite);
              legendz ->AddEntry(h_chg_reco_ratio[0][0], "Data (nominal)", "l");
              legendz ->Draw("same");
          }
          if(ibin==4){
              TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
              legendx ->SetLineColor(kWhite);
              for(int ibin3=0; ibin3<5;ibin3++){
                legendx ->AddEntry(h_chg_reco_ratio_trkeff[0][0][ibin3], trkeffvary[ibin3], "l");
              }
              legendx ->Draw("same");
          }
          if(ibin==3){
              TLegend *legendy = new TLegend(0.65,0.55,0.98,0.98);
              legendy ->SetLineColor(kWhite);
              for(int ibin3=5; ibin3<10;ibin3++){
                legendy ->AddEntry(h_chg_reco_ratio_trkeff[0][0][ibin3], trkeffvary[ibin3], "l");
              }
              legendy ->Draw("same");
          }

          TString tmp1 = cent_tag_pp[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.35,0.85, tmp1);

          tl1->Draw("same");
          
        }
      }

      ///chi2/ndf plot for varying trk eff
      TCanvas *c_chg_qg_MC_data_kappa_chi2ndf;
      c_chg_qg_MC_data_kappa_chi2ndf = new TCanvas("c_chg_qg_MC_data_kappa_chi2ndf","c_chg_qg_MC_data_kappa_chi2ndf",1800,300);
      c_chg_qg_MC_data_kappa_chi2ndf->Divide(5,1);

      for(int ibin=0;ibin<nCBins;ibin++){
        for(int ibin2=0;ibin2<4;ibin2++){
          if(ibin2==1) continue;

          h_least_chi2ndf_kappa[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp[ibin]);
          h_least_chi2ndf_kappa[ibin2]->SetBinContent(ibin+1,h_chi2ndf_kappa[ibin][ibin2]->GetBinContent(6));
          h_least_chi2ndf_kappa[ibin2]->SetBinError(ibin+1,0.);

          h_chi2ndf_kappa[ibin][ibin2]->GetXaxis()->SetTitle("trk. eff. vary (%)");
          h_chi2ndf_kappa[ibin][ibin2]->GetYaxis()->SetTitle("#chi^{2}/NDF");
          h_cosm(h_chi2ndf_kappa[ibin][ibin2]);
          h_chi2ndf_kappa[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
          h_chi2ndf_kappa[ibin][ibin2]->GetYaxis()->SetTitleSize(0.04);
          //h_chi2ndf[ibin][ibin2]->GetXaxis()->SetRangeUser(-5.,5.);
          //h_chi2ndf[ibin][ibin2]->GetYaxis()->SetRangeUser(0.8,1.2);
          h_chi2ndf_kappa[ibin][ibin2]->SetLineColor(ibin2+1); h_chi2ndf_kappa[ibin][ibin2]->SetMarkerColor(ibin2+1);
          h_chi2ndf_kappa[ibin][ibin2]->SetMarkerStyle(34); h_chi2ndf_kappa[ibin][ibin2]->SetMarkerSize(1.5);

          if(ibin==0) c_chg_qg_MC_data_kappa_chi2ndf->cd(1);
          else c_chg_qg_MC_data_kappa_chi2ndf->cd(6-ibin);   
          h_chi2ndf_kappa[ibin][ibin2]->Draw("p same");

          if(ibin==0){
              TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
              legendx ->SetLineColor(kWhite);
              legendx ->AddEntry(h_chi2ndf_kappa[0][ibin2], kappa_tag[ibin2], "lepf");
              legendx ->Draw("same");
          }
        }

        TString tmp1 = cent_tag_pp[ibin];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.35,0.85, tmp1);
      }

      TCanvas *c_chg_qg_least_chi2ndf_kappa;
      c_chg_qg_least_chi2ndf_kappa = new TCanvas("c_chg_qg_least_chi2ndf_kappa","c_chg_qg_least_chi2ndf_kappa",500,500);
      c_chg_qg_least_chi2ndf_kappa->cd(1);

      TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
      legendx ->SetLineColor(kWhite);

      for(int ibin2=0;ibin2<4;ibin2++){
        if(ibin2==1) continue;

        h_least_chi2ndf_kappa[ibin2]->GetXaxis()->SetTitle("");
        h_least_chi2ndf_kappa[ibin2]->GetYaxis()->SetTitle("#chi^{2}/NDF");
        h_cosm(h_least_chi2ndf_kappa[ibin2]);
        h_least_chi2ndf_kappa[ibin2]->SetLineColor(ibin2+1); h_least_chi2ndf_kappa[ibin2]->SetMarkerColor(ibin2+1);
        h_least_chi2ndf_kappa[ibin2]->SetMarkerStyle(34); h_least_chi2ndf_kappa[ibin2]->SetMarkerSize(1.5);
        h_least_chi2ndf_kappa[ibin2]->Draw("p same");
        legendx ->AddEntry(h_least_chi2ndf_kappa[ibin2], kappa_tag[ibin2], "lepf");
      }
      legendx ->Draw("same");
    }
/*

//////////////////spillover effect////////////////////////
  TCanvas *c_chg_spillover[ntrkbins];

  for(int ibin2=0;ibin2<1;ibin2++){
    c_chg_spillover[ibin2] = new TCanvas((TString)("c_chg_spillover"+trk[ibin2]),(TString)("c_chg_spillover"+trk[ibin2]),1200,450);

    TPad *pad_title = new TPad("pad_title", "",0.,0.92,1.,1.);
    TPad *pad_label = new TPad("pad_label", "",0.,0.,0.02,0.95);
    TPad *pad1 = new TPad("pad1", "",0.02,0.3,1.,0.98);
    pad1->Divide(4,1,0);  
    pad1->Draw();
    TPad *pad2 = new TPad("pad2", "",0.02,0.,1.,0.3);
    pad2->Divide(4,1,0);  
    pad2->Draw();    
    pad_title->Draw();
    pad_label->Draw();
    pad_title->cd();
    title_tex->DrawLatexNDC(0.8,0.3, trk_tag[ibin2]);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    pad_label->cd();
    ylabel_tex->DrawLatexNDC(0.7,0.5, y_label);
    ylabel_tex->DrawLatexNDC(0.7,0.1, diff);

    //c_chg_udg_MC_data_rebin[ibin2]->Divide(5,2,0);
    //gStyle->SetOptStat(0);
    //gStyle->SetOptFit(1);
    for(int ibin=1;ibin<nCBins;ibin++){
  
      //if(ibin==0) c_chg_udg_MC_data_rebin[ibin2]->cd(1);
      //else c_chg_udg_MC_data_rebin[ibin2]->cd(6-ibin);   
      pad1->cd(5-ibin);   

      //h_chg_ref[ibin][ibin2]->Rebin(10);
      h_chg_ref_sube0[ibin][ibin2]->Rebin(10);

      h_chg_spill[ibin][ibin2] = (TH1D*)h_chg_ref[ibin][ibin2]->Clone((TString)("h_chg_spill_"+cent[ibin]+"_"+trk[ibin2])); 
      h_chg_spill[ibin][ibin2] -> Add(h_chg_ref_sube0[ibin][ibin2],-1);

      h_chg_ref_sube0[ibin][ibin2]->GetXaxis()->SetTitle("jet charge");
      h_cosm(h_chg_ref_sube0[ibin][ibin2]);
      h_chg_ref_sube0[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
      h_chg_ref_sube0[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.82);
      h_chg_ref[ibin][ibin2]->SetLineColor(kBlue); h_chg_ref[ibin][ibin2]->SetMarkerColor(kBlue); h_chg_ref[ibin][ibin2]->SetMarkerStyle(20);
      h_chg_ref_sube0[ibin][ibin2]->SetLineColor(kRed); h_chg_ref_sube0[ibin][ibin2]->SetMarkerColor(kRed); h_chg_ref_sube0[ibin][ibin2]->SetMarkerStyle(20);    

      h_chg_ref_sube0[ibin][ibin2]->Draw("e0 same");    
      h_chg_ref[ibin][ibin2]->Draw("e0 same");      

      TString tmp1 = cent_tag_pp_MC[ibin];
      tx1 = new TLatex(); tx1->SetTextSize(.1);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);

      if(ibin==4){
      TLegend *legends = new TLegend(0.35,0.15,0.75,0.4);
      legends ->SetLineColor(kWhite);
      legends ->AddEntry(h_chg_ref[4][0], "sube all", "lep");
      legends ->AddEntry(h_chg_ref_sube0[4][0], "sube > 0 - #eta refl bkg", "lep");
      legends ->Draw("same");
      }
  
      pad2->cd(5-ibin);   

      h_chg_spill[ibin][ibin2]->GetXaxis()->SetTitle("");
      h_chg_spill[ibin][ibin2]->GetYaxis()->SetTitle("");
      h_cosm(h_chg_spill[ibin][ibin2]);
      h_chg_spill[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.5,1.5);
      h_chg_spill[ibin][ibin2]->GetYaxis()->SetRangeUser(-0.06,0.06);
      h_chg_spill[ibin][ibin2]->SetLineColor(kBlack); h_chg_spill[ibin][ibin2]->SetMarkerColor(kBlack);
      h_chg_spill[ibin][ibin2]->SetMarkerStyle(20); h_chg_spill[ibin][ibin2]->Draw("e0 same");
      tl3->Draw("same");
    }
  }
*/

    if(do_ptcut_fitting){
      TCanvas *c_chg_udg_MC_data_rebin[ntrkbins];

      for(int ibin2=3;ibin2<4;ibin2++){
        c_chg_udg_MC_data_rebin[ibin2] = new TCanvas((TString)("c_chg_udg_MC_data_rebin"+trk[ibin2]),(TString)("c_chg_udg_MC_data_rebin"+trk[ibin2]),1500,450);

        TPad *pad_title = new TPad("pad_title", "",0.,0.94,1.,1.);
        TPad *pad_label = new TPad("pad_label", "",0.,0.,0.02,0.95);
        TPad *pad1 = new TPad("pad1", "",0.02,0.3,1.,0.98);
        pad1->Divide(5,1,0);  
        pad1->Draw();
        TPad *pad2 = new TPad("pad2", "",0.02,0.,1.,0.3);
        pad2->Divide(5,1,0);  
        pad2->Draw();
        pad_title->Draw();
        pad_label->Draw();
        pad_title->cd();
        title_tex->DrawLatexNDC(0.4,0.1, data);
        title_tex->DrawLatexNDC(0.05,0.1, cms);  
        title_tex->DrawLatexNDC(0.8,0.3, trk_tag[ibin2]);  
        pad_label->cd();
        ylabel_tex->DrawLatexNDC(0.7,0.5, y_label);
        ylabel_tex->DrawLatexNDC(0.7,0.06, datafit);

        for(int ibin=0;ibin<nCBins;ibin++){
      
          if(ibin==0) pad1->cd(1);
          else pad1->cd(6-ibin);   
          h_chg_data_trkeff[ibin][ibin2][8]->Rebin(10);  
          h_chg_data_trkeff[ibin][ibin2][2]->Rebin(10);

          h_chg_data_cop[ibin][ibin2]->GetXaxis()->SetTitle("jet charge");
          h_cosm(h_chg_data_cop[ibin][ibin2]);
          h_chg_data_cop[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.8,1.8);
          //h_chg_data_cop[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.7);
          h_chg_data_cop[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.25);


          h_chg_data_cop[ibin][ibin2]->Draw("e0 same");

          //scaling according to data fractions
          if(ibin==0){
              h_chg_reco_up_cop[ibin][ibin2]->Scale(0.65*f_chg_qg_data[ibin][ibin2]->GetParameter(0));
              h_chg_reco_down_cop[ibin][ibin2]->Scale(0.35*f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
          }
          else{
              h_chg_reco_up_cop[ibin][ibin2]->Scale(0.47*f_chg_qg_data[ibin][ibin2]->GetParameter(0));
              h_chg_reco_down_cop[ibin][ibin2]->Scale(0.53*f_chg_qg_data[ibin][ibin2]->GetParameter(0));   
          }
          h_chg_reco_g_cop[ibin][ibin2]->Scale((1.-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)))*gluon_ubdbcsb[ibin]);
          h_chg_reco_others_cop[ibin][ibin2]->Scale((1.-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)))*(1-gluon_ubdbcsb[ibin]));

          h_chg_data_stk[ibin][ibin2] = new THStack((TString)("h_chg_data_stk"+cent[ibin]+"_"+trk[ibin2]),"");
          h_chg_data_stk[ibin][ibin2]->Add(h_chg_reco_others_cop[ibin][ibin2]);
          h_chg_data_stk[ibin][ibin2]->Add(h_chg_reco_down_cop[ibin][ibin2]);
          h_chg_data_stk[ibin][ibin2]->Add(h_chg_reco_up_cop[ibin][ibin2]);
          h_chg_data_stk[ibin][ibin2]->Add(h_chg_reco_g_cop[ibin][ibin2]);

          h_chg_data_stk[ibin][ibin2]->Draw("hist same");      

          h_chg_data_trkeff[ibin][ibin2][8]->Draw("hist same");  
          h_chg_data_trkeff[ibin][ibin2][2]->Draw("hist same");

          h_chg_data_cop[ibin][ibin2]->Draw("e0 same");  
      
          TString tmp1 = cent_tag_pp[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.85, tmp1);

          if(ibin==0){
            TLegend *legendx = new TLegend(0.15,0.7,0.5,0.8);
            legendx ->SetLineColor(kWhite);
            legendx ->AddEntry(h_chg_data_cop[0][ibin2], "Data", "lep");
            legendx ->Draw("same");

            TLegend *legendy = new TLegend(0.62,0.62,0.98,0.98);
            legendy ->SetLineColor(kWhite);
            legendy ->AddEntry((TObject*)0, "Fitting results", "");
            legendy ->AddEntry(h_chg_reco_g_cop[0][ibin2], "Gluon", "lepf");
            legendy ->AddEntry(h_chg_reco_up_cop[0][ibin2], "Up quark", "lepf");
            legendy ->AddEntry(h_chg_reco_down_cop[0][ibin2], "Down quark", "lepf");
            legendy ->AddEntry(h_chg_reco_others_cop[0][ibin2], "Other flavors", "lepf");
            legendy ->Draw("same");
          }

          if(do_fit==true){  
        
            if(ibin==0) pad2->cd(1);
            else pad2->cd(6-ibin);   

            h_chg_reco_ratio[ibin][ibin2] = (TH1D*)h_chg_reco_g_cop[ibin][ibin2]->Clone("h_chg_reco_ratio");
            h_chg_reco_ratio[ibin][ibin2]->Add(h_chg_reco_up_cop[ibin][ibin2]);
            h_chg_reco_ratio[ibin][ibin2]->Add(h_chg_reco_down_cop[ibin][ibin2]);
            h_chg_reco_ratio[ibin][ibin2]->Add(h_chg_reco_others_cop[ibin][ibin2]);
            h_chg_reco_ratio[ibin][ibin2]->Divide(h_chg_data_cop[ibin][ibin2]);
            h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetTitle("");
            h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetTitle("");
            //h_chg_reco_ratio[ibin][ibin2]->Rebin(20); h_chg_reco_ratio[ibin][ibin2]->Scale(1./20.);
            h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
            h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
            h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetLabelSize(0.13);
            h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->CenterTitle();
            h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->CenterTitle();
            h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetTitleOffset(0.7);
            h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
            h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
            h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
            h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetLabelSize(0.06);
            h_chg_reco_ratio[ibin][ibin2]->GetXaxis()->SetRangeUser(-1.8,1.8);
            h_chg_reco_ratio[ibin][ibin2]->GetYaxis()->SetRangeUser(0.5,1.5);
            h_chg_reco_ratio[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio[ibin][ibin2]->SetMarkerColor(kBlack);
            h_chg_reco_ratio[ibin][ibin2]->SetMarkerStyle(20); h_chg_reco_ratio[ibin][ibin2]->Draw("e0 same");
            tl1->Draw("same");
          }
        }
      } 
    }

    /////kappa final plot
    if(do_kappa_fitting){
      TCanvas *c_chg_udg_MC_data_rebin_kappa[ntrkbins];

      for(int ibin2=3;ibin2<4;ibin2++){
        c_chg_udg_MC_data_rebin_kappa[ibin2] = new TCanvas((TString)("c_chg_udg_MC_data_rebin_kappa"+trk[ibin2]),(TString)("c_chg_udg_MC_data_rebin_kappa"+trk[ibin2]),1500,450);

        TPad *pad_title = new TPad("pad_title", "",0.,0.94,1.,1.);
        TPad *pad_label = new TPad("pad_label", "",0.,0.,0.02,0.95);
        TPad *pad1 = new TPad("pad1", "",0.02,0.3,1.,0.98);
        pad1->Divide(5,1,0);  
        pad1->Draw();
        TPad *pad2 = new TPad("pad2", "",0.02,0.,1.,0.3);
        pad2->Divide(5,1,0);  
        pad2->Draw();
        pad_title->Draw();
        pad_label->Draw();
        pad_title->cd();
        title_tex->DrawLatexNDC(0.4,0.1, data);
        title_tex->DrawLatexNDC(0.05,0.1, cms);  
        title_tex->DrawLatexNDC(0.8,0.3, kappa_tag[ibin2]);  
        //title_tex->DrawLatexNDC(0.8,0.3, kappa_tag[ibin2]);  
        pad_label->cd();
        ylabel_tex->DrawLatexNDC(0.7,0.5, y_label);
        ylabel_tex->DrawLatexNDC(0.7,0.06, datafit);

        for(int ibin=0;ibin<nCBins;ibin++){
      
          if(ibin==0) pad1->cd(1);
          else pad1->cd(6-ibin);   

          h_chg_data_cop_kappa[ibin][ibin2]->Rebin(2);
          h_chg_reco_up_cop_kappa[ibin][ibin2]->Rebin(2);
          h_chg_reco_down_cop_kappa[ibin][ibin2]->Rebin(2);
          h_chg_reco_others_cop_kappa[ibin][ibin2]->Rebin(2);
          h_chg_reco_g_cop_kappa[ibin][ibin2]->Rebin(2);

          h_chg_data_cop_kappa[ibin][ibin2]->GetXaxis()->SetTitle("jet charge");
          h_cosm(h_chg_data_cop_kappa[ibin][ibin2]);
          h_chg_data_cop_kappa[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
          //h_chg_data_cop_kappa[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.7);
          h_chg_data_cop_kappa[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.4);

          h_chg_data_cop_kappa[ibin][ibin2]->Draw("e0 same");

          if(ibin==0){
              h_chg_reco_up_cop_kappa[ibin][ibin2]->Scale(0.65*f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));
              h_chg_reco_down_cop_kappa[ibin][ibin2]->Scale(0.35*f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
          }
          else{
              h_chg_reco_up_cop_kappa[ibin][ibin2]->Scale(0.47*f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));
              h_chg_reco_down_cop_kappa[ibin][ibin2]->Scale(0.53*f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));   
          }
          h_chg_reco_g_cop_kappa[ibin][ibin2]->Scale((1.-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)))*gluon_ubdbcsb[ibin]);
          h_chg_reco_others_cop_kappa[ibin][ibin2]->Scale((1.-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)))*(1-gluon_ubdbcsb[ibin]));

          h_chg_data_stk_kappa[ibin][ibin2] = new THStack((TString)("h_chg_data_stk_kappa_"+cent[ibin]+"_"+trk[ibin2]),"");
          h_chg_data_stk_kappa[ibin][ibin2]->Add(h_chg_reco_others_cop_kappa[ibin][ibin2]);
          h_chg_data_stk_kappa[ibin][ibin2]->Add(h_chg_reco_down_cop_kappa[ibin][ibin2]);
          h_chg_data_stk_kappa[ibin][ibin2]->Add(h_chg_reco_up_cop_kappa[ibin][ibin2]);
          h_chg_data_stk_kappa[ibin][ibin2]->Add(h_chg_reco_g_cop_kappa[ibin][ibin2]);

          h_chg_data_stk_kappa[ibin][ibin2]->Draw("hist same");      
          h_chg_data_cop_kappa[ibin][ibin2]->Draw("e0 same");  

          TString tmp1 = cent_tag_pp[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.85, tmp1);

          if(ibin==0){
            TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
            legendx ->SetLineColor(kWhite);
            legendx ->AddEntry(h_chg_data_cop_kappa[0][ibin2], "Data", "lep");
            legendx ->AddEntry(h_chg_reco_g_cop_kappa[0][ibin2], "Gluon", "lepf");
            legendx ->AddEntry(h_chg_reco_up_cop_kappa[0][ibin2], "Up quark", "lepf");
            legendx ->AddEntry(h_chg_reco_down_cop_kappa[0][ibin2], "Down quark", "lepf");
            legendx ->AddEntry(h_chg_reco_others_cop_kappa[0][ibin2], "Others", "lepf");
            legendx ->Draw("same");
          }

          if(do_fit==true){  
        
            if(ibin==0) pad2->cd(1);
            else pad2->cd(6-ibin);   

            h_chg_data_ratio_kappa[ibin][ibin2]->GetXaxis()->SetTitle("");
            h_chg_data_ratio_kappa[ibin][ibin2]->GetYaxis()->SetTitle("");
            h_chg_data_ratio_kappa[ibin][ibin2]->Rebin(2); h_chg_data_ratio_kappa[ibin][ibin2]->Scale(1./2.);
            h_chg_data_ratio_kappa[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
            h_chg_data_ratio_kappa[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
            h_chg_data_ratio_kappa[ibin][ibin2]->GetYaxis()->SetLabelSize(0.13);
            h_chg_data_ratio_kappa[ibin][ibin2]->GetXaxis()->CenterTitle();
            h_chg_data_ratio_kappa[ibin][ibin2]->GetYaxis()->CenterTitle();
            h_chg_data_ratio_kappa[ibin][ibin2]->GetYaxis()->SetTitleOffset(0.7);
            h_chg_data_ratio_kappa[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
            h_chg_data_ratio_kappa[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
            h_chg_data_ratio_kappa[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
            h_chg_data_ratio_kappa[ibin][ibin2]->GetXaxis()->SetLabelSize(0.06);
            h_chg_data_ratio_kappa[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
            h_chg_data_ratio_kappa[ibin][ibin2]->GetYaxis()->SetRangeUser(0.5,1.5);
            h_chg_data_ratio_kappa[ibin][ibin2]->SetLineColor(kBlack); h_chg_data_ratio_kappa[ibin][ibin2]->SetMarkerColor(kBlack);
            h_chg_data_ratio_kappa[ibin][ibin2]->SetMarkerStyle(20); h_chg_data_ratio_kappa[ibin][ibin2]->Draw("e0 same");
            tl1->Draw("same");
          }
        }
      }  
    }
  }
/*
  ////spillover in data
  TCanvas *c_chg_udg_MC_data_rebin_spillover[ntrkbins];

  for(int ibin2=0;ibin2<1;ibin2++){
    c_chg_udg_MC_data_rebin_spillover[ibin2] = new TCanvas((TString)("c_chg_udg_MC_data_rebin_spillover"+trk[ibin2]),(TString)("c_chg_udg_MC_data_rebin_spillover"+trk[ibin2]),1200,450);

    TPad *pad_title = new TPad("pad_title", "",0.,0.94,1.,1.);
    TPad *pad_label = new TPad("pad_label", "",0.,0.,0.02,0.95);
    TPad *pad1 = new TPad("pad1", "",0.02,0.3,1.,0.98);
    pad1->Divide(4,1,0);  
    pad1->Draw();
    TPad *pad2 = new TPad("pad2", "",0.02,0.,1.,0.3);
    pad2->Divide(4,1,0);  
    pad2->Draw();
    pad_title->Draw();
    pad_label->Draw();
    pad_title->cd();
    title_tex->DrawLatexNDC(0.4,0.1, data);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    title_tex->DrawLatexNDC(0.8,0.3, trk_tag[ibin2]);  
    //title_tex->DrawLatexNDC(0.8,0.3, kappa_tag[ibin2]);  
    pad_label->cd();
    ylabel_tex->DrawLatexNDC(0.7,0.5, y_label);
    ylabel_tex->DrawLatexNDC(0.7,0.06, "ratio");

    //c_chg_udg_MC_data_rebin[ibin2]->Divide(5,2,0);
    //gStyle->SetOptStat(0);
    //gStyle->SetOptFit(1);
    for(int ibin=1;ibin<nCBins;ibin++){
  
      //if(ibin==0) c_chg_udg_MC_data_rebin[ibin2]->cd(1);
      //else c_chg_udg_MC_data_rebin[ibin2]->cd(6-ibin);   
      pad1->cd(5-ibin);  

      h_chg_data_spillover[ibin][ibin2]->SetLineColor(kBlack); 
      h_chg_data_spillover[ibin][ibin2]->SetMarkerColor(kBlack); 
      h_chg_data_spillover[ibin][ibin2]->SetMarkerStyle(4);
 
      h_chg_data_cop[ibin][ibin2]->Draw("e0 same");
      h_chg_data_spillover[ibin][ibin2]->Draw("e0 same");
      h_chg_data_spilloverdn[ibin][ibin2]->Draw("e0 same");      

      TString tmp1 = cent_tag_pp[ibin];
      tx1 = new TLatex(); tx1->SetTextSize(.1);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);

      if(ibin==4){
      TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
      legendx ->SetLineColor(kWhite);
      legendx ->AddEntry(h_chg_data_cop[4][ibin2], "Data (nominal)", "lep");
      legendx ->AddEntry(h_chg_data_spillover[4][ibin2], "Data (bkg. fluctuation up)", "lep");
      legendx ->AddEntry(h_chg_data_spilloverdn[4][ibin2], "Data (bkg. fluctuation down)", "lep");
      legendx ->Draw("same");
      }

      pad2->cd(5-ibin);   

      h_chg_reco_ratio_spillover[ibin][ibin2] = (TH1D*)h_chg_data_cop[ibin][ibin2]->Clone("h_chg_reco_ratio_spillover");
      h_chg_reco_ratio_spillover[ibin][ibin2]->Divide(h_chg_data_spillover[ibin][ibin2]);
      h_chg_reco_ratio_spilloverdn[ibin][ibin2] = (TH1D*)h_chg_data_cop[ibin][ibin2]->Clone("h_chg_reco_ratio_spilloverdn");
      h_chg_reco_ratio_spilloverdn[ibin][ibin2]->Divide(h_chg_data_spilloverdn[ibin][ibin2]);
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetXaxis()->SetTitle("");
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetYaxis()->SetTitle("");
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetYaxis()->SetLabelSize(0.13);
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetXaxis()->CenterTitle();
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetYaxis()->CenterTitle();
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetYaxis()->SetTitleOffset(0.7);
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetXaxis()->SetLabelSize(0.06);
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
      h_chg_reco_ratio_spillover[ibin][ibin2]->GetYaxis()->SetRangeUser(0.8,1.2);
      h_chg_reco_ratio_spillover[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio_spillover[ibin][ibin2]->SetMarkerColor(kBlack);
      h_chg_reco_ratio_spillover[ibin][ibin2]->SetMarkerStyle(20); h_chg_reco_ratio_spillover[ibin][ibin2]->Draw("e0 same");
      h_chg_reco_ratio_spilloverdn[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio_spilloverdn[ibin][ibin2]->SetMarkerColor(kBlack);
      h_chg_reco_ratio_spilloverdn[ibin][ibin2]->SetMarkerStyle(29); h_chg_reco_ratio_spilloverdn[ibin][ibin2]->Draw("e0 same");
      tl1->Draw("same");

      if(ibin==4){
      TLegend *legendy = new TLegend(0.65,0.55,0.98,0.98);
      legendy ->SetLineColor(kWhite);
      legendy ->AddEntry(h_chg_reco_ratio_spillover[4][ibin2], "nominal / bkg. fluctuation up", "lep");
      legendy ->AddEntry(h_chg_reco_ratio_spilloverdn[4][ibin2], "nominal / bkg. fluctuation down", "lep");
      legendy ->Draw("same");
      }
    }
  }  
*/
/*  
////varying positve trk eff in data
  TCanvas *c_chg_udg_MC_data_rebin_trkeff[ntrkbins];

  for(int ibin2=0;ibin2<1;ibin2++){
    c_chg_udg_MC_data_rebin_trkeff[ibin2] = new TCanvas((TString)("c_chg_udg_MC_data_rebin_trkeff"+trk[ibin2]),(TString)("c_chg_udg_MC_data_rebin_trkeff"+trk[ibin2]),1500,450);

    TPad *pad_title = new TPad("pad_title", "",0.,0.94,1.,1.);
    TPad *pad_label = new TPad("pad_label", "",0.,0.,0.02,0.95);
    TPad *pad1 = new TPad("pad1", "",0.02,0.3,1.,0.98);
    pad1->Divide(5,1,0);  
    pad1->Draw();
    TPad *pad2 = new TPad("pad2", "",0.02,0.,1.,0.3);
    pad2->Divide(5,1,0);  
    pad2->Draw();
    pad_title->Draw();
    pad_label->Draw();
    pad_title->cd();
    title_tex->DrawLatexNDC(0.4,0.1, data);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    title_tex->DrawLatexNDC(0.8,0.3, trk_tag[ibin2]);  
    //title_tex->DrawLatexNDC(0.8,0.3, kappa_tag[ibin2]);  
    pad_label->cd();
    ylabel_tex->DrawLatexNDC(0.7,0.5, y_label);
    ylabel_tex->DrawLatexNDC(0.7,0.06, "ratio");

    //c_chg_udg_MC_data_rebin[ibin2]->Divide(5,2,0);
    //gStyle->SetOptStat(0);
    //gStyle->SetOptFit(1);
    for(int ibin=0;ibin<nCBins;ibin++){
  
      if(ibin==0) pad1->cd(1);
      else pad1->cd(6-ibin);  

      h_chg_data_trkpos[ibin][ibin2]->Rebin(10);
      h_chg_data_trkpos[ibin][ibin2]->SetLineColor(kBlack); 
      h_chg_data_trkpos[ibin][ibin2]->SetMarkerColor(kBlack); 
      h_chg_data_trkpos[ibin][ibin2]->SetMarkerStyle(4);
 
      h_chg_data_cop[ibin][ibin2]->Draw("e0 same");
      h_chg_data_trkpos[ibin][ibin2]->Draw("e0 same");

      TString tmp1 = cent_tag_pp[ibin];
      tx1 = new TLatex(); tx1->SetTextSize(.1);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);

      if(ibin==4){
      TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
      legendx ->SetLineColor(kWhite);
      legendx ->AddEntry(h_chg_data_cop[4][ibin2], "Data (nominal)", "lep");
      legendx ->AddEntry(h_chg_data_spillover[4][ibin2], "Data (pos/neg trk. eff.)", "lep");
      legendx ->Draw("same");
      }
    }
  }
*/
/*
///////JER effects
////trk eff in data
TCanvas *c_chg_udg_MC_data_rebin_jer[ntrkbins];

  for(int ibin2=0;ibin2<1;ibin2++){
    c_chg_udg_MC_data_rebin_jer[ibin2] = new TCanvas((TString)("c_chg_udg_MC_data_rebin_jer"+trk[ibin2]),(TString)("c_chg_udg_MC_data_rebin_jer"+trk[ibin2]),300,450);

    TPad *pad_title = new TPad("pad_title", "",0.,0.94,1.,1.);
    TPad *pad_label = new TPad("pad_label", "",0.,0.,0.02,0.95);
    TPad *pad1 = new TPad("pad1", "",0.02,0.3,1.,0.98);
    pad1->Divide(1,1,0);  
    pad1->Draw();
    TPad *pad2 = new TPad("pad2", "",0.02,0.,1.,0.3);
    pad2->Divide(1,1,0);  
    pad2->Draw();
    pad_title->Draw();
    pad_label->Draw();
    pad_title->cd();
    //title_tex->DrawLatexNDC(0.4,0.1, data);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    title_tex->DrawLatexNDC(0.8,0.3, trk_tag[ibin2]);  
    //title_tex->DrawLatexNDC(0.8,0.3, kappa_tag[ibin2]);  
    pad_label->cd();
    ylabel_tex->DrawLatexNDC(0.7,0.5, y_label);
    ylabel_tex->DrawLatexNDC(0.7,0.06, "ratio");

    pad1->cd(1);
    h_chg_data_jer[1][ibin2]->Rebin(10);
    h_chg_data_jer[1][ibin2]->SetLineColor(kBlack); 
    h_chg_data_jer[1][ibin2]->SetMarkerColor(kBlack); 
    h_chg_data_jer[1][ibin2]->SetMarkerStyle(4);

    h_chg_data_cop[1][ibin2]->Draw("e0 same");
    h_chg_data_jer[1][ibin2]->Draw("e0 same");

    TString tmp1 = cent_tag_pp[1];
    tx1 = new TLatex(); tx1->SetTextSize(.1);
    tx1->DrawLatexNDC(0.15,0.85, tmp1);

    TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
    legendx ->SetLineColor(kWhite);
    legendx ->AddEntry(h_chg_data_cop[1][ibin2], "Data (nominal)", "lep");
    legendx ->AddEntry(h_chg_data_spillover[1][ibin2], "Data (modified)", "lep");
    legendx ->Draw("same");

    pad2->cd(1);   
    h_chg_reco_ratio_jer[1][ibin2] = (TH1D*)h_chg_data_cop[1][ibin2]->Clone(Form("h_chg_data_cop_%d",ibin2));
    h_chg_reco_ratio_jer[1][ibin2]->Divide(h_chg_data_jer[1][ibin2]);
    h_cosm(h_chg_reco_ratio_jer[1][ibin2]);
    h_chg_reco_ratio_jer[1][ibin2]->Draw("e0 same");
    TLegend *legendy = new TLegend(0.65,0.55,0.98,0.98);
    legendy ->SetLineColor(kWhite);
    legendy ->AddEntry(h_chg_reco_ratio_spillover[1][ibin2], "nominal / modified", "lep");
    legendy ->Draw("same");
    drawline(-1.8,1.,1.8,1.);
    
  }
*/

  if(do_data){
    ////calculating final q and g values
    for(int ibin=0; ibin<nCBins; ibin++){
      for(int ibin2=0;ibin2<ntrkbins-1;ibin2++){
        if(ibin2==0){
          if(do_ptcut_fitting){
            if(ibin==0){
              up_fit[ibin][ibin2] = 0.65*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); up_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
              down_fit[ibin][ibin2] = 0.35*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); down_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
            }
            else{
              up_fit[ibin][ibin2] = 0.47*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); up_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
              down_fit[ibin][ibin2] = 0.53*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); down_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
            }
            quark_fit[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParameter(0)+(1-gluon_ubdbcsb[ibin])*(1-(f_chg_qg_data[ibin][ibin2]->GetParameter(0))); quark_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
            quark_fit[ibin][ibin2] = quark_fit[ibin][ibin2]*quark_corr_factor[ibin][ibin2];
            quark_fit[ibin][ibin2] = (quark_fit[ibin][ibin2]*data_no0trk_int[ibin][ibin2] + quark_data_0trk_int[ibin][ibin2])/data_int[ibin];
            gluon_fit[ibin][ibin2] = 1.-(quark_fit[ibin][ibin2]); gluon_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
            /*
            gluon_fit[ibin][ibin2] = gluon_ubdbcsb[ibin]*(1-(f_chg_qg_data[ibin][ibin2]->GetParameter(0))); gluon_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
            gluon_fit[ibin][ibin2] = gluon_fit[ibin][ibin2]*gluon_corr_factor[ibin][ibin2];
            quark_fit[ibin][ibin2] = 1.-(gluon_fit[ibin][ibin2]); gluon_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
            */
            gluon_fit_high_diff[ibin][ibin2] =  fabs(gluon_fit[ibin][ibin2] - gluon_fit_high[ibin][ibin2]);
            gluon_fit_low_diff[ibin][ibin2] = fabs(gluon_fit_low[ibin][ibin2] - gluon_fit[ibin][ibin2]);

            //gluon_fit_high_diff[ibin][ibin2] =  sqrt(pow(fabs(gluon_fit[ibin][ibin2] - gluon_fit_high[ibin][ibin2]),2)+pow(diff_0trk_data_MC[ibin][ibin2],2));
            //gluon_fit_low_diff[ibin][ibin2] = sqrt(pow(fabs(gluon_fit_low[ibin][ibin2] - gluon_fit[ibin][ibin2]),2)+pow(diff_0trk_data_MC[ibin][ibin2],2));
          }

          if(do_kappa_fitting){
            if(ibin==0){
              up_fit_kappa[ibin][ibin2] = 0.65*f_chg_qg_kappa[ibin][ibin2]->GetParameter(0); up_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParError(0);
              down_fit_kappa[ibin][ibin2] = 0.35*f_chg_qg_kappa[ibin][ibin2]->GetParameter(0); down_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParError(0);
            }
            else{
              up_fit_kappa[ibin][ibin2] = 0.47*f_chg_qg_kappa[ibin][ibin2]->GetParameter(0); up_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParError(0);
              down_fit_kappa[ibin][ibin2] = 0.53*f_chg_qg_kappa[ibin][ibin2]->GetParameter(0); down_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParError(0);
            }
            quark_fit_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)+(1-gluon_ubdbcsb[ibin])*(1-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0))); quark_fiterr[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParError(0);
            quark_fit_kappa[ibin][ibin2] = quark_fit_kappa[ibin][ibin2]*quark_corr_factor[ibin][ibin2];
            quark_fit_kappa[ibin][ibin2] = (quark_fit_kappa[ibin][ibin2]*data_no0trk_int[ibin][0] + quark_data_0trk_int[ibin][0])/data_int[ibin];
            gluon_fit_kappa[ibin][ibin2] = 1.-(quark_fit_kappa[ibin][ibin2]); gluon_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParError(0);

            //quark_fit[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)+(1-gluon_ubdbcsb[ibin])*(1-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0))); quark_fiterr[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParError(0);
            //quark_fit[ibin][ibin2] = quark_fit[ibin][ibin2]*quark_corr_factor[ibin][ibin2];
            //gluon_fit[ibin][ibin2] = 1.-(quark_fit[ibin][ibin2]); gluon_fiterr[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParError(0);

            gluon_fit_kappa_high_diff[ibin][ibin2] =  fabs(gluon_fit_kappa[ibin][ibin2] - gluon_fit_kappa_high[ibin][ibin2]);
            gluon_fit_kappa_low_diff[ibin][ibin2] = fabs(gluon_fit_kappa_low[ibin][ibin2] - gluon_fit_kappa[ibin][ibin2]);
          }
        }
        else{
          if(do_ptcut_fitting){
            if(ibin==0){
              up_fit[ibin][ibin2] = 0.65*(f_chg_qg_data[ibin][ibin2+1]->GetParameter(0)); up_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParError(0);
              down_fit[ibin][ibin2] = 0.35*(f_chg_qg_data[ibin][ibin2+1]->GetParameter(0)); down_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParError(0);          
            }
            else{    
              up_fit[ibin][ibin2] = 0.47*(f_chg_qg_data[ibin][ibin2+1]->GetParameter(0)); up_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParError(0);
              down_fit[ibin][ibin2] = 0.53*(f_chg_qg_data[ibin][ibin2+1]->GetParameter(0)); down_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParError(0);
            }
            quark_fit[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParameter(0)+(1-gluon_ubdbcsb[ibin])*(1-(f_chg_qg_data[ibin][ibin2+1]->GetParameter(0))); quark_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParError(0);
            quark_fit[ibin][ibin2] = quark_fit[ibin][ibin2]*quark_corr_factor[ibin][ibin2+1];
            quark_fit[ibin][ibin2] = (quark_fit[ibin][ibin2]*data_no0trk_int[ibin][ibin2+1] + quark_data_0trk_int[ibin][ibin2+1])/data_int[ibin];
            gluon_fit[ibin][ibin2] = 1.-(quark_fit[ibin][ibin2]); gluon_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParError(0);
            /* 
            gluon_fit[ibin][ibin2] = gluon_ubdbcsb[ibin]*(1-(f_chg_qg_data[ibin][ibin2+1]->GetParameter(0))); gluon_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParError(0);
            gluon_fit[ibin][ibin2] = gluon_fit[ibin][ibin2]*gluon_corr_factor[ibin][ibin2+1];
            quark_fit[ibin][ibin2] = 1.-(gluon_fit[ibin][ibin2]); gluon_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2+1]->GetParError(0);
            */ 
            //gluon_fit_high_diff[ibin][ibin2] =  sqrt(pow(fabs(gluon_fit[ibin][ibin2] - gluon_fit_high[ibin][ibin2+1]),2)+pow(diff_0trk_data_MC[ibin][ibin2+1],2));
            //gluon_fit_low_diff[ibin][ibin2] = sqrt(pow(fabs(gluon_fit_low[ibin][ibin2+1] - gluon_fit[ibin][ibin2]),2)+pow(diff_0trk_data_MC[ibin][ibin2+1],2));

            gluon_fit_high_diff[ibin][ibin2] =  fabs(gluon_fit[ibin][ibin2] - gluon_fit_high[ibin][ibin2+1]);
            gluon_fit_low_diff[ibin][ibin2] = fabs(gluon_fit_low[ibin][ibin2+1] - gluon_fit[ibin][ibin2]);
          }

          if(do_kappa_fitting){
            if(ibin==0){
              up_fit_kappa[ibin][ibin2] = 0.65*f_chg_qg_kappa[ibin][ibin2+1]->GetParameter(0); up_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParError(0);
              down_fit_kappa[ibin][ibin2] = 0.35*f_chg_qg_kappa[ibin][ibin2+1]->GetParameter(0); down_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParError(0);
            }
            else{
              up_fit_kappa[ibin][ibin2] = 0.47*f_chg_qg_kappa[ibin][ibin2+1]->GetParameter(0); up_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParError(0);
              down_fit_kappa[ibin][ibin2] = 0.53*f_chg_qg_kappa[ibin][ibin2+1]->GetParameter(0); down_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParError(0);
            }
            quark_fit_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParameter(0)+(1-gluon_ubdbcsb[ibin])*(1-(f_chg_qg_kappa[ibin][ibin2+1]->GetParameter(0))); quark_fiterr[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParError(0);
            quark_fit_kappa[ibin][ibin2] = quark_fit_kappa[ibin][ibin2]*quark_corr_factor[ibin][0];
            quark_fit_kappa[ibin][ibin2] = (quark_fit_kappa[ibin][ibin2]*data_no0trk_int[ibin][0] + quark_data_0trk_int[ibin][0])/data_int[ibin];
            gluon_fit_kappa[ibin][ibin2] = 1.-(quark_fit_kappa[ibin][ibin2]); gluon_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParError(0);

            //quark_fit[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParameter(0)+(1-gluon_ubdbcsb[ibin])*(1-(f_chg_qg_kappa[ibin][ibin2+1]->GetParameter(0))); quark_fiterr[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParError(0);
            //quark_fit[ibin][ibin2] = quark_fit[ibin][ibin2]*quark_corr_factor[ibin][ibin2+1];
            //gluon_fit[ibin][ibin2] = 1.-(quark_fit[ibin][ibin2]); gluon_fiterr[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2+1]->GetParError(0);

            gluon_fit_kappa_high_diff[ibin][ibin2] =  fabs(gluon_fit_kappa[ibin][ibin2] - gluon_fit_kappa_high[ibin][ibin2+1]);
            gluon_fit_kappa_low_diff[ibin][ibin2] = fabs(gluon_fit_kappa_low[ibin][ibin2+1] - gluon_fit_kappa[ibin][ibin2]);

          }
        }

        gluon_fitsys[ibin][ibin2] = gluon_fit[ibin][ibin2];
        quark_fitsys[ibin][ibin2] = quark_fit[ibin][ibin2];

        if(do_kappa_fitting) cout<<"gluon_withsys_kappa["<<ibin<<"]["<<ibin2<<"] = "<<gluon_fit_kappa[ibin][ibin2]<<";"<<endl;
        if(do_ptcut_fitting) cout<<"gluon_withsys["<<ibin<<"]["<<ibin2<<"] = "<<gluon_fit[ibin][ibin2]<<";"<<endl;

        //cout<<"trkeff_err["<<ibin<<"]["<<ibin2<<"] = "<<(gluon_fit_high_diff[ibin][ibin2]+gluon_fit_low_diff[ibin][ibin2])/2.<<";"<<endl;
        
        gluon_recofit_ratio[ibin][ibin2] = gluon_fitsys[ibin][ibin2]/gluon_reco_fitsys[ibin][ibin2]; 
        
        quark_ratio[ibin][ibin2] = quark_fitsys[ibin][ibin2]/quark_fitsys[0][ibin2];
        gluon_ratio[ibin][ibin2] = gluon_fitsys[ibin][ibin2]/gluon_fitsys[0][ibin2];
        
      }

      //data qg graphs
      if(do_ptcut_fitting){
          up_data[ibin] = new TGraphErrors(ntrkbins-1,x_mean,up_fit[ibin],x_mean_err,up_fiterr[ibin]);
          up_data[ibin]->SetName((TString)("up_data_cent"+cent[ibin]));
          down_data[ibin] = new TGraphErrors(ntrkbins-1,x_mean,down_fit[ibin],x_mean_err,down_fiterr[ibin]);
          down_data[ibin]->SetName((TString)("down_data_cent"+cent[ibin]));
          gluon_data[ibin] = new TGraphErrors(ntrkbins-1,x_mean,gluon_fit[ibin],x_mean_err,gluon_fiterr[ibin]);
          gluon_data[ibin]->SetName((TString)("gluon_data_cent"+cent[ibin]));
          gluon_data_asymm_sys[ibin] = new TGraphAsymmErrors(ntrkbins-1,x_mean,gluon_fit[ibin],x_mean_err,x_mean_err,gluon_fit_low_diff[ibin],gluon_fit_high_diff[ibin]);
          gluon_data_asymm_sys[ibin]->SetName((TString)("gluon_data_sys_cent"+cent[ibin]));
          quark_data[ibin] = new TGraphErrors(ntrkbins-1,x_mean,quark_fit[ibin],x_mean_err,quark_fiterr[ibin]);
          quark_data[ibin]->SetName((TString)("quark_data_cent"+cent[ibin]));
      }

      if(do_kappa_fitting){
          up_data_kappa[ibin] = new TGraphErrors(ntrkbins-1,x_mean,up_fit_kappa[ibin],x_mean_err,up_fiterr_kappa[ibin]);
          up_data_kappa[ibin]->SetName((TString)("up_data_kappa_cent"+cent[ibin]));
          down_data_kappa[ibin] = new TGraphErrors(ntrkbins-1,x_mean,down_fit_kappa[ibin],x_mean_err,down_fiterr_kappa[ibin]);
          down_data_kappa[ibin]->SetName((TString)("down_data_kappa_cent"+cent[ibin]));
          gluon_data_kappa[ibin] = new TGraphErrors(ntrkbins-1,kappa,gluon_fit_kappa[ibin],kappa_err,gluon_fiterr_kappa[ibin]);
          gluon_data_kappa[ibin]->SetName((TString)("gluon_data_kappa_cent"+cent[ibin]));
          quark_data_kappa[ibin] = new TGraphErrors(ntrkbins-1,kappa,quark_fit_kappa[ibin],kappa_err,quark_fiterr_kappa[ibin]);
          quark_data_kappa[ibin]->SetName((TString)("quark_data_kappa_cent"+cent[ibin]));
      }
    }

    if(do_ptcut_fitting){
      TCanvas *c_data_graph = new TCanvas("c_data_graph","c_data_graph",1800,350);
      c_data_graph->Divide(6,1,0);
      
      for(int ibin=0; ibin<nCBins; ibin++){
        if(ibin==0) c_data_graph->cd(1);
        else c_data_graph->cd(7-ibin);
        h_dummy->Draw("same");
        up_data[ibin]->SetTitle(" ; pT cut; fraction; ");
        up_data[ibin]->GetXaxis()->SetTitle(" ; pT cut; fraction; ");
        up_data[ibin]->SetMarkerColor(kBlue); up_data[ibin]->SetMarkerStyle(26); up_data[ibin]->SetMarkerSize(1.3);
        up_data[ibin]->Draw("PL");
        down_data[ibin]->SetMarkerColor(kBlue); down_data[ibin]->SetMarkerStyle(32); down_data[ibin]->SetMarkerSize(1.3);
        down_data[ibin]->Draw("PL");
        gluon_data[ibin]->SetMarkerColor(kRed); gluon_data[ibin]->SetMarkerStyle(20); gluon_data[ibin]->SetMarkerSize(1.2);
        gluon_data[ibin]->Draw("PL");
        quark_data[ibin]->SetMarkerColor(kBlue); quark_data[ibin]->SetMarkerStyle(4); quark_data[ibin]->SetMarkerSize(1.2);
        quark_data[ibin]->Draw("PL");

        if(ibin==0) {tl_u[ibin]->Draw("same"); tl_d[ibin]->Draw("same"); tl_g[ibin]->Draw("same"); tl_q[ibin]->Draw("same");}
        else {tl_u_data->Draw("same"); tl_d_data->Draw("same"); tl_g[ibin]->Draw("same"); tl_q[ibin]->Draw("same");}
            
        TString tmp1 = cent_tag_pp[ibin];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.15,0.85, tmp1);
        
      }
    }

    TCanvas *c_data_qg_graph = new TCanvas("c_data_qg_graph","c_data_qg_graph",1500,350);
    drawgraphPad();

    pad_titleg->cd();
    title_tex->DrawLatexNDC(0.5,0.1, data);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    pad_labelg->cd();
    chg_tex->DrawLatexNDC(0.7,0.4, fraction);

    for(int ibin=0; ibin<nCBins; ibin++){
      //if(ibin==0) c_data_qg_graph->cd(1);
      //else c_data_qg_graph->cd(6-ibin);
      if(ibin==0) pada->cd(1);
      else pada->cd(6-ibin);
      h_dummy->Draw("same");
      
      if(do_ptcut_fitting){
          gluon_data[ibin]->SetMarkerColor(kRed); gluon_data[ibin]->SetMarkerStyle(4); gluon_data[ibin]->SetMarkerSize(1.2);
          gluon_data[ibin]->Draw("PL");
          //quark_data[ibin]->SetMarkerColor(kBlue); quark_data[ibin]->SetMarkerStyle(4); quark_data[ibin]->SetMarkerSize(1.2);
          //quark_data[ibin]->Draw("PL");
      }

      if(do_kappa_fitting){
          gluon_data_kappa[ibin]->SetMarkerColor(kRed); gluon_data_kappa[ibin]->SetMarkerStyle(5); gluon_data_kappa[ibin]->SetMarkerSize(1.2);
          gluon_data_kappa[ibin]->Draw("PL");
          //quark_data_kappa[ibin]->SetMarkerColor(kBlue); quark_data_kappa[ibin]->SetMarkerStyle(5); quark_data_kappa[ibin]->SetMarkerSize(1.2);
          //quark_data_kappa[ibin]->Draw("PL");
      }

      tl_g[ibin]->Draw("same"); //tl_q[ibin]->Draw("same");
      
      TString tmp1 = cent_tag_pp[ibin];
      tx1 = new TLatex(); tx1->SetTextSize(.1);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);

      if(ibin==0){
        TLegend *legend = new TLegend(0.35,0.15,0.75,0.4);
        legend ->SetLineColor(kWhite);
        if(do_ptcut_fitting){
            //legend ->AddEntry(quark_data[ibin], "quark (pT cut)", "lepf");
            legend ->AddEntry(gluon_data[ibin], "gluon (pTcut)", "lepf");
        }
        if(do_kappa_fitting){
            //legend ->AddEntry(quark_data_kappa[ibin], "quark (kappa)", "lepf");
            legend ->AddEntry(gluon_data_kappa[ibin], "gluon (kappa)", "lepf");
        }
        legend ->Draw("same");
      }    
      if(ibin==4){
        if(do_kappa_fitting){
          TLegend *legenda = new TLegend(0.25,0.15,0.85,0.3);
          legenda ->SetLineColor(kWhite);
          legenda ->AddEntry((TObject*)0, "x1 : pT > 2 (#kappa=0.5); #kappa=0.3 (pT > 2)", "");
          legenda ->AddEntry((TObject*)0, "x2 : pT > 4 (#kappa=0.5); #kappa=0.5 (pT > 2)", "");
          legenda ->AddEntry((TObject*)0, "x3 : pT > 5 (#kappa=0.5); #kappa=0.7 (pT > 2)", "");
          legenda ->Draw("same");
        }
      }   
    }

    TCanvas *c_data_q_graph = new TCanvas("c_data_q_graph","c_data_q_graph",1800,350);
    drawgraphPad();

    pad_titleg->cd();
    title_tex->DrawLatexNDC(0.5,0.1, data);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    pad_labelg->cd();
    chg_tex->DrawLatexNDC(0.7,0.4, fraction);

    for(int ibin=0; ibin<nCBins; ibin++){
      if(ibin==0) pada->cd(1);
      else pada->cd(6-ibin);
      h_dummy2->Draw("same");

      if(do_ptcut_fitting){
          up_data[ibin]->SetMarkerColor(kAzure+2); up_data[ibin]->SetMarkerStyle(22); up_data[ibin]->SetMarkerSize(1.3);
          up_data[ibin]->Draw("PL");
          down_data[ibin]->SetMarkerColor(kGreen-2); down_data[ibin]->SetMarkerStyle(23); down_data[ibin]->SetMarkerSize(1.3);
          down_data[ibin]->Draw("PL");
      }
      if(do_kappa_fitting){
          up_data_kappa[ibin]->SetMarkerColor(kAzure+2); up_data_kappa[ibin]->SetMarkerStyle(22); up_data_kappa[ibin]->SetMarkerSize(1.3);
          up_data_kappa[ibin]->Draw("PL");
          down_data_kappa[ibin]->SetMarkerColor(kGreen-2); down_data_kappa[ibin]->SetMarkerStyle(23); down_data_kappa[ibin]->SetMarkerSize(1.3);
          down_data_kappa[ibin]->Draw("PL");
      }

      if(ibin==0) {tl_u[ibin]->Draw("same"); tl_d[ibin]->Draw("same");} //tl_ubar[ibin]->Draw("same"); tl_dbar[ibin]->Draw("same"); tl_c[ibin]->Draw("same"); tl_s[ibin]->Draw("same");tl_b[ibin]->Draw("same");}
      else {tl_u_data->Draw("same"); tl_d_data->Draw("same");} //tl_ubar[ibin]->Draw("same"); tl_dbar[ibin]->Draw("same"); tl_c[ibin]->Draw("same"); tl_s[ibin]->Draw("same");tl_b[ibin]->Draw("same");}

      if(ibin==0){
        TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
        legend ->SetLineColor(kWhite);
        legend ->AddEntry(up_data_kappa[ibin], "up", "lepf");
        legend ->AddEntry(down_data_kappa[ibin], "down", "lepf");
        legend ->Draw("same");
      }

          TString tmp1 = cent_tag_pp[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.85, tmp1);
      
    }
  }

/////////////////////final plot///////////////////////////////

  Double_t quark_ncs[nCBins][ntrkbins-1];
  Double_t gluon_ncs[nCBins][ntrkbins-1];

  quark_ncs[0][0] = 0.523287;
  gluon_ncs[0][0] = 0.476713;
  quark_ncs[0][1] = 0.485832;
  gluon_ncs[0][1] = 0.514168;
  quark_ncs[0][2] = 0.466797;
  gluon_ncs[0][2] = 0.533203;
  quark_ncs[1][0] = 0.457592;
  gluon_ncs[1][0] = 0.542408;
  quark_ncs[1][1] = 0.669312;
  gluon_ncs[1][1] = 0.330688;
  quark_ncs[1][2] = 0.739939;
  gluon_ncs[1][2] = 0.260061;
  quark_ncs[2][0] = 0.468311;
  gluon_ncs[2][0] = 0.531689;
  quark_ncs[2][1] = 0.721057;
  gluon_ncs[2][1] = 0.278943;
  quark_ncs[2][2] = 0.806671;
  gluon_ncs[2][2] = 0.193329;
  quark_ncs[3][0] = 0.527101;
  gluon_ncs[3][0] = 0.472899;
  quark_ncs[3][1] = 0.659708; 
  gluon_ncs[3][1] = 0.340292;
  quark_ncs[3][2] = 0.70566;
  gluon_ncs[3][2] = 0.29434;
  quark_ncs[4][0] = 0.556182;
  gluon_ncs[4][0] = 0.443818;
  quark_ncs[4][1] = 0.625937;
  gluon_ncs[4][1] = 0.374063;
  quark_ncs[4][2] = 0.655478;
  gluon_ncs[4][2] = 0.344522;
  quark_ncs[5][0] = 0.552878;
  gluon_ncs[5][0] = 0.447122;
  quark_ncs[5][1] = 0.692227; 
  gluon_ncs[5][1] = 0.307773;
  quark_ncs[5][2] = 0.730981;
  gluon_ncs[5][2] = 0.269019;  

  Double_t quark_unc[nCBins][ntrkbins-1];
  Double_t quark_trkup[nCBins][ntrkbins], gluon_trkup[nCBins][ntrkbins];
  Double_t quark_trkdn[nCBins][ntrkbins], gluon_trkdn[nCBins][ntrkbins];
  Double_t quark_trkpos[nCBins][ntrkbins], gluon_trkpos[nCBins][ntrkbins];
  Double_t quark_trkspillup[nCBins][ntrkbins], gluon_trkspillup[nCBins][ntrkbins];
  Double_t quark_trkspilldn[nCBins][ntrkbins], gluon_trkspilldn[nCBins][ntrkbins];
  Double_t quark_jtptup[nCBins][ntrkbins], gluon_jtptup[nCBins][ntrkbins];
  Double_t quark_jtptdn[nCBins][ntrkbins], gluon_jtptdn[nCBins][ntrkbins];
  Double_t quark_jer[nCBins][ntrkbins], gluon_jer[nCBins][ntrkbins];
  Double_t quark_syst[nCBins][ntrkbins], gluon_syst[nCBins][ntrkbins];
  Double_t quark_syst_corr[nCBins][ntrkbins], gluon_syst_corr[nCBins][ntrkbins];
  Double_t quark0trk_syst_corr[nCBins][ntrkbins], gluon0trk_syst_corr[nCBins][ntrkbins];
  Double_t quark0trk_nocorr[nCBins][ntrkbins], gluon0trk_nocorr[nCBins][ntrkbins];

  Double_t gluon_nosys[nCBins][ntrkbins-1];
  Double_t gluon_withsys[nCBins][ntrkbins-1];
  Double_t gluon_withsys_posneg[nCBins][ntrkbins-1];  
  Double_t gluon_withsys_JER[nCBins][ntrkbins-1];  
  Double_t gluon_withsys_kappa[nCBins][ntrkbins-1];
  Double_t gluon_withsys_kappa_posneg[nCBins][ntrkbins-1];
  Double_t gluon_withsys_kappa_JER[nCBins][ntrkbins-1];

  gluon_nosys[0][0] = 0.571991;
  gluon_nosys[0][1] = 0.569827;
  gluon_nosys[0][2] = 0.572606;
  gluon_nosys[1][0] = 0.565647;
  gluon_nosys[1][1] = 0.591335;
  gluon_nosys[1][2] = 0.596355;
  gluon_nosys[2][0] = 0.591641;
  gluon_nosys[2][1] = 0.566076;
  gluon_nosys[2][2] = 0.595584;
  gluon_nosys[3][0] = 0.625194;
  gluon_nosys[3][1] = 0.6034;
  gluon_nosys[3][2] = 0.594667;
  gluon_nosys[4][0] = 0.565619;
  gluon_nosys[4][1] = 0.62434;
  gluon_nosys[4][2] = 0.62779;

  gluon_withsys[0][0] = 0.584382;
  gluon_withsys[0][1] = 0.584906;
  gluon_withsys[0][2] = 0.589134;
  gluon_withsys[1][0] = 0.589512;
  gluon_withsys[1][1] = 0.590269;
  gluon_withsys[1][2] = 0.597911;
  gluon_withsys[2][0] = 0.602511;
  gluon_withsys[2][1] = 0.567025;
  gluon_withsys[2][2] = 0.593571;
  gluon_withsys[3][0] = 0.555443;
  gluon_withsys[3][1] = 0.55368;
  gluon_withsys[3][2] = 0.573965;
  gluon_withsys[4][0] = 0.544928;
  gluon_withsys[4][1] = 0.599688;
  gluon_withsys[4][2] = 0.606219;

  gluon_withsys_posneg[0][0] = 0.580448;
  gluon_withsys_posneg[0][1] = 0.585204;
  gluon_withsys_posneg[0][2] = 0.587226;
  gluon_withsys_posneg[1][0] = 0.583916;
  gluon_withsys_posneg[1][1] = 0.594779;
  gluon_withsys_posneg[1][2] = 0.600619;
  gluon_withsys_posneg[2][0] = 0.593043;
  gluon_withsys_posneg[2][1] = 0.557135;
  gluon_withsys_posneg[2][2] = 0.594748;
  gluon_withsys_posneg[3][0] = 0.535608;
  gluon_withsys_posneg[3][1] = 0.542154;
  gluon_withsys_posneg[3][2] = 0.565738;
  gluon_withsys_posneg[4][0] = 0.544139;
  gluon_withsys_posneg[4][1] = 0.587163;
  gluon_withsys_posneg[4][2] = 0.609625;

  gluon_withsys_JER[0][0] = 0.594424;
  gluon_withsys_JER[0][1] = 0.584409;
  gluon_withsys_JER[0][2] = 0.588187;
  gluon_withsys_JER[1][0] = 0.598242;
  gluon_withsys_JER[1][1] = 0.595588;
  gluon_withsys_JER[1][2] = 0.614596;
  gluon_withsys_JER[2][0] = 0.616343;
  gluon_withsys_JER[2][1] = 0.577143;
  gluon_withsys_JER[2][2] = 0.606756;
  gluon_withsys_JER[3][0] = 0.565531;
  gluon_withsys_JER[3][1] = 0.543098;
  gluon_withsys_JER[3][2] = 0.564662;
  gluon_withsys_JER[4][0] = 0.592001;
  gluon_withsys_JER[4][1] = 0.600779;
  gluon_withsys_JER[4][2] = 0.603968;

  gluon_withsys_kappa[0][0] = 0.580611;
  gluon_withsys_kappa[0][1] = 0.584382;
  gluon_withsys_kappa[0][2] = 0.592963;
  gluon_withsys_kappa[1][0] = 0.609421;
  gluon_withsys_kappa[1][1] = 0.589512;
  gluon_withsys_kappa[1][2] = 0.621667;
  gluon_withsys_kappa[2][0] = 0.557694;
  gluon_withsys_kappa[2][1] = 0.602511;
  gluon_withsys_kappa[2][2] = 0.596131;
  gluon_withsys_kappa[3][0] = 0.53002;
  gluon_withsys_kappa[3][1] = 0.555443;
  gluon_withsys_kappa[3][2] = 0.610621;
  gluon_withsys_kappa[4][0] = 0.575807;
  gluon_withsys_kappa[4][1] = 0.544928;
  gluon_withsys_kappa[4][2] = 0.548352;

  gluon_withsys_kappa_posneg[0][0] = 0.579907;
  gluon_withsys_kappa_posneg[0][1] = 0.580448;
  gluon_withsys_kappa_posneg[0][2] = 0.589738;
  gluon_withsys_kappa_posneg[1][0] = 0.615658;
  gluon_withsys_kappa_posneg[1][1] = 0.583916;
  gluon_withsys_kappa_posneg[1][2] = 0.618005;
  gluon_withsys_kappa_posneg[2][0] = 0.562723;
  gluon_withsys_kappa_posneg[2][1] = 0.593043;
  gluon_withsys_kappa_posneg[2][2] = 0.587804;
  gluon_withsys_kappa_posneg[3][0] = 0.587908;
  gluon_withsys_kappa_posneg[3][1] = 0.535608;
  gluon_withsys_kappa_posneg[3][2] = 0.601266;
  gluon_withsys_kappa_posneg[4][0] = 0.607357;
  gluon_withsys_kappa_posneg[4][1] = 0.544139;
  gluon_withsys_kappa_posneg[4][2] = 0.532832;

  gluon_withsys_kappa_JER[0][0] = 0.585971;
  gluon_withsys_kappa_JER[0][1] = 0.594424;
  gluon_withsys_kappa_JER[0][2] = 0.609343;
  gluon_withsys_kappa_JER[1][0] = 0.610972;
  gluon_withsys_kappa_JER[1][1] = 0.598242;
  gluon_withsys_kappa_JER[1][2] = 0.6435;
  gluon_withsys_kappa_JER[2][0] = 0.58198;
  gluon_withsys_kappa_JER[2][1] = 0.616343;
  gluon_withsys_kappa_JER[2][2] = 0.604644;
  gluon_withsys_kappa_JER[3][0] = 0.584166;
  gluon_withsys_kappa_JER[3][1] = 0.565531;
  gluon_withsys_kappa_JER[3][2] = 0.621031;
  gluon_withsys_kappa_JER[4][0] = 0.59299;
  gluon_withsys_kappa_JER[4][1] = 0.592001;
  gluon_withsys_kappa_JER[4][2] = 0.573069;

  Double_t gluon_relerr_ptcut_0trk[nCBins][ntrkbins-1];
  Double_t gluon_relerr_ptcut_trkunc[nCBins][ntrkbins-1];
  Double_t gluon_relerr_ptcut_posneg[nCBins][ntrkbins-1];  
  Double_t gluon_relerr_ptcut_JER[nCBins][ntrkbins-1];  

  Double_t gluon_relerr_kappa_0trk[nCBins][ntrkbins-1];
  Double_t gluon_relerr_kappa_trkunc[nCBins][ntrkbins-1];
  Double_t gluon_relerr_kappa_posneg[nCBins][ntrkbins-1];
  Double_t gluon_relerr_kappa_JER[nCBins][ntrkbins-1];


  //ptcut
  gluon_relerr_ptcut_0trk[0][0] = 0.1;   gluon_relerr_ptcut_0trk[0][1] = 0.1;   gluon_relerr_ptcut_0trk[0][2] = 0.2;
  gluon_relerr_ptcut_0trk[4][0] = 0.2;   gluon_relerr_ptcut_0trk[4][1] = 0.5;   gluon_relerr_ptcut_0trk[4][2] = 1.0;
  gluon_relerr_ptcut_0trk[3][0] = 0.4;   gluon_relerr_ptcut_0trk[3][1] = 1.0;   gluon_relerr_ptcut_0trk[3][2] = 2.0;
  gluon_relerr_ptcut_0trk[2][0] = 0.4;   gluon_relerr_ptcut_0trk[2][1] = 1.7;   gluon_relerr_ptcut_0trk[2][2] = 3.0;
  gluon_relerr_ptcut_0trk[1][0] = 0.4;   gluon_relerr_ptcut_0trk[1][1] = 2.4;   gluon_relerr_ptcut_0trk[1][2] = 4.5;

  gluon_relerr_ptcut_trkunc[0][0] = 4.0;   gluon_relerr_ptcut_trkunc[0][1] = 4.0;   gluon_relerr_ptcut_trkunc[0][2] = 4.0;
  gluon_relerr_ptcut_trkunc[4][0] = 10.;   gluon_relerr_ptcut_trkunc[4][1] = 10.;   gluon_relerr_ptcut_trkunc[4][2] = 10.;
  gluon_relerr_ptcut_trkunc[3][0] = 8.0;   gluon_relerr_ptcut_trkunc[3][1] = 8.0;   gluon_relerr_ptcut_trkunc[3][2] = 8.0;
  gluon_relerr_ptcut_trkunc[2][0] = 6.0;   gluon_relerr_ptcut_trkunc[2][1] = 8.0;   gluon_relerr_ptcut_trkunc[2][2] = 8.0;
  gluon_relerr_ptcut_trkunc[1][0] = 5.0;   gluon_relerr_ptcut_trkunc[1][1] = 5.0;   gluon_relerr_ptcut_trkunc[1][2] = 5.0;

  gluon_relerr_ptcut_posneg[0][0] = 0.5;   gluon_relerr_ptcut_posneg[0][1] = 0.5;   gluon_relerr_ptcut_posneg[0][2] = 0.5;
  gluon_relerr_ptcut_posneg[4][0] = 3.5;   gluon_relerr_ptcut_posneg[4][1] = 2.0;   gluon_relerr_ptcut_posneg[4][2] = 1.0;
  gluon_relerr_ptcut_posneg[3][0] = 3.5;   gluon_relerr_ptcut_posneg[3][1] = 2.0;   gluon_relerr_ptcut_posneg[3][2] = 2.0;
  gluon_relerr_ptcut_posneg[2][0] = 1.5;   gluon_relerr_ptcut_posneg[2][1] = 1.5;   gluon_relerr_ptcut_posneg[2][2] = 1.0;
  gluon_relerr_ptcut_posneg[1][0] = 1.0;   gluon_relerr_ptcut_posneg[1][1] = 1.0;   gluon_relerr_ptcut_posneg[1][2] = 1.0;

  gluon_relerr_ptcut_JER[0][0] = 1.0;   gluon_relerr_ptcut_JER[0][1] = 0.5;   gluon_relerr_ptcut_JER[0][2] = 0.5;
  gluon_relerr_ptcut_JER[4][0] = 2.0;   gluon_relerr_ptcut_JER[4][1] = 0.5;   gluon_relerr_ptcut_JER[4][2] = 0.5;
  gluon_relerr_ptcut_JER[3][0] = 2.0;   gluon_relerr_ptcut_JER[3][1] = 1.8;   gluon_relerr_ptcut_JER[3][2] = 1.8;
  gluon_relerr_ptcut_JER[2][0] = 2.0;   gluon_relerr_ptcut_JER[2][1] = 2.0;   gluon_relerr_ptcut_JER[2][2] = 2.0;
  gluon_relerr_ptcut_JER[1][0] = 2.0;   gluon_relerr_ptcut_JER[1][1] = 2.0;   gluon_relerr_ptcut_JER[1][2] = 2.0;

  //kappa
  gluon_relerr_kappa_0trk[0][0] = 0.1;   gluon_relerr_kappa_0trk[0][1] = 0.1;   gluon_relerr_kappa_0trk[0][2] = 0.1;
  gluon_relerr_kappa_0trk[4][0] = 0.2;   gluon_relerr_kappa_0trk[4][1] = 0.2;   gluon_relerr_kappa_0trk[4][2] = 0.2;
  gluon_relerr_kappa_0trk[3][0] = 0.4;   gluon_relerr_kappa_0trk[3][1] = 0.4;   gluon_relerr_kappa_0trk[3][2] = 0.4;
  gluon_relerr_kappa_0trk[2][0] = 0.4;   gluon_relerr_kappa_0trk[2][1] = 0.4;   gluon_relerr_kappa_0trk[2][2] = 0.4;
  gluon_relerr_kappa_0trk[1][0] = 0.4;   gluon_relerr_kappa_0trk[1][1] = 0.4;   gluon_relerr_kappa_0trk[1][2] = 0.4;

  gluon_relerr_kappa_trkunc[0][0] = 4.0;   gluon_relerr_kappa_trkunc[0][1] = 4.0;   gluon_relerr_kappa_trkunc[0][2] = 4.0;
  gluon_relerr_kappa_trkunc[4][0] = 8.0;   gluon_relerr_kappa_trkunc[4][1] = 10.;   gluon_relerr_kappa_trkunc[4][2] = 10.;
  gluon_relerr_kappa_trkunc[3][0] = 6.0;   gluon_relerr_kappa_trkunc[3][1] = 8.0;   gluon_relerr_kappa_trkunc[3][2] = 8.0;
  gluon_relerr_kappa_trkunc[2][0] = 5.0;   gluon_relerr_kappa_trkunc[2][1] = 6.0;   gluon_relerr_kappa_trkunc[2][2] = 8.0;
  gluon_relerr_kappa_trkunc[1][0] = 5.0;   gluon_relerr_kappa_trkunc[1][1] = 5.0;   gluon_relerr_kappa_trkunc[1][2] = 5.0;

  gluon_relerr_kappa_posneg[0][0] = 0.5;   gluon_relerr_kappa_posneg[0][1] = 0.5;   gluon_relerr_kappa_posneg[0][2] = 0.5;
  gluon_relerr_kappa_posneg[4][0] = 4.0;   gluon_relerr_kappa_posneg[4][1] = 3.5;   gluon_relerr_kappa_posneg[4][2] = 2.0;
  gluon_relerr_kappa_posneg[3][0] = 4.0;   gluon_relerr_kappa_posneg[3][1] = 3.5;   gluon_relerr_kappa_posneg[3][2] = 2.0;
  gluon_relerr_kappa_posneg[2][0] = 1.5;   gluon_relerr_kappa_posneg[2][1] = 1.5;   gluon_relerr_kappa_posneg[2][2] = 1.5;
  gluon_relerr_kappa_posneg[1][0] = 1.0;   gluon_relerr_kappa_posneg[1][1] = 1.0;   gluon_relerr_kappa_posneg[1][2] = 1.0;

  gluon_relerr_kappa_JER[0][0] = 1.0;   gluon_relerr_kappa_JER[0][1] = 1.0;   gluon_relerr_kappa_JER[0][2] = 1.5;
  gluon_relerr_kappa_JER[4][0] = 2.0;   gluon_relerr_kappa_JER[4][1] = 2.0;   gluon_relerr_kappa_JER[4][2] = 2.0;
  gluon_relerr_kappa_JER[3][0] = 2.0;   gluon_relerr_kappa_JER[3][1] = 2.0;   gluon_relerr_kappa_JER[3][2] = 1.8;
  gluon_relerr_kappa_JER[2][0] = 3.0;   gluon_relerr_kappa_JER[2][1] = 2.0;   gluon_relerr_kappa_JER[2][2] = 3.0;
  gluon_relerr_kappa_JER[1][0] = 2.0;   gluon_relerr_kappa_JER[1][1] = 2.0;   gluon_relerr_kappa_JER[1][2] = 3.0;

////plotting systematcs

  Double_t gluon_fit_syserr[nCBins][ntrkbins-1];

  if(do_systematics){
    for(int ibin=0; ibin<nCBins; ibin++){
      for(int ibin2=0;ibin2<ntrkbins-1;ibin2++){ 
        if(do_ptcut_fitting){
/*
          h_syst_posneg[ibin]->SetBinContent(ibin2+1,fabs(gluon_withsys_posneg[ibin][ibin2] - gluon_withsys[ibin][ibin2])/gluon_withsys[ibin][ibin2]);
          h_syst_JER[ibin]->SetBinContent(ibin2+1,fabs(gluon_withsys_JER[ibin][ibin2] - gluon_withsys[ibin][ibin2])/gluon_withsys[ibin][ibin2]);
          if(ibin2==0) h_syst_0trk[ibin]->SetBinContent(ibin2+1,diff_0trk_data_MC[ibin][ibin2]/gluon_withsys[ibin][ibin2]);
          else h_syst_0trk[ibin]->SetBinContent(ibin2+1,diff_0trk_data_MC[ibin][ibin2+1]/gluon_withsys[ibin][ibin2]);

          if(ibin != 1)h_syst_trk[ibin]->SetBinContent(ibin2+1,(gluon_fit_high_diff[ibin][ibin2]+gluon_fit_low_diff[ibin][ibin2])/(2*gluon_withsys[ibin][ibin2]));
          else h_syst_trk[ibin]->SetBinContent(ibin2+1,(gluon_fit_high_diff[ibin+1][ibin2]+gluon_fit_low_diff[ibin+1][ibin2])/(2*gluon_withsys[ibin][ibin2]));
*/
          h_syst_posneg[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_ptcut_posneg[ibin][ibin2]);
          h_syst_JER[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_ptcut_JER[ibin][ibin2]);
          h_syst_0trk[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_ptcut_0trk[ibin][ibin2]);
          h_syst_trk[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_ptcut_trkunc[ibin][ibin2]);

          gluon_fit_syserr[ibin][ibin2] = gluon_fit[ibin][ibin2]*sqrt(pow(h_syst_posneg[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_JER[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_trk[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_0trk[ibin]->GetBinContent(ibin2+1),2));
          h_syst[ibin]->SetBinContent(ibin2+1,sqrt(pow(h_syst_posneg[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_JER[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_trk[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_0trk[ibin]->GetBinContent(ibin2+1),2)));

          //cout<<ibin<<". "<<ibin2<<". "<<100.*(gluon_fit_high_diff[ibin][ibin2]+gluon_fit_low_diff[ibin][ibin2])/(2*gluon_withsys[ibin][ibin2])<<endl; 
          //cout<<ibin<<". "<<ibin2<<". "<<100.*fabs(gluon_withsys_posneg[ibin][ibin2] - gluon_withsys[ibin][ibin2])/gluon_withsys[ibin][ibin2]<<endl;
        }
        else if(do_kappa_fitting){
/*
          h_syst_posneg[ibin]->SetBinContent(ibin2+1,fabs(gluon_withsys_kappa_posneg[ibin][ibin2] - gluon_withsys_kappa[ibin][ibin2])/gluon_withsys_kappa[ibin][ibin2]);
          h_syst_JER[ibin]->SetBinContent(ibin2+1,fabs(gluon_withsys_kappa_JER[ibin][ibin2] - gluon_withsys_kappa[ibin][ibin2])/gluon_withsys_kappa[ibin][ibin2]);
          h_syst_0trk[ibin]->SetBinContent(ibin2+1,0.0025);
          h_syst_trk[ibin]->SetBinContent(ibin2+1,(gluon_fit_kappa_high_diff[ibin][ibin2]+gluon_fit_kappa_low_diff[ibin][ibin2])/(2*gluon_withsys_kappa[ibin][ibin2]));
*/
          h_syst_posneg[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_kappa_posneg[ibin][ibin2]);
          h_syst_JER[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_kappa_JER[ibin][ibin2]);
          h_syst_0trk[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_kappa_0trk[ibin][ibin2]);
          h_syst_trk[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_kappa_trkunc[ibin][ibin2]);

          gluon_fit_syserr[ibin][ibin2] = gluon_fit_kappa[ibin][ibin2]*sqrt(pow(h_syst_posneg[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_JER[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_trk[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_0trk[ibin]->GetBinContent(ibin2+1),2));
          h_syst[ibin]->SetBinContent(ibin2+1,sqrt(pow(h_syst_posneg[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_JER[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_trk[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_0trk[ibin]->GetBinContent(ibin2+1),2)));

          //cout<<ibin<<". "<<ibin2<<". "<<100.*(gluon_fit_kappa_high_diff[ibin][ibin2]+gluon_fit_kappa_low_diff[ibin][ibin2])/(2*gluon_withsys_kappa[ibin][ibin2])<<endl; 
          //cout<<ibin<<". "<<ibin2<<". "<<100.*fabs(gluon_withsys_kappa_posneg[ibin][ibin2] - gluon_withsys_kappa[ibin][ibin2])/gluon_withsys_kappa[ibin][ibin2]<<endl; 
        }
        h_syst[ibin]->SetBinError(ibin2+1,0.);
      }
    }

    //data qg graphs
    //if(do_ptcut_fitting){
      for(int ibin=0; ibin<nCBins; ibin++){
        if(do_ptcut_fitting) gluon_data_sys[ibin] = new TGraphErrors(ntrkbins-1,x_mean,gluon_fit[ibin],x_mean_err,gluon_fit_syserr[ibin]);
        if(do_kappa_fitting) gluon_data_sys[ibin] = new TGraphErrors(ntrkbins-1,kappa,gluon_fit_kappa[ibin],kappa_err,gluon_fit_syserr[ibin]);        
        gluon_data_sys[ibin]->SetName((TString)("gluon_data_sys_cent"+cent[ibin]));
      }
    //}

    TCanvas *c_syst = new TCanvas("c_syst","c_syst",1500,350);
    drawgraphPad();

    pad_titleg->cd();
    //title_tex->DrawLatexNDC(0.5,0.1, pythia);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    pad_labelg->cd();
    chg_tex->DrawLatexNDC(0.7,0.4, "relative error");

    for(int ibin=0; ibin<nCBins; ibin++){
      if(ibin==0) pada->cd(1);
      else pada->cd(6-ibin);
      h_cosm(h_syst[ibin]);
      if(do_ptcut_fitting) h_syst[ibin]->GetXaxis()->SetTitle("track pT cut");
      else if(do_kappa_fitting) h_syst[ibin]->GetXaxis()->SetTitle("kappa");

      h_syst[ibin]->GetYaxis()->SetRangeUser(0.,0.25);
      h_syst[ibin]->SetLineColor(kBlack); h_syst[ibin]->SetMarkerColor(kBlack); h_syst[ibin]->SetMarkerStyle(20);
      h_syst[ibin]->SetLineWidth(3);
      h_syst[ibin]->Draw("L same");
      h_syst_posneg[ibin]->SetLineColor(kGreen-2); h_syst_posneg[ibin]->SetMarkerColor(kGreen-2); h_syst_posneg[ibin]->SetMarkerStyle(20);
      h_syst_posneg[ibin]->SetLineWidth(3);
      h_syst_posneg[ibin]->Draw("L same");
      h_syst_JER[ibin]->SetLineColor(kRed-2); h_syst_JER[ibin]->SetMarkerColor(kRed-2); h_syst_JER[ibin]->SetMarkerStyle(20);
      h_syst_JER[ibin]->SetLineWidth(3);
      h_syst_JER[ibin]->Draw("L same");
      h_syst_0trk[ibin]->SetLineColor(kBlue-4); h_syst_0trk[ibin]->SetMarkerColor(kBlue-4); h_syst_0trk[ibin]->SetMarkerStyle(20);
      h_syst_0trk[ibin]->SetLineWidth(3);
      h_syst_0trk[ibin]->Draw("L same");
      h_syst_trk[ibin]->SetLineColor(kOrange); h_syst_trk[ibin]->SetMarkerColor(kOrange); h_syst_trk[ibin]->SetMarkerStyle(20);
      h_syst_trk[ibin]->SetLineWidth(3);
      h_syst_trk[ibin]->Draw("L same");

      TString tmp1 = cent_tag_pp[ibin];
      tx1 = new TLatex(); tx1->SetTextSize(.1);
      tx1->DrawLatexNDC(0.25,0.85, tmp1);

      if(ibin==0){
        TLegend *legend = new TLegend(0.15,0.35,0.75,0.85);
        legend ->SetLineColor(kWhite);
        legend ->AddEntry(h_syst[ibin], "total unc.", "lepf");
        legend ->AddEntry(h_syst_trk[ibin], "tracking", "lepf");
        legend ->AddEntry(h_syst_JER[ibin], "JER", "lepf");
        legend ->AddEntry(h_syst_posneg[ibin], "pos/neg trk", "lepf");
        legend ->AddEntry(h_syst_0trk[ibin], "0 trk jets", "lepf");
        legend ->Draw("same");
      }    
    }
  }


/////plotting final graphs///////////

  for(int ibin=0; ibin<nCBins; ibin++){
    for(int ibin2=0;ibin2<ntrkbins-1;ibin2++){
        //cout<<"centbin: "<<ibin<<". "<<"trkbin: "<<ibin2<<". "<<fabs(gluon_fit[ibin][ibin2] - gluon_withsys[ibin][ibin2])/gluon_withsys[ibin][ibin2]<<endl;

        gluon_fitsyserr[ibin][ibin2] = quark_fitsys[ibin][ibin2]*quark_unc[ibin][ibin2];
        quark_fitsyserr[ibin][ibin2] = quark_fitsys[ibin][ibin2]*quark_unc[ibin][ibin2];  
    }

    gluon_sys[ibin] = new TGraphErrors(ntrkbins-1,x_mean,gluon_fitsys[ibin],x_mean_err,gluon_fitsyserr[ibin]);
    gluon_sys[ibin]->SetName((TString)("gluon_fitsys_cent"+cent[ibin]));
    quark_sys[ibin] = new TGraphErrors(ntrkbins-1,x_mean,quark_fitsys[ibin],x_mean_err,quark_fitsyserr[ibin]);
    quark_sys[ibin]->SetName((TString)("quark_fitsys_cent"+cent[ibin]));

    gluon_ratio_data[ibin] = new TGraphErrors(ntrkbins-1,x_mean,gluon_ratio[ibin],x_mean_err,gluon_fitsyserr[ibin]);
    gluon_ratio_data[ibin]->SetName((TString)("gluon_ratio_data_cent"+cent[ibin]));
    quark_ratio_data[ibin] = new TGraphErrors(ntrkbins-1,x_mean,quark_ratio[ibin],x_mean_err,quark_fitsyserr[ibin]);
    quark_ratio_data[ibin]->SetName((TString)("quark_ratio_data_cent"+cent[ibin]));

    gluon_ratio_reco[ibin] = new TGraphErrors(ntrkbins-1,x_mean,gluon_recofit_ratio[ibin],x_mean_err,gluon_fitsyserr[ibin]);
    gluon_ratio_reco[ibin]->SetName((TString)("gluon_ratio_reco_cent"+cent[ibin]));
    quark_ratio_reco[ibin] = new TGraphErrors(ntrkbins-1,x_mean,quark_recofit_ratio[ibin],x_mean_err,quark_fitsyserr[ibin]);
    quark_ratio_reco[ibin]->SetName((TString)("quark_ratio_reco_cent"+cent[ibin]));

    gluon_ncs_data[ibin] = new TGraph(ntrkbins-1,x_mean,gluon_ncs[ibin]);
    gluon_ncs_data[ibin]->SetName((TString)("gluon_ncs_cent"+cent[ibin]));
    quark_ncs_data[ibin] = new TGraph(ntrkbins-1,x_mean,quark_ncs[ibin]);
    quark_ncs_data[ibin]->SetName((TString)("quark_ncs_cent"+cent[ibin]));

    gluon_trkup_gr[ibin] = new TGraph(ntrkbins-1,x_mean,gluon_trkup[ibin]);
    gluon_trkup_gr[ibin]->SetName((TString)("gluon_trkup_cent"+cent[ibin]));
    quark_trkup_gr[ibin] = new TGraph(ntrkbins-1,x_mean,quark_trkup[ibin]);
    quark_trkup_gr[ibin]->SetName((TString)("quark_trkup_cent"+cent[ibin]));

    gluon_trkdn_gr[ibin] = new TGraph(ntrkbins-1,x_mean,gluon_trkdn[ibin]);
    gluon_trkdn_gr[ibin]->SetName((TString)("gluon_trkdn_cent"+cent[ibin]));
    quark_trkdn_gr[ibin] = new TGraph(ntrkbins-1,x_mean,quark_trkdn[ibin]);
    quark_trkdn_gr[ibin]->SetName((TString)("quark_trkdn_cent"+cent[ibin])); 

    gluon_trkpos_gr[ibin] = new TGraph(ntrkbins-1,x_mean,gluon0trk_syst_corr[ibin]);
    gluon_trkpos_gr[ibin]->SetName((TString)("gluon_trkpos_cent"+cent[ibin]));
    quark_trkpos_gr[ibin] = new TGraph(ntrkbins-1,x_mean,quark_trkpos[ibin]);
    quark_trkpos_gr[ibin]->SetName((TString)("quark_trkpos_cent"+cent[ibin])); 

    gluon_spillup_gr[ibin] = new TGraph(ntrkbins-1,x_mean,gluon0trk_nocorr[ibin]);
    gluon_spillup_gr[ibin]->SetName((TString)("gluon_spillup_cent"+cent[ibin]));
    gluon_spilldn_gr[ibin] = new TGraph(ntrkbins-1,x_mean,gluon_syst_corr[ibin]);
    gluon_spilldn_gr[ibin]->SetName((TString)("gluon_spilldn_cent"+cent[ibin]));
  }

  if(do_data){
      TCanvas *c_data_syserr_graph = new TCanvas("c_data_syserr_graph","c_data_syserr_graph",1500,350);
      drawgraphPad();

      pad_titleg->cd();
      title_tex->DrawLatexNDC(0.5,0.1, data);
      title_tex->DrawLatexNDC(0.05,0.1, cms);  
      pad_labelg->cd();
      chg_tex->DrawLatexNDC(0.7,0.4, fraction);

      for(int ibin=0; ibin<nCBins; ibin++){
        if(ibin==0) pada->cd(1);
        else pada->cd(6-ibin);
        h_dummy->Draw("same");
        //quark_data[ibin]->SetMarkerColor(kBlue); quark_data[ibin]->SetMarkerStyle(4); quark_data[ibin]->SetMarkerSize(1.3);
        //quark_data[ibin]->Draw("PL");
    
        gluon_data_sys[ibin]->SetMarkerColor(kRed); gluon_data_sys[ibin]->SetMarkerStyle(20); gluon_data_sys[ibin]->SetMarkerSize(1.3);
        gluon_data_sys[ibin]->SetLineColor(kRed+2); gluon_data_sys[ibin]->SetFillColorAlpha(kRed+2,0.3);
        gluon_data_sys[ibin]->Draw("a3");

        gluon_data_sys[ibin]->SetTitle("");
        gluon_data_sys[ibin]->GetXaxis()->SetRangeUser(1.,6.);
        gluon_data_sys[ibin]->GetYaxis()->SetRangeUser(0.1,1.);
        gluon_data_sys[ibin]->GetXaxis()->SetTitle("trk p_{T} cut");
        //gluon_data_sys[ibin]->GetXaxis()->SetTitle("kappa");
        gluon_data_sys[ibin]->GetYaxis()->SetTitle("");
        gluon_data_sys[ibin]->GetYaxis()->SetNdivisions(505);
        gluon_data_sys[ibin]->GetYaxis()->SetTitleSize(0.07);
        gluon_data_sys[ibin]->GetYaxis()->SetLabelSize(0.05);
        gluon_data_sys[ibin]->GetXaxis()->CenterTitle();
        gluon_data_sys[ibin]->GetYaxis()->CenterTitle();
        gluon_data_sys[ibin]->GetYaxis()->SetTitleOffset(0.91);
        gluon_data_sys[ibin]->GetXaxis()->SetTitleOffset(0.93);
        gluon_data_sys[ibin]->GetXaxis()->SetNdivisions(505);
        gluon_data_sys[ibin]->GetXaxis()->SetTitleSize(0.05);
        gluon_data_sys[ibin]->GetXaxis()->SetLabelSize(0.05);
        gluon_data_sys[ibin]->Draw("a3 same");
    
    /*
        quark_sys[ibin]->SetMarkerColor(kBlue); quark_sys[ibin]->SetMarkerStyle(20); quark_sys[ibin]->SetMarkerSize(1.3);
        quark_sys[ibin]->SetFillStyle(3001); quark_sys[ibin]->SetFillColorAlpha(kBlue, 0.3); quark_sys[ibin]->SetFillColor(kBlue); quark_sys[ibin]->Draw("2[]");
        gluon_sys[ibin]->SetMarkerColor(kRed); gluon_sys[ibin]->SetMarkerStyle(20); gluon_sys[ibin]->SetMarkerSize(1.3);
        gluon_sys[ibin]->SetFillStyle(3001); gluon_sys[ibin]->SetFillColorAlpha(kRed, 0.3); gluon_sys[ibin]->SetFillColor(kRed); gluon_sys[ibin]->Draw("2[]");
    */

        if(do_ptcut_fitting){
          gluon_data[ibin]->SetMarkerColor(kRed); gluon_data[ibin]->SetLineColor(kRed); gluon_data[ibin]->SetMarkerStyle(4); gluon_data[ibin]->SetMarkerSize(1.3);
          gluon_data[ibin]->Draw("PL");
          tl_g[ibin]->Draw("same"); //tl_q[ibin]->Draw("same");
        }

        if(do_kappa_fitting){
          gluon_data_kappa[ibin]->SetMarkerColor(kRed); gluon_data_kappa[ibin]->SetLineColor(kRed); gluon_data_kappa[ibin]->SetMarkerStyle(4); gluon_data_kappa[ibin]->SetMarkerSize(1.3);
          gluon_data_kappa[ibin]->Draw("PL");
          tl_g_kappa[ibin]->Draw("same"); //tl_q[ibin]->Draw("same");
        } 
        
        TString tmp1 = cent_tag_pp[ibin];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.15,0.85, tmp1);
      
      }
    }

    if(do_data && do_ptcut_fitting){

    //plotting qg final graph
    //Color_t co[6] = {kBlue-4, kMagenta-2, kGreen-2, kRed-3, kOrange+1, kCyan+1};
    Color_t co[6] = {kBlack, kRed-3, kOrange+1, kGreen-2, kBlue-4, kCyan+1};
    //Color_t co_marker[6] = {kBlue+3, kMagenta+3, kGreen+3, kRed+3, kOrange+9, kCyan+4};
    Color_t co_marker[6] = {kGray+1, kRed+2, kOrange+9, kGreen+3, kBlue+2, kCyan+4};
    int markerstyle[6] = {20,47,21,29,34,34};
    //Color_t co[6] = {kRed-2,kRed-1,kRed,kRed+1,kRed+2,kRed+3};

    TMultiGraph *gluon_mg = new TMultiGraph();
    gluon_mg->SetTitle("");
    TMultiGraph *quark_mg = new TMultiGraph();
    quark_mg->SetTitle("");

    gluon_sys[1]->SetMarkerColor(co_marker[1]); gluon_sys[1]->SetMarkerStyle(markerstyle[1]); gluon_sys[1]->SetMarkerSize(1.);
    gluon_sys[1]->SetFillStyle(1001); gluon_sys[1]->SetFillColorAlpha(co[1], 0.7);
    gluon_mg->Add(gluon_sys[1]);    
    gluon_sys[2]->SetMarkerColor(co_marker[2]); gluon_sys[2]->SetMarkerStyle(markerstyle[2]); gluon_sys[2]->SetMarkerSize(1.);
    gluon_sys[2]->SetFillStyle(1001); gluon_sys[2]->SetFillColorAlpha(co[2], 0.7);
    gluon_mg->Add(gluon_sys[2]);
    gluon_sys[3]->SetMarkerColor(co_marker[3]); gluon_sys[3]->SetMarkerStyle(markerstyle[3]); gluon_sys[3]->SetMarkerSize(1.);
    gluon_sys[3]->SetFillStyle(1001); gluon_sys[3]->SetFillColorAlpha(co[3], 0.7);
    gluon_mg->Add(gluon_sys[3]);
    gluon_sys[4]->SetMarkerColor(co_marker[4]); gluon_sys[4]->SetMarkerStyle(markerstyle[4]); gluon_sys[4]->SetMarkerSize(1.);
    gluon_sys[4]->SetFillStyle(1001); gluon_sys[4]->SetFillColorAlpha(co[4], 0.7);
    gluon_mg->Add(gluon_sys[4]);

    gluon_sys[0]->SetMarkerColor(co_marker[0]); gluon_sys[0]->SetMarkerStyle(markerstyle[0]); gluon_sys[0]->SetMarkerSize(1.);
    gluon_sys[0]->SetFillStyle(1001); gluon_sys[0]->SetFillColorAlpha(co[0], 0.7);
    gluon_mg->Add(gluon_sys[0]);

    quark_sys[1]->SetMarkerColor(co_marker[1]); quark_sys[1]->SetMarkerStyle(markerstyle[1]); quark_sys[1]->SetMarkerSize(1.);
    quark_sys[1]->SetFillStyle(1001); quark_sys[1]->SetFillColorAlpha(co[1], 0.7);
    quark_mg->Add(quark_sys[1]);    
    quark_sys[2]->SetMarkerColor(co_marker[2]); quark_sys[2]->SetMarkerStyle(markerstyle[2]); quark_sys[2]->SetMarkerSize(1.);
    quark_sys[2]->SetFillStyle(1001); quark_sys[2]->SetFillColorAlpha(co[2], 0.7);
    quark_mg->Add(quark_sys[2]);
    quark_sys[3]->SetMarkerColor(co_marker[3]); quark_sys[3]->SetMarkerStyle(markerstyle[3]); quark_sys[3]->SetMarkerSize(1.);
    quark_sys[3]->SetFillStyle(1001); quark_sys[3]->SetFillColorAlpha(co[3], 0.7);
    quark_mg->Add(quark_sys[3]);
    quark_sys[4]->SetMarkerColor(co_marker[4]); quark_sys[4]->SetMarkerStyle(markerstyle[4]); quark_sys[4]->SetMarkerSize(1.);
    quark_sys[4]->SetFillStyle(1001); quark_sys[4]->SetFillColorAlpha(co[4], 0.7);
    quark_mg->Add(quark_sys[4]);

    quark_sys[0]->SetMarkerColor(co_marker[0]); quark_sys[0]->SetMarkerStyle(markerstyle[0]); quark_sys[0]->SetMarkerSize(1.);
    quark_sys[0]->SetFillStyle(1001); quark_sys[0]->SetFillColorAlpha(co[0], 0.7);
    quark_mg->Add(quark_sys[0]);


    TCanvas *c_qg_chg_graph = new TCanvas("c_qg_chg_graph","c_qg_chg_graph",650,900);
    c_qg_chg_graph->Divide(1,1,0);

    TPad *pad_title = new TPad("pad_title", "",0.,0.95,1.,1.);
    TPad *pad1 = new TPad("pad1", "",0.,0.35,1.,0.98);
    pad1->Divide(1,1);  
    pad1->Draw();
    TPad *pad2 = new TPad("pad2", "",0.,0.,1.,0.35);
    pad2->Divide(1,1);  
    pad2->Draw();
    pad_title->Draw();
    pad_title->cd();
    title_tex->DrawLatexNDC(0.3,0.1, data);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    /*
      pad1->cd(1);
      //h_dummy7->Draw("same");
      quark_mg->Draw("a3");  
      quark_mg->GetXaxis()->SetRangeUser(1.,6.);
      quark_mg->GetYaxis()->SetRangeUser(0.1,1.1);
      quark_mg->GetXaxis()->SetTitle("Track p_{T} cut (GeV)");
      //quark_mg->GetXaxis()->SetTitle("kappa");
      quark_mg->GetYaxis()->SetTitle("Quark jet fraction");
      quark_mg->GetYaxis()->SetNdivisions(505);
      quark_mg->GetYaxis()->SetTitleSize(0.07);
      quark_mg->GetYaxis()->SetLabelSize(0.05);
      quark_mg->GetXaxis()->CenterTitle();
      quark_mg->GetYaxis()->CenterTitle();
      quark_mg->GetYaxis()->SetTitleOffset(0.91);
      quark_mg->GetXaxis()->SetTitleOffset(0.93);
      quark_mg->GetXaxis()->SetNdivisions(505);
      quark_mg->GetXaxis()->SetTitleSize(0.05);
      quark_mg->GetXaxis()->SetLabelSize(0.05);
      quark_mg->Draw("a3");
      */
      pad1->cd(1);
      //h_dummy8->Draw("same");
      gluon_mg->Draw("a3");  
      gluon_mg->GetXaxis()->SetRangeUser(1.,6.);
      //gluon_mg->GetXaxis()->SetRangeUser(0.3,2.);
      gluon_mg->GetYaxis()->SetRangeUser(0.3,0.9);
      gluon_mg->GetXaxis()->SetTitle("Track p_{T} cut (GeV)");
      //gluon_mg->GetXaxis()->SetTitle("kappa");
      gluon_mg->GetYaxis()->SetTitle("Gluon jet fraction");
      gluon_mg->GetYaxis()->SetNdivisions(505);
      gluon_mg->GetYaxis()->SetTitleSize(0.07);
      gluon_mg->GetYaxis()->SetLabelSize(0.05);
      gluon_mg->GetXaxis()->CenterTitle();
      gluon_mg->GetYaxis()->CenterTitle();
      gluon_mg->GetYaxis()->SetTitleOffset(0.91);
      gluon_mg->GetXaxis()->SetTitleOffset(0.93);
      gluon_mg->GetXaxis()->SetNdivisions(505);
      gluon_mg->GetXaxis()->SetTitleSize(0.05);
      gluon_mg->GetXaxis()->SetLabelSize(0.05);
      gluon_mg->Draw("a3");


    for(int ibin=0; ibin<nCBins; ibin++){
      /*
      pad1->cd(1);
      //h_dummy7->Draw("same");
      quark_data[ibin]->SetMarkerColor(co_marker[ibin]); quark_data[ibin]->SetMarkerStyle(markerstyle[ibin]); quark_data[ibin]->SetMarkerSize(2.);
      quark_data[ibin]->Draw("P");
      */
      pad1->cd(1);
      h_dummy8->Draw("same");
      //gluon_data[ibin]->SetMarkerColor(co_marker[ibin]); gluon_data[ibin]->SetMarkerStyle(markerstyle[ibin]); gluon_data[ibin]->SetMarkerSize(2.);
      gluon_data[ibin]->Draw("P");

      if(ibin>0){
      /*
      pad2->cd(1);
      h_dummy7->Draw("same");
      quark_ratio_data[ibin]->SetMarkerColor(co_marker[ibin]); quark_ratio_data[ibin]->SetMarkerStyle(markerstyle[ibin]); quark_ratio_data[ibin]->SetMarkerSize(2.);
      quark_ratio_data[ibin]->Draw("PL");
      quark_ratio_data[ibin]->GetXaxis()->SetRangeUser(1.8,5.2);
      quark_ratio_data[ibin]->GetXaxis()->SetRangeUser(0.,1.);
      //quark_ratio_data[ibin]->GetYaxis()->SetRangeUser(0.3,1.9);
      quark_ratio_data[ibin]->Draw("PL");
      */ 
      pad2->cd(1);
      h_dummy8->Draw("same");
      //gluon_ratio_data[ibin]->SetMarkerColor(co_marker[ibin]); gluon_ratio_data[ibin]->SetMarkerStyle(markerstyle[ibin]); gluon_ratio_data[ibin]->SetMarkerSize(2.);
      gluon_ratio_data[ibin]->Draw("PL");    
      gluon_ratio_data[ibin]->GetXaxis()->SetRangeUser(1.8,5.2);
      gluon_ratio_data[ibin]->GetXaxis()->SetRangeUser(0.,1.);
      //gluon_ratio_data[ibin]->GetYaxis()->SetRangeUser(0.3,1.9);
      gluon_ratio_data[ibin]->Draw("PL");
      }
    }

    pad1->cd(1);
    TLegend *legend2 = new TLegend(0.4,0.8,0.7,0.87);
    legend2 ->SetLineColor(kWhite);
    legend2 ->AddEntry(quark_sys[0], "pp ref", "pf");
    legend2 ->Draw("same"); 

    TLegend *legend3 = new TLegend(0.15,0.65,0.55,0.8);
    legend3 ->SetLineColor(kWhite);
    legend3 ->AddEntry(quark_sys[4], "50-100% PbPb", "pf");
    legend3 ->AddEntry(quark_sys[3], "30-50% PbPb", "pf");
    legend3 ->Draw("same");

    TLegend *legend4 = new TLegend(0.55,0.65,0.85,0.8);
    legend4 ->SetLineColor(kWhite);
    legend4 ->AddEntry(quark_sys[2], "10-30% PbPb", "pf");
    legend4 ->AddEntry(quark_sys[1], "0-10% PbPb", "pf");
    legend4 ->Draw("same");

    pad2->cd(1);
    TLegend *legend6 = new TLegend(0.15,0.65,0.5,0.85);
    legend6 ->SetLineColor(kWhite);
    legend6 ->AddEntry(quark_sys[4], "50-100% PbPb", "p");
    legend6 ->AddEntry(quark_sys[3], "30-50% PbPb", "p");
    legend6 ->Draw("same");

    TLegend *legend7 = new TLegend(0.5,0.65,0.85,0.85);
    legend7 ->SetLineColor(kWhite);
    legend7 ->AddEntry(quark_sys[2], "10-30% PbPb", "p");
    legend7 ->AddEntry(quark_sys[1], "0-10% PbPb", "p");
    legend7 ->Draw("same");

  }

  if(draw_tracking){
  /////tracking closures

    TF1 *f_inc = new TF1("f_inc","pol0",2.,10.);
    TF1 *f_pos = new TF1("f_pos","pol0",2.,10.);
    TF1 *f_neg = new TF1("f_neg","pol0",2.,10.);

    TCanvas *c_trk_pos = new TCanvas("c_trk_pos","c_trk_pos",1500,600);

      drawPad();
      pad_titlex->cd();
      //title_tex->DrawLatexNDC(0.4,0.1, data);
      title_tex->DrawLatexNDC(0.05,0.1, cms);
      pad_labelx->cd();
      ylabel_tex->DrawLatexNDC(0.7,0.6, "Entries");
      ylabel_tex->DrawLatexNDC(0.7,0.1, "Reco / Gen");

      for(int ibin=0;ibin<nCBins;ibin++){
    
        if(ibin==0) padx->cd(1);
        else padx->cd(6-ibin);   

        padx->SetMargin(10.,0.,0.,0.);      
        h_trkpt_MC_full[ibin]->Rebin(5);
        h_trkpt_MC_pos[ibin]->Rebin(5);
        h_trkpt_MC_neg[ibin]->Rebin(5);

        h_genpt_MC_full[ibin]->Rebin(5);
        h_genpt_MC_pos[ibin]->Rebin(5);
        h_genpt_MC_neg[ibin]->Rebin(5);

        h_trkpt_MC_full[ibin]->GetXaxis()->SetTitle("trk p_{T}");
        h_cosm(h_trkpt_MC_full[ibin]);
        h_trkpt_MC_full[ibin]->GetXaxis()->SetRangeUser(2.,10.);
        //h_trkpt_MC_full[ibin]->GetYaxis()->SetRangeUser(0.,0.22);

        h_trkpt_MC_full[ibin]->SetLineColor(kBlack); h_trkpt_MC_full[ibin]->SetMarkerColor(kBlack); h_trkpt_MC_full[ibin]->SetMarkerStyle(20);
        h_trkpt_MC_pos[ibin]->SetLineColor(kBlue); h_trkpt_MC_pos[ibin]->SetMarkerColor(kBlue); h_trkpt_MC_pos[ibin]->SetMarkerStyle(20);    
        h_trkpt_MC_neg[ibin]->SetLineColor(kRed); h_trkpt_MC_neg[ibin]->SetMarkerColor(kRed); h_trkpt_MC_neg[ibin]->SetMarkerStyle(20); 

        h_trkpt_MC_full[ibin]->Draw("e0 same");
        h_trkpt_MC_pos[ibin]->Draw("e0 same");
        h_trkpt_MC_neg[ibin]->Draw("e0 same");

        TString tmp1 = cent_tag_pp_MC[ibin];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.15,0.85, tmp1);
        
        if(ibin==0){
          TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
          legendx ->SetLineColor(kWhite);
          legendx ->AddEntry(h_trkpt_MC_full[ibin], "all tracks", "lep");
          legendx ->AddEntry(h_trkpt_MC_pos[ibin], "positive tracks", "lepf");
          legendx ->AddEntry(h_trkpt_MC_neg[ibin], "negative tracks", "lepf");
          legendx ->Draw("same");
        }
        
        if(ibin==0) pady->cd(1);
        else pady->cd(6-ibin);   

        h_trkpt_MC_full_ratio[ibin] = (TH1D*)h_trkpt_MC_full[ibin]->Clone((TString)("h_trkpt_MC_full"+cent[ibin]));
        h_trkpt_MC_full_ratio[ibin]->Divide(h_genpt_MC_full[ibin]);
        h_trkpt_MC_pos_ratio[ibin] = (TH1D*)h_trkpt_MC_pos[ibin]->Clone((TString)("h_trkpt_MC_pos"+cent[ibin]));
        h_trkpt_MC_pos_ratio[ibin]->Divide(h_genpt_MC_pos[ibin]);
        h_trkpt_MC_neg_ratio[ibin] = (TH1D*)h_trkpt_MC_neg[ibin]->Clone((TString)("h_trkpt_MC_neg"+cent[ibin]));
        h_trkpt_MC_neg_ratio[ibin]->Divide(h_genpt_MC_neg[ibin]);
        h_trkpt_MC_full_ratio[ibin]->GetXaxis()->SetTitle("");
        h_trkpt_MC_full_ratio[ibin]->GetYaxis()->SetTitle("");
        h_cosm(h_trkpt_MC_full_ratio[ibin]);
        h_trkpt_MC_full_ratio[ibin]->GetXaxis()->SetRangeUser(2.,10.);
        h_trkpt_MC_full_ratio[ibin]->GetYaxis()->SetRangeUser(0.8,1.2);
        h_trkpt_MC_full_ratio[ibin]->SetLineColor(kBlack); h_trkpt_MC_full_ratio[ibin]->SetMarkerColor(kBlack);
        h_trkpt_MC_full_ratio[ibin]->SetMarkerStyle(21); h_trkpt_MC_full_ratio[ibin]->Draw("e0 same");
        h_trkpt_MC_pos_ratio[ibin]->SetLineColor(kBlue); h_trkpt_MC_pos_ratio[ibin]->SetMarkerColor(kBlue);
        h_trkpt_MC_pos_ratio[ibin]->SetMarkerStyle(4); h_trkpt_MC_pos_ratio[ibin]->Draw("e0 same");
        h_trkpt_MC_neg_ratio[ibin]->SetLineColor(kRed); h_trkpt_MC_neg_ratio[ibin]->SetMarkerColor(kRed);
        h_trkpt_MC_neg_ratio[ibin]->SetMarkerStyle(4); h_trkpt_MC_neg_ratio[ibin]->Draw("e0 same");
        h_trkpt_MC_full_ratio[ibin]->Fit(f_inc,"Q M R","sames",5.,10.);
        h_trkpt_MC_pos_ratio[ibin]->Fit(f_pos,"Q M R","sames",5.,10.);
        h_trkpt_MC_neg_ratio[ibin]->Fit(f_neg,"Q M R","sames",5.,10.);
        drawline(2.,0.95,10.,0.95);
        drawline(2.,1.,10.,1.);
        drawline(2.,1.05,10.,1.05);

        cout<<ibin<<"  : "<<"trk pos diff:  "<<100*(f_pos->GetParameter(0)-f_inc->GetParameter(0))<<"%"<<endl;
        cout<<ibin<<"  : "<<"trk neg diff:  "<<100*(f_neg->GetParameter(0)-f_inc->GetParameter(0))<<"%"<<endl;

        if(ibin==0){
          TLegend *legendy = new TLegend(0.65,0.55,0.98,0.98);
          legendy ->SetLineColor(kWhite);
          legendy ->AddEntry(h_trkpt_MC_full_ratio[ibin], "all tracks", "lep");
          legendy ->AddEntry(h_trkpt_MC_pos_ratio[ibin], "positive tracks", "lepf");
          legendy ->AddEntry(h_trkpt_MC_neg_ratio[ibin], "negative tracks", "lepf");
          legendy ->Draw("same");
        }
      }

    TCanvas *c_trk_sig_bkg = new TCanvas("c_trk_sig_bkg","c_trk_sig_bkg",1500,300);
    c_trk_sig_bkg->Divide(5,1,0);

    for(int ibin=0;ibin<nCBins;ibin++){  
        if(ibin==0) c_trk_sig_bkg->cd(1);
        else c_trk_sig_bkg->cd(6-ibin);   
        
        h_trkpt_MC[ibin][0]->Rebin(10);
        h_bkgtrkpt_MC[ibin][0]->Rebin(10);
        h_trkpt_data[ibin][0]->Rebin(10);
        h_bkgtrkpt_data[ibin][0]->Rebin(10);  
        
        h_bkgtrkpt_MC[ibin][0]->Scale(1./h_trkpt_MC[ibin][0]->Integral());
        h_trkpt_MC[ibin][0]->Scale(1./h_trkpt_MC[ibin][0]->Integral());
        h_bkgtrkpt_data[ibin][0]->Scale(1./h_trkpt_data[ibin][0]->Integral());
        h_trkpt_data[ibin][0]->Scale(1./h_trkpt_data[ibin][0]->Integral());

        //gPad->SetLogy();
        h_trkpt_MC[ibin][0]->GetXaxis()->SetTitle("trk p_{T}");
        h_cosm(h_trkpt_MC[ibin][0]);
        h_trkpt_MC[ibin][0]->GetXaxis()->SetRangeUser(1.,10.);
        h_trkpt_MC[ibin][0]->GetYaxis()->SetRangeUser(0.,0.55);

        h_trkpt_MC[ibin][0]->SetLineColor(kBlack); h_trkpt_MC[ibin][0]->SetMarkerColor(kBlack); h_trkpt_MC[ibin][0]->SetMarkerStyle(20);
        h_bkgtrkpt_MC[ibin][0]->SetLineColor(kBlack); h_bkgtrkpt_MC[ibin][0]->SetMarkerColor(kBlack); h_bkgtrkpt_MC[ibin][0]->SetMarkerStyle(4);      

        h_trkpt_data[ibin][0]->SetLineColor(kMagenta); h_trkpt_data[ibin][0]->SetMarkerColor(kMagenta); h_trkpt_data[ibin][0]->SetMarkerStyle(20);
        h_bkgtrkpt_data[ibin][0]->SetLineColor(kMagenta); h_bkgtrkpt_data[ibin][0]->SetMarkerColor(kMagenta); h_bkgtrkpt_data[ibin][0]->SetMarkerStyle(4);      

        h_trkpt_MC[ibin][0]->Draw("e0 same");
        h_bkgtrkpt_MC[ibin][0]->Draw("e0 same");

        h_trkpt_data[ibin][0]->Draw("e0 same");
        h_bkgtrkpt_data[ibin][0]->Draw("e0 same");

        TString tmp1 = cent_tag_pp_MC[ibin];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.15,0.85, tmp1);

        if(ibin==0){
          TLegend *legendy = new TLegend(0.65,0.55,0.98,0.98);
          legendy ->SetLineColor(kWhite);
          legendy ->AddEntry(h_trkpt_MC[ibin][0], "raw signal (MC)", "lep");
          legendy ->AddEntry(h_bkgtrkpt_MC[ibin][0], "#eta refl bkg (MC)", "lep");
          legendy ->AddEntry(h_trkpt_data[ibin][0], "raw signal (data)", "lep");
          legendy ->AddEntry(h_bkgtrkpt_data[ibin][0], "#eta refl bkg (data)", "lep");
          legendy ->Draw("same");
        }
    }

    TCanvas *c_trk_sig_bkg_kappa_MC = new TCanvas("c_trk_sig_bkg_kappa_MC","c_trk_sig_bkg_kappa_MC",1500,600);
    c_trk_sig_bkg_kappa_MC->Divide(5,2,0);

    for(int ibin=0;ibin<nCBins;ibin++){  
        for(int ibin2=0;ibin2<ntrkbins;ibin2++){  
          if(ibin==0) c_trk_sig_bkg_kappa_MC->cd(1);
          else c_trk_sig_bkg_kappa_MC->cd(6-ibin);   
          
          h_trkk_MC[ibin][ibin2]->Rebin(10);
          h_bkgtrkk_MC[ibin][ibin2]->Rebin(10);

          h_bkgtrkk_MC[ibin][ibin2]->Scale(1./h_trkk_MC[ibin][ibin2]->Integral());
          h_trkk_MC[ibin][ibin2]->Scale(1./h_trkk_MC[ibin][ibin2]->Integral());

          //gPad->SetLogy();
          h_trkk_MC[ibin][ibin2]->GetXaxis()->SetTitle("trk p_{T}");
          h_cosm(h_trkk_MC[ibin][ibin2]);
          h_trkk_MC[ibin][ibin2]->GetXaxis()->SetRangeUser(2.,10.);
          h_trkk_MC[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.3);

          h_trkk_MC[ibin][ibin2]->SetLineColor(ibin2+1); h_trkk_MC[ibin][ibin2]->SetMarkerColor(ibin2+1); h_trkk_MC[ibin][ibin2]->SetMarkerStyle(20);
          h_bkgtrkk_MC[ibin][ibin2]->SetLineColor(ibin2+1); h_bkgtrkk_MC[ibin][ibin2]->SetMarkerColor(ibin2+1); h_bkgtrkk_MC[ibin][ibin2]->SetMarkerStyle(4);      

          h_trkk_MC[ibin][ibin2]->Draw("e0 same");
          h_bkgtrkk_MC[ibin][ibin2]->Draw("e0 same");

          TString tmp1 = cent_tag_pp_MC[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.85, tmp1);

          if(ibin==0) c_trk_sig_bkg_kappa_MC->cd(6);
          else c_trk_sig_bkg_kappa_MC->cd(11-ibin);   

          //gPad->SetLogy();
          h_trkk_MC_sub[ibin][ibin2] = (TH1D*)h_trkk_MC[ibin][ibin2]->Clone((TString)("h_trkk_MC_sub"+cent[ibin]+"_"+kappa_str[ibin2]));
          h_trkk_MC_sub[ibin][ibin2]->Add(h_bkgtrkk_MC[ibin][ibin2],-1);
          h_trkk_MC_sub[ibin][ibin2]->GetXaxis()->SetRangeUser(2.,10.);
          h_trkk_MC_sub[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.3);
          h_trkk_MC_sub[ibin][ibin2]->Draw("e0 same");
        }
    }

    TCanvas *c_trk_sig_bkg_kappa_data = new TCanvas("c_trk_sig_bkg_kappa_data","c_trk_sig_bkg_kappa_data",1500,600);
    c_trk_sig_bkg_kappa_data->Divide(5,2,0);

    for(int ibin=0;ibin<nCBins;ibin++){  
        for(int ibin2=0;ibin2<ntrkbins;ibin2++){  
          if(ibin==0) c_trk_sig_bkg_kappa_data->cd(1);
          else c_trk_sig_bkg_kappa_data->cd(6-ibin);   
          
          h_trkk_data[ibin][ibin2]->Rebin(10);
          h_bkgtrkk_data[ibin][ibin2]->Rebin(10);

          h_bkgtrkk_data[ibin][ibin2]->Scale(1./h_trkk_data[ibin][ibin2]->Integral());
          h_trkk_data[ibin][ibin2]->Scale(1./h_trkk_data[ibin][ibin2]->Integral());

          //gPad->SetLogy();
          h_trkk_data[ibin][ibin2]->GetXaxis()->SetTitle("trk p_{T}");
          h_cosm(h_trkk_data[ibin][ibin2]);
          h_trkk_data[ibin][ibin2]->GetXaxis()->SetRangeUser(2.,20.);
          h_trkk_data[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.4);

          h_trkk_data[ibin][ibin2]->SetLineColor(ibin2+1); h_trkk_data[ibin][ibin2]->SetMarkerColor(ibin2+1); h_trkk_data[ibin][ibin2]->SetMarkerStyle(20);
          h_bkgtrkk_data[ibin][ibin2]->SetLineColor(ibin2+1); h_bkgtrkk_data[ibin][ibin2]->SetMarkerColor(ibin2+1); h_bkgtrkk_data[ibin][ibin2]->SetMarkerStyle(4);      

          h_trkk_data[ibin][ibin2]->Draw("e0 same");
          h_bkgtrkk_data[ibin][ibin2]->Draw("e0 same");

          TString tmp1 = cent_tag_pp[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.85, tmp1);

          if(ibin==0){
            TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
            legendx ->SetLineColor(kWhite);
            legendx ->AddEntry(h_trkk_data[ibin][0], "Raw signal", "lep");
            legendx ->AddEntry(h_bkgtrkk_data[ibin][0], "#eta refl bkg", "lep");
            legendx ->Draw("same");
          }

          if(ibin==0) c_trk_sig_bkg_kappa_data->cd(6);
          else c_trk_sig_bkg_kappa_data->cd(11-ibin);   

          //gPad->SetLogy();
          h_trkk_data_sub[ibin][ibin2] = (TH1D*)h_trkk_data[ibin][ibin2]->Clone((TString)("h_trkk_data_sub"+cent[ibin]+"_"+kappa_str[ibin2]));
          h_trkk_data_sub[ibin][ibin2]->Add(h_bkgtrkk_data[ibin][ibin2],-1);
          h_trkk_data_sub[ibin][ibin2]->GetXaxis()->SetRangeUser(2.,20.);
          h_trkk_data_sub[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.4);
          h_trkk_data_sub[ibin][ibin2]->Draw("e0 same");
          if(ibin==0) cout<<ibin2<<"  "<<h_trkk_data_sub[ibin][ibin2]->GetMean()<<endl;
          
          TLegend *legendy = new TLegend(0.15,0.55,0.98,0.98);
          legendy ->SetLineColor(kWhite);
          legendy ->AddEntry(h_trkk_data[ibin][0], Form("#kappa = 1. (#mu = %.1f)",h_trkk_data[ibin][0]->GetMean()), "lep");
          legendy ->AddEntry(h_trkk_data[ibin][1], Form("#kappa = 1.5 (#mu = %.1f)",h_trkk_data[ibin][1]->GetMean()), "lep");
          legendy ->AddEntry(h_trkk_data[ibin][2], Form("#kappa = 2. (#mu = %.1f)",h_trkk_data[ibin][2]->GetMean()), "lep");
          legendy ->AddEntry(h_trkk_data[ibin][3], Form("#kappa = 2.5 (#mu = %.1f)",h_trkk_data[ibin][3]->GetMean()), "lep");
          legendy ->Draw("same");
          
        }
    }

    TCanvas *c_trk_sig_bkg_pt_data = new TCanvas("c_trk_sig_bkg_pt_data","c_trk_sig_bkg_pt_data",1500,600);
    c_trk_sig_bkg_pt_data->Divide(5,2,0);

    for(int ibin=0;ibin<nCBins;ibin++){  
        for(int ibin2=1;ibin2<ntrkbins;ibin2++){  
          if(ibin==0) c_trk_sig_bkg_pt_data->cd(1);
          else c_trk_sig_bkg_pt_data->cd(6-ibin);   
          
          h_trkpt_data[ibin][ibin2]->Rebin(10);
          h_bkgtrkpt_data[ibin][ibin2]->Rebin(10);

          double integral_trk = h_trkpt_data[ibin][ibin2]->Integral();
          h_bkgtrkpt_data[ibin][ibin2]->Scale(1./integral_trk);
          h_trkpt_data[ibin][ibin2]->Scale(1./integral_trk);

          //gPad->SetLogy();
          h_trkpt_data[ibin][ibin2]->GetXaxis()->SetTitle("trk p_{T}");
          h_cosm(h_trkpt_data[ibin][ibin2]);
          h_trkpt_data[ibin][ibin2]->GetXaxis()->SetRangeUser(2.,20.);
          h_trkpt_data[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.5);

          h_trkpt_data[ibin][ibin2]->SetLineColor(ibin2+1); h_trkpt_data[ibin][ibin2]->SetMarkerColor(ibin2+1); h_trkpt_data[ibin][ibin2]->SetMarkerStyle(20);
          h_bkgtrkpt_data[ibin][ibin2]->SetLineColor(ibin2+1); h_bkgtrkpt_data[ibin][ibin2]->SetMarkerColor(ibin2+1); h_bkgtrkpt_data[ibin][ibin2]->SetMarkerStyle(4);      
   
          if(ibin2==1){  
            h_trkpt_data[ibin][ibin2]->Draw("e0 same");
            h_bkgtrkpt_data[ibin][ibin2]->Draw("e0 same");
          }

          TString tmp1 = cent_tag_pp[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.85, tmp1);

          if(ibin==0){
            TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
            legendx ->SetLineColor(kWhite);
            legendx ->AddEntry(h_trkpt_data[ibin][1], "Raw signal", "lep");
            legendx ->AddEntry(h_bkgtrkpt_data[ibin][1], "#eta refl bkg", "lep");
            legendx ->Draw("same");
          }

          if(ibin==0) c_trk_sig_bkg_pt_data->cd(6);
          else c_trk_sig_bkg_pt_data->cd(11-ibin);   

          //gPad->SetLogy();
          h_trkpt_data_sub[ibin][ibin2] = (TH1D*)h_trkpt_data[ibin][ibin2]->Clone((TString)("h_trkpt_data_sub"+cent[ibin]+"_"+kappa_str[ibin2]));
          h_trkpt_data_sub[ibin][ibin2]->Add(h_bkgtrkpt_data[ibin][ibin2],-1);
          h_trkpt_data_sub[ibin][ibin2]->GetXaxis()->SetRangeUser(2.,20.);
          h_trkpt_data_sub[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.5);
          if(ibin2==1) h_trkpt_data_sub[ibin][ibin2]->Draw("e0 same");
          if(ibin==0) cout<<ibin2<<"  "<<h_trkpt_data_sub[ibin][ibin2]->GetMean()<<endl;
          TLegend *legendy = new TLegend(0.15,0.55,0.98,0.98);
          legendy ->SetLineColor(kWhite);
          legendy ->AddEntry((TObject*)0, Form("trk p_{T} > 2 GeV (#mu = %.1f)",h_trkpt_data[ibin][1]->GetMean()), "");
          legendy ->AddEntry((TObject*)0, Form("trk p_{T} > 4 GeV (#mu = %.1f)",h_trkpt_data[ibin][2]->GetMean()), "");
          legendy ->AddEntry((TObject*)0, Form("trk p_{T} > 5 GeV (#mu = %.1f)",h_trkpt_data[ibin][3]->GetMean()), "");
          legendy ->Draw("same");
        }
    }

    TCanvas *c_ntrk_jet = new TCanvas("c_ntrk_jet","c_ntrk_jet",1500,300);
    c_ntrk_jet->Divide(5,1,0);

    for(int ibin=0;ibin<nCBins;ibin++){  
        if(ibin==0) c_ntrk_jet->cd(1);
        else c_ntrk_jet->cd(6-ibin);   
                
        h_bkgtrk_data[ibin][0]->Scale(1./h_trk_data[ibin][0]->Integral());
        h_trk_data[ibin][0]->Scale(1./h_trk_data[ibin][0]->Integral());
        h_bkgtrk_data[ibin][1]->Scale(1./h_trk_data[ibin][1]->Integral());
        h_trk_data[ibin][1]->Scale(1./h_trk_data[ibin][1]->Integral());
        h_bkgtrk_data[ibin][2]->Scale(1./h_trk_data[ibin][2]->Integral());
        h_trk_data[ibin][2]->Scale(1./h_trk_data[ibin][2]->Integral());
        h_bkgtrk_data[ibin][3]->Scale(1./h_trk_data[ibin][3]->Integral());
        h_trk_data[ibin][3]->Scale(1./h_trk_data[ibin][3]->Integral());

        h_trk_data[ibin][1]->GetXaxis()->SetTitle("# trks");
        h_trk_data[ibin][1]->GetYaxis()->SetTitle("Counts");
        h_cosm(h_trk_data[ibin][1]);  
        h_trk_data[ibin][1]->GetXaxis()->SetRangeUser(0.,40.);
        h_trk_data[ibin][1]->GetYaxis()->SetRangeUser(0.,0.3);

        h_trk_data[ibin][1]->SetLineColor(kGreen); h_trk_data[ibin][1]->SetMarkerColor(kGreen); h_trk_data[ibin][1]->SetMarkerStyle(20);
        h_trk_data[ibin][1]->Draw("e0 same");
        h_trk_data[ibin][0]->SetLineColor(kBlack); h_trk_data[ibin][0]->SetMarkerColor(kBlack); h_trk_data[ibin][0]->SetMarkerStyle(20);
        h_trk_data[ibin][0]->Draw("e0 same");
        h_trk_data[ibin][2]->SetLineColor(kRed); h_trk_data[ibin][2]->SetMarkerColor(kRed); h_trk_data[ibin][2]->SetMarkerStyle(20);
        h_trk_data[ibin][2]->Draw("e0 same");
        h_trk_data[ibin][3]->SetLineColor(kBlue); h_trk_data[ibin][3]->SetMarkerColor(kBlue); h_trk_data[ibin][3]->SetMarkerStyle(20);
        h_trk_data[ibin][3]->Draw("e0 same");

        TString tmp1 = cent_tag_pp[ibin];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.15,0.85, tmp1);

        //cout<<ibin<<"  "<<h_trk_data[ibin][0]->GetBinContent(1)/h_trk_data[ibin][0]->Integral()<<endl;

        if(ibin==0){
          TLegend *legendy = new TLegend(0.35,0.55,0.98,0.98);
          legendy ->SetLineColor(kWhite);
          legendy ->AddEntry(h_trk_data[ibin][0], "jets (#trks pT>1)", "lep");
          legendy ->AddEntry(h_trk_data[ibin][1], "jets (#trks pT>2)", "lep");
          legendy ->AddEntry(h_trk_data[ibin][2], "jets (#trks pT>4)", "lep");
          legendy ->AddEntry(h_trk_data[ibin][3], "jets (#trks pT>5)", "lep");
          legendy ->Draw("same");
        }
    }
  }

  if(!do_data){
    //fraction of jets with 0 tracks
    TCanvas *c_q_0trk_graph = new TCanvas("c_q_0trk_graph","c_q_0trk_graph",1500,350);
    drawgraphPad();

    pad_titleg->cd();
    title_tex->DrawLatexNDC(0.5,0.1, pythia);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    pad_labelg->cd();
    chg_tex->DrawLatexNDC(0.7,0.4, fraction);

    for(int ibin=0; ibin<nCBins; ibin++){
      if(ibin==0) pada->cd(1);
      else pada->cd(6-ibin);
      h_dummy->Draw("same");
      q_reco_0trk[ibin]->SetMarkerColor(kBlue); q_reco_0trk[ibin]->SetMarkerStyle(5); q_reco_0trk[ibin]->SetMarkerSize(1.2);
      q_reco_0trk[ibin]->Draw("PL");  
      tl_q[ibin]->Draw("same");

      TString tmp1 = cent_tag_pp_MC[ibin];
      tx1 = new TLatex(); tx1->SetTextSize(.1);
      tx1->DrawLatexNDC(0.25,0.85, tmp1);

      if(ibin==0){
        TLegend *legend = new TLegend(0.35,0.15,0.75,0.4);
        legend ->SetLineColor(kWhite);
        legend ->AddEntry(q_reco_0trk[ibin], "quark jet fraction", "lepf");
        legend ->AddEntry((TObject*)0, "(0 trks above pTcut)", "");
        legend ->Draw("same");
      }    
    }

    TCanvas *c_inc_0trk = new TCanvas("c_inc_0trk","c_inc_0trk",1500,350);
    drawgraphPad();

    pad_titleg->cd();
    //title_tex->DrawLatexNDC(0.5,0.1, pythia);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    pad_labelg->cd();
    chg_tex->DrawLatexNDC(0.7,0.4, fraction);

    for(int ibin=0; ibin<nCBins; ibin++){
      if(ibin==0) pada->cd(1);
      else pada->cd(6-ibin);
      h_cosm(h_0trks_data[ibin]);
      h_0trks_data[ibin]->GetXaxis()->SetTitle("track pT cut");
      h_0trks_data[ibin]->GetYaxis()->SetRangeUser(0.,0.12);
      h_0trks_data[ibin]->SetLineColor(kGreen-2); h_0trks_data[ibin]->SetMarkerColor(kGreen-2); h_0trks_data[ibin]->SetMarkerStyle(20);
      h_0trks_data[ibin]->Draw("same");
      h_0trks_MC[ibin]->SetLineColor(kRed-2); h_0trks_MC[ibin]->SetMarkerColor(kRed-2); h_0trks_MC[ibin]->SetMarkerStyle(20);
      h_0trks_MC[ibin]->Draw("same");

      TString tmp1 = cent_tag_pp[ibin];
      tx1 = new TLatex(); tx1->SetTextSize(.1);
      tx1->DrawLatexNDC(0.25,0.85, tmp1);

      if(ibin==0){
        TLegend *legend = new TLegend(0.35,0.15,0.75,0.4);
        legend ->SetLineColor(kWhite);
        legend ->AddEntry(h_0trks_data[ibin], "frac of jets with 0 trks (data)", "lepf");
        legend ->AddEntry(h_0trks_MC[ibin], "frac of jets with 0 trks (MC)", "lepf");
        legend ->Draw("same");
      }    

    }

    TCanvas *c_corr_factor = new TCanvas("c_corr_factor","c_corr_factor",1500,350);
    drawgraphPad();

    pad_titleg->cd();
    title_tex->DrawLatexNDC(0.5,0.1, pythia);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    pad_labelg->cd();
    chg_tex->DrawLatexNDC(0.7,0.4, "quark corr factor");

    for(int ibin=0; ibin<nCBins; ibin++){
      if(ibin==0) pada->cd(1);
      else pada->cd(6-ibin);
      h_cosm(h_corr_factor[ibin]);
      h_corr_factor[ibin]->GetXaxis()->SetTitle("track pT cut");
      h_corr_factor[ibin]->GetYaxis()->SetRangeUser(0.85,1.);
      h_corr_factor[ibin]->SetLineColor(kBlue); h_corr_factor[ibin]->SetMarkerColor(kBlue); h_corr_factor[ibin]->SetMarkerStyle(20);
      h_corr_factor[ibin]->Draw("same");

      TString tmp1 = cent_tag_pp_MC[ibin];
      tx1 = new TLatex(); tx1->SetTextSize(.1);
      tx1->DrawLatexNDC(0.25,0.85, tmp1);
  /*
      if(ibin==0){
        TLegend *legend = new TLegend(0.35,0.15,0.75,0.4);
        legend ->SetLineColor(kWhite);
        legend ->AddEntry(h_0trks_data[ibin], "frac of jets with 0 trks (data)", "lepf");
        legend ->AddEntry(h_0trks_MC[ibin], "frac of jets with 0 trks (MC)", "lepf");
        legend ->Draw("same");
      }    
  */
    }
  }

  if(do_mean_sigma){
    TCanvas *c_chg_avg_MC[ntrkbins];

    for(int ibin2=0;ibin2<ntrkbins;ibin2++){
      c_chg_avg_MC[ibin2] = new TCanvas((TString)("c_chg_avg_MC"+trk[ibin2]),(TString)("c_chg_avg_MC"+trk[ibin2]),1500,600);
      c_chg_avg_MC[ibin2]->Divide(5,2,0);
      //c_chg_avg_MC[ibin2]->Divide(5,1,0);
      gStyle->SetOptStat(0);
      gStyle->SetOptFit(0);
      for(int ibin=0;ibin<nCBins;ibin++){
       for(int ibin3=0;ibin3<njtptbins;ibin3++){  
        c_chg_avg_MC[ibin2]->cd(ibin3+1);
        //if(ibin==0) c_chg_avg_MC[ibin2]->cd(1);   
        //else c_chg_avg_MC[ibin2]->cd(6-ibin);   
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->GetXaxis()->SetTitle("jet chg");
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->GetYaxis()->SetTitle("# of jets");
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->GetYaxis()->SetNdivisions(505);
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->GetYaxis()->SetTitleSize(0.05);
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->GetYaxis()->SetLabelSize(0.05);
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->GetXaxis()->CenterTitle();
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->GetXaxis()->SetNdivisions(505);
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->GetXaxis()->SetTitleSize(0.05);
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->GetXaxis()->SetLabelSize(0.05);
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->GetXaxis()->SetRangeUser(-2.,2.);
        //h_chg_reco_g_avg[ibin][ibin3][ibin2]->GetYaxis()->SetRangeUser(0.,1.);
        h_chg_reco_avg[ibin][ibin3][ibin2]->SetLineColor(kBlack); h_chg_reco_avg[ibin][ibin3][ibin2]->SetMarkerColor(kBlack); h_chg_reco_avg[ibin][ibin3][ibin2]->SetMarkerStyle(2);
        h_chg_reco_q_avg[ibin][ibin3][ibin2]->SetLineColor(kBlue); h_chg_reco_q_avg[ibin][ibin3][ibin2]->SetMarkerColor(kBlue); h_chg_reco_q_avg[ibin][ibin3][ibin2]->SetMarkerStyle(2);    
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->SetLineColor(kRed); h_chg_reco_g_avg[ibin][ibin3][ibin2]->SetMarkerColor(kRed); h_chg_reco_g_avg[ibin][ibin3][ibin2]->SetMarkerStyle(2);
        h_chg_reco_up_avg[ibin][ibin3][ibin2]->SetLineColor(kBlue); h_chg_reco_up_avg[ibin][ibin3][ibin2]->SetMarkerColor(kBlue); h_chg_reco_up_avg[ibin][ibin3][ibin2]->SetMarkerStyle(2);
        h_chg_reco_down_avg[ibin][ibin3][ibin2]->SetLineColor(kGreen-2); h_chg_reco_down_avg[ibin][ibin3][ibin2]->SetMarkerColor(kGreen-2); h_chg_reco_down_avg[ibin][ibin3][ibin2]->SetMarkerStyle(2);

        //h_chg_reco_avg[ibin][ibin3][ibin2]->Draw("e0 same");
        //h_chg_reco_q_avg[ibin][ibin3][ibin2]->Draw("e0 same");
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->Draw("e0 same");
        h_chg_reco_up_avg[ibin][ibin3][ibin2]->Draw("e0 same");
        h_chg_reco_down_avg[ibin][ibin3][ibin2]->Draw("e0 same");
        
        h_chg_reco_avg[ibin][ibin3][ibin2]->Fit(f_chg_reco_avg[ibin][ibin3][ibin2],"Q M R","sames",-2.,2.);
        //h_chg_data_avg[ibin][ibin3][ibin2]->Fit(f_chg_data_avg[ibin][ibin3][ibin2],"Q M R","sames",-2.,2.);
        //h_chg_reco_q_avg[ibin][ibin3][ibin2]->Fit(f_chg_reco_q_avg[ibin][ibin3][ibin2],"Q M R","sames",-2.,2.);
        h_chg_reco_g_avg[ibin][ibin3][ibin2]->Fit(f_chg_reco_g_avg[ibin][ibin3][ibin2],"Q M R","sames",-2.,2.);
        h_chg_reco_up_avg[ibin][ibin3][ibin2]->Fit(f_chg_reco_up_avg[ibin][ibin3][ibin2],"Q M R","sames",-2.,2.);
        h_chg_reco_down_avg[ibin][ibin3][ibin2]->Fit(f_chg_reco_down_avg[ibin][ibin3][ibin2],"Q M R","sames",-2.,2.);
       }
       
          TString tmp1 = cent_tag_pp_MC[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.25,0.85, tmp1);
       
     }
   }

    Double_t inclusive_data_mean[njtptbins], inclusive_mean[njtptbins], quark_mean[njtptbins], gluon_mean[njtptbins], up_mean[njtptbins], down_mean[njtptbins];
    Double_t inclusive_data_res[njtptbins], inclusive_res[njtptbins], quark_res[njtptbins], gluon_res[njtptbins], up_res[njtptbins], down_res[njtptbins];

    TGraph *inclusive_avg[nCBins][ntrkbins]; TGraph *quark_avg[nCBins][ntrkbins]; TGraph *gluon_avg[nCBins][ntrkbins]; TGraph *up_avg[nCBins][ntrkbins]; TGraph *down_avg[nCBins][ntrkbins];
    TGraph *inclusive_sig[nCBins][ntrkbins]; TGraph *quark_sig[nCBins][ntrkbins]; TGraph *gluon_sig[nCBins][ntrkbins]; TGraph *up_sig[nCBins][ntrkbins]; TGraph *down_sig[nCBins][ntrkbins];
    TGraphErrors *gluon_avg_sig[nCBins][ntrkbins];
    TGraph *inclusive_data_avg[nCBins][ntrkbins]; TGraph *inclusive_MC_avg[nCBins][ntrkbins];

    for(int ibin=0; ibin<nCBins; ibin++){
      for(int ibin2=0;ibin2<ntrkbins;ibin2++){
        for(int ibin3=0;ibin3<njtptbins;ibin3++){
          //inclusive_data_mean[ibin3] = f_chg_data_avg[ibin][ibin3][ibin2]->GetParameter(1);
          inclusive_mean[ibin3] = f_chg_reco_avg[ibin][ibin3][ibin2]->GetParameter(1);
          quark_mean[ibin3] = f_chg_reco_q_avg[ibin][ibin3][ibin2]->GetParameter(1);
          gluon_mean[ibin3] = f_chg_reco_g_avg[ibin][ibin3][ibin2]->GetParameter(1);
          up_mean[ibin3] = f_chg_reco_up_avg[ibin][ibin3][ibin2]->GetParameter(1);
          down_mean[ibin3] = f_chg_reco_down_avg[ibin][ibin3][ibin2]->GetParameter(1);

          //inclusive_data_res[ibin3] = f_chg_data_avg[ibin][ibin3][ibin2]->GetParameter(2)/2.;
          inclusive_res[ibin3] = f_chg_reco_avg[ibin][ibin3][ibin2]->GetParameter(2)/2.;
          quark_res[ibin3] = f_chg_reco_q_avg[ibin][ibin3][ibin2]->GetParameter(2)/2.;
          gluon_res[ibin3] = f_chg_reco_g_avg[ibin][ibin3][ibin2]->GetParameter(2)/2.;
          up_res[ibin3] = f_chg_reco_up_avg[ibin][ibin3][ibin2]->GetParameter(2)/2.;
          down_res[ibin3] = f_chg_reco_down_avg[ibin][ibin3][ibin2]->GetParameter(2)/2.;
        }

      inclusive_MC_avg[ibin][ibin2] = new TGraph(njtptbins,jtpt_bound,inclusive_mean);
      inclusive_MC_avg[ibin][ibin2]->SetName((TString)("inclusive_MC_mean_"+cent[ibin]+"_"+trk[ibin2]));
      //inclusive_data_avg[ibin][ibin2] = new TGraph(njtptbins,jtpt_bound,inclusive_data_mean);
      //inclusive_data_avg[ibin][ibin2]->SetName((TString)("inclusive_data_mean_"+cent[ibin]+"_"+trk[ibin2]));

      inclusive_avg[ibin][ibin2] = new TGraph(njtptbins,jtpt_bound,inclusive_mean);
      inclusive_avg[ibin][ibin2]->SetName((TString)("inclusive_mean_"+cent[ibin]+"_"+trk[ibin2]));
      quark_avg[ibin][ibin2] = new TGraph(njtptbins,jtpt_bound,quark_mean);
      quark_avg[ibin][ibin2]->SetName((TString)("quark_mean_"+cent[ibin]+"_"+trk[ibin2]));
      gluon_avg[ibin][ibin2] = new TGraph(njtptbins,jtpt_bound,gluon_mean);
      gluon_avg[ibin][ibin2]->SetName((TString)("gluon_mean_"+cent[ibin]+"_"+trk[ibin2]));
      up_avg[ibin][ibin2] = new TGraph(njtptbins,jtpt_bound,up_mean);
      up_avg[ibin][ibin2]->SetName((TString)("up_mean_"+cent[ibin]+"_"+trk[ibin2]));
      down_avg[ibin][ibin2] = new TGraph(njtptbins,jtpt_bound,down_mean);
      down_avg[ibin][ibin2]->SetName((TString)("down_mean_"+cent[ibin]+"_"+trk[ibin2]));

      inclusive_sig[ibin][ibin2] = new TGraph(njtptbins,jtpt_bound,inclusive_res);
      inclusive_sig[ibin][ibin2]->SetName((TString)("inclusive_res_"+cent[ibin]+"_"+trk[ibin2]));
      quark_sig[ibin][ibin2] = new TGraph(njtptbins,jtpt_bound,quark_res);
      quark_sig[ibin][ibin2]->SetName((TString)("quark_res_"+cent[ibin]+"_"+trk[ibin2]));
      gluon_sig[ibin][ibin2] = new TGraph(njtptbins,jtpt_bound, gluon_res);
      gluon_sig[ibin][ibin2]->SetName((TString)("gluon_res_"+cent[ibin]+"_"+trk[ibin2]));
      up_sig[ibin][ibin2] = new TGraph(njtptbins,jtpt_bound,up_res);
      up_sig[ibin][ibin2]->SetName((TString)("up_res_"+cent[ibin]+"_"+trk[ibin2]));
      down_sig[ibin][ibin2] = new TGraph(njtptbins,jtpt_bound,down_res);
      down_sig[ibin][ibin2]->SetName((TString)("down_res_"+cent[ibin]+"_"+trk[ibin2]));

      gluon_avg_sig[ibin][ibin2] = new TGraphErrors(njtptbins,jtpt_bound, jtpt_bound_err, gluon_mean, gluon_res);
      }
    }

    TCanvas *c_MC_avg_graph = new TCanvas("c_MC_avg_graph","c_MC_avg_graph",1500,500);
    c_MC_avg_graph->Divide(3,1);
    
    for(int ibin=0; ibin<1; ibin++){
     for(int ibin2=0; ibin2<ntrkbins; ibin2++){
      if(ibin2==1) continue;
      if(ibin2==0)c_MC_avg_graph->cd(ibin2+1);
      else c_MC_avg_graph->cd(ibin2);
      //h_dummy3->Draw("same");
      //inclusive_avg[ibin][ibin2]->SetMarkerColor(kBlack); inclusive_avg[ibin][ibin2]->SetMarkerStyle(20); inclusive_avg[ibin][ibin2]->SetMarkerSize(1.3);
      //inclusive_avg[ibin][ibin2]->Draw("PL");
      //quark_avg[ibin][ibin2]->SetMarkerColor(kBlue); quark_avg[ibin][ibin2]->SetMarkerStyle(20); quark_avg[ibin][ibin2]->SetMarkerSize(1.3);
      //quark_avg[ibin][ibin2]->Draw("PL");
      gluon_avg_sig[ibin][ibin2]->SetMarkerColor(kRed-2); gluon_avg_sig[ibin][ibin2]->SetMarkerStyle(kRed-2); gluon_avg_sig[ibin][ibin2]->SetMarkerSize(1.);
      gluon_avg_sig[ibin][ibin2]->SetFillStyle(1001); gluon_avg_sig[ibin][ibin2]->SetFillColorAlpha(kRed-2, 0.3);
      gluon_avg_sig[ibin][ibin2]->Draw("a3");  
      gluon_avg_sig[ibin][ibin2]->SetTitle("");
      gluon_avg_sig[ibin][ibin2]->GetXaxis()->SetRangeUser(100.,220.);
      //gluon_avg_sig[ibin][ibin2]->GetYaxis()->SetRangeUser(-0.7,0.7);
      gluon_avg_sig[ibin][ibin2]->GetXaxis()->SetTitle("jet p_{T} (GeV)");
      gluon_avg_sig[ibin][ibin2]->GetYaxis()->SetTitle("<jet charge>");
      gluon_avg_sig[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
      gluon_avg_sig[ibin][ibin2]->GetYaxis()->SetTitleSize(0.07);
      gluon_avg_sig[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
      gluon_avg_sig[ibin][ibin2]->GetXaxis()->CenterTitle();
      gluon_avg_sig[ibin][ibin2]->GetYaxis()->CenterTitle();
      gluon_avg_sig[ibin][ibin2]->GetYaxis()->SetTitleOffset(0.91);
      gluon_avg_sig[ibin][ibin2]->GetXaxis()->SetTitleOffset(0.93);
      gluon_avg_sig[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
      gluon_avg_sig[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
      gluon_avg_sig[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
      gluon_avg_sig[ibin][ibin2]->Draw("a3");
      gluon_avg[ibin][ibin2]->SetMarkerColor(kRed); gluon_avg[ibin][ibin2]->SetMarkerStyle(20); gluon_avg[ibin][ibin2]->SetMarkerSize(1.3);
      gluon_avg[ibin][ibin2]->Draw("P");
      up_avg[ibin][ibin2]->SetMarkerColor(kBlue); up_avg[ibin][ibin2]->SetMarkerStyle(20); up_avg[ibin][ibin2]->SetMarkerSize(1.3);
      up_avg[ibin][ibin2]->Draw("P");
      down_avg[ibin][ibin2]->SetMarkerColor(kGreen-2); down_avg[ibin][ibin2]->SetMarkerStyle(20); down_avg[ibin][ibin2]->SetMarkerSize(1.3);
      down_avg[ibin][ibin2]->Draw("P");
        
        TString tmp = trk_tag[ibin2];
        tx = new TLatex(); tx->SetTextSize(.08);
        tx->DrawLatexNDC(0.15,0.8, tmp);

        if(ibin2==0){
            TString tmp1 = cent_tag_pp_MC[ibin];
            tx1 = new TLatex(); tx1->SetTextSize(.1);
            tx1->DrawLatexNDC(0.35,0.9, tmp1);

            TLegend *legendy = new TLegend(0.35,0.55,0.98,0.98);
            legendy ->SetLineColor(kWhite);
            legendy ->AddEntry(gluon_avg[ibin][3], "gluon mean", "lep");
            legendy ->AddEntry(up_avg[ibin][3], "up mean", "lep");
            legendy ->AddEntry(down_avg[ibin][3], "down mean", "lep");
            legendy ->AddEntry(gluon_avg_sig[ibin][3], "gluon template width", "f");
            legendy ->Draw("same");
        }
     }
    }

    TCanvas *c_MC_res_graph = new TCanvas("c_MC_res_graph","c_MC_res_graph",1500,500);
    c_MC_res_graph->Divide(3,1,0);
    
    for(int ibin=0; ibin<1; ibin++){
     for(int ibin2=0; ibin2<ntrkbins; ibin2++){
      if(ibin2==1) continue;
      if(ibin2==0)c_MC_res_graph->cd(ibin2+1);
      else c_MC_res_graph->cd(ibin2);
      //c_MC_res_graph->cd(ibin2+1);
      h_dummy4->Draw("same");
      inclusive_sig[ibin][ibin2]->SetMarkerColor(kBlack); inclusive_sig[ibin][ibin2]->SetMarkerStyle(20); inclusive_sig[ibin][ibin2]->SetMarkerSize(1.3);
      //inclusive_sig[ibin][ibin2]->Draw("PL");
      quark_sig[ibin][ibin2]->SetMarkerColor(kBlue); quark_sig[ibin][ibin2]->SetMarkerStyle(20); quark_sig[ibin][ibin2]->SetMarkerSize(1.3);
      //quark_sig[ibin][ibin2]->Draw("PL");
      gluon_sig[ibin][ibin2]->SetMarkerColor(kRed); gluon_sig[ibin][ibin2]->SetMarkerStyle(20); gluon_sig[ibin][ibin2]->SetMarkerSize(1.3);
      gluon_sig[ibin][ibin2]->Draw("PL");
      up_sig[ibin][ibin2]->SetMarkerColor(kBlue); up_sig[ibin][ibin2]->SetMarkerStyle(20); up_sig[ibin][ibin2]->SetMarkerSize(1.3);
      up_sig[ibin][ibin2]->Draw("PL");
      down_sig[ibin][ibin2]->SetMarkerColor(kGreen-2); down_sig[ibin][ibin2]->SetMarkerStyle(20); down_sig[ibin][ibin2]->SetMarkerSize(1.3);
      down_sig[ibin][ibin2]->Draw("PL");
        
        TString tmp = trk_tag[ibin2];
        tx = new TLatex(); tx->SetTextSize(.1);
        tx->DrawLatexNDC(0.15,0.85, tmp);
     }
    }

    TCanvas *c_MC_data_avg_graph = new TCanvas("c_MC_data_avg_graph","c_MC_data_avg_graph",1500,300);
    c_MC_data_avg_graph->Divide(5,1,0);
    
    for(int ibin=0; ibin<nCBins; ibin++){
       for(int ibin2=0; ibin2<1; ibin2++){
          if(ibin==0)c_MC_data_avg_graph->cd(1);
          else c_MC_data_avg_graph->cd(6-ibin);
          h_chg_data_avg[ibin][ibin2]->SetTitle("");
          h_chg_data_avg[ibin][ibin2]->Rebin(10);
          h_chg_data_avg[ibin][ibin2]->GetXaxis()->SetRangeUser(120.,300.);
          h_chg_data_avg[ibin][ibin2]->GetYaxis()->SetRangeUser(-0.02,0.07);
          h_chg_data_avg[ibin][ibin2]->GetXaxis()->SetTitle("jet p_{T} (GeV)");
          h_chg_data_avg[ibin][ibin2]->GetYaxis()->SetTitle("<jet charge>");
          h_chg_data_avg[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
          h_chg_data_avg[ibin][ibin2]->GetYaxis()->SetTitleSize(0.07);
          h_chg_data_avg[ibin][ibin2]->GetYaxis()->SetLabelSize(0.05);
          h_chg_data_avg[ibin][ibin2]->GetXaxis()->CenterTitle();
          h_chg_data_avg[ibin][ibin2]->GetYaxis()->CenterTitle();
          h_chg_data_avg[ibin][ibin2]->GetYaxis()->SetTitleOffset(0.91);
          h_chg_data_avg[ibin][ibin2]->GetXaxis()->SetTitleOffset(0.93);
          h_chg_data_avg[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
          h_chg_data_avg[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
          h_chg_data_avg[ibin][ibin2]->GetXaxis()->SetLabelSize(0.05);
          inclusive_MC_avg[ibin][ibin2]->SetMarkerColor(kBlack); inclusive_MC_avg[ibin][ibin2]->SetMarkerStyle(20); inclusive_MC_avg[ibin][ibin2]->SetMarkerSize(1.3);
          if(ibin==0)
          {  h_chg_data_avg[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_data_avg[ibin][ibin2]->SetMarkerStyle(4); h_chg_data_avg[ibin][ibin2]->SetMarkerSize(1.3);}
          else if(ibin==4)
          {  h_chg_data_avg[ibin][ibin2]->SetMarkerColor(kGreen-2); h_chg_data_avg[ibin][ibin2]->SetMarkerStyle(20); h_chg_data_avg[ibin][ibin2]->SetMarkerSize(1.3);}
          else 
          {  h_chg_data_avg[ibin][ibin2]->SetMarkerColor(kGreen-2); h_chg_data_avg[ibin][ibin2]->SetMarkerStyle(20); h_chg_data_avg[ibin][ibin2]->SetMarkerSize(1.3);}
          h_chg_data_avg[ibin][ibin2]->Draw("e0"); //inclusive_MC_avg[ibin][ibin2]->Draw("P");   
          //h_chg_data_avg[0][ibin2]->Draw("e0 same");
          //h_chg_data_avg[4][ibin2]->Draw("e0 same");
          drawline(120.,0.,300.,0.);

          TString tmp1 = cent_tag_pp[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.35,0.9, tmp1);

          if(ibin==0){
            TLegend *legendy = new TLegend(0.35,0.55,0.98,0.98);
            legendy ->SetLineColor(kWhite);
            //legendy ->AddEntry(h_chg_data_avg[0][ibin2], "pp ref", "lep");
            //legendy ->AddEntry(h_chg_data_avg[4][ibin2], "50-100% PbPb", "lep");
            //legendy ->AddEntry(h_chg_data_avg[ibin][ibin2], "data", "lep");
            //legendy ->Draw("same");
          }
       }
    }
  }
}