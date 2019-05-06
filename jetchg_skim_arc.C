#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TRandom.h"
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
#include <TCut.h>
#include "TDatime.h"
#include <vector>
#include "TCanvas.h"
//#include "jffcorr/nCScorr.h"
#include "TrkCorr_July22_Iterative_pp_eta2p4/getTrkCorr.h"
#include "xiaoTrkCorr/xiaoTrkCorr.h"

bool ispp;

const int nCBins = (ispp==true) ? 1 : 4;
const int nptBins = 55;
const int ntrkBins = 3;
const int nkbins = 3;
const int netaBins = 15;
const int trkeffbins = 11;

using namespace std;

int mypbin, mycbin, myptbin, myrefptbin, myetabin;

char saythis[500];

TDatime* date = new TDatime();

TString cent[5] = {"0","1","2","3","4"};
TString trk[4] = {"0p7","2","4","5"};

TString pt[56] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55"};

int jt_nbins = 55;
Double_t jt_bin_bounds[56] = {50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280.,290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.,510.,520.,530.,540.,550.,560.,570.,580.,590.,600.};

Double_t eta_bounds[16] = {-1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5};

///pp MC Pythia6 pthat weights
double pthatbins[9] = {50,80,100,120,170,220,280,370,9999};
double xsecs[9] = {3.778E-03, 4.412E-04, 1.511E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
double pthatEntries[8] ={766451, 292063, 91200, 468748, 447938, 259208, 234447, 50972};

///pp MC Pythia8 pthat weights
double pthatbins_Pythia8[9] = {50,80,100,120,170,220,280,370,9999};
double xsecs_Pythia8[9] = {3.778E-03, 4.412E-04, 1.511E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
double pthatEntries_Pythia8[8] ={175696, 145099, 160445, 258400, 189793, 196579, 54724, 12064};

///PbPb MC Pythia6+Hydjet pthat weights
double pthatbins_PH[9] = {15,30,50,80,120,170,220,280,9999};
double xsecs_PH[9] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 0};
double pthatEntries_PH[8] ={0, 0, 0, 2571566, 2.85082e+06, 2.68057e+06, 2.89138e+06, 950344};

float CBins[5] = {-1., 20., 60., 100., 200.};
float kappa[4] = {0.3, 0.5, 0.7};

/*
//ptcut
double trk_eff_scale_pt2[5] = {0.99,1.,0.98,1.,1.};
double trk_eff_scale_pt4[5] = {0.99,0.99,1.,1.,1.};
double trk_eff_scale_pt5[5] = {0.99,0.99,1.,1.,1.};

//kappa
double trk_eff_scale_k0p3[5] = {0.99,1.,0.98,1.,1.};
double trk_eff_scale_k0p5[5] = {0.99,1.,0.98,1.,1.};
double trk_eff_scale_k0p7[5] = {0.975,0.975,0.98,0.97,1.};
*/
 
int n_trk[ntrkBins], n_gen[ntrkBins], n_bkggen[ntrkBins], n_bkgtrk[ntrkBins], n_eff_gen[ntrkBins], n_sube0_gen[ntrkBins], n_eff_sube0_gen[ntrkBins], n_sube0_gengen[ntrkBins];
double trk_chg_ptcut_trkeff[ntrkBins][trkeffbins], trk_bkg_ptcut_trkeff[ntrkBins][trkeffbins], trk_chg_kappa_trkeff[ntrkBins][trkeffbins], trk_bkg_kappa_trkeff[ntrkBins][trkeffbins];
double trk_chg_ptcut[ntrkBins], trk_bkg_ptcut[ntrkBins], trk_chg_kappa[ntrkBins], trk_bkg_kappa[ntrkBins];
double gen_chg_ptcut[ntrkBins], gen_bkg_ptcut[ntrkBins], gen_chg_ptcut_eff[ntrkBins], gen_bkg_ptcut_eff[ntrkBins];
double gen_chg_ptcut_sube0[ntrkBins], gen_bkg_ptcut_sube0[ntrkBins], gen_chg_ptcut_sube0_bkgadd[ntrkBins], gen_chg_ptcut_sube0_eff[ntrkBins], gengen_chg_ptcut_sube0[ntrkBins];
double gen_chg_kappa[ntrkBins], gen_chg_kappa_sube0[ntrkBins], gen_chg_kappa_sube0_eff[ntrkBins], gen_chg_kappa_eff[ntrkBins], gengen_chg_kappa_sube0[ntrkBins];
double trk_correction, trk_eff, jet_corrpt;

float vz, pthat, pthat_weight, weight_vz, weight_cen; 
int hiBin;
vector<float> *calo_jtpt=0, *calo_corrpt=0, *calo_jteta=0, *calo_jtphi=0, *calo_refpt=0, *calo_refeta=0, *calo_refphi=0; 
vector<int> *calo_refparton_flavor=0, *trkChg=0, *genChg=0, *gensube=0, *genpdg=0;
vector<float> *trkPt=0, *trkEta=0, *trkPhi=0, *genPt=0, *genEta=0, *genPhi=0;

//cymbal tune centrality reweighting
double fcent_cymbal(double centrality, TF1* fcent1){ 
  return (centrality < 194) ? fcent1->Eval(centrality) : 1;
}

TF1 *f_res_gauss = new TF1("f_res_gauss", "exp(-((x-1.)^2)/(2*([0]^2)))", 0., 2.);

void jetchg_skim_arc(bool ispp=1, bool isdata=1, bool do_eff_studies=0){

  //nCScorr *corrpt = new nCScorr(ispp,true);

  TrkCorr* pptrkCorr = new TrkCorr("TrkCorr_July22_Iterative_pp_eta2p4/");
  xiaoTrkCorr* PbPbtrkCorr = new xiaoTrkCorr("xiaoTrkCorr/eta_symmetry_cymbalCorr_FineBin.root");

  f_res_gauss->FixParameter(0,0.05);

  gRandom->SetSeed(0);

//// defining histos and profiles

  TH1D *h_hibin = new TH1D("h_hibin","",200,0.,200.);
  h_hibin->Sumw2();
  TH1D *h_hibin_nrw = new TH1D("h_hibin_nrw","",200,0.,200.);
  h_hibin_nrw->Sumw2();
  TH2F *h_hibin_cbin = new TH2F("h_hibin_cbin","",200,0.,200.,7,-1.,6.);
  h_hibin_cbin->Sumw2();

  TH1D *h_vz[nCBins];
  TH1D *h_vz_nrw[nCBins];
  TH1D *h_pthat[nCBins];

  TH1D *h_ntrk[nCBins];
  TH1D *h_nbkgtrk[nCBins];

  TH1D *h_reco_full[nCBins];
  TH1D *h_reco_full_q[nCBins];
  TH1D *h_reco_full_g[nCBins];

  TH1D *h_reco_corr[nCBins];
  TH1D *h_reco_corr_q[nCBins];
  TH1D *h_reco_corr_g[nCBins];
  TH1D *h_reco_corr_u[nCBins];

  TH1D *h_reco_corr_no0trks[nCBins][ntrkBins];
  TH1D *h_reco_corr_no0trks_q[nCBins][ntrkBins];
  TH1D *h_reco_corr_no0trks_g[nCBins][ntrkBins];

  TH1D *h_reco_corr_0trks[nCBins][ntrkBins];
  TH1D *h_reco_corr_0trks_q[nCBins][ntrkBins];
  TH1D *h_reco_corr_0trks_up[nCBins][ntrkBins];
  TH1D *h_reco_corr_0trks_down[nCBins][ntrkBins];
  TH1D *h_reco_corr_0trks_others[nCBins][ntrkBins];
  TH1D *h_reco_corr_0trks_g[nCBins][ntrkBins];

  TH1D *h_gen_full_no0trks[nCBins][ntrkBins];
  TH1D *h_gen_full_no0trks_q[nCBins][ntrkBins];
  TH1D *h_gen_full_no0trks_g[nCBins][ntrkBins];

  TH1D *h_gen_full_0trks[nCBins][ntrkBins];
  TH1D *h_gen_full_0trks_q[nCBins][ntrkBins];
  TH1D *h_gen_full_0trks_up[nCBins][ntrkBins];
  TH1D *h_gen_full_0trks_down[nCBins][ntrkBins];
  TH1D *h_gen_full_0trks_others[nCBins][ntrkBins];
  TH1D *h_gen_full_0trks_g[nCBins][ntrkBins];

  TH1D *h_reco_corr_up[nCBins];
  TH1D *h_reco_corr_upbar[nCBins];
  TH1D *h_reco_corr_d[nCBins];
  TH1D *h_reco_corr_dbar[nCBins];
  TH1D *h_reco_corr_c[nCBins];
  TH1D *h_reco_corr_s[nCBins];
  TH1D *h_reco_corr_b[nCBins];

  TH1D *h_gen_full[nCBins];
  TH1D *h_gen_full_q[nCBins];
  TH1D *h_gen_full_g[nCBins];  
  TH1D *h_gen_full_u[nCBins];

  TH1D *h_gen_full_c[nCBins];
  TH1D *h_gen_full_s[nCBins];
  TH1D *h_gen_full_b[nCBins];

  TH1D *h_gen_full_up[nCBins];
  TH1D *h_gen_full_upbar[nCBins];
  TH1D *h_gen_full_d[nCBins];
  TH1D *h_gen_full_dbar[nCBins];

  TH2F *h_chg_refpt[nCBins][ntrkBins];
  TH2F *h_chg_refpt_nobkgsub[nCBins][ntrkBins];
  TH2F *h_chg_refpt_eff[nCBins][ntrkBins];
  TH2F *h_chg_refpt_sube0[nCBins][ntrkBins];
  TH2F *h_chg_refpt_sube0_bkgadd[nCBins][ntrkBins];
  TH2F *h_chg_refpt_q[nCBins][ntrkBins];  
  TH2F *h_chg_refpt_g[nCBins][ntrkBins];  
  TH2F *h_chg_refpt_u[nCBins][ntrkBins];  

  TH2F *h_chg_corrpt[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_q[nCBins][ntrkBins];  
  TH2F *h_chg_corrpt_g[nCBins][ntrkBins];  
  TH2F *h_chg_corrpt_u[nCBins][ntrkBins];

  TH2F *h_chg_corrpt_nobkgsub[nCBins][ntrkBins];

  TH2F *h_chg_corrpt_trkeff[nCBins][ntrkBins][trkeffbins];
  TH2F *h_chg_corrpt_kappa_trkeff[nCBins][ntrkBins][trkeffbins];

  TH2F *h_chg_corrpt_kappa[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_q_kappa[nCBins][ntrkBins];  
  TH2F *h_chg_corrpt_g_kappa[nCBins][ntrkBins];  
  TH2F *h_chg_corrpt_u_kappa[nCBins][ntrkBins];

  TH2F *h_chg_refpt_bkg[nCBins][ntrkBins];
  TH2F *h_chg_refpt_bkg_q[nCBins][ntrkBins];
  TH2F *h_chg_refpt_bkg_g[nCBins][ntrkBins];

  TH2F *h_chg_refpt_subenon0[nCBins];

  TH2F *h_chg_corrpt_bkg[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_bkg_q[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_bkg_g[nCBins][ntrkBins];

  TH2F *h_chg_corrpt_bkg_kappa[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_bkg_q_kappa[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_bkg_g_kappa[nCBins][ntrkBins];

  TH2F *h_chg_corrpt_mix[nCBins][ntrkBins];

  TH2F *h_chg_refpt_up[nCBins][ntrkBins];  
  TH2F *h_chg_refpt_upbar[nCBins][ntrkBins];
  TH2F *h_chg_refpt_d[nCBins][ntrkBins];
  TH2F *h_chg_refpt_dbar[nCBins][ntrkBins];

  TH2F *h_chg_refpt_c[nCBins][ntrkBins];  
  TH2F *h_chg_refpt_s[nCBins][ntrkBins];  
  TH2F *h_chg_refpt_b[nCBins][ntrkBins];  

  TH2F *h_chg_corrpt_up[nCBins][ntrkBins];  
  TH2F *h_chg_corrpt_upbar[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_d[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_dbar[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_c[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_s[nCBins][ntrkBins]; 
  TH2F *h_chg_corrpt_b[nCBins][ntrkBins];

  TH2F *h_chg_sig_bkg[nCBins][ntrkBins];  
  TH2F *h_chg_sig_ntrk[nCBins][ntrkBins];
  TH2F *h_chg_bkg_nbkg[nCBins][ntrkBins];
  TH2F *h_chg_sig_ntrk_pos[nCBins][ntrkBins];
  TH2F *h_chg_bkg_nbkg_pos[nCBins][ntrkBins];
  TH2F *h_chg_sig_ntrk_neg[nCBins][ntrkBins];
  TH2F *h_chg_bkg_nbkg_neg[nCBins][ntrkBins];

  TH2F *h_chg_corrpt_up_kappa[nCBins][ntrkBins];  
  TH2F *h_chg_corrpt_upbar_kappa[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_d_kappa[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_dbar_kappa[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_c_kappa[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_s_kappa[nCBins][ntrkBins]; 
  TH2F *h_chg_corrpt_b_kappa[nCBins][ntrkBins];

  TH2F *h_chg_corrpt_highest[nCBins];
  TH2F *h_chg_corrpt_highest_2[nCBins];

  TH2F *h_trk_refpt[nCBins];
  TH2F *h_trk_refpt_q[nCBins];  
  TH2F *h_trk_refpt_g[nCBins];  
  TH2F *h_trk_refpt_u[nCBins];

  TH2F *h_bkgtrk_refpt[nCBins];
  TH2F *h_sube0trk_refpt[nCBins];

  TH2F *h_trk_corrpt[nCBins][ntrkBins];
  TH2F *h_bkgtrk_corrpt[nCBins][ntrkBins];
  TH2F *h_trk_corrpt_q[nCBins];  
  TH2F *h_trk_corrpt_g[nCBins];  
  TH2F *h_trk_corrpt_u[nCBins];

  TH1D *h_eta_full[nCBins];
  TH1D *h_eta_full_q[nCBins];
  TH1D *h_eta_full_g[nCBins];
  TH1D *h_eta_full_u[nCBins];

  TH1D *h_eta_full_up[nCBins];
  TH1D *h_eta_full_upbar[nCBins];
  TH1D *h_eta_full_d[nCBins];
  TH1D *h_eta_full_dbar[nCBins];

  TH2F *h_jt_gen_cloure[nCBins];  
  TH2F *h_jt_gen_q_cloure[nCBins];  
  TH2F *h_jt_gen_g_cloure[nCBins];  

  TH1D *h_sampled_bkg[nCBins];

  TH2F *h_chg_recopt_genpt[nCBins][ntrkBins];
  TH2F *h_chg_recoptsube0_genptsube0[nCBins][ntrkBins];
  TH2F *h_chg_recopt_recoptsube0[nCBins][ntrkBins];
  TH2F *h_chg_recopt_genptsube0[nCBins][ntrkBins];
  TH2F *h_chg_effgenpt_genptsube0[nCBins][ntrkBins];
  TH2F *h_chg_recorecopt_gengensube0pt[nCBins][ntrkBins];

  TH2F *h_chg_recoptsube0_genptsube0_up[nCBins][ntrkBins];
  TH2F *h_chg_recopt_recoptsube0_up[nCBins][ntrkBins];
  TH2F *h_chg_recopt_genptsube0_up[nCBins][ntrkBins];
  TH2F *h_chg_effgenpt_genptsube0_up[nCBins][ntrkBins];
  TH2F *h_chg_recorecopt_gengensube0pt_up[nCBins][ntrkBins];
  TH2F *h_chg_recoptsube0_genptsube0_d[nCBins][ntrkBins];
  TH2F *h_chg_recopt_recoptsube0_d[nCBins][ntrkBins];
  TH2F *h_chg_recopt_genptsube0_d[nCBins][ntrkBins];
  TH2F *h_chg_effgenpt_genptsube0_d[nCBins][ntrkBins];
  TH2F *h_chg_recorecopt_gengensube0pt_d[nCBins][ntrkBins];
  TH2F *h_chg_recoptsube0_genptsube0_g_oth[nCBins][ntrkBins];
  TH2F *h_chg_recopt_recoptsube0_g_oth[nCBins][ntrkBins];
  TH2F *h_chg_recopt_genptsube0_g_oth[nCBins][ntrkBins];
  TH2F *h_chg_effgenpt_genptsube0_g_oth[nCBins][ntrkBins];
  TH2F *h_chg_recorecopt_gengensube0pt_g_oth[nCBins][ntrkBins];

  TH1D *h_chg_recopt[nCBins][ntrkBins];
  TH1D *h_chg_recopt_datasample[nCBins][ntrkBins];
  TH1D *h_chg_genpt[nCBins][ntrkBins];
  TH1D *h_chg_eff_genpt[nCBins][ntrkBins];
  TH1D *h_chg_recopt_sube0[nCBins][ntrkBins];
  TH1D *h_chg_genpt_sube0[nCBins][ntrkBins];
  TH1D *h_chg_gengensube0pt[nCBins][ntrkBins];

  TH1D *h_chg_recopt_bkg[nCBins][ntrkBins];
  TH1D *h_chg_genpt_bkg[nCBins][ntrkBins];

  TH1D *h_chg_recopt_up[nCBins][ntrkBins];
  TH1D *h_chg_recopt_up_datasample[nCBins][ntrkBins];
  TH1D *h_chg_genpt_up[nCBins][ntrkBins];
  TH1D *h_chg_eff_genpt_up[nCBins][ntrkBins];
  TH1D *h_chg_recopt_sube0_up[nCBins][ntrkBins];
  TH1D *h_chg_genpt_sube0_up[nCBins][ntrkBins];
  TH1D *h_chg_gengensube0pt_up[nCBins][ntrkBins];

  TH1D *h_chg_recopt_d[nCBins][ntrkBins];
  TH1D *h_chg_recopt_d_datasample[nCBins][ntrkBins];
  TH1D *h_chg_genpt_d[nCBins][ntrkBins];
  TH1D *h_chg_eff_genpt_d[nCBins][ntrkBins];
  TH1D *h_chg_recopt_sube0_d[nCBins][ntrkBins];
  TH1D *h_chg_genpt_sube0_d[nCBins][ntrkBins];
  TH1D *h_chg_gengensube0pt_d[nCBins][ntrkBins];

  TH1D *h_chg_recopt_g_oth[nCBins][ntrkBins];
  TH1D *h_chg_recopt_g_oth_datasample[nCBins][ntrkBins];
  TH1D *h_chg_genpt_g_oth[nCBins][ntrkBins];
  TH1D *h_chg_eff_genpt_g_oth[nCBins][ntrkBins];
  TH1D *h_chg_recopt_sube0_g_oth[nCBins][ntrkBins];
  TH1D *h_chg_genpt_sube0_g_oth[nCBins][ntrkBins];
  TH1D *h_chg_gengensube0pt_g_oth[nCBins][ntrkBins];

//kappa
  TH2F *h_chg_recopt_genpt_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recoptsube0_genptsube0_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recopt_recoptsube0_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recopt_genptsube0_kappa[nCBins][ntrkBins];
  TH2F *h_chg_effgenpt_genptsube0_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recorecopt_gengensube0pt_kappa[nCBins][ntrkBins];

  TH2F *h_chg_recoptsube0_genptsube0_up_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recopt_recoptsube0_up_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recopt_genptsube0_up_kappa[nCBins][ntrkBins];
  TH2F *h_chg_effgenpt_genptsube0_up_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recorecopt_gengensube0pt_up_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recoptsube0_genptsube0_d_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recopt_recoptsube0_d_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recopt_genptsube0_d_kappa[nCBins][ntrkBins];
  TH2F *h_chg_effgenpt_genptsube0_d_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recorecopt_gengensube0pt_d_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recoptsube0_genptsube0_g_oth_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recopt_recoptsube0_g_oth_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recopt_genptsube0_g_oth_kappa[nCBins][ntrkBins];
  TH2F *h_chg_effgenpt_genptsube0_g_oth_kappa[nCBins][ntrkBins];
  TH2F *h_chg_recorecopt_gengensube0pt_g_oth_kappa[nCBins][ntrkBins];

  TH1D *h_chg_recopt_kappa[nCBins][ntrkBins];
  TH1D *h_chg_recopt_kappa_datasample[nCBins][ntrkBins];
  TH1D *h_chg_genpt_kappa[nCBins][ntrkBins];
  TH1D *h_chg_eff_genpt_kappa[nCBins][ntrkBins];
  TH1D *h_chg_recopt_sube0_kappa[nCBins][ntrkBins];
  TH1D *h_chg_genpt_sube0_kappa[nCBins][ntrkBins];
  TH1D *h_chg_gengensube0pt_kappa[nCBins][ntrkBins];

  TH1D *h_chg_recopt_up_kappa[nCBins][ntrkBins];
  TH1D *h_chg_recopt_up_kappa_datasample[nCBins][ntrkBins];
  TH1D *h_chg_genpt_up_kappa[nCBins][ntrkBins];
  TH1D *h_chg_eff_genpt_up_kappa[nCBins][ntrkBins];
  TH1D *h_chg_recopt_sube0_up_kappa[nCBins][ntrkBins];
  TH1D *h_chg_genpt_sube0_up_kappa[nCBins][ntrkBins];
  TH1D *h_chg_gengensube0pt_up_kappa[nCBins][ntrkBins];

  TH1D *h_chg_recopt_d_kappa[nCBins][ntrkBins];
  TH1D *h_chg_recopt_d_kappa_datasample[nCBins][ntrkBins];
  TH1D *h_chg_genpt_d_kappa[nCBins][ntrkBins];
  TH1D *h_chg_eff_genpt_d_kappa[nCBins][ntrkBins];
  TH1D *h_chg_recopt_sube0_d_kappa[nCBins][ntrkBins];
  TH1D *h_chg_genpt_sube0_d_kappa[nCBins][ntrkBins];
  TH1D *h_chg_gengensube0pt_d_kappa[nCBins][ntrkBins];

  TH1D *h_chg_recopt_g_oth_kappa[nCBins][ntrkBins];
  TH1D *h_chg_recopt_g_oth_kappa_datasample[nCBins][ntrkBins];
  TH1D *h_chg_genpt_g_oth_kappa[nCBins][ntrkBins];
  TH1D *h_chg_eff_genpt_g_oth_kappa[nCBins][ntrkBins];
  TH1D *h_chg_recopt_sube0_g_oth_kappa[nCBins][ntrkBins];
  TH1D *h_chg_genpt_sube0_g_oth_kappa[nCBins][ntrkBins];
  TH1D *h_chg_gengensube0pt_g_oth_kappa[nCBins][ntrkBins];

  TH1D *h_chg_genpt_sube0_bkgadd[nCBins][ntrkBins];
  TH1D *h_chg_recopt_sube0_bkgadd[nCBins][ntrkBins];
  TH1D *h_chg_genpt_sube0_bkgadd_kappa[nCBins][ntrkBins];
  TH1D *h_chg_recopt_sube0_bkgadd_kappa[nCBins][ntrkBins];

  TH1D *h_trkpt_reco[nCBins]; 
  TH1D *h_trkpt_pos_reco[nCBins]; 
  TH1D *h_trkpt_neg_reco[nCBins]; 

  TH1D *h_trkpt_gen[nCBins]; 
  TH1D *h_trkpt_pos_gen[nCBins]; 
  TH1D *h_trkpt_neg_gen[nCBins]; 

  for (int ibin=0;ibin<nCBins;ibin++){

    sprintf(saythis,"h_ntrk_cent%d",ibin);
    h_ntrk[ibin] = new TH1D(saythis, "", 100,0,100);
    h_ntrk[ibin]->Sumw2();

    sprintf(saythis,"h_nbkgtrk_cent%d",ibin);
    h_nbkgtrk[ibin] = new TH1D(saythis, "", 100,0,100);
    h_nbkgtrk[ibin]->Sumw2();

    sprintf(saythis,"h_eta_full_cent%d",ibin);
    h_eta_full[ibin] = new TH1D(saythis, "", netaBins-1,eta_bounds);
    h_eta_full[ibin]->Sumw2();

    sprintf(saythis,"h_eta_full_q_cent%d",ibin);
    h_eta_full_q[ibin] = new TH1D(saythis,"",netaBins-1,eta_bounds);
    h_eta_full_q[ibin]->Sumw2();

    sprintf(saythis,"h_eta_full_g_cent%d",ibin);
    h_eta_full_g[ibin] = new TH1D(saythis,"",netaBins-1,eta_bounds);
    h_eta_full_g[ibin]->Sumw2();

    sprintf(saythis,"h_eta_full_u_cent%d",ibin);
    h_eta_full_u[ibin] = new TH1D(saythis,"",netaBins-1,eta_bounds);
    h_eta_full_u[ibin]->Sumw2();

    sprintf(saythis,"h_eta_full_up_cent%d",ibin);
    h_eta_full_up[ibin] = new TH1D(saythis,"",netaBins-1,eta_bounds);
    h_eta_full_up[ibin]->Sumw2();

    sprintf(saythis,"h_eta_full_upbar_cent%d",ibin);
    h_eta_full_upbar[ibin] = new TH1D(saythis,"",netaBins-1,eta_bounds);
    h_eta_full_upbar[ibin]->Sumw2();

    sprintf(saythis,"h_eta_full_d_cent%d",ibin);
    h_eta_full_d[ibin] = new TH1D(saythis,"",netaBins-1,eta_bounds);
    h_eta_full_d[ibin]->Sumw2();

    sprintf(saythis,"h_eta_full_dbar_cent%d",ibin);
    h_eta_full_dbar[ibin] = new TH1D(saythis,"",netaBins-1,eta_bounds);
    h_eta_full_dbar[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_cent%d",ibin);
    h_gen_full[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_gen_full[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_q_cent%d",ibin);
    h_gen_full_q[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_gen_full_q[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_g_cent%d",ibin);
    h_gen_full_g[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_gen_full_g[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_u_cent%d",ibin);
    h_gen_full_u[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_gen_full_u[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_up_cent%d",ibin);
    h_gen_full_up[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_gen_full_up[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_upbar_cent%d",ibin);
    h_gen_full_upbar[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_gen_full_upbar[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_d_cent%d",ibin);
    h_gen_full_d[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_gen_full_d[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_dbar_cent%d",ibin);
    h_gen_full_dbar[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_gen_full_dbar[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_c_cent%d",ibin);
    h_gen_full_c[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_gen_full_c[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_s_cent%d",ibin);
    h_gen_full_s[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_gen_full_s[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_b_cent%d",ibin);
    h_gen_full_b[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_gen_full_b[ibin]->Sumw2();

    sprintf(saythis,"h_pthat_cent%d",ibin);
    h_pthat[ibin] = new TH1D(saythis, "",50,0.,500.);
    h_pthat[ibin]->Sumw2();

    sprintf(saythis,"h_vz_cent%d",ibin);
    h_vz[ibin] = new TH1D(saythis, "", 30,-15.,15.);
    h_vz[ibin]->Sumw2();

    sprintf(saythis,"h_vz_nrw_cent%d",ibin);
    h_vz_nrw[ibin] = new TH1D(saythis, "", 30,-15.,15.);
    h_vz_nrw[ibin]->Sumw2();

    sprintf(saythis,"h_reco_full_cent%d",ibin);
    h_reco_full[ibin] = new TH1D(saythis, "", jt_nbins,jt_bin_bounds);
    h_reco_full[ibin]->Sumw2();

    sprintf(saythis,"h_reco_full_q_cent%d",ibin);
    h_reco_full_q[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_reco_full_q[ibin]->Sumw2();

    sprintf(saythis,"h_reco_full_g_cent%d",ibin);
    h_reco_full_g[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_reco_full_g[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_cent%d",ibin);
    h_reco_corr[ibin] = new TH1D(saythis, "", jt_nbins,jt_bin_bounds);
    h_reco_corr[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_q_cent%d",ibin);
    h_reco_corr_q[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_reco_corr_q[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_g_cent%d",ibin);
    h_reco_corr_g[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_reco_corr_g[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_u_cent%d",ibin);
    h_reco_corr_u[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_reco_corr_u[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_up_cent%d",ibin);
    h_reco_corr_up[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_reco_corr_up[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_upbar_cent%d",ibin);
    h_reco_corr_upbar[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_reco_corr_upbar[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_d_cent%d",ibin);
    h_reco_corr_d[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_reco_corr_d[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_dbar_cent%d",ibin);
    h_reco_corr_dbar[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_reco_corr_dbar[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_c_cent%d",ibin);
    h_reco_corr_c[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_reco_corr_c[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_s_cent%d",ibin);
    h_reco_corr_s[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_reco_corr_s[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_b_cent%d",ibin);
    h_reco_corr_b[ibin] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
    h_reco_corr_b[ibin]->Sumw2();

    sprintf(saythis,"h_trk_refpt_cent%d",ibin);
    h_trk_refpt[ibin] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,40,0.,40.);
    h_trk_refpt[ibin]->Sumw2();

    sprintf(saythis,"h_trk_refpt_q_cent%d",ibin);
    h_trk_refpt_q[ibin] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,40,0.,40.);
    h_trk_refpt_q[ibin]->Sumw2();

    sprintf(saythis,"h_trk_refpt_g_cent%d",ibin);
    h_trk_refpt_g[ibin] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,40,0.,40.);
    h_trk_refpt_g[ibin]->Sumw2();

    sprintf(saythis,"h_trk_refpt_u_cent%d",ibin);
    h_trk_refpt_u[ibin] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,40,0.,40.);
    h_trk_refpt_u[ibin]->Sumw2();

    sprintf(saythis,"h_bkgtrk_refpt_cent%d",ibin);
    h_bkgtrk_refpt[ibin] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,40,0.,40.);
    h_bkgtrk_refpt[ibin]->Sumw2();

    sprintf(saythis,"h_sube0trk_refpt_cent%d",ibin);
    h_sube0trk_refpt[ibin] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,40,0.,40.);
    h_sube0trk_refpt[ibin]->Sumw2();

    sprintf(saythis,"h_sampled_bkg_cent%d",ibin);
    h_sampled_bkg[ibin] = new TH1D(saythis,"",500,-2.555,2.445);
    h_sampled_bkg[ibin]->Sumw2();

    sprintf(saythis,"h_chg_refpt_subenon0_cent%d",ibin);
    h_chg_refpt_subenon0[ibin] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
    h_chg_refpt_subenon0[ibin]->Sumw2();

    for (int ibin2=0;ibin2<ntrkBins;ibin2++){
      sprintf(saythis,"h_trk_corrpt_cent%d_trk%d",ibin,ibin2);
      h_trk_corrpt[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,40,0.,40.);
      h_trk_corrpt[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_bkgtrk_corrpt_cent%d_trk%d",ibin,ibin2);
      h_bkgtrk_corrpt[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,40,0.,40.);
      h_bkgtrk_corrpt[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_reco_corr_no0trks_cent%d_trk%d",ibin,ibin2);
      h_reco_corr_no0trks[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_reco_corr_no0trks[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_reco_corr_no0trks_q_cent%d_trk%d",ibin,ibin2);
      h_reco_corr_no0trks_q[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_reco_corr_no0trks_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_reco_corr_no0trks_g_cent%d_trk%d",ibin,ibin2);
      h_reco_corr_no0trks_g[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_reco_corr_no0trks_g[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_reco_corr_0trks_cent%d_trk%d",ibin,ibin2);
      h_reco_corr_0trks[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_reco_corr_0trks[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_reco_corr_0trks_q_cent%d_trk%d",ibin,ibin2);
      h_reco_corr_0trks_q[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_reco_corr_0trks_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_reco_corr_0trks_g_cent%d_trk%d",ibin,ibin2);
      h_reco_corr_0trks_g[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_reco_corr_0trks_g[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_reco_corr_0trks_up_cent%d_trk%d",ibin,ibin2);
      h_reco_corr_0trks_up[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_reco_corr_0trks_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_reco_corr_0trks_down_cent%d_trk%d",ibin,ibin2);
      h_reco_corr_0trks_down[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_reco_corr_0trks_down[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_reco_corr_0trks_others_cent%d_trk%d",ibin,ibin2);
      h_reco_corr_0trks_others[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_reco_corr_0trks_others[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_gen_full_no0trks_cent%d_trk%d",ibin,ibin2);
      h_gen_full_no0trks[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_gen_full_no0trks[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_gen_full_no0trks_q_cent%d_trk%d",ibin,ibin2);
      h_gen_full_no0trks_q[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_gen_full_no0trks_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_gen_full_no0trks_g_cent%d_trk%d",ibin,ibin2);
      h_gen_full_no0trks_g[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_gen_full_no0trks_g[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_gen_full_0trks_cent%d_trk%d",ibin,ibin2);
      h_gen_full_0trks[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_gen_full_0trks[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_gen_full_0trks_q_cent%d_trk%d",ibin,ibin2);
      h_gen_full_0trks_q[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_gen_full_0trks_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_gen_full_0trks_g_cent%d_trk%d",ibin,ibin2);
      h_gen_full_0trks_g[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_gen_full_0trks_g[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_gen_full_0trks_up_cent%d_trk%d",ibin,ibin2);
      h_gen_full_0trks_up[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_gen_full_0trks_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_gen_full_0trks_down_cent%d_trk%d",ibin,ibin2);
      h_gen_full_0trks_down[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_gen_full_0trks_down[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_gen_full_0trks_others_cent%d_trk%d",ibin,ibin2);
      h_gen_full_0trks_others[ibin][ibin2] = new TH1D(saythis,"",jt_nbins,jt_bin_bounds);
      h_gen_full_0trks_others[ibin][ibin2]->Sumw2();


 ///unfolding histograms
      sprintf(saythis,"h_chg_recopt_genpt_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_genpt[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_genpt[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recoptsube0_genptsube0_cent%d_trk%d",ibin,ibin2);
      h_chg_recoptsube0_genptsube0[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recoptsube0_genptsube0[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_recoptsube0_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_recoptsube0[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_recoptsube0[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_genptsube0_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_genptsube0[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_genptsube0[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_effgenpt_genptsube0_cent%d_trk%d",ibin,ibin2);
      h_chg_effgenpt_genptsube0[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_effgenpt_genptsube0[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recorecopt_gengensube0pt_cent%d_trk%d",ibin,ibin2);
      h_chg_recorecopt_gengensube0pt[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recorecopt_gengensube0pt[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recoptsube0_genptsube0_up_cent%d_trk%d",ibin,ibin2);
      h_chg_recoptsube0_genptsube0_up[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recoptsube0_genptsube0_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_recoptsube0_up_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_recoptsube0_up[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_recoptsube0_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_genptsube0_up_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_genptsube0_up[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_genptsube0_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_effgenpt_genptsube0_up_cent%d_trk%d",ibin,ibin2);
      h_chg_effgenpt_genptsube0_up[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_effgenpt_genptsube0_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recorecopt_gengensube0pt_up_cent%d_trk%d",ibin,ibin2);
      h_chg_recorecopt_gengensube0pt_up[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recorecopt_gengensube0pt_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recoptsube0_genptsube0_d_cent%d_trk%d",ibin,ibin2);
      h_chg_recoptsube0_genptsube0_d[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recoptsube0_genptsube0_d[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_recoptsube0_d_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_recoptsube0_d[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_recoptsube0_d[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_genptsube0_d_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_genptsube0_d[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_genptsube0_d[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_effgenpt_genptsube0_d_cent%d_trk%d",ibin,ibin2);
      h_chg_effgenpt_genptsube0_d[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_effgenpt_genptsube0_d[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recorecopt_gengensube0pt_d_cent%d_trk%d",ibin,ibin2);
      h_chg_recorecopt_gengensube0pt_d[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recorecopt_gengensube0pt_d[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recoptsube0_genptsube0_g_oth_cent%d_trk%d",ibin,ibin2);
      h_chg_recoptsube0_genptsube0_g_oth[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recoptsube0_genptsube0_g_oth[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_recoptsube0_g_oth_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_recoptsube0_g_oth[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_recoptsube0_g_oth[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_genptsube0_g_oth_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_genptsube0_g_oth[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_genptsube0_g_oth[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_effgenpt_genptsube0_g_oth_cent%d_trk%d",ibin,ibin2);
      h_chg_effgenpt_genptsube0_g_oth[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_effgenpt_genptsube0_g_oth[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recorecopt_gengensube0pt_g_oth_cent%d_trk%d",ibin,ibin2);
      h_chg_recorecopt_gengensube0pt_g_oth[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recorecopt_gengensube0pt_g_oth[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_datasample_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_datasample[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_datasample[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_eff_genpt_cent%d_trk%d",ibin,ibin2);
      h_chg_eff_genpt[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_eff_genpt[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_bkg_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_bkg[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_bkg[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_bkg_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_bkg[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_bkg[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_sube0_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_sube0[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_sube0[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_sube0_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_sube0[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_sube0[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_gengensube0pt_cent%d_trk%d",ibin,ibin2);
      h_chg_gengensube0pt[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_gengensube0pt[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_sube0_bkgadd_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_sube0_bkgadd[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_sube0_bkgadd[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_sube0_bkgadd_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_sube0_bkgadd[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_sube0_bkgadd[ibin][ibin2]->Sumw2();

 ///up
      sprintf(saythis,"h_chg_recopt_up_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_up[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_up_datasample_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_up_datasample[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_up_datasample[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_up_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_up[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_eff_genpt_up_cent%d_trk%d",ibin,ibin2);
      h_chg_eff_genpt_up[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_eff_genpt_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_sube0_up_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_sube0_up[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_sube0_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_sube0_up_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_sube0_up[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_sube0_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_gengensube0pt_up_cent%d_trk%d",ibin,ibin2);
      h_chg_gengensube0pt_up[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_gengensube0pt_up[ibin][ibin2]->Sumw2();

 ///down
      sprintf(saythis,"h_chg_recopt_d_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_d[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_d[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_d_datasample_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_d_datasample[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_d_datasample[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_d_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_d[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_d[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_eff_genpt_d_cent%d_trk%d",ibin,ibin2);
      h_chg_eff_genpt_d[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_eff_genpt_d[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_sube0_d_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_sube0_d[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_sube0_d[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_sube0_d_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_sube0_d[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_sube0_d[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_gengensube0pt_d_cent%d_trk%d",ibin,ibin2);
      h_chg_gengensube0pt_d[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_gengensube0pt_d[ibin][ibin2]->Sumw2();

 ///g_oth
      sprintf(saythis,"h_chg_recopt_g_oth_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_g_oth[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_g_oth[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_g_oth_datasample_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_g_oth_datasample[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_g_oth_datasample[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_g_oth_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_g_oth[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_g_oth[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_eff_genpt_g_oth_cent%d_trk%d",ibin,ibin2);
      h_chg_eff_genpt_g_oth[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_eff_genpt_g_oth[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_sube0_g_oth_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_sube0_g_oth[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_sube0_g_oth[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_sube0_g_oth_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_sube0_g_oth[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_sube0_g_oth[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_gengensube0pt_g_oth_cent%d_trk%d",ibin,ibin2);
      h_chg_gengensube0pt_g_oth[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_gengensube0pt_g_oth[ibin][ibin2]->Sumw2();

///////////////////////////////////kappa histograms//////////////////////////////////////////////////
 ///unfolding histograms
      sprintf(saythis,"h_chg_recopt_genpt_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_genpt_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_genpt_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recoptsube0_genptsube0_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recoptsube0_genptsube0_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recoptsube0_genptsube0_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_recoptsube0_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_recoptsube0_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_recoptsube0_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_genptsube0_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_genptsube0_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_genptsube0_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_effgenpt_genptsube0_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_effgenpt_genptsube0_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_effgenpt_genptsube0_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recorecopt_gengensube0pt_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recorecopt_gengensube0pt_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recorecopt_gengensube0pt_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recoptsube0_genptsube0_up_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recoptsube0_genptsube0_up_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recoptsube0_genptsube0_up_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_recoptsube0_up_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_recoptsube0_up_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_recoptsube0_up_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_genptsube0_up_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_genptsube0_up_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_genptsube0_up_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_effgenpt_genptsube0_up_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_effgenpt_genptsube0_up_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_effgenpt_genptsube0_up_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recorecopt_gengensube0pt_up_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recorecopt_gengensube0pt_up_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recorecopt_gengensube0pt_up_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recoptsube0_genptsube0_d_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recoptsube0_genptsube0_d_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recoptsube0_genptsube0_d_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_recoptsube0_d_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_recoptsube0_d_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_recoptsube0_d_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_genptsube0_d_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_genptsube0_d_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_genptsube0_d_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_effgenpt_genptsube0_d_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_effgenpt_genptsube0_d_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_effgenpt_genptsube0_d_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recorecopt_gengensube0pt_d_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recorecopt_gengensube0pt_d_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recorecopt_gengensube0pt_d_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recoptsube0_genptsube0_g_oth_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recoptsube0_genptsube0_g_oth_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recoptsube0_genptsube0_g_oth_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_recoptsube0_g_oth_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_recoptsube0_g_oth_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_recoptsube0_g_oth_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_genptsube0_g_oth_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_genptsube0_g_oth_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recopt_genptsube0_g_oth_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_effgenpt_genptsube0_g_oth_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_effgenpt_genptsube0_g_oth_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_effgenpt_genptsube0_g_oth_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recorecopt_gengensube0pt_g_oth_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recorecopt_gengensube0pt_g_oth_kappa[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,500,-2.555,2.445);
      h_chg_recorecopt_gengensube0pt_g_oth_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_kappa_datasample_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_kappa_datasample[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_kappa_datasample[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_eff_genpt_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_eff_genpt_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_eff_genpt_kappa[ibin][ibin2]->Sumw2();
/*
      sprintf(saythis,"h_chg_recopt_bkg_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_bkg_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_bkg_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_bkg_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_bkg_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_bkg_kappa[ibin][ibin2]->Sumw2();
*/
      sprintf(saythis,"h_chg_recopt_sube0_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_sube0_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_sube0_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_sube0_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_sube0_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_sube0_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_gengensube0pt_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_gengensube0pt_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_gengensube0pt_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_sube0_bkgadd_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_sube0_bkgadd_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_sube0_bkgadd_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_sube0_bkgadd_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_sube0_bkgadd_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_sube0_bkgadd_kappa[ibin][ibin2]->Sumw2();

 ///up
      sprintf(saythis,"h_chg_recopt_up_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_up_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_up_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_up_kappa_datasample_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_up_kappa_datasample[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_up_kappa_datasample[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_up_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_up_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_up_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_eff_genpt_up_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_eff_genpt_up_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_eff_genpt_up_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_sube0_up_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_sube0_up_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_sube0_up_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_sube0_up_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_sube0_up_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_sube0_up_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_gengensube0pt_up_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_gengensube0pt_up_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_gengensube0pt_up_kappa[ibin][ibin2]->Sumw2();

 ///down
      sprintf(saythis,"h_chg_recopt_d_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_d_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_d_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_d_kappa_datasample_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_d_kappa_datasample[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_d_kappa_datasample[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_d_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_d_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_d_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_eff_genpt_d_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_eff_genpt_d_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_eff_genpt_d_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_sube0_d_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_sube0_d_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_sube0_d_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_sube0_d_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_sube0_d_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_sube0_d_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_gengensube0pt_d_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_gengensube0pt_d_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_gengensube0pt_d_kappa[ibin][ibin2]->Sumw2();

 ///g_oth
      sprintf(saythis,"h_chg_recopt_g_oth_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_g_oth_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_g_oth_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_g_oth_kappa_datasample_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_g_oth_kappa_datasample[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_g_oth_kappa_datasample[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_g_oth_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_g_oth_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_g_oth_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_eff_genpt_g_oth_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_eff_genpt_g_oth_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_eff_genpt_g_oth_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_recopt_sube0_g_oth_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_recopt_sube0_g_oth_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_recopt_sube0_g_oth_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_genpt_sube0_g_oth_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_genpt_sube0_g_oth_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_genpt_sube0_g_oth_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_gengensube0pt_g_oth_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_gengensube0pt_g_oth_kappa[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
      h_chg_gengensube0pt_g_oth_kappa[ibin][ibin2]->Sumw2();
    }

    sprintf(saythis,"h_trk_corrpt_q_cent%d",ibin);
    h_trk_corrpt_q[ibin] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,40,0.,40.);
    h_trk_corrpt_q[ibin]->Sumw2();

    sprintf(saythis,"h_trk_corrpt_g_cent%d",ibin);
    h_trk_corrpt_g[ibin] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,40,0.,40.);
    h_trk_corrpt_g[ibin]->Sumw2();

    sprintf(saythis,"h_trk_corrpt_u_cent%d",ibin);
    h_trk_corrpt_u[ibin] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,40,0.,40.);
    h_trk_corrpt_u[ibin]->Sumw2();

    sprintf(saythis,"h_chg_corrpt_highest_cent%d",ibin);
    h_chg_corrpt_highest[ibin] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
    h_chg_corrpt_highest[ibin]->Sumw2();    

    sprintf(saythis,"h_chg_corrpt_highest_2_cent%d",ibin);
    h_chg_corrpt_highest_2[ibin] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
    h_chg_corrpt_highest_2[ibin]->Sumw2();    

    for (int ibin2=0;ibin2<ntrkBins;ibin2++){

      sprintf(saythis,"h_chg_refpt_eff_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_eff[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_eff[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_sube0_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_sube0[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_sube0[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_sube0_bkgadd_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_sube0_bkgadd[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_sube0_bkgadd[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_nobkgsub_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_nobkgsub[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_nobkgsub[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_q_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_q[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_g_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_g[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_g[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_bkg_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_bkg[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_bkg[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_bkg_q_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_bkg_q[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_bkg_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_bkg_g_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_bkg_g[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_bkg_g[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_u_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_u[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_u[ibin][ibin2]->Sumw2();
  
      for(int ibin3=0; ibin3<trkeffbins; ibin3++){
          sprintf(saythis,"h_chg_corrpt_trkeff_cent%d_trk%d_eff%d",ibin,ibin2,ibin3);
          h_chg_corrpt_trkeff[ibin][ibin2][ibin3] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
          h_chg_corrpt_trkeff[ibin][ibin2][ibin3]->Sumw2();

          sprintf(saythis,"h_chg_corrpt_kappa_trkeff_cent%d_trk%d_eff%d",ibin,ibin2,ibin3);
          h_chg_corrpt_kappa_trkeff[ibin][ibin2][ibin3] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
          h_chg_corrpt_kappa_trkeff[ibin][ibin2][ibin3]->Sumw2();
      }

      sprintf(saythis,"h_chg_sig_bkg_cent%d_trk%d",ibin,ibin2);
      h_chg_sig_bkg[ibin][ibin2] = new TH2F(saythis,"",40.,0.,40.,100,0.,40.);
      h_chg_sig_bkg[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_sig_ntrk_cent%d_trk%d",ibin,ibin2);
      h_chg_sig_ntrk[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,40,0.,40.);
      h_chg_sig_ntrk[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_bkg_nbkg_cent%d_trk%d",ibin,ibin2);
      h_chg_bkg_nbkg[ibin][ibin2] = new TH2F(saythis,"",500,-2.555,2.445,40,0.,40.);
      h_chg_bkg_nbkg[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_sig_ntrk_pos_cent%d_trk%d",ibin,ibin2);
      h_chg_sig_ntrk_pos[ibin][ibin2] = new TH2F(saythis,"",40,0.,40.,40,0.,40.);
      h_chg_sig_ntrk_pos[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_bkg_nbkg_pos_cent%d_trk%d",ibin,ibin2);
      h_chg_bkg_nbkg_pos[ibin][ibin2] = new TH2F(saythis,"",40,0.,40.,40,0.,40.);
      h_chg_bkg_nbkg_pos[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_sig_ntrk_neg_cent%d_trk%d",ibin,ibin2);
      h_chg_sig_ntrk_neg[ibin][ibin2] = new TH2F(saythis,"",40,-40.,0.,40,0.,40.);
      h_chg_sig_ntrk_neg[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_bkg_nbkg_neg_cent%d_trk%d",ibin,ibin2);
      h_chg_bkg_nbkg_neg[ibin][ibin2] = new TH2F(saythis,"",40,-40.,0.,40,0.,40.);
      h_chg_bkg_nbkg_neg[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_nobkgsub_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_nobkgsub[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_nobkgsub[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_q_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_q[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_g_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_g[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_g[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_mix_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_mix[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_mix[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_bkg_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_bkg[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_bkg[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_bkg_q_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_bkg_q[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_bkg_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_bkg_g_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_bkg_g[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_bkg_g[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_u_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_u[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_u[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_up_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_up[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_upbar_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_upbar[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_upbar[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_d_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_d[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_d[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_dbar_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_dbar[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_dbar[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_c_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_c[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_c[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_s_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_s[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_s[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_b_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_b[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_b[ibin][ibin2]->Sumw2();

      //////kappa histos
      sprintf(saythis,"h_chg_corrpt_kappa_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_kappa_q_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_q_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_q_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_kappa_g_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_g_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_g_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_kappa_bkg_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_bkg_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_bkg_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_kappa_bkg_q_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_bkg_q_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_bkg_q_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_kappa_bkg_g_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_bkg_g_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_bkg_g_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_kappa_u_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_u_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_u_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_kappa_up_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_up_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_up_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_kappa_upbar_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_upbar_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_upbar_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_kappa_d_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_d_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_d_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_kappa_dbar_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_dbar_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_dbar_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_kappa_c_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_c_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_c_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_kappa_s_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_s_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_s_kappa[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_kappa_b_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_b_kappa[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_corrpt_b_kappa[ibin][ibin2]->Sumw2();
      //////////////////

      sprintf(saythis,"h_chg_refpt_up_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_up[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_upbar_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_upbar[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_upbar[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_d_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_d[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_d[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_dbar_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_dbar[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_dbar[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_c_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_c[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_c[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_s_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_s[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_s[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_b_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_b[ibin][ibin2] = new TH2F(saythis,"",jt_nbins,jt_bin_bounds,500,-2.555,2.445);
      h_chg_refpt_b[ibin][ibin2]->Sumw2();
    }

    sprintf(saythis,"h_jt_gen_closure_cent%d",ibin);
    h_jt_gen_cloure[ibin] = new TH2F(saythis, "", 40, 100., 500., 100, 0., 2.);
    h_jt_gen_cloure[ibin]->Sumw2();

    sprintf(saythis,"h_jt_gen_q_closure_cent%d",ibin);
    h_jt_gen_q_cloure[ibin] = new TH2F(saythis, "", 40, 100., 500., 100, 0., 2.);
    h_jt_gen_q_cloure[ibin]->Sumw2();

    sprintf(saythis,"h_jt_gen_g_closure_cent%d",ibin);
    h_jt_gen_g_cloure[ibin] = new TH2F(saythis, "", 40, 100., 500., 100, 0., 2.);
    h_jt_gen_g_cloure[ibin]->Sumw2();

    //f_res_gauss[ibin] = new TF1((TString)("f_gauss"+cent[ibin]), "exp(-((x-1.)^2)/(2*([0]^2)))", 0., 2.);
    //f_res_gauss[ibin]->FixParameter(0,0.0255);

    sprintf(saythis,"h_trkpt_reco_cent%d",ibin);
    h_trkpt_reco[ibin] = new TH1D(saythis, "", 20,0.,10.);
    h_trkpt_reco[ibin]->Sumw2();

    sprintf(saythis,"h_trkpt_pos_reco_cent%d",ibin);
    h_trkpt_pos_reco[ibin] = new TH1D(saythis, "", 20,0.,10.);
    h_trkpt_pos_reco[ibin]->Sumw2();

    sprintf(saythis,"h_trkpt_neg_reco_cent%d",ibin);
    h_trkpt_neg_reco[ibin] = new TH1D(saythis, "", 20,0.,10.);
    h_trkpt_neg_reco[ibin]->Sumw2();    

    sprintf(saythis,"h_trkpt_gen_cent%d",ibin);
    h_trkpt_gen[ibin] = new TH1D(saythis, "", 20,0.,10.);
    h_trkpt_gen[ibin]->Sumw2();

    sprintf(saythis,"h_trkpt_pos_gen_cent%d",ibin);
    h_trkpt_pos_gen[ibin] = new TH1D(saythis, "", 20,0.,10.);
    h_trkpt_pos_gen[ibin]->Sumw2();

    sprintf(saythis,"h_trkpt_neg_gen_cent%d",ibin);
    h_trkpt_neg_gen[ibin] = new TH1D(saythis, "", 20,0.,10.);
    h_trkpt_neg_gen[ibin]->Sumw2();    

  }

  int total_n_jets = 0;

  ///// CUTS ///////
  const double etamaxcut = 1.5;
  const double etamincut = 0.5;
  const double pTmincut = 50.;
  const double pTmaxcut = 600.;
  const double refpTmincut = 50.;

  ///////////////// centrality reweighting ///////////////////////

  //cymbal
  TF1* f_cent= new TF1("f_cent","[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4+[7]*exp([5]+[6]*x)",0,180);  
  f_cent->SetParameters(4.40810, -7.75301e-02, 4.91953e-04, -1.34961e-06, 1.44407e-09, -160, 1, 3.68078e-15);

  ///////////////// vz reweighting ///////////////////////

  //Pythia6
  TF1 *fvz_pp = new TF1("fvz_pp","gaus",-15,15); 
  fvz_pp->SetParameter(0,1.10477);
  fvz_pp->SetParameter(1,2.52738);
  fvz_pp->SetParameter(2,1.30296e1);

  //Pythia8 vz reweighting
  TF1 *fvz_pythia8 = new TF1("fvz_rw","pol4",-15.,15.); 
  fvz_pythia8->SetParameters(9.22005e-01,-1.42259e-02,3.47708e-03,-2.73684e-05,-4.09970e-06);

  //cymbal
  TF1 *fvz_cymbal = new TF1("fvz","pol6",-15,15); 
  fvz_cymbal->SetParameters(1.18472,-0.132675,0.00857998,-0.000326085,-1.48786e-06,4.68665e-07,-7.32942e-09 );
    
  //pyquen 
  TF1 *fvz_pyquen = new TF1("fvz","pol4",-15,15);   
  fvz_pyquen->SetParameters(1.27257,8.45372e-03,-8.32117e-03,-3.41576e-05,1.60058e-05);

  ///skims
  TFile *my_file;
  int nfiles;
  if(!ispp && !isdata) nfiles = 9;
  else nfiles = 2;

  for(int file=1; file<nfiles; file++){

    TFile *my_file;  
    if(isdata){
      if(ispp) my_file = TFile::Open("skims/arc_skims/ppdata_jetcharge_skim_Mar6.root");
      else my_file = TFile::Open("skims/arc_skims/PbPbdata2015_jetchargeskims_Mar18.root");
    }
    else{
      if(ispp) my_file = TFile::Open("skims/arc_skims/Pythia_jetcharge_skim_Mar6.root"); 
      //else my_file = TFile::Open("skims/arc_skims/PbPbMC_gensmear_jetchg_Feb28.root");    
      else my_file = TFile::Open(Form("skims/arc_skims/PbPbMC_P+H_jetchargeskims_Apr12/HydJet_added_%d.root",file));
    }

    cout<<"got file "<<file<<" out of "<<nfiles<<endl;

    TTree *inp_tree = (TTree*)my_file->Get("mixing_tree");

    inp_tree->SetBranchAddress("vz",&vz);
    if(!isdata) inp_tree->SetBranchAddress("pthat",&pthat);
    if(!ispp) inp_tree->SetBranchAddress("hiBin",&hiBin);

    inp_tree->SetBranchAddress("calo_jtpt",&calo_jtpt);
    inp_tree->SetBranchAddress("calo_corrpt",&calo_corrpt);
    inp_tree->SetBranchAddress("calo_jteta",&calo_jteta);
    inp_tree->SetBranchAddress("calo_jtphi",&calo_jtphi);
      
    if(!isdata){
      inp_tree->SetBranchAddress("calo_refpt",&calo_refpt);
      inp_tree->SetBranchAddress("calo_refeta",&calo_refeta);
      inp_tree->SetBranchAddress("calo_refphi",&calo_refphi);
      inp_tree->SetBranchAddress("calo_refparton_flavor",&calo_refparton_flavor);
    }

    inp_tree->SetBranchAddress("trkChg",&trkChg);
    inp_tree->SetBranchAddress("trkPt",&trkPt);
    inp_tree->SetBranchAddress("trkEta",&trkEta);
    inp_tree->SetBranchAddress("trkPhi",&trkPhi);

    if(!isdata){
      inp_tree->SetBranchAddress("pt",&genPt);
      inp_tree->SetBranchAddress("eta",&genEta);
      inp_tree->SetBranchAddress("phi",&genPhi);
      inp_tree->SetBranchAddress("chg",&genChg);
      inp_tree->SetBranchAddress("sube",&gensube);
      inp_tree->SetBranchAddress("pdg",&genpdg);
    }

    cout << "Retrieved tree from input file!" << endl;
    Long64_t n_evts = inp_tree->GetEntriesFast();
    cout<<"n evts: "<<n_evts<<endl;
    int evi_frac=0;

    for (int evi = 0; evi < n_evts; evi++){

      inp_tree->GetEntry(evi);
      evi_frac = 100*evi/n_evts;
      if(evi%100000==0) cout<<"evi frac: "<<evi_frac<<"%"<<endl;

      bool fill_data_sample = false;
      //if(gRandom->Rndm() >= 0.5) fill_data_sample = true;

      if (fabs(vz) > 15.) continue;

      //pythia6 
      if(ispp) weight_vz = 1./fvz_pp->Eval(vz);
      //Pythia6+Hydjet
      else weight_vz = fvz_cymbal->Eval(vz);

      if(isdata) weight_vz = 1.;  

      ///// pthat weight
      if(!isdata && pthat < 80.) continue;

      int ibin=0;
      if(ispp){
        while(pthat>pthatbins[ibin+1]) ibin++;
        pthat_weight = (xsecs[ibin]-xsecs[ibin+1])/pthatEntries[ibin];
      }
      else{
        while(pthat>pthatbins_PH[ibin+1]) ibin++;
        pthat_weight = (xsecs_PH[ibin]-xsecs_PH[ibin+1])/pthatEntries_PH[ibin];
      }

      if(isdata) pthat_weight = 1.;

      //// centrality bin and weight 
      if(ispp) hiBin = 1;

      double weight_cen = fcent_cymbal(hiBin,f_cent);
      if(ispp || isdata) weight_cen = 1.;

      for (int cbin = 0; cbin < nCBins; cbin++){ 
        if (hiBin > CBins[cbin] && hiBin <= CBins[cbin+1]){
          mycbin = cbin; 
        }
      }

      //jet loop
      for (int jet = 0; jet < (int) calo_jtpt->size(); jet++){ 

        //if(fabs(calo_jteta->at(jet)) >= etamaxcut || fabs(calo_jteta->at(jet)) <= etamincut) continue ; 
        if(fabs(calo_jteta->at(jet)) >= etamaxcut) continue ; 
        
        jet_corrpt = calo_corrpt->at(jet);        
        //jet_corrpt = jet_corrpt*f_res_gauss->GetRandom();

        if(jet_corrpt <= pTmincut || jet_corrpt >= pTmaxcut) continue;

        if(!isdata){
            if (calo_refpt->at(jet) <= refpTmincut) continue;
        }

        for(int ibin2=0; ibin2<ntrkBins; ibin2++){          
            trk_chg_ptcut[ibin2]=0.;
            trk_bkg_ptcut[ibin2]=0.;

            trk_chg_kappa[ibin2]=0.;
            trk_bkg_kappa[ibin2]=0.;

            n_trk[ibin2]=0;

            for(int effbin=0; effbin<trkeffbins; effbin++){
                trk_chg_ptcut_trkeff[ibin2][effbin]=0.;
                trk_chg_kappa_trkeff[ibin2][effbin]=0.;
            }
        }

        //reco track loop 
        for (int trk = 0; trk < (int) trkPt->size(); trk++){
            h_trkpt_reco[mycbin]->Fill(trkPt->at(trk),weight_vz*weight_cen*pthat_weight);
            if(trkChg->at(trk) == 1) h_trkpt_pos_reco[mycbin]->Fill(trkPt->at(trk),weight_vz*weight_cen*pthat_weight);
            else if (trkChg->at(trk) == -1) h_trkpt_neg_reco[mycbin]->Fill(trkPt->at(trk),weight_vz*weight_cen*pthat_weight);

            double dr = sqrt(pow(trkEta->at(trk) - calo_jteta->at(jet), 2) + pow(acos(cos(calo_jtphi->at(jet) - trkPhi->at(trk))),2));
            double refl_dr = sqrt(pow(trkEta->at(trk) + calo_jteta->at(jet), 2) + pow(acos(cos(calo_jtphi->at(jet) - trkPhi->at(trk))),2));

            if(dr <= 0.4){ 
                if(trkPt->at(trk) >= 2.){
                  n_trk[0]++;
                  trk_chg_ptcut[0] += pow(trkPt->at(trk),0.5)*trkChg->at(trk);    
                  
                  trk_chg_kappa[0] += pow(trkPt->at(trk),0.3)*trkChg->at(trk);    
                  trk_chg_kappa[1] += pow(trkPt->at(trk),0.5)*trkChg->at(trk);    
                  trk_chg_kappa[2] += pow(trkPt->at(trk),0.7)*trkChg->at(trk);    
                }
                if(trkPt->at(trk) >= 4.){
                    n_trk[1]++;
                    trk_chg_ptcut[1] += pow(trkPt->at(trk),0.5)*trkChg->at(trk);    
                }
                if(trkPt->at(trk) >= 5.){ 
                    n_trk[2]++;        
                    trk_chg_ptcut[2] += pow(trkPt->at(trk),0.5)*trkChg->at(trk);    
                }
            }
            else if(refl_dr <= 0.4){
                if(trkPt->at(trk) >= 2.){
                  trk_bkg_ptcut[0] += pow(trkPt->at(trk),0.5)*trkChg->at(trk);                    

                  trk_bkg_kappa[0] += pow(trkPt->at(trk),0.3)*trkChg->at(trk);
                  trk_bkg_kappa[1] += pow(trkPt->at(trk),0.5)*trkChg->at(trk);
                  trk_bkg_kappa[2] += pow(trkPt->at(trk),0.7)*trkChg->at(trk);
                }
                if(trkPt->at(trk) >= 4.){
                  trk_bkg_ptcut[1] += pow(trkPt->at(trk),0.5)*trkChg->at(trk);                    
                }
                if(trkPt->at(trk) >= 5.){
                  trk_bkg_ptcut[2] += pow(trkPt->at(trk),0.5)*trkChg->at(trk);                    
                }
            }
        }//trk loop 

        if(!isdata){
          for(int ibin2=0; ibin2<ntrkBins; ibin2++){          
              gen_chg_ptcut[ibin2]=0.;
              gen_bkg_ptcut[ibin2]=0.;
              gen_chg_ptcut_sube0[ibin2]=0.;
              gen_chg_ptcut_eff[ibin2]=0.;
              gen_chg_ptcut_sube0_eff[ibin2]=0.;
              gengen_chg_ptcut_sube0[ibin2]=0.;

              gen_chg_kappa[ibin2]=0.;
              gen_chg_kappa_sube0[ibin2]=0.;
              gen_chg_kappa_eff[ibin2]=0.;
              gen_chg_kappa_sube0_eff[ibin2]=0.;
              gengen_chg_kappa_sube0[ibin2]=0.;

              n_gen[ibin2]=0;
              n_eff_gen[ibin2]=0;
              n_sube0_gen[ibin2]=0;
              n_eff_sube0_gen[ibin2]=0;
              n_sube0_gengen[ibin2]=0;
          }

          //gen track loop 
          for (int gen = 0; gen < (int) genPt->size(); gen++){
              bool no_track=0;
              if(genPt->at(gen) > 400.) continue; 

              h_trkpt_gen[mycbin]->Fill(genPt->at(gen),weight_vz*weight_cen*pthat_weight);
              if(genChg->at(gen) == 1) h_trkpt_pos_gen[mycbin]->Fill(genPt->at(gen),weight_vz*weight_cen*pthat_weight);
              else if (genChg->at(gen) == -1) h_trkpt_neg_gen[mycbin]->Fill(genPt->at(gen),weight_vz*weight_cen*pthat_weight);

  //gen jet loop
              double gen_dr = sqrt(pow(genEta->at(gen) - calo_refeta->at(jet), 2) + pow(acos(cos(calo_refphi->at(jet) - genPhi->at(gen))),2));

              if(gen_dr <= 0.4){
                if(gensube->at(gen) == 0) {
                    if(genPt->at(gen) >= 2.){
                        n_sube0_gengen[0]++;
                        gengen_chg_ptcut_sube0[0] += pow(genPt->at(gen),0.5)*genChg->at(gen);     

                        gengen_chg_kappa_sube0[0] += pow(genPt->at(gen),0.3)*genChg->at(gen);    
                        gengen_chg_kappa_sube0[1] += pow(genPt->at(gen),0.5)*genChg->at(gen);    
                        gengen_chg_kappa_sube0[2] += pow(genPt->at(gen),0.7)*genChg->at(gen);    
                    }
                    if(genPt->at(gen) >= 4.){
                        n_sube0_gengen[1]++;
                        gengen_chg_ptcut_sube0[1] += pow(genPt->at(gen),0.5)*genChg->at(gen);    
                    }
                    if(genPt->at(gen) >= 5.){ 
                        n_sube0_gengen[2]++;        
                        gengen_chg_ptcut_sube0[2] += pow(genPt->at(gen),0.5)*genChg->at(gen);    
                    }
                }
              }


              double dr = sqrt(pow(genEta->at(gen) - calo_jteta->at(jet), 2) + pow(acos(cos(calo_jtphi->at(jet) - genPhi->at(gen))),2));
              double refl_dr = sqrt(pow(genEta->at(gen) + calo_jteta->at(jet), 2) + pow(acos(cos(calo_jtphi->at(jet) - genPhi->at(gen))),2));

              float rmin = 999.;  
              if(ispp) trk_correction = pptrkCorr->getTrkCorr(genPt->at(gen),genEta->at(gen),genPhi->at(gen),1,rmin); 
              else trk_correction = PbPbtrkCorr->getTrkCorr(genPt->at(gen),genEta->at(gen),genPhi->at(gen),hiBin); 
              //trk_correction = 1.04*trk_correction;
              trk_eff = 1./trk_correction;
              if(gRandom->Rndm() > trk_eff) no_track = 1; 

              if(dr <= 0.4){

                  if(genPt->at(gen) >= 2.){
                      n_gen[0]++;
                      gen_chg_ptcut[0] += pow(genPt->at(gen),0.5)*genChg->at(gen);    

                      gen_chg_kappa[0] += pow(genPt->at(gen),0.3)*genChg->at(gen);    
                      gen_chg_kappa[1] += pow(genPt->at(gen),0.5)*genChg->at(gen);    
                      gen_chg_kappa[2] += pow(genPt->at(gen),0.7)*genChg->at(gen);    

                      if(no_track == 0) {
                          n_eff_gen[0]++;
                          gen_chg_ptcut_eff[0] += pow(genPt->at(gen),0.5)*genChg->at(gen);    

                          gen_chg_kappa_eff[0] += pow(genPt->at(gen),0.3)*genChg->at(gen);    
                          gen_chg_kappa_eff[1] += pow(genPt->at(gen),0.5)*genChg->at(gen);    
                          gen_chg_kappa_eff[2] += pow(genPt->at(gen),0.7)*genChg->at(gen);    
                      }

                      if(gensube->at(gen) == 0) {
                          n_sube0_gen[0]++;
                          gen_chg_ptcut_sube0[0] += pow(genPt->at(gen),0.5)*genChg->at(gen);

                          gen_chg_kappa_sube0[0] += pow(genPt->at(gen),0.3)*genChg->at(gen);    
                          gen_chg_kappa_sube0[1] += pow(genPt->at(gen),0.5)*genChg->at(gen);    
                          gen_chg_kappa_sube0[2] += pow(genPt->at(gen),0.7)*genChg->at(gen);    

                          if(no_track == 0){
                              n_eff_sube0_gen[0]++;
                              gen_chg_ptcut_sube0_eff[0] += pow(genPt->at(gen),0.5)*genChg->at(gen);             

                              gen_chg_kappa_sube0_eff[0] += pow(genPt->at(gen),0.3)*genChg->at(gen);    
                              gen_chg_kappa_sube0_eff[1] += pow(genPt->at(gen),0.5)*genChg->at(gen);    
                              gen_chg_kappa_sube0_eff[2] += pow(genPt->at(gen),0.7)*genChg->at(gen);    
                          }
                      }                           
                  }
                  if(genPt->at(gen) >= 4.){
                      n_gen[1]++;
                      gen_chg_ptcut[1] += pow(genPt->at(gen),0.5)*genChg->at(gen);    

                      if(no_track == 0) {
                          n_eff_gen[1]++;
                          gen_chg_ptcut_eff[1] += pow(genPt->at(gen),0.5)*genChg->at(gen);    
                      }

                      if(gensube->at(gen) == 0) {
                          n_sube0_gen[1]++;
                          gen_chg_ptcut_sube0[1] += pow(genPt->at(gen),0.5)*genChg->at(gen);
                          if(no_track == 0){
                              n_eff_sube0_gen[1]++;
                              gen_chg_ptcut_sube0_eff[1] += pow(genPt->at(gen),0.5)*genChg->at(gen);             
                          }
                      }
                  }
                  if(genPt->at(gen) >= 5.){ 
                      n_gen[2]++;        
                      gen_chg_ptcut[2] += pow(genPt->at(gen),0.5)*genChg->at(gen);    

                      if(no_track == 0) {
                          n_eff_gen[2]++;
                          gen_chg_ptcut_eff[2] += pow(genPt->at(gen),0.5)*genChg->at(gen);    
                      }

                      if(gensube->at(gen) == 0) {
                          n_sube0_gen[2]++;
                          gen_chg_ptcut_sube0[2] += pow(genPt->at(gen),0.5)*genChg->at(gen);
                          if(no_track == 0){
                              n_eff_sube0_gen[2]++;
                              gen_chg_ptcut_sube0_eff[2] += pow(genPt->at(gen),0.5)*genChg->at(gen);             
                          }
                      }
                  }
              }
              else if(refl_dr <= 0.4){
                  if(genPt->at(gen) >= 2.){
                    gen_bkg_ptcut[0] += pow(genPt->at(gen),0.5)*genChg->at(gen);
                  }   
                  if(genPt->at(gen) >= 4.){
                    gen_bkg_ptcut[1] += pow(genPt->at(gen),0.5)*genChg->at(gen);
                  }   
                  if(genPt->at(gen) >= 5.){
                    gen_bkg_ptcut[2] += pow(genPt->at(gen),0.5)*genChg->at(gen);
                  }   
              }
          }//gen track loop 
        }

        for(int ibin2=0; ibin2<ntrkBins; ibin2++){
            trk_chg_ptcut[ibin2] = trk_chg_ptcut[ibin2]/pow(jet_corrpt,0.5);
            trk_bkg_ptcut[ibin2] = trk_bkg_ptcut[ibin2]/pow(jet_corrpt,0.5);

            trk_chg_kappa[ibin2] = trk_chg_kappa[ibin2]/pow(jet_corrpt,kappa[ibin2]);
            trk_bkg_kappa[ibin2] = trk_bkg_kappa[ibin2]/pow(jet_corrpt,kappa[ibin2]);

            if(!isdata) {
              gen_chg_ptcut[ibin2] = gen_chg_ptcut[ibin2]/pow(jet_corrpt,0.5);
              gen_chg_ptcut_eff[ibin2] = gen_chg_ptcut_eff[ibin2]/pow(jet_corrpt,0.5);
              gen_chg_ptcut_sube0[ibin2] = gen_chg_ptcut_sube0[ibin2]/pow(jet_corrpt,0.5);
              gen_chg_ptcut_sube0_eff[ibin2] = gen_chg_ptcut_sube0_eff[ibin2]/pow(jet_corrpt,0.5);
              gen_bkg_ptcut[ibin2] = gen_bkg_ptcut[ibin2]/pow(jet_corrpt,0.5);  
              gengen_chg_ptcut_sube0[ibin2] = gengen_chg_ptcut_sube0[ibin2]/pow(calo_refpt->at(jet),0.5);

              gen_chg_kappa[ibin2] = gen_chg_kappa[ibin2]/pow(jet_corrpt,kappa[ibin2]);
              gen_chg_kappa_eff[ibin2] = gen_chg_kappa_eff[ibin2]/pow(jet_corrpt,kappa[ibin2]);
              gen_chg_kappa_sube0[ibin2] = gen_chg_kappa_sube0[ibin2]/pow(jet_corrpt,kappa[ibin2]);
              gen_chg_kappa_sube0_eff[ibin2] = gen_chg_kappa_sube0_eff[ibin2]/pow(jet_corrpt,kappa[ibin2]);
              gengen_chg_kappa_sube0[ibin2] = gengen_chg_kappa_sube0[ibin2]/pow(calo_refpt->at(jet),kappa[ibin2]);
            }
        }

        /////filling histos

        if(fill_data_sample == true){
          for(int ibin2=0;ibin2<ntrkBins;ibin2++){
            if(jet_corrpt >= 120. && n_trk[ibin2] > 0){
              h_chg_recopt_datasample[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
              
              if(calo_refparton_flavor->at(jet)==2) h_chg_recopt_up_datasample[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
              else if(calo_refparton_flavor->at(jet)==1) h_chg_recopt_d_datasample[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
              else h_chg_recopt_g_oth_datasample[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
            }
            if(jet_corrpt >= 120. && n_trk[0] > 0){
              h_chg_recopt_kappa_datasample[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
              
              if(calo_refparton_flavor->at(jet)==2) h_chg_recopt_up_kappa_datasample[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
              else if(calo_refparton_flavor->at(jet)==1) h_chg_recopt_d_kappa_datasample[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
              else h_chg_recopt_g_oth_kappa_datasample[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
            }  
          }
        }
        else{
        if(!isdata){      
          h_jt_gen_cloure[mycbin]->Fill(calo_refpt->at(jet),jet_corrpt/calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
          if(calo_refparton_flavor->at(jet) == -999){
            if(calo_refpt->at(jet)>=120.) h_gen_full_u[mycbin]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
            if(jet_corrpt>=120.) {
              h_reco_corr_u[mycbin]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
              h_eta_full_u[mycbin]->Fill(calo_jteta->at(jet),pthat_weight*weight_cen*weight_vz);
            }
            continue;
          }      

          if(calo_refpt->at(jet)>=120.){
            h_trk_refpt[mycbin]->Fill(calo_refpt->at(jet),n_gen[0],pthat_weight*weight_cen*weight_vz);
            h_bkgtrk_refpt[mycbin]->Fill(calo_refpt->at(jet),n_bkggen[0],pthat_weight*weight_cen*weight_vz);
            h_sube0trk_refpt[mycbin]->Fill(calo_refpt->at(jet),n_sube0_gen[0],pthat_weight*weight_cen*weight_vz);

            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              if(n_gen[ibin2]>0) h_chg_refpt[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
              h_chg_refpt_bkg[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_bkg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);

              if(n_trk[ibin2]>0) h_gen_full_no0trks[mycbin][ibin2]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
              else if(n_trk[ibin2]==0) h_gen_full_0trks[mycbin][ibin2]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);         

              h_chg_refpt_subenon0[mycbin]->Fill(calo_refpt->at(jet),gen_chg_ptcut[1] - gen_chg_ptcut_sube0[0],pthat_weight*weight_cen*weight_vz);

              if(n_eff_gen[ibin2]>0) h_chg_refpt_eff[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_chg_ptcut_eff[ibin2],pthat_weight*weight_cen*weight_vz);
              if(n_gen[ibin2]>0) h_chg_refpt_sube0[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_chg_ptcut_sube0[ibin2],pthat_weight*weight_cen*weight_vz);
              if(n_gen[ibin2]>0) h_chg_refpt_sube0_bkgadd[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_chg_ptcut_sube0_bkgadd[ibin2],pthat_weight*weight_cen*weight_vz);

              if(n_sube0_gengen[ibin2] > 0) {
                h_chg_gengensube0pt[mycbin][ibin2]->Fill(gengen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);          
              
                if(calo_refparton_flavor->at(jet)==2) h_chg_gengensube0pt_up[mycbin][ibin2]->Fill(gengen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                else if(calo_refparton_flavor->at(jet)==1) h_chg_gengensube0pt_d[mycbin][ibin2]->Fill(gengen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                else h_chg_gengensube0pt_g_oth[mycbin][ibin2]->Fill(gengen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
              }
              if(n_sube0_gengen[0] > 0) {
                h_chg_gengensube0pt_kappa[mycbin][ibin2]->Fill(gengen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);          
              
                if(calo_refparton_flavor->at(jet)==2) h_chg_gengensube0pt_up_kappa[mycbin][ibin2]->Fill(gengen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                else if(calo_refparton_flavor->at(jet)==1) h_chg_gengensube0pt_d_kappa[mycbin][ibin2]->Fill(gengen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                else h_chg_gengensube0pt_g_oth_kappa[mycbin][ibin2]->Fill(gengen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
              }

            }
            h_gen_full[mycbin]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
          }
        }  
      
        if(jet_corrpt>=120.){

          h_ntrk[mycbin]->Fill(n_trk[0],pthat_weight*weight_cen*weight_vz);
          h_nbkgtrk[mycbin]->Fill(n_bkgtrk[0],pthat_weight*weight_cen*weight_vz);

          for(int ibin2=0;ibin2<ntrkBins;ibin2++){
            h_trk_corrpt[mycbin][ibin2]->Fill(jet_corrpt,n_trk[ibin2],pthat_weight*weight_cen*weight_vz);
            h_bkgtrk_corrpt[mycbin][ibin2]->Fill(jet_corrpt,n_bkgtrk[ibin2],pthat_weight*weight_cen*weight_vz); 

            if(n_trk[ibin2]>0) h_chg_corrpt[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
            //if(n_sube0_gen[ibin2]>0) h_chg_corrpt[mycbin][ibin2]->Fill(jet_corrpt,gen_chg_ptcut_sube0[ibin2],pthat_weight*weight_cen*weight_vz);        
            h_chg_corrpt_bkg[mycbin][ibin2]->Fill(jet_corrpt,trk_bkg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);

            if(n_trk[0]>0) h_chg_corrpt_kappa[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
            h_chg_corrpt_bkg_kappa[mycbin][ibin2]->Fill(jet_corrpt,trk_bkg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);
     
            if(n_trk[ibin2]>0) h_reco_corr_no0trks[mycbin][ibin2]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
            else if(n_trk[ibin2]==0) h_reco_corr_0trks[mycbin][ibin2]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);         

            if(isdata && do_eff_studies==1){    
              for(int ibin3=0;ibin3<trkeffbins;ibin3++){
                if(n_trk[ibin2]>0) h_chg_corrpt_trkeff[mycbin][ibin2][ibin3]->Fill(jet_corrpt,trk_chg_ptcut_trkeff[ibin2][ibin3],pthat_weight*weight_cen*weight_vz);
                if(n_trk[0]>0) h_chg_corrpt_kappa_trkeff[mycbin][ibin2][ibin3]->Fill(jet_corrpt,trk_chg_kappa_trkeff[ibin2][ibin3],pthat_weight*weight_cen*weight_vz);
              }
            }

            if(!isdata){

//ptcut
                if(n_gen[ibin2] > 0 && n_trk[ibin2] > 0) h_chg_recopt_genpt[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],gen_chg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_eff_sube0_gen[ibin2] > 0 && n_trk[ibin2] > 0) h_chg_recopt_recoptsube0[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],gen_chg_ptcut_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_eff_sube0_gen[ibin2] > 0 && n_sube0_gen[ibin2] > 0) h_chg_recoptsube0_genptsube0[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0_eff[ibin2],gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_trk[ibin2] > 0 && n_sube0_gen[ibin2] > 0) h_chg_recopt_genptsube0[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_eff_gen[ibin2] > 0 && n_sube0_gen[ibin2] > 0) h_chg_effgenpt_genptsube0[mycbin][ibin2]->Fill(gen_chg_ptcut_eff[ibin2],gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                if(calo_refpt->at(jet)>=120. && n_trk[ibin2] > 0 && n_sube0_gengen[ibin2] > 0) h_chg_recorecopt_gengensube0pt[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],gengen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);

                if(calo_refparton_flavor->at(jet)==2){
                    if(n_eff_sube0_gen[ibin2] > 0 && n_trk[ibin2] > 0) h_chg_recopt_recoptsube0_up[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],gen_chg_ptcut_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_sube0_gen[ibin2] > 0 && n_sube0_gen[ibin2] > 0) h_chg_recoptsube0_genptsube0_up[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0_eff[ibin2],gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_trk[ibin2] > 0 && n_sube0_gen[ibin2] > 0) h_chg_recopt_genptsube0_up[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_gen[ibin2] > 0 && n_sube0_gen[ibin2] > 0) h_chg_effgenpt_genptsube0_up[mycbin][ibin2]->Fill(gen_chg_ptcut_eff[ibin2],gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(calo_refpt->at(jet)>=120. && n_trk[ibin2] > 0 && n_sube0_gengen[ibin2] > 0) h_chg_recorecopt_gengensube0pt_up[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],gengen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);

                    if(n_gen[ibin2] > 0) h_chg_genpt_up[mycbin][ibin2]->Fill(gen_chg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_gen[ibin2] > 0) h_chg_eff_genpt_up[mycbin][ibin2]->Fill(gen_chg_ptcut_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_sube0_gen[ibin2] > 0) h_chg_genpt_sube0_up[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_sube0_gen[ibin2] > 0) h_chg_recopt_sube0_up[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_trk[ibin2] > 0) h_chg_recopt_up[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
                }
                else if(calo_refparton_flavor->at(jet)==1){
                    if(n_eff_sube0_gen[ibin2] > 0 && n_trk[ibin2] > 0) h_chg_recopt_recoptsube0_d[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],gen_chg_ptcut_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_sube0_gen[ibin2] > 0 && n_sube0_gen[ibin2] > 0) h_chg_recoptsube0_genptsube0_d[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0_eff[ibin2],gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_trk[ibin2] > 0 && n_sube0_gen[ibin2] > 0) h_chg_recopt_genptsube0_d[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_gen[ibin2] > 0 && n_sube0_gen[ibin2] > 0) h_chg_effgenpt_genptsube0_d[mycbin][ibin2]->Fill(gen_chg_ptcut_eff[ibin2],gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(calo_refpt->at(jet)>=120. && n_trk[ibin2] > 0 && n_sube0_gengen[ibin2] > 0) h_chg_recorecopt_gengensube0pt_d[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],gengen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);

                    if(n_gen[ibin2] > 0) h_chg_genpt_d[mycbin][ibin2]->Fill(gen_chg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_gen[ibin2] > 0) h_chg_eff_genpt_d[mycbin][ibin2]->Fill(gen_chg_ptcut_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_sube0_gen[ibin2] > 0) h_chg_genpt_sube0_d[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_sube0_gen[ibin2] > 0) h_chg_recopt_sube0_d[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_trk[ibin2] > 0) h_chg_recopt_d[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
                }
                else{
                    if(n_eff_sube0_gen[ibin2] > 0 && n_trk[ibin2] > 0) h_chg_recopt_recoptsube0_g_oth[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],gen_chg_ptcut_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_sube0_gen[ibin2] > 0 && n_sube0_gen[ibin2] > 0) h_chg_recoptsube0_genptsube0_g_oth[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0_eff[ibin2],gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_trk[ibin2] > 0 && n_sube0_gen[ibin2] > 0) h_chg_recopt_genptsube0_g_oth[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_gen[ibin2] > 0 && n_sube0_gen[ibin2] > 0) h_chg_effgenpt_genptsube0_g_oth[mycbin][ibin2]->Fill(gen_chg_ptcut_eff[ibin2],gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(calo_refpt->at(jet)>=120. && n_trk[ibin2] > 0 && n_sube0_gengen[ibin2] > 0) h_chg_recorecopt_gengensube0pt_g_oth[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],gengen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);

                    if(n_gen[ibin2] > 0) h_chg_genpt_g_oth[mycbin][ibin2]->Fill(gen_chg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_gen[ibin2] > 0) h_chg_eff_genpt_g_oth[mycbin][ibin2]->Fill(gen_chg_ptcut_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_sube0_gen[ibin2] > 0) h_chg_genpt_sube0_g_oth[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_sube0_gen[ibin2] > 0) h_chg_recopt_sube0_g_oth[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_trk[ibin2] > 0) h_chg_recopt_g_oth[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
                }                   

                if(n_gen[ibin2] > 0) h_chg_genpt[mycbin][ibin2]->Fill(gen_chg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_eff_gen[ibin2] > 0) h_chg_eff_genpt[mycbin][ibin2]->Fill(gen_chg_ptcut_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_sube0_gen[ibin2] > 0) h_chg_genpt_sube0[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_eff_sube0_gen[ibin2] > 0) h_chg_recopt_sube0[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);

                if(n_sube0_gen[ibin2] > 0) h_chg_genpt_sube0_bkgadd[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0[ibin2]/*+h_chg_genpt_bkg_prior[mycbin][1]->GetRandom()*/,weight_vz*weight_cen*pthat_weight);
                if(n_eff_sube0_gen[ibin2] > 0) h_chg_recopt_sube0_bkgadd[mycbin][ibin2]->Fill(gen_chg_ptcut_sube0_eff[ibin2]/*+h_chg_recopt_bkg_prior[mycbin][1]->GetRandom()*/,weight_vz*weight_cen*pthat_weight);
                h_chg_genpt_bkg[mycbin][ibin2]->Fill(gen_bkg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
                h_chg_recopt_bkg[mycbin][ibin2]->Fill(trk_bkg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);

//kappa
                if(n_gen[0] > 0 && n_trk[0] > 0) h_chg_recopt_genpt_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],gen_chg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_eff_sube0_gen[0] > 0 && n_trk[0] > 0) h_chg_recopt_recoptsube0_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],gen_chg_kappa_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_eff_sube0_gen[0] > 0 && n_sube0_gen[0] > 0) h_chg_recoptsube0_genptsube0_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0_eff[ibin2],gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_trk[0] > 0 && n_sube0_gen[0] > 0) h_chg_recopt_genptsube0_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_eff_gen[0] > 0 && n_sube0_gen[0] > 0) h_chg_effgenpt_genptsube0_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_eff[ibin2],gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                if(calo_refpt->at(jet)>=120. && n_trk[0] > 0 && n_sube0_gengen[0] > 0) h_chg_recorecopt_gengensube0pt_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],gengen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);

                if(calo_refparton_flavor->at(jet)==2){
                    if(n_eff_sube0_gen[0] > 0 && n_trk[0] > 0) h_chg_recopt_recoptsube0_up_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],gen_chg_kappa_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_sube0_gen[0] > 0 && n_sube0_gen[0] > 0) h_chg_recoptsube0_genptsube0_up_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0_eff[ibin2],gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_trk[0] > 0 && n_sube0_gen[0] > 0) h_chg_recopt_genptsube0_up_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_gen[0] > 0 && n_sube0_gen[0] > 0) h_chg_effgenpt_genptsube0_up_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_eff[ibin2],gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(calo_refpt->at(jet)>=120. && n_trk[0] > 0 && n_sube0_gengen[0] > 0) h_chg_recorecopt_gengensube0pt_up_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],gengen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);

                    if(n_gen[0] > 0) h_chg_genpt_up_kappa[mycbin][ibin2]->Fill(gen_chg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_gen[0] > 0) h_chg_eff_genpt_up_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_sube0_gen[0] > 0) h_chg_genpt_sube0_up_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_sube0_gen[0] > 0) h_chg_recopt_sube0_up_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_trk[0] > 0) h_chg_recopt_up_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
                }
                else if(calo_refparton_flavor->at(jet)==1){
                    if(n_eff_sube0_gen[0] > 0 && n_trk[0] > 0) h_chg_recopt_recoptsube0_d_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],gen_chg_kappa_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_sube0_gen[0] > 0 && n_sube0_gen[0] > 0) h_chg_recoptsube0_genptsube0_d_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0_eff[ibin2],gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_trk[0] > 0 && n_sube0_gen[0] > 0) h_chg_recopt_genptsube0_d_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_gen[0] > 0 && n_sube0_gen[0] > 0) h_chg_effgenpt_genptsube0_d_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_eff[ibin2],gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(calo_refpt->at(jet)>=120. && n_trk[0] > 0 && n_sube0_gengen[0] > 0) h_chg_recorecopt_gengensube0pt_d_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],gengen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);

                    if(n_gen[0] > 0) h_chg_genpt_d_kappa[mycbin][ibin2]->Fill(gen_chg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_gen[0] > 0) h_chg_eff_genpt_d_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_sube0_gen[0] > 0) h_chg_genpt_sube0_d_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_sube0_gen[0] > 0) h_chg_recopt_sube0_d_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_trk[0] > 0) h_chg_recopt_d_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
                }
                else{
                    if(n_eff_sube0_gen[0] > 0 && n_trk[0] > 0) h_chg_recopt_recoptsube0_g_oth_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],gen_chg_kappa_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_sube0_gen[0] > 0 && n_sube0_gen[0] > 0) h_chg_recoptsube0_genptsube0_g_oth_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0_eff[ibin2],gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_trk[0] > 0 && n_sube0_gen[0] > 0) h_chg_recopt_genptsube0_g_oth_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_gen[0] > 0 && n_sube0_gen[0] > 0) h_chg_effgenpt_genptsube0_g_oth_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_eff[ibin2],gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(calo_refpt->at(jet)>=120. && n_trk[0] > 0 && n_sube0_gengen[0] > 0) h_chg_recorecopt_gengensube0pt_g_oth_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],gengen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);

                    if(n_gen[0] > 0) h_chg_genpt_g_oth_kappa[mycbin][ibin2]->Fill(gen_chg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_gen[0] > 0) h_chg_eff_genpt_g_oth_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_sube0_gen[0] > 0) h_chg_genpt_sube0_g_oth_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_eff_sube0_gen[0] > 0) h_chg_recopt_sube0_g_oth_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                    if(n_trk[0] > 0) h_chg_recopt_g_oth_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
                }                   

                if(n_gen[0] > 0) h_chg_genpt_kappa[mycbin][ibin2]->Fill(gen_chg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_eff_gen[0] > 0) h_chg_eff_genpt_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_eff[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_sube0_gen[0] > 0) h_chg_genpt_sube0_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0[ibin2],weight_vz*weight_cen*pthat_weight);
                if(n_eff_sube0_gen[0] > 0) h_chg_recopt_sube0_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0_eff[ibin2],weight_vz*weight_cen*pthat_weight);

                if(n_sube0_gen[0] > 0) h_chg_genpt_sube0_bkgadd_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0[ibin2]/*+h_chg_genpt_bkg_prior_kappa[mycbin][1]->GetRandom()*/,weight_vz*weight_cen*pthat_weight);
                if(n_eff_sube0_gen[0] > 0) h_chg_recopt_sube0_bkgadd_kappa[mycbin][ibin2]->Fill(gen_chg_kappa_sube0_eff[ibin2]/*+h_chg_recopt_bkg_prior_kappa[mycbin][1]->GetRandom()*/,weight_vz*weight_cen*pthat_weight);
                //h_chg_genpt_bkg_kappa[mycbin][ibin2]->Fill(gen_bkg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
                //h_chg_recopt_bkg_kappa[mycbin][ibin2]->Fill(trk_bkg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
            }
            if(n_trk[ibin2] > 0) h_chg_recopt[mycbin][ibin2]->Fill(trk_chg_ptcut[ibin2],weight_vz*weight_cen*pthat_weight);
            if(n_trk[0] > 0) h_chg_recopt_kappa[mycbin][ibin2]->Fill(trk_chg_kappa[ibin2],weight_vz*weight_cen*pthat_weight);
          }      

          h_eta_full[mycbin]->Fill(calo_jteta->at(jet),pthat_weight*weight_cen*weight_vz);
          h_reco_corr[mycbin]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
        }

        if(!isdata){    
          if(calo_refparton_flavor->at(jet) == -999) continue;
          else if (fabs(calo_refparton_flavor->at(jet)) == 21){
            h_jt_gen_g_cloure[mycbin]->Fill(calo_refpt->at(jet),jet_corrpt/calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);

            if(calo_refpt->at(jet)>=100.){
              h_trk_refpt_g[mycbin]->Fill(calo_refpt->at(jet),n_gen[0],pthat_weight*weight_cen*weight_vz);

              for(int ibin2=0;ibin2<ntrkBins;ibin2++){
                if(n_gen[ibin2]>0) h_chg_refpt_g[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
                h_chg_refpt_bkg_g[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_bkg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
                if(n_trk[ibin2]>0) h_gen_full_no0trks_g[mycbin][ibin2]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
                else if(n_trk[ibin2]==0) h_gen_full_0trks_g[mycbin][ibin2]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);         
              }

              h_gen_full_g[mycbin]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
            }
            if(jet_corrpt>=100.){
              h_trk_corrpt_g[mycbin]->Fill(jet_corrpt,n_trk[0],pthat_weight*weight_cen*weight_vz);

              for(int ibin2=0;ibin2<ntrkBins;ibin2++){
                //if(n_trk[ibin2]>0) h_chg_corrpt_g[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
                if(n_sube0_gen[ibin2]>0) h_chg_corrpt_g[mycbin][ibin2]->Fill(jet_corrpt,gen_chg_ptcut_sube0[ibin2],pthat_weight*weight_cen*weight_vz);        
                h_chg_corrpt_bkg_g[mycbin][ibin2]->Fill(jet_corrpt,trk_bkg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);

                if(n_trk[0]>0) h_chg_corrpt_g_kappa[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
                h_chg_corrpt_bkg_g_kappa[mycbin][ibin2]->Fill(jet_corrpt,trk_bkg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);

                if(n_trk[ibin2]>0) h_reco_corr_no0trks_g[mycbin][ibin2]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
                else if(n_trk[ibin2]==0) h_reco_corr_0trks_g[mycbin][ibin2]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);         
              }

              h_reco_corr_g[mycbin]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
              h_eta_full_g[mycbin]->Fill(calo_jteta->at(jet),pthat_weight*weight_cen*weight_vz);
            }
          }   
          else {
            h_jt_gen_q_cloure[mycbin]->Fill(calo_refpt->at(jet),jet_corrpt/calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
            
            //gen
            if(calo_refpt->at(jet)>=100.){
              h_trk_refpt_q[mycbin]->Fill(calo_refpt->at(jet),n_gen[0],pthat_weight*weight_cen*weight_vz);

              for(int ibin2=0;ibin2<ntrkBins;ibin2++){
                if(n_gen[ibin2]>0) h_chg_refpt_q[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
                h_chg_refpt_bkg_q[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_bkg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
                if(n_trk[ibin2]>0) h_gen_full_no0trks_q[mycbin][ibin2]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
                else if(n_trk[ibin2]==0) h_gen_full_0trks_q[mycbin][ibin2]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);         
              }

              if(calo_refparton_flavor->at(jet)==1){
                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
                  if(n_gen[ibin2]>0) h_chg_refpt_d[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
                  if(n_trk[ibin2]==0) h_gen_full_0trks_down[mycbin][ibin2]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
                }
                h_gen_full_d[mycbin]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
              }
              else if(calo_refparton_flavor->at(jet)==-1){
                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
                  h_chg_refpt_dbar[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
                }
                h_gen_full_dbar[mycbin]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
              }
              else if(calo_refparton_flavor->at(jet)==2){
                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
                  if(n_gen[ibin2]>0) h_chg_refpt_up[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
                  if(n_trk[ibin2]==0) h_gen_full_0trks_up[mycbin][ibin2]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
                }
                h_gen_full_up[mycbin]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
              }
              else if(calo_refparton_flavor->at(jet)==-2){
                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
                  if(n_gen[ibin2]>0) h_chg_refpt_upbar[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
                }
                h_gen_full_upbar[mycbin]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
              }
              else if(fabs(calo_refparton_flavor->at(jet))==3){
                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
                  if(n_gen[ibin2]>0) h_chg_refpt_s[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
                }
                h_gen_full_s[mycbin]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
              }          
              else if(fabs(calo_refparton_flavor->at(jet))==4){       
                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
                  if(n_gen[ibin2]>0) h_chg_refpt_c[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
                }
                h_gen_full_c[mycbin]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
              }
              else if(fabs(calo_refparton_flavor->at(jet))==5){
                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
                  if(n_gen[ibin2]>0) h_chg_refpt_b[mycbin][ibin2]->Fill(calo_refpt->at(jet),gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
                }
                h_gen_full_b[mycbin]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);
              }

              h_gen_full_q[mycbin]->Fill(calo_refpt->at(jet),pthat_weight*weight_cen*weight_vz);

            }

            //reco
            if(jet_corrpt>=100.){
              h_trk_corrpt_q[mycbin]->Fill(jet_corrpt,n_trk[0],pthat_weight*weight_cen*weight_vz);

              for(int ibin2=0;ibin2<ntrkBins;ibin2++){
  //              if(n_trk[ibin2]>0) h_chg_corrpt_q[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
                if(n_sube0_gen[ibin2]>0) h_chg_corrpt_q[mycbin][ibin2]->Fill(jet_corrpt,gen_chg_ptcut_sube0[ibin2],pthat_weight*weight_cen*weight_vz);        
                h_chg_corrpt_bkg_q[mycbin][ibin2]->Fill(jet_corrpt,trk_bkg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);

                if(n_trk[0]>0) h_chg_corrpt_q_kappa[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
                h_chg_corrpt_bkg_q_kappa[mycbin][ibin2]->Fill(jet_corrpt,trk_bkg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);
                if(n_trk[ibin2]>0) h_reco_corr_no0trks_q[mycbin][ibin2]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);       
                else if(n_trk[ibin2]==0) h_reco_corr_0trks_q[mycbin][ibin2]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);   
              }

              if(calo_refparton_flavor->at(jet)==1){

                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
  //                if(n_trk[ibin2]>0) h_chg_corrpt_d[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
                  if(n_sube0_gen[ibin2]>0) h_chg_corrpt_d[mycbin][ibin2]->Fill(jet_corrpt,gen_chg_ptcut_sube0[ibin2],pthat_weight*weight_cen*weight_vz);        
                  if(n_trk[0]>0) h_chg_corrpt_d_kappa[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
                  if(n_trk[ibin2]==0) h_reco_corr_0trks_down[mycbin][ibin2]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
                }

                h_reco_corr_d[mycbin]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
                h_eta_full_d[mycbin]->Fill(calo_jteta->at(jet),pthat_weight*weight_cen*weight_vz);
              }
              else if(calo_refparton_flavor->at(jet)==-1){

                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
  //                if(n_trk[ibin2]>0) h_chg_corrpt_dbar[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
                  if(n_sube0_gen[ibin2]>0) h_chg_corrpt_dbar[mycbin][ibin2]->Fill(jet_corrpt,gen_chg_ptcut_sube0[ibin2],pthat_weight*weight_cen*weight_vz);        
                  if(n_trk[0]>0) h_chg_corrpt_dbar_kappa[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
                }

                h_reco_corr_dbar[mycbin]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
                h_eta_full_dbar[mycbin]->Fill(calo_jteta->at(jet),pthat_weight*weight_cen*weight_vz);
              }
              else if(calo_refparton_flavor->at(jet)==2){

                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
  //                if(n_trk[ibin2]>0) h_chg_corrpt_up[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
                  if(n_sube0_gen[ibin2]>0) h_chg_corrpt_up[mycbin][ibin2]->Fill(jet_corrpt,gen_chg_ptcut_sube0[ibin2],pthat_weight*weight_cen*weight_vz);        
                  if(n_trk[0]>0) h_chg_corrpt_up_kappa[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);  
                  if(n_trk[ibin2]==0) h_reco_corr_0trks_up[mycbin][ibin2]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);      
                }

                h_reco_corr_up[mycbin]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
                h_eta_full_up[mycbin]->Fill(calo_jteta->at(jet),pthat_weight*weight_cen*weight_vz);
              }
              else if(calo_refparton_flavor->at(jet)==-2){

                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
  //                if(n_trk[ibin2]>0) h_chg_corrpt_upbar[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
                  if(n_sube0_gen[ibin2]>0) h_chg_corrpt_upbar[mycbin][ibin2]->Fill(jet_corrpt,gen_chg_ptcut_sube0[ibin2],pthat_weight*weight_cen*weight_vz);        
                  if(n_trk[0]>0) h_chg_corrpt_upbar_kappa[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
                }

                h_reco_corr_upbar[mycbin]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
                h_eta_full_upbar[mycbin]->Fill(calo_jteta->at(jet),pthat_weight*weight_cen*weight_vz);
              }
              else if(fabs(calo_refparton_flavor->at(jet))==3){

                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
  //                if(n_trk[ibin2]>0) h_chg_corrpt_s[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
                  if(n_sube0_gen[ibin2]>0) h_chg_corrpt_s[mycbin][ibin2]->Fill(jet_corrpt,gen_chg_ptcut_sube0[ibin2],pthat_weight*weight_cen*weight_vz);        
                  if(n_trk[0]>0) h_chg_corrpt_s_kappa[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
                }

                h_reco_corr_s[mycbin]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
              }          
              else if(fabs(calo_refparton_flavor->at(jet))==4){

                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
  //                if(n_trk[ibin2]>0) h_chg_corrpt_c[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
                  if(n_sube0_gen[ibin2]>0) h_chg_corrpt_c[mycbin][ibin2]->Fill(jet_corrpt,gen_chg_ptcut_sube0[ibin2],pthat_weight*weight_cen*weight_vz);        
                  if(n_trk[0]>0) h_chg_corrpt_c_kappa[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
                }

                h_reco_corr_c[mycbin]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
              }
              else if(fabs(calo_refparton_flavor->at(jet))==5){

                for(int ibin2=0;ibin2<ntrkBins;ibin2++){
  //                if(n_trk[ibin2]>0) h_chg_corrpt_b[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
                  if(n_sube0_gen[ibin2]>0) h_chg_corrpt_b[mycbin][ibin2]->Fill(jet_corrpt,gen_chg_ptcut_sube0[ibin2],pthat_weight*weight_cen*weight_vz);        
                  if(n_trk[0]>0) h_chg_corrpt_b_kappa[mycbin][ibin2]->Fill(jet_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
                }

                h_reco_corr_b[mycbin]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
              }

              h_reco_corr_q[mycbin]->Fill(jet_corrpt,pthat_weight*weight_cen*weight_vz);
              h_eta_full_q[mycbin]->Fill(calo_jteta->at(jet),pthat_weight*weight_cen*weight_vz);
            }
          }
        }
        }//unfolding sample
      }//end of jet loop
      h_hibin_nrw->Fill(hiBin);
      h_hibin->Fill(hiBin,pthat_weight*weight_cen*weight_vz);
      h_vz[mycbin]->Fill(vz,pthat_weight*weight_cen*weight_vz);
      h_vz_nrw[mycbin]->Fill(vz,pthat_weight*weight_cen);
      h_pthat[mycbin]->Fill(pthat,pthat_weight*weight_cen*weight_vz);
      h_hibin_cbin->Fill(hiBin,mycbin,pthat_weight*weight_cen*weight_vz);

    }//event loop
  }//file loop

  TFile *closure_histos;

  if(isdata){
    if(ispp) closure_histos = new TFile(Form("ppdata_fulleta_jetchg_trkeffconst_no0trkjets_%d.root",date->GetDate()), "RECREATE");
    else closure_histos = new TFile(Form("PbPbdata_fulleta_jetchg_trkeffconst_no0trkjets_%d.root",date->GetDate()), "RECREATE");
  }
  else{
    if(ispp) closure_histos = new TFile(Form("Pythia6_fulleta_trkeffconst_split_jetchg_no0trkjets_%d.root",date->GetDate()), "RECREATE");
    else closure_histos = new TFile(Form("P+H_fulleta_trkeffconst_split_jetchg_no0trkjets_%d.root",date->GetDate()), "RECREATE");
  }  

  closure_histos->cd();

  for(int ibin=0;ibin<nCBins;ibin++){

    h_pthat[ibin]->Write();
    h_vz[ibin]->Write();
    h_vz_nrw[ibin]->Write();

    h_eta_full[ibin]->Write();
    h_reco_corr[ibin]->Write();

    h_chg_corrpt_highest[ibin]->Write();
    h_chg_corrpt_highest_2[ibin]->Write();

    h_ntrk[ibin]->Write();
    h_nbkgtrk[ibin]->Write();

    h_sampled_bkg[ibin]->Write();
    h_chg_refpt_subenon0[ibin]->Write();

   if(!isdata){ 
      h_trkpt_reco[ibin]->Write();
      h_trkpt_neg_reco[ibin]->Write();
      h_trkpt_pos_reco[ibin]->Write();

      h_trkpt_gen[ibin]->Write();
      h_trkpt_neg_gen[ibin]->Write();
      h_trkpt_pos_gen[ibin]->Write();  

      h_gen_full[ibin]->Write();
      h_gen_full_u[ibin]->Write();
      h_gen_full_q[ibin]->Write();
      h_gen_full_g[ibin]->Write();

      h_gen_full_c[ibin]->Write();
      h_gen_full_s[ibin]->Write();
      h_gen_full_b[ibin]->Write();

      h_gen_full_up[ibin]->Write();
      h_gen_full_upbar[ibin]->Write();
      h_gen_full_d[ibin]->Write();
      h_gen_full_dbar[ibin]->Write();

      h_eta_full_u[ibin]->Write();
      h_eta_full_q[ibin]->Write();
      h_eta_full_g[ibin]->Write();

      h_eta_full_up[ibin]->Write();
      h_eta_full_upbar[ibin]->Write();
      h_eta_full_d[ibin]->Write();
      h_eta_full_dbar[ibin]->Write();

      h_reco_corr_u[ibin]->Write();
      h_reco_corr_q[ibin]->Write();
      h_reco_corr_g[ibin]->Write();

      h_reco_corr_up[ibin]->Write();
      h_reco_corr_upbar[ibin]->Write();
      h_reco_corr_d[ibin]->Write();
      h_reco_corr_dbar[ibin]->Write();

      h_reco_corr_c[ibin]->Write();
      h_reco_corr_s[ibin]->Write();
      h_reco_corr_b[ibin]->Write();

      h_trk_refpt[ibin]->Write();
      h_trk_refpt_q[ibin]->Write();
      h_trk_refpt_g[ibin]->Write();
      h_trk_refpt_u[ibin]->Write();

      h_bkgtrk_refpt[ibin]->Write();
      h_sube0trk_refpt[ibin]->Write();

      h_trk_corrpt_q[ibin]->Write();
      h_trk_corrpt_g[ibin]->Write();
      h_trk_corrpt_u[ibin]->Write();

      h_jt_gen_cloure[ibin]->Write();
      h_jt_gen_q_cloure[ibin]->Write();
      h_jt_gen_g_cloure[ibin]->Write();

      for(int ibin2=0;ibin2<ntrkBins;ibin2++){
        h_chg_refpt[ibin][ibin2]->Write();
        h_chg_refpt_nobkgsub[ibin][ibin2]->Write();
        h_chg_refpt_eff[ibin][ibin2]->Write();
        h_chg_refpt_sube0[ibin][ibin2]->Write();
        h_chg_refpt_sube0_bkgadd[ibin][ibin2]->Write();
        h_chg_refpt_q[ibin][ibin2]->Write();
        h_chg_refpt_g[ibin][ibin2]->Write();
        h_chg_refpt_u[ibin][ibin2]->Write(); 

        h_chg_corrpt_q[ibin][ibin2]->Write();
        h_chg_corrpt_g[ibin][ibin2]->Write();
        h_chg_corrpt_u[ibin][ibin2]->Write();

        h_chg_corrpt_q_kappa[ibin][ibin2]->Write();
        h_chg_corrpt_g_kappa[ibin][ibin2]->Write();
        
        h_chg_refpt_bkg[ibin][ibin2]->Write();
        h_chg_refpt_bkg_q[ibin][ibin2]->Write();
        h_chg_refpt_bkg_g[ibin][ibin2]->Write();

        h_chg_corrpt_bkg_q[ibin][ibin2]->Write();
        h_chg_corrpt_bkg_g[ibin][ibin2]->Write();

        h_chg_corrpt_bkg_q_kappa[ibin][ibin2]->Write();
        h_chg_corrpt_bkg_g_kappa[ibin][ibin2]->Write();

        h_chg_refpt_up[ibin][ibin2]->Write();
        h_chg_refpt_upbar[ibin][ibin2]->Write();
        h_chg_refpt_d[ibin][ibin2]->Write();
        h_chg_refpt_dbar[ibin][ibin2]->Write();

        h_chg_refpt_c[ibin][ibin2]->Write();
        h_chg_refpt_s[ibin][ibin2]->Write();
        h_chg_refpt_b[ibin][ibin2]->Write();

        h_chg_corrpt_up[ibin][ibin2]->Write();
        h_chg_corrpt_upbar[ibin][ibin2]->Write();
        h_chg_corrpt_d[ibin][ibin2]->Write();
        h_chg_corrpt_dbar[ibin][ibin2]->Write();

        h_chg_corrpt_c[ibin][ibin2]->Write();
        h_chg_corrpt_s[ibin][ibin2]->Write();
        h_chg_corrpt_b[ibin][ibin2]->Write();

        h_chg_corrpt_up_kappa[ibin][ibin2]->Write();
        h_chg_corrpt_upbar_kappa[ibin][ibin2]->Write();
        h_chg_corrpt_d_kappa[ibin][ibin2]->Write();
        h_chg_corrpt_dbar_kappa[ibin][ibin2]->Write();

        h_chg_corrpt_c_kappa[ibin][ibin2]->Write();
        h_chg_corrpt_s_kappa[ibin][ibin2]->Write();
        h_chg_corrpt_b_kappa[ibin][ibin2]->Write();

        h_gen_full_no0trks[ibin][ibin2]->Write();
        h_gen_full_no0trks_q[ibin][ibin2]->Write();
        h_gen_full_no0trks_g[ibin][ibin2]->Write();

        h_reco_corr_no0trks_q[ibin][ibin2]->Write();
        h_reco_corr_no0trks_g[ibin][ibin2]->Write();

        h_reco_corr_0trks_q[ibin][ibin2]->Write();
        h_reco_corr_0trks_g[ibin][ibin2]->Write();
        h_reco_corr_0trks_up[ibin][ibin2]->Write();
        h_reco_corr_0trks_down[ibin][ibin2]->Write();

        h_gen_full_0trks[ibin][ibin2]->Write();
        h_gen_full_0trks_q[ibin][ibin2]->Write();
        h_gen_full_0trks_g[ibin][ibin2]->Write();
        h_gen_full_0trks_up[ibin][ibin2]->Write();
        h_gen_full_0trks_down[ibin][ibin2]->Write();

///unfolding histos
        h_chg_recopt_genpt[ibin][ibin2]->Write();
        h_chg_recoptsube0_genptsube0[ibin][ibin2]->Write();
        h_chg_recopt_recoptsube0[ibin][ibin2]->Write();
        h_chg_recopt_genptsube0[ibin][ibin2]->Write();
        h_chg_effgenpt_genptsube0[ibin][ibin2]->Write();
        h_chg_recorecopt_gengensube0pt[ibin][ibin2]->Write();

        h_chg_recoptsube0_genptsube0_up[ibin][ibin2]->Write();
        h_chg_recopt_recoptsube0_up[ibin][ibin2]->Write();
        h_chg_recopt_genptsube0_up[ibin][ibin2]->Write();
        h_chg_effgenpt_genptsube0_up[ibin][ibin2]->Write();
        h_chg_recorecopt_gengensube0pt_up[ibin][ibin2]->Write();

        h_chg_recoptsube0_genptsube0_d[ibin][ibin2]->Write();
        h_chg_recopt_recoptsube0_d[ibin][ibin2]->Write();
        h_chg_recopt_genptsube0_d[ibin][ibin2]->Write();
        h_chg_effgenpt_genptsube0_d[ibin][ibin2]->Write();
        h_chg_recorecopt_gengensube0pt_d[ibin][ibin2]->Write();

        h_chg_recoptsube0_genptsube0_g_oth[ibin][ibin2]->Write();
        h_chg_recopt_recoptsube0_g_oth[ibin][ibin2]->Write();
        h_chg_recopt_genptsube0_g_oth[ibin][ibin2]->Write();
        h_chg_effgenpt_genptsube0_g_oth[ibin][ibin2]->Write();
        h_chg_recorecopt_gengensube0pt_g_oth[ibin][ibin2]->Write();

        h_chg_genpt[ibin][ibin2]->Write();
        h_chg_genpt_bkg[ibin][ibin2]->Write();
        h_chg_eff_genpt[ibin][ibin2]->Write();
        h_chg_genpt_sube0[ibin][ibin2]->Write();
        h_chg_gengensube0pt[ibin][ibin2]->Write();
        h_chg_recopt_sube0[ibin][ibin2]->Write();
        h_chg_recopt_bkg[ibin][ibin2]->Write();

        h_chg_genpt_up[ibin][ibin2]->Write();
        h_chg_eff_genpt_up[ibin][ibin2]->Write();
        h_chg_genpt_sube0_up[ibin][ibin2]->Write();
        h_chg_recopt_sube0_up[ibin][ibin2]->Write();
        h_chg_recopt_up[ibin][ibin2]->Write();
        h_chg_gengensube0pt_up[ibin][ibin2]->Write();

        h_chg_recopt_datasample[ibin][ibin2]->Write();
        h_chg_recopt_d_datasample[ibin][ibin2]->Write();
        h_chg_recopt_up_datasample[ibin][ibin2]->Write();
        h_chg_recopt_g_oth_datasample[ibin][ibin2]->Write();

        h_chg_genpt_d[ibin][ibin2]->Write();
        h_chg_eff_genpt_d[ibin][ibin2]->Write();
        h_chg_genpt_sube0_d[ibin][ibin2]->Write();
        h_chg_recopt_sube0_d[ibin][ibin2]->Write();
        h_chg_recopt_d[ibin][ibin2]->Write();
        h_chg_gengensube0pt_d[ibin][ibin2]->Write();

        h_chg_genpt_g_oth[ibin][ibin2]->Write();
        h_chg_eff_genpt_g_oth[ibin][ibin2]->Write();
        h_chg_genpt_sube0_g_oth[ibin][ibin2]->Write();
        h_chg_recopt_sube0_g_oth[ibin][ibin2]->Write();
        h_chg_recopt_g_oth[ibin][ibin2]->Write();
        h_chg_gengensube0pt_g_oth[ibin][ibin2]->Write();

        h_chg_genpt_sube0_bkgadd[ibin][ibin2]->Write();
        h_chg_recopt_sube0_bkgadd[ibin][ibin2]->Write(); 

//kappa 
///unfolding histos
        h_chg_recopt_genpt_kappa[ibin][ibin2]->Write();
        h_chg_recoptsube0_genptsube0_kappa[ibin][ibin2]->Write();
        h_chg_recopt_recoptsube0_kappa[ibin][ibin2]->Write();
        h_chg_recopt_genptsube0_kappa[ibin][ibin2]->Write();
        h_chg_effgenpt_genptsube0_kappa[ibin][ibin2]->Write();
        h_chg_recorecopt_gengensube0pt_kappa[ibin][ibin2]->Write();

        h_chg_recoptsube0_genptsube0_up_kappa[ibin][ibin2]->Write();
        h_chg_recopt_recoptsube0_up_kappa[ibin][ibin2]->Write();
        h_chg_recopt_genptsube0_up_kappa[ibin][ibin2]->Write();
        h_chg_effgenpt_genptsube0_up_kappa[ibin][ibin2]->Write();
        h_chg_recorecopt_gengensube0pt_up_kappa[ibin][ibin2]->Write();

        h_chg_recoptsube0_genptsube0_d_kappa[ibin][ibin2]->Write();
        h_chg_recopt_recoptsube0_d_kappa[ibin][ibin2]->Write();
        h_chg_recopt_genptsube0_d_kappa[ibin][ibin2]->Write();
        h_chg_effgenpt_genptsube0_d_kappa[ibin][ibin2]->Write();
        h_chg_recorecopt_gengensube0pt_d_kappa[ibin][ibin2]->Write();

        h_chg_recoptsube0_genptsube0_g_oth_kappa[ibin][ibin2]->Write();
        h_chg_recopt_recoptsube0_g_oth_kappa[ibin][ibin2]->Write();
        h_chg_recopt_genptsube0_g_oth_kappa[ibin][ibin2]->Write();
        h_chg_effgenpt_genptsube0_g_oth_kappa[ibin][ibin2]->Write();
        h_chg_recorecopt_gengensube0pt_g_oth_kappa[ibin][ibin2]->Write();

        h_chg_genpt_kappa[ibin][ibin2]->Write();
        //h_chg_genpt_bkg_kappa[ibin][ibin2]->Write();
        h_chg_eff_genpt_kappa[ibin][ibin2]->Write();
        h_chg_genpt_sube0_kappa[ibin][ibin2]->Write();
        h_chg_gengensube0pt_kappa[ibin][ibin2]->Write();
        h_chg_recopt_sube0_kappa[ibin][ibin2]->Write();
        h_chg_recopt_kappa_datasample[ibin][ibin2]->Write();
        //h_chg_recopt_bkg_kappa[ibin][ibin2]->Write();

        h_chg_genpt_up_kappa[ibin][ibin2]->Write();
        h_chg_eff_genpt_up_kappa[ibin][ibin2]->Write();
        h_chg_genpt_sube0_up_kappa[ibin][ibin2]->Write();
        h_chg_recopt_sube0_up_kappa[ibin][ibin2]->Write();
        h_chg_recopt_up_kappa[ibin][ibin2]->Write();
        h_chg_recopt_up_kappa_datasample[ibin][ibin2]->Write();
        h_chg_gengensube0pt_up_kappa[ibin][ibin2]->Write();

        h_chg_genpt_d_kappa[ibin][ibin2]->Write();
        h_chg_eff_genpt_d_kappa[ibin][ibin2]->Write();
        h_chg_genpt_sube0_d_kappa[ibin][ibin2]->Write();
        h_chg_recopt_sube0_d_kappa[ibin][ibin2]->Write();
        h_chg_recopt_d_kappa[ibin][ibin2]->Write();
        h_chg_recopt_d_kappa_datasample[ibin][ibin2]->Write();
        h_chg_gengensube0pt_d_kappa[ibin][ibin2]->Write();

        h_chg_genpt_g_oth_kappa[ibin][ibin2]->Write();
        h_chg_eff_genpt_g_oth_kappa[ibin][ibin2]->Write();
        h_chg_genpt_sube0_g_oth_kappa[ibin][ibin2]->Write();
        h_chg_recopt_sube0_g_oth_kappa[ibin][ibin2]->Write();
        h_chg_recopt_g_oth_kappa[ibin][ibin2]->Write();
        h_chg_recopt_g_oth_kappa_datasample[ibin][ibin2]->Write();
        h_chg_gengensube0pt_g_oth_kappa[ibin][ibin2]->Write();

        h_chg_genpt_sube0_bkgadd_kappa[ibin][ibin2]->Write();
        h_chg_recopt_sube0_bkgadd_kappa[ibin][ibin2]->Write();            
      } 
    }

    for(int ibin2=0;ibin2<ntrkBins;ibin2++){
      h_reco_corr_no0trks[ibin][ibin2]->Write();
      h_reco_corr_0trks[ibin][ibin2]->Write();

      h_trk_corrpt[ibin][ibin2]->Write();
      h_bkgtrk_corrpt[ibin][ibin2]->Write();

      h_chg_corrpt[ibin][ibin2]->Write();
      h_chg_corrpt_nobkgsub[ibin][ibin2]->Write();
      h_chg_corrpt_bkg[ibin][ibin2]->Write();

      h_chg_recopt[ibin][ibin2]->Write();
      h_chg_recopt_kappa[ibin][ibin2]->Write();

      if(isdata && do_eff_studies){
        for(int ibin3=0; ibin3<trkeffbins;ibin3++){
          h_chg_corrpt_trkeff[ibin][ibin2][ibin3]->Write();
          h_chg_corrpt_kappa_trkeff[ibin][ibin2][ibin3]->Write();
        }
      }

      h_chg_corrpt_kappa[ibin][ibin2]->Write();
      h_chg_corrpt_bkg_kappa[ibin][ibin2]->Write();

      h_chg_corrpt_mix[ibin][ibin2]->Write();
    }
  }
  h_hibin_nrw->Write();
  h_hibin->Write();
  h_hibin_cbin->Write();

  closure_histos->Close();
}