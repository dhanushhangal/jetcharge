#include <iostream>
#include "TFile.h"
#include "TRandom.h"
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
#include <TCut.h>
#include "TDatime.h"
#include <vector>
#include "TCanvas.h"
//#include "pp_jetchg_Oct19.h"
#include "PbPb_jetchg_Oct19.h"
#include "jffcorr/nCScorr.h"

const int nCBins = 4;
const int nptBins = 55;
const int ntrkBins = 4;
const int nkbins = 4;
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

double trkeffvary[trkeffbins] = {0.95,0.96,0.97,0.98,0.99,1.01,1.02,1.03,1.04,1.05};

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

float kappa[5] = {0.3, 0.5, 0.7, 1., 2.};

double data_MC_MB = 0.032416;
double bkg_scale[4] = {1-(0.077202947*data_MC_MB), 1-(0.10725959*data_MC_MB), 1-(0.13315407*data_MC_MB), 1-(0.17509690*data_MC_MB)};

//ptcut
double trk_eff_scale_pt2[5] = {0.99,1.,0.98,1.,1.};
double trk_eff_scale_pt4[5] = {0.99,0.99,1.,1.,1.};
double trk_eff_scale_pt5[5] = {0.99,0.99,1.,1.,1.};

//kappa
double trk_eff_scale_k0p3[5] = {0.99,1.,0.98,1.,1.};
double trk_eff_scale_k0p5[5] = {0.99,1.,0.98,1.,1.};
double trk_eff_scale_k0p7[5] = {0.975,0.975,0.98,0.97,1.};

/*
//ptcut
double trk_eff_scale_pt2[5] = {1.,1.,1.,1.,1.};
double trk_eff_scale_pt4[5] = {1.,1.,1.,1.,1.};
double trk_eff_scale_pt5[5] = {1.,1.,1.,1.,1.};

//kappa
double trk_eff_scale_k0p3[5] = {1.,1.,1.,1.,1.};
double trk_eff_scale_k0p5[5] = {1.,1.,1.,1.,1.};
double trk_eff_scale_k0p7[5] = {1.,1.,1.,1.,1.};
*/
double vz, calo_jtpt, calo_corrpt, calo_refpt, calo_jteta, calo_jtphi, calo_refeta, hiBin, pthat, pthat_weight, weight_vz, weight_cen; 
int refparton_flavor, nPF_2, nCS_2, n_sube0_gen, n_bkg_sube0_gen;
int n_trk[ntrkBins], n_gen[ntrkBins], n_bkggen[ntrkBins], n_bkgtrk[ntrkBins];
double trk_chg_ptcut_trkeff[ntrkBins][trkeffbins], trk_bkg_ptcut_trkeff[ntrkBins][trkeffbins], trk_chg_kappa_trkeff[ntrkBins][trkeffbins], trk_bkg_kappa_trkeff[ntrkBins][trkeffbins];
double trk_chg_ptcut[ntrkBins], trk_bkg_ptcut[ntrkBins], trk_chg_kappa[ntrkBins], trk_bkg_kappa[ntrkBins];
double gen_chg_ptcut[ntrkBins], gen_bkg_ptcut[ntrkBins];
double trk_chg_pt2_highest, trk_chg_pt2_highest_2;

//cymbal tune centrality reweighting
double fcent_cymbal(double centrality, TF1* fcent1){ 
  return (centrality < 194) ? fcent1->Eval(centrality) : 1;
}

//////////////jet resolution file///////////
TFile *f_res_param = TFile::Open("f_res_param_Dec1.root");
TF1 *f_res[nCBins];

TF1 *f_res_gauss = new TF1("f_res_gauss", "exp(-((x-1.)^2)/(2*([0]^2)))", 0., 2.);

void jetchg_kappa(bool ispp=1, bool isdata=1){

  nCScorr *corrpt = new nCScorr(ispp,true);

  f_res_gauss->FixParameter(0,0.05);

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
  TH2F *h_chg_refpt_q[nCBins][ntrkBins];  
  TH2F *h_chg_refpt_g[nCBins][ntrkBins];  
  TH2F *h_chg_refpt_u[nCBins][ntrkBins];  

  TH2F *h_chg_corrpt[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_q[nCBins][ntrkBins];  
  TH2F *h_chg_corrpt_g[nCBins][ntrkBins];  
  TH2F *h_chg_corrpt_u[nCBins][ntrkBins];

  TH2F *h_chg_corrpt_trkeff[nCBins][ntrkBins][trkeffbins];
  TH2F *h_chg_corrpt_kappa_trkeff[nCBins][ntrkBins][trkeffbins];

  TH2F *h_chg_corrpt_kappa[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_q_kappa[nCBins][ntrkBins];  
  TH2F *h_chg_corrpt_g_kappa[nCBins][ntrkBins];  
  TH2F *h_chg_corrpt_u_kappa[nCBins][ntrkBins];

  TH2F *h_chg_refpt_bkg[nCBins][ntrkBins];
  TH2F *h_chg_refpt_bkg_q[nCBins][ntrkBins];
  TH2F *h_chg_refpt_bkg_g[nCBins][ntrkBins];

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
    
  ///skims

  TFile *my_file;  
  if(isdata){
    if(ispp) my_file = TFile::Open("skims/preapp/ppdata_jetcharge_Oct19.root");
    //else my_file = TFile::Open("skims/preapp/PbPbdata_kappa_jetchg_Sep10.root");
    else my_file = TFile::Open("skims/preapp/PbPbdata_jetcharge_Oct19.root");
    //else my_file = TFile::Open("skims/preapp/PbPbdata_JER_jetchg_kappa_Sep14.root");
  }
  else{
    if(ispp) my_file = TFile::Open("skims/preapp/ppMC_jetcharge_Oct19.root"); 
    //else my_file = TFile::Open("skims/preapp/PbPbMC_kappa_jetchg_Sep10.root");
    else my_file = TFile::Open("skims/preapp/PbPbMC_jetcharge_Oct19.root");
  }

  cout<<"got file"<<endl;

  TTree *inp_tree = (TTree*)my_file->Get("unzipMixTree");

  //PbPbchg_Sep10 *my_primary = new PbPbchg_Sep10(inp_tree);
  PbPb_jetchg_Oct19 *my_primary = new PbPb_jetchg_Oct19(inp_tree);
  //pp_jetchg_Oct19 *my_primary = new pp_jetchg_Oct19(inp_tree);

  cout << "Successfully retrieved tree from input file!" << endl;
  Long64_t n_jets = my_primary->fChain->GetEntriesFast();
  cout<<"n jets: "<<n_jets<<endl;

  //// Loop over all reco jets ////

  for (int jet = 0; jet < n_jets; jet++){ 
  //for (int jet = 0; jet < 100; jet++){

    //// vz and weight
    vz = my_primary->vz;
    if (fabs(vz) > 15.) continue;

    if (jet%100000==0) cout<<jet<<endl;
    my_primary->fChain->GetEntry(jet);

    calo_jteta = my_primary->jteta;  
    if(fabs(calo_jteta) >= etamaxcut || fabs(calo_jteta) <= etamincut) continue ;   

    calo_jtpt = my_primary->jtpt;
    calo_corrpt = my_primary->corrpt;
    
    //JER smearing
    //calo_corrpt = calo_corrpt*(f_res_gauss->GetRandom());
    
    if (calo_corrpt <= pTmincut || calo_corrpt >= pTmaxcut) continue;

    calo_jtphi = my_primary->jtphi;

    if(!isdata){
        calo_refpt = my_primary->refpt;
        if (calo_refpt <= refpTmincut) continue;

        refparton_flavor = my_primary->refparton_flavor;

    }

    //// centrality bin and weight 
    if(ispp) hiBin = 1;
    //else hiBin = my_primary->hiBin;
    
    double weight_cen = fcent_cymbal(hiBin,f_cent);
    if(ispp || isdata) weight_cen = 1.;

    for (int cbin = 0; cbin < nCBins; cbin++){ 
      if (hiBin > CBins[cbin] && hiBin <= CBins[cbin+1]){
        mycbin = cbin; 
      }
    }

    n_trk[0] = my_primary->n_trkpt2; n_trk[1] = my_primary->n_trkpt3; n_trk[2] = my_primary->n_trkpt4; n_trk[3] = my_primary->n_trkpt5;
    n_bkgtrk[0] = my_primary->n_bkgtrk;
    //if(n_trk[0]==0) continue;

    trk_chg_kappa[0] = -999.; trk_chg_kappa[1] = -999.; trk_chg_kappa[2] = -999.; trk_chg_kappa[3] = -999.;

    //kappa
    if(n_trk[0] > 0) {
        trk_chg_kappa[0] = trk_eff_scale_k0p3[mycbin]*(1.001*my_primary->trk_totchg_pt2_k0p3_pos + my_primary->trk_totchg_pt2_k0p3_neg)/pow(calo_corrpt,0.3);    
        trk_chg_kappa[1] = trk_eff_scale_k0p5[mycbin]*(1.001*my_primary->trk_totchg_pt2_k0p5_pos + my_primary->trk_totchg_pt2_k0p5_neg)/pow(calo_corrpt,0.5);    
        trk_chg_kappa[2] = trk_eff_scale_k0p5[mycbin]*(1.001*my_primary->trk_totchg_pt2_k0p5_pos + my_primary->trk_totchg_pt2_k0p5_neg)/pow(calo_corrpt,0.5);    
        trk_chg_kappa[3] = trk_eff_scale_k0p7[mycbin]*(1.001*my_primary->trk_totchg_pt2_k0p7_pos + my_primary->trk_totchg_pt2_k0p7_neg)/pow(calo_corrpt,0.7);    
    }

    trk_bkg_kappa[0] = trk_eff_scale_k0p3[mycbin]*(1.001*my_primary->trk_bkgchg_pt2_k0p3_pos + my_primary->trk_bkgchg_pt2_k0p3_neg)/pow(calo_corrpt,0.3);    
    trk_bkg_kappa[1] = trk_eff_scale_k0p5[mycbin]*(1.001*my_primary->trk_bkgchg_pt2_k0p5_pos + my_primary->trk_bkgchg_pt2_k0p5_neg)/pow(calo_corrpt,0.5);    
    trk_bkg_kappa[2] = trk_eff_scale_k0p5[mycbin]*(1.001*my_primary->trk_bkgchg_pt2_k0p5_pos + my_primary->trk_bkgchg_pt2_k0p5_neg)/pow(calo_corrpt,0.5);    
    trk_bkg_kappa[3] = trk_eff_scale_k0p7[mycbin]*(1.001*my_primary->trk_bkgchg_pt2_k0p7_pos + my_primary->trk_bkgchg_pt2_k0p7_neg)/pow(calo_corrpt,0.7);    

    for(int ibin2=0; ibin2<ntrkBins; ibin2++){
      trk_chg_kappa[ibin2] = trk_chg_kappa[ibin2] - trk_bkg_kappa[ibin2];
    }

    for(int effbin=0; effbin<trkeffbins; effbin++){

        trk_chg_kappa_trkeff[0][effbin] = -999.; trk_chg_kappa_trkeff[1][effbin] = -999.; trk_chg_kappa_trkeff[2][effbin] = -999.; trk_chg_kappa_trkeff[3][effbin] = -999.;

        if(n_trk[0] > 0) {
            trk_chg_kappa_trkeff[0][effbin] = (trk_eff_scale_k0p3[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_totchg_pt2_k0p3_pos + my_primary->trk_totchg_pt2_k0p3_neg)/pow(calo_corrpt,0.3);    
            trk_chg_kappa_trkeff[1][effbin] = (trk_eff_scale_k0p5[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_totchg_pt2_k0p5_pos + my_primary->trk_totchg_pt2_k0p5_neg)/pow(calo_corrpt,0.5);    
            trk_chg_kappa_trkeff[2][effbin] = (trk_eff_scale_k0p5[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_totchg_pt2_k0p5_pos + my_primary->trk_totchg_pt2_k0p5_neg)/pow(calo_corrpt,0.5);    
            trk_chg_kappa_trkeff[3][effbin] = (trk_eff_scale_k0p7[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_totchg_pt2_k0p7_pos + my_primary->trk_totchg_pt2_k0p7_neg)/pow(calo_corrpt,0.7);    
        }

        trk_bkg_kappa_trkeff[0][effbin] = (trk_eff_scale_k0p3[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_bkgchg_pt2_k0p3_pos + my_primary->trk_bkgchg_pt2_k0p3_neg)/pow(calo_corrpt,0.3);    
        trk_bkg_kappa_trkeff[1][effbin] = (trk_eff_scale_k0p5[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_bkgchg_pt2_k0p5_pos + my_primary->trk_bkgchg_pt2_k0p5_neg)/pow(calo_corrpt,0.5);    
        trk_bkg_kappa_trkeff[2][effbin] = (trk_eff_scale_k0p5[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_bkgchg_pt2_k0p5_pos + my_primary->trk_bkgchg_pt2_k0p5_neg)/pow(calo_corrpt,0.5);    
        trk_bkg_kappa_trkeff[3][effbin] = (trk_eff_scale_k0p7[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_bkgchg_pt2_k0p7_pos + my_primary->trk_bkgchg_pt2_k0p7_neg)/pow(calo_corrpt,0.7); 

        for(int trkbin=0; trkbin<ntrkBins; trkbin++){
          trk_chg_kappa_trkeff[trkbin][effbin] = trk_chg_kappa_trkeff[trkbin][effbin] - trk_bkg_kappa_trkeff[trkbin][effbin];
        }
    }

    trk_chg_ptcut[0] = -999.; trk_chg_ptcut[1] = -999.; trk_chg_ptcut[2] = -999.; trk_chg_ptcut[3] = -999.;
    
    //ptcut
    if(n_trk[0] > 0) trk_chg_ptcut[0] = trk_eff_scale_pt2[mycbin]*(1.001*my_primary->trk_totchg_pt2_k0p5_pos + my_primary->trk_totchg_pt2_k0p5_neg)/pow(calo_corrpt,0.5);    
    if(n_trk[2] > 0) trk_chg_ptcut[1] = trk_eff_scale_pt4[mycbin]*(1.001*my_primary->trk_totchg_pt4_k0p5_pos + my_primary->trk_totchg_pt4_k0p5_neg)/pow(calo_corrpt,0.5);    
    if(n_trk[2] > 0) trk_chg_ptcut[2] = trk_eff_scale_pt4[mycbin]*(1.001*my_primary->trk_totchg_pt4_k0p5_pos + my_primary->trk_totchg_pt4_k0p5_neg)/pow(calo_corrpt,0.5);    
    if(n_trk[3] > 0) trk_chg_ptcut[3] = trk_eff_scale_pt5[mycbin]*(1.001*my_primary->trk_totchg_pt5_k0p5_pos + my_primary->trk_totchg_pt5_k0p5_neg)/pow(calo_corrpt,0.5);    

    trk_bkg_ptcut[0] = trk_eff_scale_pt2[mycbin]*(1.001*my_primary->trk_bkgchg_pt2_k0p5_pos + my_primary->trk_bkgchg_pt2_k0p5_neg)/pow(calo_corrpt,0.5);    
    trk_bkg_ptcut[1] = trk_eff_scale_pt4[mycbin]*(1.001*my_primary->trk_bkgchg_pt4_k0p5_pos + my_primary->trk_bkgchg_pt4_k0p5_neg)/pow(calo_corrpt,0.5);    
    trk_bkg_ptcut[2] = trk_eff_scale_pt4[mycbin]*(1.001*my_primary->trk_bkgchg_pt4_k0p5_pos + my_primary->trk_bkgchg_pt4_k0p5_neg)/pow(calo_corrpt,0.5);    
    trk_bkg_ptcut[3] = trk_eff_scale_pt5[mycbin]*(1.001*my_primary->trk_bkgchg_pt5_k0p5_pos + my_primary->trk_bkgchg_pt5_k0p5_neg)/pow(calo_corrpt,0.5);    

    for(int ibin2=0; ibin2<ntrkBins; ibin2++){
      trk_chg_ptcut[ibin2] = trk_chg_ptcut[ibin2] - trk_bkg_ptcut[ibin2];
    }

    for(int effbin=0; effbin<trkeffbins; effbin++){

        trk_chg_ptcut_trkeff[0][effbin] = -999.; trk_chg_ptcut_trkeff[1][effbin] = -999.; trk_chg_ptcut_trkeff[2][effbin] = -999.; trk_chg_ptcut_trkeff[3][effbin] = -999.;

        if(n_trk[0] > 0) trk_chg_ptcut_trkeff[0][effbin] = (trk_eff_scale_pt2[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_totchg_pt2_k0p5_pos + my_primary->trk_totchg_pt2_k0p5_neg)/pow(calo_corrpt,0.5);    
        if(n_trk[2] > 0) trk_chg_ptcut_trkeff[1][effbin] = (trk_eff_scale_pt4[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_totchg_pt4_k0p5_pos + my_primary->trk_totchg_pt4_k0p5_neg)/pow(calo_corrpt,0.5);    
        if(n_trk[2] > 0) trk_chg_ptcut_trkeff[2][effbin] = (trk_eff_scale_pt4[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_totchg_pt4_k0p5_pos + my_primary->trk_totchg_pt4_k0p5_neg)/pow(calo_corrpt,0.5);    
        if(n_trk[3] > 0) trk_chg_ptcut_trkeff[3][effbin] = (trk_eff_scale_pt5[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_totchg_pt5_k0p5_pos + my_primary->trk_totchg_pt5_k0p5_neg)/pow(calo_corrpt,0.5);    

        trk_bkg_ptcut_trkeff[0][effbin] = (trk_eff_scale_pt2[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_bkgchg_pt2_k0p5_pos + my_primary->trk_bkgchg_pt2_k0p5_neg)/pow(calo_corrpt,0.5);    
        trk_bkg_ptcut_trkeff[1][effbin] = (trk_eff_scale_pt4[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_bkgchg_pt4_k0p5_pos + my_primary->trk_bkgchg_pt4_k0p5_neg)/pow(calo_corrpt,0.5);    
        trk_bkg_ptcut_trkeff[2][effbin] = (trk_eff_scale_pt4[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_bkgchg_pt4_k0p5_pos + my_primary->trk_bkgchg_pt4_k0p5_neg)/pow(calo_corrpt,0.5);    
        trk_bkg_ptcut_trkeff[3][effbin] = (trk_eff_scale_pt5[mycbin]-0.05+(effbin*0.01))*(1.001*my_primary->trk_bkgchg_pt5_k0p5_pos + my_primary->trk_bkgchg_pt5_k0p5_neg)/pow(calo_corrpt,0.5); 

        for(int trkbin=0; trkbin<ntrkBins; trkbin++){
          trk_chg_ptcut_trkeff[trkbin][effbin] = trk_chg_ptcut_trkeff[trkbin][effbin] - trk_bkg_ptcut_trkeff[trkbin][effbin];
        }
    }
   
/*
    //nPF_2 = my_primary->nPFcand2_id1;
    nPF_2 = my_primary->npf_id1_pt2;
    if(nPF_2 > 20 && calo_jtpt <= 120.) calo_corrpt = corrpt->getCorrection(ispp, 20, hiBin, calo_jtpt, calo_jteta);
    else calo_corrpt = corrpt->getCorrection(ispp, nPF_2, hiBin, calo_jtpt, calo_jteta);

    //nCS_2 = my_primary->nCScand2_id145;
    nCS_2 = my_primary->ncs_id145_pt2;
    if(nCS_2 > 20 && calo_jtpt <= 120.) calo_corrpt = corrpt->getCorrection(ispp, 20, hiBin, calo_jtpt, calo_jteta);
    else calo_corrpt = corrpt->getCorrection(ispp, nCS_2, hiBin, calo_jtpt, calo_jteta);
*/
//    if (calo_corrpt <= pTmincut || calo_corrpt >= pTmaxcut) continue;

    //pythia8 
    //if(ispp) weight_vz = fvz_pythia8->Eval(vz);
    //pythia6 
    if(ispp) weight_vz = 1./fvz_pp->Eval(vz);
    //Pythia6+Hydjet
    else weight_vz = fvz_cymbal->Eval(vz);

    if(isdata) weight_vz = 1.;  

    ///// pthat weight
    pthat = my_primary->pthat;
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

    if(!isdata){
      for (int ptbin = 0; ptbin < nptBins; ptbin++){
        if (calo_refpt > jt_bin_bounds[ptbin] && calo_refpt <= jt_bin_bounds[ptbin+1]){       
          myrefptbin = ptbin; 
        }
      }
    }

    /////filling histos

    if(!isdata){      
      h_jt_gen_cloure[mycbin]->Fill(calo_refpt,calo_corrpt/calo_refpt,pthat_weight*weight_cen*weight_vz);
      if(refparton_flavor == -999){
        h_trk_refpt_u[mycbin]->Fill(calo_refpt,n_gen[0],pthat_weight*weight_cen*weight_vz);
        h_chg_refpt_u[mycbin][0]->Fill(calo_refpt,gen_chg_ptcut[0],pthat_weight*weight_cen*weight_vz);
        if(calo_refpt>=120.) h_gen_full_u[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
        if(calo_corrpt>=120.) {
          h_reco_corr_u[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
          h_eta_full_u[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
        }
        continue;
      }      
      if(calo_refpt>=100.){
        h_trk_refpt[mycbin]->Fill(calo_refpt,n_gen[0],pthat_weight*weight_cen*weight_vz);
        h_bkgtrk_refpt[mycbin]->Fill(calo_refpt,n_bkggen[0],pthat_weight*weight_cen*weight_vz);
        h_sube0trk_refpt[mycbin]->Fill(calo_refpt,n_gen[0]-n_sube0_gen,pthat_weight*weight_cen*weight_vz);

        for(int ibin2=0;ibin2<ntrkBins;ibin2++){
          h_chg_refpt[mycbin][ibin2]->Fill(calo_refpt,gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
          h_chg_refpt_bkg[mycbin][ibin2]->Fill(calo_refpt,gen_bkg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
          if(n_trk[ibin2]>0) h_gen_full_no0trks[mycbin][ibin2]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
          else if(n_trk[ibin2]==0) h_gen_full_0trks[mycbin][ibin2]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);         
        }
        h_gen_full[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
      }
    }  
    
    if(calo_corrpt>=100.){

      h_ntrk[mycbin]->Fill(n_trk[0],pthat_weight*weight_cen*weight_vz);
      h_nbkgtrk[mycbin]->Fill(n_bkgtrk[0],pthat_weight*weight_cen*weight_vz);

      for(int ibin2=0;ibin2<ntrkBins;ibin2++){
        h_trk_corrpt[mycbin][ibin2]->Fill(calo_corrpt,n_trk[ibin2],pthat_weight*weight_cen*weight_vz);
        h_bkgtrk_corrpt[mycbin][ibin2]->Fill(calo_corrpt,n_bkgtrk[ibin2],pthat_weight*weight_cen*weight_vz); 

        h_chg_corrpt[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
        h_chg_corrpt_bkg[mycbin][ibin2]->Fill(calo_corrpt,trk_bkg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);

        h_chg_corrpt_kappa[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
        h_chg_corrpt_bkg_kappa[mycbin][ibin2]->Fill(calo_corrpt,trk_bkg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);
 
        if(n_trk[ibin2]>0) h_reco_corr_no0trks[mycbin][ibin2]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
        else if(n_trk[ibin2]==0) h_reco_corr_0trks[mycbin][ibin2]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);         

        if(isdata){    
          for(int ibin3=0;ibin3<trkeffbins;ibin3++){
            h_chg_corrpt_trkeff[mycbin][ibin2][ibin3]->Fill(calo_corrpt,trk_chg_ptcut_trkeff[ibin2][ibin3],pthat_weight*weight_cen*weight_vz);
            h_chg_corrpt_kappa_trkeff[mycbin][ibin2][ibin3]->Fill(calo_corrpt,trk_chg_kappa_trkeff[ibin2][ibin3],pthat_weight*weight_cen*weight_vz);
          }
        }
      }      

      h_eta_full[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
      h_reco_corr[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);

      h_chg_corrpt_highest[mycbin]->Fill(calo_corrpt,trk_chg_pt2_highest,pthat_weight*weight_cen*weight_vz);
      h_chg_corrpt_highest_2[mycbin]->Fill(calo_corrpt,trk_chg_pt2_highest_2,pthat_weight*weight_cen*weight_vz);

      h_hibin_nrw->Fill(hiBin);
      h_hibin->Fill(hiBin,pthat_weight*weight_cen*weight_vz);
    }

    h_pthat[mycbin]->Fill(pthat,pthat_weight*weight_cen*weight_vz);
    h_vz[mycbin]->Fill(vz,pthat_weight*weight_cen*weight_vz);
    h_vz_nrw[mycbin]->Fill(vz,pthat_weight*weight_cen);
    h_hibin_cbin->Fill(hiBin,mycbin,pthat_weight*weight_cen*weight_vz);

    if(!isdata){    
      if(refparton_flavor == -999) continue;
      else if (fabs(refparton_flavor) == 21){
        h_jt_gen_g_cloure[mycbin]->Fill(calo_refpt,calo_corrpt/calo_refpt,pthat_weight*weight_cen*weight_vz);

        if(calo_refpt>=100.){
          h_trk_refpt_g[mycbin]->Fill(calo_refpt,n_gen[0],pthat_weight*weight_cen*weight_vz);

          for(int ibin2=0;ibin2<ntrkBins;ibin2++){
            h_chg_refpt_g[mycbin][ibin2]->Fill(calo_refpt,gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
            h_chg_refpt_bkg_g[mycbin][ibin2]->Fill(calo_refpt,gen_bkg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
            if(n_trk[ibin2]>0) h_gen_full_no0trks_g[mycbin][ibin2]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
            else if(n_trk[ibin2]==0) h_gen_full_0trks_g[mycbin][ibin2]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);         
          }

          h_gen_full_g[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
        }
        if(calo_corrpt>=100.){
          h_trk_corrpt_g[mycbin]->Fill(calo_corrpt,n_trk[0],pthat_weight*weight_cen*weight_vz);

          for(int ibin2=0;ibin2<ntrkBins;ibin2++){
            h_chg_corrpt_g[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
            h_chg_corrpt_bkg_g[mycbin][ibin2]->Fill(calo_corrpt,trk_bkg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);

            h_chg_corrpt_g_kappa[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
            h_chg_corrpt_bkg_g_kappa[mycbin][ibin2]->Fill(calo_corrpt,trk_bkg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);

            if(n_trk[ibin2]>0) h_reco_corr_no0trks_g[mycbin][ibin2]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
            else if(n_trk[ibin2]==0) h_reco_corr_0trks_g[mycbin][ibin2]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);         
          }

          h_reco_corr_g[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
          h_eta_full_g[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
        }
      }   
      else {
        h_jt_gen_q_cloure[mycbin]->Fill(calo_refpt,calo_corrpt/calo_refpt,pthat_weight*weight_cen*weight_vz);
        
        //gen
        if(calo_refpt>=100.){
          h_trk_refpt_q[mycbin]->Fill(calo_refpt,n_gen[0],pthat_weight*weight_cen*weight_vz);

          for(int ibin2=0;ibin2<ntrkBins;ibin2++){
            h_chg_refpt_q[mycbin][ibin2]->Fill(calo_refpt,gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
            h_chg_refpt_bkg_q[mycbin][ibin2]->Fill(calo_refpt,gen_bkg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
            if(n_trk[ibin2]>0) h_gen_full_no0trks_q[mycbin][ibin2]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
            else if(n_trk[ibin2]==0) h_gen_full_0trks_q[mycbin][ibin2]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);         
          }

          if(refparton_flavor==1){
            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_refpt_d[mycbin][ibin2]->Fill(calo_refpt,gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
              if(n_trk[ibin2]==0) h_gen_full_0trks_down[mycbin][ibin2]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
            }
            h_gen_full_d[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
          }
          else if(refparton_flavor==-1){
            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_refpt_dbar[mycbin][ibin2]->Fill(calo_refpt,gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
            }
            h_gen_full_dbar[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
          }
          else if(refparton_flavor==2){
            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_refpt_up[mycbin][ibin2]->Fill(calo_refpt,gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
              if(n_trk[ibin2]==0) h_gen_full_0trks_up[mycbin][ibin2]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
            }
            h_gen_full_up[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
          }
          else if(refparton_flavor==-2){
            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_refpt_upbar[mycbin][ibin2]->Fill(calo_refpt,gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
            }
            h_gen_full_upbar[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
          }
          else if(fabs(refparton_flavor)==3){
            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_refpt_s[mycbin][ibin2]->Fill(calo_refpt,gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
            }
            h_gen_full_s[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
          }          
          else if(fabs(refparton_flavor)==4){       
            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_refpt_c[mycbin][ibin2]->Fill(calo_refpt,gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
            }
            h_gen_full_c[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
          }
          else if(fabs(refparton_flavor)==5){
            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_refpt_b[mycbin][ibin2]->Fill(calo_refpt,gen_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
            }
            h_gen_full_b[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
          }

          h_gen_full_q[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);

        }

        //reco
        if(calo_corrpt>=100.){
          h_trk_corrpt_q[mycbin]->Fill(calo_corrpt,n_trk[0],pthat_weight*weight_cen*weight_vz);

          for(int ibin2=0;ibin2<ntrkBins;ibin2++){
            h_chg_corrpt_q[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
            h_chg_corrpt_bkg_q[mycbin][ibin2]->Fill(calo_corrpt,trk_bkg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);

            h_chg_corrpt_q_kappa[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
            h_chg_corrpt_bkg_q_kappa[mycbin][ibin2]->Fill(calo_corrpt,trk_bkg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);
            if(n_trk[ibin2]>0) h_reco_corr_no0trks_q[mycbin][ibin2]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);       
            else if(n_trk[ibin2]==0) h_reco_corr_0trks_q[mycbin][ibin2]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);   
          }

          if(refparton_flavor==1){

            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_corrpt_d[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);        
              h_chg_corrpt_d_kappa[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
              if(n_trk[ibin2]==0) h_reco_corr_0trks_down[mycbin][ibin2]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
            }

            h_reco_corr_d[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
            h_eta_full_d[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
          }
          else if(refparton_flavor==-1){

            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_corrpt_dbar[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
              h_chg_corrpt_dbar_kappa[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
            }

            h_reco_corr_dbar[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
            h_eta_full_dbar[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
          }
          else if(refparton_flavor==2){

            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_corrpt_up[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
              h_chg_corrpt_up_kappa[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);  
              if(n_trk[ibin2]==0) h_reco_corr_0trks_up[mycbin][ibin2]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);      
            }

            h_reco_corr_up[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
            h_eta_full_up[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
          }
          else if(refparton_flavor==-2){

            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_corrpt_upbar[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
              h_chg_corrpt_upbar_kappa[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
            }

            h_reco_corr_upbar[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
            h_eta_full_upbar[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
          }
          else if(fabs(refparton_flavor)==3){

            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_corrpt_s[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
              h_chg_corrpt_s_kappa[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
            }

            h_reco_corr_s[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
          }          
          else if(fabs(refparton_flavor)==4){

            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_corrpt_c[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
              h_chg_corrpt_c_kappa[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
            }

            h_reco_corr_c[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
          }
          else if(fabs(refparton_flavor)==5){

            for(int ibin2=0;ibin2<ntrkBins;ibin2++){
              h_chg_corrpt_b[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_ptcut[ibin2],pthat_weight*weight_cen*weight_vz);
              h_chg_corrpt_b_kappa[mycbin][ibin2]->Fill(calo_corrpt,trk_chg_kappa[ibin2],pthat_weight*weight_cen*weight_vz);        
            }

            h_reco_corr_b[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
          }

          h_reco_corr_q[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
          h_eta_full_q[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
        }
      }
    }

  }//end of jet loop

  TFile *closure_histos;

  if(isdata){
    if(ispp) closure_histos = new TFile(Form("ppdata_jetchg_posneg_trkeffvary_no0trkjets_%d.root",date->GetDate()), "RECREATE");
    else closure_histos = new TFile(Form("PbPbdata_jetchg_posneg_trkeffvary_no0trkjets_%d.root",date->GetDate()), "RECREATE");
  }
  else{
    if(ispp) closure_histos = new TFile(Form("Pythia6_jetchg_no0trkjets_%d.root",date->GetDate()), "RECREATE");
    else closure_histos = new TFile(Form("P+H_jetchg_trkeffconst_no0trkjets_%d.root",date->GetDate()), "RECREATE");
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

    if(!isdata){ 
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
      } 
    }

    for(int ibin2=0;ibin2<ntrkBins;ibin2++){
      h_reco_corr_no0trks[ibin][ibin2]->Write();
      h_reco_corr_0trks[ibin][ibin2]->Write();

      h_trk_corrpt[ibin][ibin2]->Write();
      h_bkgtrk_corrpt[ibin][ibin2]->Write();

      h_chg_corrpt[ibin][ibin2]->Write();
      h_chg_corrpt_bkg[ibin][ibin2]->Write();

      if(isdata){
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