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
#include <vector>
#include "TCanvas.h"
#include "skims/pp_chg_tree.h"
#include "skims/PbPb_chg_tree.h"
#include "jffcorr/nCScorr.h"

const int nCBins = 1;
const int nptBins = 55;
const int ntrkBins = 4;
const int nkbins = 4;
const int netaBins = 15;

using namespace std;

int mypbin, mycbin, myptbin, myrefptbin, myetabin;

char saythis[500];

TString cent[4] = {"0","1","2","3"};
TString trk[4] = {"0p7","2","4","5"};

TString pt[56] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55"};

int jt_nbins = 55;
Double_t jt_bin_bounds[56] = {50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280.,290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.,510.,520.,530.,540.,550.,560.,570.,580.,590.,600.};

Double_t eta_bounds[16] = {-1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5};
/*
///pp MC Pythia8 pthat weights
double pthatbins_Pythia8[9] = {50,80,100,120,170,220,280,370,9999};
double xsecs_Pythia8[9] = {3.778E-03, 4.412E-04, 1.511E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
double pthatEntries_Pythia8[8] ={175696, 145099, 160445, 258400, 189793, 196579, 54724, 12064};
*/

///pp MC Pythia6 pthat weights
double pthatbins[9] = {50,80,100,120,170,220,280,370,9999};
double xsecs[9] = {3.778E-03, 4.412E-04, 1.511E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
double pthatEntries[8] ={766451, 292063, 91200, 468748, 447938, 259208, 234447, 50972};

/*
///PbPb MC pthat weights
double pthatbins_P+H[9] = {15,30,50,80,120,170,220,280,9999};
double xsecs_P+H[9] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 0};
double pthatEntries_P+H[8] ={0, 0, 0, 2571566, 2.85082e+06, 2.68057e+06, 2.89138e+06, 950344};
*/
float CBins[5] = {0, 20, 60, 100, 200};

double vz, calo_jtpt, calo_corrpt, calo_refpt, calo_jteta, calo_jtphi, hiBin, pthat, pthat_weight, weight_vz, weight_cen, gen_chg, trk_chg, gen_chg_pt2, trk_chg_pt2, gen_chg_pt4, trk_chg_pt4, gen_chg_pt5, trk_chg_pt5;
double gen_chg_k0p3, gen_chg_k0p5, gen_chg_k0p7, gen_chg_k1, trk_chg_k0p3, trk_chg_k0p5, trk_chg_k0p7, trk_chg_k1;
double gen_bkg, trk_bkg, gen_bkg_pt2, trk_bkg_pt2, gen_bkg_pt4, trk_bkg_pt4, gen_bkg_pt5, trk_bkg_pt5;
int refparton_flavor, n_trk, n_gen, nPF_2, nCS_2, n_bkgtrk;

//cymbal tune centrality reweighting
double fcent_cymbal(double centrality, TF1* fcent1){ 
  return (centrality < 194) ? fcent1->Eval(centrality) : 1;
}

void jetchg(bool ispp=1, bool isdata=1){

  nCScorr *corrpt = new nCScorr(ispp,false);

//// defining histos and profiles

  TH1D *h_hibin = new TH1D("h_hibin","",200,0.,200.);
  TH1D *h_hibin_nrw = new TH1D("h_hibin_nrw","",200,0.,200.);
  TH1D *h_vz[nCBins];
  TH1D *h_vz_nrw[nCBins];

  TH1D *h_pthat[nCBins];

  TH1D *h_reco_full[nCBins];
  TH1D *h_reco_full_q[nCBins];
  TH1D *h_reco_full_g[nCBins];

  TH1D *h_reco_corr[nCBins];
  TH1D *h_reco_corr_q[nCBins];
  TH1D *h_reco_corr_g[nCBins];
  TH1D *h_reco_corr_u[nCBins];

  TH1D *h_reco_corr_up[nCBins];
  TH1D *h_reco_corr_upbar[nCBins];
  TH1D *h_reco_corr_d[nCBins];
  TH1D *h_reco_corr_dbar[nCBins];

  TH1D *h_gen_full[nCBins];
  TH1D *h_gen_full_q[nCBins];
  TH1D *h_gen_full_g[nCBins];  
  TH1D *h_gen_full_u[nCBins];

  TH2F *h_chg_refpt[nCBins][ntrkBins];
  TH2F *h_chg_refpt_q[nCBins][ntrkBins];  
  TH2F *h_chg_refpt_g[nCBins][ntrkBins];  
  TH2F *h_chg_refpt_u[nCBins][ntrkBins];  

  TH2F *h_chg_corrpt[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_q[nCBins][ntrkBins];  
  TH2F *h_chg_corrpt_g[nCBins][ntrkBins];  
  TH2F *h_chg_corrpt_u[nCBins][ntrkBins];

  TH2F *h_chg_refpt_bkg[nCBins][ntrkBins];
  TH2F *h_chg_refpt_bkg_q[nCBins][ntrkBins];
  TH2F *h_chg_refpt_bkg_g[nCBins][ntrkBins];

  TH2F *h_chg_corrpt_bkg[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_bkg_q[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_bkg_g[nCBins][ntrkBins];

  TH2F *h_chg_corrpt_up[nCBins][ntrkBins];  
  TH2F *h_chg_corrpt_upbar[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_d[nCBins][ntrkBins];
  TH2F *h_chg_corrpt_dbar[nCBins][ntrkBins];

  TH2F *h_trk_refpt[nCBins];
  TH2F *h_trk_refpt_q[nCBins];  
  TH2F *h_trk_refpt_g[nCBins];  
  TH2F *h_trk_refpt_u[nCBins];

  TH2F *h_trk_corrpt[nCBins];
  TH2F *h_bkgtrk_corrpt[nCBins];
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

  for (int ibin=0;ibin<nCBins;ibin++){

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
    h_gen_full[ibin] = new TH1D(saythis,"",jt_nbins-1,jt_bin_bounds);
    h_gen_full[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_q_cent%d",ibin);
    h_gen_full_q[ibin] = new TH1D(saythis,"",jt_nbins-1,jt_bin_bounds);
    h_gen_full_q[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_g_cent%d",ibin);
    h_gen_full_g[ibin] = new TH1D(saythis,"",jt_nbins-1,jt_bin_bounds);
    h_gen_full_g[ibin]->Sumw2();

    sprintf(saythis,"h_gen_full_u_cent%d",ibin);
    h_gen_full_u[ibin] = new TH1D(saythis,"",jt_nbins-1,jt_bin_bounds);
    h_gen_full_u[ibin]->Sumw2();

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
    h_reco_full[ibin] = new TH1D(saythis, "", jt_nbins-1,jt_bin_bounds);
    h_reco_full[ibin]->Sumw2();

    sprintf(saythis,"h_reco_full_q_cent%d",ibin);
    h_reco_full_q[ibin] = new TH1D(saythis,"",jt_nbins-1,jt_bin_bounds);
    h_reco_full_q[ibin]->Sumw2();

    sprintf(saythis,"h_reco_full_g_cent%d",ibin);
    h_reco_full_g[ibin] = new TH1D(saythis,"",jt_nbins-1,jt_bin_bounds);
    h_reco_full_g[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_cent%d",ibin);
    h_reco_corr[ibin] = new TH1D(saythis, "", jt_nbins-1,jt_bin_bounds);
    h_reco_corr[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_q_cent%d",ibin);
    h_reco_corr_q[ibin] = new TH1D(saythis,"",jt_nbins-1,jt_bin_bounds);
    h_reco_corr_q[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_g_cent%d",ibin);
    h_reco_corr_g[ibin] = new TH1D(saythis,"",jt_nbins-1,jt_bin_bounds);
    h_reco_corr_g[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_u_cent%d",ibin);
    h_reco_corr_u[ibin] = new TH1D(saythis,"",jt_nbins-1,jt_bin_bounds);
    h_reco_corr_u[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_up_cent%d",ibin);
    h_reco_corr_up[ibin] = new TH1D(saythis,"",jt_nbins-1,jt_bin_bounds);
    h_reco_corr_up[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_upbar_cent%d",ibin);
    h_reco_corr_upbar[ibin] = new TH1D(saythis,"",jt_nbins-1,jt_bin_bounds);
    h_reco_corr_upbar[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_d_cent%d",ibin);
    h_reco_corr_d[ibin] = new TH1D(saythis,"",jt_nbins-1,jt_bin_bounds);
    h_reco_corr_d[ibin]->Sumw2();

    sprintf(saythis,"h_reco_corr_dbar_cent%d",ibin);
    h_reco_corr_dbar[ibin] = new TH1D(saythis,"",jt_nbins-1,jt_bin_bounds);
    h_reco_corr_dbar[ibin]->Sumw2();

      sprintf(saythis,"h_trk_refpt_cent%d",ibin);
      h_trk_refpt[ibin] = new TH2F(saythis,"",600,0.,600.,40,0.,40.);
      h_trk_refpt[ibin]->Sumw2();

      sprintf(saythis,"h_trk_refpt_q_cent%d",ibin);
      h_trk_refpt_q[ibin] = new TH2F(saythis,"",600,0.,600.,40,0.,40.);
      h_trk_refpt_q[ibin]->Sumw2();

      sprintf(saythis,"h_trk_refpt_g_cent%d",ibin);
      h_trk_refpt_g[ibin] = new TH2F(saythis,"",600,0.,600.,40,0.,40.);
      h_trk_refpt_g[ibin]->Sumw2();

      sprintf(saythis,"h_trk_refpt_u_cent%d",ibin);
      h_trk_refpt_u[ibin] = new TH2F(saythis,"",600,0.,600.,40,0.,40.);
      h_trk_refpt_u[ibin]->Sumw2();

      sprintf(saythis,"h_trk_corrpt_cent%d",ibin);
      h_trk_corrpt[ibin] = new TH2F(saythis,"",600,0.,600.,40,0.,40.);
      h_trk_corrpt[ibin]->Sumw2();

      sprintf(saythis,"h_bkgtrk_corrpt_cent%d",ibin);
      h_bkgtrk_corrpt[ibin] = new TH2F(saythis,"",600,0.,600.,40,0.,40.);
      h_bkgtrk_corrpt[ibin]->Sumw2();

      sprintf(saythis,"h_trk_corrpt_q_cent%d",ibin);
      h_trk_corrpt_q[ibin] = new TH2F(saythis,"",600,0.,600.,40,0.,40.);
      h_trk_corrpt_q[ibin]->Sumw2();

      sprintf(saythis,"h_trk_corrpt_g_cent%d",ibin);
      h_trk_corrpt_g[ibin] = new TH2F(saythis,"",600,0.,600.,40,0.,40.);
      h_trk_corrpt_g[ibin]->Sumw2();

      sprintf(saythis,"h_trk_corrpt_u_cent%d",ibin);
      h_trk_corrpt_u[ibin] = new TH2F(saythis,"",600,0.,600.,40,0.,40.);
      h_trk_corrpt_u[ibin]->Sumw2();

    for (int ibin2=0;ibin2<ntrkBins;ibin2++){

      sprintf(saythis,"h_chg_refpt_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_refpt[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_q_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_q[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_refpt_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_g_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_g[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_refpt_g[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_bkg_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_bkg[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_refpt_bkg[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_bkg_q_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_bkg_q[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_refpt_bkg_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_bkg_g_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_bkg_g[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_refpt_bkg_g[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_refpt_u_cent%d_trk%d",ibin,ibin2);
      h_chg_refpt_u[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_refpt_u[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_corrpt[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_q_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_q[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_corrpt_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_g_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_g[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_corrpt_g[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_bkg_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_bkg[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_corrpt_bkg[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_bkg_q_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_bkg_q[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_corrpt_bkg_q[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_bkg_g_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_bkg_g[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_corrpt_bkg_g[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_u_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_u[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_corrpt_u[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_up_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_up[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_corrpt_up[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_upbar_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_upbar[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_corrpt_upbar[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_d_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_d[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_corrpt_d[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_chg_corrpt_dbar_cent%d_trk%d",ibin,ibin2);
      h_chg_corrpt_dbar[ibin][ibin2] = new TH2F(saythis,"",600,0.,600.,250,-2.5,2.5);
      h_chg_corrpt_dbar[ibin][ibin2]->Sumw2();
    }
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
   if(ispp) my_file = TFile::Open("skims/ppdata_bkg_jetcharge_Apr10.root");
   else my_file = TFile::Open("skims/PbPbdata_jetcharge_bkg_ptcut_eta0p5_1p6_Apr5.root");
  }
  else{
    if(ispp) my_file = TFile::Open("skims/Pythia6_bkg_jetcharge_Apr10.root"); 
    else my_file = TFile::Open("skims/P+H_jetcharge_bkg_ptcut_eta0p5_1p6_Apr5.root");
  } 

  cout<<"got file"<<endl;

  TTree *inp_tree = (TTree*)my_file->Get("unzipMixTree");
  pp_chg_tree *my_primary = new pp_chg_tree(inp_tree);
  //PbPb_chg_tree *my_primary = new PbPb_chg_tree(inp_tree);
  cout << "Successfully retrieved tree from input file!" << endl;
  Long64_t n_jets = my_primary->fChain->GetEntriesFast();

  //// Loop over all reco jets ////

  for (int jet = 0; jet < n_jets; jet++){ 
  //for (int jet = 0; jet < 20000; jet++){

    if (jet%100000==0) cout<<jet<<endl;

    my_primary->fChain->GetEntry(jet);

    calo_jteta = my_primary->jteta;  
    if(fabs(calo_jteta) >= etamaxcut || fabs(calo_jteta) <= etamincut) continue ;  

    calo_jtpt = my_primary->jtpt;
    //calo_corrpt = my_primary->corrpt;

    calo_jtphi = my_primary->jtphi;

    if(!isdata){
        calo_refpt = my_primary->refpt;
        if (calo_refpt <= refpTmincut) continue;

        refparton_flavor = my_primary->refparton_flavor;

        n_gen = my_primary->n_gen;

        gen_chg = my_primary->gen_totchg_k0p5;
        gen_chg_pt2 = my_primary->gen_totchg_pt2;
        gen_chg_pt4 = my_primary->gen_totchg_pt4;
        gen_chg_pt5 = my_primary->gen_totchg_pt5;

        gen_bkg = my_primary->gen_bkgchg_k0p5;
        gen_bkg_pt2 = my_primary->gen_bkgchg_pt2;
        gen_bkg_pt4 = my_primary->gen_bkgchg_pt4;
        gen_bkg_pt5 = my_primary->gen_bkgchg_pt5;

        gen_chg = gen_chg - gen_bkg;
        gen_chg_pt2 = gen_chg_pt2 - gen_bkg_pt2;
        gen_chg_pt4 = gen_chg_pt4 - gen_bkg_pt4;
        gen_chg_pt5 = gen_chg_pt5 - gen_bkg_pt5;
    }

    n_trk = my_primary->n_trk;

    //trk_chg = my_primary->trk_totchg;
    trk_chg = my_primary->trk_totchg_k0p5;
    trk_chg_pt2 = my_primary->trk_totchg_pt2;
    trk_chg_pt4 = my_primary->trk_totchg_pt4;
    trk_chg_pt5 = my_primary->trk_totchg_pt5;

    n_bkgtrk = my_primary->n_bkgtrk;

    trk_bkg = my_primary->trk_bkgchg_k0p5;
    trk_bkg_pt2 = my_primary->trk_bkgchg_pt2;
    trk_bkg_pt4 = my_primary->trk_bkgchg_pt4;
    trk_bkg_pt5 = my_primary->trk_bkgchg_pt5;

    trk_chg = trk_chg - trk_bkg;
    trk_chg_pt2 = trk_chg_pt2 - trk_bkg_pt2;
    trk_chg_pt4 = trk_chg_pt4 - trk_bkg_pt4;
    trk_chg_pt5 = trk_chg_pt5 - trk_bkg_pt5;

    //// centrality bin and weight 
    if(ispp) hiBin = 1;
    //else hiBin = my_primary->hiBin;

    double weight_cen = fcent_cymbal(hiBin,f_cent);
    if(ispp || isdata) weight_cen = 1.;

    nPF_2 = my_primary->nPFcand2_id1;
    if(nPF_2 > 20 && calo_jtpt <= 120.) calo_corrpt = corrpt->getCorrection(ispp, 20, hiBin, calo_jtpt, calo_jteta);
    else calo_corrpt = corrpt->getCorrection(ispp, nPF_2, hiBin, calo_jtpt, calo_jteta);

    //nCS_2 = my_primary->nCScand2_id145;
    //if(nCS_2 > 20 && calo_jtpt <= 120.) calo_corrpt = corrpt->getCorrection(ispp, 20, hiBin, calo_jtpt, calo_jteta);
    //else calo_corrpt = corrpt->getCorrection(ispp, nCS_2, hiBin, calo_jtpt, calo_jteta);
        
    if (calo_corrpt <= pTmincut || calo_corrpt >= pTmaxcut) continue;

    //// vz and weight
    vz = my_primary->vz;
    if (fabs(vz) > 15.) continue;

    //pythia8 
    //weight_vz = fvz_pythia8->Eval(vz);
    //pythia6 
    if(ispp) weight_vz = 1./fvz_pp->Eval(vz);
    //cymbal
    else weight_vz = fvz_cymbal->Eval(vz);

    if(isdata) weight_vz = 1.;  
  
    for (int cbin = 0; cbin < nCBins; cbin++){ 
      if (hiBin > CBins[cbin] && hiBin <= CBins[cbin+1]){
        mycbin = cbin; 
      }
    }
    if(ispp) mycbin = 0;
    
    ///// pthat weight

    pthat = my_primary->pthat;

    int ibin=0;
    while(pthat>pthatbins[ibin+1]) ibin++;
    pthat_weight = (xsecs[ibin]-xsecs[ibin+1])/pthatEntries[ibin];
    //pthat_weight = my_primary->weight;
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
      if(refparton_flavor == -999){
        h_trk_refpt_u[mycbin]->Fill(calo_refpt,n_gen,pthat_weight*weight_cen*weight_vz);
        h_chg_refpt_u[mycbin][0]->Fill(calo_refpt,gen_chg,pthat_weight*weight_cen*weight_vz);
        if(calo_refpt>=120.) h_gen_full_u[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
        if(calo_corrpt>=120.) {
          h_reco_corr_u[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
          h_eta_full_u[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
        }
        continue;
      }      
      else if(calo_refpt>=120.){
        h_trk_refpt[mycbin]->Fill(calo_refpt,n_gen,pthat_weight*weight_cen*weight_vz);

        h_chg_refpt[mycbin][0]->Fill(calo_refpt,gen_chg,pthat_weight*weight_cen*weight_vz);
        h_chg_refpt[mycbin][1]->Fill(calo_refpt,gen_chg_pt2,pthat_weight*weight_cen*weight_vz);
        h_chg_refpt[mycbin][2]->Fill(calo_refpt,gen_chg_pt4,pthat_weight*weight_cen*weight_vz);
        h_chg_refpt[mycbin][3]->Fill(calo_refpt,gen_chg_pt5,pthat_weight*weight_cen*weight_vz);

        h_chg_refpt_bkg[mycbin][0]->Fill(calo_refpt,gen_bkg,pthat_weight*weight_cen*weight_vz);
        h_chg_refpt_bkg[mycbin][1]->Fill(calo_refpt,gen_bkg_pt2,pthat_weight*weight_cen*weight_vz);
        h_chg_refpt_bkg[mycbin][2]->Fill(calo_refpt,gen_bkg_pt4,pthat_weight*weight_cen*weight_vz);
        h_chg_refpt_bkg[mycbin][3]->Fill(calo_refpt,gen_bkg_pt5,pthat_weight*weight_cen*weight_vz);

        h_gen_full[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
      }
    }  
    
    if(calo_corrpt>=120.){
      h_trk_corrpt[mycbin]->Fill(calo_corrpt,n_trk,pthat_weight*weight_cen*weight_vz);
      h_bkgtrk_corrpt[mycbin]->Fill(calo_corrpt,n_bkgtrk,pthat_weight*weight_cen*weight_vz);

      h_chg_corrpt[mycbin][0]->Fill(calo_corrpt,trk_chg,pthat_weight*weight_cen*weight_vz);
      h_chg_corrpt[mycbin][1]->Fill(calo_corrpt,trk_chg_pt2,pthat_weight*weight_cen*weight_vz);
      h_chg_corrpt[mycbin][2]->Fill(calo_corrpt,trk_chg_pt4,pthat_weight*weight_cen*weight_vz);
      h_chg_corrpt[mycbin][3]->Fill(calo_corrpt,trk_chg_pt5,pthat_weight*weight_cen*weight_vz);

      h_chg_corrpt_bkg[mycbin][0]->Fill(calo_corrpt,trk_bkg,pthat_weight*weight_cen*weight_vz);
      h_chg_corrpt_bkg[mycbin][1]->Fill(calo_corrpt,trk_bkg_pt2,pthat_weight*weight_cen*weight_vz);
      h_chg_corrpt_bkg[mycbin][2]->Fill(calo_corrpt,trk_bkg_pt4,pthat_weight*weight_cen*weight_vz);
      h_chg_corrpt_bkg[mycbin][3]->Fill(calo_corrpt,trk_bkg_pt5,pthat_weight*weight_cen*weight_vz);

      h_reco_corr[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
      h_eta_full[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
    }
    h_pthat[mycbin]->Fill(pthat,pthat_weight*weight_cen*weight_vz);
    h_vz[mycbin]->Fill(vz,pthat_weight*weight_cen*weight_vz);
    h_vz_nrw[mycbin]->Fill(vz,pthat_weight*weight_cen);
    h_hibin->Fill(hiBin,pthat_weight*weight_cen*weight_vz);

    if(!isdata){    
      if(refparton_flavor == -999) continue;
      else if (fabs(refparton_flavor) == 21){
        if(calo_refpt>=120.){
          h_trk_refpt_g[mycbin]->Fill(calo_refpt,n_gen,pthat_weight*weight_cen*weight_vz);

          h_chg_refpt_g[mycbin][0]->Fill(calo_refpt,gen_chg,pthat_weight*weight_cen*weight_vz);
          h_chg_refpt_g[mycbin][1]->Fill(calo_refpt,gen_chg_pt2,pthat_weight*weight_cen*weight_vz);          
          h_chg_refpt_g[mycbin][2]->Fill(calo_refpt,gen_chg_pt4,pthat_weight*weight_cen*weight_vz); 
          h_chg_refpt_g[mycbin][3]->Fill(calo_refpt,gen_chg_pt5,pthat_weight*weight_cen*weight_vz);                   

          h_chg_refpt_bkg_g[mycbin][0]->Fill(calo_refpt,gen_bkg,pthat_weight*weight_cen*weight_vz);
          h_chg_refpt_bkg_g[mycbin][1]->Fill(calo_refpt,gen_bkg_pt2,pthat_weight*weight_cen*weight_vz);
          h_chg_refpt_bkg_g[mycbin][2]->Fill(calo_refpt,gen_bkg_pt4,pthat_weight*weight_cen*weight_vz);
          h_chg_refpt_bkg_g[mycbin][3]->Fill(calo_refpt,gen_bkg_pt5,pthat_weight*weight_cen*weight_vz);

          h_gen_full_g[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);
        }
        if(calo_corrpt>=120.){
          h_trk_corrpt_g[mycbin]->Fill(calo_corrpt,n_trk,pthat_weight*weight_cen*weight_vz);
    
          h_chg_corrpt_g[mycbin][0]->Fill(calo_corrpt,trk_chg,pthat_weight*weight_cen*weight_vz);
          h_chg_corrpt_g[mycbin][1]->Fill(calo_corrpt,trk_chg_pt2,pthat_weight*weight_cen*weight_vz);
          h_chg_corrpt_g[mycbin][2]->Fill(calo_corrpt,trk_chg_pt4,pthat_weight*weight_cen*weight_vz);
          h_chg_corrpt_g[mycbin][3]->Fill(calo_corrpt,trk_chg_pt5,pthat_weight*weight_cen*weight_vz);

          h_chg_corrpt_bkg_g[mycbin][0]->Fill(calo_corrpt,trk_bkg,pthat_weight*weight_cen*weight_vz); 
          h_chg_corrpt_bkg_g[mycbin][1]->Fill(calo_corrpt,trk_bkg_pt2,pthat_weight*weight_cen*weight_vz);
          h_chg_corrpt_bkg_g[mycbin][2]->Fill(calo_corrpt,trk_bkg_pt4,pthat_weight*weight_cen*weight_vz);
          h_chg_corrpt_bkg_g[mycbin][3]->Fill(calo_corrpt,trk_bkg_pt5,pthat_weight*weight_cen*weight_vz);

          h_reco_corr_g[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
          h_eta_full_g[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
        }
      }   
      else {
        if(calo_refpt>=120.){
          h_trk_refpt_q[mycbin]->Fill(calo_refpt,n_gen,pthat_weight*weight_cen*weight_vz);

          h_chg_refpt_q[mycbin][0]->Fill(calo_refpt,gen_chg,pthat_weight*weight_cen*weight_vz);
          h_chg_refpt_q[mycbin][1]->Fill(calo_refpt,gen_chg_pt2,pthat_weight*weight_cen*weight_vz);
          h_chg_refpt_q[mycbin][2]->Fill(calo_refpt,gen_chg_pt4,pthat_weight*weight_cen*weight_vz);
          h_chg_refpt_q[mycbin][3]->Fill(calo_refpt,gen_chg_pt5,pthat_weight*weight_cen*weight_vz);

          h_chg_refpt_bkg_q[mycbin][0]->Fill(calo_refpt,gen_bkg,pthat_weight*weight_cen*weight_vz);
          h_chg_refpt_bkg_q[mycbin][1]->Fill(calo_refpt,gen_bkg_pt2,pthat_weight*weight_cen*weight_vz);
          h_chg_refpt_bkg_q[mycbin][2]->Fill(calo_refpt,gen_bkg_pt4,pthat_weight*weight_cen*weight_vz);
          h_chg_refpt_bkg_q[mycbin][3]->Fill(calo_refpt,gen_bkg_pt5,pthat_weight*weight_cen*weight_vz);

          h_gen_full_q[mycbin]->Fill(calo_refpt,pthat_weight*weight_cen*weight_vz);

        }
        if(calo_corrpt>=120.){
          h_trk_corrpt_q[mycbin]->Fill(calo_corrpt,n_trk,pthat_weight*weight_cen*weight_vz);

          h_chg_corrpt_q[mycbin][0]->Fill(calo_corrpt,trk_chg,pthat_weight*weight_cen*weight_vz);
          h_chg_corrpt_q[mycbin][1]->Fill(calo_corrpt,trk_chg_pt2,pthat_weight*weight_cen*weight_vz);
          h_chg_corrpt_q[mycbin][2]->Fill(calo_corrpt,trk_chg_pt4,pthat_weight*weight_cen*weight_vz);
          h_chg_corrpt_q[mycbin][3]->Fill(calo_corrpt,trk_chg_pt5,pthat_weight*weight_cen*weight_vz);

          h_chg_corrpt_bkg_q[mycbin][0]->Fill(calo_corrpt,trk_bkg,pthat_weight*weight_cen*weight_vz);
          h_chg_corrpt_bkg_q[mycbin][1]->Fill(calo_corrpt,trk_bkg_pt2,pthat_weight*weight_cen*weight_vz);
          h_chg_corrpt_bkg_q[mycbin][2]->Fill(calo_corrpt,trk_bkg_pt4,pthat_weight*weight_cen*weight_vz);
          h_chg_corrpt_bkg_q[mycbin][3]->Fill(calo_corrpt,trk_bkg_pt5,pthat_weight*weight_cen*weight_vz);

          if(refparton_flavor==1){
            h_chg_corrpt_d[mycbin][0]->Fill(calo_corrpt,trk_chg,pthat_weight*weight_cen*weight_vz);
            h_chg_corrpt_d[mycbin][1]->Fill(calo_corrpt,trk_chg_pt2,pthat_weight*weight_cen*weight_vz);
            h_chg_corrpt_d[mycbin][2]->Fill(calo_corrpt,trk_chg_pt4,pthat_weight*weight_cen*weight_vz); 
            h_chg_corrpt_d[mycbin][3]->Fill(calo_corrpt,trk_chg_pt5,pthat_weight*weight_cen*weight_vz);                    

            h_reco_corr_d[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
            h_eta_full_d[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
          }
          else if(refparton_flavor==-1){
            h_chg_corrpt_dbar[mycbin][0]->Fill(calo_corrpt,trk_chg,pthat_weight*weight_cen*weight_vz);
            h_chg_corrpt_dbar[mycbin][1]->Fill(calo_corrpt,trk_chg_pt2,pthat_weight*weight_cen*weight_vz);
            h_chg_corrpt_dbar[mycbin][2]->Fill(calo_corrpt,trk_chg_pt4,pthat_weight*weight_cen*weight_vz);
            h_chg_corrpt_dbar[mycbin][3]->Fill(calo_corrpt,trk_chg_pt5,pthat_weight*weight_cen*weight_vz);

            h_reco_corr_dbar[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
            h_eta_full_dbar[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
          }
          else if(refparton_flavor==2){
            h_chg_corrpt_up[mycbin][0]->Fill(calo_corrpt,trk_chg,pthat_weight*weight_cen*weight_vz);
            h_chg_corrpt_up[mycbin][1]->Fill(calo_corrpt,trk_chg_pt2,pthat_weight*weight_cen*weight_vz);
            h_chg_corrpt_up[mycbin][2]->Fill(calo_corrpt,trk_chg_pt4,pthat_weight*weight_cen*weight_vz);
            h_chg_corrpt_up[mycbin][3]->Fill(calo_corrpt,trk_chg_pt5,pthat_weight*weight_cen*weight_vz);

            h_reco_corr_up[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
            h_eta_full_up[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
          }
          else if(refparton_flavor==-2){
            h_chg_corrpt_upbar[mycbin][0]->Fill(calo_corrpt,trk_chg,pthat_weight*weight_cen*weight_vz);
            h_chg_corrpt_upbar[mycbin][1]->Fill(calo_corrpt,trk_chg_pt2,pthat_weight*weight_cen*weight_vz);
            h_chg_corrpt_upbar[mycbin][2]->Fill(calo_corrpt,trk_chg_pt4,pthat_weight*weight_cen*weight_vz); 
            h_chg_corrpt_upbar[mycbin][3]->Fill(calo_corrpt,trk_chg_pt5,pthat_weight*weight_cen*weight_vz);

            h_reco_corr_upbar[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
            h_eta_full_upbar[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
          }
          h_reco_corr_q[mycbin]->Fill(calo_corrpt,pthat_weight*weight_cen*weight_vz);
          h_eta_full_q[mycbin]->Fill(calo_jteta,pthat_weight*weight_cen*weight_vz);
        }
      }
    }

  }//end of jet loop

  TFile *closure_histos;

  if(isdata){
    if(ispp) closure_histos = new TFile("ppdata_jetchg_bkgsub_eta0p5_1p5_Apr26.root", "RECREATE");
    else closure_histos = new TFile("PbPbdata_jetchg_bkgsub_eta0p5_1p5_Apr26.root", "RECREATE");
  }
  else{
    if(ispp) closure_histos = new TFile("Pythia6_jetchg_bkgsub_eta0p5_1p5_Apr24.root", "RECREATE");
    else closure_histos = new TFile("P+H_jetchg_cymbal_bkgsub_eta0p5_1p5_Apr24.root", "RECREATE");
  }  

  closure_histos->cd();

  for(int ibin=0;ibin<nCBins;ibin++){

  if(!isdata){ 
    h_gen_full[ibin]->Write();
    h_gen_full_u[ibin]->Write();
    h_gen_full_q[ibin]->Write();
    h_gen_full_g[ibin]->Write();

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

    h_trk_refpt[ibin]->Write();
    h_trk_refpt_q[ibin]->Write();
    h_trk_refpt_g[ibin]->Write();
    h_trk_refpt_u[ibin]->Write();

    h_trk_corrpt_q[ibin]->Write();
    h_trk_corrpt_g[ibin]->Write();
    h_trk_corrpt_u[ibin]->Write();

    for(int ibin2=0;ibin2<ntrkBins;ibin2++){
      h_chg_refpt[ibin][ibin2]->Write();
      h_chg_refpt_q[ibin][ibin2]->Write();
      h_chg_refpt_g[ibin][ibin2]->Write();
      h_chg_refpt_u[ibin][ibin2]->Write(); 

      h_chg_corrpt_q[ibin][ibin2]->Write();
      h_chg_corrpt_g[ibin][ibin2]->Write();
      h_chg_corrpt_u[ibin][ibin2]->Write();

      h_chg_refpt_bkg[ibin][ibin2]->Write();
      h_chg_refpt_bkg_q[ibin][ibin2]->Write();
      h_chg_refpt_bkg_g[ibin][ibin2]->Write();

      h_chg_corrpt_bkg_q[ibin][ibin2]->Write();
      h_chg_corrpt_bkg_g[ibin][ibin2]->Write();

      h_chg_corrpt_up[ibin][ibin2]->Write();
      h_chg_corrpt_upbar[ibin][ibin2]->Write();
      h_chg_corrpt_d[ibin][ibin2]->Write();
      h_chg_corrpt_dbar[ibin][ibin2]->Write();
    } 
  }
    h_pthat[ibin]->Write();
    h_vz[ibin]->Write();
    h_vz_nrw[ibin]->Write();

    h_eta_full[ibin]->Write();
    h_reco_corr[ibin]->Write();
    h_trk_corrpt[ibin]->Write();
    h_bkgtrk_corrpt[ibin]->Write();

    for(int ibin2=0;ibin2<ntrkBins;ibin2++){
      h_chg_corrpt[ibin][ibin2]->Write();
      h_chg_corrpt_bkg[ibin][ibin2]->Write();
    }
  }
  h_hibin->Write();

  closure_histos->Close();
}