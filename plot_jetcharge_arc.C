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
#include "fitting_templates_arc.h"

#define nCBins 5
#define nptBins 3
#define ntrkbins 3
#define nkbins 3
#define njtptbins 5
#define trkeffbins 11

const Double_t pt_low = 120.;
const Double_t pt_high = 600.;

bool do_fit = true;
bool do_systematics = true;
bool do_ptcut_fitting = true;
bool do_kappa_fitting = true;
bool do_mean_sigma = false;
bool do_pyquen = false;
bool do_eff_studies = true;

char saythis[500];

using namespace std;

Double_t fit_limit_ptcut[ntrkbins] = {1.,1.,1.};
Double_t fit_limit_kappa[ntrkbins] = {1.5,1.,0.8};

TString cent[5] = {"0","1","2","3","4"};
TString trk[3] = {"0","1","2"};
TString kbin[3] = {"0","1","2"};
TString pt[56] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55"};

TString trkeff[11] = {"0","1","2","3","4","5","6","7","8","9","10"};
TString trkeffvary[trkeffbins] = {"trk eff -5%","trk eff -4%","trk eff -3%","trk eff -2%","trk eff -1%","trk eff nom.","trk eff +1%","trk eff +2%","trk eff +3%","trk eff +4%","trk eff +5%"};

TString trk_pt[4] = {"","2","4","5"};
TString kappa_str[4] = {"1","1p5","2","2p5"};

TString hi_mc_tag[] = {"PYTHIA+HYDJET","","PYQUEN (Coll)","","PYQUEN (Wide)"};

TString cent_tag_pp_data[] = {"pp ref","0-10% PbPb","10-30% PbPb","30-50% PbPb","50-100% PbPb"};
//TString cent_tag_pp_data[] = {"pp ref","0-10%","10-30%","30-50%","50-100%"};
TString cent_tag_pp_MC[] = {"PYTHIA ref","0-10% P+H","10-30% P+H","30-50% P+H","50-100% P+H"};
TString x_label[] = {"x1","x2","x3"};
TString x_label_kappa[] = {"pT > 2 (#kappa=0.5); #kappa=0.5 (pT > 2)","pT > 4 (#kappa=0.5); #kappa=1 (pT > 2)","pT > 5 (#kappa=0.5); #kappa=2 (pT > 2)"};
TString cent_tag[] = {"0-10%","10-30%","30-50%","50-100%","70-100%"};
TString trk_tag[] = {"p_{T}^{trk} > 2 GeV, #kappa = 0.5","p_{T}^{trk} > 4 GeV, #kappa = 0.5", "p_{T}^{trk} > 5 GeV, #kappa = 0.5"};
TString kappa_tag[] = {"#kappa = 0.3 , p_{T}^{trk} > 2 GeV","#kappa = 0.5 , p_{T}^{trk} > 2 GeV", "#kappa = 0.7 , p_{T}^{trk} > 2 GeV"};

Double_t jtpt_bound[njtptbins+1] = {120.,130.,140.,160.,190.,500.};
Double_t jtpt_bound_err[njtptbins+1] = {0.,0.,0.,0.,0.,0.};
Double_t x_mean[ntrkbins] = {2.,4.,5.};
Double_t x_mean_err[ntrkbins] = {0.,0.,0.};
Double_t kappa[ntrkbins] = {0.3,0.5,0.7};
Double_t kappa_err[ntrkbins] = {0.,0.,0.};
Double_t cent_mean[nCBins] = {1.,2.,3.,4.,5.};
Double_t cent_mean_err[nCBins] = {0.1,0.1,0.1,0.1,0.1};

Double_t trkpt_bounds[ntrkbins+1] = {1.,3.,4.5,6.};
Double_t kappa_bounds[ntrkbins+1] = {0.2,0.4,0.6,0.8};

Double_t diff_0trk_data_MC[nCBins][ntrkbins];
Double_t gluon[nCBins],quark[nCBins],up[nCBins],upbar[nCBins],down[nCBins],downbar[nCBins],gluon_ubar_dbar[nCBins], charm[nCBins], strange[nCBins], bottom[nCBins];
Double_t gluon_gen[nCBins],quark_gen[nCBins], down_gen[nCBins], up_gen[nCBins], gluon_ubar_dbar_gen[nCBins];
Double_t ubdbcsb[nCBins],gluon_ubdbcsb[nCBins], u_ud[nCBins], d_ud[nCBins];
Double_t gluon_corr_factor[nCBins][ntrkbins], quark_corr_factor[nCBins][ntrkbins];
Double_t gluon_gen_int[nCBins][ntrkbins], quark_gen_int[nCBins][ntrkbins];
Double_t gluon_reco_int[nCBins][ntrkbins], quark_reco_int[nCBins][ntrkbins];
Double_t gluon_no0trk_int[nCBins][ntrkbins], quark_no0trk_int[nCBins][ntrkbins], quark_0trk_int[nCBins][ntrkbins], gluon_0trk_int[nCBins][ntrkbins]; Double_t gluon_no0trk_gen_int[nCBins][ntrkbins], quark_no0trk_gen_int[nCBins][ntrkbins];
Double_t quark_0trk_int_gen[nCBins][ntrkbins], quark_0trk_int_gr[nCBins][ntrkbins];
Double_t data_MC_0trk_factor[nCBins][ntrkbins], data_int[nCBins], quark_data_0trk_int[nCBins][ntrkbins], quark_mc_0trk_int[nCBins][ntrkbins], gluon_data_0trk_int[nCBins][ntrkbins]; 
Double_t data_0trk_int[nCBins][ntrkbins], data_no0trk_int[nCBins][ntrkbins], mc_0trk_int[nCBins][ntrkbins], mc_no0trk_int[nCBins][ntrkbins];
/*
TFile *jetcharge_pp_MC = TFile::Open("Pythia6_jetchg_recogen_no0trkjets_20190328.root");
//TFile *jetcharge_pp_MC = TFile::Open("Pythia6_jetchg_recogen_no0trkjets_20190320.root");
TFile *jetcharge_pp_data = TFile::Open("ppdata_jetchg_trkeffconst_recogen_no0trkjets_20190321.root");
//TFile *jetcharge_pp_data = TFile::Open("ppdata_fulleta_jetchg_trkeffconst_recogen_no0trkjets_20190329.root");
TFile *jetcharge_PbPb_MC = TFile::Open("P+H_jetchg_recogen_trkeffconst_no0trkjets_20190407.root");
//TFile *jetcharge_PbPb_MC = TFile::Open("P+H_jetchg_recogen_trkeffconst_no0trkjets_20190320.root");
TFile *jetcharge_PbPb_data = TFile::Open("PbPbdata_jetchg_trkeffconst_recogen_no0trkjets_20190408.root");
//TFile *jetcharge_PbPb_data = TFile::Open("PbPbdata_fulleta_jetchg_trkeffconst_recogen_no0trkjets_20190329.root");

TFile *jetcharge_unfolded_data = TFile::Open("data_unfolded_jetchg_trkeffconst_20190408.root");
*/
//full eta
//TFile *jetcharge_pp_MC = TFile::Open("Pythia6_fulleta_trkeffconst_jetchg_no0trkjets_20190415.root");
//TFile *jetcharge_PbPb_MC = TFile::Open("P+H_fulleta_trkeffconst_jetchg_no0trkjets_20190415.root");
TFile *jetcharge_pp_MC = TFile::Open("Pythia6_fulleta_trkeffconst_split_jetchg_no0trkjets_20190501.root");
TFile *jetcharge_PbPb_MC = TFile::Open("P+H_fulleta_trkeffconst_split_jetchg_no0trkjets_20190501.root");
TFile *jetcharge_pp_data = TFile::Open("ppdata_fulleta_jetchg_trkeffconst_recogen_no0trkjets_20190415.root");
TFile *jetcharge_PbPb_data = TFile::Open("PbPbdata_fulleta_jetchg_trkeffconst_no0trkjets_20190415.root");

//TFile *jetcharge_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffconst_kappa_effgengen_20190416.root");
//TFile *jetcharge_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffconst_kappa_recogen_20190415.root");

//TFile *jetcharge_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffconst_rescaled_kappa_recogen_20190423.root");
//TFile *jetcharge_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffconst_kappa_recogen_20190425.root");
TFile *jetcharge_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffnom_kappa_recogen_20190502.root");
//TFile *jetcharge_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffconst_rescaled_kappa_effgengen_20190418.root");
//TFile *jetcharge_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffconst_udrescaled_kappa_recogen_20190427.root");

TFile *jetcharge_trkinc_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffinc_kappa_recogen_20190416.root");
TFile *jetcharge_trkdec_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffdec_kappa_recogen_20190416.root");

//TFile *jetcharge_trkinc_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffdec_kappa_recogen_20190502.root");
//TFile *jetcharge_trkdec_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffdec_kappa_recogen_20190502.root");
TFile *jetcharge_trkposneg_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffconst_svd_kappa_recogen_20190416.root");
TFile *jetcharge_JER_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffconst_jer_kappa_recogen_20190416.root");
TFile *jetcharge_glusamp_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffnom_glusamp_kappa_recogen_20190502.root");

/*
TFile *jetcharge_trkinc0p99_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffvary0p99_kappa_effgengen_20190418.root");
TFile *jetcharge_trkinc0p98_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffvary0p98_kappa_effgengen_20190418.root");
TFile *jetcharge_trkinc0p97_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffvary0p97_kappa_effgengen_20190418.root");
TFile *jetcharge_trkinc0p96_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffvary0p96_kappa_effgengen_20190418.root");

TFile *jetcharge_trkinc1p01_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffvary1p01_kappa_effgengen_20190418.root");
TFile *jetcharge_trkinc1p02_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffvary1p02_kappa_effgengen_20190418.root");
TFile *jetcharge_trkinc1p03_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffvary1p03_kappa_effgengen_20190418.root");
TFile *jetcharge_trkinc1p04_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffvary1p04_kappa_effgengen_20190418.root");
*/

TFile *jetcharge_trkinc0p99_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeff0p99_kappa_recogen_20190425.root");
TFile *jetcharge_trkinc0p98_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeff0p98_kappa_recogen_20190425.root");
TFile *jetcharge_trkinc0p97_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeff0p97_kappa_recogen_20190425.root");
TFile *jetcharge_trkinc0p96_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeff0p96_kappa_recogen_20190425.root");

TFile *jetcharge_trkinc1p01_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeff1p01_kappa_recogen_20190425.root");
TFile *jetcharge_trkinc1p02_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeff1p02_kappa_recogen_20190425.root");
TFile *jetcharge_trkinc1p03_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeff1p03_kappa_recogen_20190425.root");
TFile *jetcharge_trkinc1p04_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeff1p04_kappa_recogen_20190425.root");

TFile *jetcharge_gengen_unfolded_data = TFile::Open("data_gengen_fulleta_unfolded_trk_jetchg_trkeffconst_kappa_recogen_20190503.root");

/*
TFile *jetcharge_trkinc_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffconst_kappa_rebin5_20190417.root");
TFile *jetcharge_trkdec_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffconst_kappa_rebin5_20190417.root");
TFile *jetcharge_trkposneg_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffconst_kappa_rebin5_20190417.root");
TFile *jetcharge_JER_unfolded_data = TFile::Open("data_fulleta_unfolded_trk_jetchg_trkeffconst_kappa_rebin5_20190417.root");
*/
TH2F *h_chg_refpt[nCBins][ntrkbins];
TH2F *h_chg_refpt_q[nCBins][ntrkbins];
TH2F *h_chg_refpt_g[nCBins][ntrkbins];
TH2F *h_chg_refpt_up[nCBins][ntrkbins];
TH2F *h_chg_refpt_down[nCBins][ntrkbins];
TH2F *h_chg_refpt_upb[nCBins][ntrkbins];
TH2F *h_chg_refpt_downb[nCBins][ntrkbins];
TH2F *h_chg_refpt_c[nCBins][ntrkbins];
TH2F *h_chg_refpt_s[nCBins][ntrkbins];
TH2F *h_chg_refpt_b[nCBins][ntrkbins];

TH2F *h_chg_refpt_sube0[nCBins][ntrkbins];
TH2F *h_chg_refpt_subenon0[nCBins][ntrkbins];

TH2F *h_chg_corrpt_data[nCBins][ntrkbins];

TH2F *h_chg_corrpt[nCBins][ntrkbins];
TH2F *h_chg_corrpt_q[nCBins][ntrkbins];
TH2F *h_chg_corrpt_g[nCBins][ntrkbins];
TH2F *h_chg_corrpt_up[nCBins][ntrkbins];
TH2F *h_chg_corrpt_down[nCBins][ntrkbins];
TH2F *h_chg_corrpt_upb[nCBins][ntrkbins];
TH2F *h_chg_corrpt_downb[nCBins][ntrkbins];
TH2F *h_chg_corrpt_c[nCBins][ntrkbins];
TH2F *h_chg_corrpt_s[nCBins][ntrkbins];
TH2F *h_chg_corrpt_b[nCBins][ntrkbins];

TH2F *h_chg_corrpt_g_kappa[nCBins][ntrkbins];
TH2F *h_chg_corrpt_upb_kappa[nCBins][ntrkbins];
TH2F *h_chg_corrpt_downb_kappa[nCBins][ntrkbins];
TH2F *h_chg_corrpt_c_kappa[nCBins][ntrkbins];
TH2F *h_chg_corrpt_s_kappa[nCBins][ntrkbins];
TH2F *h_chg_corrpt_b_kappa[nCBins][ntrkbins];

TH2F *h_chg_corrpt_data_trkeff[nCBins][ntrkbins][trkeffbins];
TH1D *h_chg_data_trkeff[nCBins][ntrkbins][trkeffbins];
TH2F *h_chg_corrpt_data_kappa_trkeff[nCBins][ntrkbins][trkeffbins];
TH1D *h_chg_data_kappa_trkeff[nCBins][ntrkbins][trkeffbins];

TH1D *h_chg_data_glusamp[nCBins][ntrkbins];

TH1D *h_chg_ref[nCBins][ntrkbins];
TH1D *h_chg_ref_q[nCBins][ntrkbins];
TH1D *h_chg_ref_g[nCBins][ntrkbins];
TH1D *h_chg_ref_up[nCBins][ntrkbins];
TH1D *h_chg_ref_down[nCBins][ntrkbins];
TH1D *h_chg_ref_upb[nCBins][ntrkbins];
TH1D *h_chg_ref_downb[nCBins][ntrkbins];
TH1D *h_chg_ref_c[nCBins][ntrkbins];
TH1D *h_chg_ref_s[nCBins][ntrkbins];
TH1D *h_chg_ref_b[nCBins][ntrkbins];

TH1D *h_chg_ref_sube0[nCBins][ntrkbins];
TH1D *h_chg_ref_subenon0[nCBins][ntrkbins];

//TH1D *h_chg_ref_gubdb[nCBins][ntrkbins];
//TH1D *h_chg_reco_ud_scaled[nCBins][ntrkbins];
//TH1D *h_chg_reco_gubdb[nCBins][ntrkbins];

TH1D *h_chg_data[nCBins][ntrkbins];
TH1D *h_chg_data_trkinc[nCBins][ntrkbins];
TH1D *h_chg_data_trkdec[nCBins][ntrkbins];
TH1D *h_chg_data_svd[nCBins][ntrkbins];
TH1D *h_chg_data_jer[nCBins][ntrkbins];
TH1D *h_chg_data_unfolded[nCBins][ntrkbins];

TH1D *h_chg_data_gengen[nCBins][ntrkbins];

TH1D *h_chg_data_trkrel[nCBins][ntrkbins];

TH1D *h_chg_data_kappa[nCBins][ntrkbins];
TH1D *h_chg_data_kappa_trkinc[nCBins][ntrkbins];
TH1D *h_chg_data_kappa_trkdec[nCBins][ntrkbins];
TH1D *h_chg_data_kappa_trkposneg[nCBins][ntrkbins];
TH1D *h_chg_data_kappa_jer[nCBins][ntrkbins];
TH1D *h_chg_data_kappa_unfolded[nCBins][ntrkbins];

TH1D *h_chg_reco[nCBins][ntrkbins];
TH1D *h_chg_reco_q[nCBins][ntrkbins];
TH1D *h_chg_reco_g[nCBins][ntrkbins];
//TH1D *h_chg_reco_up[nCBins][ntrkbins];
//TH1D *h_chg_reco_down[nCBins][ntrkbins];
TH1D *h_chg_reco_upb[nCBins][ntrkbins];
TH1D *h_chg_reco_downb[nCBins][ntrkbins];
TH1D *h_chg_reco_c[nCBins][ntrkbins];
TH1D *h_chg_reco_s[nCBins][ntrkbins];
TH1D *h_chg_reco_b[nCBins][ntrkbins];
TH1D *h_chg_reco_others[nCBins][ntrkbins];

TH1D *h_chg_gengensube0[nCBins][ntrkbins];
TH1D *h_chg_gengensube0_up[nCBins][ntrkbins];
TH1D *h_chg_gengensube0_d[nCBins][ntrkbins];

TH1D *h_chg_reco_kappa[nCBins][ntrkbins];
TH1D *h_chg_reco_q_kappa[nCBins][ntrkbins];
TH1D *h_chg_reco_g_kappa[nCBins][ntrkbins];
TH1D *h_chg_reco_up_kappa[nCBins][ntrkbins];
TH1D *h_chg_reco_down_kappa[nCBins][ntrkbins];
TH1D *h_chg_reco_upb_kappa[nCBins][ntrkbins];
TH1D *h_chg_reco_downb_kappa[nCBins][ntrkbins];
TH1D *h_chg_reco_c_kappa[nCBins][ntrkbins];
TH1D *h_chg_reco_s_kappa[nCBins][ntrkbins];
TH1D *h_chg_reco_b_kappa[nCBins][ntrkbins];
TH1D *h_chg_reco_others_kappa[nCBins][ntrkbins];

TH1D *h_ref_MC[nCBins];
TH1D *h_ref_MC_q[nCBins];
TH1D *h_ref_MC_g[nCBins];
TH1D *h_ref_MC_up[nCBins];
TH1D *h_ref_MC_d[nCBins];
TH1D *h_ref_MC_upbar[nCBins];
TH1D *h_ref_MC_dbar[nCBins];

TH1D *h_gen_full_0trks[nCBins][ntrkbins];
TH1D *h_gen_full_0trks_q[nCBins][ntrkbins];
TH1D *h_gen_full_0trks_g[nCBins][ntrkbins];

TH1D *h_reco_MC[nCBins];
TH1D *h_reco_MC_q[nCBins];
TH1D *h_reco_MC_g[nCBins];
TH1D *h_reco_MC_up[nCBins];
TH1D *h_reco_MC_d[nCBins];
TH1D *h_reco_MC_upbar[nCBins];
TH1D *h_reco_MC_dbar[nCBins];
TH1D *h_reco_MC_c[nCBins];
TH1D *h_reco_MC_s[nCBins];
TH1D *h_reco_MC_b[nCBins];
TH1D *h_reco_data[nCBins];
TH1D *h_reco_corr_0trks[nCBins][ntrkbins];
TH1D *h_reco_corr_0trks_q[nCBins][ntrkbins];
TH1D *h_reco_corr_0trks_g[nCBins][ntrkbins];
TH1D *h_reco_data_0trks[nCBins][ntrkbins];

TH1D *h_eta_MC[nCBins];
TH1D *h_eta_MC_g[nCBins];
TH1D *h_eta_MC_up[nCBins];
TH1D *h_eta_MC_d[nCBins];
TH1D *h_eta_data[nCBins];

TH1D *h_chg_reco_cop[nCBins][ntrkbins];
TH1D *h_chg_reco_up_cop[nCBins][ntrkbins];
TH1D *h_chg_reco_down_cop[nCBins][ntrkbins];
TH1D *h_chg_reco_gubdb_cop[nCBins][ntrkbins];
TH1D *h_chg_data_cop[nCBins][ntrkbins];
TH1D *h_chg_data_cop_sysm[nCBins][ntrkbins];
TH1D *h_chg_reco_g_cop[nCBins][ntrkbins];
TH1D *h_chg_reco_others_cop[nCBins][ntrkbins];

TH1D *h_chg_reco_kappa_cop[nCBins][ntrkbins];
TH1D *h_chg_reco_up_kappa_cop[nCBins][ntrkbins];
TH1D *h_chg_reco_down_kappa_cop[nCBins][ntrkbins];
TH1D *h_chg_reco_gubdb_kappa_cop[nCBins][ntrkbins];
TH1D *h_chg_data_kappa_cop[nCBins][ntrkbins];
TH1D *h_chg_reco_g_kappa_cop[nCBins][ntrkbins];
TH1D *h_chg_reco_others_kappa_cop[nCBins][ntrkbins];
TH1D *h_chg_data_kappa_cop_sysm[nCBins][ntrkbins];

TH1D *h_chg_reco_gubdb_afterfit[nCBins][ntrkbins];
TH1D *h_chg_reco_fit_total[nCBins][ntrkbins];
TH1D *h_chg_reco_gubdb_afterfit_kappa[nCBins][ntrkbins];
TH1D *h_chg_reco_fit_total_kappa[nCBins][ntrkbins];

TH1D* h_jetchg_avg_data[ntrkbins];
TH1D* h_jetchg_avg_data_trkinc[ntrkbins];
TH1D* h_jetchg_avg_data_trkdec[ntrkbins];
TH1D* h_jetchg_avg_MC[ntrkbins];
TH1D* h_jetchg_avg_data_pyquen[ntrkbins];
TH1D* h_jetchg_avg_data_sysm[ntrkbins];
TH1D* h_jetchg_rms_data[ntrkbins];
TH1D* h_jetchg_rms_data_sysm[ntrkbins];
TH1D* h_jetchg_rms_data_trkinc[ntrkbins];
TH1D* h_jetchg_rms_data_trkdec[ntrkbins];
TH1D* h_jetchg_rms_data_svd[ntrkbins];
TH1D* h_jetchg_rms_data_jer[ntrkbins];
TH1D* h_jetchg_rms_data_pyquen[ntrkbins];
TH1D* h_jetchg_rms_data_pyquencoll[ntrkbins];
TH1D* h_jetchg_rms_MC[ntrkbins];

//ratio histos
TH1D *h_chg_reco_ratio_2param[nCBins][ntrkbins];
TH1D *h_chg_reco_ratio_2param_kappa[nCBins][ntrkbins];
TH1D *h_chg_reco_ratio_trkeff[nCBins][ntrkbins][trkeffbins];
TH1D *h_chg_reco_ratio_kappa_trkeff[nCBins][ntrkbins][trkeffbins];
TH1D *h_chg_reco_ratio[nCBins][ntrkbins];

//functions
TF1 *f_chg_qg_data[nCBins][ntrkbins];
TF1 *f_chg_qg_kappa[nCBins][ntrkbins];

TF1 *f_chg_ud_data[nCBins][ntrkbins];

TGraphErrors *up_data[nCBins];
TGraphErrors *down_data[nCBins];
TGraphErrors *gluon_data[nCBins];
TGraphErrors *quark_data[nCBins];
//TGraphAsymmErrors *gluon_data_asymm_sys[nCBins];
TGraphErrors *gluon_data_asymm_sys[nCBins];

TGraphErrors *gluon_data_old[nCBins];
TGraphErrors *gluon_data_arc[nCBins];

TGraphErrors *up_data_kappa[nCBins];
TGraphErrors *down_data_kappa[nCBins];
TGraphErrors *gluon_data_kappa[nCBins];
TGraphErrors *quark_data_kappa[nCBins];
//TGraphAsymmErrors *gluon_data_asymm_sys[nCBins];
TGraphErrors *gluon_data_kappa_sysm[nCBins];

THStack *h_chg_data_stk[nCBins][ntrkbins];
THStack *h_chg_data_stk_kappa[nCBins][ntrkbins];
THStack *h_chg_reco_stk[nCBins][ntrkbins];
THStack *h_chg_reco_stk_kappa[nCBins][ntrkbins];

TH1D *h_trkpt_reco[nCBins]; 
TH1D *h_trkpt_pos_reco[nCBins]; 
TH1D *h_trkpt_neg_reco[nCBins]; 

TH1D *h_trkpt_gen[nCBins]; 
TH1D *h_trkpt_pos_gen[nCBins]; 
TH1D *h_trkpt_neg_gen[nCBins]; 

TH1D *h_trkpt_reco_gen[nCBins]; 
TH1D *h_trkpt_pos_reco_gen[nCBins]; 
TH1D *h_trkpt_neg_reco_gen[nCBins]; 

TH1D *h_syst[nCBins];  
TH1D *h_syst_MCstat[nCBins];  
TH1D *h_stat[nCBins];  
TH1D *h_syst_trk[nCBins];  
TH1D *h_syst_posneg[nCBins];  
TH1D *h_syst_unf[nCBins];  
TH1D *h_syst_JER[nCBins];  
TH1D *h_syst_0trk[nCBins]; 

TH1D *h_syst_kappa[nCBins];  
TH1D *h_syst_MCstat_kappa[nCBins];  
TH1D *h_stat_kappa[nCBins];  
TH1D *h_syst_trk_kappa[nCBins];  
TH1D *h_syst_posneg_kappa[nCBins];  
TH1D *h_syst_unf_kappa[nCBins];  
TH1D *h_syst_JER_kappa[nCBins];  
TH1D *h_syst_0trk_kappa[nCBins]; 

Double_t up_fit_high[nCBins][ntrkbins], down_fit_high[nCBins][ntrkbins], gluon_fit_high[nCBins][ntrkbins], quark_fit_high[nCBins][ntrkbins]; 
Double_t up_fit_low[nCBins][ntrkbins], down_fit_low[nCBins][ntrkbins], gluon_fit_low[nCBins][ntrkbins],quark_fit_low[nCBins][ntrkbins]; 
Double_t up_fit_posneg[nCBins][ntrkbins], down_fit_posneg[nCBins][ntrkbins], gluon_fit_posneg[nCBins][ntrkbins], quark_fit_posneg[nCBins][ntrkbins]; 
Double_t up_fit_jer[nCBins][ntrkbins], down_fit_jer[nCBins][ntrkbins], gluon_fit_jer[nCBins][ntrkbins], quark_fit_jer[nCBins][ntrkbins]; 
Double_t up_fit_glusamp[nCBins][ntrkbins], down_fit_glusamp[nCBins][ntrkbins], gluon_fit_glusamp[nCBins][ntrkbins], quark_fit_glusamp[nCBins][ntrkbins]; 

Double_t up_fit_high_kappa[nCBins][ntrkbins], down_fit_high_kappa[nCBins][ntrkbins], gluon_fit_high_kappa[nCBins][ntrkbins], quark_fit_high_kappa[nCBins][ntrkbins]; 
Double_t up_fit_low_kappa[nCBins][ntrkbins], down_fit_low_kappa[nCBins][ntrkbins], gluon_fit_low_kappa[nCBins][ntrkbins],quark_fit_low_kappa[nCBins][ntrkbins]; 
Double_t up_fit_posneg_kappa[nCBins][ntrkbins], down_fit_posneg_kappa[nCBins][ntrkbins], gluon_fit_posneg_kappa[nCBins][ntrkbins], quark_fit_posneg_kappa[nCBins][ntrkbins]; 
Double_t up_fit_kappa_jer[nCBins][ntrkbins], down_fit_kappa_jer[nCBins][ntrkbins], gluon_fit_kappa_jer[nCBins][ntrkbins], quark_fit_kappa_jer[nCBins][ntrkbins]; 

void abs_histos(TH1D *h_old, TH1D *h_new){
  for(int ibin=1; ibin<h_old->GetNbinsX()+1; ibin++){
    if(ibin <= h_old->GetNbinsX()/2) h_new->SetBinContent(ibin,0.);
    else if(ibin > 1+(h_old->GetNbinsX()/2)) h_new->SetBinContent(ibin,h_old->GetBinContent(ibin)+h_old->GetBinContent(h_old->GetNbinsX()-ibin));
  }  
}

double bincontent[500];

void diff_histos(TH1D *h_old, TH1D *h_new){
  for(int ibin=1; ibin<h_old->GetNbinsX()+1; ibin++){
    bincontent[ibin] = h_old->GetBinContent(ibin);
    h_new->SetBinContent(ibin,h_old->GetBinContent(ibin) - h_old->GetBinContent(h_old->GetNbinsX()+1 - ibin));
  }  
}

void inv_histos(TH1D *h_old, TH1D *h_new){
  for(int ibin=1; ibin<h_old->GetNbinsX()+1; ibin++){
    bincontent[ibin] = h_old->GetBinContent(ibin);
    h_new->SetBinContent(ibin,-1.*h_old->GetBinContent(ibin));
  }  
}

void get_histos(bool do_data = false){
    for(int ibin=0;ibin<nCBins;ibin++){
        if(ibin==0){
            h_reco_MC[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_reco_corr_cent"+cent[ibin]));
            h_reco_data[ibin] = (TH1D*)jetcharge_pp_data->Get((TString)("h_reco_corr_cent"+cent[ibin]));
            h_reco_MC_q[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_reco_corr_q_cent"+cent[ibin]));
            h_reco_MC_g[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_reco_corr_g_cent"+cent[ibin]));
            h_reco_MC_up[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_reco_corr_up_cent"+cent[ibin]));
            h_reco_MC_d[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_reco_corr_d_cent"+cent[ibin]));
            h_reco_MC_upbar[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_reco_corr_upbar_cent"+cent[ibin]));
            h_reco_MC_dbar[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_reco_corr_dbar_cent"+cent[ibin]));
            h_reco_MC_c[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_reco_corr_c_cent"+cent[ibin]));
            h_reco_MC_s[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_reco_corr_s_cent"+cent[ibin]));
            h_reco_MC_b[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_reco_corr_b_cent"+cent[ibin]));

            h_eta_MC[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_eta_full_cent"+cent[ibin]));
            h_eta_data[ibin] = (TH1D*)jetcharge_pp_data->Get((TString)("h_eta_full_cent"+cent[ibin]));
            h_eta_MC_g[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_eta_full_g_cent"+cent[ibin]));
            h_eta_MC_up[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_eta_full_up_cent"+cent[ibin]));
            h_eta_MC_d[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_eta_full_d_cent"+cent[ibin]));

            h_ref_MC[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_gen_full_cent"+cent[ibin]));
            h_ref_MC_q[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_gen_full_q_cent"+cent[ibin]));
            h_ref_MC_g[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_gen_full_g_cent"+cent[ibin]));
            h_ref_MC_up[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_gen_full_up_cent"+cent[ibin]));
            h_ref_MC_d[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_gen_full_d_cent"+cent[ibin]));
            h_ref_MC_upbar[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_gen_full_upbar_cent"+cent[ibin]));
            h_ref_MC_dbar[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_gen_full_dbar_cent"+cent[ibin]));

            h_trkpt_reco[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_trkpt_reco_cent"+cent[ibin]));
            h_trkpt_neg_reco[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_trkpt_neg_reco_cent"+cent[ibin]));
            h_trkpt_pos_reco[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_trkpt_pos_reco_cent"+cent[ibin]));

            h_trkpt_gen[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_trkpt_gen_cent"+cent[ibin]));
            h_trkpt_neg_gen[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_trkpt_neg_gen_cent"+cent[ibin]));
            h_trkpt_pos_gen[ibin] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_trkpt_pos_gen_cent"+cent[ibin]));
        }
        else{
            h_reco_MC[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_reco_corr_cent"+cent[ibin-1]));       
            h_reco_data[ibin] = (TH1D*)jetcharge_PbPb_data->Get((TString)("h_reco_corr_cent"+cent[ibin-1]));
            h_reco_MC_q[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_reco_corr_q_cent"+cent[ibin-1]));
            h_reco_MC_g[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_reco_corr_g_cent"+cent[ibin-1]));
            h_reco_MC_up[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_reco_corr_up_cent"+cent[ibin-1]));
            h_reco_MC_d[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_reco_corr_d_cent"+cent[ibin-1]));
            h_reco_MC_upbar[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_reco_corr_upbar_cent"+cent[ibin-1]));
            h_reco_MC_dbar[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_reco_corr_dbar_cent"+cent[ibin-1]));
            h_reco_MC_c[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_reco_corr_c_cent"+cent[ibin-1]));
            h_reco_MC_s[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_reco_corr_s_cent"+cent[ibin-1]));
            h_reco_MC_b[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_reco_corr_b_cent"+cent[ibin-1]));

            h_eta_MC[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_eta_full_cent"+cent[ibin-1]));
            h_eta_data[ibin] = (TH1D*)jetcharge_PbPb_data->Get((TString)("h_eta_full_cent"+cent[ibin-1]));
            h_eta_MC_g[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_eta_full_g_cent"+cent[ibin-1]));
            h_eta_MC_up[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_eta_full_up_cent"+cent[ibin-1]));
            h_eta_MC_d[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_eta_full_d_cent"+cent[ibin-1]));

            h_ref_MC[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_gen_full_cent"+cent[ibin-1]));
            h_ref_MC_q[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_gen_full_q_cent"+cent[ibin-1]));
            h_ref_MC_g[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_gen_full_g_cent"+cent[ibin-1]));
            h_ref_MC_up[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_gen_full_up_cent"+cent[ibin-1]));
            h_ref_MC_d[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_gen_full_d_cent"+cent[ibin-1]));
            h_ref_MC_upbar[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_gen_full_upbar_cent"+cent[ibin-1]));
            h_ref_MC_dbar[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_gen_full_dbar_cent"+cent[ibin-1]));

            h_trkpt_reco[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_trkpt_reco_cent"+cent[ibin-1]));
            h_trkpt_neg_reco[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_trkpt_neg_reco_cent"+cent[ibin-1]));
            h_trkpt_pos_reco[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_trkpt_pos_reco_cent"+cent[ibin-1]));

            h_trkpt_gen[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_trkpt_gen_cent"+cent[ibin-1]));
            h_trkpt_neg_gen[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_trkpt_neg_gen_cent"+cent[ibin-1]));
            h_trkpt_pos_gen[ibin] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_trkpt_pos_gen_cent"+cent[ibin-1]));
        }

        h_eta_MC[ibin] -> Scale(1./h_eta_MC[ibin]->Integral());
        h_eta_data[ibin] -> Scale(1./h_eta_data[ibin]->Integral());
        h_eta_MC_g[ibin] -> Scale(1./h_eta_MC_g[ibin]->Integral());
        h_eta_MC_up[ibin] -> Scale(1./h_eta_MC_up[ibin]->Integral());
        h_eta_MC_d[ibin] -> Scale(1./h_eta_MC_d[ibin]->Integral());
    }  

    for(int ibin2=0;ibin2<ntrkbins;ibin2++){
        for(int ibin=0;ibin<nCBins;ibin++){
            if(ibin==0){
              h_reco_corr_0trks[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_reco_corr_0trks_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_reco_corr_0trks_q[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_reco_corr_0trks_q_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_reco_corr_0trks_g[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_reco_corr_0trks_g_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_reco_data_0trks[ibin][ibin2] = (TH1D*)jetcharge_pp_data->Get((TString)("h_reco_corr_0trks_cent"+cent[ibin]+"_trk"+trk[ibin2]));

              h_gen_full_0trks[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_gen_full_0trks_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_gen_full_0trks_q[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_gen_full_0trks_q_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_gen_full_0trks_g[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_gen_full_0trks_g_cent"+cent[ibin]+"_trk"+trk[ibin2]));

              h_chg_corrpt_data[ibin][ibin2] = (TH2F*)jetcharge_pp_data->Get((TString)("h_chg_corrpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 

              for(int ibin3=0; ibin3<trkeffbins;ibin3++){
                  h_chg_corrpt_data_trkeff[ibin][ibin2][ibin3] = (TH2F*)jetcharge_pp_data->Get((TString)("h_chg_corrpt_trkeff_cent"+cent[ibin]+"_trk"+trk[ibin2]+"_eff"+trkeff[ibin3]));           
                  //if(do_kappa_fitting) h_chg_corrpt_data_kappa_trkeff[ibin][ibin2][ibin3] = (TH2F*)closure_histos_pp_data->Get((TString)("h_chg_corrpt_kappa_trkeff_cent"+cent[ibin]+"_trk"+trk[ibin2]+"_eff"+trkeff[ibin3]));           
              }

              h_chg_refpt[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_refpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
              h_chg_refpt_q[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_refpt_q_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_refpt_g[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_refpt_g_cent"+cent[ibin]+"_trk"+trk[ibin2]));

              h_chg_refpt_up[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_refpt_up_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_refpt_down[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_refpt_d_cent"+cent[ibin]+"_trk"+trk[ibin2]));      
              h_chg_refpt_upb[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_refpt_upbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_refpt_downb[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_refpt_dbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));

              h_chg_refpt_c[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_refpt_c_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_refpt_s[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_refpt_s_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_refpt_b[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_refpt_b_cent"+cent[ibin]+"_trk"+trk[ibin2]));

              h_chg_corrpt[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
              h_chg_corrpt[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_cent"+cent[ibin]+"_trk"+trk[ibin2])); 
              h_chg_corrpt_q[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_q_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_corrpt_g[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_g_cent"+cent[ibin]+"_trk"+trk[ibin2]));

              h_chg_corrpt_g_kappa[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_kappa_g_cent"+cent[ibin]+"_trk"+trk[ibin2]));

              h_chg_corrpt_up[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_up_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_corrpt_down[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_d_cent"+cent[ibin]+"_trk"+trk[ibin2]));      
              h_chg_corrpt_upb[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_upbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_corrpt_downb[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_dbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));

              h_chg_corrpt_c[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_c_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_corrpt_s[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_s_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_corrpt_b[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_b_cent"+cent[ibin]+"_trk"+trk[ibin2]));

              h_chg_corrpt_upb_kappa[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_upbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_corrpt_downb_kappa[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_dbar_cent"+cent[ibin]+"_trk"+trk[ibin2]));

              h_chg_corrpt_c_kappa[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_kappa_c_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_corrpt_s_kappa[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_kappa_s_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_corrpt_b_kappa[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get((TString)("h_chg_corrpt_kappa_b_cent"+cent[ibin]+"_trk"+trk[ibin2]));
            }
            else{
              h_reco_corr_0trks[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_reco_corr_0trks_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_reco_corr_0trks_q[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_reco_corr_0trks_q_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_reco_corr_0trks_g[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_reco_corr_0trks_g_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_reco_data_0trks[ibin][ibin2] = (TH1D*)jetcharge_PbPb_data->Get((TString)("h_reco_corr_0trks_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

              h_gen_full_0trks[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_gen_full_0trks_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_gen_full_0trks_q[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_gen_full_0trks_q_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_gen_full_0trks_g[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_gen_full_0trks_g_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

              h_chg_corrpt_data[ibin][ibin2] = (TH2F*)jetcharge_PbPb_data->Get((TString)("h_chg_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 

              for(int ibin3=0; ibin3<trkeffbins;ibin3++){
                  h_chg_corrpt_data_trkeff[ibin][ibin2][ibin3] = (TH2F*)jetcharge_PbPb_data->Get((TString)("h_chg_corrpt_trkeff_cent"+cent[ibin-1]+"_trk"+trk[ibin2]+"_eff"+trkeff[ibin3]));           
                  //if(do_kappa_fitting) h_chg_corrpt_data_kappa_trkeff[ibin][ibin2][ibin3] = (TH2F*)closure_histos_PbPb_data->Get((TString)("h_chg_corrpt_kappa_trkeff_cent"+cent[ibin-1]+"_trk"+trk[ibin2]+"_eff"+trkeff[ibin3]));           
              }

              h_chg_refpt[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_refpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
              h_chg_refpt_q[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_refpt_q_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_refpt_g[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_refpt_g_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

              h_chg_refpt_up[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_refpt_up_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_refpt_down[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_refpt_d_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));      
              h_chg_refpt_upb[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_refpt_upbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_refpt_downb[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_refpt_dbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

              h_chg_refpt_c[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_refpt_c_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_refpt_s[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_refpt_s_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_refpt_b[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_refpt_b_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

              h_chg_refpt_sube0[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_refpt_sube0_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
              //h_chg_refpt_subenon0[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_refpt_subenon0_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 

              h_chg_corrpt[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_cent"+cent[ibin-1]+"_trk"+trk[ibin2])); 
              h_chg_corrpt_q[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_q_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_corrpt_g[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_g_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

              h_chg_corrpt_g_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_g_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

              h_chg_corrpt_up[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_up_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_corrpt_down[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_d_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));      
              h_chg_corrpt_upb[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_upbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_corrpt_downb[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_dbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

              h_chg_corrpt_c[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_c_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_corrpt_s[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_s_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_corrpt_b[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_b_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));              

              h_chg_corrpt_upb_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_upbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_corrpt_downb_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_dbar_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

              h_chg_corrpt_c_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_c_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_corrpt_s_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_s_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_corrpt_b_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get((TString)("h_chg_corrpt_kappa_b_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));              
            }

            //h_chg_data[ibin][ibin2] = h_chg_corrpt_data[ibin][ibin2]-> ProjectionY(Form("h_chg_data_%d_%d",ibin,ibin2),h_chg_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_data[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_data[ibin][ibin2] = (TH1D*)jetcharge_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            //if(ibin==0) h_chg_data[ibin][ibin2] = (TH1D*)jetcharge_pp_data->Get((TString)("h_chg_recopt_cent"+cent[ibin]+"_trk"+trk[ibin2]));
            //else h_chg_data[ibin][ibin2] = (TH1D*)jetcharge_PbPb_data->Get((TString)("h_chg_recopt_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
            //h_chg_data[ibin][ibin2]->Rebin(10);
            h_chg_data[ibin][ibin2]->Scale(1./h_chg_data[ibin][ibin2]->Integral()); 

            h_chg_data_trkinc[ibin][ibin2] = (TH1D*)jetcharge_trkinc_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_trkinc[ibin][ibin2]->Scale(1./h_chg_data_trkinc[ibin][ibin2]->Integral()); 
            h_chg_data_trkdec[ibin][ibin2] = (TH1D*)jetcharge_trkdec_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_trkdec[ibin][ibin2]->Scale(1./h_chg_data_trkdec[ibin][ibin2]->Integral()); 
            h_chg_data_svd[ibin][ibin2] = (TH1D*)jetcharge_trkposneg_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_svd[ibin][ibin2]->Scale(1./h_chg_data_svd[ibin][ibin2]->Integral()); 
            h_chg_data_jer[ibin][ibin2] = (TH1D*)jetcharge_JER_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_jer[ibin][ibin2]->Scale(1./h_chg_data_jer[ibin][ibin2]->Integral()); 
            h_chg_data_glusamp[ibin][ibin2] = (TH1D*)jetcharge_glusamp_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_glusamp[ibin][ibin2]->Scale(1./h_chg_data_glusamp[ibin][ibin2]->Integral()); 

            h_chg_data_gengen[ibin][ibin2] = (TH1D*)jetcharge_gengen_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));          
            h_chg_data_gengen[ibin][ibin2]->Scale(1./h_chg_data_gengen[ibin][ibin2]->Integral());

//kappa
            h_chg_data_kappa[ibin][ibin2] = (TH1D*)jetcharge_unfolded_data->Get((TString)("h_jetcharge_response_data_kappa_cent"+cent[ibin]+"_trk"+trk[ibin2]));           

            //h_chg_data_kappa[ibin][ibin2] = h_chg_corrpt_data_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_data_kappa_%d_%d",ibin,ibin2),h_chg_corrpt_data_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_data_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            //if(ibin==0) h_chg_data_kappa[ibin][ibin2] = (TH1D*)jetcharge_pp_data->Get((TString)("h_chg_recopt_kappa_cent"+cent[ibin]+"_trk"+trk[ibin2]));
            //else h_chg_data_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_data->Get((TString)("h_chg_recopt_kappa_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
            //h_chg_data_kappa[ibin][ibin2]->Rebin(10);
            h_chg_data_kappa[ibin][ibin2]->Scale(1./h_chg_data_kappa[ibin][ibin2]->Integral()); 

            h_chg_data_kappa_trkinc[ibin][ibin2] = (TH1D*)jetcharge_trkinc_unfolded_data->Get((TString)("h_jetcharge_response_data_kappa_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_kappa_trkinc[ibin][ibin2]->Scale(1./h_chg_data_kappa_trkinc[ibin][ibin2]->Integral()); 
            h_chg_data_kappa_trkdec[ibin][ibin2] = (TH1D*)jetcharge_trkdec_unfolded_data->Get((TString)("h_jetcharge_response_data_kappa_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_kappa_trkdec[ibin][ibin2]->Scale(1./h_chg_data_kappa_trkdec[ibin][ibin2]->Integral()); 
            h_chg_data_kappa_trkposneg[ibin][ibin2] = (TH1D*)jetcharge_trkposneg_unfolded_data->Get((TString)("h_jetcharge_response_data_kappa_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_kappa_trkposneg[ibin][ibin2]->Scale(1./h_chg_data_kappa_trkposneg[ibin][ibin2]->Integral()); 
            h_chg_data_kappa_jer[ibin][ibin2] = (TH1D*)jetcharge_JER_unfolded_data->Get((TString)("h_jetcharge_response_data_kappa_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_kappa_jer[ibin][ibin2]->Scale(1./h_chg_data_kappa_jer[ibin][ibin2]->Integral()); 
/*
            for(int ibin3=0; ibin3<6;ibin3++){        
                h_chg_data_trkeff[ibin][ibin2][ibin3] = (TH1D*)jetcharge_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
                h_chg_data_trkeff[ibin][ibin2][ibin3]->Scale(1./h_chg_data_trkeff[ibin][ibin2][ibin3]->Integral()); 
            }
*/
            h_chg_data_trkeff[ibin][ibin2][0] = (TH1D*)jetcharge_trkdec_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_trkeff[ibin][ibin2][0]->Scale(1./h_chg_data_trkeff[ibin][ibin2][0]->Integral()); 
            h_chg_data_trkeff[ibin][ibin2][1] = (TH1D*)jetcharge_trkinc0p96_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_trkeff[ibin][ibin2][1]->Scale(1./h_chg_data_trkeff[ibin][ibin2][1]->Integral()); 
            h_chg_data_trkeff[ibin][ibin2][2] = (TH1D*)jetcharge_trkinc0p97_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_trkeff[ibin][ibin2][2]->Scale(1./h_chg_data_trkeff[ibin][ibin2][2]->Integral()); 
            h_chg_data_trkeff[ibin][ibin2][3] = (TH1D*)jetcharge_trkinc0p98_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_trkeff[ibin][ibin2][3]->Scale(1./h_chg_data_trkeff[ibin][ibin2][3]->Integral()); 
            h_chg_data_trkeff[ibin][ibin2][4] = (TH1D*)jetcharge_trkinc0p99_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_trkeff[ibin][ibin2][4]->Scale(1./h_chg_data_trkeff[ibin][ibin2][4]->Integral());
            h_chg_data_trkeff[ibin][ibin2][5] = (TH1D*)jetcharge_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_trkeff[ibin][ibin2][5]->Scale(1./h_chg_data_trkeff[ibin][ibin2][5]->Integral());
            h_chg_data_trkeff[ibin][ibin2][6] = (TH1D*)jetcharge_trkinc1p01_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_trkeff[ibin][ibin2][6]->Scale(1./h_chg_data_trkeff[ibin][ibin2][6]->Integral()); 
            h_chg_data_trkeff[ibin][ibin2][7] = (TH1D*)jetcharge_trkinc1p02_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_trkeff[ibin][ibin2][7]->Scale(1./h_chg_data_trkeff[ibin][ibin2][7]->Integral()); 
            h_chg_data_trkeff[ibin][ibin2][8] = (TH1D*)jetcharge_trkinc1p03_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_trkeff[ibin][ibin2][8]->Scale(1./h_chg_data_trkeff[ibin][ibin2][8]->Integral()); 
            h_chg_data_trkeff[ibin][ibin2][9] = (TH1D*)jetcharge_trkinc1p04_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_trkeff[ibin][ibin2][9]->Scale(1./h_chg_data_trkeff[ibin][ibin2][9]->Integral()); 
            h_chg_data_trkeff[ibin][ibin2][10] = (TH1D*)jetcharge_trkinc_unfolded_data->Get((TString)("h_jetcharge_response_data_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_data_trkeff[ibin][ibin2][10]->Scale(1./h_chg_data_trkeff[ibin][ibin2][10]->Integral()); 

            h_chg_reco[ibin][ibin2] = (TH1D*)jetcharge_unfolded_data->Get((TString)("h_jetcharge_response_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            //if(ibin==0) h_chg_reco[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_chg_genpt_sube0_cent"+cent[ibin]+"_trk"+trk[ibin2]));
            //else h_chg_reco[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_chg_genpt_sube0_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
            h_chg_reco[ibin][ibin2]->Rebin(10);
            h_chg_reco[ibin][ibin2]->Scale(1./h_chg_reco[ibin][ibin2]->Integral());

            if(ibin==0) h_chg_gengensube0[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_chg_gengensube0pt_cent"+cent[ibin]+"_trk"+trk[ibin2]));
            else {
                h_chg_gengensube0[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_chg_gengensube0pt_g_oth_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
                h_chg_gengensube0_up[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_chg_gengensube0pt_up_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
                h_chg_gengensube0_d[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_chg_gengensube0pt_d_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

                h_chg_gengensube0[ibin][ibin2]->Scale(1./h_chg_gengensube0[ibin][ibin2]->Integral());            
                h_chg_gengensube0_up[ibin][ibin2]->Scale(1./h_chg_gengensube0_up[ibin][ibin2]->Integral());            
                h_chg_gengensube0_d[ibin][ibin2]->Scale(1./h_chg_gengensube0_d[ibin][ibin2]->Integral());            

                h_chg_gengensube0[ibin][ibin2]->Scale(0.762);
                h_chg_gengensube0_up[ibin][ibin2]->Scale(0.113);
                h_chg_gengensube0_d[ibin][ibin2]->Scale(0.127);
                h_chg_gengensube0[ibin][ibin2]->Add(h_chg_gengensube0_up[ibin][ibin2]);
                h_chg_gengensube0[ibin][ibin2]->Add(h_chg_gengensube0_d[ibin][ibin2]);
            }

            h_chg_gengensube0[ibin][ibin2]->Scale(1./h_chg_gengensube0[ibin][ibin2]->Integral());

            h_chg_reco_kappa[ibin][ibin2] = (TH1D*)jetcharge_unfolded_data->Get((TString)("h_jetcharge_response_kappa_cent"+cent[ibin]+"_trk"+trk[ibin2]));           
            h_chg_reco_kappa[ibin][ibin2]->Scale(1./h_chg_reco_kappa[ibin][ibin2]->Integral());

            h_chg_reco_q[ibin][ibin2] = h_chg_corrpt_q[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_q_%d_%d",ibin,ibin2),h_chg_corrpt_q[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_q[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_reco_g[ibin][ibin2] = h_chg_corrpt_g[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_g_%d_%d",ibin,ibin2),h_chg_corrpt_g[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_g[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            //if(ibin==0) h_chg_reco_up[ibin][ibin2] = h_chg_corrpt_up[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_up_%d_%d",ibin,ibin2),h_chg_corrpt_up[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_up[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            //if(ibin==0) h_chg_reco_down[ibin][ibin2] = h_chg_corrpt_down[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_down_%d_%d",ibin,ibin2),h_chg_corrpt_down[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_down[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_reco_upb[ibin][ibin2] = h_chg_corrpt_upb[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_upb_%d_%d",ibin,ibin2),h_chg_corrpt_upb[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_upb[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_reco_downb[ibin][ibin2] = h_chg_corrpt_downb[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_downb_%d_%d",ibin,ibin2),h_chg_corrpt_downb[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_downb[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_reco_c[ibin][ibin2] = h_chg_corrpt_c[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_c_%d_%d",ibin,ibin2),h_chg_corrpt_c[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_c[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_reco_s[ibin][ibin2] = h_chg_corrpt_s[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_s_%d_%d",ibin,ibin2),h_chg_corrpt_s[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_s[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_reco_b[ibin][ibin2] = h_chg_corrpt_b[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_b_%d_%d",ibin,ibin2),h_chg_corrpt_b[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_b[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");

            if(ibin==0){
              h_chg_reco_up[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_chg_genpt_sube0_up_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_reco_down[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_chg_genpt_sube0_d_cent"+cent[ibin]+"_trk"+trk[ibin2]));      
              h_chg_reco_gubdb[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_chg_genpt_sube0_g_oth_cent"+cent[ibin]+"_trk"+trk[ibin2]));

              h_chg_reco_up_kappa[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_chg_genpt_sube0_up_kappa_cent"+cent[ibin]+"_trk"+trk[ibin2]));
              h_chg_reco_down_kappa[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_chg_genpt_sube0_d_kappa_cent"+cent[ibin]+"_trk"+trk[ibin2]));      
              h_chg_reco_gubdb_kappa[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get((TString)("h_chg_genpt_sube0_g_oth_kappa_cent"+cent[ibin]+"_trk"+trk[ibin2]));
            }
            if(ibin>0){
              h_chg_reco_up[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_chg_genpt_sube0_up_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_reco_down[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_chg_genpt_sube0_d_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));      
              h_chg_reco_gubdb[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_chg_genpt_sube0_g_oth_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));

              h_chg_reco_up_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_chg_genpt_sube0_up_kappa_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
              h_chg_reco_down_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_chg_genpt_sube0_d_kappa_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));      
              h_chg_reco_gubdb_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get((TString)("h_chg_genpt_sube0_g_oth_kappa_cent"+cent[ibin-1]+"_trk"+trk[ibin2]));
            }

            h_chg_reco_g[ibin][ibin2]->Rebin(10);
            h_chg_reco_g[ibin][ibin2]->Scale(1./h_chg_reco_g[ibin][ibin2]->Integral());

            h_chg_ref[ibin][ibin2] = h_chg_refpt[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_%d_%d",ibin,ibin2),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_ref_q[ibin][ibin2] = h_chg_refpt_q[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_q_%d_%d",ibin,ibin2),h_chg_refpt_q[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_q[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_ref_g[ibin][ibin2] = h_chg_refpt_g[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_g_%d_%d",ibin,ibin2),h_chg_refpt_g[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_g[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_ref_up[ibin][ibin2] = h_chg_refpt_up[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_up_%d_%d",ibin,ibin2),h_chg_refpt_up[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_up[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_ref_down[ibin][ibin2] = h_chg_refpt_down[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_down_%d_%d",ibin,ibin2),h_chg_refpt_down[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_down[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_ref_upb[ibin][ibin2] = h_chg_refpt_upb[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_upb_%d_%d",ibin,ibin2),h_chg_refpt_upb[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_upb[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_ref_downb[ibin][ibin2] = h_chg_refpt_downb[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_downb_%d_%d",ibin,ibin2),h_chg_refpt_downb[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_downb[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_ref_c[ibin][ibin2] = h_chg_refpt_c[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_c_%d_%d",ibin,ibin2),h_chg_refpt_c[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_c[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_ref_s[ibin][ibin2] = h_chg_refpt_s[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_s_%d_%d",ibin,ibin2),h_chg_refpt_s[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_s[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_ref_b[ibin][ibin2] = h_chg_refpt_b[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_b_%d_%d",ibin,ibin2),h_chg_refpt_b[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt_b[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");

            if(ibin>0){       
              h_chg_ref_sube0[ibin][ibin2] = h_chg_refpt_sube0[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_sube0_%d_%d",ibin,ibin2),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
              //h_chg_ref_subenon0[ibin][ibin2] = h_chg_refpt_subenon0[ibin][ibin2]-> ProjectionY(Form("h_chg_ref_subenon0_%d_%d",ibin,ibin2),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_refpt[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            }
/*
            //defining gubdb templates
            if(ibin==0){
              h_chg_reco_gubdb[ibin][ibin2] = (TH1D*)h_chg_reco_g[ibin][ibin2]->Clone(Form("h_chg_reco_gubdb_%d_%d",ibin,ibin2));
              h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_upb[ibin][ibin2]);
              h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_downb[ibin][ibin2]);
              h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_c[ibin][ibin2]);
              h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_s[ibin][ibin2]);
              h_chg_reco_gubdb[ibin][ibin2]->Add(h_chg_reco_b[ibin][ibin2]);
            }
*/    
            h_chg_reco_gubdb[ibin][ibin2]->Rebin(10);
            h_chg_reco_gubdb[ibin][ibin2]->Scale(1./h_chg_reco_gubdb[ibin][ibin2]->Integral());

            h_chg_reco_gubdb_kappa[ibin][ibin2]->Rebin(10);
            h_chg_reco_gubdb_kappa[ibin][ibin2]->Scale(1./h_chg_reco_gubdb_kappa[ibin][ibin2]->Integral());
     
            h_chg_ref_gubdb[ibin][ibin2] = (TH1D*)h_chg_ref_g[ibin][ibin2]->Clone(Form("h_chg_ref_gubdb_%d_%d",ibin,ibin2));
            h_chg_ref_gubdb[ibin][ibin2]->Add(h_chg_ref_upb[ibin][ibin2]);
            h_chg_ref_gubdb[ibin][ibin2]->Add(h_chg_ref_downb[ibin][ibin2]);
            h_chg_ref_gubdb[ibin][ibin2]->Add(h_chg_ref_c[ibin][ibin2]);
            h_chg_ref_gubdb[ibin][ibin2]->Add(h_chg_ref_s[ibin][ibin2]);
            h_chg_ref_gubdb[ibin][ibin2]->Add(h_chg_ref_b[ibin][ibin2]);

            h_chg_reco_others[ibin][ibin2] = (TH1D*)h_chg_reco_upb[ibin][ibin2]->Clone(Form("h_chg_reco_others_%d_%d",ibin,ibin2));
            h_chg_reco_others[ibin][ibin2]->Add(h_chg_reco_downb[ibin][ibin2]);
            h_chg_reco_others[ibin][ibin2]->Add(h_chg_reco_c[ibin][ibin2]);
            h_chg_reco_others[ibin][ibin2]->Add(h_chg_reco_s[ibin][ibin2]);
            h_chg_reco_others[ibin][ibin2]->Add(h_chg_reco_b[ibin][ibin2]);

            h_chg_reco_others[ibin][ibin2]->Rebin(10);
            h_chg_reco_others[ibin][ibin2]->Scale(1./h_chg_reco_others[ibin][ibin2]->Integral());

            h_chg_reco_g_kappa[ibin][ibin2] = h_chg_corrpt_g_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_g_kappa_%d_%d",ibin,ibin2),h_chg_corrpt_g_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_g_kappa[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_reco_g_kappa[ibin][ibin2]->Rebin(10);
            h_chg_reco_g_kappa[ibin][ibin2]->Scale(1./h_chg_reco_g_kappa[ibin][ibin2]->Integral());

            h_chg_reco_upb_kappa[ibin][ibin2] = h_chg_corrpt_upb_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_upb_kappa_%d_%d",ibin,ibin2),h_chg_corrpt_upb[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_upb[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_reco_downb_kappa[ibin][ibin2] = h_chg_corrpt_downb_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_downb_kappa_%d_%d",ibin,ibin2),h_chg_corrpt_downb[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_downb[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_reco_c_kappa[ibin][ibin2] = h_chg_corrpt_c_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_c_kappa_%d_%d",ibin,ibin2),h_chg_corrpt_c[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_c[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_reco_s_kappa[ibin][ibin2] = h_chg_corrpt_s_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_s_kappa_%d_%d",ibin,ibin2),h_chg_corrpt_s[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_s[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");
            h_chg_reco_b_kappa[ibin][ibin2] = h_chg_corrpt_b_kappa[ibin][ibin2]-> ProjectionY(Form("h_chg_reco_b_kappa_%d_%d",ibin,ibin2),h_chg_corrpt_b[ibin][ibin2]->GetXaxis()->FindBin(pt_low),h_chg_corrpt_b[ibin][ibin2]->GetXaxis()->FindBin(pt_high),"");

            h_chg_reco_others_kappa[ibin][ibin2] = (TH1D*)h_chg_reco_upb_kappa[ibin][ibin2]->Clone(Form("h_chg_reco_others_kappa_%d_%d",ibin,ibin2));
            h_chg_reco_others_kappa[ibin][ibin2]->Add(h_chg_reco_downb_kappa[ibin][ibin2]);
            h_chg_reco_others_kappa[ibin][ibin2]->Add(h_chg_reco_c_kappa[ibin][ibin2]);
            h_chg_reco_others_kappa[ibin][ibin2]->Add(h_chg_reco_s_kappa[ibin][ibin2]);
            h_chg_reco_others_kappa[ibin][ibin2]->Add(h_chg_reco_b_kappa[ibin][ibin2]);

            h_chg_reco_others_kappa[ibin][ibin2]->Rebin(10);
            h_chg_reco_others_kappa[ibin][ibin2]->Scale(1./h_chg_reco_others_kappa[ibin][ibin2]->Integral());

            if(do_data){
                //scaling up and down to match data
                h_chg_reco_up[ibin][ibin2]->Scale(1./(h_chg_reco_up[ibin][ibin2]->Integral()));
                h_chg_reco_down[ibin][ibin2]->Scale(1./(h_chg_reco_down[ibin][ibin2]->Integral()));

                h_chg_reco_ud_scaled[ibin][ibin2] = (TH1D*)h_chg_reco_down[ibin][ibin2]->Clone((TString)("h_chg_reco_ud_scaled_"+cent[ibin]+"_"+trk[ibin2]));
                if(ibin==0) h_chg_reco_up[ibin][ibin2]->Scale(1.896);
                else h_chg_reco_up[ibin][ibin2]->Scale(0.882);
                h_chg_reco_ud_scaled[ibin][ibin2]->Add(h_chg_reco_up[ibin][ibin2]); 
                h_chg_reco_ud_scaled[ibin][ibin2]->Rebin(10);
                h_chg_reco_ud_scaled[ibin][ibin2]->Scale(1./h_chg_reco_ud_scaled[ibin][ibin2]->Integral());

                if(ibin==0) h_chg_reco_up[ibin][ibin2]->Scale(1./1.896);
                else h_chg_reco_up[ibin][ibin2]->Scale(1./0.882);
 
                //scaling up and down to match data
                h_chg_reco_up_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_up_kappa[ibin][ibin2]->Integral()));
                h_chg_reco_down_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_down_kappa[ibin][ibin2]->Integral()));

                h_chg_reco_ud_scaled_kappa[ibin][ibin2] = (TH1D*)h_chg_reco_down_kappa[ibin][ibin2]->Clone((TString)("h_chg_reco_ud_scaled_kappa_"+cent[ibin]+"_"+trk[ibin2]));
                if(ibin==0) h_chg_reco_up_kappa[ibin][ibin2]->Scale(1.896);
                else h_chg_reco_up_kappa[ibin][ibin2]->Scale(0.882);
                h_chg_reco_ud_scaled_kappa[ibin][ibin2]->Add(h_chg_reco_up_kappa[ibin][ibin2]); 
                h_chg_reco_ud_scaled_kappa[ibin][ibin2]->Rebin(10);
                h_chg_reco_ud_scaled_kappa[ibin][ibin2]->Scale(1./h_chg_reco_ud_scaled_kappa[ibin][ibin2]->Integral());

                if(ibin==0) h_chg_reco_up_kappa[ibin][ibin2]->Scale(1./1.896);
                else h_chg_reco_up_kappa[ibin][ibin2]->Scale(1./0.882);

            }
            else{
                //MC 
                h_chg_reco_up[ibin][ibin2]->Scale(1./(h_chg_reco[ibin][ibin2]->Integral()));
                h_chg_reco_down[ibin][ibin2]->Scale(1./(h_chg_reco[ibin][ibin2]->Integral()));

                h_chg_reco_ud_scaled[ibin][ibin2] = (TH1D*)h_chg_reco_down[ibin][ibin2]->Clone((TString)("h_chg_reco_ud_scaled_"+cent[ibin]+"_"+trk[ibin2]));
                h_chg_reco_ud_scaled[ibin][ibin2]->Add(h_chg_reco_up[ibin][ibin2]); 
                h_chg_reco_ud_scaled[ibin][ibin2]->Rebin(10);
                h_chg_reco_ud_scaled[ibin][ibin2]->Scale(1./h_chg_reco_ud_scaled[ibin][ibin2]->Integral());

                h_chg_reco_up[ibin][ibin2]->Scale(1./(h_chg_reco_up[ibin][ibin2]->Integral()));
                h_chg_reco_down[ibin][ibin2]->Scale(1./(h_chg_reco_down[ibin][ibin2]->Integral()));

                //MC 
                h_chg_reco_up_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_kappa[ibin][ibin2]->Integral()));
                h_chg_reco_down_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_kappa[ibin][ibin2]->Integral()));

                h_chg_reco_ud_scaled_kappa[ibin][ibin2] = (TH1D*)h_chg_reco_down_kappa[ibin][ibin2]->Clone((TString)("h_chg_reco_ud_scaled_kappa_"+cent[ibin]+"_"+trk[ibin2]));
                h_chg_reco_ud_scaled_kappa[ibin][ibin2]->Add(h_chg_reco_up_kappa[ibin][ibin2]); 
                h_chg_reco_ud_scaled_kappa[ibin][ibin2]->Rebin(10);
                h_chg_reco_ud_scaled_kappa[ibin][ibin2]->Scale(1./h_chg_reco_ud_scaled_kappa[ibin][ibin2]->Integral());

                h_chg_reco_up_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_up_kappa[ibin][ibin2]->Integral()));
                h_chg_reco_down_kappa[ibin][ibin2]->Scale(1./(h_chg_reco_down_kappa[ibin][ibin2]->Integral()));
            }

            h_chg_reco_up[ibin][ibin2]->Rebin(10);
            h_chg_reco_up[ibin][ibin2]->Scale(1./h_chg_reco_up[ibin][ibin2]->Integral());
            h_chg_reco_down[ibin][ibin2]->Rebin(10);
            h_chg_reco_down[ibin][ibin2]->Scale(1./h_chg_reco_down[ibin][ibin2]->Integral());

            h_chg_reco_up_kappa[ibin][ibin2]->Rebin(10);
            h_chg_reco_up_kappa[ibin][ibin2]->Scale(1./h_chg_reco_up_kappa[ibin][ibin2]->Integral());
            h_chg_reco_down_kappa[ibin][ibin2]->Rebin(10);
            h_chg_reco_down_kappa[ibin][ibin2]->Scale(1./h_chg_reco_down_kappa[ibin][ibin2]->Integral());

            //replicas for plotting 
            sprintf(saythis,"h_chg_data_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_data_cop[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
            h_chg_data_cop[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_data_cop_sysm_cent%d_trk%d",ibin,ibin2);
            h_chg_data_cop_sysm[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
            h_chg_data_cop_sysm[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_data_cop_sysm_kappa_cent%d_trk%d",ibin,ibin2);
            h_chg_data_kappa_cop_sysm[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
            h_chg_data_kappa_cop_sysm[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_data_trkrel_cent%d_trk%d",ibin,ibin2);
            h_chg_data_trkrel[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
            h_chg_data_trkrel[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_reco_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_reco_cop[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
            h_chg_reco_cop[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_reco_up_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_reco_up_cop[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
            h_chg_reco_up_cop[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_reco_down_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_reco_down_cop[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
            h_chg_reco_down_cop[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_reco_gubdb_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_reco_gubdb_cop[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
            h_chg_reco_gubdb_cop[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_reco_g_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_reco_g_cop[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
            h_chg_reco_g_cop[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_reco_others_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_reco_others_cop[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
            h_chg_reco_others_cop[ibin][ibin2]->Sumw2();

            //replicas for plotting 
            sprintf(saythis,"h_chg_data_kappa_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_data_kappa_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
            h_chg_data_kappa_cop[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_reco_kappa_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_reco_kappa_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
            h_chg_reco_kappa_cop[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_reco_up_kappa_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_reco_up_kappa_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
            h_chg_reco_up_kappa_cop[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_reco_down_kappa_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_reco_down_kappa_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
            h_chg_reco_down_kappa_cop[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_reco_gubdb_kappa_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_reco_gubdb_kappa_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
            h_chg_reco_gubdb_kappa_cop[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_reco_g_kappa_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_reco_g_kappa_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
            h_chg_reco_g_kappa_cop[ibin][ibin2]->Sumw2();
            sprintf(saythis,"h_chg_reco_others_kappa_cop_cent%d_trk%d",ibin,ibin2);
            h_chg_reco_others_kappa_cop[ibin][ibin2] = new TH1D(saythis,"",500,-2.555,2.445);
            h_chg_reco_others_kappa_cop[ibin][ibin2]->Sumw2();

            for(int ibin3=1;ibin3<(h_chg_data[ibin][ibin2]->GetXaxis()->GetNbins())+1;ibin3++){
              if(do_ptcut_fitting){
                  h_chg_data_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_data[ibin][ibin2]->GetBinContent(ibin3));
                  h_chg_data_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_data[ibin][ibin2]->GetBinError(ibin3));
                  if(h_chg_data_trkeff[ibin][ibin2][5]->GetBinContent(ibin3) > 0){
                      if(ibin==0) h_chg_data_trkrel[ibin][ibin2]->SetBinContent(ibin3,fabs(h_chg_data_trkeff[ibin][ibin2][1]->GetBinContent(ibin3) - h_chg_data_trkeff[ibin][ibin2][9]->GetBinContent(ibin3))/h_chg_data[ibin][ibin2]->GetBinContent(ibin3)/2.);
                      else h_chg_data_trkrel[ibin][ibin2]->SetBinContent(ibin3,fabs(h_chg_data_trkeff[ibin][ibin2][0]->GetBinContent(ibin3) - h_chg_data_trkeff[ibin][ibin2][10]->GetBinContent(ibin3))/h_chg_data[ibin][ibin2]->GetBinContent(ibin3)/2.);
                      h_chg_data_trkrel[ibin][ibin2]->SetBinError(ibin3,0);
                  }
                  else{
                      h_chg_data_trkrel[ibin][ibin2]->SetBinContent(ibin3,0);
                  }

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

                  h_chg_data_cop[ibin][ibin2]->SetLineColor(kBlack); h_chg_data_cop[ibin][ibin2]->SetFillColor(kGray); h_chg_data_cop[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_data_cop[ibin][ibin2]->SetMarkerStyle(1);
                  h_chg_reco_cop[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_cop[ibin][ibin2]->SetFillColor(kGray); h_chg_reco_cop[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_reco_cop[ibin][ibin2]->SetMarkerStyle(20);
                  h_chg_reco_up_cop[ibin][ibin2]->SetLineColor(kRed+2); h_chg_reco_up_cop[ibin][ibin2]->SetFillColor(kRed+2); h_chg_reco_up_cop[ibin][ibin2]->SetMarkerColor(kRed+2); h_chg_reco_up_cop[ibin][ibin2]->SetMarkerStyle(0);    
                  h_chg_reco_down_cop[ibin][ibin2]->SetLineColor(kYellow-7); h_chg_reco_down_cop[ibin][ibin2]->SetFillColorAlpha(kYellow-7,0.7); h_chg_reco_down_cop[ibin][ibin2]->SetMarkerColor(kYellow-7); h_chg_reco_down_cop[ibin][ibin2]->SetMarkerStyle(0); 
                  h_chg_reco_g_cop[ibin][ibin2]->SetLineColor(kAzure+1); h_chg_reco_g_cop[ibin][ibin2]->SetFillColorAlpha(kAzure+1,0.7); h_chg_reco_g_cop[ibin][ibin2]->SetMarkerColor(kAzure+1); h_chg_reco_g_cop[ibin][ibin2]->SetMarkerStyle(0);      
                  h_chg_reco_others_cop[ibin][ibin2]->SetLineColor(kGreen-2); h_chg_reco_others_cop[ibin][ibin2]->SetFillColorAlpha(kGreen-2,0.7); h_chg_reco_others_cop[ibin][ibin2]->SetMarkerColor(kGreen-2); h_chg_reco_others_cop[ibin][ibin2]->SetMarkerStyle(0);
/*
                  h_chg_reco_gubdb_cop[ibin][ibin2]->Rebin(2);
                  h_chg_reco_g_cop[ibin][ibin2]->Rebin(2);
                  h_chg_reco_others_cop[ibin][ibin2]->Rebin(2);
                  h_chg_reco_down_cop[ibin][ibin2]->Rebin(2);
                  h_chg_reco_cop[ibin][ibin2]->Rebin(2);
                  h_chg_data_cop[ibin][ibin2]->Rebin(2);
                  h_chg_reco_up_cop[ibin][ibin2]->Rebin(2);
*/

                  h_chg_data_kappa_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_data_kappa[ibin][ibin2]->GetBinContent(ibin3));
                  h_chg_data_kappa_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_data_kappa[ibin][ibin2]->GetBinError(ibin3));
                  h_chg_reco_up_kappa_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_up_kappa[ibin][ibin2]->GetBinContent(ibin3));
                  h_chg_reco_up_kappa_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_up_kappa[ibin][ibin2]->GetBinError(ibin3));
                  h_chg_reco_kappa_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_kappa[ibin][ibin2]->GetBinContent(ibin3));
                  h_chg_reco_kappa_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_kappa[ibin][ibin2]->GetBinError(ibin3));
                  h_chg_reco_down_kappa_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_down_kappa[ibin][ibin2]->GetBinContent(ibin3));
                  h_chg_reco_down_kappa_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_down_kappa[ibin][ibin2]->GetBinError(ibin3));
                  h_chg_reco_g_kappa_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_g_kappa[ibin][ibin2]->GetBinContent(ibin3));
                  h_chg_reco_g_kappa_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_g_kappa[ibin][ibin2]->GetBinError(ibin3));
                  h_chg_reco_gubdb_kappa_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_gubdb_kappa[ibin][ibin2]->GetBinContent(ibin3));
                  h_chg_reco_gubdb_kappa_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_gubdb_kappa[ibin][ibin2]->GetBinError(ibin3));
                  h_chg_reco_others_kappa_cop[ibin][ibin2]->SetBinContent(ibin3,h_chg_reco_others_kappa[ibin][ibin2]->GetBinContent(ibin3));
                  h_chg_reco_others_kappa_cop[ibin][ibin2]->SetBinError(ibin3,h_chg_reco_others_kappa[ibin][ibin2]->GetBinError(ibin3));

                  h_chg_data_kappa_cop[ibin][ibin2]->SetLineColor(kBlack); h_chg_data_kappa_cop[ibin][ibin2]->SetFillColor(kGray); h_chg_data_kappa_cop[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_data_kappa_cop[ibin][ibin2]->SetMarkerStyle(1);
                  h_chg_reco_kappa_cop[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_kappa_cop[ibin][ibin2]->SetFillColor(kGray); h_chg_reco_kappa_cop[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_reco_kappa_cop[ibin][ibin2]->SetMarkerStyle(20);
                  h_chg_reco_up_kappa_cop[ibin][ibin2]->SetLineColor(kRed+2); h_chg_reco_up_kappa_cop[ibin][ibin2]->SetFillColor(kRed+2); h_chg_reco_up_kappa_cop[ibin][ibin2]->SetMarkerColor(kRed+2); h_chg_reco_up_kappa_cop[ibin][ibin2]->SetMarkerStyle(0);    
                  h_chg_reco_down_kappa_cop[ibin][ibin2]->SetLineColor(kYellow-7); h_chg_reco_down_kappa_cop[ibin][ibin2]->SetFillColorAlpha(kYellow-7,0.7); h_chg_reco_down_kappa_cop[ibin][ibin2]->SetMarkerColor(kYellow-7); h_chg_reco_down_kappa_cop[ibin][ibin2]->SetMarkerStyle(0); 
                  h_chg_reco_g_kappa_cop[ibin][ibin2]->SetLineColor(kAzure+1); h_chg_reco_g_kappa_cop[ibin][ibin2]->SetFillColorAlpha(kAzure+1,0.7); h_chg_reco_g_kappa_cop[ibin][ibin2]->SetMarkerColor(kAzure+1); h_chg_reco_g_kappa_cop[ibin][ibin2]->SetMarkerStyle(0);      
                  h_chg_reco_others_kappa_cop[ibin][ibin2]->SetLineColor(kGreen-2); h_chg_reco_others_kappa_cop[ibin][ibin2]->SetFillColorAlpha(kGreen-2,0.7); h_chg_reco_others_kappa_cop[ibin][ibin2]->SetMarkerColor(kGreen-2); h_chg_reco_others_kappa_cop[ibin][ibin2]->SetMarkerStyle(0);

              }
            }

        }

    }



}

TLine *tl;

void drawline(double x1, double y1, double x2, double y2){
    tl = new TLine(x1,y1,x2,y2);  
    tl->SetLineStyle(2);
    tl->Draw("same");
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

void drawgraphPad(){

  pad_titleg = new TPad("pad_titleg", "",0.,0.93,1.,1.);
  pad_labelg = new TPad("pad_labelg", "",0.,0.,0.02,0.95);
  pada = new TPad("pada", "",0.02,0.,1.,0.98);
  pada->Divide(5,1,0);  
  pada->Draw();
  pad_titleg->Draw();
  pad_labelg->Draw();

}

void plot_jetcharge_arc(bool do_data = false){

    get_histos(do_data);

    //plotting texts
    char g_frac[50],g_error[50],chi2ndf[50], up_frac[50], down_frac[50];

    TString tmp;
    TString pythia = "MC";
    TString ppdata = "Data";
    TString data = "pp 27.4 pb^{-1} (5.02 TeV)  PbPb 404 #mub^{-1} (5.02 TeV)";
    //TString data = "pp 27.4 pb^{-1} (5.02 TeV)";
    TString cms = "CMS Preliminary";
    TString gen_tex = "MC Gen";
    TString chg_label = "1/N_{jets}";
    TString ratio = "Data / MC";
    TString fraction = "gluon like jet fraction";
    TString gluon_like_fraction = "gluon like jet fraction";
    TString y_label =  "1/N_{jets} (dN/dQ) [1/e]";
    TString datafit =  "Data / Fit";
    TString tag_ratio =  "Ratio";
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

    TLatex *tx, *tx1, *tx2, *tx3, *tx4, *tx5;

    TLine *tl_q_reco[nCBins]; TLine *tl_g_reco[nCBins]; TLine *tl_u_reco[nCBins]; TLine *tl_d_reco[nCBins];
    TLine *tl_q_ref[nCBins]; TLine *tl_g_ref[nCBins]; TLine *tl_u_ref[nCBins]; TLine *tl_d_ref[nCBins];
    TLine *tl_q_ref_kappa[nCBins]; TLine *tl_g_ref_kappa[nCBins]; TLine *tl_u_ref_kappa[nCBins]; TLine *tl_d_ref_kappa[nCBins];
    
    //direct frac from MC 
    Double_t gluon_reco[nCBins],quark_reco[nCBins],up_reco[nCBins],upbar_reco[nCBins],down_reco[nCBins],downbar_reco[nCBins],gluon_ubar_dbar_reco[nCBins], charm_reco[nCBins], strange_reco[nCBins], bottom_reco[nCBins];
    Double_t gluon_gen[nCBins],quark_gen[nCBins], down_gen[nCBins], up_gen[nCBins], gluon_ubar_dbar_gen[nCBins], data_int[nCBins], mc_int[nCBins];
    Double_t ubdbcsb_reco[nCBins],gluon_ubdbcsb[nCBins], u_ud[nCBins], d_ud[nCBins], gluon_corr_factor[nCBins], quark_corr_factor[nCBins];

    for(int ibin=0;ibin<nCBins;ibin++){
        int low_bin = h_reco_MC[ibin]->FindBin(pt_low); int high_bin = h_reco_MC[ibin]->FindBin(pt_high); 

        gluon_gen[ibin] = h_ref_MC_g[ibin]->Integral(low_bin,high_bin)/h_ref_MC[ibin]->Integral(low_bin,high_bin);
        quark_gen[ibin] = h_ref_MC_q[ibin]->Integral(low_bin,high_bin)/h_ref_MC[ibin]->Integral(low_bin,high_bin);
        up_gen[ibin] = h_ref_MC_up[ibin]->Integral(low_bin,high_bin)/h_ref_MC[ibin]->Integral(low_bin,high_bin);
        down_gen[ibin] = h_ref_MC_d[ibin]->Integral(low_bin,high_bin)/h_ref_MC[ibin]->Integral(low_bin,high_bin);
        gluon_ubar_dbar_gen[ibin] = (h_ref_MC_g[ibin]->Integral(low_bin,high_bin)+h_ref_MC_upbar[ibin]->Integral(low_bin,high_bin)+h_ref_MC_dbar[ibin]->Integral(low_bin,high_bin))/h_ref_MC[ibin]->Integral(low_bin,high_bin);

        gluon_reco[ibin] = h_reco_MC_g[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
        quark_reco[ibin] = h_reco_MC_q[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
        up_reco[ibin] = h_reco_MC_up[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
        down_reco[ibin] = h_reco_MC_d[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
        upbar_reco[ibin] = h_reco_MC_upbar[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
        downbar_reco[ibin] = h_reco_MC_dbar[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
        charm_reco[ibin] = h_reco_MC_c[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
        strange_reco[ibin] = h_reco_MC_s[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
        bottom_reco[ibin] = h_reco_MC_b[ibin]->Integral(low_bin,high_bin)/h_reco_MC[ibin]->Integral(low_bin,high_bin);
        gluon_ubar_dbar_reco[ibin] = (h_reco_MC_g[ibin]->Integral(low_bin,high_bin)+h_reco_MC_upbar[ibin]->Integral(low_bin,high_bin)+h_reco_MC_dbar[ibin]->Integral(low_bin,high_bin)+h_reco_MC_c[ibin]->Integral(low_bin,high_bin)+h_reco_MC_s[ibin]->Integral(low_bin,high_bin)+h_reco_MC_b[ibin]->Integral(low_bin,high_bin))/h_reco_MC[ibin]->Integral(low_bin,high_bin);

        data_int[ibin] = h_reco_data[ibin]->Integral(low_bin,high_bin);
        mc_int[ibin] = h_reco_MC[ibin]->Integral(low_bin,high_bin);

        ubdbcsb_reco[ibin] = upbar_reco[ibin]+downbar_reco[ibin]+charm_reco[ibin]+strange_reco[ibin]+bottom_reco[ibin];
        gluon_ubdbcsb[ibin] = gluon_reco[ibin]/(gluon_reco[ibin]+ubdbcsb_reco[ibin]);
        u_ud[ibin] = up_reco[ibin]/(up_reco[ibin]+down_reco[ibin]);
        d_ud[ibin] = down_reco[ibin]/(up_reco[ibin]+down_reco[ibin]);

        gluon_corr_factor[ibin] = gluon_gen[ibin]/gluon_reco[ibin];
        quark_corr_factor[ibin] = quark_gen[ibin]/quark_reco[ibin]; 

        cout<<"cent bin: "<<ibin<<endl;
        cout<<"up: "<<up_reco[ibin]<<"  down:  "<<down_reco[ibin]<<"  g_ub_db:  "<<gluon_ubar_dbar_reco[ibin]<<endl;
        //cout<<"charm: "<<charm_reco[ibin]<<"  strange:  "<<strange_reco[ibin]<<"  bottom:  "<<bottom_reco[ibin]<<endl;
        cout<<"total: "<<gluon_ubar_dbar_reco[ibin]+up_reco[ibin]+down_reco[ibin]<<endl;
        cout<<"gluongen / gluonreco"<<gluon_corr_factor[ibin]<<".  quarkgen / quarkreco"<<quark_corr_factor[ibin]<<endl;

        //0 trk jets 
        for(int ibin2=0; ibin2<ntrkbins; ibin2++){

            quark_0trk_int[ibin][ibin2] = h_gen_full_0trks_q[ibin][ibin2]->Integral(low_bin,high_bin)/h_gen_full_0trks[ibin][ibin2]->Integral(low_bin,high_bin);
            gluon_0trk_int[ibin][ibin2] = h_gen_full_0trks_g[ibin][ibin2]->Integral(low_bin,high_bin)/h_gen_full_0trks[ibin][ibin2]->Integral(low_bin,high_bin);

            data_0trk_int[ibin][ibin2] = h_reco_data_0trks[ibin][ibin2]->Integral(low_bin,high_bin);
            data_no0trk_int[ibin][ibin2] = (data_int[ibin] - data_0trk_int[ibin][ibin2])/data_int[ibin];
            mc_0trk_int[ibin][ibin2] = h_reco_corr_0trks[ibin][ibin2]->Integral(low_bin,high_bin);
            mc_no0trk_int[ibin][ibin2] = (mc_int[ibin] - mc_0trk_int[ibin][ibin2])/mc_int[ibin];
            quark_data_0trk_int[ibin][ibin2] = quark_0trk_int[ibin][ibin2]*data_0trk_int[ibin][ibin2]/data_int[ibin];
            quark_mc_0trk_int[ibin][ibin2] = quark_0trk_int[ibin][ibin2]*mc_0trk_int[ibin][ibin2]/mc_int[ibin];
            gluon_data_0trk_int[ibin][ibin2] = (1-quark_0trk_int[ibin][ibin2])*data_0trk_int[ibin][ibin2]/data_int[ibin];
        }

        tl_g_ref[ibin] = new TLine(1.,gluon_gen[ibin],6.,gluon_gen[ibin]); tl_g_ref[ibin]->SetLineStyle(2); tl_g_ref[ibin]->SetLineColor(kRed);
        tl_g_ref_kappa[ibin] = new TLine(0.,gluon_gen[ibin],1.,gluon_gen[ibin]); tl_g_ref_kappa[ibin]->SetLineStyle(2); tl_g_ref_kappa[ibin]->SetLineColor(kRed);

        tl_g_reco[ibin] = new TLine(1.,gluon_reco[ibin],6.,gluon_reco[ibin]); tl_g_reco[ibin]->SetLineStyle(2); tl_g_reco[ibin]->SetLineColor(kRed);
    }

    //setsystematics

    Double_t trkeff_up_sys, trkeff_dn_sys, trkeff_pos_sys, jes_up_sys, jes_dn_sys;
    Double_t trkeff_up_kappa_sys, trkeff_dn_kappa_sys, trkeff_pos_kappa_sys;
    Double_t trkeff_sys, jes_res_sys, total_err;
    Double_t trkeff_kappa_sys, jes_res_kappa_sys, total_kappa_err;
    Double_t unf_err_ptcut, unf_err_kappa;

    for(int ibin = 0; ibin < nCBins; ibin++){
      for(int ibin2 = 0; ibin2 < ntrkbins; ibin2++){
        for(int ibin3 = 1; ibin3 < h_chg_data[0][0]->GetNbinsX()+1; ibin3++){
          if(do_ptcut_fitting && do_systematics){
            //pt cut
            trkeff_up_sys = fabs(h_chg_data[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_trkinc[ibin][ibin2]->GetBinContent(ibin3));
            trkeff_dn_sys = fabs(h_chg_data[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_trkdec[ibin][ibin2]->GetBinContent(ibin3));
            trkeff_sys = 0.;
            trkeff_sys = (trkeff_up_sys+trkeff_dn_sys)/2.;
            trkeff_pos_sys = 0.;   
            trkeff_pos_sys = fabs(h_chg_data[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_svd[ibin][ibin2]->GetBinContent(ibin3));       
            jes_res_sys = 0.;
            jes_res_sys = fabs(h_chg_data[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_jer[ibin][ibin2]->GetBinContent(ibin3));
            unf_err_ptcut = 0.;
            unf_err_ptcut = fabs(h_chg_data[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_glusamp[ibin][ibin2]->GetBinContent(ibin3));
            total_err = sqrt(pow(unf_err_ptcut,2)+pow(trkeff_pos_sys,2)+pow(trkeff_sys,2)+pow(jes_res_sys,2)+pow(h_chg_data[ibin][ibin2]->GetBinError(ibin3),2));
            h_chg_data[ibin][ibin2]->SetBinError(ibin3, total_err);
            h_chg_data_glusamp[ibin][ibin2]->SetBinError(ibin3, total_err);
            //h_chg_data_trkinc[ibin][ibin2]->SetBinError(ibin3, total_err);
            //h_chg_data_trkdec[ibin][ibin2]->SetBinError(ibin3, total_err);
            //h_chg_data_svd[ibin][ibin2]->SetBinError(ibin3, total_err);
            for(int ibin4=0; ibin4<trkeffbins;ibin4++){
                h_chg_data_trkeff[ibin][ibin2][ibin4]->SetBinError(ibin3, total_err);
            }
          }
          if(do_kappa_fitting && do_systematics){
            //kappa
            trkeff_up_kappa_sys = fabs(h_chg_data_kappa[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_kappa_trkinc[ibin][ibin2]->GetBinContent(ibin3));
            trkeff_dn_kappa_sys = fabs(h_chg_data_kappa[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_kappa_trkdec[ibin][ibin2]->GetBinContent(ibin3));
            trkeff_kappa_sys = 0.;
            trkeff_kappa_sys = (trkeff_up_kappa_sys+trkeff_dn_kappa_sys)/2.;
            trkeff_pos_kappa_sys = 0.;   
            trkeff_pos_kappa_sys = fabs(h_chg_data_kappa[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_kappa_trkposneg[ibin][ibin2]->GetBinContent(ibin3));       
            jes_res_kappa_sys = 0.;
            //jes_res_sys = 0.02*h_chg_data[ibin][ibin2]->GetBinContent(ibin3);
            jes_res_kappa_sys = fabs(h_chg_data_kappa[ibin][ibin2]->GetBinContent(ibin3) - h_chg_data_kappa_jer[ibin][ibin2]->GetBinContent(ibin3));
            total_kappa_err = sqrt(pow(trkeff_pos_kappa_sys,2)+pow(trkeff_kappa_sys,2)+pow(jes_res_kappa_sys,2)+pow(h_chg_data_kappa[ibin][ibin2]->GetBinError(ibin3),2));
            h_chg_data_kappa[ibin][ibin2]->SetBinError(ibin3, total_kappa_err);
            //h_chg_data_kappa_trkinc[ibin][ibin2]->SetBinError(ibin3, total_kappa_err);
            //h_chg_data_kappa_trkdec[ibin][ibin2]->SetBinError(ibin3, total_kappa_err);
            //h_chg_data_kappa_trkposneg[ibin][ibin2]->SetBinError(ibin3, total_kappa_err);
          }
        }
      }
    }

///updated systematics for plotting
    for(int ibin=0; ibin<nCBins; ibin++){
        for(int ibin2=0; ibin2<ntrkbins; ibin2++){
            for(int ibin3=1;ibin3<(h_chg_data[ibin][ibin2]->GetXaxis()->GetNbins())+1;ibin3++){
                if(do_ptcut_fitting){
                    h_chg_data_cop_sysm[ibin][ibin2]->SetBinContent(ibin3,h_chg_data[ibin][ibin2]->GetBinContent(ibin3));
                    h_chg_data_cop_sysm[ibin][ibin2]->SetBinError(ibin3,h_chg_data[ibin][ibin2]->GetBinError(ibin3));
                    h_chg_data_cop_sysm[ibin][ibin2]->SetLineColor(kBlack); h_chg_data_cop_sysm[ibin][ibin2]->SetFillColorAlpha(kBlack,0.4); h_chg_data_cop_sysm[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_data_cop_sysm[ibin][ibin2]->SetMarkerStyle(20); h_chg_data_cop_sysm[ibin][ibin2]->SetMarkerSize(0.5);
                }
                if(do_kappa_fitting){
                    h_chg_data_kappa_cop_sysm[ibin][ibin2]->SetBinContent(ibin3,h_chg_data_kappa[ibin][ibin2]->GetBinContent(ibin3));
                    h_chg_data_kappa_cop_sysm[ibin][ibin2]->SetBinError(ibin3,h_chg_data_kappa[ibin][ibin2]->GetBinError(ibin3));
                    h_chg_data_kappa_cop_sysm[ibin][ibin2]->SetLineColor(kBlack); h_chg_data_kappa_cop_sysm[ibin][ibin2]->SetFillColorAlpha(kBlack,0.4); h_chg_data_kappa_cop_sysm[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_data_kappa_cop_sysm[ibin][ibin2]->SetMarkerStyle(20); h_chg_data_kappa_cop_sysm[ibin][ibin2]->SetMarkerSize(0.5);
                }
            }
        }
    }

    //templates
    //reco

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

    ////up/down templates
    f_chg_ud_data[0][0] = new TF1("f_chg_ud_data_cent0_trkpt2",f_ud_reco_cent0_pt2,-2.,2.,2);
    f_chg_ud_data[1][0] = new TF1("f_chg_ud_data_cent1_trkpt2",f_ud_reco_cent1_pt2,-2.,2.,2);
    f_chg_ud_data[2][0] = new TF1("f_chg_ud_data_cent2_trkpt2",f_ud_reco_cent2_pt2,-2.,2.,2);
    f_chg_ud_data[3][0] = new TF1("f_chg_ud_data_cent3_trkpt2",f_ud_reco_cent3_pt2,-2.,2.,2);
    f_chg_ud_data[4][0] = new TF1("f_chg_ud_data_cent4_trkpt2",f_ud_reco_cent4_pt2,-2.,2.,2);

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
        h_eta_MC_g[ibin]->GetYaxis()->SetRangeUser(0.,0.12);
        h_eta_MC_g[ibin]->SetLineColor(kRed-2); h_eta_MC_g[ibin]->SetMarkerColor(kRed-2);
        h_eta_MC_g[ibin]->SetMarkerStyle(4); h_eta_MC_g[ibin]->Draw("e0 same");
        h_eta_MC[ibin]->SetLineColor(kBlack); h_eta_MC[ibin]->SetMarkerColor(kBlack);
        h_eta_MC[ibin]->SetMarkerStyle(20); h_eta_MC[ibin]->Draw("e0 same");
        h_eta_MC_up[ibin]->SetLineColor(kBlue-4); h_eta_MC_up[ibin]->SetMarkerColor(kBlue-4);
        h_eta_MC_up[ibin]->SetMarkerStyle(22); h_eta_MC_up[ibin]->Draw("e0 same");
        h_eta_MC_d[ibin]->SetLineColor(kOrange+1); h_eta_MC_d[ibin]->SetMarkerColor(kOrange+1);
        h_eta_MC_d[ibin]->SetMarkerStyle(23); h_eta_MC_d[ibin]->Draw("e0 same");
        h_eta_data[ibin]->SetLineColor(kGreen-2); h_eta_data[ibin]->SetMarkerColor(kGreen-2);
        h_eta_data[ibin]->SetMarkerStyle(5); h_eta_data[ibin]->Draw("e0 same");

        if(ibin==0){
          TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
          legend ->SetLineColor(kWhite);
          legend ->AddEntry(h_eta_data[ibin], "inclusive jets (data)", "lepf");
          legend ->AddEntry(h_eta_MC[ibin], "inclusive jets (MC)", "lepf");
          legend ->AddEntry(h_eta_MC_g[ibin], "gluon jets (MC)", "lepf");
          legend ->AddEntry(h_eta_MC_up[ibin], "up quark jets (MC)", "lepf");
          legend ->AddEntry(h_eta_MC_d[ibin], "down quark jets (MC)", "lepf");
          legend ->Draw("same");
        } 

        TString tmp1 = cent_tag_pp_MC[ibin];
        tx1 = new TLatex(); tx1->SetTextSize(.12);
        tx1->DrawLatexNDC(0.35,0.85, tmp1);
    }

    if(!do_data){
    ///two parameter ptcut fitting
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
            h_chg_reco[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
            h_chg_reco[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);

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
          
            h_chg_reco[ibin][ibin2]->Fit(f_chg_qg_data[ibin][ibin2],"Q N M R","sames",-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
          
            if(ibin==0) c_chg_qg_MC_reco[ibin2]->cd(6);
            else c_chg_qg_MC_reco[ibin2]->cd(11-ibin);   

            h_chg_reco_fit_total[ibin][ibin2] = (TH1D*)h_chg_reco_ud_scaled[ibin][ibin2]->Clone("h_chg_reco_fit_total");
            h_chg_reco_fit_total[ibin][ibin2]->Scale(f_chg_qg_data[ibin][ibin2]->GetParameter(0));
            
            h_chg_reco_gubdb_afterfit[ibin][ibin2] = (TH1D*)h_chg_reco_gubdb[ibin][ibin2]->Clone("h_chg_reco_gubdb_afterfit");
            h_chg_reco_gubdb_afterfit[ibin][ibin2]->Scale(1. - f_chg_qg_data[ibin][ibin2]->GetParameter(0));
            h_chg_reco_fit_total[ibin][ibin2]->Add(h_chg_reco_gubdb_afterfit[ibin][ibin2]); 

            h_chg_reco_ratio_2param[ibin][ibin2] = (TH1D*)h_chg_reco[ibin][ibin2]->Clone("h_chg_reco_ratio");
            h_chg_reco_ratio_2param[ibin][ibin2]->Divide(h_chg_reco_fit_total[ibin][ibin2]);
            //h_chg_reco_ratio_2param[ibin][ibin2]->Rebin(10); h_chg_reco_ratio_2param[ibin][ibin2]->Scale(1./10.);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetTitle("reco jet charge");
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetTitle("Reco / Fit");
            h_cosm(h_chg_reco_ratio_2param[ibin][ibin2]);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetRangeUser(0.8,1.2);
            h_chg_reco_ratio_2param[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerColor(kBlack);
            h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerStyle(1); h_chg_reco_ratio_2param[ibin][ibin2]->Draw("e0 same");
            drawline(-2.,1.,2.,1.);
            sprintf(g_frac, "g+others (MC/Fit) %.3f",(1.-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)))/gluon_ubar_dbar_reco[ibin]);
            sprintf(g_error, "#pm %.3f",f_chg_qg_data[ibin][ibin2]->GetParError(0));
            sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_qg_data[ibin][ibin2]->GetChisquare()/f_chg_qg_data[ibin][ibin2]->GetNDF());
            sprintf(up_frac, "up+down (MC/Fit) %.3f",(f_chg_qg_data[ibin][ibin2]->GetParameter(0)/(up_reco[ibin]+down_reco[ibin])));
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

        TCanvas *c_chg_udg_MC_MC_rebin[ntrkbins];

        for(int ibin2=0;ibin2<ntrkbins;ibin2++){
            c_chg_udg_MC_MC_rebin[ibin2] = new TCanvas((TString)("c_chg_udg_MC_MC_rebin"+trk[ibin2]),(TString)("c_chg_udg_MC_MC_rebin"+trk[ibin2]),1500,450);

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
              h_chg_reco_cop[ibin][ibin2]->GetXaxis()->SetTitle("jet charge");
              h_cosm(h_chg_reco_cop[ibin][ibin2]);
              h_chg_reco_cop[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
              h_chg_reco_cop[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.7);
              h_chg_reco_cop[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.2);

              h_chg_reco_cop[ibin][ibin2]->Draw("e0 same");

              h_chg_reco_up_cop[ibin][ibin2]->Scale(0.65*f_chg_qg_data[ibin][ibin2]->GetParameter(0));
              h_chg_reco_down_cop[ibin][ibin2]->Scale(0.35*f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 

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
                legendx ->AddEntry(h_chg_reco_cop[0][ibin2], "MC (Gen)", "lep");
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

                h_chg_reco_ratio_2param[ibin][ibin2]->Draw("same");
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
                h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetLabelSize(0.12);
                h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
                h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetRangeUser(0.5,1.5);
                h_chg_reco_ratio_2param[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerColor(kBlack);
                h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerStyle(21); h_chg_reco_ratio_2param[ibin][ibin2]->Draw("e0 same");

                drawline(-fit_limit_ptcut[ibin2],1.,fit_limit_ptcut[ibin2],1.);
              }
              
            }
         }
      }

    ///two parameter kappa fitting
      if(do_kappa_fitting){
        TCanvas *c_chg_qg_MC_reco_kappa[ntrkbins];
        for(int ibin2=0;ibin2<ntrkbins;ibin2++){
          c_chg_qg_MC_reco_kappa[ibin2] = new TCanvas((TString)("c_chg_qg_MC_reco_kappa"+trk[ibin2]),(TString)("c_chg_qg_MC_reco_kappa"+trk[ibin2]),1800,600);
          c_chg_qg_MC_reco_kappa[ibin2]->Divide(5,2,0);
          gStyle->SetOptStat(0);
          for(int ibin=0;ibin<nCBins;ibin++){
            if(ibin==0) c_chg_qg_MC_reco_kappa[ibin2]->cd(1);
            else c_chg_qg_MC_reco_kappa[ibin2]->cd(6-ibin);   
            h_chg_reco_kappa[ibin][ibin2]->GetXaxis()->SetTitle("reco jet chg");
            h_chg_reco_kappa[ibin][ibin2]->GetYaxis()->SetTitle("# of jets");
            h_cosm(h_chg_reco_kappa[ibin][ibin2]);
            h_chg_reco_kappa[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_kappa[ibin2],fit_limit_kappa[ibin2]);
            h_chg_reco_kappa[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);

            h_chg_reco_kappa[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_kappa[ibin][ibin2]->SetMarkerColor(kBlack); h_chg_reco_kappa[ibin][ibin2]->SetMarkerStyle(2);
            h_chg_reco_gubdb_kappa[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_gubdb_kappa[ibin][ibin2]->SetMarkerColor(kRed); h_chg_reco_gubdb_kappa[ibin][ibin2]->SetMarkerStyle(2);
            h_chg_reco_ud_scaled_kappa[ibin][ibin2]->SetLineColor(kBlue); h_chg_reco_ud_scaled_kappa[ibin][ibin2]->SetMarkerColor(kBlue); h_chg_reco_ud_scaled_kappa[ibin][ibin2]->SetMarkerStyle(2);

            h_chg_reco_kappa[ibin][ibin2]->Draw("e0 same");
            h_chg_reco_ud_scaled_kappa[ibin][ibin2]->Draw("e0 same");
            h_chg_reco_gubdb_kappa[ibin][ibin2]->Draw("e0 same"); 

            TString tmp1 = cent_tag_pp_MC[ibin];
            tx1 = new TLatex(); tx1->SetTextSize(.1);
            tx1->DrawLatexNDC(0.15,0.85, tmp1);
          
            TString tmp = kappa_tag[ibin2];
            tx = new TLatex(); tx->SetTextSize(.08);
            tx->DrawLatexNDC(0.15,0.65, tmp);
          
            h_chg_reco_kappa[ibin][ibin2]->Fit(f_chg_qg_kappa[ibin][ibin2],"Q N M R","sames",-fit_limit_kappa[ibin2],fit_limit_kappa[ibin2]);
          
            if(ibin==0) c_chg_qg_MC_reco_kappa[ibin2]->cd(6);
            else c_chg_qg_MC_reco_kappa[ibin2]->cd(11-ibin);   

            h_chg_reco_fit_total_kappa[ibin][ibin2] = (TH1D*)h_chg_reco_ud_scaled_kappa[ibin][ibin2]->Clone("h_chg_reco_fit_total_kappa");
            h_chg_reco_fit_total_kappa[ibin][ibin2]->Scale(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));
            
            h_chg_reco_gubdb_afterfit_kappa[ibin][ibin2] = (TH1D*)h_chg_reco_gubdb_kappa[ibin][ibin2]->Clone("h_chg_reco_gubdb_afterfit_kappa");
            h_chg_reco_gubdb_afterfit_kappa[ibin][ibin2]->Scale(1. - f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));
            h_chg_reco_fit_total_kappa[ibin][ibin2]->Add(h_chg_reco_gubdb_afterfit_kappa[ibin][ibin2]); 

            h_chg_reco_ratio_2param_kappa[ibin][ibin2] = (TH1D*)h_chg_reco_kappa[ibin][ibin2]->Clone("h_chg_reco_ratio_kappa");
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->Divide(h_chg_reco_fit_total_kappa[ibin][ibin2]);
            //h_chg_reco_ratio_2param[ibin][ibin2]->Rebin(10); h_chg_reco_ratio_2param[ibin][ibin2]->Scale(1./10.);
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetXaxis()->SetTitle("reco jet charge");
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetYaxis()->SetTitle("Reco / Fit");
            h_cosm(h_chg_reco_ratio_2param_kappa[ibin][ibin2]);
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_kappa[ibin2],fit_limit_kappa[ibin2]);
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetYaxis()->SetRangeUser(0.8,1.2);
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio_2param_kappa[ibin][ibin2]->SetMarkerColor(kBlack);
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->SetMarkerStyle(1); h_chg_reco_ratio_2param_kappa[ibin][ibin2]->Draw("e0 same");
            drawline(-2.,1.,2.,1.);
            sprintf(g_frac, "g+others (MC/Fit) %.3f",(1.-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)))/gluon_ubar_dbar_reco[ibin]);
            sprintf(g_error, "#pm %.3f",f_chg_qg_kappa[ibin][ibin2]->GetParError(0));
            sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_qg_kappa[ibin][ibin2]->GetChisquare()/f_chg_qg_kappa[ibin][ibin2]->GetNDF());
            sprintf(up_frac, "up+down (MC/Fit) %.3f",(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)/(up_reco[ibin]+down_reco[ibin])));
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
/*
        TCanvas *c_chg_udg_MC_MC_rebin[ntrkbins];

        for(int ibin2=0;ibin2<ntrkbins;ibin2++){
            c_chg_udg_MC_MC_rebin[ibin2] = new TCanvas((TString)("c_chg_udg_MC_MC_rebin"+trk[ibin2]),(TString)("c_chg_udg_MC_MC_rebin"+trk[ibin2]),1500,450);

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
              h_chg_reco_cop[ibin][ibin2]->GetXaxis()->SetTitle("jet charge");
              h_cosm(h_chg_reco_cop[ibin][ibin2]);
              h_chg_reco_cop[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
              h_chg_reco_cop[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.7);
              h_chg_reco_cop[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.2);

              h_chg_reco_cop[ibin][ibin2]->Draw("e0 same");

              h_chg_reco_up_cop[ibin][ibin2]->Scale(0.65*f_chg_qg_data[ibin][ibin2]->GetParameter(0));
              h_chg_reco_down_cop[ibin][ibin2]->Scale(0.35*f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 

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
                legendx ->AddEntry(h_chg_reco_cop[0][ibin2], "MC (Gen)", "lep");
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

                h_chg_reco_ratio_2param[ibin][ibin2]->Draw("same");
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
                h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetLabelSize(0.12);
                h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
                h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetRangeUser(0.5,1.5);
                h_chg_reco_ratio_2param[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerColor(kBlack);
                h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerStyle(21); h_chg_reco_ratio_2param[ibin][ibin2]->Draw("e0 same");

                drawline(-fit_limit_ptcut[ibin2],1.,fit_limit_ptcut[ibin2],1.);
              }
              
            }
         }
*/
      }

    }

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
            h_chg_data[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
            h_chg_data[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);
            h_chg_data[ibin][ibin2]->SetLineColor(kGreen); h_chg_data[ibin][ibin2]->SetMarkerColor(kGreen); h_chg_data[ibin][ibin2]->SetMarkerStyle(2);
          
            h_chg_data[ibin][ibin2]->Draw("e0 same");
            h_chg_reco_ud_scaled[ibin][ibin2]->Draw("e0 same");
            h_chg_reco_gubdb[ibin][ibin2]->Draw("e0 same");

            //h_chg_data_trkinc[ibin][ibin2]->Draw("hist same");
            //h_chg_data_trkdec[ibin][ibin2]->Draw("hist same");
            
            TString tmp1 = cent_tag_pp_data[ibin];
            tx1 = new TLatex(); tx1->SetTextSize(.1);
            tx1->DrawLatexNDC(0.15,0.85, tmp1);
          
            TString tmp = trk_tag[ibin2];
            tx = new TLatex(); tx->SetTextSize(.1);
            tx->DrawLatexNDC(0.15,0.65, tmp);

            h_chg_data[ibin][ibin2]->Fit(f_chg_qg_data[ibin][ibin2],"Q N M R","sames",-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
        
            if(ibin==0) c_chg_qg_MC_data[ibin2]->cd(6);
            else c_chg_qg_MC_data[ibin2]->cd(11-ibin);   

            h_chg_reco_fit_total[ibin][ibin2] = (TH1D*)h_chg_reco_ud_scaled[ibin][ibin2]->Clone("h_chg_reco_fit_total");
            h_chg_reco_fit_total[ibin][ibin2]->Scale(f_chg_qg_data[ibin][ibin2]->GetParameter(0));
            
            h_chg_reco_gubdb_afterfit[ibin][ibin2] = (TH1D*)h_chg_reco_gubdb[ibin][ibin2]->Clone("h_chg_reco_gubdb_afterfit");
            h_chg_reco_gubdb_afterfit[ibin][ibin2]->Scale(1. - f_chg_qg_data[ibin][ibin2]->GetParameter(0));
            h_chg_reco_fit_total[ibin][ibin2]->Add(h_chg_reco_gubdb_afterfit[ibin][ibin2]); 

            h_chg_reco_ratio_2param[ibin][ibin2] = (TH1D*)h_chg_data[ibin][ibin2]->Clone("h_chg_data_ratio");
            h_chg_reco_ratio_2param[ibin][ibin2]->Divide(h_chg_reco_fit_total[ibin][ibin2]); 

            //h_chg_reco_ratio_2param[ibin][ibin2]->Rebin(10); h_chg_reco_ratio_2param[ibin][ibin2]->Scale(1./10.);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetTitle("data jet charge");
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetTitle("Data / Fit");
            h_cosm(h_chg_reco_ratio_2param[ibin][ibin2]);  
            h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
            h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetRangeUser(0.5,1.5);
            h_chg_reco_ratio_2param[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerColor(kRed);
            h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerStyle(21); h_chg_reco_ratio_2param[ibin][ibin2]->Draw("e0 same");
            //tl1->Draw("same");
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
 
        //final plot
        TCanvas *c_chg_udg_MC_data_rebin[ntrkbins];

        for(int ibin2=0;ibin2<ntrkbins;ibin2++){
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
              h_chg_data_cop_sysm[ibin][ibin2]->GetXaxis()->SetTitle("jet charge");
              h_cosm(h_chg_data_cop_sysm[ibin][ibin2]);
              h_chg_data_cop_sysm[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
              h_chg_data_cop_sysm[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.7);
              h_chg_data_cop_sysm[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.18);

              //h_chg_data_cop_sysm[ibin][ibin2]->Draw("e0 same");
              h_chg_data_cop_sysm[ibin][ibin2]->Draw("e2 same");

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

              //h_chg_data_trkinc[ibin][ibin2]->Draw("hist same");  
              //h_chg_data_trkdec[ibin][ibin2]->Draw("hist same");

              h_chg_data_cop_sysm[ibin][ibin2]->Draw("e2 same");  

              TString tmp1 = cent_tag_pp_data[ibin];
              tx1 = new TLatex(); tx1->SetTextSize(.1);
              tx1->DrawLatexNDC(0.15,0.85, tmp1);

              if(ibin==0){
                TLegend *legendx = new TLegend(0.15,0.7,0.5,0.8);
                legendx ->SetLineColor(kWhite);
                legendx ->AddEntry(h_chg_data_cop_sysm[0][ibin2], "Data", "lepf");
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

                h_chg_reco_ratio_2param[ibin][ibin2]->Draw("same");
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
                h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetLabelSize(0.12);
                h_chg_reco_ratio_2param[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
                h_chg_reco_ratio_2param[ibin][ibin2]->GetYaxis()->SetRangeUser(0.5,1.5);
                h_chg_reco_ratio_2param[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerColor(kBlack);
                h_chg_reco_ratio_2param[ibin][ibin2]->SetMarkerStyle(21); h_chg_reco_ratio_2param[ibin][ibin2]->Draw("e0 same");

                drawline(-fit_limit_ptcut[ibin2],1.,fit_limit_ptcut[ibin2],1.);
              }
              
            }
        }

        //systematics for q and g
        for(int ibin2=0;ibin2<ntrkbins;ibin2++){
           for(int ibin=0;ibin<nCBins;ibin++){    
              h_chg_data_trkinc[ibin][ibin2]->Fit(f_chg_qg_data[ibin][ibin2],"Q N M R","sames",-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);

              if(ibin==0){
                up_fit_high[ibin][ibin2] = 0.65*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
                down_fit_high[ibin][ibin2] = 0.35*(f_chg_qg_data[ibin][ibin2]->GetParameter(0));
              }
              else{
                up_fit_high[ibin][ibin2] = 0.47*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
                down_fit_high[ibin][ibin2] = 0.53*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
              }
              quark_fit_high[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParameter(0)+((1.-gluon_ubdbcsb[ibin])*(1.-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)))); 
              quark_fit_high[ibin][ibin2] = quark_fit_high[ibin][ibin2]*quark_corr_factor[ibin];
              quark_fit_high[ibin][ibin2] = (quark_fit_high[ibin][ibin2]*data_no0trk_int[ibin][ibin2] + quark_data_0trk_int[ibin][ibin2]);
              gluon_fit_high[ibin][ibin2] = 1.-(quark_fit_high[ibin][ibin2]); 

              h_chg_data_trkdec[ibin][ibin2]->Fit(f_chg_qg_data[ibin][ibin2],"Q N M R","sames",-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);

              if(ibin==0){
                up_fit_low[ibin][ibin2] = 0.65*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
                down_fit_low[ibin][ibin2] = 0.35*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
              }
              else{
                up_fit_low[ibin][ibin2] = 0.47*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
                down_fit_low[ibin][ibin2] = 0.53*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
              }
              quark_fit_low[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParameter(0)+((1.-gluon_ubdbcsb[ibin])*(1.-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)))); 
              quark_fit_low[ibin][ibin2] = quark_fit_low[ibin][ibin2]*quark_corr_factor[ibin];
              quark_fit_low[ibin][ibin2] = (quark_fit_low[ibin][ibin2]*data_no0trk_int[ibin][ibin2] + quark_data_0trk_int[ibin][ibin2]);
              gluon_fit_low[ibin][ibin2] = 1.-(quark_fit_low[ibin][ibin2]);

            //trkposneg
              h_chg_data_svd[ibin][ibin2]->Fit(f_chg_qg_data[ibin][ibin2],"Q N M R","sames",-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);

              if(ibin==0){
                up_fit_posneg[ibin][ibin2] = 0.65*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
                down_fit_posneg[ibin][ibin2] = 0.35*(f_chg_qg_data[ibin][ibin2]->GetParameter(0));
              }
              else{
                up_fit_posneg[ibin][ibin2] = 0.47*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
                down_fit_posneg[ibin][ibin2] = 0.53*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
              }
              quark_fit_posneg[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParameter(0)+((1.-gluon_ubdbcsb[ibin])*(1.-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)))); 
              quark_fit_posneg[ibin][ibin2] = quark_fit_posneg[ibin][ibin2]*quark_corr_factor[ibin];
              quark_fit_posneg[ibin][ibin2] = (quark_fit_posneg[ibin][ibin2]*data_no0trk_int[ibin][ibin2] + quark_data_0trk_int[ibin][ibin2]);
              gluon_fit_posneg[ibin][ibin2] = 1.-(quark_fit_posneg[ibin][ibin2]);


            //JER
              h_chg_data_jer[ibin][ibin2]->Fit(f_chg_qg_data[ibin][ibin2],"Q N M R","sames",-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);

              if(ibin==0){
                up_fit_jer[ibin][ibin2] = 0.65*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
                down_fit_jer[ibin][ibin2] = 0.35*(f_chg_qg_data[ibin][ibin2]->GetParameter(0));
              }
              else{
                up_fit_jer[ibin][ibin2] = 0.47*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
                down_fit_jer[ibin][ibin2] = 0.53*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
              }
              quark_fit_jer[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParameter(0)+((1.-gluon_ubdbcsb[ibin])*(1.-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)))); 
              quark_fit_jer[ibin][ibin2] = quark_fit_jer[ibin][ibin2]*quark_corr_factor[ibin];
              quark_fit_jer[ibin][ibin2] = (quark_fit_jer[ibin][ibin2]*data_no0trk_int[ibin][ibin2] + quark_data_0trk_int[ibin][ibin2]);
              gluon_fit_jer[ibin][ibin2] = 1.-(quark_fit_jer[ibin][ibin2]);

            //glu sample
              h_chg_data_glusamp[ibin][ibin2]->Fit(f_chg_qg_data[ibin][ibin2],"Q N M R","sames",-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);

              if(ibin==0){
                up_fit_glusamp[ibin][ibin2] = 0.65*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
                down_fit_glusamp[ibin][ibin2] = 0.35*(f_chg_qg_data[ibin][ibin2]->GetParameter(0));
              }
              else{
                up_fit_glusamp[ibin][ibin2] = 0.47*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
                down_fit_glusamp[ibin][ibin2] = 0.53*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); 
              }
              quark_fit_glusamp[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParameter(0)+((1.-gluon_ubdbcsb[ibin])*(1.-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)))); 
              quark_fit_glusamp[ibin][ibin2] = quark_fit_glusamp[ibin][ibin2]*quark_corr_factor[ibin];
              quark_fit_glusamp[ibin][ibin2] = (quark_fit_glusamp[ibin][ibin2]*data_no0trk_int[ibin][ibin2] + quark_data_0trk_int[ibin][ibin2]);
              gluon_fit_glusamp[ibin][ibin2] = 1.-(quark_fit_glusamp[ibin][ibin2]);

              h_chg_data[ibin][ibin2]->Fit(f_chg_qg_data[ibin][ibin2],"Q N M R","sames",-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
            }
         }
      }  

      if(do_kappa_fitting){
        //////2 parameter fitting data
        TCanvas *c_chg_qg_MC_data_kappa[ntrkbins];

        for(int ibin2=0;ibin2<ntrkbins;ibin2++){
          c_chg_qg_MC_data_kappa[ibin2] = new TCanvas((TString)("c_chg_qg_MC_data_kappa"+trk[ibin2]),(TString)("c_chg_qg_MC_data_kappa"+trk[ibin2]),1800,600);
          c_chg_qg_MC_data_kappa[ibin2]->Divide(5,2,0);
          //gStyle->SetOptStat(0);
          gStyle->SetOptFit(1);
          for(int ibin=0;ibin<nCBins;ibin++){
        
            if(ibin==0) c_chg_qg_MC_data_kappa[ibin2]->cd(1);
            else c_chg_qg_MC_data_kappa[ibin2]->cd(6-ibin);   
            h_chg_data_kappa[ibin][ibin2]->GetXaxis()->SetTitle("data jet chg");
            h_chg_data_kappa[ibin][ibin2]->GetYaxis()->SetTitle("# of jets");
            h_cosm(h_chg_data_kappa[ibin][ibin2]);
            h_chg_data_kappa[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_kappa[ibin2],fit_limit_kappa[ibin2]);
            h_chg_data_kappa[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.15);
            h_chg_data_kappa[ibin][ibin2]->SetLineColor(kGreen); h_chg_data_kappa[ibin][ibin2]->SetMarkerColor(kGreen); h_chg_data_kappa[ibin][ibin2]->SetMarkerStyle(2);
          
            h_chg_data_kappa[ibin][ibin2]->Draw("e0 same");
            h_chg_reco_ud_scaled_kappa[ibin][ibin2]->Draw("e0 same");
            h_chg_reco_gubdb_kappa[ibin][ibin2]->Draw("e0 same");

            //h_chg_data_kappa_trkinc[ibin][ibin2]->Draw("hist same");
            //h_chg_data_kappa_trkdec[ibin][ibin2]->Draw("hist same");
            
            TString tmp1 = cent_tag_pp_data[ibin];
            tx1 = new TLatex(); tx1->SetTextSize(.1);
            tx1->DrawLatexNDC(0.15,0.85, tmp1);
          
            TString tmp = trk_tag[ibin2];
            tx = new TLatex(); tx->SetTextSize(.1);
            tx->DrawLatexNDC(0.15,0.65, tmp);

            h_chg_data_kappa[ibin][ibin2]->Fit(f_chg_qg_kappa[ibin][ibin2],"Q N M R","sames",-fit_limit_kappa[ibin2],fit_limit_kappa[ibin2]);
        
            if(ibin==0) c_chg_qg_MC_data_kappa[ibin2]->cd(6);
            else c_chg_qg_MC_data_kappa[ibin2]->cd(11-ibin);   

            h_chg_reco_fit_total_kappa[ibin][ibin2] = (TH1D*)h_chg_reco_ud_scaled_kappa[ibin][ibin2]->Clone("h_chg_reco_fit_total");
            h_chg_reco_fit_total_kappa[ibin][ibin2]->Scale(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));
            
            h_chg_reco_gubdb_afterfit_kappa[ibin][ibin2] = (TH1D*)h_chg_reco_gubdb_kappa[ibin][ibin2]->Clone("h_chg_reco_gubdb_afterfit");
            h_chg_reco_gubdb_afterfit_kappa[ibin][ibin2]->Scale(1. - f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));
            h_chg_reco_fit_total_kappa[ibin][ibin2]->Add(h_chg_reco_gubdb_afterfit_kappa[ibin][ibin2]); 

            h_chg_reco_ratio_2param_kappa[ibin][ibin2] = (TH1D*)h_chg_data_kappa[ibin][ibin2]->Clone("h_chg_data_ratio");
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->Divide(h_chg_reco_fit_total_kappa[ibin][ibin2]); 

            //h_chg_reco_ratio_2param_kappa[ibin][ibin2]->Rebin(10); h_chg_reco_ratio_2param_kappa[ibin][ibin2]->Scale(1./10.);
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetXaxis()->SetTitle("data jet charge");
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetYaxis()->SetTitle("Data / Fit");
            h_cosm(h_chg_reco_ratio_2param_kappa[ibin][ibin2]);  
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_kappa[ibin2],fit_limit_kappa[ibin2]);
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetYaxis()->SetRangeUser(0.5,1.5);
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->SetLineColor(kRed); h_chg_reco_ratio_2param_kappa[ibin][ibin2]->SetMarkerColor(kRed);
            h_chg_reco_ratio_2param_kappa[ibin][ibin2]->SetMarkerStyle(21); h_chg_reco_ratio_2param_kappa[ibin][ibin2]->Draw("e0 same");
            //tl1->Draw("same");
            sprintf(g_frac, "g+others %.3f",(1.-f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)));
            sprintf(g_error, "#pm %.3f",f_chg_qg_kappa[ibin][ibin2]->GetParError(0));
            sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_qg_kappa[ibin][ibin2]->GetChisquare()/f_chg_qg_kappa[ibin][ibin2]->GetNDF());
            sprintf(up_frac, "up+down %.3f",f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));
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

        //final plot
        TCanvas *c_chg_udg_MC_data_rebin_kappa[ntrkbins];

        for(int ibin2=0;ibin2<ntrkbins;ibin2++){
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
            pad_label->cd();
            ylabel_tex->DrawLatexNDC(0.7,0.5, y_label);
            ylabel_tex->DrawLatexNDC(0.7,0.06, datafit);

            for(int ibin=0;ibin<nCBins;ibin++){

              if(ibin==0) pad1->cd(1);
              else pad1->cd(6-ibin);   
              h_chg_data_kappa_cop_sysm[ibin][ibin2]->GetXaxis()->SetTitle("jet charge");
              h_cosm(h_chg_data_kappa_cop_sysm[ibin][ibin2]);
              h_chg_data_kappa_cop_sysm[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_kappa[ibin2],fit_limit_kappa[ibin2]);
              h_chg_data_kappa_cop_sysm[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.7);
              h_chg_data_kappa_cop_sysm[ibin][ibin2]->GetYaxis()->SetRangeUser(0.,0.18);

              h_chg_data_kappa_cop_sysm[ibin][ibin2]->Draw("e2 same");

              //scaling according to data fractions
              if(ibin==0){
                  h_chg_reco_up_kappa_cop[ibin][ibin2]->Scale(0.65*f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));
                  h_chg_reco_down_kappa_cop[ibin][ibin2]->Scale(0.35*f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
              }
              else{
                  h_chg_reco_up_kappa_cop[ibin][ibin2]->Scale(0.47*f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));
                  h_chg_reco_down_kappa_cop[ibin][ibin2]->Scale(0.53*f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));   
              }
              h_chg_reco_g_kappa_cop[ibin][ibin2]->Scale((1.-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)))*gluon_ubdbcsb[ibin]);
              h_chg_reco_others_kappa_cop[ibin][ibin2]->Scale((1.-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)))*(1-gluon_ubdbcsb[ibin]));

              h_chg_data_stk_kappa[ibin][ibin2] = new THStack((TString)("h_chg_data_stk_kappa"+cent[ibin]+"_"+trk[ibin2]),"");
              h_chg_data_stk_kappa[ibin][ibin2]->Add(h_chg_reco_others_kappa_cop[ibin][ibin2]);
              h_chg_data_stk_kappa[ibin][ibin2]->Add(h_chg_reco_down_kappa_cop[ibin][ibin2]);
              h_chg_data_stk_kappa[ibin][ibin2]->Add(h_chg_reco_up_kappa_cop[ibin][ibin2]);
              h_chg_data_stk_kappa[ibin][ibin2]->Add(h_chg_reco_g_kappa_cop[ibin][ibin2]);

              h_chg_data_stk_kappa[ibin][ibin2]->Draw("hist same");      

              h_chg_data_kappa_cop_sysm[ibin][ibin2]->Draw("e2 same");  

              TString tmp1 = cent_tag_pp_data[ibin];
              tx1 = new TLatex(); tx1->SetTextSize(.1);
              tx1->DrawLatexNDC(0.15,0.85, tmp1);

              if(ibin==0){
                TLegend *legendx = new TLegend(0.15,0.7,0.5,0.8);
                legendx ->SetLineColor(kWhite);
                legendx ->AddEntry(h_chg_data_kappa_cop_sysm[0][ibin2], "Data", "lepf");
                legendx ->Draw("same");

                TLegend *legendy = new TLegend(0.62,0.62,0.98,0.98);
                legendy ->SetLineColor(kWhite);
                legendy ->AddEntry((TObject*)0, "Fitting results", "");
                legendy ->AddEntry(h_chg_reco_g_kappa_cop[0][ibin2], "Gluon", "lepf");
                legendy ->AddEntry(h_chg_reco_up_kappa_cop[0][ibin2], "Up quark", "lepf");
                legendy ->AddEntry(h_chg_reco_down_kappa_cop[0][ibin2], "Down quark", "lepf");
                legendy ->AddEntry(h_chg_reco_others_kappa_cop[0][ibin2], "Other flavors", "lepf");
                legendy ->Draw("same");
              }

              if(do_fit==true){  

                if(ibin==0) pad2->cd(1);
                else pad2->cd(6-ibin);   

                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->Draw("same");
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetXaxis()->SetTitle("");
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetYaxis()->SetTitle("");
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetYaxis()->SetTitleSize(0.05);
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetYaxis()->SetLabelSize(0.13);
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetXaxis()->CenterTitle();
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetYaxis()->CenterTitle();
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetYaxis()->SetTitleOffset(0.7);
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetYaxis()->SetNdivisions(505);
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetXaxis()->SetNdivisions(505);
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetXaxis()->SetTitleSize(0.05);
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetXaxis()->SetLabelSize(0.12);
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetXaxis()->SetRangeUser(-fit_limit_kappa[ibin2],fit_limit_kappa[ibin2]);
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->GetYaxis()->SetRangeUser(0.5,1.5);
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->SetLineColor(kBlack); h_chg_reco_ratio_2param_kappa[ibin][ibin2]->SetMarkerColor(kBlack);
                h_chg_reco_ratio_2param_kappa[ibin][ibin2]->SetMarkerStyle(21); h_chg_reco_ratio_2param_kappa[ibin][ibin2]->Draw("e0 same");

                drawline(-fit_limit_kappa[ibin2],1.,fit_limit_kappa[ibin2],1.);
              }
              
            }
        }


        //systematics   
        for(int ibin2=0;ibin2<ntrkbins;ibin2++){
           for(int ibin=0;ibin<nCBins;ibin++){    
              h_chg_data_kappa_trkinc[ibin][ibin2]->Fit(f_chg_qg_kappa[ibin][ibin2],"Q N M R","sames",-fit_limit_kappa[ibin2],fit_limit_kappa[ibin2]);

              if(ibin==0){
                up_fit_high_kappa[ibin][ibin2] = 0.65*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
                down_fit_high_kappa[ibin][ibin2] = 0.35*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));
              }
              else{
                up_fit_high_kappa[ibin][ibin2] = 0.47*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
                down_fit_high_kappa[ibin][ibin2] = 0.53*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
              }
              quark_fit_high_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)+((1.-gluon_ubdbcsb[ibin])*(1.-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)))); 
              quark_fit_high_kappa[ibin][ibin2] = quark_fit_high_kappa[ibin][ibin2]*quark_corr_factor[ibin];
              quark_fit_high_kappa[ibin][ibin2] = (quark_fit_high_kappa[ibin][ibin2]*data_no0trk_int[ibin][0] + quark_data_0trk_int[ibin][0]);
              gluon_fit_high_kappa[ibin][ibin2] = 1.-(quark_fit_high_kappa[ibin][ibin2]); 

              h_chg_data_kappa_trkdec[ibin][ibin2]->Fit(f_chg_qg_kappa[ibin][ibin2],"Q N M R","sames",-fit_limit_kappa[ibin2],fit_limit_kappa[ibin2]);

              if(ibin==0){
                up_fit_low_kappa[ibin][ibin2] = 0.65*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
                down_fit_low_kappa[ibin][ibin2] = 0.35*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
              }
              else{
                up_fit_low_kappa[ibin][ibin2] = 0.47*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
                down_fit_low_kappa[ibin][ibin2] = 0.53*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
              }
              quark_fit_low_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)+((1.-gluon_ubdbcsb[ibin])*(1.-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)))); 
              quark_fit_low_kappa[ibin][ibin2] = quark_fit_low_kappa[ibin][ibin2]*quark_corr_factor[ibin];
              quark_fit_low_kappa[ibin][ibin2] = (quark_fit_low_kappa[ibin][ibin2]*data_no0trk_int[ibin][0] + quark_data_0trk_int[ibin][0]);
              gluon_fit_low_kappa[ibin][ibin2] = 1.-(quark_fit_low_kappa[ibin][ibin2]);

            //trkposneg
              h_chg_data_kappa_trkposneg[ibin][ibin2]->Fit(f_chg_qg_kappa[ibin][ibin2],"Q N M R","sames",-fit_limit_kappa[ibin2],fit_limit_kappa[ibin2]);

              if(ibin==0){
                up_fit_posneg_kappa[ibin][ibin2] = 0.65*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
                down_fit_posneg_kappa[ibin][ibin2] = 0.35*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));
              }
              else{
                up_fit_posneg_kappa[ibin][ibin2] = 0.47*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
                down_fit_posneg_kappa[ibin][ibin2] = 0.53*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
              }
              quark_fit_posneg_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)+((1.-gluon_ubdbcsb[ibin])*(1.-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)))); 
              quark_fit_posneg_kappa[ibin][ibin2] = quark_fit_posneg_kappa[ibin][ibin2]*quark_corr_factor[ibin];
              quark_fit_posneg_kappa[ibin][ibin2] = (quark_fit_posneg_kappa[ibin][ibin2]*data_no0trk_int[ibin][0] + quark_data_0trk_int[ibin][0]);
              gluon_fit_posneg_kappa[ibin][ibin2] = 1.-(quark_fit_posneg_kappa[ibin][ibin2]);


            //JER
              h_chg_data_kappa_jer[ibin][ibin2]->Fit(f_chg_qg_kappa[ibin][ibin2],"Q N M R","sames",-fit_limit_kappa[ibin2],fit_limit_kappa[ibin2]);

              if(ibin==0){
                up_fit_kappa_jer[ibin][ibin2] = 0.65*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
                down_fit_kappa_jer[ibin][ibin2] = 0.35*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0));
              }
              else{
                up_fit_kappa_jer[ibin][ibin2] = 0.47*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
                down_fit_kappa_jer[ibin][ibin2] = 0.53*(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)); 
              }
              quark_fit_kappa_jer[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)+((1.-gluon_ubdbcsb[ibin])*(1.-(f_chg_qg_kappa[ibin][ibin2]->GetParameter(0)))); 
              quark_fit_kappa_jer[ibin][ibin2] = quark_fit_kappa_jer[ibin][ibin2]*quark_corr_factor[ibin];
              quark_fit_kappa_jer[ibin][ibin2] = (quark_fit_kappa_jer[ibin][ibin2]*data_no0trk_int[ibin][0] + quark_data_0trk_int[ibin][0]);
              gluon_fit_kappa_jer[ibin][ibin2] = 1.-(quark_fit_kappa_jer[ibin][ibin2]);

              h_chg_data_kappa[ibin][ibin2]->Fit(f_chg_qg_kappa[ibin][ibin2],"Q N M R","sames",-fit_limit_kappa[ibin2],fit_limit_kappa[ibin2]);
            }
         }
      }
    }

    Double_t up_fit[nCBins][ntrkbins], up_fiterr[nCBins][ntrkbins], down_fit[nCBins][ntrkbins], down_fiterr[nCBins][ntrkbins];
    Double_t gluon_fit[nCBins][ntrkbins], gluon_fiterr[nCBins][ntrkbins], quark_fit[nCBins][ntrkbins], quark_fiterr[nCBins][ntrkbins];
    Double_t gluon_fit_old[nCBins][ntrkbins], gluon_fiterr_old[nCBins][ntrkbins], gluon_fit_arc[nCBins][ntrkbins], gluon_fiterr_arc[nCBins][ntrkbins]; 
    Double_t gluon_fit_high_diff[nCBins][ntrkbins], gluon_fit_low_diff[nCBins][ntrkbins], gluon_fit_posneg_diff[nCBins][ntrkbins], gluon_fit_jer_diff[nCBins][ntrkbins]; 
    Double_t gluon_fit_glusamp_diff[nCBins][ntrkbins];
    Double_t gluon_fit_sysm[nCBins][ntrkbins]; 

    Double_t up_fit_kappa[nCBins][ntrkbins], up_fiterr_kappa[nCBins][ntrkbins], down_fit_kappa[nCBins][ntrkbins], down_fiterr_kappa[nCBins][ntrkbins];
    Double_t gluon_fit_kappa[nCBins][ntrkbins], gluon_fiterr_kappa[nCBins][ntrkbins];
    Double_t quark_fit_kappa[nCBins][ntrkbins], quark_fiterr_kappa[nCBins][ntrkbins];
    Double_t gluon_fit_old_kappa[nCBins][ntrkbins], gluon_fiterr_old_kappa[nCBins][ntrkbins]; 
    Double_t gluon_fit_arc_kappa[nCBins][ntrkbins], gluon_fiterr_arc_kappa[nCBins][ntrkbins]; 
    Double_t gluon_fit_high_diff_kappa[nCBins][ntrkbins], gluon_fit_low_diff_kappa[nCBins][ntrkbins];
    Double_t gluon_fit_posneg_diff_kappa[nCBins][ntrkbins], gluon_fit_jer_diff_kappa[nCBins][ntrkbins]; 
    Double_t gluon_fit_sysm_kappa[nCBins][ntrkbins];

    gluon_fit_old[0][0] = 0.59;
    gluon_fit_old[0][1] = 0.59;
    gluon_fit_old[0][2] = 0.59;
    gluon_fit_old[1][0] = 0.59;
    gluon_fit_old[1][1] = 0.59;
    gluon_fit_old[1][2] = 0.60;
    gluon_fit_old[2][0] = 0.60;
    gluon_fit_old[2][1] = 0.58;
    gluon_fit_old[2][2] = 0.60;
    gluon_fit_old[3][0] = 0.55;
    gluon_fit_old[3][1] = 0.55;
    gluon_fit_old[3][2] = 0.58;
    gluon_fit_old[4][0] = 0.55;
    gluon_fit_old[4][1] = 0.60;
    gluon_fit_old[4][2] = 0.61;

    gluon_fiterr_old[0][0] = 0.005;
    gluon_fiterr_old[0][1] = 0.005;
    gluon_fiterr_old[0][2] = 0.005;
    gluon_fiterr_old[1][0] = 0.01;
    gluon_fiterr_old[1][1] = 0.01;
    gluon_fiterr_old[1][2] = 0.01;
    gluon_fiterr_old[2][0] = 0.012;
    gluon_fiterr_old[2][1] = 0.012;
    gluon_fiterr_old[2][2] = 0.012;
    gluon_fiterr_old[3][0] = 0.025;
    gluon_fiterr_old[3][1] = 0.025;
    gluon_fiterr_old[3][2] = 0.025;
    gluon_fiterr_old[4][0] = 0.05;
    gluon_fiterr_old[4][1] = 0.05;
    gluon_fiterr_old[4][2] = 0.05;

    gluon_fit_arc[0][0] = 0.58;
    gluon_fit_arc[0][1] = 0.58;
    gluon_fit_arc[0][2] = 0.58;
    gluon_fit_arc[1][0] = 0.58;
    gluon_fit_arc[1][1] = 0.55;
    gluon_fit_arc[1][2] = 0.58;
    gluon_fit_arc[2][0] = 0.56;
    gluon_fit_arc[2][1] = 0.58;
    gluon_fit_arc[2][2] = 0.60;
    gluon_fit_arc[3][0] = 0.52;
    gluon_fit_arc[3][1] = 0.55;
    gluon_fit_arc[3][2] = 0.54;
    gluon_fit_arc[4][0] = 0.60;
    gluon_fit_arc[4][1] = 0.57;
    gluon_fit_arc[4][2] = 0.62;

    gluon_fiterr_arc[0][0] = 0.005;
    gluon_fiterr_arc[0][1] = 0.005;
    gluon_fiterr_arc[0][2] = 0.005;
    gluon_fiterr_arc[1][0] = 0.01;
    gluon_fiterr_arc[1][1] = 0.01;
    gluon_fiterr_arc[1][2] = 0.01;
    gluon_fiterr_arc[2][0] = 0.012;
    gluon_fiterr_arc[2][1] = 0.012;
    gluon_fiterr_arc[2][2] = 0.012;
    gluon_fiterr_arc[3][0] = 0.025;
    gluon_fiterr_arc[3][1] = 0.025;
    gluon_fiterr_arc[3][2] = 0.025;
    gluon_fiterr_arc[4][0] = 0.05;
    gluon_fiterr_arc[4][1] = 0.05;
    gluon_fiterr_arc[4][2] = 0.05;

    ////calculating final q and g values
    for(int ibin=0; ibin<nCBins; ibin++){
      for(int ibin2=0;ibin2<ntrkbins;ibin2++){
        
          if(do_ptcut_fitting){
            if(ibin==0){
              up_fit[ibin][ibin2] = 0.65*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); up_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
              down_fit[ibin][ibin2] = 0.35*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); down_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
            }
            else{
              up_fit[ibin][ibin2] = 0.47*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); up_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
              down_fit[ibin][ibin2] = 0.53*(f_chg_qg_data[ibin][ibin2]->GetParameter(0)); down_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
            }
            quark_fit[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParameter(0)+((1.-gluon_ubdbcsb[ibin])*(1.-(f_chg_qg_data[ibin][ibin2]->GetParameter(0)))); quark_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);
            quark_fit[ibin][ibin2] = quark_fit[ibin][ibin2]*quark_corr_factor[ibin];
            if(do_data) quark_fit[ibin][ibin2] = (quark_fit[ibin][ibin2]*data_no0trk_int[ibin][ibin2] + quark_data_0trk_int[ibin][ibin2]);
            else quark_fit[ibin][ibin2] = (quark_fit[ibin][ibin2]*mc_no0trk_int[ibin][ibin2] + quark_mc_0trk_int[ibin][ibin2]);
            gluon_fit[ibin][ibin2] = 1.-(quark_fit[ibin][ibin2]); gluon_fiterr[ibin][ibin2] = f_chg_qg_data[ibin][ibin2]->GetParError(0);

            gluon_fit_high_diff[ibin][ibin2] =  fabs(gluon_fit[ibin][ibin2] - gluon_fit_high[ibin][ibin2]);
            gluon_fit_high_diff[ibin][ibin2] +=  fabs(gluon_fit[ibin][ibin2] - gluon_fit_low[ibin][ibin2]);
            gluon_fit_high_diff[ibin][ibin2] = gluon_fit_high_diff[ibin][ibin2]/2.;

            gluon_fit_posneg_diff[ibin][ibin2] =  fabs(gluon_fit[ibin][ibin2] - gluon_fit_posneg[ibin][ibin2]);

            gluon_fit_jer_diff[ibin][ibin2] =  fabs(gluon_fit[ibin][ibin2] - gluon_fit_jer[ibin][ibin2]);

            gluon_fit_glusamp_diff[ibin][ibin2] =  fabs(gluon_fit[ibin][ibin2] - gluon_fit_glusamp[ibin][ibin2]);

            gluon_fit_sysm[ibin][ibin2] = sqrt(pow(gluon_fit_high_diff[ibin][ibin2],2)+pow(gluon_fit_posneg_diff[ibin][ibin2],2)+pow(gluon_fit_jer_diff[ibin][ibin2],2));  
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
            quark_fit_kappa[ibin][ibin2] = quark_fit_kappa[ibin][ibin2]*quark_corr_factor[ibin];
            quark_fit_kappa[ibin][ibin2] = (quark_fit_kappa[ibin][ibin2]*data_no0trk_int[ibin][0] + quark_data_0trk_int[ibin][0]);
            gluon_fit_kappa[ibin][ibin2] = 1.-(quark_fit_kappa[ibin][ibin2]); gluon_fiterr_kappa[ibin][ibin2] = f_chg_qg_kappa[ibin][ibin2]->GetParError(0);

            gluon_fit_high_diff_kappa[ibin][ibin2] =  fabs(gluon_fit_kappa[ibin][ibin2] - gluon_fit_high_kappa[ibin][ibin2]);
            gluon_fit_high_diff_kappa[ibin][ibin2] +=  fabs(gluon_fit_kappa[ibin][ibin2] - gluon_fit_low_kappa[ibin][ibin2]);
            gluon_fit_high_diff_kappa[ibin][ibin2] = gluon_fit_high_diff_kappa[ibin][ibin2]/2.;

            gluon_fit_posneg_diff_kappa[ibin][ibin2] =  fabs(gluon_fit_kappa[ibin][ibin2] - gluon_fit_posneg_kappa[ibin][ibin2]);

            gluon_fit_jer_diff_kappa[ibin][ibin2] =  fabs(gluon_fit_kappa[ibin][ibin2] - gluon_fit_kappa_jer[ibin][ibin2]);

            gluon_fit_sysm_kappa[ibin][ibin2] = sqrt(pow(gluon_fit_high_diff_kappa[ibin][ibin2],2)+pow(gluon_fit_posneg_diff_kappa[ibin][ibin2],2)+pow(gluon_fit_jer_diff_kappa[ibin][ibin2],2));  
            //gluon_fit_sysm_kappa[ibin][ibin2] = gluon_fit_sysm[ibin][ibin2];
          }                
      }

      //data qg graphs
      if(do_ptcut_fitting){
          up_data[ibin] = new TGraphErrors(ntrkbins,x_mean,up_fit[ibin],x_mean_err,up_fiterr[ibin]);
          up_data[ibin]->SetName((TString)("up_data_cent"+cent[ibin]));
          down_data[ibin] = new TGraphErrors(ntrkbins,x_mean,down_fit[ibin],x_mean_err,down_fiterr[ibin]);
          down_data[ibin]->SetName((TString)("down_data_cent"+cent[ibin]));
          gluon_data[ibin] = new TGraphErrors(ntrkbins,x_mean,gluon_fit[ibin],x_mean_err,gluon_fiterr[ibin]);
          gluon_data[ibin]->SetName((TString)("gluon_data_cent"+cent[ibin]));

          gluon_data_old[ibin] = new TGraphErrors(ntrkbins,x_mean,gluon_fit_old[ibin],x_mean_err,gluon_fiterr_old[ibin]);
          gluon_data_old[ibin]->SetName((TString)("gluon_data_old_cent"+cent[ibin]));
          gluon_data_arc[ibin] = new TGraphErrors(ntrkbins,x_mean,gluon_fit_arc[ibin],x_mean_err,gluon_fiterr_arc[ibin]);
          gluon_data_arc[ibin]->SetName((TString)("gluon_data_arc_cent"+cent[ibin]));

          //gluon_data_asymm_sys[ibin] = new TGraphAsymmErrors(ntrkbins,x_mean,gluon_fit[ibin],x_mean_err,x_mean_err,gluon_fit_sysm[ibin],gluon_fit_sysm[ibin]);
          quark_data[ibin] = new TGraphErrors(ntrkbins,x_mean,quark_fit[ibin],x_mean_err,quark_fiterr[ibin]);
          quark_data[ibin]->SetName((TString)("quark_data_cent"+cent[ibin]));
      }

      if(do_kappa_fitting){
          up_data_kappa[ibin] = new TGraphErrors(ntrkbins,kappa,up_fit_kappa[ibin],kappa_err,up_fiterr_kappa[ibin]);
          up_data_kappa[ibin]->SetName((TString)("up_data_kappa_cent"+cent[ibin]));
          down_data_kappa[ibin] = new TGraphErrors(ntrkbins,kappa,down_fit_kappa[ibin],kappa_err,down_fiterr_kappa[ibin]);
          down_data_kappa[ibin]->SetName((TString)("down_data_kappa_cent"+cent[ibin]));
          gluon_data_kappa[ibin] = new TGraphErrors(ntrkbins,kappa,gluon_fit_kappa[ibin],kappa_err,gluon_fiterr_kappa[ibin]);
          gluon_data_kappa[ibin]->SetName((TString)("gluon_data_kappa_cent"+cent[ibin]));
          //gluon_data_kappa_sysm[ibin] = new TGraphErrors(ntrkbins,kappa,gluon_fit_kappa[ibin],kappa_err,gluon_fit_sysm_kappa[ibin]);
          //gluon_data_kappa_sysm[ibin]->SetName((TString)("gluon_data_kappa_cent"+cent[ibin]));
          quark_data_kappa[ibin] = new TGraphErrors(ntrkbins,kappa,quark_fit_kappa[ibin],kappa_err,quark_fiterr_kappa[ibin]);
          quark_data_kappa[ibin]->SetName((TString)("quark_data_kappa_cent"+cent[ibin]));
      }

    }
  
    TCanvas *c_data_q_graph = new TCanvas("c_data_q_graph","c_data_q_graph",1500,350);
    drawgraphPad();
    gStyle->SetOptStat(0);

    pad_titleg->cd();
    title_tex->DrawLatexNDC(0.5,0.1, data);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    pad_labelg->cd();
    chg_tex->DrawLatexNDC(0.7,0.4, fraction);

    TH1D *h_data_ud_fit_results = new TH1D("h_data_ud_fit_results","",5,1.,6.);
    h_cosm(h_data_ud_fit_results);
    h_data_ud_fit_results->GetXaxis()->SetTitle("track p_{T} cut");
    h_data_ud_fit_results->GetXaxis()->SetRangeUser(1.,6.);
    h_data_ud_fit_results->GetYaxis()->SetTitle("");
    h_data_ud_fit_results->GetYaxis()->SetRangeUser(0.,1.);

    for(int ibin=0; ibin<nCBins; ibin++){
      if(ibin==0) pada->cd(1);
      else pada->cd(6-ibin);

      h_data_ud_fit_results->Draw("same");

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
       
      //if(ibin==0) {tl_u[ibin]->Draw("same"); tl_d[ibin]->Draw("same");} //tl_ubar[ibin]->Draw("same"); tl_dbar[ibin]->Draw("same"); tl_c[ibin]->Draw("same"); tl_s[ibin]->Draw("same");tl_b[ibin]->Draw("same");}
      //else {tl_u_data->Draw("same"); tl_d_data->Draw("same");} //tl_ubar[ibin]->Draw("same"); tl_dbar[ibin]->Draw("same"); tl_c[ibin]->Draw("same"); tl_s[ibin]->Draw("same");tl_b[ibin]->Draw("same");}

      if(ibin==0){
        TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
        legend ->SetLineColor(kWhite);
        legend ->AddEntry(up_data[ibin], "up", "lepf");
        legend ->AddEntry(down_data[ibin], "down", "lepf");
        legend ->Draw("same");
      }

      TString tmp1;
      if(do_data) tmp1 = cent_tag_pp_data[ibin];
      else tmp1 = cent_tag_pp_MC[ibin];
      tx1 = new TLatex(); tx1->SetTextSize(.1);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);
      
    }

    ////trkeff studies
    //////vary trk eff fitting
    
    TH1D* h_chi2ndf[nCBins][ntrkbins];
    TH1D* h_least_chi2ndf[ntrkbins];

    for(int ibin2=0; ibin2<ntrkbins; ibin2++){
        h_chi2ndf[0][ibin2] = new TH1D(Form("h_chi2ndf_%d_%d",0,ibin2),"",11,-5.,5.);
        h_chi2ndf[1][ibin2] = new TH1D(Form("h_chi2ndf_%d_%d",1,ibin2),"",11,-5.,5.);
        h_chi2ndf[2][ibin2] = new TH1D(Form("h_chi2ndf_%d_%d",2,ibin2),"",11,-5.,5.);
        h_chi2ndf[3][ibin2] = new TH1D(Form("h_chi2ndf_%d_%d",3,ibin2),"",11,-5.,5.);
        h_chi2ndf[4][ibin2] = new TH1D(Form("h_chi2ndf_%d_%d",4,ibin2),"",11,-5.,5.);

        h_least_chi2ndf[ibin2] = new TH1D(Form("h_least_chi2ndf_%d",ibin2),"",5,0.,5.);
    } 

    Double_t chi2NDF[nCBins][ntrkbins][trkeffbins+1];

    //////two parameter fitting

    if(do_data && do_eff_studies){
      TCanvas *c_chg_qg_MC_data_trkeff[ntrkbins];

      for(int ibin2=0;ibin2<ntrkbins;ibin2++){
        c_chg_qg_MC_data_trkeff[ibin2] = new TCanvas((TString)("c_chg_qg_MC_data_trkeff"+trk[ibin2]),(TString)("c_chg_qg_MC_data_trkeff"+trk[ibin2]),1500,300);
        c_chg_qg_MC_data_trkeff[ibin2]->Divide(5,1,0);
        for(int ibin=0;ibin<nCBins;ibin++){
      
          if(ibin==0) c_chg_qg_MC_data_trkeff[ibin2]->cd(1);
          else c_chg_qg_MC_data_trkeff[ibin2]->cd(6-ibin);   

          if(do_fit==true){  
            for(int ibin3=0; ibin3<trkeffbins;ibin3++){
               h_chg_data_trkeff[ibin][ibin2][ibin3]->Fit(f_chg_qg_data[ibin][ibin2],"Q N M R","sames",-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);

               h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3] = (TH1D*)h_chg_data_trkeff[ibin][ibin2][ibin3]->Clone("h_chg_data_ratio_trkeff");
               h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->Divide(f_chg_qg_data[ibin][ibin2]);
            
               chi2NDF[ibin][ibin2][ibin3] = f_chg_qg_data[ibin][ibin2]->GetChisquare()/f_chg_qg_data[ibin][ibin2]->GetNDF();
               h_chi2ndf[ibin][ibin2]->SetBinContent(ibin3+1,chi2NDF[ibin][ibin2][ibin3]);
            }

            h_chg_data[ibin][ibin2]->Fit(f_chg_qg_data[ibin][ibin2],"Q N M R","sames",-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);

            if(ibin==0) c_chg_qg_MC_data_trkeff[ibin2]->cd(1);
            else c_chg_qg_MC_data_trkeff[ibin2]->cd(6-ibin);   

            for(int ibin3=0; ibin3<trkeffbins;ibin3++){
                if(ibin3<5){
                  h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->SetLineColor(kCyan-5+ibin3); h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->SetMarkerColor(kCyan-5+ibin3);
                }
                else{
                  h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->SetLineColor(kGreen-5+ibin3); h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->SetMarkerColor(kGreen-5+ibin3);
                }
                h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->GetYaxis()->SetRangeUser(0.8,1.2);
                h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->GetXaxis()->SetRangeUser(-fit_limit_ptcut[ibin2],fit_limit_ptcut[ibin2]);
                h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->SetMarkerStyle(20); h_chg_reco_ratio_trkeff[ibin][ibin2][ibin3]->Draw("HIST L same");        
            }
            h_chg_reco_ratio_2param[ibin][ibin2]->Draw("HIST L same");
           
            if(ibin==0){
                TLegend *legendz = new TLegend(0.65,0.55,0.98,0.98);
                legendz ->SetLineColor(kWhite);
                legendz ->AddEntry(h_chg_reco_ratio_2param[0][0], "Data (nominal)", "l");
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

            TString tmp1 = cent_tag_pp_data[ibin];
            tx1 = new TLatex(); tx1->SetTextSize(.1);
            tx1->DrawLatexNDC(0.35,0.85, tmp1);

            drawline(-2.,1.,2.,1.);
          }
        }
      }

      ///chi2/ndf plot for varying trk eff
      TCanvas *c_chg_qg_MC_data_chi2ndf;
      c_chg_qg_MC_data_chi2ndf = new TCanvas("c_chg_qg_MC_data_chi2ndf","c_chg_qg_MC_data_chi2ndf",1800,300);
      c_chg_qg_MC_data_chi2ndf->Divide(5,1);

      for(int ibin=0;ibin<nCBins;ibin++){        
        for(int ibin2=0;ibin2<ntrkbins;ibin2++){
          //gStyle->SetOptStat(0);
          //gStyle->SetOptFit(1);
          if(ibin==0){
            h_least_chi2ndf[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp_data[ibin]);
            h_least_chi2ndf[ibin2]->SetBinContent(ibin+1,h_chi2ndf[ibin][ibin2]->GetBinContent(5));
            h_least_chi2ndf[ibin2]->SetBinError(ibin+1,0.);
          }
          else{
            h_least_chi2ndf[ibin2]->GetXaxis()->SetBinLabel(6-ibin,cent_tag_pp_data[ibin]);
            h_least_chi2ndf[ibin2]->SetBinContent(6-ibin,h_chi2ndf[ibin][ibin2]->GetBinContent(5));
            h_least_chi2ndf[ibin2]->SetBinError(6-ibin,0.);            
          }

          h_chi2ndf[ibin][ibin2]->GetXaxis()->SetTitle("trk. eff. vary (%)");
          h_chi2ndf[ibin][ibin2]->GetYaxis()->SetTitle("#chi^{2}/NDF");
          h_cosm(h_chi2ndf[ibin][ibin2]);
          h_chi2ndf[ibin][ibin2]->SetLineColor(ibin2+1); h_chi2ndf[ibin][ibin2]->SetMarkerColor(ibin2+1);
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
          TString tmp1 = cent_tag_pp_data[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.75, tmp1);
      }

      TCanvas *c_chg_qg_least_chi2ndf;
      c_chg_qg_least_chi2ndf = new TCanvas("c_chg_qg_least_chi2ndf","c_chg_qg_least_chi2ndf",500,500);
      c_chg_qg_least_chi2ndf->cd(1);

      TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
      legendx ->SetLineColor(kWhite);

      for(int ibin2=0;ibin2<ntrkbins;ibin2++){
        h_least_chi2ndf[ibin2]->GetXaxis()->SetTitle("");
        h_least_chi2ndf[ibin2]->GetYaxis()->SetTitle("#chi^{2}/NDF");
        h_cosm(h_least_chi2ndf[ibin2]);
        h_least_chi2ndf[ibin2]->GetYaxis()->SetRangeUser(0.,2.);
        h_least_chi2ndf[ibin2]->SetLineColor(ibin2+1); h_least_chi2ndf[ibin2]->SetMarkerColor(ibin2+1);
        h_least_chi2ndf[ibin2]->SetMarkerStyle(34); h_least_chi2ndf[ibin2]->SetMarkerSize(1.5);
        h_least_chi2ndf[ibin2]->Draw("p same");
        legendx ->AddEntry(h_least_chi2ndf[ibin2], trk_tag[ibin2], "lepf");
      }
      legendx ->Draw("same");

    }

    ///jet charge avg and rms

      TH1D *h_chg_data_abs[nCBins][ntrkbins];
      TH1D *h_chg_data_abs_trkinc[nCBins][ntrkbins];
      TH1D *h_chg_data_abs_trkdec[nCBins][ntrkbins];

      TH1D *h_chg_reco_abs[nCBins][ntrkbins];
      //TH1D *h_chg_reco_up_abs[nCBins][ntrkbins];
      //TH1D *h_chg_reco_down_abs[nCBins][ntrkbins];
      //TH1D *h_chg_reco_gubdb_abs[nCBins][ntrkbins];

      Double_t trkeff_rms_err[nCBins][ntrkbins], jer_rms_err[nCBins][ntrkbins], svd_rms_err[nCBins][ntrkbins], rms_err[nCBins][ntrkbins];
      Double_t trkeff_avg_err[nCBins][ntrkbins], jer_avg_err[nCBins][ntrkbins], svd_avg_err[nCBins][ntrkbins], avg_err[nCBins][ntrkbins];

      for(int ibin2=0; ibin2<ntrkbins; ibin2++){
          h_jetchg_avg_data[ibin2] = new TH1D(Form("h_jetchg_avg_data_%d",ibin2),"",4,0.,4.);
          h_jetchg_avg_data_trkinc[ibin2] = new TH1D(Form("h_jetchg_avg_data_trkinc_%d",ibin2),"",4,0.,4.);
          h_jetchg_avg_data_trkdec[ibin2] = new TH1D(Form("h_jetchg_avg_data_trkdec_%d",ibin2),"",4,0.,4.);
          h_jetchg_avg_MC[ibin2] = new TH1D(Form("h_jetchg_avg_MC_%d",ibin2),"",4,0.,4.);
          h_jetchg_avg_data_pyquen[ibin2] = new TH1D(Form("h_jetchg_avg_data_pyquen_%d",ibin2),"",4,0.,4.);
          h_jetchg_avg_data_sysm[ibin2] = new TH1D(Form("h_jetchg_avg_data_sysm_%d",ibin2),"",4,0.,4.);

          h_jetchg_rms_data[ibin2] = new TH1D(Form("h_jetchg_rms_data_%d",ibin2),"",4,0.,4.);
          h_jetchg_rms_data_sysm[ibin2] = new TH1D(Form("h_jetchg_rms_data_sysm_%d",ibin2),"",4,0.,4.);
          h_jetchg_rms_data_trkinc[ibin2] = new TH1D(Form("h_jetchg_rms_data_trkinc_%d",ibin2),"",4,0.,4.);
          h_jetchg_rms_data_trkdec[ibin2] = new TH1D(Form("h_jetchg_rms_data_trkdec_%d",ibin2),"",4,0.,4.);
          h_jetchg_rms_data_svd[ibin2] = new TH1D(Form("h_jetchg_rms_data_svd_%d",ibin2),"",4,0.,4.);
          h_jetchg_rms_data_jer[ibin2] = new TH1D(Form("h_jetchg_rms_data_jer_%d",ibin2),"",4,0.,4.);
          h_jetchg_rms_data_pyquen[ibin2] = new TH1D(Form("h_jetchg_rms_data_pyquen_%d",ibin2),"",4,0.,4.);
          h_jetchg_rms_data_pyquencoll[ibin2] = new TH1D(Form("h_jetchg_rms_data_pyquencoll_%d",ibin2),"",4,0.,4.);          
          h_jetchg_rms_MC[ibin2] = new TH1D(Form("h_jetchg_rms_MC_%d",ibin2),"",4,0.,4.);
      }

        //h_chg_reco_up[0][0]->Scale(up_reco[0]);
        //h_chg_reco_down[0][0]->Scale(down_reco[0]);
      h_chg_reco_gubdb[4][0]->Scale(1./h_chg_reco_gubdb[4][0]->Integral());

      for(int ibin=0;ibin<nCBins;ibin++){        
        for(int ibin2=0;ibin2<ntrkbins;ibin2++){
          h_chg_data_abs[ibin][ibin2] = (TH1D*)h_chg_data[ibin][ibin2]->Clone(Form("h_chg_data_abs_%d_%d",ibin,ibin2));
          //abs_histos(h_chg_data[ibin][ibin2], h_chg_data_abs[ibin][ibin2]); 
          diff_histos(h_chg_data[ibin][ibin2], h_chg_data_abs[ibin][ibin2]); 
          h_chg_data_abs_trkinc[ibin][ibin2] = (TH1D*)h_chg_data_trkinc[ibin][ibin2]->Clone(Form("h_chg_data_abs_trkinc_%d_%d",ibin,ibin2));
          abs_histos(h_chg_data_trkinc[ibin][ibin2], h_chg_data_abs_trkinc[ibin][ibin2]); 
          h_chg_data_abs_trkdec[ibin][ibin2] = (TH1D*)h_chg_data_trkdec[ibin][ibin2]->Clone(Form("h_chg_data_abs_trkdec_%d_%d",ibin,ibin2));
          abs_histos(h_chg_data_trkdec[ibin][ibin2], h_chg_data_abs_trkdec[ibin][ibin2]); 

          h_chg_reco_abs[ibin][ibin2] = (TH1D*)h_chg_reco[ibin][ibin2]->Clone(Form("h_chg_reco_abs_%d_%d",ibin,ibin2));
          //abs_histos(h_chg_reco[ibin][ibin2], h_chg_reco_abs[ibin][ibin2]); 
          diff_histos(h_chg_reco[ibin][ibin2], h_chg_reco_abs[ibin][ibin2]); 
          //h_chg_reco_up[ibin][ibin2]->Rebin(10); h_chg_reco_up[ibin][ibin2]->Scale(1./10.);
          h_chg_reco_up_abs[ibin][ibin2] = (TH1D*)h_chg_reco_up[ibin][ibin2]->Clone(Form("h_chg_reco_up_abs_%d_%d",ibin,ibin2));
          diff_histos(h_chg_reco_up[ibin][ibin2], h_chg_reco_up_abs[ibin][ibin2]); 
          //h_chg_reco_down[ibin][ibin2]->Rebin(10); h_chg_reco_down[ibin][ibin2]->Scale(1./10.);
          h_chg_reco_down_abs[ibin][ibin2] = (TH1D*)h_chg_reco_down[ibin][ibin2]->Clone(Form("h_chg_reco_down_abs_%d_%d",ibin,ibin2));
          diff_histos(h_chg_reco_down[ibin][ibin2], h_chg_reco_down_abs[ibin][ibin2]); 
          //inv_histos(h_chg_reco_down_abs[ibin][ibin2], h_chg_reco_down_abs[ibin][ibin2]);
          //h_chg_reco_gubdb[ibin][ibin2]->Scale(1./10.);
          h_chg_reco_gubdb_abs[ibin][ibin2] = (TH1D*)h_chg_reco_gubdb[ibin][ibin2]->Clone(Form("h_chg_reco_gubdb_abs_%d_%d",ibin,ibin2));
          //abs_histos(h_chg_reco_gubdb[ibin][ibin2], h_chg_reco_gubdb_abs[ibin][ibin2]); 
          diff_histos(h_chg_reco_gubdb[ibin][ibin2], h_chg_reco_gubdb_abs[ibin][ibin2]);           
          if(ibin==0){
            h_jetchg_avg_data[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp_data[ibin]);
            h_jetchg_avg_data[ibin2]->SetBinContent(ibin+1,0.987*h_chg_data_abs[ibin][ibin2]->GetMean());
            h_jetchg_avg_data[ibin2]->SetBinError(ibin+1,h_chg_data_abs[ibin][ibin2]->GetMeanError());

            h_jetchg_avg_data_trkinc[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp_data[ibin]);
            h_jetchg_avg_data_trkinc[ibin2]->SetBinContent(ibin+1,0.987*h_chg_data_abs_trkinc[ibin][ibin2]->GetMean());
            h_jetchg_avg_data_trkinc[ibin2]->SetBinError(ibin+1,h_chg_data_abs_trkinc[ibin][ibin2]->GetMeanError());

            h_jetchg_avg_data_trkdec[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp_data[ibin]);
            h_jetchg_avg_data_trkdec[ibin2]->SetBinContent(ibin+1,0.987*h_chg_data_abs_trkdec[ibin][ibin2]->GetMean());
            h_jetchg_avg_data_trkdec[ibin2]->SetBinError(ibin+1,h_chg_data_abs_trkdec[ibin][ibin2]->GetMeanError());
/*
            h_jetchg_rms_data[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp_data[ibin]);
            h_jetchg_rms_data[ibin2]->SetBinContent(ibin+1,h_chg_data[ibin][ibin2]->GetRMS());
            h_jetchg_rms_data[ibin2]->SetBinError(ibin+1,h_chg_data[ibin][ibin2]->GetRMSError());

            h_jetchg_rms_data_trkinc[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp_data[ibin]);
            h_jetchg_rms_data_trkinc[ibin2]->SetBinContent(ibin+1,h_chg_data_trkinc[ibin][ibin2]->GetRMS());
            h_jetchg_rms_data_trkinc[ibin2]->SetBinError(ibin+1,h_chg_data_trkinc[ibin][ibin2]->GetRMSError());

            h_jetchg_rms_data_trkdec[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp_data[ibin]);
            h_jetchg_rms_data_trkdec[ibin2]->SetBinContent(ibin+1,h_chg_data_trkdec[ibin][ibin2]->GetRMS());
            h_jetchg_rms_data_trkdec[ibin2]->SetBinError(ibin+1,h_chg_data_trkdec[ibin][ibin2]->GetRMSError());

            h_jetchg_rms_data_svd[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp_data[ibin]);
            h_jetchg_rms_data_svd[ibin2]->SetBinContent(ibin+1,h_chg_data_svd[ibin][ibin2]->GetRMS());
            h_jetchg_rms_data_svd[ibin2]->SetBinError(ibin+1,h_chg_data_svd[ibin][ibin2]->GetRMSError());

            h_jetchg_rms_data_jer[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp_data[ibin]);
            h_jetchg_rms_data_jer[ibin2]->SetBinContent(ibin+1,h_chg_data_jer[ibin][ibin2]->GetRMS());
            h_jetchg_rms_data_jer[ibin2]->SetBinError(ibin+1,h_chg_data_jer[ibin][ibin2]->GetRMSError());

            h_jetchg_rms_MC[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp_MC[ibin]);
            h_jetchg_rms_MC[ibin2]->SetBinContent(ibin+1,h_chg_reco[ibin][ibin2]->GetRMS());
            h_jetchg_rms_MC[ibin2]->SetBinError(ibin+1,h_chg_reco[ibin][ibin2]->GetRMSError());
*/
            h_jetchg_avg_MC[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp_data[ibin]);
/*
            cout<<"h_chg_reco[ibin][ibin2]->GetMean(): "<<h_chg_reco[ibin][ibin2]->GetMean()<<endl;
            cout<<"h_chg_reco_up[ibin][ibin2]->GetMean(): "<<h_chg_reco_up[ibin][ibin2]->GetMean()<<endl;
            cout<<"h_chg_reco_down[ibin][ibin2]->GetMean(): "<<h_chg_reco_down[ibin][ibin2]->GetMean()<<endl;
            cout<<"h_chg_reco_gubdb[ibin][ibin2]->GetMean(): "<<h_chg_reco_gubdb[ibin][ibin2]->GetMean()<<endl;
            cout<<"h_chg_reco[ibin][ibin2]->GetRMS(): "<<h_chg_reco[ibin][ibin2]->GetRMS()<<endl;
            cout<<"h_chg_reco_up[ibin][ibin2]->GetRMS(): "<<h_chg_reco_up[ibin][ibin2]->GetRMS()<<endl;
            cout<<"h_chg_reco_down[ibin][ibin2]->GetRMS(): "<<h_chg_reco_down[ibin][ibin2]->GetRMS()<<endl;
            cout<<"h_chg_reco_gubdb[ibin][ibin2]->GetRMS(): "<<h_chg_reco_gubdb[ibin][ibin2]->GetRMS()<<endl;
cout<<"..............................................."<<endl;
            cout<<"h_chg_reco_abs[ibin][ibin2]->GetMean(): "<<h_chg_reco_abs[ibin][ibin2]->GetMean()<<endl;
            cout<<"h_chg_reco_up_abs[ibin][ibin2]->GetMean(): "<<h_chg_reco_up_abs[ibin][ibin2]->GetMean()<<endl;
            cout<<"h_chg_reco_down_abs[ibin][ibin2]->GetMean(): "<<h_chg_reco_down_abs[ibin][ibin2]->GetMean()<<endl;
            cout<<"h_chg_reco_gubdb_abs[ibin][ibin2]->GetMean(): "<<h_chg_reco_gubdb_abs[ibin][ibin2]->GetMean()<<endl;
cout<<"..............................................."<<endl;
cout<<"..............................................."<<endl;
*/
            h_jetchg_avg_MC[ibin2]->SetBinContent(ibin+1,h_chg_reco[ibin][ibin2]->GetMean());
            h_jetchg_avg_MC[ibin2]->SetBinError(ibin+1,h_chg_reco[ibin][ibin2]->GetMeanError());
          }
          else{
            h_jetchg_avg_data[ibin2]->GetXaxis()->SetBinLabel(5-ibin,cent_tag_pp_data[ibin]);
            h_jetchg_avg_data[ibin2]->SetBinContent(5-ibin,h_chg_data[ibin][ibin2]->GetMean());
            h_jetchg_avg_data[ibin2]->SetBinError(5-ibin,h_chg_data[ibin][ibin2]->GetMeanError());

            h_jetchg_avg_data_trkinc[ibin2]->GetXaxis()->SetBinLabel(5-ibin,cent_tag_pp_data[ibin]);
            h_jetchg_avg_data_trkinc[ibin2]->SetBinContent(5-ibin,h_chg_data_abs_trkinc[ibin][ibin2]->GetMean());
            h_jetchg_avg_data_trkinc[ibin2]->SetBinError(5-ibin,h_chg_data_abs_trkinc[ibin][ibin2]->GetMeanError());

            h_jetchg_avg_data_trkdec[ibin2]->GetXaxis()->SetBinLabel(5-ibin,cent_tag_pp_data[ibin]);
            h_jetchg_avg_data_trkdec[ibin2]->SetBinContent(5-ibin,h_chg_data_abs_trkdec[ibin][ibin2]->GetMean());
            h_jetchg_avg_data_trkdec[ibin2]->SetBinError(5-ibin,h_chg_data_abs_trkdec[ibin][ibin2]->GetMeanError());

            h_jetchg_rms_data[ibin2]->GetXaxis()->SetBinLabel(5-ibin,cent_tag_pp_data[ibin]);
            h_jetchg_rms_data[ibin2]->SetBinContent(5-ibin,h_chg_data[ibin][ibin2]->GetRMS());
            h_jetchg_rms_data[ibin2]->SetBinError(5-ibin,h_chg_data[ibin][ibin2]->GetRMSError());

            h_jetchg_rms_data_trkinc[ibin2]->GetXaxis()->SetBinLabel(5-ibin,cent_tag_pp_data[ibin]);
            h_jetchg_rms_data_trkinc[ibin2]->SetBinContent(5-ibin,h_chg_data_trkinc[ibin][ibin2]->GetRMS());
            h_jetchg_rms_data_trkinc[ibin2]->SetBinError(5-ibin,h_chg_data_trkinc[ibin][ibin2]->GetRMSError());

            h_jetchg_rms_data_trkdec[ibin2]->GetXaxis()->SetBinLabel(5-ibin,cent_tag_pp_data[ibin]);
            h_jetchg_rms_data_trkdec[ibin2]->SetBinContent(5-ibin,h_chg_data_trkdec[ibin][ibin2]->GetRMS());
            h_jetchg_rms_data_trkdec[ibin2]->SetBinError(5-ibin,h_chg_data_trkdec[ibin][ibin2]->GetRMSError());

            h_jetchg_rms_data_svd[ibin2]->GetXaxis()->SetBinLabel(5-ibin,cent_tag_pp_data[ibin]);
            h_jetchg_rms_data_svd[ibin2]->SetBinContent(5-ibin,h_chg_data_svd[ibin][ibin2]->GetRMS());
            h_jetchg_rms_data_svd[ibin2]->SetBinError(5-ibin,h_chg_data_svd[ibin][ibin2]->GetRMSError());

            h_jetchg_rms_data_jer[ibin2]->GetXaxis()->SetBinLabel(5-ibin,cent_tag_pp_data[ibin]);
            h_jetchg_rms_data_jer[ibin2]->SetBinContent(5-ibin,h_chg_data_jer[ibin][ibin2]->GetRMS());
            h_jetchg_rms_data_jer[ibin2]->SetBinError(5-ibin,h_chg_data_jer[ibin][ibin2]->GetRMSError());

            h_jetchg_rms_MC[ibin2]->GetXaxis()->SetBinLabel(5-ibin,cent_tag_pp_MC[ibin]);
            h_jetchg_rms_MC[ibin2]->SetBinContent(5-ibin,h_chg_reco[ibin][ibin2]->GetRMS());
            h_jetchg_rms_MC[ibin2]->SetBinError(5-ibin,h_chg_reco[ibin][ibin2]->GetRMSError());

            h_jetchg_avg_MC[ibin2]->GetXaxis()->SetBinLabel(5-ibin,cent_tag_pp_data[ibin]);
            h_jetchg_avg_MC[ibin2]->SetBinContent(5-ibin,h_chg_reco[ibin][ibin2]->GetMean());
            h_jetchg_avg_MC[ibin2]->SetBinError(5-ibin,h_chg_reco[ibin][ibin2]->GetMeanError());

            h_jetchg_avg_data_pyquen[ibin2]->SetBinContent(ibin,0.);
            h_jetchg_avg_data_pyquen[ibin2]->SetBinError(ibin,0.004);  

            h_jetchg_rms_data_pyquen[ibin2]->SetBinContent(ibin,0.444);
            h_jetchg_rms_data_pyquen[ibin2]->SetBinError(ibin,0.0025);  

            h_jetchg_rms_data_pyquencoll[ibin2]->SetBinContent(ibin,0.426);
            h_jetchg_rms_data_pyquencoll[ibin2]->SetBinError(ibin,0.0025);  
          }

            trkeff_rms_err[ibin][ibin2] = fabs(h_jetchg_rms_data_trkinc[ibin2]->GetBinContent(ibin) - h_jetchg_rms_data_trkdec[ibin2]->GetBinContent(ibin))/2.;
            jer_rms_err[ibin][ibin2] = fabs(h_jetchg_rms_data_jer[ibin2]->GetBinContent(ibin) - h_jetchg_rms_data[ibin2]->GetBinContent(ibin));
            svd_rms_err[ibin][ibin2] = fabs(h_jetchg_rms_data_svd[ibin2]->GetBinContent(ibin) - h_jetchg_rms_data[ibin2]->GetBinContent(ibin));      
            if(ibin>2) rms_err[ibin][ibin2] = sqrt(pow(trkeff_rms_err[ibin][ibin2],2)+pow(jer_rms_err[ibin][ibin2],2)+pow(svd_rms_err[ibin][ibin2],2));
            //if(ibin==3) rms_err[ibin][ibin2] = 0.008;
            else rms_err[ibin][ibin2] = 0.005; 

            avg_err[ibin][ibin2] = fabs(h_jetchg_avg_data_trkinc[ibin2]->GetBinContent(ibin) - h_jetchg_avg_data_trkdec[ibin2]->GetBinContent(ibin))/2.;

            if(ibin==0){
                //h_jetchg_rms_data_sysm[ibin2]->GetXaxis()->SetBinLabel(ibin+1,cent_tag_pp_data[ibin]);
                //h_jetchg_rms_data_sysm[ibin2]->SetBinContent(1+ibin,h_chg_data[ibin][ibin2]->GetRMS());
                //h_jetchg_rms_data_sysm[ibin2]->SetBinError(1+ibin,rms_err[ibin][ibin2]);
            }
            else{
                h_jetchg_avg_data_sysm[ibin2]->GetXaxis()->SetBinLabel(5-ibin,cent_tag_pp_data[ibin]);
                h_jetchg_avg_data_sysm[ibin2]->SetBinContent(5-ibin,h_chg_data[ibin][ibin2]->GetMean());
                h_jetchg_avg_data_sysm[ibin2]->SetBinError(5-ibin,avg_err[ibin][ibin2]);

                h_jetchg_rms_data_sysm[ibin2]->GetXaxis()->SetBinLabel(5-ibin,cent_tag_pp_data[ibin]);
                h_jetchg_rms_data_sysm[ibin2]->SetBinContent(5-ibin,h_chg_data[ibin][ibin2]->GetRMS());
                h_jetchg_rms_data_sysm[ibin2]->SetBinError(5-ibin,rms_err[ibin][ibin2]);                
            }
        }
      }

      TCanvas *c_chg_abs_jetchg;
      c_chg_abs_jetchg = new TCanvas("c_chg_abs_jetchg","c_chg_abs_jetchg",1500,400);
      c_chg_abs_jetchg->Divide(5,1,0);
      gStyle->SetOptStat(0);
      gStyle->SetOptFit(0);

      for(int ibin=0; ibin<nCBins; ibin++){

          if(ibin==0) c_chg_abs_jetchg->cd(1); 
          else c_chg_abs_jetchg->cd(6-ibin); 

          h_chg_data_abs[ibin][0]->SetMarkerColor(kBlack); h_chg_data_abs[ibin][0]->SetMarkerStyle(20);    
          h_chg_data_abs[ibin][0]->GetXaxis()->SetRangeUser(-1.5,1.5);
          h_chg_data_abs[ibin][0]->GetYaxis()->SetRangeUser(-0.2,0.2);
          h_chg_data_abs[ibin][0]->Draw("e0 same");

/*  
          h_chg_reco_abs[ibin][0]->SetMarkerColor(kBlack); h_chg_reco_abs[ibin][0]->SetMarkerStyle(20);    
          h_chg_reco_abs[ibin][0]->GetXaxis()->SetRangeUser(-1.5,1.5);
          h_chg_reco_abs[ibin][0]->GetYaxis()->SetRangeUser(-0.2,0.2);
          h_chg_reco_abs[ibin][0]->Draw("e0 same");
*/
          h_chg_reco_up_abs[ibin][0]->Rebin(10);
          h_chg_reco_down_abs[ibin][0]->Rebin(10);
          //h_chg_reco_gubdb_abs[ibin][0]->Rebin(10);
          h_chg_reco_up_abs[ibin][0]->SetMarkerColor(kAzure-2); h_chg_reco_up_abs[ibin][0]->SetMarkerStyle(20); h_chg_reco_up_abs[ibin][0]->Draw("e0 same");
          h_chg_reco_down_abs[ibin][0]->SetMarkerColor(kGreen-2); h_chg_reco_down_abs[ibin][0]->SetMarkerStyle(20); h_chg_reco_down_abs[ibin][0]->Draw("e0 same");
          h_chg_reco_gubdb_abs[ibin][0]->SetMarkerColor(kRed-2); h_chg_reco_gubdb_abs[ibin][0]->SetMarkerStyle(20); h_chg_reco_gubdb_abs[ibin][0]->Draw("e0 same");
          f_chg_ud_data[ibin][0]->SetParameter(1,7.47768e-01);
          h_chg_data_abs[ibin][0]->Fit(f_chg_ud_data[ibin][0],"Q M R","sames",-fit_limit_ptcut[0],fit_limit_ptcut[0]);

            sprintf(g_frac, "g+others %.3f",f_chg_ud_data[ibin][0]->GetParameter(1));
            sprintf(chi2ndf, "chi2/NDF %.3f",f_chg_ud_data[ibin][0]->GetChisquare()/f_chg_ud_data[ibin][0]->GetNDF());
            sprintf(up_frac, "up %.3f",(f_chg_ud_data[ibin][0]->GetParameter(0)));
            sprintf(down_frac, "down %.3f",(1. - (f_chg_ud_data[ibin][0]->GetParameter(0)+f_chg_ud_data[ibin][0]->GetParameter(1))));

            tx1 = new TLatex(); tx1->SetTextSize(.07);
            tx1->DrawLatexNDC(0.15,0.9, g_frac);
            tx3 = new TLatex(); tx3->SetTextSize(.06);
            tx3->DrawLatexNDC(0.15,0.7, down_frac);
            tx4 = new TLatex(); tx4->SetTextSize(.07);
            tx4->DrawLatexNDC(0.15,0.8, up_frac);
            tx2 = new TLatex(); tx2->SetTextSize(.07);
            tx2->DrawLatexNDC(0.3,0.2, chi2ndf);

            TString tmp1 = cent_tag_pp_data[ibin];
            tx1 = new TLatex(); tx1->SetTextSize(.07);
            tx1->DrawLatexNDC(0.65,0.85, tmp1);
      }

      TCanvas *c_chg_avg_jetchg;
      c_chg_avg_jetchg = new TCanvas("c_chg_avg_jetchg","c_chg_avg_jetchg",1000,500);
      c_chg_avg_jetchg->Divide(2,1);

      c_chg_avg_jetchg->cd(1);
      TLegend *legendx = new TLegend(0.65,0.55,0.98,0.98);
      legendx ->SetLineColor(kWhite);

      for(int ibin2=0;ibin2<1;ibin2++){

        h_jetchg_avg_data[ibin2]->GetXaxis()->SetTitle("");
        h_jetchg_avg_data[ibin2]->GetYaxis()->SetTitle("Average Jet Charge [e]");
        h_cosm(h_jetchg_avg_data[ibin2]);
        h_jetchg_avg_data[ibin2]->GetYaxis()->SetRangeUser(-0.5,0.5);
        h_jetchg_avg_data[ibin2]->GetYaxis()->SetLabelSize(0.04);
        h_jetchg_avg_data[ibin2]->GetYaxis()->SetTitleOffset(0.85);
        h_jetchg_avg_data[ibin2]->SetLineColor(ibin2+1); h_jetchg_avg_data[ibin2]->SetMarkerColor(ibin2+1);
        h_jetchg_avg_data[ibin2]->SetMarkerStyle(34); h_jetchg_avg_data[ibin2]->SetMarkerSize(1.5);
        h_jetchg_avg_data_sysm[ibin2]->SetLineColor(ibin2+1); h_jetchg_avg_data_sysm[ibin2]->SetMarkerColor(ibin2+1); h_jetchg_avg_data_sysm[ibin2]->SetFillColorAlpha(kAzure+2,0.3);
        h_jetchg_avg_data_sysm[ibin2]->SetMarkerStyle(1); h_jetchg_avg_data_sysm[ibin2]->SetMarkerSize(1.);
        h_jetchg_avg_data_pyquen[ibin2]->SetLineColor(kRed-2); h_jetchg_avg_data_pyquen[ibin2]->SetMarkerColor(kRed-2); h_jetchg_avg_data_pyquen[ibin2]->SetLineStyle(2);
        h_jetchg_avg_data_pyquen[ibin2]->SetFillColorAlpha(kRed-2,0.3); h_jetchg_avg_data_pyquen[ibin2]->SetMarkerStyle(34); h_jetchg_avg_data_pyquen[ibin2]->SetMarkerSize(0.);
        h_jetchg_avg_MC[ibin2]->SetLineColor(kGreen-2); h_jetchg_avg_MC[ibin2]->SetMarkerColor(kGreen-2); h_jetchg_avg_MC[ibin2]->SetLineStyle(2);
        h_jetchg_avg_MC[ibin2]->SetMarkerStyle(34); h_jetchg_avg_MC[ibin2]->SetMarkerSize(1.);
        //h_jetchg_avg_data[ibin2]->Draw("p same");
        h_jetchg_avg_data_sysm[ibin2]->Draw("e2 same");
        //h_jetchg_avg_data_pyquen[ibin2]->Draw("e2 same");
        //h_jetchg_avg_MC[ibin2]->Draw("l same");        

        legendx ->AddEntry((TObject*)0, trk_tag[ibin2], "");
      }
      legendx ->Draw("same");

      TCanvas *c_chg_rms_jetchg;
      c_chg_rms_jetchg = new TCanvas("c_chg_rms_jetchg","c_chg_rms_jetchg",500,500);

      c_chg_rms_jetchg->cd(1);
      for(int ibin2=0;ibin2<1;ibin2++){

        h_jetchg_rms_data[ibin2]->GetXaxis()->SetTitle("");
        h_jetchg_rms_data[ibin2]->GetYaxis()->SetTitle("Jet Charge Standard Deviation [e]");
        h_cosm(h_jetchg_rms_data[ibin2]);
        h_jetchg_rms_data[ibin2]->GetYaxis()->SetRangeUser(0.3,0.5);
        h_jetchg_rms_data[ibin2]->GetYaxis()->SetLabelSize(0.04);
        h_jetchg_rms_data[ibin2]->GetYaxis()->SetTitleOffset(0.85);
        h_jetchg_rms_data[ibin2]->SetLineColor(ibin2+1); h_jetchg_rms_data[ibin2]->SetMarkerColor(ibin2+1);
        h_jetchg_rms_data[ibin2]->SetMarkerStyle(34); h_jetchg_rms_data[ibin2]->SetMarkerSize(1.5);
        h_jetchg_rms_data_sysm[ibin2]->SetLineColor(ibin2+1); h_jetchg_rms_data_sysm[ibin2]->SetMarkerColor(ibin2+1); h_jetchg_rms_data_sysm[ibin2]->SetFillColorAlpha(kAzure+2,0.3);
        h_jetchg_rms_data_sysm[ibin2]->SetMarkerStyle(34); h_jetchg_rms_data_sysm[ibin2]->SetMarkerSize(1.);
        h_jetchg_rms_data_pyquen[ibin2]->SetLineColor(kRed-2); h_jetchg_rms_data_pyquen[ibin2]->SetMarkerColor(kRed-2); h_jetchg_rms_data_pyquen[ibin2]->SetLineStyle(2);
        h_jetchg_rms_data_pyquen[ibin2]->SetFillColorAlpha(kRed-2,0.3); h_jetchg_rms_data_pyquen[ibin2]->SetMarkerStyle(34); h_jetchg_rms_data_pyquen[ibin2]->SetMarkerSize(0.);
        h_jetchg_rms_data_pyquencoll[ibin2]->SetLineColor(kRed-2); h_jetchg_rms_data_pyquencoll[ibin2]->SetMarkerColor(kRed-2); h_jetchg_rms_data_pyquencoll[ibin2]->SetLineStyle(2);
        h_jetchg_rms_data_pyquencoll[ibin2]->SetFillColorAlpha(kRed-2,0.3); h_jetchg_rms_data_pyquencoll[ibin2]->SetMarkerStyle(34); h_jetchg_rms_data_pyquencoll[ibin2]->SetMarkerSize(0.);
        h_jetchg_rms_MC[ibin2]->SetLineColor(kGreen-2); h_jetchg_rms_MC[ibin2]->SetMarkerColor(kGreen-2); h_jetchg_rms_MC[ibin2]->SetLineStyle(2);
        h_jetchg_rms_MC[ibin2]->SetMarkerStyle(34); h_jetchg_rms_MC[ibin2]->SetMarkerSize(1.);
        h_jetchg_rms_data[ibin2]->Draw("p same");
        h_jetchg_rms_data_sysm[ibin2]->Draw("e2 same");
        h_jetchg_rms_data_pyquen[ibin2]->Draw("e2 same");
        h_jetchg_rms_data_pyquencoll[ibin2]->Draw("e2 same");
        h_jetchg_rms_MC[ibin2]->Draw("l same");        
      }

      legendx ->Draw("same");
      c_chg_rms_jetchg->cd(1);
        TLegend *legendrms = new TLegend(0.62,0.62,0.98,0.98);
        legendrms->SetLineColor(kWhite);
        legendrms->AddEntry(h_jetchg_rms_data_sysm[0], "Data", "lepf");
        legendrms->AddEntry(h_jetchg_rms_MC[0], "Pythia+Hydjet", "lepf");
        legendrms->AddEntry(h_jetchg_rms_data_pyquen[0], "PYQUEN (Wide)", "f");
        legendrms->AddEntry(h_jetchg_rms_data_pyquencoll[0], "PYQUEN (Coll)", "f");
        legendrms->Draw("same");


///trkpt closure
    TCanvas *c_trkpt_MC;

    c_trkpt_MC = new TCanvas("c_trkpt_MC","c_trkpt_MC",1800,600);
    c_trkpt_MC->Divide(5,2,0);
    gStyle->SetOptStat(0);
    for(int ibin=0;ibin<nCBins;ibin++){
        if(ibin==0) c_trkpt_MC->cd(1); 
        else c_trkpt_MC->cd(6-ibin); 
        h_cosm(h_trkpt_reco[ibin]);
        h_trkpt_reco[ibin]->GetXaxis()->SetRangeUser(2.,10.); 
        h_trkpt_reco[ibin]->GetXaxis()->SetTitle("track p_{T}");
        h_trkpt_reco[ibin]->GetXaxis()->CenterTitle();
        h_trkpt_reco[ibin]->GetYaxis()->SetTitle("Entries");
        h_trkpt_reco[ibin]->GetYaxis()->CenterTitle();
        h_trkpt_reco[ibin]->SetLineColor(kBlack); h_trkpt_reco[ibin]->SetMarkerColor(kBlack); h_trkpt_reco[ibin]->SetMarkerStyle(20);
        h_trkpt_pos_reco[ibin]->SetLineColor(kBlue-2); h_trkpt_pos_reco[ibin]->SetMarkerColor(kBlue-2); h_trkpt_pos_reco[ibin]->SetMarkerStyle(5);
        h_trkpt_neg_reco[ibin]->SetLineColor(kRed-2); h_trkpt_neg_reco[ibin]->SetMarkerColor(kRed-2); h_trkpt_neg_reco[ibin]->SetMarkerStyle(5);

        h_trkpt_reco[ibin]->Draw("same");
        h_trkpt_neg_reco[ibin]->Draw("same");
        h_trkpt_pos_reco[ibin]->Draw("same");

        TString tmp1 = cent_tag_pp_MC[ibin];
        tx1 = new TLatex(); tx1->SetTextSize(.09);
        tx1->DrawLatexNDC(0.35,0.85, tmp1);

        if(ibin==0){
            TLegend *legendx = new TLegend(0.15,0.7,0.5,0.8);
            legendx ->SetLineColor(kWhite);
            legendx ->AddEntry(h_trkpt_reco[0], "all tracks (reco)", "lep");
            legendx ->AddEntry(h_trkpt_pos_reco[0], "positive tracks (reco)", "lep");
            legendx ->AddEntry(h_trkpt_neg_reco[0], "negative tracks (reco)", "lep");            
            legendx ->Draw("same");
        }


        if(ibin==0) c_trkpt_MC->cd(6); 
        else c_trkpt_MC->cd(11-ibin);

        h_trkpt_reco_gen[ibin] = (TH1D*)h_trkpt_reco[ibin]->Clone(Form("h_trkpt_reco_gen_%d",ibin)); 
        h_trkpt_pos_reco_gen[ibin] = (TH1D*)h_trkpt_pos_reco[ibin]->Clone(Form("h_trkpt_pos_reco_gen_%d",ibin)); 
        h_trkpt_neg_reco_gen[ibin] = (TH1D*)h_trkpt_neg_reco[ibin]->Clone(Form("h_trkpt_neg_reco_gen_%d",ibin)); 

        h_trkpt_reco_gen[ibin]->Divide(h_trkpt_gen[ibin]);
        h_trkpt_pos_reco_gen[ibin]->Divide(h_trkpt_pos_gen[ibin]);
        h_trkpt_neg_reco_gen[ibin]->Divide(h_trkpt_neg_gen[ibin]);

        h_cosm(h_trkpt_reco_gen[ibin]);

        h_trkpt_reco_gen[ibin]->GetYaxis()->SetRangeUser(0.5,1.); 
/*
        h_trkpt_reco_gen[ibin]->GetXaxis()->SetTitle("track p_{T}");
        h_trkpt_reco_gen[ibin]->GetXaxis()->CenterTitle();
*/
        h_trkpt_reco_gen[ibin]->GetYaxis()->SetTitle("Reco / Gen");
        h_trkpt_reco_gen[ibin]->GetYaxis()->CenterTitle();
        h_trkpt_reco_gen[ibin]->SetLineColor(kBlack); h_trkpt_reco_gen[ibin]->SetMarkerColor(kBlack); h_trkpt_reco_gen[ibin]->SetMarkerStyle(4);
        h_trkpt_pos_reco_gen[ibin]->SetLineColor(kBlue-2); h_trkpt_pos_reco_gen[ibin]->SetMarkerColor(kBlue-2); h_trkpt_pos_reco_gen[ibin]->SetMarkerStyle(4);
        h_trkpt_neg_reco_gen[ibin]->SetLineColor(kRed-2); h_trkpt_neg_reco_gen[ibin]->SetMarkerColor(kRed-2); h_trkpt_neg_reco_gen[ibin]->SetMarkerStyle(4);

        h_trkpt_reco_gen[ibin]->Draw("same");
        h_trkpt_neg_reco_gen[ibin]->Draw("same");
        h_trkpt_pos_reco_gen[ibin]->Draw("same");

        if(ibin==0){
            TLegend *legendy = new TLegend(0.15,0.7,0.5,0.8);
            legendy ->SetLineColor(kWhite);
            legendy ->AddEntry(h_trkpt_reco_gen[0], "all tracks", "lep");
            legendy ->AddEntry(h_trkpt_pos_reco_gen[0], "positive tracks", "lep");
            legendy ->AddEntry(h_trkpt_neg_reco_gen[0], "negative tracks", "lep");            
            legendy ->Draw("same");
        }
    }

///relative trk eff. error

    TCanvas *c_chg_trk_eff = new TCanvas("c_chg_trk_eff","c_chg_trk_eff",1500,300);
    c_chg_trk_eff->Divide(5,1,0);
    gStyle->SetOptStat(0);

    for(int ibin=0; ibin<nCBins; ibin++){

        if(ibin==0) c_chg_trk_eff->cd(1);
        else c_chg_trk_eff->cd(6-ibin);

        //h_chg_data_trkrel[ibin][0]->Scale(1./h_chg_data_trkrel[ibin][0]->Integral());

        h_chg_data_trkrel[ibin][0]->GetXaxis()->SetTitle("jetcharge");
        h_chg_data_trkrel[ibin][0]->GetYaxis()->SetTitle("rel. err");
        h_chg_data_trkrel[ibin][0]->GetXaxis()->SetTitleSize(0.05);
        h_chg_data_trkrel[ibin][0]->GetXaxis()->SetLabelSize(0.05);
        h_chg_data_trkrel[ibin][0]->GetXaxis()->CenterTitle();

        h_chg_data_trkrel[ibin][0]->GetXaxis()->SetRangeUser(-1.,1.);
        h_chg_data_trkrel[ibin][0]->GetYaxis()->SetRangeUser(0.,0.11);
        h_chg_data_trkrel[ibin][0]->SetLineColor(kRed); h_chg_data_trkrel[ibin][0]->SetLineStyle(1); h_chg_data_trkrel[ibin][0]->SetMarkerColor(kRed); h_chg_data_trkrel[ibin][0]->SetMarkerStyle(20);
        h_chg_data_trkrel[ibin][1]->SetLineColor(kBlue); h_chg_data_trkrel[ibin][1]->SetLineStyle(1); h_chg_data_trkrel[ibin][1]->SetMarkerColor(kBlue); h_chg_data_trkrel[ibin][1]->SetMarkerStyle(20);
        h_chg_data_trkrel[ibin][2]->SetLineColor(kGreen-2); h_chg_data_trkrel[ibin][2]->SetLineStyle(1); h_chg_data_trkrel[ibin][2]->SetMarkerColor(kGreen-2); h_chg_data_trkrel[ibin][2]->SetMarkerStyle(20);

        h_chg_data_trkrel[ibin][0]->Draw("hist same");
        h_chg_data_trkrel[ibin][1]->Draw("hist same");
        h_chg_data_trkrel[ibin][2]->Draw("hist same");

        TString tmp1 = cent_tag_pp_data[ibin];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.35,0.85, tmp1);

        if(ibin==1){
            TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
            legend ->SetLineColor(kWhite);
            legend ->AddEntry(h_chg_data_trkrel[ibin][0], "trk p_{T} > 2 GeV", "l");
            legend ->AddEntry(h_chg_data_trkrel[ibin][1], "trk p_{T} > 4 GeV", "l");
            legend ->AddEntry(h_chg_data_trkrel[ibin][2], "trk p_{T} > 5 GeV", "l");
            legend ->Draw("same");
        }
    }

//systematic uncertainties
    for(int ibin=0; ibin<nCBins; ibin++){
        sprintf(saythis,"h_syst_cent%d",ibin);
        h_syst[ibin] = new TH1D(saythis,"",3,trkpt_bounds);
        h_syst[ibin]->Sumw2();

        sprintf(saythis,"h_stat_cent%d",ibin);
        h_stat[ibin] = new TH1D(saythis,"",3,trkpt_bounds);
        h_stat[ibin]->Sumw2();

        sprintf(saythis,"h_syst_MCstat_cent%d",ibin);
        h_syst_MCstat[ibin] = new TH1D(saythis,"",3,trkpt_bounds);
        h_syst_MCstat[ibin]->Sumw2();

        sprintf(saythis,"h_syst_JER_cent%d",ibin);
        h_syst_JER[ibin] = new TH1D(saythis,"",3,trkpt_bounds);
        h_syst_JER[ibin]->Sumw2();

        sprintf(saythis,"h_syst_trk_cent%d",ibin);
        h_syst_trk[ibin] = new TH1D(saythis,"",3,trkpt_bounds);
        h_syst_trk[ibin]->Sumw2();

        sprintf(saythis,"h_syst_unf_cent%d",ibin);
        h_syst_unf[ibin] = new TH1D(saythis,"",3,trkpt_bounds);
        h_syst_unf[ibin]->Sumw2();

        sprintf(saythis,"h_syst_posneg_cent%d",ibin);
        h_syst_posneg[ibin] = new TH1D(saythis,"",3,trkpt_bounds);
        h_syst_posneg[ibin]->Sumw2();

        sprintf(saythis,"h_syst_0trk_cent%d",ibin);
        h_syst_0trk[ibin] = new TH1D(saythis,"",3,trkpt_bounds);
        h_syst_0trk[ibin]->Sumw2();

//kappa
        sprintf(saythis,"h_syst_kappa_cent%d",ibin);
        h_syst_kappa[ibin] = new TH1D(saythis,"",3,kappa_bounds);
        h_syst_kappa[ibin]->Sumw2();

        sprintf(saythis,"h_stat_kappa_cent%d",ibin);
        h_stat_kappa[ibin] = new TH1D(saythis,"",3,kappa_bounds);
        h_stat_kappa[ibin]->Sumw2();

        sprintf(saythis,"h_syst_MCstat_kappa_cent%d",ibin);
        h_syst_MCstat_kappa[ibin] = new TH1D(saythis,"",3,kappa_bounds);
        h_syst_MCstat_kappa[ibin]->Sumw2();

        sprintf(saythis,"h_syst_JER_kappa_cent%d",ibin);
        h_syst_JER_kappa[ibin] = new TH1D(saythis,"",3,kappa_bounds);
        h_syst_JER_kappa[ibin]->Sumw2();

        sprintf(saythis,"h_syst_trk_kappa_cent%d",ibin);
        h_syst_trk_kappa[ibin] = new TH1D(saythis,"",3,kappa_bounds);
        h_syst_trk_kappa[ibin]->Sumw2();

        sprintf(saythis,"h_syst_unf_kappa_cent%d",ibin);
        h_syst_unf_kappa[ibin] = new TH1D(saythis,"",3,kappa_bounds);
        h_syst_unf_kappa[ibin]->Sumw2();

        sprintf(saythis,"h_syst_posneg_kappa_cent%d",ibin);
        h_syst_posneg_kappa[ibin] = new TH1D(saythis,"",3,kappa_bounds);
        h_syst_posneg_kappa[ibin]->Sumw2();

        sprintf(saythis,"h_syst_0trk_kappa_cent%d",ibin);
        h_syst_0trk_kappa[ibin] = new TH1D(saythis,"",3,kappa_bounds);
        h_syst_0trk_kappa[ibin]->Sumw2();
    }


    Double_t gluon_relerr_ptcut_MCstat[nCBins][ntrkbins];
    Double_t gluon_relerr_ptcut_0trk[nCBins][ntrkbins];
    Double_t gluon_relerr_ptcut_trkunc[nCBins][ntrkbins];
    Double_t gluon_relerr_ptcut_posneg[nCBins][ntrkbins];  
    Double_t gluon_relerr_ptcut_JER[nCBins][ntrkbins];  
    Double_t gluon_relerr_ptcut_unf[nCBins][ntrkbins];  

    Double_t gluon_relerr_kappa_MCstat[nCBins][ntrkbins];
    Double_t gluon_relerr_kappa_0trk[nCBins][ntrkbins];
    Double_t gluon_relerr_kappa_trkunc[nCBins][ntrkbins];
    Double_t gluon_relerr_kappa_posneg[nCBins][ntrkbins];
    Double_t gluon_relerr_kappa_JER[nCBins][ntrkbins];
    Double_t gluon_relerr_kappa_unf[nCBins][ntrkbins];  

    Double_t gluon_stat_kappa[nCBins][ntrkbins];
    Double_t gluon_stat_ptcut[nCBins][ntrkbins];

    //ptcut
    gluon_relerr_ptcut_MCstat[0][0] = 1.;   gluon_relerr_ptcut_MCstat[0][1] = 1.;   gluon_relerr_ptcut_MCstat[0][2] = 1.;
    gluon_relerr_ptcut_MCstat[4][0] = 1.8;   gluon_relerr_ptcut_MCstat[4][1] = 1.8;   gluon_relerr_ptcut_MCstat[4][2] = 1.8;
    gluon_relerr_ptcut_MCstat[3][0] = 1.8;   gluon_relerr_ptcut_MCstat[3][1] = 1.8;   gluon_relerr_ptcut_MCstat[3][2] = 1.8;
    gluon_relerr_ptcut_MCstat[2][0] = 1.8;   gluon_relerr_ptcut_MCstat[2][1] = 1.8;   gluon_relerr_ptcut_MCstat[2][2] = 1.8;
    gluon_relerr_ptcut_MCstat[1][0] = 1.8;   gluon_relerr_ptcut_MCstat[1][1] = 1.8;   gluon_relerr_ptcut_MCstat[1][2] = 1.8;

    gluon_relerr_ptcut_0trk[0][0] = 0.1;   gluon_relerr_ptcut_0trk[0][1] = 0.1;   gluon_relerr_ptcut_0trk[0][2] = 0.2;
    gluon_relerr_ptcut_0trk[4][0] = 0.2;   gluon_relerr_ptcut_0trk[4][1] = 0.5;   gluon_relerr_ptcut_0trk[4][2] = 1.0;
    gluon_relerr_ptcut_0trk[3][0] = 0.4;   gluon_relerr_ptcut_0trk[3][1] = 1.0;   gluon_relerr_ptcut_0trk[3][2] = 2.0;
    gluon_relerr_ptcut_0trk[2][0] = 0.4;   gluon_relerr_ptcut_0trk[2][1] = 1.7;   gluon_relerr_ptcut_0trk[2][2] = 3.0;
    gluon_relerr_ptcut_0trk[1][0] = 0.4;   gluon_relerr_ptcut_0trk[1][1] = 2.4;   gluon_relerr_ptcut_0trk[1][2] = 4.5;

    gluon_relerr_ptcut_JER[0][0] = 1.0;   gluon_relerr_ptcut_JER[0][1] = 1.0;   gluon_relerr_ptcut_JER[0][2] = 1.0;
    gluon_relerr_ptcut_JER[4][0] = 2.0;   gluon_relerr_ptcut_JER[4][1] = 1.0;   gluon_relerr_ptcut_JER[4][2] = 1.0;
    gluon_relerr_ptcut_JER[3][0] = 2.0;   gluon_relerr_ptcut_JER[3][1] = 2.0;   gluon_relerr_ptcut_JER[3][2] = 2.0;
    gluon_relerr_ptcut_JER[2][0] = 2.0;   gluon_relerr_ptcut_JER[2][1] = 2.0;   gluon_relerr_ptcut_JER[2][2] = 2.0;
    gluon_relerr_ptcut_JER[1][0] = 2.0;   gluon_relerr_ptcut_JER[1][1] = 2.0;   gluon_relerr_ptcut_JER[1][2] = 2.0;

    gluon_relerr_ptcut_posneg[0][0] = 0.5;   gluon_relerr_ptcut_posneg[0][1] = 0.5;   gluon_relerr_ptcut_posneg[0][2] = 0.5;
    gluon_relerr_ptcut_posneg[4][0] = 1.;   gluon_relerr_ptcut_posneg[4][1] = 1.0;   gluon_relerr_ptcut_posneg[4][2] = 1.0;
    gluon_relerr_ptcut_posneg[3][0] = 1.;   gluon_relerr_ptcut_posneg[3][1] = 1.0;   gluon_relerr_ptcut_posneg[3][2] = 1.0;
    gluon_relerr_ptcut_posneg[2][0] = 1.;   gluon_relerr_ptcut_posneg[2][1] = 0.5;   gluon_relerr_ptcut_posneg[2][2] = 0.5;
    gluon_relerr_ptcut_posneg[1][0] = 1.;   gluon_relerr_ptcut_posneg[1][1] = 0.5;   gluon_relerr_ptcut_posneg[1][2] = 0.5;

    gluon_relerr_ptcut_unf[0][0] = 5.;   gluon_relerr_ptcut_unf[0][1] = 4.;   gluon_relerr_ptcut_unf[0][2] = 3.5;
    gluon_relerr_ptcut_unf[4][0] = 6.5;   gluon_relerr_ptcut_unf[4][1] = 5.5;   gluon_relerr_ptcut_unf[4][2] = 5.;
    gluon_relerr_ptcut_unf[3][0] = 6.5;   gluon_relerr_ptcut_unf[3][1] = 5.5;   gluon_relerr_ptcut_unf[3][2] = 5.;
    gluon_relerr_ptcut_unf[2][0] = 6.5;   gluon_relerr_ptcut_unf[2][1] = 5.5;   gluon_relerr_ptcut_unf[2][2] = 5.;
    gluon_relerr_ptcut_unf[1][0] = 6.5;   gluon_relerr_ptcut_unf[1][1] = 5.5;   gluon_relerr_ptcut_unf[1][2] = 5.;

    gluon_relerr_ptcut_trkunc[0][0] = 1.;   gluon_relerr_ptcut_trkunc[0][1] = 1.;   gluon_relerr_ptcut_trkunc[0][2] = 1.;
    gluon_relerr_ptcut_trkunc[4][0] = 2.;   gluon_relerr_ptcut_trkunc[4][1] = 2.;   gluon_relerr_ptcut_trkunc[4][2] = 2.;
    gluon_relerr_ptcut_trkunc[3][0] = 2.;   gluon_relerr_ptcut_trkunc[3][1] = 2.;   gluon_relerr_ptcut_trkunc[3][2] = 2.;
    gluon_relerr_ptcut_trkunc[2][0] = 2.;   gluon_relerr_ptcut_trkunc[2][1] = 2.;   gluon_relerr_ptcut_trkunc[2][2] = 2.;
    gluon_relerr_ptcut_trkunc[1][0] = 2.;   gluon_relerr_ptcut_trkunc[1][1] = 2.;   gluon_relerr_ptcut_trkunc[1][2] = 2.;

    gluon_stat_ptcut[0][0] = 3.;   gluon_stat_ptcut[0][1] = 3.;   gluon_stat_ptcut[0][2] = 3.;
    gluon_stat_ptcut[4][0] = 6.;   gluon_stat_ptcut[4][1] = 6.;   gluon_stat_ptcut[4][2] = 6.;
    gluon_stat_ptcut[3][0] = 6.;   gluon_stat_ptcut[3][1] = 6.;   gluon_stat_ptcut[3][2] = 6.;
    gluon_stat_ptcut[2][0] = 6.;   gluon_stat_ptcut[2][1] = 6.;   gluon_stat_ptcut[2][2] = 6.;
    gluon_stat_ptcut[1][0] = 6.;   gluon_stat_ptcut[1][1] = 6.;   gluon_stat_ptcut[1][2] = 6.;

    //kappa
    gluon_relerr_kappa_MCstat[0][0] = 1.;   gluon_relerr_kappa_MCstat[0][1] = 1.;   gluon_relerr_kappa_MCstat[0][2] = 1.;
    gluon_relerr_kappa_MCstat[4][0] = 1.8;   gluon_relerr_kappa_MCstat[4][1] = 1.8;   gluon_relerr_kappa_MCstat[4][2] = 1.8;
    gluon_relerr_kappa_MCstat[3][0] = 1.8;   gluon_relerr_kappa_MCstat[3][1] = 1.8;   gluon_relerr_kappa_MCstat[3][2] = 1.8;
    gluon_relerr_kappa_MCstat[2][0] = 1.8;   gluon_relerr_kappa_MCstat[2][1] = 1.8;   gluon_relerr_kappa_MCstat[2][2] = 1.8;
    gluon_relerr_kappa_MCstat[1][0] = 1.8;   gluon_relerr_kappa_MCstat[1][1] = 1.8;   gluon_relerr_kappa_MCstat[1][2] = 1.8;

    gluon_relerr_kappa_0trk[0][0] = 0.1;   gluon_relerr_kappa_0trk[0][1] = 0.1;   gluon_relerr_kappa_0trk[0][2] = 0.1;
    gluon_relerr_kappa_0trk[4][0] = 0.2;   gluon_relerr_kappa_0trk[4][1] = 0.2;   gluon_relerr_kappa_0trk[4][2] = 0.2;
    gluon_relerr_kappa_0trk[3][0] = 0.4;   gluon_relerr_kappa_0trk[3][1] = 0.4;   gluon_relerr_kappa_0trk[3][2] = 0.4;
    gluon_relerr_kappa_0trk[2][0] = 0.4;   gluon_relerr_kappa_0trk[2][1] = 0.4;   gluon_relerr_kappa_0trk[2][2] = 0.4;
    gluon_relerr_kappa_0trk[1][0] = 0.4;   gluon_relerr_kappa_0trk[1][1] = 0.4;   gluon_relerr_kappa_0trk[1][2] = 0.4;

    gluon_relerr_kappa_JER[0][0] = 1.0;   gluon_relerr_kappa_JER[0][1] = 1.0;   gluon_relerr_kappa_JER[0][2] = 1.5;
    gluon_relerr_kappa_JER[4][0] = 2.0;   gluon_relerr_kappa_JER[4][1] = 2.0;   gluon_relerr_kappa_JER[4][2] = 2.0;
    gluon_relerr_kappa_JER[3][0] = 2.0;   gluon_relerr_kappa_JER[3][1] = 2.0;   gluon_relerr_kappa_JER[3][2] = 2.0;
    gluon_relerr_kappa_JER[2][0] = 3.0;   gluon_relerr_kappa_JER[2][1] = 2.0;   gluon_relerr_kappa_JER[2][2] = 3.0;
    gluon_relerr_kappa_JER[1][0] = 2.0;   gluon_relerr_kappa_JER[1][1] = 2.0;   gluon_relerr_kappa_JER[1][2] = 3.0;

    gluon_relerr_kappa_posneg[0][0] = 1.;   gluon_relerr_kappa_posneg[0][1] = 0.5;   gluon_relerr_kappa_posneg[0][2] = 0.5;
    gluon_relerr_kappa_posneg[4][0] = 2.0;   gluon_relerr_kappa_posneg[4][1] = 1.;   gluon_relerr_kappa_posneg[4][2] = 1.;
    gluon_relerr_kappa_posneg[3][0] = 2.0;   gluon_relerr_kappa_posneg[3][1] = 1.;   gluon_relerr_kappa_posneg[3][2] = 1.;
    gluon_relerr_kappa_posneg[2][0] = 1.5;   gluon_relerr_kappa_posneg[2][1] = 1.;   gluon_relerr_kappa_posneg[2][2] = 1.;
    gluon_relerr_kappa_posneg[1][0] = 1.5;   gluon_relerr_kappa_posneg[1][1] = 1.;   gluon_relerr_kappa_posneg[1][2] = 1.;

    gluon_relerr_kappa_trkunc[0][0] = 1.;   gluon_relerr_kappa_trkunc[0][1] = 1.;   gluon_relerr_kappa_trkunc[0][2] = 1.;
    gluon_relerr_kappa_trkunc[4][0] = 2.;   gluon_relerr_kappa_trkunc[4][1] = 2.;   gluon_relerr_kappa_trkunc[4][2] = 2.;
    gluon_relerr_kappa_trkunc[3][0] = 2.;   gluon_relerr_kappa_trkunc[3][1] = 2.;   gluon_relerr_kappa_trkunc[3][2] = 2.;
    gluon_relerr_kappa_trkunc[2][0] = 2.;   gluon_relerr_kappa_trkunc[2][1] = 2.;   gluon_relerr_kappa_trkunc[2][2] = 2.;
    gluon_relerr_kappa_trkunc[1][0] = 2.;   gluon_relerr_kappa_trkunc[1][1] = 2.;   gluon_relerr_kappa_trkunc[1][2] = 2.;

    gluon_relerr_kappa_unf[0][0] = 6.;   gluon_relerr_kappa_unf[0][1] = 5.;   gluon_relerr_kappa_unf[0][2] = 4.;
    gluon_relerr_kappa_unf[4][0] = 7.5;   gluon_relerr_kappa_unf[4][1] = 6.5;   gluon_relerr_kappa_unf[4][2] = 5.;
    gluon_relerr_kappa_unf[3][0] = 7.5;   gluon_relerr_kappa_unf[3][1] = 6.5;   gluon_relerr_kappa_unf[3][2] = 5.;
    gluon_relerr_kappa_unf[2][0] = 7.5;   gluon_relerr_kappa_unf[2][1] = 6.5;   gluon_relerr_kappa_unf[2][2] = 5.;
    gluon_relerr_kappa_unf[1][0] = 7.5;   gluon_relerr_kappa_unf[1][1] = 6.5;   gluon_relerr_kappa_unf[1][2] = 5.;

    gluon_stat_kappa[0][0] = 3.;   gluon_stat_kappa[0][1] = 3.;   gluon_stat_kappa[0][2] = 3.;
    gluon_stat_kappa[4][0] = 6.;   gluon_stat_kappa[4][1] = 6.;   gluon_stat_kappa[4][2] = 6.;
    gluon_stat_kappa[3][0] = 6.;   gluon_stat_kappa[3][1] = 6.;   gluon_stat_kappa[3][2] = 6.;
    gluon_stat_kappa[2][0] = 6.;   gluon_stat_kappa[2][1] = 6.;   gluon_stat_kappa[2][2] = 6.;
    gluon_stat_kappa[1][0] = 6.;   gluon_stat_kappa[1][1] = 6.;   gluon_stat_kappa[1][2] = 6.;

    for(int ibin=0; ibin<nCBins; ibin++){
      for(int ibin2=0;ibin2<ntrkbins;ibin2++){ 
        if(do_ptcut_fitting){
/*
          h_syst_posneg[ibin]->SetBinContent(ibin2+1,fabs(gluon_withsys_posneg[ibin][ibin2] - gluon_withsys[ibin][ibin2])/gluon_withsys[ibin][ibin2]);
          h_syst_JER[ibin]->SetBinContent(ibin2+1,fabs(gluon_withsys_JER[ibin][ibin2] - gluon_withsys[ibin][ibin2])/gluon_withsys[ibin][ibin2]);
          if(ibin2==0) h_syst_0trk[ibin]->SetBinContent(ibin2+1,diff_0trk_data_MC[ibin][ibin2]/gluon_withsys[ibin][ibin2]);
          else h_syst_0trk[ibin]->SetBinContent(ibin2+1,diff_0trk_data_MC[ibin][ibin2+1]/gluon_withsys[ibin][ibin2]);

          if(ibin != 1)h_syst_trk[ibin]->SetBinContent(ibin2+1,(gluon_fit_high_diff[ibin][ibin2]+gluon_fit_low_diff[ibin][ibin2])/(2*gluon_withsys[ibin][ibin2]));
          else h_syst_trk[ibin]->SetBinContent(ibin2+1,(gluon_fit_high_diff[ibin+1][ibin2]+gluon_fit_low_diff[ibin+1][ibin2])/(2*gluon_withsys[ibin][ibin2]));
*/
          h_stat[ibin]->SetBinContent(ibin2+1,0.01*gluon_stat_ptcut[ibin][ibin2]);
          //h_stat[ibin]->SetBinContent(ibin2+1,gluon_fiterr[ibin][ibin2]/gluon_fit[ibin][ibin2]);
          h_syst_MCstat[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_ptcut_MCstat[ibin][ibin2]/gluon_fit[ibin][ibin2]);
          h_syst_posneg[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_ptcut_posneg[ibin][ibin2]);
          h_syst_JER[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_ptcut_JER[ibin][ibin2]);
          h_syst_0trk[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_ptcut_0trk[ibin][ibin2]);
          //h_syst_trk[ibin]->SetBinContent(ibin2+1,gluon_fit_high_diff[ibin][ibin2]/gluon_fit[ibin][ibin2]);
          h_syst_trk[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_ptcut_trkunc[ibin][ibin2]);
          //h_syst_unf[ibin]->SetBinContent(ibin2+1,gluon_fit_glusamp_diff[ibin][ibin2]/gluon_fit[ibin][ibin2]);
          h_syst_unf[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_ptcut_unf[ibin][ibin2]);
          //h_syst_trk[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_ptcut_trkunc[ibin][ibin2]);

          h_syst[ibin]->SetBinContent(ibin2+1,sqrt(pow(h_syst_MCstat[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_posneg[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_JER[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_trk[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_0trk[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_unf[ibin]->GetBinContent(ibin2+1),2)));
          gluon_fit_sysm[ibin][ibin2] = gluon_fit[ibin][ibin2]*h_syst[ibin]->GetBinContent(ibin2+1);
        }
        if(do_kappa_fitting){
/*
          h_syst_posneg_kappa[ibin]->SetBinContent(ibin2+1,fabs(gluon_withsys_kappa_posneg_kappa[ibin][ibin2] - gluon_withsys_kappa_kappa[ibin][ibin2])/gluon_withsys_kappa_kappa[ibin][ibin2]);
          h_syst_JER_kappa[ibin]->SetBinContent(ibin2+1,fabs(gluon_withsys_kappa_JER_kappa[ibin][ibin2] - gluon_withsys_kappa_kappa[ibin][ibin2])/gluon_withsys_kappa_kappa[ibin][ibin2]);
          h_syst_0trk_kappa[ibin]->SetBinContent(ibin2+1,0.0025);
          h_syst_trk_kappa[ibin]->SetBinContent(ibin2+1,(gluon_fit_kappa_high_diff_kappa[ibin][ibin2]+gluon_fit_kappa_low_diff_kappa[ibin][ibin2])/(2*gluon_withsys_kappa_kappa[ibin][ibin2]));
*/
          h_stat_kappa[ibin]->SetBinContent(ibin2+1,0.01*gluon_stat_kappa[ibin][ibin2]);
          //h_stat_kappa[ibin]->SetBinContent(ibin2+1,gluon_fiterr_kappa[ibin][ibin2]/gluon_fit_kappa[ibin][ibin2]);
          h_syst_MCstat_kappa[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_kappa_MCstat[ibin][ibin2]/gluon_fit_kappa[ibin][ibin2]);
          h_syst_posneg_kappa[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_kappa_posneg[ibin][ibin2]);
          //h_syst_posneg_kappa[ibin]->SetBinContent(ibin2+1,gluon_fit_posneg_diff_kappa[ibin][ibin2]/gluon_fit_kappa[ibin][ibin2]);
          h_syst_JER_kappa[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_kappa_JER[ibin][ibin2]);
          h_syst_0trk_kappa[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_kappa_0trk[ibin][ibin2]);
          h_syst_trk_kappa[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_kappa_trkunc[ibin][ibin2]);
          //h_syst_trk_kappa[ibin]->SetBinContent(ibin2+1,gluon_fit_high_diff_kappa[ibin][ibin2]/gluon_fit_kappa[ibin][ibin2]);
          h_syst_unf_kappa[ibin]->SetBinContent(ibin2+1,0.01*gluon_relerr_kappa_unf[ibin][ibin2]);

          h_syst_kappa[ibin]->SetBinContent(ibin2+1,sqrt(pow(h_syst_MCstat_kappa[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_posneg_kappa[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_JER_kappa[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_trk_kappa[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_0trk_kappa[ibin]->GetBinContent(ibin2+1),2)+pow(h_syst_unf_kappa[ibin]->GetBinContent(ibin2+1),2)));
          gluon_fit_sysm_kappa[ibin][ibin2] = gluon_fit_kappa[ibin][ibin2]*h_syst_kappa[ibin]->GetBinContent(ibin2+1);
        }
      }
    }

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
      if(do_ptcut_fitting) h_syst[ibin]->GetXaxis()->SetTitle("track p_{T} cut");
      else if(do_kappa_fitting) h_syst[ibin]->GetXaxis()->SetTitle("kappa");

      h_syst[ibin]->GetYaxis()->SetRangeUser(0.,0.2);
      h_syst[ibin]->SetLineColor(kBlack); h_syst[ibin]->SetMarkerColor(kBlack); h_syst[ibin]->SetMarkerStyle(20);
      h_syst[ibin]->SetLineWidth(2);
      h_syst[ibin]->Draw("L same");
      h_stat[ibin]->SetLineColor(kBlack); h_stat[ibin]->SetMarkerColor(kBlack); h_stat[ibin]->SetMarkerStyle(20);
      h_stat[ibin]->SetLineWidth(2); h_stat[ibin]->SetLineStyle(2);
      h_stat[ibin]->Draw("L same");
      h_syst_MCstat[ibin]->SetLineColor(kOrange+1); h_syst_MCstat[ibin]->SetMarkerColor(kOrange+1); h_syst_MCstat[ibin]->SetMarkerStyle(20);
      h_syst_MCstat[ibin]->SetLineWidth(2);
      h_syst_MCstat[ibin]->Draw("L same");
      h_syst_posneg[ibin]->SetLineColor(kGreen-2); h_syst_posneg[ibin]->SetMarkerColor(kGreen-2); h_syst_posneg[ibin]->SetMarkerStyle(20);
      h_syst_posneg[ibin]->SetLineWidth(2);
      h_syst_posneg[ibin]->Draw("L same");
      h_syst_JER[ibin]->SetLineColor(kRed-2); h_syst_JER[ibin]->SetMarkerColor(kRed-2); h_syst_JER[ibin]->SetMarkerStyle(20);
      h_syst_JER[ibin]->SetLineWidth(2);
      h_syst_JER[ibin]->Draw("L same");
      h_syst_0trk[ibin]->SetLineColor(kBlue-4); h_syst_0trk[ibin]->SetMarkerColor(kBlue-4); h_syst_0trk[ibin]->SetMarkerStyle(20);
      h_syst_0trk[ibin]->SetLineWidth(2);
      h_syst_0trk[ibin]->Draw("L same");
      h_syst_trk[ibin]->SetLineColor(kOrange); h_syst_trk[ibin]->SetMarkerColor(kOrange); h_syst_trk[ibin]->SetMarkerStyle(20);
      h_syst_trk[ibin]->SetLineWidth(2);
      h_syst_trk[ibin]->Draw("L same");
      h_syst_unf[ibin]->SetLineColor(kAzure+8); h_syst_unf[ibin]->SetMarkerColor(kAzure+8); h_syst_unf[ibin]->SetMarkerStyle(20);
      h_syst_unf[ibin]->SetLineWidth(2);
      h_syst_unf[ibin]->Draw("L same");

      TString tmp1 = cent_tag_pp_data[ibin];
      tx1 = new TLatex(); tx1->SetTextSize(.1);
      tx1->DrawLatexNDC(0.25,0.9, tmp1);

      if(ibin==0){
        TLegend *legend = new TLegend(0.15,0.35,0.75,0.85);
        legend ->SetLineColor(kWhite);
        legend ->AddEntry(h_syst[ibin], "total syst. unc.", "l");
        legend ->AddEntry(h_syst_trk[ibin], "tracking", "l");
        legend ->AddEntry(h_syst_JER[ibin], "JER", "l");
        legend ->AddEntry(h_syst_posneg[ibin], "pos/neg trk", "l");
        legend ->AddEntry(h_syst_0trk[ibin], "0 trk jets", "l");
        legend ->AddEntry(h_syst_MCstat[ibin], "MC stats", "l");
        legend ->AddEntry(h_syst_unf[ibin], "Unfolding", "l");
        legend ->Draw("same");

        TLegend *legendr = new TLegend(0.55,0.55,0.75,0.85);
        legendr->SetLineColor(kWhite);
        legendr ->AddEntry(h_stat[ibin], "stat. unc.", "l");
        legendr ->Draw("same");
      }    
    }

//kappa
    TCanvas *c_syst_kappa = new TCanvas("c_syst_kappa","c_syst_kappa",1500,350);
    drawgraphPad();

    pad_titleg->cd();
    //title_tex->DrawLatexNDC(0.5,0.1, pythia);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    pad_labelg->cd();
    chg_tex->DrawLatexNDC(0.7,0.4, "relative error");

    for(int ibin=0; ibin<nCBins; ibin++){
      if(ibin==0) pada->cd(1);
      else pada->cd(6-ibin);
      h_cosm(h_syst_kappa[ibin]);
      h_syst_kappa[ibin]->GetXaxis()->SetTitle("kappa");

      h_syst_kappa[ibin]->GetYaxis()->SetRangeUser(0.,0.2);
      h_syst_kappa[ibin]->SetLineColor(kBlack); h_syst_kappa[ibin]->SetMarkerColor(kBlack); h_syst_kappa[ibin]->SetMarkerStyle(20);
      h_syst_kappa[ibin]->SetLineWidth(2);
      h_syst_kappa[ibin]->Draw("L same");
      h_stat_kappa[ibin]->SetLineColor(kBlack); h_stat_kappa[ibin]->SetMarkerColor(kBlack); h_stat_kappa[ibin]->SetMarkerStyle(20);
      h_stat_kappa[ibin]->SetLineWidth(2); h_stat_kappa[ibin]->SetLineStyle(2);
      h_stat_kappa[ibin]->Draw("L same");
      h_syst_MCstat_kappa[ibin]->SetLineColor(kOrange+1); h_syst_MCstat_kappa[ibin]->SetMarkerColor(kOrange+1); h_syst_MCstat_kappa[ibin]->SetMarkerStyle(20);
      h_syst_MCstat_kappa[ibin]->SetLineWidth(2);
      h_syst_MCstat_kappa[ibin]->Draw("L same");
      h_syst_posneg_kappa[ibin]->SetLineColor(kGreen-2); h_syst_posneg_kappa[ibin]->SetMarkerColor(kGreen-2); h_syst_posneg_kappa[ibin]->SetMarkerStyle(20);
      h_syst_posneg_kappa[ibin]->SetLineWidth(2);
      h_syst_posneg_kappa[ibin]->Draw("L same");
      h_syst_JER_kappa[ibin]->SetLineColor(kRed-2); h_syst_JER_kappa[ibin]->SetMarkerColor(kRed-2); h_syst_JER_kappa[ibin]->SetMarkerStyle(20);
      h_syst_JER_kappa[ibin]->SetLineWidth(2);
      h_syst_JER_kappa[ibin]->Draw("L same");
      h_syst_0trk_kappa[ibin]->SetLineColor(kBlue-4); h_syst_0trk_kappa[ibin]->SetMarkerColor(kBlue-4); h_syst_0trk_kappa[ibin]->SetMarkerStyle(20);
      h_syst_0trk_kappa[ibin]->SetLineWidth(2);
      h_syst_0trk_kappa[ibin]->Draw("L same");
      h_syst_trk_kappa[ibin]->SetLineColor(kOrange); h_syst_trk_kappa[ibin]->SetMarkerColor(kOrange); h_syst_trk_kappa[ibin]->SetMarkerStyle(20);
      h_syst_trk_kappa[ibin]->SetLineWidth(2);
      h_syst_trk_kappa[ibin]->Draw("L same");
      h_syst_unf_kappa[ibin]->SetLineColor(kAzure+8); h_syst_unf_kappa[ibin]->SetMarkerColor(kAzure+8); h_syst_unf_kappa[ibin]->SetMarkerStyle(20);
      h_syst_unf_kappa[ibin]->SetLineWidth(2);
      h_syst_unf_kappa[ibin]->Draw("L same");

      TString tmp1 = cent_tag_pp_data[ibin];
      tx1 = new TLatex(); tx1->SetTextSize(.1);
      tx1->DrawLatexNDC(0.25,0.9, tmp1);

      if(ibin==0){
        TLegend *legend = new TLegend(0.15,0.35,0.75,0.85);
        legend ->SetLineColor(kWhite);
        legend ->AddEntry(h_syst_kappa[ibin], "total syst. unc.", "l");
        legend ->AddEntry(h_syst_trk_kappa[ibin], "tracking", "l");
        legend ->AddEntry(h_syst_JER_kappa[ibin], "JER", "l");
        legend ->AddEntry(h_syst_posneg_kappa[ibin], "pos/neg trk", "l");
        legend ->AddEntry(h_syst_0trk_kappa[ibin], "0 trk jets", "l");
        legend ->AddEntry(h_syst_MCstat_kappa[ibin], "MC stats", "l");
        legend ->AddEntry(h_syst_unf_kappa[ibin], "Unfolding", "l");
        legend ->Draw("same");

        TLegend *legendr = new TLegend(0.55,0.55,0.75,0.85);
        legendr->SetLineColor(kWhite);
        legendr ->AddEntry(h_stat_kappa[ibin], "stat. unc.", "l");
        legendr ->Draw("same");
      }    
    }

    for(int ibin=0; ibin<nCBins; ibin++){
        for(int ibin2=0; ibin2<ntrkbins; ibin2++){
            gluon_data_asymm_sys[ibin] = new TGraphErrors(ntrkbins,x_mean,gluon_fit[ibin],x_mean_err,gluon_fit_sysm[ibin]);
            gluon_data_asymm_sys[ibin]->SetName((TString)("gluon_data_sys_cent"+cent[ibin]));
            gluon_data_kappa_sysm[ibin] = new TGraphErrors(ntrkbins,kappa,gluon_fit_kappa[ibin],kappa_err,gluon_fit_sysm_kappa[ibin]);
            gluon_data_kappa_sysm[ibin]->SetName((TString)("gluon_data_sys_kappa_cent"+cent[ibin]));
        }
    }

//final plots

    TCanvas *c_data_qg_graph = new TCanvas("c_data_qg_graph","c_data_qg_graph",1500,350);
    drawgraphPad();
    gStyle->SetOptStat(0);

    pad_titleg->cd();
    title_tex->DrawLatexNDC(0.5,0.1, data);
    title_tex->DrawLatexNDC(0.05,0.1, cms);  
    pad_labelg->cd();
    chg_tex->DrawLatexNDC(0.7,0.4, fraction);

    for(int ibin=0; ibin<nCBins; ibin++){
      if(ibin==0) pada->cd(1);
      else pada->cd(6-ibin);
      h_data_ud_fit_results->Draw("same");
      
      if(do_ptcut_fitting){
          gluon_data[ibin]->SetMarkerColor(kRed); gluon_data[ibin]->SetLineColor(kRed); gluon_data[ibin]->SetMarkerStyle(4); gluon_data[ibin]->SetMarkerSize(1.2);
          gluon_data[ibin]->Draw("PL");
          gluon_data_asymm_sys[ibin]->SetMarkerColor(kRed); gluon_data_asymm_sys[ibin]->SetLineColor(kRed); gluon_data_asymm_sys[ibin]->SetFillColorAlpha(kRed,0.3); gluon_data_asymm_sys[ibin]->SetMarkerStyle(4); gluon_data_asymm_sys[ibin]->SetMarkerSize(1.2);
          gluon_data_asymm_sys[ibin]->Draw("3 same");
          gluon_data_old[ibin]->SetMarkerColor(kRed); gluon_data_old[ibin]->SetLineColor(kRed); gluon_data_old[ibin]->SetMarkerStyle(4); gluon_data_old[ibin]->SetMarkerSize(1.2);
          gluon_data_old[ibin]->SetLineStyle(2); //gluon_data_old[ibin]->Draw("PL");
          gluon_data_arc[ibin]->SetMarkerColor(kBlue-2); gluon_data_arc[ibin]->SetLineColor(kBlue-2); gluon_data_arc[ibin]->SetMarkerStyle(4); gluon_data_arc[ibin]->SetMarkerSize(1.2);
          //gluon_data_arc[ibin]->Draw("PL");
          //quark_data[ibin]->SetMarkerColor(kBlue); quark_data[ibin]->SetMarkerStyle(4); quark_data[ibin]->SetMarkerSize(1.2);
          //quark_data[ibin]->Draw("PL");
      }

      tl_g_ref[ibin]->Draw("same"); //tl_q[ibin]->Draw("same");
      
      TString tmp1 = cent_tag_pp_data[ibin];
      tx1 = new TLatex(); tx1->SetTextSize(.1);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);

      if(ibin==0){
        TLegend *legend = new TLegend(0.35,0.15,0.75,0.4);
        legend ->SetLineColor(kWhite);
        legend ->AddEntry(gluon_data_asymm_sys[ibin], "Data", "lepf");
        legend ->AddEntry(gluon_data_old[ibin], "MC", "l");
        legend ->Draw("same");
      }      
    }

    if(do_kappa_fitting){
        TCanvas *c_data_qg_graph_kappa = new TCanvas("c_data_qg_graph_kappa","c_data_qg_graph_kappa",1500,350);
        drawgraphPad();
        gStyle->SetOptStat(0);

        pad_titleg->cd();
        title_tex->DrawLatexNDC(0.5,0.1, data);
        title_tex->DrawLatexNDC(0.05,0.1, cms);  
        pad_labelg->cd();
        chg_tex->DrawLatexNDC(0.7,0.4, fraction);

        TH1D *h_data_ud_fit_results_kappa = new TH1D("h_data_ud_fit_results_kappa","",5,0.2,0.8);
        h_cosm(h_data_ud_fit_results_kappa);
        h_data_ud_fit_results_kappa->GetXaxis()->SetTitle("kappa");
        h_data_ud_fit_results_kappa->GetXaxis()->SetRangeUser(0.2,0.8);
        h_data_ud_fit_results_kappa->GetYaxis()->SetTitle("");
        h_data_ud_fit_results_kappa->GetYaxis()->SetRangeUser(0.,1.);

        for(int ibin=0; ibin<nCBins; ibin++){
          if(ibin==0) pada->cd(1);
          else pada->cd(6-ibin);
          h_data_ud_fit_results_kappa->Draw("same");
          
          gluon_data_kappa[ibin]->SetMarkerColor(kRed); gluon_data_kappa[ibin]->SetLineColor(kRed); gluon_data_kappa[ibin]->SetMarkerStyle(5); gluon_data_kappa[ibin]->SetMarkerSize(1.2);
          gluon_data_kappa[ibin]->Draw("PL");
          gluon_data_kappa_sysm[ibin]->SetMarkerColor(kRed); gluon_data_kappa_sysm[ibin]->SetLineColor(kRed); gluon_data_kappa_sysm[ibin]->SetFillColorAlpha(kRed,0.3); gluon_data_kappa_sysm[ibin]->SetMarkerStyle(5); gluon_data_kappa_sysm[ibin]->SetMarkerSize(1.2);
          gluon_data_kappa_sysm[ibin]->Draw("3 same");
          //quark_data_kappa[ibin]->SetMarkerColor(kBlue); quark_data_kappa[ibin]->SetMarkerStyle(5); quark_data_kappa[ibin]->SetMarkerSize(1.2);
          //quark_data_kappa[ibin]->Draw("PL");

          tl_g_ref_kappa[ibin]->Draw("same"); //tl_q[ibin]->Draw("same");
          
          TString tmp1 = cent_tag_pp_data[ibin];
          tx1 = new TLatex(); tx1->SetTextSize(.1);
          tx1->DrawLatexNDC(0.15,0.85, tmp1);

          if(ibin==0){
            TLegend *legend = new TLegend(0.35,0.15,0.75,0.4);
            legend ->SetLineColor(kWhite);
            legend ->AddEntry(gluon_data_kappa_sysm[ibin], "Data", "lepf");
            legend ->AddEntry(gluon_data_old[ibin], "MC", "l");
            legend ->Draw("same");
          }      
        }
    }















}

