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
#include "TDatime.h"

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
//#include "RooUnfoldIds.h"
#endif

#define nCBins 5
#define nptBins 3
#define ntrkbins 3
#define nkbins 3
#define ngenbins 50

TDatime* date = new TDatime();

char saythis[500];

TString cent_tag_pp_MC[] = {"PYTHIA ref","0-10% P+H","10-30% P+H","30-50% P+H","50-100% P+H"};
TString cent_tag_pp_data[] = {"pp ref","0-10% PbPb","10-30% PbPb","30-50% PbPb","50-100% PbPb"};
TString trk_tag[] = {"#kappa = 0.5 , p_{T}^{trk} > 2 GeV","#kappa = 0.5 , p_{T}^{trk} > 4 GeV", "#kappa = 0.5 , p_{T}^{trk} > 5 GeV"};
TString kappa_tag[] = {"#kappa = 0.3 , p_{T}^{trk} > 2 GeV","#kappa = 0.5 , p_{T}^{trk} > 2 GeV", "#kappa = 0.7 , p_{T}^{trk} > 2 GeV"};

void unfold_jetcharge_fulleta(){

//ptcut and kappa
    TFile *jetcharge_pp_MC = TFile::Open("../../../Documents/research/jetcharge/Pythia6_fulleta_trkeffconst_jetchg_no0trkjets_20190415.root");
    //TFile *jetcharge_pp_MC_tounfold = TFile::Open("../../../Documents/research/jetcharge/Pythia6_fulleta_trkeffconst_split_jetchg_no0trkjets_20190429.root");
    TFile *jetcharge_pp_data = TFile::Open("../../../Documents/research/jetcharge/ppdata_fulleta_jetchg_trkeffconst_recogen_no0trkjets_20190415.root");
    //TFile *jetcharge_pp_data = TFile::Open("../../../Documents/research/jetcharge/Herwig_fulleta_trkeffconst_jetchg_no0trkjets_20190429.root");
    TFile *jetcharge_PbPb_MC = TFile::Open("../../../Documents/research/jetcharge/P+H_fulleta_trkeffconst_jetchg_no0trkjets_20190415.root");
    //TFile *jetcharge_PbPb_MC_tounfold = TFile::Open("../../../Documents/research/jetcharge/P+H_fulleta_trkeffconst_split_jetchg_no0trkjets_20190430.root");
    TFile *jetcharge_PbPb_data = TFile::Open("../../../Documents/research/jetcharge/PbPbdata_fulleta_jetchg_trkeffconst_no0trkjets_20190415.root");

    TFile *jetcharge_pp_MC_vary = TFile::Open("../../../Documents/research/jetcharge/Pythia6_fulleta_trkeffconst_jetchg_no0trkjets_20190415.root");
    TFile *jetcharge_PbPb_MC_vary = TFile::Open("../../../Documents/research/jetcharge/P+H_fulleta_trkeffconst_jetchg_no0trkjets_20190415.root");

/*
    TFile *jetcharge_pp_MC_vary = TFile::Open("../../../Documents/research/jetcharge/Pythia6_fulleta_trkeffdec_jetchg_no0trkjets_20190416.root");
    TFile *jetcharge_PbPb_MC_vary = TFile::Open("../../../Documents/research/jetcharge/P+H_fulleta_trkeffdec_jetchg_no0trkjets_20190416.root");
*/
//ptcut
    TH2F *h_chg_recoptsube0_genptsube0[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_recoptsube0[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_genptsube0[nCBins][ntrkbins];    
    TH2F *h_chg_effgenpt_genptsube0[nCBins][ntrkbins];    
    TH2F *h_chg_effgenpt_genptsube0_trkvary[nCBins][ntrkbins];    
    TH2F *h_chg_recorecopt_gengensube0pt[nCBins][ntrkbins];    

    TH2F *h_chg_recoptsube0_genptsube0_up[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_recoptsube0_up[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_genptsube0_up[nCBins][ntrkbins];    
    TH2F *h_chg_effgenpt_genptsube0_up[nCBins][ntrkbins];    
    TH2F *h_chg_effgenpt_genptsube0_trkvary_up[nCBins][ntrkbins];    
    TH2F *h_chg_recorecopt_gengensube0pt_up[nCBins][ntrkbins];    

    TH2F *h_chg_recoptsube0_genptsube0_d[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_recoptsube0_d[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_genptsube0_d[nCBins][ntrkbins];    
    TH2F *h_chg_effgenpt_genptsube0_d[nCBins][ntrkbins];    
    TH2F *h_chg_effgenpt_genptsube0_trkvary_d[nCBins][ntrkbins];    
    TH2F *h_chg_recorecopt_gengensube0pt_d[nCBins][ntrkbins];    

    TH2F *h_chg_recoptsube0_genptsube0_g_oth[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_recoptsube0_g_oth[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_genptsube0_g_oth[nCBins][ntrkbins];    
    TH2F *h_chg_effgenpt_genptsube0_g_oth[nCBins][ntrkbins];    
    TH2F *h_chg_effgenpt_genptsube0_trkvary_g_oth[nCBins][ntrkbins];    
    TH2F *h_chg_recorecopt_gengensube0pt_g_oth[nCBins][ntrkbins];    

//kappa
    TH2F *h_chg_recopt_genpt_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recoptsube0_genptsube0_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_recoptsube0_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_genptsube0_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_effgenpt_genptsube0_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recorecopt_gengensube0pt_kappa[nCBins][ntrkbins];    

    TH2F *h_chg_recopt_genpt_up_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recoptsube0_genptsube0_up_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_recoptsube0_up_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_genptsube0_up_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_effgenpt_genptsube0_up_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recorecopt_gengensube0pt_up_kappa[nCBins][ntrkbins];    

    TH2F *h_chg_recopt_genpt_d_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recoptsube0_genptsube0_d_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_recoptsube0_d_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_genptsube0_d_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_effgenpt_genptsube0_d_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recorecopt_gengensube0pt_d_kappa[nCBins][ntrkbins];    

    TH2F *h_chg_recopt_genpt_g_oth_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recoptsube0_genptsube0_g_oth_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_recoptsube0_g_oth_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recopt_genptsube0_g_oth_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_effgenpt_genptsube0_g_oth_kappa[nCBins][ntrkbins];    
    TH2F *h_chg_recorecopt_gengensube0pt_g_oth_kappa[nCBins][ntrkbins];

//ptcut
    TH1D *h_chg_data[nCBins][ntrkbins];    
    TH1D *h_chg_data_trkunfold[nCBins][ntrkbins];    
    TH1D *h_chg_data_bkgunfold[nCBins][ntrkbins];    

    TH1D *h_chg_reco_prior[nCBins][ntrkbins][ngenbins];
    TH1D *h_chg_reco_unfolded[nCBins][ntrkbins];
    TH1D *h_chg_reco_unfolded_forbkg[nCBins][ntrkbins];

    TH1D *h_chg_reco_unfolded_rand[nCBins][ntrkbins];
    TH1D *h_chg_data_unfolded_rand[nCBins][ntrkbins];

    TH1D *h_chg_reco_unfolded_forbkg_rand[nCBins][ntrkbins];
    TH1D *h_chg_data_unfolded_forbkg_rand[nCBins][ntrkbins];

    TH1D *h_chg_reco_full_sube0[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_sube0[nCBins][ntrkbins];

    TH1D *h_chg_reco_full_sube0_up[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_sube0_up[nCBins][ntrkbins];
    TH1D *h_chg_reco_full_sube0_d[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_sube0_d[nCBins][ntrkbins];
    TH1D *h_chg_reco_full_sube0_g_oth[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_sube0_g_oth[nCBins][ntrkbins];

    TH1D *h_chg_reco[nCBins][ntrkbins];
    TH1D *h_chg_gen[nCBins][ntrkbins];

    TH1D *h_chg_reco_forbkg[nCBins][ntrkbins];
    TH1D *h_chg_gen_forbkg[nCBins][ntrkbins];

    TH1D *h_chg_reco_forbkg_prior[nCBins][ntrkbins][ngenbins];

    TH1D *h_chg_reco_full[nCBins][ntrkbins];
    TH1D *h_chg_reco_full_tounfold[nCBins][ntrkbins];    
    TH1D *h_chg_gen_full_forbkg[nCBins][ntrkbins];

    TH1D *h_chg_reco_full_up_tounfold[nCBins][ntrkbins];    
    TH1D *h_chg_reco_full_d_tounfold[nCBins][ntrkbins];    
    TH1D *h_chg_reco_full_g_oth_tounfold[nCBins][ntrkbins];    

    TH1D *h_chg_reco_full_trkvary[nCBins][ntrkbins];
    TH1D *h_chg_reco_full_trkvary_up[nCBins][ntrkbins];
    TH1D *h_chg_reco_full_trkvary_d[nCBins][ntrkbins];
    TH1D *h_chg_reco_full_trkvary_g_oth[nCBins][ntrkbins];

    TH1D *h_chg_reco_full_up[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_forbkg_up[nCBins][ntrkbins];
    TH1D *h_chg_reco_full_d[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_forbkg_d[nCBins][ntrkbins];
    TH1D *h_chg_reco_full_g_oth[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_forbkg_g_oth[nCBins][ntrkbins];

    TH1D *h_chg_eff_gen_full[nCBins][ntrkbins];
    TH1D *h_chg_eff_gen_full_up[nCBins][ntrkbins];
    TH1D *h_chg_eff_gen_full_d[nCBins][ntrkbins];
    TH1D *h_chg_eff_gen_full_g_oth[nCBins][ntrkbins];

    TH1D *h_chg_eff_gen_full_sube0[nCBins][ntrkbins];
    TH1D *h_chg_eff_gen_full_sube0_up[nCBins][ntrkbins];
    TH1D *h_chg_eff_gen_full_sube0_d[nCBins][ntrkbins];
    TH1D *h_chg_eff_gen_full_sube0_g_oth[nCBins][ntrkbins];

    TH1D *h_chg_gengensube0_full[nCBins][ntrkbins];
    TH1D *h_chg_gengensube0_full_up[nCBins][ntrkbins];
    TH1D *h_chg_gengensube0_full_d[nCBins][ntrkbins];
    TH1D *h_chg_gengensube0_full_g_oth[nCBins][ntrkbins];

    TH1D *h_chg_reco_gen_full[nCBins][ntrkbins];
    TH1D *h_chg_reco_unfolded_gen[nCBins][ntrkbins];
    TH1D *h_chg_data_unfolded_gen[nCBins][ntrkbins];
    TH1D *h_chg_reco_unfolded_rand_gen[nCBins][ntrkbins];
    TH1D *h_chg_data_unfolded_rand_gen[nCBins][ntrkbins];
    TH1D *h_chg_mc_data[nCBins][ntrkbins];

    TH1D *h_chg_reco_gen_full_forbkg[nCBins][ntrkbins];
    TH1D *h_chg_reco_unfolded_gen_forbkg[nCBins][ntrkbins];
    TH1D *h_chg_reco_unfolded_rand_gen_forbkg[nCBins][ntrkbins];

//kappa
    TH1D *h_chg_data_kappa[nCBins][ntrkbins];    
    TH1D *h_chg_data_trkunfold_kappa[nCBins][ntrkbins];    
    TH1D *h_chg_data_bkgunfold_kappa[nCBins][ntrkbins];    

    TH1D *h_chg_reco_unfolded_kappa[nCBins][ntrkbins];
    TH1D *h_chg_reco_unfolded_forbkg_kappa[nCBins][ntrkbins];

    TH1D *h_chg_reco_full_sube0_kappa[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_sube0_kappa[nCBins][ntrkbins];

    TH1D *h_chg_reco_full_sube0_up_kappa[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_sube0_up_kappa[nCBins][ntrkbins];
    TH1D *h_chg_reco_full_sube0_d_kappa[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_sube0_d_kappa[nCBins][ntrkbins];
    TH1D *h_chg_reco_full_sube0_g_oth_kappa[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_sube0_g_oth_kappa[nCBins][ntrkbins];

    TH1D *h_chg_reco_kappa[nCBins][ntrkbins];
    TH1D *h_chg_gen_kappa[nCBins][ntrkbins];

    TH1D *h_chg_reco_forbkg_kappa[nCBins][ntrkbins];
    TH1D *h_chg_gen_forbkg_kappa[nCBins][ntrkbins];

    TH1D *h_chg_reco_full_kappa[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_forbkg_kappa[nCBins][ntrkbins];

    TH1D *h_chg_reco_full_up_kappa[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_forbkg_up_kappa[nCBins][ntrkbins];
    TH1D *h_chg_reco_full_d_kappa[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_forbkg_d_kappa[nCBins][ntrkbins];
    TH1D *h_chg_reco_full_g_oth_kappa[nCBins][ntrkbins];
    TH1D *h_chg_gen_full_forbkg_g_oth_kappa[nCBins][ntrkbins];

    TH1D *h_chg_eff_gen_full_kappa[nCBins][ntrkbins];
    TH1D *h_chg_eff_gen_full_up_kappa[nCBins][ntrkbins];
    TH1D *h_chg_eff_gen_full_d_kappa[nCBins][ntrkbins];
    TH1D *h_chg_eff_gen_full_g_oth_kappa[nCBins][ntrkbins];

    TH1D *h_chg_eff_gen_full_sube0_kappa[nCBins][ntrkbins];
    TH1D *h_chg_eff_gen_full_sube0_up_kappa[nCBins][ntrkbins];
    TH1D *h_chg_eff_gen_full_sube0_d_kappa[nCBins][ntrkbins];
    TH1D *h_chg_eff_gen_full_sube0_g_oth_kappa[nCBins][ntrkbins];

    TH1D *h_chg_gengensube0_full_kappa[nCBins][ntrkbins];
    TH1D *h_chg_gengensube0_full_up_kappa[nCBins][ntrkbins];
    TH1D *h_chg_gengensube0_full_d_kappa[nCBins][ntrkbins];
    TH1D *h_chg_gengensube0_full_g_oth_kappa[nCBins][ntrkbins];

    TH1D *h_chg_reco_gen_full_kappa[nCBins][ntrkbins];
    TH1D *h_chg_reco_unfolded_gen_kappa[nCBins][ntrkbins];
    TH1D *h_chg_data_unfolded_gen_kappa[nCBins][ntrkbins];
    TH1D *h_chg_mc_data_kappa[nCBins][ntrkbins];

    TH1D *h_chg_reco_gen_full_forbkg_kappa[nCBins][ntrkbins];
    TH1D *h_chg_reco_unfolded_gen_forbkg_kappa[nCBins][ntrkbins];
    TH1D *h_chg_data_unfolded_gen_forbkg[nCBins][ntrkbins];
    TH1D *h_chg_reco_unfolded_rand_gen_forbkg_kappa[nCBins][ntrkbins];

//gen
//    double up = 0.126;
//    double down = 0.143;
//    double g_oth = 0.731;

//reco

    double up = 0.763;
    double down = 1.379;
    double g_oth = 1.;
    //double g_oth = 1.184;
    //double g_oth = 0.631;

/*
    double up = 0.13;
    double down = 0.13;
    double g_oth = 1.32;    
*/
/*
    double up = 0.;
    double down = 0.;
    double g_oth = 1.;
*/
    for (int ibin=0; ibin<nCBins; ibin++){

	    for (int ibin2=0; ibin2<ntrkbins; ibin2++){
			sprintf(saythis,"h_chg_reco_unfolded_rand_cent%d_trk%d",ibin,ibin2);
			h_chg_reco_unfolded_rand[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
			h_chg_reco_unfolded_rand[ibin][ibin2]->Sumw2();

			sprintf(saythis,"h_chg_data_unfolded_rand_cent%d_trk%d",ibin,ibin2);
			h_chg_data_unfolded_rand[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
			h_chg_data_unfolded_rand[ibin][ibin2]->Sumw2();

			sprintf(saythis,"h_chg_reco_unfolded_forbkg_rand_cent%d_trk%d",ibin,ibin2);
			h_chg_reco_unfolded_forbkg_rand[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
			h_chg_reco_unfolded_forbkg_rand[ibin][ibin2]->Sumw2();

			sprintf(saythis,"h_chg_data_unfolded_forbkg_rand_cent%d_trk%d",ibin,ibin2);
			h_chg_data_unfolded_forbkg_rand[ibin][ibin2] = new TH1D(saythis,"",50,-2.555,2.445);
			h_chg_data_unfolded_forbkg_rand[ibin][ibin2]->Sumw2();
		}

        if(ibin==0){
	        for (int ibin2=0; ibin2<ntrkbins; ibin2++){
//ptcut
                h_chg_recoptsube0_genptsube0[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get(Form("h_chg_recopt_genpt_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_genpt_cent%d_trk%d",ibin,ibin2));
/*	  
                h_chg_recoptsube0_genptsube0[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get(Form("h_chg_recoptsube0_genptsube0_g_oth_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_genpt_cent%d_trk%d",ibin,ibin2));
                h_chg_recoptsube0_genptsube0_up[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get(Form("h_chg_recoptsube0_genptsube0_up_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_genpt_up_cent%d_trk%d",ibin,ibin2));
                h_chg_recoptsube0_genptsube0_d[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get(Form("h_chg_recoptsube0_genptsube0_d_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_genpt_d_cent%d_trk%d",ibin,ibin2));
                h_chg_recoptsube0_genptsube0_up[ibin][ibin2]->Scale(0.13);
                h_chg_recoptsube0_genptsube0_d[ibin][ibin2]->Scale(0.135);
                h_chg_recoptsube0_genptsube0[ibin][ibin2]->Scale(1.32);
                h_chg_recoptsube0_genptsube0[ibin][ibin2]->Add(h_chg_recoptsube0_genptsube0_d[ibin][ibin2]);
                h_chg_recoptsube0_genptsube0[ibin][ibin2]->Add(h_chg_recoptsube0_genptsube0_up[ibin][ibin2]);
*/
                h_chg_effgenpt_genptsube0[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get(Form("h_chg_effgenpt_genptsube0_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_cent%d_trk%d",ibin,ibin2));
                //h_chg_effgenpt_genptsube0_trkvary[ibin][ibin2] = (TH2F*)jetcharge_pp_MC_vary->Get(Form("h_chg_effgenpt_genptsube0_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_vary_cent%d_trk%d",ibin,ibin2));
                h_chg_recorecopt_gengensube0pt[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get(Form("h_chg_recorecopt_gengensube0pt_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recorecopt_gengensube0pt_cent%d_trk%d",ibin,ibin2));
                h_chg_data[ibin][ibin2] = (TH1D*)jetcharge_pp_data->Get(Form("h_chg_recopt_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_data_cent%d_trk%d",ibin,ibin2)); 

                h_chg_reco_full[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_recopt_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_%d_%d",ibin,ibin2));           	        
/*
                h_chg_reco_full[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_recopt_g_oth_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_%d_%d",ibin,ibin2));       
                h_chg_reco_full_up[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_recopt_up_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_up_%d_%d",ibin,ibin2));       
                h_chg_reco_full_d[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_recopt_d_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_d_%d_%d",ibin,ibin2));       
                h_chg_reco_full_up[ibin][ibin2]->Scale(0.13);
                h_chg_reco_full_d[ibin][ibin2]->Scale(0.135);
                h_chg_reco_full[ibin][ibin2]->Scale(1.32);
                h_chg_reco_full[ibin][ibin2]->Add(h_chg_reco_full_d[ibin][ibin2]);
                h_chg_reco_full[ibin][ibin2]->Add(h_chg_reco_full_up[ibin][ibin2]);
*/
                h_chg_reco_full_tounfold[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_recopt_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_%d_%d",ibin,ibin2));       
                //h_chg_reco_full_sube0[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_recopt_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_%d_%d",ibin,ibin2));       
                //h_chg_reco_full_trkvary[ibin][ibin2] = (TH1D*)jetcharge_pp_MC_vary->Get(Form("h_chg_eff_genpt_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_%d_%d",ibin,ibin2));       

                h_chg_gen_full_sube0[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_genpt_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_genpt_%d_%d",ibin,ibin2));       
/*
	            h_chg_gen_full_sube0[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_genpt_g_oth_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_genpt_%d_%d",ibin,ibin2));       
                h_chg_gen_full_sube0_up[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_genpt_up_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_genpt_up_%d_%d",ibin,ibin2));       
                h_chg_gen_full_sube0_d[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_genpt_d_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_genpt_d_%d_%d",ibin,ibin2));       
                h_chg_gen_full_sube0_up[ibin][ibin2]->Scale(0.13);
                h_chg_gen_full_sube0_d[ibin][ibin2]->Scale(0.135);
                h_chg_gen_full_sube0[ibin][ibin2]->Scale(1.32);
                h_chg_gen_full_sube0[ibin][ibin2]->Add(h_chg_gen_full_sube0_d[ibin][ibin2]);
                h_chg_gen_full_sube0[ibin][ibin2]->Add(h_chg_gen_full_sube0_up[ibin][ibin2]);
*/
                h_chg_eff_gen_full_sube0[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_eff_genpt_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_eff_genpt_%d_%d",ibin,ibin2));       
                h_chg_gengensube0_full[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_gengensube0pt_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_gengensube0pt_%d_%d",ibin,ibin2));       
                h_chg_eff_gen_full_sube0[ibin][ibin2]->Rebin(10);

                h_chg_gen_full_sube0_g_oth[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_genpt_sube0_g_oth_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_genpt_sube0_g_oth_%d_%d",ibin,ibin2));

//kappa
                h_chg_recoptsube0_genptsube0_kappa[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get(Form("h_chg_recopt_genpt_kappa_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_genpt_kappa_cent%d_trk%d",ibin,ibin2));
                //h_chg_recopt_genpt_kappa[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get(Form("h_chg_recopt_genpt_kappa_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_genpt_kappa_cent%d_trk%d",ibin,ibin2));
                h_chg_effgenpt_genptsube0_kappa[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get(Form("h_chg_effgenpt_genptsube0_kappa_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_kappa_cent%d_trk%d",ibin,ibin2));
                h_chg_recorecopt_gengensube0pt_kappa[ibin][ibin2] = (TH2F*)jetcharge_pp_MC->Get(Form("h_chg_recorecopt_gengensube0pt_kappa_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recorecopt_gengensube0pt_kappa_cent%d_trk%d",ibin,ibin2));
                h_chg_data_kappa[ibin][ibin2] = (TH1D*)jetcharge_pp_data->Get(Form("h_chg_recopt_kappa_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_data_kappa_cent%d_trk%d",ibin,ibin2)); 

                h_chg_reco_full_kappa[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_recopt_kappa_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_kappa_%d_%d",ibin,ibin2));       
                //h_chg_reco_full_sube0_kappa[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_recopt_kappa_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_recopt_kappa_%d_%d",ibin,ibin2));       
                h_chg_gen_full_sube0_kappa[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_genpt_kappa_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_genpt_kappa_%d_%d",ibin,ibin2));       
                h_chg_eff_gen_full_sube0_kappa[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_eff_genpt_kappa_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_eff_genpt_kappa_%d_%d",ibin,ibin2));       
                h_chg_gengensube0_full_kappa[ibin][ibin2] = (TH1D*)jetcharge_pp_MC->Get(Form("h_chg_gengensube0pt_kappa_cent%d_trk%d",ibin,ibin2))->Clone(Form("h_chg_gengensube0pt_kappa_%d_%d",ibin,ibin2));       
                h_chg_eff_gen_full_sube0_kappa[ibin][ibin2]->Rebin(10);
	        }
	    }
	    else{
	        for (int ibin2=0; ibin2<ntrkbins; ibin2++){
//up
	            h_chg_recoptsube0_genptsube0_up[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recoptsube0_genptsube0_up_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_genpt_up_cent%d_trk%d",ibin,ibin2));       
	            h_chg_recopt_recoptsube0_up[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_recoptsube0_up_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_recoptsube0_up_%d_%d",ibin,ibin2));       
	            h_chg_recopt_genptsube0_up[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_genptsube0_up_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_genptsube0_up_%d_%d",ibin,ibin2));       
                h_chg_effgenpt_genptsube0_up[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_effgenpt_genptsube0_up_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_up_cent%d_trk%d",ibin,ibin2));
                //h_chg_effgenpt_genptsube0_trkvary_up[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC_vary->Get(Form("h_chg_effgenpt_genptsube0_up_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_up_trkvary_cent%d_trk%d",ibin,ibin2));
	
                //h_chg_reco_full_trkvary_up[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC_vary->Get(Form("h_chg_eff_genpt_up_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_up_trkvary_%d_%d",ibin,ibin2));       
	            h_chg_reco_full_up[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_up_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_up_%d_%d",ibin,ibin2));       
                h_chg_reco_full_up_tounfold[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_up_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_up_datasample_%d_%d",ibin,ibin2));       
                h_chg_eff_gen_full_up[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_eff_genpt_up_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_eff_genpt_up_%d_%d",ibin,ibin2));       

    	        h_chg_reco_full_sube0_up[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_sube0_up_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_sube0_up_%d_%d",ibin,ibin2));
	            h_chg_gen_full_sube0_up[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_genpt_sube0_up_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_genpt_sube0_up_%d_%d",ibin,ibin2));

//down
	            h_chg_recoptsube0_genptsube0_d[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recoptsube0_genptsube0_d_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_genpt_d_cent%d_trk%d",ibin,ibin2));       
	            h_chg_recopt_recoptsube0_d[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_recoptsube0_d_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_recoptsube0_d_%d_%d",ibin,ibin2));       
	            h_chg_recopt_genptsube0_d[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_genptsube0_d_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_genptsube0_d_%d_%d",ibin,ibin2));       
                h_chg_effgenpt_genptsube0_d[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_effgenpt_genptsube0_d_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_d_cent%d_trk%d",ibin,ibin2));
                //h_chg_effgenpt_genptsube0_trkvary_d[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC_vary->Get(Form("h_chg_effgenpt_genptsube0_d_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_d_trkvary_cent%d_trk%d",ibin,ibin2));

                //h_chg_reco_full_trkvary_d[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC_vary->Get(Form("h_chg_eff_genpt_d_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_d_trkvary_%d_%d",ibin,ibin2));       
	            h_chg_reco_full_d[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_d_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_d_%d_%d",ibin,ibin2));       
                h_chg_reco_full_d_tounfold[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_d_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_d_datasample_%d_%d",ibin,ibin2));       
                h_chg_eff_gen_full_d[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_eff_genpt_d_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_eff_genpt_d_%d_%d",ibin,ibin2));       

    	        h_chg_reco_full_sube0_d[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_sube0_d_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_sube0_d_%d_%d",ibin,ibin2));
	            h_chg_gen_full_sube0_d[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_genpt_sube0_d_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_genpt_sube0_d_%d_%d",ibin,ibin2));

//g+others
	            h_chg_recoptsube0_genptsube0_g_oth[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recoptsube0_genptsube0_g_oth_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_genpt_g_oth_cent%d_trk%d",ibin,ibin2));       
	            h_chg_recopt_recoptsube0_g_oth[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_recoptsube0_g_oth_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_recoptsube0_g_oth_%d_%d",ibin,ibin2));       
	            h_chg_recopt_genptsube0_g_oth[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_genptsube0_g_oth_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_genptsube0_g_oth_%d_%d",ibin,ibin2));       
                h_chg_effgenpt_genptsube0_g_oth[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_effgenpt_genptsube0_g_oth_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_g_oth_cent%d_trk%d",ibin,ibin2));
                //h_chg_effgenpt_genptsube0_trkvary_g_oth[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC_vary->Get(Form("h_chg_effgenpt_genptsube0_g_oth_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_g_oth_trkvary_cent%d_trk%d",ibin,ibin2));

                //h_chg_reco_full_trkvary_g_oth[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC_vary->Get(Form("h_chg_eff_genpt_g_oth_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_g_oth_trkvary_%d_%d",ibin,ibin2));       
	            h_chg_reco_full_g_oth[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_g_oth_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_g_oth_%d_%d",ibin,ibin2));       
                h_chg_reco_full_g_oth_tounfold[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_g_oth_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_g_oth_datasample_%d_%d",ibin,ibin2));       
                h_chg_eff_gen_full_g_oth[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_eff_genpt_g_oth_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_eff_genpt_g_oth_%d_%d",ibin,ibin2));       

    	        h_chg_reco_full_sube0_g_oth[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_sube0_g_oth_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_sube0_g_oth_%d_%d",ibin,ibin2));
	            h_chg_gen_full_sube0_g_oth[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_genpt_sube0_g_oth_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_genpt_sube0_g_oth_%d_%d",ibin,ibin2));

//inclusive
	            //h_chg_recoptsube0_genptsube0[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recoptsube0_genptsube0_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recoptsube0_genptsube0_cent%d_trk%d",ibin,ibin2));       

                // h_chg_recoptsube0_genptsube0_up[ibin][ibin2]->Scale(1./h_chg_recopt_genpt_up[ibin][ibin2]->Integral());
                // h_chg_recoptsube0_genptsube0_d[ibin][ibin2]->Scale(1./h_chg_recopt_genpt_d[ibin][ibin2]->Integral());
                // h_chg_recoptsube0_genptsube0_g_oth[ibin][ibin2]->Scale(1./h_chg_recopt_genpt_g_oth[ibin][ibin2]->Integral());
                h_chg_recoptsube0_genptsube0_up[ibin][ibin2]->Scale(up);
                h_chg_recoptsube0_genptsube0_d[ibin][ibin2]->Scale(down);
                h_chg_recoptsube0_genptsube0_g_oth[ibin][ibin2]->Scale(g_oth);
                h_chg_recoptsube0_genptsube0[ibin][ibin2] = (TH2F*)h_chg_recoptsube0_genptsube0_up[ibin][ibin2]->Clone(Form("h_chg_recoptsube0_genptsube0_prior_cent%d_trk%d",ibin,ibin2));
                h_chg_recoptsube0_genptsube0[ibin][ibin2]->Add(h_chg_recoptsube0_genptsube0_d[ibin][ibin2]);
                h_chg_recoptsube0_genptsube0[ibin][ibin2]->Add(h_chg_recoptsube0_genptsube0_g_oth[ibin][ibin2]);

	            //h_chg_recopt_recoptsube0[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_recoptsube0_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_recoptsube0_%d_%d",ibin,ibin2));       

                // h_chg_recopt_recoptsube0_up[ibin][ibin2]->Scale(1./h_chg_recopt_recoptsube0_up[ibin][ibin2]->Integral());
                // h_chg_recopt_recoptsube0_d[ibin][ibin2]->Scale(1./h_chg_recopt_recoptsube0_d[ibin][ibin2]->Integral());
                // h_chg_recopt_recoptsube0_g_oth[ibin][ibin2]->Scale(1./h_chg_recopt_recoptsube0_g_oth[ibin][ibin2]->Integral());
                h_chg_recopt_recoptsube0_up[ibin][ibin2]->Scale(up);
                h_chg_recopt_recoptsube0_d[ibin][ibin2]->Scale(down);
                h_chg_recopt_recoptsube0_g_oth[ibin][ibin2]->Scale(g_oth);
                h_chg_recopt_recoptsube0[ibin][ibin2] = (TH2F*)h_chg_recopt_recoptsube0_up[ibin][ibin2]->Clone(Form("h_chg_recopt_recoptsube0_prior_cent%d_trk%d",ibin,ibin2));
                h_chg_recopt_recoptsube0[ibin][ibin2]->Add(h_chg_recopt_recoptsube0_d[ibin][ibin2]);
                h_chg_recopt_recoptsube0[ibin][ibin2]->Add(h_chg_recopt_recoptsube0_g_oth[ibin][ibin2]);

	            h_chg_recopt_recoptsube0[ibin][ibin2]->RebinX(10);
	            h_chg_recopt_recoptsube0[ibin][ibin2]->RebinY(10);

	            //h_chg_recopt_genptsube0[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_genptsube0_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_genptsube0_%d_%d",ibin,ibin2));       
                //h_chg_recopt_genptsube0[ibin][ibin2]->Scale(1./h_chg_recopt_genptsube0[ibin][ibin2]->Integral());

                // h_chg_recopt_genptsube0_up[ibin][ibin2]->Scale(1./h_chg_recopt_genptsube0_up[ibin][ibin2]->Integral());
                // h_chg_recopt_genptsube0_d[ibin][ibin2]->Scale(1./h_chg_recopt_genptsube0_d[ibin][ibin2]->Integral());
                // h_chg_recopt_genptsube0_g_oth[ibin][ibin2]->Scale(1./h_chg_recopt_genptsube0_g_oth[ibin][ibin2]->Integral());
                h_chg_recopt_genptsube0_up[ibin][ibin2]->Scale(up);
                h_chg_recopt_genptsube0_d[ibin][ibin2]->Scale(down);
                h_chg_recopt_genptsube0_g_oth[ibin][ibin2]->Scale(g_oth);
                h_chg_recopt_genptsube0[ibin][ibin2] = (TH2F*)h_chg_recopt_genptsube0_up[ibin][ibin2]->Clone(Form("h_chg_recopt_genptsube0_prior_cent%d_trk%d",ibin,ibin2));
                h_chg_recopt_genptsube0[ibin][ibin2]->Add(h_chg_recopt_genptsube0_d[ibin][ibin2]);
                h_chg_recopt_genptsube0[ibin][ibin2]->Add(h_chg_recopt_genptsube0_g_oth[ibin][ibin2]);

                h_chg_recopt_genptsube0[ibin][ibin2]->RebinX(10);
                h_chg_recopt_genptsube0[ibin][ibin2]->RebinY(10);

                //h_chg_effgenpt_genptsube0[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_effgenpt_genptsube0_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_cent%d_trk%d",ibin,ibin2));

                // h_chg_effgenpt_genptsube0_up[ibin][ibin2]->Scale(1./h_chg_effgenpt_genptsube0_up[ibin][ibin2]->Integral());
                // h_chg_effgenpt_genptsube0_d[ibin][ibin2]->Scale(1./h_chg_effgenpt_genptsube0_d[ibin][ibin2]->Integral());
                // h_chg_effgenpt_genptsube0_g_oth[ibin][ibin2]->Scale(1./h_chg_effgenpt_genptsube0_g_oth[ibin][ibin2]->Integral());
                h_chg_effgenpt_genptsube0_up[ibin][ibin2]->Scale(up);
                h_chg_effgenpt_genptsube0_d[ibin][ibin2]->Scale(down);
                h_chg_effgenpt_genptsube0_g_oth[ibin][ibin2]->Scale(g_oth);
                h_chg_effgenpt_genptsube0[ibin][ibin2] = (TH2F*)h_chg_effgenpt_genptsube0_up[ibin][ibin2]->Clone(Form("h_chg_effgenpt_genptsube0_prior_cent%d_trk%d",ibin,ibin2));
                h_chg_effgenpt_genptsube0[ibin][ibin2]->Add(h_chg_effgenpt_genptsube0_d[ibin][ibin2]);
                h_chg_effgenpt_genptsube0[ibin][ibin2]->Add(h_chg_effgenpt_genptsube0_g_oth[ibin][ibin2]);

/*
                //h_chg_effgenpt_genptsube0_trkvary[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC_vary->Get(Form("h_chg_effgenpt_genptsube0_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_trkvary_cent%d_trk%d",ibin,ibin2));

                // h_chg_effgenpt_genptsube0_up[ibin][ibin2]->Scale(1./h_chg_effgenpt_genptsube0_up[ibin][ibin2]->Integral());
                // h_chg_effgenpt_genptsube0_d[ibin][ibin2]->Scale(1./h_chg_effgenpt_genptsube0_d[ibin][ibin2]->Integral());
                // h_chg_effgenpt_genptsube0_g_oth[ibin][ibin2]->Scale(1./h_chg_effgenpt_genptsube0_g_oth[ibin][ibin2]->Integral());
                h_chg_effgenpt_genptsube0_trkvary_up[ibin][ibin2]->Scale(up);
                h_chg_effgenpt_genptsube0_trkvary_d[ibin][ibin2]->Scale(down);
                h_chg_effgenpt_genptsube0_trkvary_g_oth[ibin][ibin2]->Scale(g_oth);
                h_chg_effgenpt_genptsube0_trkvary[ibin][ibin2] = (TH2F*)h_chg_effgenpt_genptsube0_trkvary_up[ibin][ibin2]->Clone(Form("h_chg_effgenpt_genptsube0_prior_trkvary_cent%d_trk%d",ibin,ibin2));
                h_chg_effgenpt_genptsube0_trkvary[ibin][ibin2]->Add(h_chg_effgenpt_genptsube0_trkvary_d[ibin][ibin2]);
                h_chg_effgenpt_genptsube0_trkvary[ibin][ibin2]->Add(h_chg_effgenpt_genptsube0_trkvary_g_oth[ibin][ibin2]);
*/
                h_chg_recorecopt_gengensube0pt[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recorecopt_gengensube0pt_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recorecopt_gengensube0pt_cent%d_trk%d",ibin,ibin2));

	            //h_chg_reco_full[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_%d_%d",ibin,ibin2));       
                //h_chg_reco_full[ibin][ibin2]->Scale(1./h_chg_reco_full[ibin][ibin2]->Integral());

                // h_chg_reco_full_up[ibin][ibin2]->Scale(1./h_chg_reco_full_up[ibin][ibin2]->Integral());
                // h_chg_reco_full_d[ibin][ibin2]->Scale(1./h_chg_reco_full_d[ibin][ibin2]->Integral());
                // h_chg_reco_full_g_oth[ibin][ibin2]->Scale(1./h_chg_reco_full_g_oth[ibin][ibin2]->Integral());
                //cout<<"before: "<<(h_chg_reco_full_up[ibin][ibin2]->Integral()+h_chg_reco_full_d[ibin][ibin2]->Integral()+h_chg_reco_full_g_oth[ibin][ibin2]->Integral())/h_chg_reco_full[ibin][ibin2]->Integral()<<endl;
                h_chg_reco_full_up[ibin][ibin2]->Scale(up);
                h_chg_reco_full_d[ibin][ibin2]->Scale(down);
                h_chg_reco_full_g_oth[ibin][ibin2]->Scale(g_oth);
                //cout<<"after: "<<(h_chg_reco_full_up[ibin][ibin2]->Integral()+h_chg_reco_full_d[ibin][ibin2]->Integral()+h_chg_reco_full_g_oth[ibin][ibin2]->Integral())/h_chg_reco_full[ibin][ibin2]->Integral()<<endl;
                h_chg_reco_full[ibin][ibin2] = (TH1D*)h_chg_reco_full_up[ibin][ibin2]->Clone(Form("h_chg_reco_full_cent%d_trk%d",ibin,ibin2));
                h_chg_reco_full[ibin][ibin2]->Add(h_chg_reco_full_d[ibin][ibin2]);
                h_chg_reco_full[ibin][ibin2]->Add(h_chg_reco_full_g_oth[ibin][ibin2]);

                //h_chg_reco_full_tounfold[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC_tounfold->Get(Form("h_chg_recopt_datasample_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_datasample_%d_%d",ibin,ibin2));       
                //h_chg_reco_full_tounfold[ibin][ibin2]->Scale(1./h_chg_reco_full_tounfold[ibin][ibin2]->Integral());

                // h_chg_reco_full_up_tounfold[ibin][ibin2]->Scale(1./h_chg_reco_full_up_tounfold[ibin][ibin2]->Integral());
                // h_chg_reco_full_d_tounfold[ibin][ibin2]->Scale(1./h_chg_reco_full_d[ibin][ibin2]->Integral());
                // h_chg_reco_full_g_oth_tounfold[ibin][ibin2]->Scale(1./h_chg_reco_full_g_oth_tounfold[ibin][ibin2]->Integral());
                //cout<<"before: "<<(h_chg_reco_full_up[ibin][ibin2]->Integral()+h_chg_reco_full_d[ibin][ibin2]->Integral()+h_chg_reco_full_g_oth[ibin][ibin2]->Integral())/h_chg_reco_full[ibin][ibin2]->Integral()<<endl;
                h_chg_reco_full_up_tounfold[ibin][ibin2]->Scale(0.763);
                h_chg_reco_full_d_tounfold[ibin][ibin2]->Scale(1.379);
                h_chg_reco_full_g_oth_tounfold[ibin][ibin2]->Scale(1.);
                //cout<<"after: "<<(h_chg_reco_full_up[ibin][ibin2]->Integral()+h_chg_reco_full_d[ibin][ibin2]->Integral()+h_chg_reco_full_g_oth[ibin][ibin2]->Integral())/h_chg_reco_full[ibin][ibin2]->Integral()<<endl;
                h_chg_reco_full_tounfold[ibin][ibin2] = (TH1D*)h_chg_reco_full_up_tounfold[ibin][ibin2]->Clone(Form("h_chg_reco_full_tounfold_cent%d_trk%d",ibin,ibin2));
                h_chg_reco_full_tounfold[ibin][ibin2]->Add(h_chg_reco_full_d_tounfold[ibin][ibin2]);
                h_chg_reco_full_tounfold[ibin][ibin2]->Add(h_chg_reco_full_g_oth_tounfold[ibin][ibin2]);

/*
                //h_chg_reco_full_trkvary[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC_vary->Get(Form("h_chg_recopt_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_trkvary_%d_%d",ibin,ibin2));       
                //h_chg_reco_full_trkvary[ibin][ibin2]->Scale(1./h_chg_reco_full_trkvary[ibin][ibin2]->Integral());

                // h_chg_reco_full_trkvary_up[ibin][ibin2]->Scale(1./h_chg_reco_full_trkvary_up[ibin][ibin2]->Integral());
                // h_chg_reco_full_trkvary_d[ibin][ibin2]->Scale(1./h_chg_reco_full_trkvary_d[ibin][ibin2]->Integral());
                // h_chg_reco_full_trkvary_g_oth[ibin][ibin2]->Scale(1./h_chg_reco_full_trkvary_g_oth[ibin][ibin2]->Integral());
                h_chg_reco_full_trkvary_up[ibin][ibin2]->Scale(up);
                h_chg_reco_full_trkvary_d[ibin][ibin2]->Scale(down);
                h_chg_reco_full_trkvary_g_oth[ibin][ibin2]->Scale(g_oth);
                h_chg_reco_full_trkvary[ibin][ibin2] = (TH1D*)h_chg_reco_full_trkvary_up[ibin][ibin2]->Clone(Form("h_chg_reco_full_trkvary_cent%d_trk%d",ibin,ibin2));
                h_chg_reco_full_trkvary[ibin][ibin2]->Add(h_chg_reco_full_trkvary_d[ibin][ibin2]);
                h_chg_reco_full_trkvary[ibin][ibin2]->Add(h_chg_reco_full_trkvary_g_oth[ibin][ibin2]);
*/
	            //h_chg_reco_full_sube0[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_sube0_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_sube0_%d_%d",ibin,ibin2));       
                //h_chg_reco_full_sube0[ibin][ibin2]->Scale(1./h_chg_reco_full_sube0[ibin][ibin2]->Integral());

                //h_chg_reco_full_sube0_up[ibin][ibin2]->Scale(1./h_chg_reco_full_sube0_up[ibin][ibin2]->Integral());
                //h_chg_reco_full_sube0_d[ibin][ibin2]->Scale(1./h_chg_reco_full_sube0_d[ibin][ibin2]->Integral());
                //h_chg_reco_full_sube0_g_oth[ibin][ibin2]->Scale(1./h_chg_reco_full_sube0_g_oth[ibin][ibin2]->Integral());
                h_chg_reco_full_sube0_up[ibin][ibin2]->Scale(up);
                h_chg_reco_full_sube0_d[ibin][ibin2]->Scale(down);
                h_chg_reco_full_sube0_g_oth[ibin][ibin2]->Scale(g_oth);
                h_chg_reco_full_sube0[ibin][ibin2] = (TH1D*)h_chg_reco_full_sube0_up[ibin][ibin2]->Clone(Form("h_chg_reco_full_sube0_cent%d_trk%d",ibin,ibin2));
                h_chg_reco_full_sube0[ibin][ibin2]->Add(h_chg_reco_full_sube0_d[ibin][ibin2]);
                h_chg_reco_full_sube0[ibin][ibin2]->Add(h_chg_reco_full_sube0_g_oth[ibin][ibin2]);
 
                h_chg_reco_full_sube0[ibin][ibin2]->Rebin(10);

	            //h_chg_gen_full_sube0[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_genpt_sube0_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_genpt_sube0_%d_%d",ibin,ibin2));
                //h_chg_gen_full_sube0[ibin][ibin2]->Scale(1./h_chg_gen_full_sube0[ibin][ibin2]->Integral());

                // h_chg_gen_full_sube0_up[ibin][ibin2]->Scale(1./h_chg_gen_full_sube0_up[ibin][ibin2]->Integral());
                // h_chg_gen_full_sube0_d[ibin][ibin2]->Scale(1./h_chg_gen_full_sube0_d[ibin][ibin2]->Integral());
                // h_chg_gen_full_sube0_g_oth[ibin][ibin2]->Scale(1./h_chg_gen_full_sube0_g_oth[ibin][ibin2]->Integral());
                h_chg_gen_full_sube0_up[ibin][ibin2]->Scale(up);
                h_chg_gen_full_sube0_d[ibin][ibin2]->Scale(down);
                h_chg_gen_full_sube0_g_oth[ibin][ibin2]->Scale(g_oth);
                h_chg_gen_full_sube0[ibin][ibin2] = (TH1D*)h_chg_gen_full_sube0_up[ibin][ibin2]->Clone(Form("h_chg_gen_full_sube0_cent%d_trk%d",ibin,ibin2));
                h_chg_gen_full_sube0[ibin][ibin2]->Add(h_chg_gen_full_sube0_d[ibin][ibin2]);
                h_chg_gen_full_sube0[ibin][ibin2]->Add(h_chg_gen_full_sube0_g_oth[ibin][ibin2]);

                //h_chg_eff_gen_full[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_eff_genpt_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_eff_genpt_%d_%d",ibin,ibin2));       

                // h_chg_eff_gen_full_up[ibin][ibin2]->Scale(1./h_chg_eff_gen_full_up[ibin][ibin2]->Integral());
                // h_chg_eff_gen_full_d[ibin][ibin2]->Scale(1./h_chg_eff_gen_full_d[ibin][ibin2]->Integral());
                // h_chg_eff_gen_full_g_oth[ibin][ibin2]->Scale(1./h_chg_eff_gen_full_g_oth[ibin][ibin2]->Integral());
                h_chg_eff_gen_full_up[ibin][ibin2]->Scale(up);
                h_chg_eff_gen_full_d[ibin][ibin2]->Scale(down);
                h_chg_eff_gen_full_g_oth[ibin][ibin2]->Scale(g_oth);
                h_chg_eff_gen_full[ibin][ibin2] = (TH1D*)h_chg_eff_gen_full_up[ibin][ibin2]->Clone(Form("h_chg_eff_gen_full_cent%d_trk%d",ibin,ibin2));
                h_chg_eff_gen_full[ibin][ibin2]->Add(h_chg_eff_gen_full_d[ibin][ibin2]);
                h_chg_eff_gen_full[ibin][ibin2]->Add(h_chg_eff_gen_full_g_oth[ibin][ibin2]);

                h_chg_eff_gen_full[ibin][ibin2]->Rebin(10);
                //cout<<"h_chg_eff_gen_full[ibin][ibin2]->"<<h_chg_eff_gen_full[ibin][ibin2]->Integral()<<endl;

                h_chg_gengensube0_full[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_gengensube0pt_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_gengensube0pt_%d_%d",ibin,ibin2));       

//////////kappa
//up
                h_chg_recoptsube0_genptsube0_up_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recoptsube0_genptsube0_up_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recoptsube0_genptsube0_up_kappa_cent%d_trk%d",ibin,ibin2));       
                h_chg_recopt_recoptsube0_up_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_recoptsube0_up_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_recoptsube0_up_kappa_%d_%d",ibin,ibin2));       
                h_chg_recopt_genptsube0_up_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_genptsube0_up_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_genptsube0_up_kappa_%d_%d",ibin,ibin2));       
                h_chg_effgenpt_genptsube0_up_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_effgenpt_genptsube0_up_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_up_kappa_cent%d_trk%d",ibin,ibin2));
    
                h_chg_reco_full_up_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_up_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_up_kappa_%d_%d",ibin,ibin2));       
                h_chg_eff_gen_full_up_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_eff_genpt_up_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_eff_genpt_up_kappa_%d_%d",ibin,ibin2));       

                h_chg_reco_full_sube0_up_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_sube0_up_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_sube0_up_kappa_%d_%d",ibin,ibin2));
                h_chg_gen_full_sube0_up_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_genpt_sube0_up_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_genpt_sube0_up_kappa_%d_%d",ibin,ibin2));

//down
                h_chg_recoptsube0_genptsube0_d_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recoptsube0_genptsube0_d_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recoptsube0_genptsube0_d_kappa_cent%d_trk%d",ibin,ibin2));       
                h_chg_recopt_recoptsube0_d_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_recoptsube0_d_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_recoptsube0_d_kappa_%d_%d",ibin,ibin2));       
                h_chg_recopt_genptsube0_d_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_genptsube0_d_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_genptsube0_d_kappa_%d_%d",ibin,ibin2));       
                h_chg_effgenpt_genptsube0_d_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_effgenpt_genptsube0_d_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_d_kappa_cent%d_trk%d",ibin,ibin2));

                h_chg_reco_full_d_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_d_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_d_kappa_%d_%d",ibin,ibin2));       
                h_chg_eff_gen_full_d_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_eff_genpt_d_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_eff_genpt_d_kappa_%d_%d",ibin,ibin2));       

                h_chg_reco_full_sube0_d_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_sube0_d_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_sube0_d_kappa_%d_%d",ibin,ibin2));
                h_chg_gen_full_sube0_d_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_genpt_sube0_d_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_genpt_sube0_d_kappa_%d_%d",ibin,ibin2));

//g+others
                h_chg_recoptsube0_genptsube0_g_oth_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recoptsube0_genptsube0_g_oth_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recoptsube0_genptsube0_g_oth_kappa_cent%d_trk%d",ibin,ibin2));       
                h_chg_recopt_recoptsube0_g_oth_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_recoptsube0_g_oth_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_recoptsube0_g_oth_kappa_%d_%d",ibin,ibin2));       
                h_chg_recopt_genptsube0_g_oth_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_genptsube0_g_oth_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_genptsube0_g_oth_kappa_%d_%d",ibin,ibin2));       
                h_chg_effgenpt_genptsube0_g_oth_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_effgenpt_genptsube0_g_oth_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_g_oth_kappa_cent%d_trk%d",ibin,ibin2));

                h_chg_reco_full_g_oth_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_g_oth_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_g_oth_kappa_%d_%d",ibin,ibin2));       
                h_chg_eff_gen_full_g_oth_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_eff_genpt_g_oth_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_eff_genpt_g_oth_kappa_%d_%d",ibin,ibin2));       

                h_chg_reco_full_sube0_g_oth_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_sube0_g_oth_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_sube0_g_oth_kappa_%d_%d",ibin,ibin2));
                h_chg_gen_full_sube0_g_oth_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_genpt_sube0_g_oth_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_genpt_sube0_g_oth_kappa_%d_%d",ibin,ibin2));

//inclusive
                //h_chg_recoptsube0_genptsube0_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recoptsube0_genptsube0_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recoptsube0_genptsube0_kappa_cent%d_trk%d",ibin,ibin2));       

                // h_chg_recoptsube0_genptsube0_up_kappa[ibin][ibin2]->Scale(1./h_chg_recoptsube0_genptsube0_up_kappa[ibin][ibin2]->Integral());
                // h_chg_recoptsube0_genptsube0_d_kappa[ibin][ibin2]->Scale(1./h_chg_recoptsube0_genptsube0_d_kappa[ibin][ibin2]->Integral());
                // h_chg_recoptsube0_genptsube0_g_oth_kappa[ibin][ibin2]->Scale(1./h_chg_recoptsube0_genptsube0_g_oth_kappa[ibin][ibin2]->Integral());
                h_chg_recoptsube0_genptsube0_up_kappa[ibin][ibin2]->Scale(up);
                h_chg_recoptsube0_genptsube0_d_kappa[ibin][ibin2]->Scale(down);
                h_chg_recoptsube0_genptsube0_g_oth_kappa[ibin][ibin2]->Scale(g_oth);
                h_chg_recoptsube0_genptsube0_kappa[ibin][ibin2] = (TH2F*)h_chg_recoptsube0_genptsube0_up_kappa[ibin][ibin2]->Clone(Form("h_chg_recoptsube0_genptsube0_prior_kappa_cent%d_trk%d",ibin,ibin2));
                h_chg_recoptsube0_genptsube0_kappa[ibin][ibin2]->Add(h_chg_recoptsube0_genptsube0_d_kappa[ibin][ibin2]);
                h_chg_recoptsube0_genptsube0_kappa[ibin][ibin2]->Add(h_chg_recoptsube0_genptsube0_g_oth_kappa[ibin][ibin2]);

                //h_chg_recopt_recoptsube0_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_recoptsube0_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_recoptsube0_kappa_%d_%d",ibin,ibin2));       

                // h_chg_recopt_recoptsube0_up_kappa[ibin][ibin2]->Scale(1./h_chg_recopt_recoptsube0_up_kappa[ibin][ibin2]->Integral());
                // h_chg_recopt_recoptsube0_d_kappa[ibin][ibin2]->Scale(1./h_chg_recopt_recoptsube0_d_kappa[ibin][ibin2]->Integral());
                // h_chg_recopt_recoptsube0_g_oth_kappa[ibin][ibin2]->Scale(1./h_chg_recopt_recoptsube0_g_oth_kappa[ibin][ibin2]->Integral());
                h_chg_recopt_recoptsube0_up_kappa[ibin][ibin2]->Scale(up);
                h_chg_recopt_recoptsube0_d_kappa[ibin][ibin2]->Scale(down);
                h_chg_recopt_recoptsube0_g_oth_kappa[ibin][ibin2]->Scale(g_oth);
                h_chg_recopt_recoptsube0_kappa[ibin][ibin2] = (TH2F*)h_chg_recopt_recoptsube0_up_kappa[ibin][ibin2]->Clone(Form("h_chg_recopt_recoptsube0_prior_kappa_cent%d_trk%d",ibin,ibin2));
                h_chg_recopt_recoptsube0_kappa[ibin][ibin2]->Add(h_chg_recopt_recoptsube0_d_kappa[ibin][ibin2]);
                h_chg_recopt_recoptsube0_kappa[ibin][ibin2]->Add(h_chg_recopt_recoptsube0_g_oth_kappa[ibin][ibin2]);

                h_chg_recopt_recoptsube0_kappa[ibin][ibin2]->RebinX(10);
                h_chg_recopt_recoptsube0_kappa[ibin][ibin2]->RebinY(10);

                //h_chg_recopt_genptsube0_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_genptsube0_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_genptsube0_kappa_%d_%d",ibin,ibin2));       
                //h_chg_recopt_genptsube0_kappa[ibin][ibin2]->Scale(1./h_chg_recopt_genptsube0_kappa[ibin][ibin2]->Integral());

                // h_chg_recopt_genptsube0_up_kappa[ibin][ibin2]->Scale(1./h_chg_recopt_genptsube0_up_kappa[ibin][ibin2]->Integral());
                // h_chg_recopt_genptsube0_d_kappa[ibin][ibin2]->Scale(1./h_chg_recopt_genptsube0_d_kappa[ibin][ibin2]->Integral());
                // h_chg_recopt_genptsube0_g_oth_kappa[ibin][ibin2]->Scale(1./h_chg_recopt_genptsube0_g_oth_kappa[ibin][ibin2]->Integral());
                h_chg_recopt_genptsube0_up_kappa[ibin][ibin2]->Scale(up);
                h_chg_recopt_genptsube0_d_kappa[ibin][ibin2]->Scale(down);
                h_chg_recopt_genptsube0_g_oth_kappa[ibin][ibin2]->Scale(g_oth);
                h_chg_recopt_genptsube0_kappa[ibin][ibin2] = (TH2F*)h_chg_recopt_genptsube0_up_kappa[ibin][ibin2]->Clone(Form("h_chg_recopt_genptsube0_prior_kappa_cent%d_trk%d",ibin,ibin2));
                h_chg_recopt_genptsube0_kappa[ibin][ibin2]->Add(h_chg_recopt_genptsube0_d_kappa[ibin][ibin2]);
                h_chg_recopt_genptsube0_kappa[ibin][ibin2]->Add(h_chg_recopt_genptsube0_g_oth_kappa[ibin][ibin2]);

                h_chg_recopt_genptsube0_kappa[ibin][ibin2]->RebinX(10);
                h_chg_recopt_genptsube0_kappa[ibin][ibin2]->RebinY(10);

                //h_chg_effgenpt_genptsube0_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_effgenpt_genptsube0_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_effgenpt_genptsube0_kappa_cent%d_trk%d",ibin,ibin2));

                // h_chg_effgenpt_genptsube0_up_kappa[ibin][ibin2]->Scale(1./h_chg_effgenpt_genptsube0_up_kappa[ibin][ibin2]->Integral());
                // h_chg_effgenpt_genptsube0_d_kappa[ibin][ibin2]->Scale(1./h_chg_effgenpt_genptsube0_d_kappa[ibin][ibin2]->Integral());
                // h_chg_effgenpt_genptsube0_g_oth_kappa[ibin][ibin2]->Scale(1./h_chg_effgenpt_genptsube0_g_oth_kappa[ibin][ibin2]->Integral());
                h_chg_effgenpt_genptsube0_up_kappa[ibin][ibin2]->Scale(up);
                h_chg_effgenpt_genptsube0_d_kappa[ibin][ibin2]->Scale(down);
                h_chg_effgenpt_genptsube0_g_oth_kappa[ibin][ibin2]->Scale(g_oth);
                h_chg_effgenpt_genptsube0_kappa[ibin][ibin2] = (TH2F*)h_chg_effgenpt_genptsube0_up_kappa[ibin][ibin2]->Clone(Form("h_chg_effgenpt_genptsube0_prior_kappa_cent%d_trk%d",ibin,ibin2));
                h_chg_effgenpt_genptsube0_kappa[ibin][ibin2]->Add(h_chg_effgenpt_genptsube0_d_kappa[ibin][ibin2]);
                h_chg_effgenpt_genptsube0_kappa[ibin][ibin2]->Add(h_chg_effgenpt_genptsube0_g_oth_kappa[ibin][ibin2]);


                h_chg_recorecopt_gengensube0pt_kappa[ibin][ibin2] = (TH2F*)jetcharge_PbPb_MC->Get(Form("h_chg_recorecopt_gengensube0pt_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recorecopt_gengensube0pt_kappa_cent%d_trk%d",ibin,ibin2));

                //h_chg_reco_full_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_kappa_%d_%d",ibin,ibin2));       
                //h_chg_reco_full_kappa[ibin][ibin2]->Scale(1./h_chg_reco_full_kappa[ibin][ibin2]->Integral());

                // h_chg_reco_full_up_kappa[ibin][ibin2]->Scale(1./h_chg_reco_full_up_kappa[ibin][ibin2]->Integral());
                // h_chg_reco_full_d_kappa[ibin][ibin2]->Scale(1./h_chg_reco_full_d_kappa[ibin][ibin2]->Integral());
                // h_chg_reco_full_g_oth_kappa[ibin][ibin2]->Scale(1./h_chg_reco_full_g_oth_kappa[ibin][ibin2]->Integral());
                h_chg_reco_full_up_kappa[ibin][ibin2]->Scale(up);
                h_chg_reco_full_d_kappa[ibin][ibin2]->Scale(down);
                h_chg_reco_full_g_oth_kappa[ibin][ibin2]->Scale(g_oth);
                h_chg_reco_full_kappa[ibin][ibin2] = (TH1D*)h_chg_reco_full_up_kappa[ibin][ibin2]->Clone(Form("h_chg_reco_full_kappa_cent%d_trk%d",ibin,ibin2));
                h_chg_reco_full_kappa[ibin][ibin2]->Add(h_chg_reco_full_d_kappa[ibin][ibin2]);
                h_chg_reco_full_kappa[ibin][ibin2]->Add(h_chg_reco_full_g_oth_kappa[ibin][ibin2]);

                h_chg_reco_full_sube0_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_recopt_sube0_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_recopt_sube0_kappa_%d_%d",ibin,ibin2));       
                //h_chg_reco_full_sube0_kappa[ibin][ibin2]->Scale(1./h_chg_reco_full_sube0_kappa[ibin][ibin2]->Integral());
/*                
                h_chg_reco_full_sube0_up_kappa[ibin][ibin2]->Scale(1./h_chg_reco_full_sube0_up_kappa[ibin][ibin2]->Integral());
                h_chg_reco_full_sube0_d_kappa[ibin][ibin2]->Scale(1./h_chg_reco_full_sube0_d_kappa[ibin][ibin2]->Integral());
                h_chg_reco_full_sube0_g_oth_kappa[ibin][ibin2]->Scale(1./h_chg_reco_full_sube0_g_oth_kappa[ibin][ibin2]->Integral());
                h_chg_reco_full_sube0_up_kappa[ibin][ibin2]->Scale(up);
                h_chg_reco_full_sube0_d_kappa[ibin][ibin2]->Scale(down);
                h_chg_reco_full_sube0_g_oth_kappa[ibin][ibin2]->Scale(g_oth);
                h_chg_reco_full_sube0_kappa[ibin][ibin2] = (TH1D*)h_chg_reco_full_sube0_up_kappa[ibin][ibin2]->Clone(Form("h_chg_reco_full_sube0_kappa_cent%d_trk%d",ibin,ibin2));
                h_chg_reco_full_sube0_kappa[ibin][ibin2]->Add(h_chg_reco_full_sube0_d_kappa[ibin][ibin2]);
                h_chg_reco_full_sube0_kappa[ibin][ibin2]->Add(h_chg_reco_full_sube0_g_oth_kappa[ibin][ibin2]);
*/ 
                //h_chg_gen_full_sube0_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_genpt_sube0_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_genpt_sube0_kappa_%d_%d",ibin,ibin2));
                //h_chg_gen_full_sube0_kappa[ibin][ibin2]->Scale(1./h_chg_gen_full_sube0_kappa[ibin][ibin2]->Integral());

                // h_chg_gen_full_sube0_up_kappa[ibin][ibin2]->Scale(1./h_chg_gen_full_sube0_up_kappa[ibin][ibin2]->Integral());
                // h_chg_gen_full_sube0_d_kappa[ibin][ibin2]->Scale(1./h_chg_gen_full_sube0_d_kappa[ibin][ibin2]->Integral());
                // h_chg_gen_full_sube0_g_oth_kappa[ibin][ibin2]->Scale(1./h_chg_gen_full_sube0_g_oth_kappa[ibin][ibin2]->Integral());
                h_chg_gen_full_sube0_up_kappa[ibin][ibin2]->Scale(up);
                h_chg_gen_full_sube0_d_kappa[ibin][ibin2]->Scale(down);
                h_chg_gen_full_sube0_g_oth_kappa[ibin][ibin2]->Scale(g_oth);
                h_chg_gen_full_sube0_kappa[ibin][ibin2] = (TH1D*)h_chg_gen_full_sube0_up_kappa[ibin][ibin2]->Clone(Form("h_chg_gen_full_sube0_kappa_cent%d_trk%d",ibin,ibin2));
                h_chg_gen_full_sube0_kappa[ibin][ibin2]->Add(h_chg_gen_full_sube0_d_kappa[ibin][ibin2]);
                h_chg_gen_full_sube0_kappa[ibin][ibin2]->Add(h_chg_gen_full_sube0_g_oth_kappa[ibin][ibin2]);

                //h_chg_eff_gen_full_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_eff_genpt_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_eff_genpt_kappa_%d_%d",ibin,ibin2));       
                
                // h_chg_eff_gen_full_up_kappa[ibin][ibin2]->Scale(1./h_chg_eff_gen_full_up_kappa[ibin][ibin2]->Integral());
                // h_chg_eff_gen_full_d_kappa[ibin][ibin2]->Scale(1./h_chg_eff_gen_full_d_kappa[ibin][ibin2]->Integral());
                // h_chg_eff_gen_full_g_oth_kappa[ibin][ibin2]->Scale(1./h_chg_eff_gen_full_g_oth_kappa[ibin][ibin2]->Integral());
                h_chg_eff_gen_full_up_kappa[ibin][ibin2]->Scale(up);
                h_chg_eff_gen_full_d_kappa[ibin][ibin2]->Scale(down);
                h_chg_eff_gen_full_g_oth_kappa[ibin][ibin2]->Scale(g_oth);
                h_chg_eff_gen_full_kappa[ibin][ibin2] = (TH1D*)h_chg_eff_gen_full_up_kappa[ibin][ibin2]->Clone(Form("h_chg_eff_gen_full_kappa_cent%d_trk%d",ibin,ibin2));
                h_chg_eff_gen_full_kappa[ibin][ibin2]->Add(h_chg_eff_gen_full_d_kappa[ibin][ibin2]);
                h_chg_eff_gen_full_kappa[ibin][ibin2]->Add(h_chg_eff_gen_full_g_oth_kappa[ibin][ibin2]);
                h_chg_eff_gen_full_kappa[ibin][ibin2]->Rebin(10);

                h_chg_gengensube0_full_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_MC->Get(Form("h_chg_gengensube0pt_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_gengensube0pt_kappa_%d_%d",ibin,ibin2));       

//data
                h_chg_data_kappa[ibin][ibin2] = (TH1D*)jetcharge_PbPb_data->Get(Form("h_chg_recopt_kappa_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_data_kappa_cent%d_trk%d",ibin,ibin2)); 
	            h_chg_data[ibin][ibin2] = (TH1D*)jetcharge_PbPb_data->Get(Form("h_chg_recopt_cent%d_trk%d",ibin-1,ibin2))->Clone(Form("h_chg_data_cent%d_trk%d",ibin,ibin2)); 

	        }	    	
	    }
			
	    for (int ibin2=0; ibin2<ntrkbins; ibin2++){
            h_chg_data[ibin][ibin2]->Rebin(10);  
            h_chg_data_kappa[ibin][ibin2]->Rebin(10);  

            //h_chg_reco_full_trkvary[ibin][ibin2]->Rebin(10);

	        h_chg_recoptsube0_genptsube0[ibin][ibin2]->RebinX(10);
	        h_chg_recoptsube0_genptsube0[ibin][ibin2]->RebinY(10);

            h_chg_reco_full[ibin][ibin2]->Rebin(10);
            h_chg_reco_full_tounfold[ibin][ibin2]->Rebin(10);
            h_chg_gen_full_sube0[ibin][ibin2]->Rebin(10); 

            h_chg_effgenpt_genptsube0[ibin][ibin2]->RebinX(10);
            h_chg_effgenpt_genptsube0[ibin][ibin2]->RebinY(10);

            //h_chg_effgenpt_genptsube0_trkvary[ibin][ibin2]->RebinX(10);
            //h_chg_effgenpt_genptsube0_trkvary[ibin][ibin2]->RebinY(10);

            h_chg_recoptsube0_genptsube0_kappa[ibin][ibin2]->RebinX(10);
            h_chg_recoptsube0_genptsube0_kappa[ibin][ibin2]->RebinY(10);

            h_chg_reco_full_kappa[ibin][ibin2]->Rebin(10);

            //h_chg_reco_full_sube0_kappa[ibin][ibin2]->Rebin(10);
            h_chg_gen_full_sube0_kappa[ibin][ibin2]->Rebin(10); 

            h_chg_effgenpt_genptsube0_kappa[ibin][ibin2]->RebinX(10);
            h_chg_effgenpt_genptsube0_kappa[ibin][ibin2]->RebinY(10);
/*
            if(ibin>0){
                //h_chg_reco_full[ibin][ibin2]->Scale(100./h_chg_reco_full[ibin][ibin2]->Integral());
                //h_chg_reco_full_sube0[ibin][ibin2]->Scale(100./h_chg_reco_full_sube0[ibin][ibin2]->Integral());

		        for (int ibin3=0; ibin3<h_chg_recopt_recoptsube0[ibin][ibin2]->GetNbinsX(); ibin3++){
		            h_chg_reco_forbkg_prior[ibin][ibin2][ibin3] = h_chg_recopt_recoptsube0[ibin][ibin2]->ProjectionY(Form("h_chg_reco_forbkg_%d_%d_%d",ibin,ibin2,ibin3),ibin3+1,ibin3+1,"");
    		        for (int ibin4=0; ibin4<h_chg_reco_full[ibin][ibin2]->GetBinContent(ibin3+1); ibin4++){       	
	                    if(h_chg_reco_forbkg_prior[ibin][ibin2][ibin3]->Integral() > 0.) h_chg_reco_unfolded_forbkg_rand[ibin][ibin2]->Fill(h_chg_reco_forbkg_prior[ibin][ibin2][ibin3]->GetRandom());
    		        }
    		        for (int ibin4=0; ibin4<h_chg_data[ibin][ibin2]->GetBinContent(ibin3+1); ibin4++){       	
	                    if(h_chg_reco_forbkg_prior[ibin][ibin2][ibin3]->Integral() > 0.) h_chg_data_unfolded_forbkg_rand[ibin][ibin2]->Fill(h_chg_reco_forbkg_prior[ibin][ibin2][ibin3]->GetRandom());
    		        }
		        }
            }

            //h_chg_reco_full_sube0[ibin][ibin2]->Scale(100./h_chg_reco_full_sube0[ibin][ibin2]->Integral());
            //h_chg_gen_full_sube0[ibin][ibin2]->Scale(100./h_chg_gen_full_sube0[ibin][ibin2]->Integral());

	        for (int ibin3=0; ibin3<h_chg_recopt_genpt[ibin][ibin2]->GetNbinsX(); ibin3++){
	            h_chg_reco_prior[ibin][ibin2][ibin3] = h_chg_recopt_genpt[ibin][ibin2]->ProjectionY(Form("h_chg_reco_prior_%d_%d_%d",ibin,ibin2,ibin3),ibin3+1,ibin3+1,"");

		        if(ibin==0){
			        for (int ibin4=0; ibin4<h_chg_reco_full_sube0[ibin][ibin2]->GetBinContent(ibin3+1); ibin4++){
	                    if(h_chg_reco_prior[ibin][ibin2][ibin3]->Integral() > 0.) h_chg_reco_unfolded_rand[ibin][ibin2]->Fill(h_chg_reco_prior[ibin][ibin2][ibin3]->GetRandom());     		        	
			        }
			        for (int ibin4=0; ibin4<h_chg_data[ibin][ibin2]->GetBinContent(ibin3+1); ibin4++){
	                    if(h_chg_reco_prior[ibin][ibin2][ibin3]->Integral() > 0.) h_chg_data_unfolded_rand[ibin][ibin2]->Fill(h_chg_reco_prior[ibin][ibin2][ibin3]->GetRandom());     		        	
			        }
			    }
                else{
			        for (int ibin4=0; ibin4<h_chg_reco_unfolded_forbkg_rand[ibin][ibin2]->GetBinContent(ibin3+1); ibin4++){
	                    if(h_chg_reco_prior[ibin][ibin2][ibin3]->Integral() > 0.) h_chg_reco_unfolded_rand[ibin][ibin2]->Fill(h_chg_reco_prior[ibin][ibin2][ibin3]->GetRandom());     		        	
			        }
			        for (int ibin4=0; ibin4<h_chg_data_unfolded_forbkg_rand[ibin][ibin2]->GetBinContent(ibin3+1); ibin4++){
	                    if(h_chg_reco_prior[ibin][ibin2][ibin3]->Integral() > 0.) h_chg_data_unfolded_rand[ibin][ibin2]->Fill(h_chg_reco_prior[ibin][ibin2][ibin3]->GetRandom());     		        	
			        }
			    }
	        }
*/
	    }
    }
/*
    TCanvas *c_check = new TCanvas("c_check","c_check",1200,400);
    c_check->Divide(4,1);
    
    for(int ibin=1; ibin<nCBins; ibin++){
        cout<<ibin<<endl;
        const int ibin2=0;
        if(ibin==0) c_check->cd(1);
        else c_check->cd(5-ibin);
        h_chg_reco_full[ibin][0]->GetYaxis()->SetTitle("reco / eff gen");
        h_chg_reco_full[ibin][0]->GetXaxis()->SetRangeUser(-1.5,1.5);
        h_chg_reco_full[ibin][0]->SetLineColor(kRed-2);
        h_chg_reco_full[ibin][0]->Divide(h_chg_eff_gen_full[ibin][0]);
        //h_chg_reco_full[ibin][0]->Rebin(10);
        //h_chg_reco_full[ibin][0]->Scale(1./10.);
        h_chg_reco_full[ibin][0]->Draw("e1");    
    
        TString tmp1 = cent_tag_pp_MC[ibin];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.35,0.85, tmp1);
    }
    //h_chg_eff_gen_full[3][0]->Draw("e1 same");    
*/
    TLatex *tx1;
    TLatex *tx;

    ///RooUnfold

    for(int ibin=0; ibin<nCBins; ibin++){
       for(int ibin2=0; ibin2<ntrkbins; ibin2++){
		    if(ibin==0){
//unfold for tracking
//reco to gen
//ptcut/
/*
		        RooUnfoldResponse h_jetcharge_response (h_chg_reco_full[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_recoptsube0_genptsube0[ibin][ibin2],Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2));
		        RooUnfoldBayes unfold (&h_jetcharge_response, h_chg_reco_full_tounfold[ibin][ibin2], 1);
                //RooUnfoldSvd unfold (&h_jetcharge_response, h_chg_reco_full_sube0[ibin][ibin2], 1);
		        h_chg_reco_unfolded[ibin][ibin2] = (TH1D*)unfold.Hreco();

                RooUnfoldResponse h_jetcharge_response_data (h_chg_reco_full[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_recoptsube0_genptsube0[ibin][ibin2],Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold_data (&h_jetcharge_response_data, h_chg_data[ibin][ibin2], 2);
                //RooUnfoldSvd unfold_data (&h_jetcharge_response_data, h_chg_data[ibin][ibin2], 3);
                h_chg_data_trkunfold[ibin][ibin2] = (TH1D*)unfold_data.Hreco();             

//kappa
                RooUnfoldResponse h_jetcharge_response_kappa (h_chg_reco_full_kappa[ibin][ibin2], h_chg_gen_full_sube0_kappa[ibin][ibin2], h_chg_recoptsube0_genptsube0_kappa[ibin][ibin2],Form("h_jetcharge_response_kappa_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_kappa_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold_kappa (&h_jetcharge_response_kappa, h_chg_reco_full_kappa[ibin][ibin2], 1);
                //RooUnfoldSvd unfold (&h_jetcharge_response_kappa, h_chg_reco_full_sube0_kappa[ibin][ibin2], 1);
                h_chg_reco_unfolded_kappa[ibin][ibin2] = (TH1D*)unfold_kappa.Hreco();

                RooUnfoldResponse h_jetcharge_response_data_kappa (h_chg_reco_full_kappa[ibin][ibin2], h_chg_gen_full_sube0_kappa[ibin][ibin2], h_chg_recoptsube0_genptsube0_kappa[ibin][ibin2],Form("h_jetcharge_response_data_kappa_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_data_kappa_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold_data_kappa (&h_jetcharge_response_data_kappa, h_chg_data_kappa[ibin][ibin2], 2);
                //RooUnfoldSvd unfold_data_kappa (&h_jetcharge_response_data_kappa, h_chg_data_kappa[ibin][ibin2], 3);
                h_chg_data_trkunfold_kappa[ibin][ibin2] = (TH1D*)unfold_data_kappa.Hreco();
*/

//eff gen to gen
//ptcut

                RooUnfoldResponse h_jetcharge_response (h_chg_eff_gen_full_sube0[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_effgenpt_genptsube0[ibin][ibin2],Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2));
                //RooUnfoldResponse h_jetcharge_response (h_chg_reco_full_trkvary[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_effgenpt_genptsube0_trkvary[ibin][ibin2],Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold (&h_jetcharge_response, h_chg_eff_gen_full_sube0[ibin][ibin2], 1);
                h_chg_reco_unfolded[ibin][ibin2] = (TH1D*)unfold.Hreco();

                RooUnfoldResponse h_jetcharge_response_data (h_chg_eff_gen_full_sube0[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_effgenpt_genptsube0[ibin][ibin2],Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2));
                //RooUnfoldResponse h_jetcharge_response_data (h_chg_reco_full_trkvary[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_effgenpt_genptsube0_trkvary[ibin][ibin2],Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold_data (&h_jetcharge_response_data, h_chg_data[ibin][ibin2], 2);
                h_chg_data_trkunfold[ibin][ibin2] = (TH1D*)unfold_data.Hreco();             

//kappa
                RooUnfoldResponse h_jetcharge_response_kappa (h_chg_eff_gen_full_sube0_kappa[ibin][ibin2], h_chg_gen_full_sube0_kappa[ibin][ibin2], h_chg_effgenpt_genptsube0_kappa[ibin][ibin2],Form("h_jetcharge_response_kappa_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_kappa_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold_kappa (&h_jetcharge_response_kappa, h_chg_eff_gen_full_sube0_kappa[ibin][ibin2], 1);
                h_chg_reco_unfolded_kappa[ibin][ibin2] = (TH1D*)unfold_kappa.Hreco();

                RooUnfoldResponse h_jetcharge_response_data_kappa (h_chg_eff_gen_full_sube0_kappa[ibin][ibin2], h_chg_gen_full_sube0_kappa[ibin][ibin2], h_chg_effgenpt_genptsube0_kappa[ibin][ibin2],Form("h_jetcharge_response_data_kappa_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_data_kappa_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold_data_kappa (&h_jetcharge_response_data_kappa, h_chg_data_kappa[ibin][ibin2], 2);
                h_chg_data_trkunfold_kappa[ibin][ibin2] = (TH1D*)unfold_data_kappa.Hreco();

/*
//unfold for tracking and jet reco
                RooUnfoldResponse h_jetcharge_jetreco_response (h_chg_reco_full_sube0[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_recopt_genpt[ibin][ibin2],Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold (&h_jetcharge_response, h_chg_reco_full_sube0[ibin][ibin2], 1);
                h_chg_reco_unfolded[ibin][ibin2] = (TH1D*)unfold.Hreco();
*/
		    }
		    else{
/////////////////////////////
////////one step/////////////
/////////////////////////////

//reco to gen sube0
/*
//ptcut
		        RooUnfoldResponse h_jetcharge_response (h_chg_reco_full[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_recopt_genptsube0[ibin][ibin2],Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2));
		        RooUnfoldBayes unfold (&h_jetcharge_response, h_chg_reco_full_tounfold[ibin][ibin2], 1);
                //RooUnfoldSvd unfold (&h_jetcharge_response, h_chg_reco_full[ibin][ibin2], 1);
		        h_chg_reco_unfolded[ibin][ibin2] = (TH1D*)unfold.Hreco();

		        RooUnfoldResponse h_jetcharge_response_data (h_chg_reco_full[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_recopt_genptsube0[ibin][ibin2],Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2));
		        RooUnfoldBayes unfold_data (&h_jetcharge_response_data, h_chg_data[ibin][ibin2], 2);
                //RooUnfoldSvd unfold_data (&h_jetcharge_response_data, h_chg_data[ibin][ibin2], 3);
		        h_chg_data_trkunfold[ibin][ibin2] = (TH1D*)unfold_data.Hreco();		    	

//kappa
                RooUnfoldResponse h_jetcharge_response_kappa (h_chg_reco_full_kappa[ibin][ibin2], h_chg_gen_full_sube0_kappa[ibin][ibin2], h_chg_recopt_genptsube0_kappa[ibin][ibin2],Form("h_jetcharge_response_kappa_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_kappa_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold_kappa (&h_jetcharge_response_kappa, h_chg_reco_full_kappa[ibin][ibin2], 1);
                //RooUnfoldSvd unfold_kappa (&h_jetcharge_response_kappa, h_chg_reco_full_kappa[ibin][ibin2], 1);
                h_chg_reco_unfolded_kappa[ibin][ibin2] = (TH1D*)unfold_kappa.Hreco();

                RooUnfoldResponse h_jetcharge_response_data_kappa (h_chg_reco_full_kappa[ibin][ibin2], h_chg_gen_full_sube0_kappa[ibin][ibin2], h_chg_recopt_genptsube0_kappa[ibin][ibin2],Form("h_jetcharge_response_data_kappa_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_data_kappa_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold_data_kappa (&h_jetcharge_response_data_kappa, h_chg_data_kappa[ibin][ibin2], 2);
                //RooUnfoldSvd unfold_data_kappa (&h_jetcharge_response_data_kappa, h_chg_data_kappa[ibin][ibin2], 3);
                h_chg_data_trkunfold_kappa[ibin][ibin2] = (TH1D*)unfold_data_kappa.Hreco();
*/
//eff gen to gen sube0

//ptcut
                RooUnfoldResponse h_jetcharge_response (h_chg_eff_gen_full[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_effgenpt_genptsube0[ibin][ibin2],Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2));
                //RooUnfoldResponse h_jetcharge_response (h_chg_reco_full_trkvary[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_effgenpt_genptsube0_trkvary[ibin][ibin2],Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold (&h_jetcharge_response, h_chg_eff_gen_full[ibin][ibin2], 1);
                h_chg_reco_unfolded[ibin][ibin2] = (TH1D*)unfold.Hreco();

                RooUnfoldResponse h_jetcharge_response_data (h_chg_eff_gen_full[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_effgenpt_genptsube0[ibin][ibin2],Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2));
                //RooUnfoldResponse h_jetcharge_response_data (h_chg_reco_full_trkvary[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_effgenpt_genptsube0_trkvary[ibin][ibin2],Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold_data (&h_jetcharge_response_data, h_chg_data[ibin][ibin2], 2);
                h_chg_data_trkunfold[ibin][ibin2] = (TH1D*)unfold_data.Hreco();             

//kappa
                RooUnfoldResponse h_jetcharge_response_kappa (h_chg_eff_gen_full_kappa[ibin][ibin2], h_chg_gen_full_sube0_kappa[ibin][ibin2], h_chg_effgenpt_genptsube0_kappa[ibin][ibin2],Form("h_jetcharge_response_kappa_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_kappa_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold_kappa (&h_jetcharge_response_kappa, h_chg_eff_gen_full_kappa[ibin][ibin2], 1);
                h_chg_reco_unfolded_kappa[ibin][ibin2] = (TH1D*)unfold_kappa.Hreco();

                RooUnfoldResponse h_jetcharge_response_data_kappa (h_chg_eff_gen_full_kappa[ibin][ibin2], h_chg_gen_full_sube0_kappa[ibin][ibin2], h_chg_effgenpt_genptsube0_kappa[ibin][ibin2],Form("h_jetcharge_response_data_kappa_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_data_kappa_cent%d_trk%d",ibin,ibin2));
                RooUnfoldBayes unfold_data_kappa (&h_jetcharge_response_data_kappa, h_chg_data_kappa[ibin][ibin2], 2);
                h_chg_data_trkunfold_kappa[ibin][ibin2] = (TH1D*)unfold_data_kappa.Hreco();

/////////////////////////////
////////twosteps/////////////
/////////////////////////////
/*
//reco to reco sube0 
		        RooUnfoldResponse h_jetcharge_response_forbkg (h_chg_reco_full[ibin][ibin2], h_chg_reco_full_sube0[ibin][ibin2], h_chg_recopt_recoptsube0[ibin][ibin2],Form("h_jetcharge_response_forbkg_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_forbkg_cent%d_trk%d",ibin,ibin2));
		        RooUnfoldBayes unfold_forbkg (&h_jetcharge_response_forbkg, h_chg_reco_full[ibin][ibin2], 1);
		        h_chg_reco_unfolded_forbkg[ibin][ibin2] = (TH1D*)unfold_forbkg.Hreco();

		        RooUnfoldResponse h_jetcharge_response_forbkg_data (h_chg_reco_full[ibin][ibin2], h_chg_reco_full_sube0[ibin][ibin2], h_chg_recopt_recoptsube0[ibin][ibin2],Form("h_jetcharge_response_forbkg_data_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_forbkg_data_cent%d_trk%d",ibin,ibin2));
		        RooUnfoldBayes unfold_forbkg_data (&h_jetcharge_response_forbkg_data, h_chg_data[ibin][ibin2], 2);
		        h_chg_data_bkgunfold[ibin][ibin2] = (TH1D*)unfold_forbkg_data.Hreco();

//reco sube0 to gen sube0
		        RooUnfoldResponse h_jetcharge_response (h_chg_reco_unfolded_forbkg[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_recoptsube0_genptsube0[ibin][ibin2],Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2));
		        //RooUnfoldResponse h_jetcharge_response (h_chg_reco_full_sube0[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_recopt_genpt[ibin][ibin2],Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_cent%d_trk%d",ibin,ibin2));
		        RooUnfoldBayes unfold (&h_jetcharge_response, h_chg_reco_unfolded_forbkg[ibin][ibin2], 1);
		        h_chg_reco_unfolded[ibin][ibin2] = (TH1D*)unfold.Hreco();

		        RooUnfoldResponse h_jetcharge_response_data (h_chg_reco_unfolded_forbkg[ibin][ibin2], h_chg_gen_full_sube0[ibin][ibin2], h_chg_recoptsube0_genptsube0[ibin][ibin2],Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2),Form("h_jetcharge_response_data_cent%d_trk%d",ibin,ibin2));
		        RooUnfoldBayes unfold_data (&h_jetcharge_response_data, h_chg_data_bkgunfold[ibin][ibin2], 2);
		        h_chg_data_trkunfold[ibin][ibin2] = (TH1D*)unfold_data.Hreco();		    	
*/
		    }
		}
    }

	//TH1D* hReco= (TH1D*) unfold.Hreco();

	//unfold.PrintTable (cout, h_chg_reco_full[0][0]);
	//unfold_forbkg.PrintTable (cout, h_chg_gen_forbkg[1][1]);
   
    TLine *tl1 = new TLine(-1.2,1.,1.2,1.);
    tl1->SetLineStyle(2);

    TCanvas *c_chg_reco_unfolded_fortrk = new TCanvas("c_chg_reco_unfolded_fortrk","c_chg_reco_unfolded_fortrk",1500,600);
    c_chg_reco_unfolded_fortrk->Divide(5,2,0);
    gStyle->SetOptStat(0);

    for(int ibin=0; ibin<nCBins; ibin++){

	    if(ibin==0) c_chg_reco_unfolded_fortrk->cd(1);
	    else c_chg_reco_unfolded_fortrk->cd(6-ibin);

	    //h_chg_reco_full_sube0[ibin][0]->Scale(1./h_chg_reco_full_sube0[ibin][0]->Integral());
        h_chg_reco_full[ibin][0]->Scale(1./h_chg_reco_full[ibin][0]->Integral());
	    h_chg_gen_full_sube0[ibin][0]->Scale(1./h_chg_gen_full_sube0[ibin][0]->Integral());
	    h_chg_reco_unfolded[ibin][0]->Scale(1./h_chg_reco_unfolded[ibin][0]->Integral());
	    h_chg_data_trkunfold[ibin][0]->Scale(1./h_chg_data_trkunfold[ibin][0]->Integral());
	    h_chg_data[ibin][0]->Scale(1./h_chg_data[ibin][0]->Integral());
	    h_chg_data_unfolded_rand[ibin][0]->Scale(1./h_chg_data_unfolded_rand[ibin][0]->Integral());
	    h_chg_reco_unfolded_rand[ibin][0]->Scale(1./h_chg_reco_unfolded_rand[ibin][0]->Integral());

        h_chg_reco_unfolded[ibin][0]->GetYaxis()->SetTitle("1/N_{jets}");
        h_chg_reco_unfolded[ibin][0]->GetYaxis()->SetTitleSize(0.06);
        h_chg_reco_unfolded[ibin][0]->GetYaxis()->CenterTitle();
        h_chg_reco_unfolded[ibin][0]->GetYaxis()->SetTitleOffset(0.8);
	    h_chg_reco_unfolded[ibin][0]->GetXaxis()->SetRangeUser(-1.2,1.2);
	    h_chg_reco_unfolded[ibin][0]->GetYaxis()->SetRangeUser(0.,0.13);
        h_chg_reco_full[ibin][0]->SetLineColor(kBlue-2); h_chg_reco_full[ibin][0]->SetMarkerColor(kBlue-2); h_chg_reco_full[ibin][0]->SetMarkerStyle(20);
	    //h_chg_reco_full_sube0[ibin][0]->SetLineColor(kBlue-2); h_chg_reco_full_sube0[ibin][0]->SetMarkerColor(kBlue-2); h_chg_reco_full_sube0[ibin][0]->SetMarkerStyle(20);
	    h_chg_gen_full_sube0[ibin][0]->SetLineColor(kGreen-2); h_chg_gen_full_sube0[ibin][0]->SetMarkerColor(kGreen-2); h_chg_gen_full_sube0[ibin][0]->SetMarkerStyle(5);
	    h_chg_reco_unfolded[ibin][0]->SetLineColor(kRed-2); h_chg_reco_unfolded[ibin][0]->SetMarkerColor(kRed-2); h_chg_reco_unfolded[ibin][0]->SetMarkerStyle(20);
	    h_chg_data_trkunfold[ibin][0]->SetLineColor(kBlack); h_chg_data_trkunfold[ibin][0]->SetMarkerColor(kBlack); h_chg_data_trkunfold[ibin][0]->SetMarkerStyle(4);
	    h_chg_data[ibin][0]->SetLineColor(kBlack); h_chg_data[ibin][0]->SetMarkerColor(kBlack); h_chg_data[ibin][0]->SetMarkerStyle(34);
	    h_chg_data_unfolded_rand[ibin][0]->SetLineColor(kBlack); h_chg_data_unfolded_rand[ibin][0]->SetMarkerColor(kBlack); h_chg_data_unfolded_rand[ibin][0]->SetMarkerStyle(5);
	    h_chg_reco_unfolded_rand[ibin][0]->SetLineColor(kBlue-2); h_chg_reco_unfolded_rand[ibin][0]->SetMarkerColor(kBlue-2); h_chg_reco_unfolded_rand[ibin][0]->SetMarkerStyle(21);

        h_chg_reco_unfolded[ibin][0]->Draw("same");
        h_chg_gen_full_sube0[ibin][0]->Draw("same");
	    //h_chg_reco_full_sube0[ibin][0]->Draw("same");
        h_chg_reco_full[ibin][0]->Draw("same");
	    h_chg_data_trkunfold[ibin][0]->Draw("same");
	    h_chg_data[ibin][0]->Draw("same");
	    //h_chg_data_unfolded_rand[ibin][0]->Draw("same");
	    //h_chg_reco_unfolded_rand[ibin][0]->Draw("same");

	    TString tmp1 = cent_tag_pp_MC[ibin];
	    tx1 = new TLatex(); tx1->SetTextSize(.1);
	    tx1->DrawLatexNDC(0.35,0.85, tmp1);

        if(ibin==1){
		    TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
			legend ->SetLineColor(kWhite);
			legend ->AddEntry((TObject*)0, kappa_tag[0], "");
			legend ->AddEntry(h_chg_reco_full[ibin][0], "Reco", "lepf");
			legend ->AddEntry(h_chg_gen_full_sube0[ibin][0], "Gen", "lepf");
			legend ->AddEntry(h_chg_reco_unfolded[ibin][0], "Reco Unfolded", "lepf");
			legend ->AddEntry(h_chg_data[ibin][0], "Data", "lepf");
			legend ->AddEntry(h_chg_data_trkunfold[ibin][0], "Data Unfolded", "lepf");
			legend ->Draw("same");
		}

	    if(ibin==0) c_chg_reco_unfolded_fortrk->cd(6);
	    else c_chg_reco_unfolded_fortrk->cd(11-ibin);

	    h_chg_reco_gen_full[ibin][0] = (TH1D*)h_chg_reco_full[ibin][0]->Clone(Form("h_chg_reco_gen_full_%d_1",ibin));
	    h_chg_reco_unfolded_gen[ibin][0] = (TH1D*)h_chg_reco_unfolded[ibin][0]->Clone(Form("h_chg_reco_unfolded_%d_1",ibin));
	    h_chg_data_unfolded_gen[ibin][0] = (TH1D*)h_chg_data_trkunfold[ibin][0]->Clone(Form("h_chg_data_unfolded_%d_1",ibin));
	    h_chg_reco_unfolded_rand_gen[ibin][0] = (TH1D*)h_chg_reco_unfolded_rand[ibin][0]->Clone(Form("h_chg_reco_unfolded_rand_%d_1",ibin));
	    h_chg_data_unfolded_rand_gen[ibin][0] = (TH1D*)h_chg_data_unfolded_rand[ibin][0]->Clone(Form("h_chg_data_unfolded_rand_%d_1",ibin));

	    h_chg_reco_gen_full[ibin][0]->Divide(h_chg_gen_full_sube0[ibin][0]);
	    h_chg_reco_unfolded_gen[ibin][0]->Divide(h_chg_gen_full_sube0[ibin][0]);
	    h_chg_data_unfolded_gen[ibin][0]->Divide(h_chg_gen_full_sube0[ibin][0]);
	    h_chg_reco_unfolded_rand_gen[ibin][0]->Divide(h_chg_gen_full_sube0[ibin][0]);
	    h_chg_data_unfolded_rand_gen[ibin][0]->Divide(h_chg_data[ibin][0]);

	    h_chg_reco_unfolded_gen[ibin][0]->GetXaxis()->SetTitle("jetcharge");
	    h_chg_reco_unfolded_gen[ibin][0]->GetXaxis()->SetTitleSize(0.05);
	    h_chg_reco_unfolded_gen[ibin][0]->GetXaxis()->SetLabelSize(0.05);
	    h_chg_reco_unfolded_gen[ibin][0]->GetXaxis()->CenterTitle();

        h_chg_reco_unfolded_gen[ibin][0]->GetYaxis()->SetTitle("Ratio");
        h_chg_reco_unfolded_gen[ibin][0]->GetYaxis()->SetTitleSize(0.06);
        h_chg_reco_unfolded_gen[ibin][0]->GetYaxis()->SetTitleOffset(0.8);
        h_chg_reco_unfolded_gen[ibin][0]->GetYaxis()->SetLabelSize(0.05);
        h_chg_reco_unfolded_gen[ibin][0]->GetYaxis()->CenterTitle();

	    //h_chg_reco_full_sube0[ibin][0]->GetXaxis()->SetRangeUser(-1.5,1.5);
	    h_chg_reco_unfolded_gen[ibin][0]->GetYaxis()->SetRangeUser(0.5,1.5);
	    //h_chg_reco_full_sube0[ibin][0]->SetLineColor(kBlue-2); h_chg_reco_full_sube0[ibin][0]->SetMarkerColor(kBlue-2); h_chg_reco_full_sube0[ibin][0]->SetMarkerStyle(4);
	    //h_chg_gen_full_sube0[ibin][0]->SetLineColor(kGreen-2); h_chg_gen_full_sube0[ibin][0]->SetMarkerColor(kGreen-2); h_chg_gen_full_sube0[ibin][0]->SetMarkerStyle(4);
	    //h_chg_reco_unfolded[ibin][0]->SetLineColor(kBlue-2); h_chg_reco_unfolded[ibin][0]->SetMarkerColor(kBlue-2); h_chg_reco_unfolded[ibin][0]->SetMarkerStyle(2);

        h_chg_reco_unfolded_gen[ibin][0]->Draw("same");
	    h_chg_reco_gen_full[ibin][0]->Draw("same");
	    h_chg_data_unfolded_gen[ibin][0]->Draw("same");
	    //h_chg_reco_unfolded_rand_gen[ibin][0]->Draw("same");
	    //h_chg_data_unfolded_rand_gen[ibin][0]->Draw("same");

        tl1->Draw("same");  

        if(ibin==1){
		    TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
			legend ->SetLineColor(kWhite);
			legend ->AddEntry(h_chg_reco_gen_full[ibin][0], "Reco / Gen", "lepf");
			legend ->AddEntry(h_chg_reco_unfolded_gen[ibin][0], "Reco Unfolded / Gen", "lepf");
			//legend ->AddEntry(h_chg_data_unfolded_gen[ibin][0], "Data Unfolded / Data", "lepf");
			//legend ->AddEntry(h_chg_reco_unfolded_rand_gen[ibin][0], "Reco Unfolded rand / Gen", "lepf");
			//legend ->AddEntry(h_chg_data_unfolded_rand_gen[ibin][0], "Data Unfolded rand / Data", "lepf");
			legend ->Draw("same");
		}
	}


//kappa
    TCanvas *c_chg_reco_unfolded_fortrk_kappa = new TCanvas("c_chg_reco_unfolded_fortrk_kappa","c_chg_reco_unfolded_fortrk_kappa",1500,600);
    c_chg_reco_unfolded_fortrk_kappa->Divide(5,2,0);
    gStyle->SetOptStat(0);

    for(int ibin=0; ibin<nCBins; ibin++){

        if(ibin==0) c_chg_reco_unfolded_fortrk_kappa->cd(1);
        else c_chg_reco_unfolded_fortrk_kappa->cd(6-ibin);

        h_chg_reco_full_kappa[ibin][0]->Scale(1./h_chg_reco_full_kappa[ibin][0]->Integral());
        //h_chg_reco_full_sube0_kappa[ibin][0]->Scale(1./h_chg_reco_full_sube0_kappa[ibin][0]->Integral());
        h_chg_gen_full_sube0_kappa[ibin][0]->Scale(1./h_chg_gen_full_sube0_kappa[ibin][0]->Integral());
        h_chg_reco_unfolded_kappa[ibin][0]->Scale(1./h_chg_reco_unfolded_kappa[ibin][0]->Integral());
        h_chg_data_trkunfold_kappa[ibin][0]->Scale(1./h_chg_data_trkunfold_kappa[ibin][0]->Integral());
        h_chg_data_kappa[ibin][0]->Scale(1./h_chg_data_kappa[ibin][0]->Integral());

        h_chg_reco_unfolded_kappa[ibin][0]->GetXaxis()->SetRangeUser(-1.5,1.5);
        h_chg_reco_unfolded_kappa[ibin][0]->GetYaxis()->SetRangeUser(0.,0.13);
        //h_chg_reco_full_sube0_kappa[ibin][0]->SetLineColor(kBlue-2); h_chg_reco_full_sube0_kappa[ibin][0]->SetMarkerColor(kBlue-2); h_chg_reco_full_sube0_kappa[ibin][0]->SetMarkerStyle(20);
        h_chg_reco_full_kappa[ibin][0]->SetLineColor(kBlue-2); h_chg_reco_full_kappa[ibin][0]->SetMarkerColor(kBlue-2); h_chg_reco_full_kappa[ibin][0]->SetMarkerStyle(20);
        h_chg_gen_full_sube0_kappa[ibin][0]->SetLineColor(kGreen-2); h_chg_gen_full_sube0_kappa[ibin][0]->SetMarkerColor(kGreen-2); h_chg_gen_full_sube0_kappa[ibin][0]->SetMarkerStyle(20);
        h_chg_reco_unfolded_kappa[ibin][0]->SetLineColor(kRed-2); h_chg_reco_unfolded_kappa[ibin][0]->SetMarkerColor(kRed-2); h_chg_reco_unfolded_kappa[ibin][0]->SetMarkerStyle(20);
        h_chg_data_trkunfold_kappa[ibin][0]->SetLineColor(kBlack); h_chg_data_trkunfold_kappa[ibin][0]->SetMarkerColor(kBlack); h_chg_data_trkunfold_kappa[ibin][0]->SetMarkerStyle(4);
        h_chg_data_kappa[ibin][0]->SetLineColor(kBlack); h_chg_data_kappa[ibin][0]->SetMarkerColor(kBlack); h_chg_data_kappa[ibin][0]->SetMarkerStyle(34);

        h_chg_reco_unfolded_kappa[ibin][0]->Draw("same");
        h_chg_reco_full_kappa[ibin][0]->Draw("same");
        //h_chg_reco_full_sube0_kappa[ibin][0]->Draw("same");
        h_chg_gen_full_sube0_kappa[ibin][0]->Draw("same");
        h_chg_data_trkunfold_kappa[ibin][0]->Draw("same");
        h_chg_data_kappa[ibin][0]->Draw("same");

        TString tmp1 = cent_tag_pp_MC[ibin];
        tx1 = new TLatex(); tx1->SetTextSize(.1);
        tx1->DrawLatexNDC(0.35,0.85, tmp1);

        if(ibin==1){
            TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
            legend ->SetLineColor(kWhite);
            legend ->AddEntry((TObject*)0, kappa_tag[0], "");
            //legend ->AddEntry(h_chg_reco_full_sube0_kappa[ibin][0], "Reco", "lepf");
            legend ->AddEntry(h_chg_reco_full_kappa[ibin][0], "Reco", "lepf");
            legend ->AddEntry(h_chg_gen_full_sube0_kappa[ibin][0], "Gen", "lepf");
            legend ->AddEntry(h_chg_reco_unfolded_kappa[ibin][0], "Reco Unfolded", "lepf");
            legend ->AddEntry(h_chg_data_kappa[ibin][0], "Data", "lepf");
            legend ->AddEntry(h_chg_data_trkunfold_kappa[ibin][0], "Data Unfolded", "lepf");
            legend ->Draw("same");
        }

        if(ibin==0) c_chg_reco_unfolded_fortrk_kappa->cd(6);
        else c_chg_reco_unfolded_fortrk_kappa->cd(11-ibin);

        h_chg_reco_gen_full_kappa[ibin][0] = (TH1D*)h_chg_reco_full_kappa[ibin][0]->Clone(Form("h_chg_reco_gen_full_%d_1",ibin));
        h_chg_reco_unfolded_gen_kappa[ibin][0] = (TH1D*)h_chg_reco_unfolded_kappa[ibin][0]->Clone(Form("h_chg_reco_unfolded_%d_1",ibin));
        h_chg_data_unfolded_gen_kappa[ibin][0] = (TH1D*)h_chg_data_trkunfold_kappa[ibin][0]->Clone(Form("h_chg_data_unfolded_%d_1",ibin));

        h_chg_reco_gen_full_kappa[ibin][0]->Divide(h_chg_gen_full_sube0_kappa[ibin][0]);
        h_chg_reco_unfolded_gen_kappa[ibin][0]->Divide(h_chg_gen_full_sube0_kappa[ibin][0]);
        h_chg_data_unfolded_gen_kappa[ibin][0]->Divide(h_chg_gen_full_sube0_kappa[ibin][0]);

        h_chg_reco_gen_full_kappa[ibin][0]->GetXaxis()->SetTitle("jetcharge");
        h_chg_reco_gen_full_kappa[ibin][0]->GetXaxis()->SetTitleSize(0.05);
        h_chg_reco_gen_full_kappa[ibin][0]->GetXaxis()->SetLabelSize(0.05);
        h_chg_reco_gen_full_kappa[ibin][0]->GetXaxis()->CenterTitle();

        //h_chg_reco_full_sube0_kappa[ibin][0]->GetXaxis()->SetRangeUser(-1.5,1.5);
        h_chg_reco_unfolded_gen_kappa[ibin][0]->GetYaxis()->SetRangeUser(0.,2.);
        //h_chg_reco_full_sube0_kappa[ibin][0]->SetLineColor(kBlue-2); h_chg_reco_full_sube0_kappa[ibin][0]->SetMarkerColor(kBlue-2); h_chg_reco_full_sube0_kappa[ibin][0]->SetMarkerStyle(4);
        //h_chg_gen_full_sube0_kappa[ibin][0]->SetLineColor(kGreen-2); h_chg_gen_full_sube0_kappa[ibin][0]->SetMarkerColor(kGreen-2); h_chg_gen_full_sube0_kappa[ibin][0]->SetMarkerStyle(4);
        //h_chg_reco_unfolded_kappa[ibin][0]->SetLineColor(kBlue-2); h_chg_reco_unfolded_kappa[ibin][0]->SetMarkerColor(kBlue-2); h_chg_reco_unfolded_kappa[ibin][0]->SetMarkerStyle(2);

        h_chg_reco_unfolded_gen_kappa[ibin][0]->Draw("same");
        h_chg_reco_gen_full_kappa[ibin][0]->Draw("same");
        h_chg_data_unfolded_gen_kappa[ibin][0]->Draw("same");

        if(ibin==1){
            TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
            legend ->SetLineColor(kWhite);
            legend ->AddEntry(h_chg_reco_gen_full_kappa[ibin][0], "Reco / Gen", "lepf");
            legend ->AddEntry(h_chg_reco_unfolded_gen_kappa[ibin][0], "Reco Unfolded / Gen", "lepf");
            legend ->AddEntry(h_chg_data_unfolded_gen_kappa[ibin][0], "Data Unfolded / Data", "lepf");
            legend ->Draw("same");
        }
    }
/*
///unfolding for bkg
    TCanvas *c_chg_reco_unfolded_forbkg = new TCanvas("c_chg_reco_unfolded_forbkg","c_chg_reco_unfolded_forbkg",1200,600);
    c_chg_reco_unfolded_forbkg->Divide(4,2,0);
    gStyle->SetOptStat(0);

    for(int ibin=1; ibin<nCBins; ibin++){

	    if(ibin==0) c_chg_reco_unfolded_forbkg->cd(1);
	    else c_chg_reco_unfolded_forbkg->cd(5-ibin);

	    //h_chg_reco_full[ibin][0]->Scale(1./h_chg_reco_full[ibin][0]->Integral());
	    h_chg_reco_full_sube0[ibin][0]->Scale(1./h_chg_reco_full_sube0[ibin][0]->Integral());
	    h_chg_reco_unfolded_forbkg[ibin][0]->Scale(1./h_chg_reco_unfolded_forbkg[ibin][0]->Integral());
	    h_chg_data_bkgunfold[ibin][0]->Scale(1./h_chg_data_bkgunfold[ibin][0]->Integral());
	    //h_chg_reco_unfolded_forbkg_rand[ibin][0]->Scale(1./h_chg_reco_unfolded_forbkg_rand[ibin][0]->Integral());
	    //h_chg_data_unfolded_forbkg_rand[ibin][0]->Scale(1./h_chg_data_unfolded_forbkg_rand[ibin][0]->Integral());

	    h_chg_reco_full[ibin][0]->GetXaxis()->SetRangeUser(-1.2,1.2);
	    h_chg_reco_full[ibin][0]->GetYaxis()->SetRangeUser(0.,0.13);
	    h_chg_reco_full[ibin][0]->SetLineColor(kBlue-2); h_chg_reco_full[ibin][0]->SetMarkerColor(kBlue-2); h_chg_reco_full[ibin][0]->SetMarkerStyle(20);
	    h_chg_reco_full_sube0[ibin][0]->SetLineColor(kGreen-2); h_chg_reco_full_sube0[ibin][0]->SetMarkerColor(kGreen-2); h_chg_reco_full_sube0[ibin][0]->SetMarkerStyle(5);
	    h_chg_reco_unfolded_forbkg[ibin][0]->SetLineColor(kRed-2); h_chg_reco_unfolded_forbkg[ibin][0]->SetMarkerColor(kRed-2); h_chg_reco_unfolded_forbkg[ibin][0]->SetMarkerStyle(20);
	    //h_chg_reco_unfolded_forbkg_rand[ibin][0]->SetLineColor(kRed-2); h_chg_reco_unfolded_forbkg_rand[ibin][0]->SetMarkerColor(kRed-2); h_chg_reco_unfolded_forbkg_rand[ibin][0]->SetMarkerStyle(5);

	    h_chg_reco_full[ibin][0]->Draw("same");
	    h_chg_reco_full_sube0[ibin][0]->Draw("same");
	    h_chg_reco_unfolded_forbkg[ibin][0]->Draw("same");
        //h_chg_data[ibin][0]->Draw("same");
	    //h_chg_data_bkgunfold[ibin][0]->Draw("same");
	    //h_chg_reco_unfolded_forbkg_rand[ibin][0]->Draw("same");

	    TString tmp1 = cent_tag_pp_MC[ibin];
	    tx1 = new TLatex(); tx1->SetTextSize(.1);
	    tx1->DrawLatexNDC(0.35,0.85, tmp1);

        if(ibin==1){
		    TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
			legend ->SetLineColor(kWhite);
			legend ->AddEntry((TObject*)0, kappa_tag[0], "");
			legend ->AddEntry(h_chg_reco_full[ibin][0], "Reco", "lepf");
			legend ->AddEntry(h_chg_reco_full_sube0[ibin][0], "Reco sube0", "lepf");
			legend ->AddEntry(h_chg_reco_unfolded_forbkg[ibin][0], "Reco Unfolded", "lepf");
            //legend ->AddEntry(h_chg_data[ibin][0], "Data", "lepf");
            //legend ->AddEntry(h_chg_data_bkgunfold[ibin][0], "Data Unfolded", "lepf");
			//legend ->AddEntry(h_chg_reco_unfolded_forbkg_rand[ibin][0], "Reco Unfolded rand", "lepf");
			legend ->Draw("same");
		}

	    if(ibin==0) c_chg_reco_unfolded_forbkg->cd(6);
	    else c_chg_reco_unfolded_forbkg->cd(9-ibin);

	    h_chg_reco_gen_full_forbkg[ibin][0] = (TH1D*)h_chg_reco_full[ibin][0]->Clone(Form("h_chg_reco_gen_full_forbkg_%d_1",ibin));
	    h_chg_reco_unfolded_gen_forbkg[ibin][0] = (TH1D*)h_chg_reco_unfolded_forbkg[ibin][0]->Clone(Form("h_chg_reco_unfolded_forbkg_%d_1",ibin));
        h_chg_data_unfolded_gen_forbkg[ibin][0] = (TH1D*)h_chg_data_bkgunfold[ibin][0]->Clone(Form("h_chg_data_unfolded_forbkg_%d_1",ibin));
	    //h_chg_reco_unfolded_rand_gen_forbkg[ibin][0] = (TH1D*)h_chg_reco_unfolded_forbkg_rand[ibin][0]->Clone(Form("h_chg_reco_unfolded_forbkg_rand_%d_1",ibin));

	    h_chg_reco_gen_full_forbkg[ibin][0]->Divide(h_chg_reco_full_sube0[ibin][0]);
	    h_chg_reco_unfolded_gen_forbkg[ibin][0]->Divide(h_chg_reco_full_sube0[ibin][0]);
        h_chg_data_unfolded_gen_forbkg[ibin][0]->Divide(h_chg_reco_full_sube0[ibin][0]);
	    //h_chg_reco_unfolded_rand_gen_forbkg[ibin][0]->Divide(h_chg_reco_full_sube0[ibin][0]);

	    h_chg_reco_gen_full_forbkg[ibin][0]->GetXaxis()->SetTitle("jetcharge");
	    h_chg_reco_gen_full_forbkg[ibin][0]->GetXaxis()->SetTitleSize(0.05);
	    h_chg_reco_gen_full_forbkg[ibin][0]->GetXaxis()->SetLabelSize(0.05);
	    h_chg_reco_gen_full_forbkg[ibin][0]->GetXaxis()->CenterTitle();

	    //h_chg_reco_full_sube0[ibin][0]->GetXaxis()->SetRangeUser(-1.5,1.5);
	    h_chg_reco_gen_full_forbkg[ibin][0]->GetYaxis()->SetRangeUser(0.,2.);
	    //h_chg_reco_full_sube0[ibin][0]->SetLineColor(kBlue-2); h_chg_reco_full_sube0[ibin][0]->SetMarkerColor(kBlue-2); h_chg_reco_full_sube0[ibin][0]->SetMarkerStyle(4);
	    //h_chg_gen_full_sube0[ibin][0]->SetLineColor(kGreen-2); h_chg_gen_full_sube0[ibin][0]->SetMarkerColor(kGreen-2); h_chg_gen_full_sube0[ibin][0]->SetMarkerStyle(4);
	    //h_chg_reco_unfolded[ibin][0]->SetLineColor(kBlue-2); h_chg_reco_unfolded[ibin][0]->SetMarkerColor(kBlue-2); h_chg_reco_unfolded[ibin][0]->SetMarkerStyle(2);

	    h_chg_reco_gen_full_forbkg[ibin][0]->Draw("same");
	    h_chg_reco_unfolded_gen_forbkg[ibin][0]->Draw("same");
        //h_chg_data_unfolded_gen_forbkg[ibin][0]->Draw("same");
	    //h_chg_reco_unfolded_rand_gen_forbkg[ibin][0]->Draw("same");

        if(ibin==1){
		    TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
			legend ->SetLineColor(kWhite);
			legend ->AddEntry(h_chg_reco_gen_full_forbkg[ibin][0], "Reco / Reco sube0", "lepf");
			legend ->AddEntry(h_chg_reco_unfolded_gen_forbkg[ibin][0], "Reco Unfolded / Reco sube0", "lepf");
            //legend ->AddEntry(h_chg_reco_unfolded_gen_forbkg[ibin][0], "Data Unfolded / Reco sube0", "lepf");
			//legend ->AddEntry(h_chg_reco_unfolded_rand_gen_forbkg[ibin][0], "Reco Unfolded rand / Reco sube0", "lepf");
			legend ->Draw("same");
		}
	}

///unfolding for MC
    TCanvas *c_chg_reco_unfolded = new TCanvas("c_chg_reco_unfolded","c_chg_reco_unfolded",1500,300);
    c_chg_reco_unfolded->Divide(5,1,0);
    gStyle->SetOptStat(0);

    for(int ibin=0; ibin<nCBins; ibin++){

	    if(ibin==0) c_chg_reco_unfolded->cd(1);
	    else c_chg_reco_unfolded->cd(6-ibin);

	    if(ibin==0){
	      //h_chg_data_bkgunfold[ibin][0] = (TH1D*)h_chg_data_trkunfold[ibin][0];	
	    }

	    h_chg_reco_full[ibin][0]->Scale(1./h_chg_reco_full[ibin][0]->Integral());
	    h_chg_reco_unfolded[ibin][0]->Scale(1./h_chg_reco_unfolded[ibin][0]->Integral());

	    h_chg_reco_full[ibin][0]->GetXaxis()->SetTitle("jetcharge");
	    h_chg_reco_full[ibin][0]->GetXaxis()->SetTitleSize(0.05);
	    h_chg_reco_full[ibin][0]->GetXaxis()->SetLabelSize(0.05);
	    h_chg_reco_full[ibin][0]->GetXaxis()->CenterTitle();
	    //h_chg_reco_full[ibin][0]->GetYaxis()->SetTitle("1/N_{jets} (dN/dQ) [1/e]");
	    //h_chg_reco_full[ibin][0]->GetYaxis()->CenterTitle();

	    h_chg_reco_full[ibin][0]->GetXaxis()->SetRangeUser(-1.5,1.5);
	    h_chg_reco_full[ibin][0]->GetYaxis()->SetRangeUser(0.,0.13);
	    h_chg_reco_full[ibin][0]->SetLineColor(kRed); h_chg_reco_full[ibin][0]->SetLineStyle(1); h_chg_reco_full[ibin][0]->SetMarkerColor(kRed); h_chg_reco_full[ibin][0]->SetMarkerStyle(20);
	    h_chg_reco_unfolded[ibin][0]->SetLineColor(kBlue); h_chg_reco_unfolded[ibin][0]->SetMarkerColor(kBlue); h_chg_reco_unfolded[ibin][0]->SetMarkerStyle(20);

	    h_chg_reco_full[ibin][0]->Draw("hist same");
	    h_chg_reco_unfolded[ibin][0]->Draw("hist same");

	    TString tmp1 = cent_tag_pp_MC[ibin];
	    tx1 = new TLatex(); tx1->SetTextSize(.1);
	    tx1->DrawLatexNDC(0.35,0.85, tmp1);

        if(ibin==1){
		    TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
			legend ->SetLineColor(kWhite);
			legend ->AddEntry((TObject*)0, kappa_tag[0], "");
			legend ->AddEntry(h_chg_reco_full[ibin][0], "Reco", "l");
			legend ->AddEntry(h_chg_reco_unfolded[ibin][0], "Reco Unfolded", "l");
			legend ->Draw("same");
		}
	}
*/
/*
///unfolding for data
    TCanvas *c_chg_data_unfolded = new TCanvas("c_chg_data_unfolded","c_chg_data_unfolded",1500,300);
    c_chg_data_unfolded->Divide(5,1,0);
    gStyle->SetOptStat(0);

    for(int ibin=0; ibin<nCBins; ibin++){

	    if(ibin==0) c_chg_data_unfolded->cd(1);
	    else c_chg_data_unfolded->cd(6-ibin);

	    if(ibin==0){
	      //h_chg_data_bkgunfold[ibin][0] = (TH1D*)h_chg_data_trkunfold[ibin][0];	
	    }

	    h_chg_data[ibin][0]->Scale(1./h_chg_data[ibin][0]->Integral());
	    h_chg_data_trkunfold[ibin][0]->Scale(1./h_chg_data_trkunfold[ibin][0]->Integral());

	    h_chg_data[ibin][0]->GetXaxis()->SetTitle("jetcharge");
        h_chg_data[ibin][0]->GetYaxis()->SetTitle("1/N_{jets} (dN/dQ) [1/e]");
	    h_chg_data[ibin][0]->GetXaxis()->SetTitleSize(0.05);
	    h_chg_data[ibin][0]->GetXaxis()->SetLabelSize(0.05);
	    h_chg_data[ibin][0]->GetXaxis()->CenterTitle();
	    //h_chg_data[ibin][0]->GetYaxis()->SetTitle("1/N_{jets} (dN/dQ) [1/e]");
	    //h_chg_data[ibin][0]->GetYaxis()->CenterTitle();

	    h_chg_data[ibin][0]->GetXaxis()->SetRangeUser(-1.2,1.2);
	    h_chg_data[ibin][0]->GetYaxis()->SetRangeUser(0.,0.13);
	    h_chg_data[ibin][0]->SetLineColor(kRed); h_chg_data[ibin][0]->SetLineStyle(1); h_chg_data[ibin][0]->SetMarkerColor(kRed); h_chg_data[ibin][0]->SetMarkerStyle(20);
	    h_chg_data_trkunfold[ibin][0]->SetLineColor(kBlue); h_chg_data_trkunfold[ibin][0]->SetMarkerColor(kBlue); h_chg_data_trkunfold[ibin][0]->SetMarkerStyle(20);

	    h_chg_data[ibin][0]->Draw("hist same");
	    h_chg_data_trkunfold[ibin][0]->Draw("hist same");

	    TString tmp1 = cent_tag_pp_data[ibin];
	    tx1 = new TLatex(); tx1->SetTextSize(.1);
	    tx1->DrawLatexNDC(0.35,0.85, tmp1);

        if(ibin==1){
		    TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
			legend ->SetLineColor(kWhite);
			legend ->AddEntry((TObject*)0, kappa_tag[0], "");
			legend ->AddEntry(h_chg_data[ibin][0], "Detector level Data", "l");
			legend ->AddEntry(h_chg_data_trkunfold[ibin][0], "Unfolded Data", "l");
			legend ->Draw("same");
		}
	}
/*
    TCanvas *c_chg_recopt_genpt = new TCanvas("c_chg_recopt_genpt","c_chg_recopt_genpt",1500,300);
    c_chg_recopt_genpt->Divide(5,1);
    gStyle->SetOptStat(0);

    for(int ibin=0; ibin<nCBins; ibin++){
        if(ibin==0) c_chg_recopt_genpt->cd(1);
        else c_chg_recopt_genpt->cd(6-ibin);
        gPad->SetLogz();

	    h_chg_recopt_genpt[ibin][1]->GetXaxis()->SetRangeUser(-2.,2.);
	    h_chg_recopt_genpt[ibin][1]->GetYaxis()->SetRangeUser(-2.,2.);
	    h_chg_recopt_genpt[ibin][1]->GetYaxis()->SetTitleSize(0.05);
	    h_chg_recopt_genpt[ibin][1]->GetXaxis()->SetTitleSize(0.05);
	    h_chg_recopt_genpt[ibin][1]->GetXaxis()->SetTitle("recosube0 jetcharge");
	    h_chg_recopt_genpt[ibin][1]->GetYaxis()->SetTitle("gensube0 jetcharge");
	    h_chg_recopt_genpt[ibin][1]->Draw("COLZ");

	    TString tmp1 = cent_tag_pp_MC[ibin];
	    tx1 = new TLatex(); tx1->SetTextSize(.1);
	    tx1->DrawLatexNDC(0.35,0.85, tmp1);

	}
*/
    TCanvas *c_chg_recopt_recoptsube0 = new TCanvas("c_chg_recopt_recoptsube0","c_chg_recopt_recoptsube0",1200,300);
    c_chg_recopt_recoptsube0->Divide(4,1);
    gStyle->SetOptStat(0);

    for(int ibin=1; ibin<nCBins; ibin++){
        if(ibin==0) c_chg_recopt_recoptsube0->cd(1);
        else c_chg_recopt_recoptsube0->cd(5-ibin);
        //gPad->SetLogz();  

	    h_chg_recopt_recoptsube0_kappa[ibin][0]->GetXaxis()->SetRangeUser(-1.2,1.2);
	    h_chg_recopt_recoptsube0_kappa[ibin][0]->GetYaxis()->SetRangeUser(-1.2,1.2);
	    h_chg_recopt_recoptsube0_kappa[ibin][0]->GetYaxis()->SetTitleSize(0.06);
	    h_chg_recopt_recoptsube0_kappa[ibin][0]->GetXaxis()->SetTitleSize(0.06);
        h_chg_recopt_recoptsube0_kappa[ibin][0]->GetYaxis()->SetTitleOffset(0.75);
        h_chg_recopt_recoptsube0_kappa[ibin][0]->GetXaxis()->SetTitleOffset(0.75);
        h_chg_recopt_recoptsube0_kappa[ibin][0]->GetYaxis()->CenterTitle();
        h_chg_recopt_recoptsube0_kappa[ibin][0]->GetXaxis()->CenterTitle();
        h_chg_recopt_recoptsube0_kappa[ibin][0]->GetYaxis()->SetLabelSize(0.05);
        h_chg_recopt_recoptsube0_kappa[ibin][0]->GetXaxis()->SetLabelSize(0.05);
	    h_chg_recopt_recoptsube0_kappa[ibin][0]->GetXaxis()->SetTitle("Reco jet charge (Signal)");
	    h_chg_recopt_recoptsube0_kappa[ibin][0]->GetYaxis()->SetTitle("Reco jet charge (Signal+Bkg)");
	    h_chg_recopt_recoptsube0_kappa[ibin][0]->Draw("COLZ");

	    TString tmp1 = cent_tag_pp_MC[ibin];
	    tx1 = new TLatex(); tx1->SetTextSize(.09);
	    tx1->DrawLatexNDC(0.25,0.91, tmp1);

        if(ibin==1){
            TString tmp = kappa_tag[0];
            tx = new TLatex(); tx->SetTextSize(.08);
            tx->DrawLatexNDC(0.15,0.83, tmp);
        }
	}

    TCanvas *c_chg_recopt_genptsube0 = new TCanvas("c_chg_recopt_genptsube0","c_chg_recopt_genptsube0",1500,300);
    c_chg_recopt_genptsube0->Divide(5,1);
    gStyle->SetOptStat(0);

    for(int ibin=0; ibin<nCBins; ibin++){
        if(ibin==0) c_chg_recopt_genptsube0->cd(1);
        else c_chg_recopt_genptsube0->cd(6-ibin);
        //gPad->SetLogz();  

	    h_chg_recoptsube0_genptsube0_kappa[ibin][0]->GetXaxis()->SetRangeUser(-1.2,1.2);
	    h_chg_recoptsube0_genptsube0_kappa[ibin][0]->GetYaxis()->SetRangeUser(-1.2,1.2);
	    h_chg_recoptsube0_genptsube0_kappa[ibin][0]->GetYaxis()->SetTitleSize(0.06);
	    h_chg_recoptsube0_genptsube0_kappa[ibin][0]->GetXaxis()->SetTitleSize(0.06);
        h_chg_recoptsube0_genptsube0_kappa[ibin][0]->GetYaxis()->SetTitleOffset(0.75);
        h_chg_recoptsube0_genptsube0_kappa[ibin][0]->GetXaxis()->SetTitleOffset(0.75);
        h_chg_recoptsube0_genptsube0_kappa[ibin][0]->GetYaxis()->CenterTitle();
        h_chg_recoptsube0_genptsube0_kappa[ibin][0]->GetXaxis()->CenterTitle();
        h_chg_recoptsube0_genptsube0_kappa[ibin][0]->GetYaxis()->SetLabelSize(0.05);
        h_chg_recoptsube0_genptsube0_kappa[ibin][0]->GetXaxis()->SetLabelSize(0.05);
	    h_chg_recoptsube0_genptsube0_kappa[ibin][0]->GetXaxis()->SetTitle("Reco jet charge (signal)");
	    h_chg_recoptsube0_genptsube0_kappa[ibin][0]->GetYaxis()->SetTitle("Gen jet charge (signal)");
	    h_chg_recoptsube0_genptsube0_kappa[ibin][0]->Draw("COLZ");

	    TString tmp1 = cent_tag_pp_MC[ibin];
	    tx1 = new TLatex(); tx1->SetTextSize(.09);
	    tx1->DrawLatexNDC(0.25,0.91, tmp1);

        if(ibin==1){
            TString tmp = kappa_tag[0];
            tx = new TLatex(); tx->SetTextSize(.08);
            tx->DrawLatexNDC(0.15,0.83, tmp);
        }
	}
/*
    TCanvas *c_chg_data_MC = new TCanvas("c_chg_data_MC","c_chg_data_MC",1500,600);
    c_chg_data_MC->Divide(5,2,0);
    gStyle->SetOptStat(0);

    for(int ibin=0; ibin<nCBins; ibin++){

	    if(ibin==0) c_chg_data_MC->cd(1);
	    else c_chg_data_MC->cd(6-ibin);

	    h_chg_reco_unfolded[ibin][0]->Scale(1./h_chg_reco_unfolded[ibin][0]->Integral());
	    h_chg_data_trkunfold[ibin][0]->Scale(1./h_chg_data_trkunfold[ibin][0]->Integral());

	    h_chg_reco_unfolded[ibin][0]->GetXaxis()->SetRangeUser(-2.,2.);
	    h_chg_reco_unfolded[ibin][0]->GetYaxis()->SetRangeUser(0.,0.14);
	    h_chg_reco_unfolded[ibin][0]->SetLineColor(kBlue-2); h_chg_reco_unfolded[ibin][0]->SetMarkerColor(kBlue-2); h_chg_reco_unfolded[ibin][0]->SetMarkerStyle(20);
	    h_chg_data_trkunfold[ibin][0]->SetLineColor(kBlack); h_chg_data_trkunfold[ibin][0]->SetMarkerColor(kBlack); h_chg_data_trkunfold[ibin][0]->SetMarkerStyle(5);

	    h_chg_reco_unfolded[ibin][0]->Draw("same");
	    h_chg_data_trkunfold[ibin][0]->Draw("same");

	    TString tmp1 = cent_tag_pp_MC[ibin];
	    tx1 = new TLatex(); tx1->SetTextSize(.1);
	    tx1->DrawLatexNDC(0.35,0.85, tmp1);

        if(ibin==1){
		    TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
			legend ->SetLineColor(kWhite);
			legend ->AddEntry((TObject*)0, kappa_tag[0], "");
			legend ->AddEntry(h_chg_reco_unfolded[ibin][0], "Reco", "lepf");
			legend ->AddEntry(h_chg_data_trkunfold[ibin][0], "Data", "lepf");
			legend ->Draw("same");
		}

	    if(ibin==0) c_chg_data_MC->cd(6);
	    else c_chg_data_MC->cd(11-ibin);

	    h_chg_mc_data[ibin][0] = (TH1D*)h_chg_reco_unfolded[ibin][0]->Clone(Form("h_chg_mc_data_%d_1",ibin));
	    h_chg_mc_data[ibin][0]->Divide(h_chg_data_trkunfold[ibin][0]);

	    h_chg_mc_data[ibin][0]->GetXaxis()->SetTitle("jetcharge");
	    h_chg_mc_data[ibin][0]->GetXaxis()->SetTitleSize(0.05);
	    h_chg_mc_data[ibin][0]->GetXaxis()->SetLabelSize(0.05);
	    h_chg_mc_data[ibin][0]->GetXaxis()->CenterTitle();

        h_chg_mc_data[ibin][0]->GetXaxis()->SetRangeUser(-2.,2.);
	    h_chg_mc_data[ibin][0]->GetYaxis()->SetRangeUser(0.5,1.5);
	    h_chg_mc_data[ibin][0]->Draw("same");

        if(ibin==1){
		    TLegend *legend = new TLegend(0.25,0.75,0.75,0.9);
			legend ->SetLineColor(kWhite);
			legend ->AddEntry(h_chg_mc_data[ibin][0], "Reco / Data", "lepf");
			legend ->Draw("same");
		}
	}
*/
	TFile *unfolded_histos = new TFile(Form("data_fulleta_unfolded_trk_jetchg_trkeffconst_kappa_recogen_%d.root",date->GetDate()), "RECREATE");
	unfolded_histos->cd();

    for(int ibin=0; ibin<nCBins; ibin++){
       for(int ibin2=0; ibin2<ntrkbins; ibin2++){
    	   h_chg_data_trkunfold[ibin][ibin2]->Scale(1./h_chg_data_trkunfold[ibin][ibin2]->Integral());
           h_chg_data_trkunfold[ibin][ibin2]->Write();
    	   h_chg_reco_unfolded[ibin][ibin2]->Scale(1./h_chg_reco_unfolded[ibin][ibin2]->Integral());
           h_chg_reco_unfolded[ibin][ibin2]->Write();
           h_chg_data_trkunfold_kappa[ibin][ibin2]->Scale(1./h_chg_data_trkunfold_kappa[ibin][ibin2]->Integral());
           h_chg_data_trkunfold_kappa[ibin][ibin2]->Write();
           h_chg_reco_unfolded_kappa[ibin][ibin2]->Scale(1./h_chg_reco_unfolded_kappa[ibin][ibin2]->Integral());
           h_chg_reco_unfolded_kappa[ibin][ibin2]->Write();
       }
    }
    
	unfolded_histos->Close();

}