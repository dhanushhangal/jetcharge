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

TString cent[5] = {"0","1","2","3","4"};
TString trk[4] = {"0","1","2","3"};
TString kbin[4] = {"0","1","2","3"};
TString pt[56] = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55"};

TString cent_tag[] = {"0-10%","10-30%","30-50%", "50-100%"};
TString trk_tag[] = {"p_{T}^{trk} > 0.7 GeV","p_{T}^{trk} > 2 GeV","p_{T}^{trk} > 4 GeV", "p_{T}^{trk} > 5 GeV"};


void plot_charge(){

  //TFile *closure_histos_pp_MC_kscan = TFile::Open("/home/dhanush/Documents/jet_chg/Pythia6_jetchg_Mar31.root");
  TFile *closure_histos_pp_MC = TFile::Open("Pythia6_jetchg_bkgsub_eta1p5_Apr24.root");
  TFile *closure_histos_PbPb_MC = TFile::Open("P+H_jetchg_cymbal_bkgsub_eta1p5_Apr24.root");
  //TFile *closure_histos_pp_data = TFile::Open("/home/dhanush/Documents/jet_chg/ppdata_jetchg_bkgsub_eta0p5_1p5_Apr16.root");
  //TFile *closure_histos_PbPb_data = TFile::Open("/home/dhanush/Documents/jet_chg/PbPbdata_jetchg_bkgsub_eta0p5_1p5_Apr16.root");

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

  TH1D *h_reco_up_inc_ratio[nCBins];
  TH1D *h_reco_upbar_inc_ratio[nCBins];  
  TH1D *h_reco_d_inc_ratio[nCBins];
  TH1D *h_reco_dbar_inc_ratio[nCBins];
  TH1D *h_reco_gluon_inc_ratio[nCBins];

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

      h_eta_MC[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_cent"+cent[ibin]));       
      h_eta_MC_q[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_q_cent"+cent[ibin]))->Clone((TString)("h_eta_q_cent"+cent[ibin]));
      h_eta_MC_g[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_g_cent"+cent[ibin]))->Clone((TString)("h_eta_g_cent"+cent[ibin]));
      //h_eta_data[ibin] = (TH1D*)closure_histos_pp_data->Get((TString)("h_eta_full_cent"+cent[ibin]));
      h_eta_MC_up[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_up_cent"+cent[ibin]));
      h_eta_MC_upbar[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_upbar_cent"+cent[ibin]));
      h_eta_MC_d[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_d_cent"+cent[ibin]));
      h_eta_MC_dbar[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_eta_full_dbar_cent"+cent[ibin]));

      h_reco_MC[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_cent"+cent[ibin]));
      //h_reco_data[ibin] = (TH1D*)closure_histos_pp_data->Get((TString)("h_reco_corr_cent"+cent[ibin]));
      h_reco_MC_q[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_q_cent"+cent[ibin]));
      h_reco_MC_g[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_g_cent"+cent[ibin]));
      h_reco_MC_up[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_up_cent"+cent[ibin]));
      h_reco_MC_upbar[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_upbar_cent"+cent[ibin]));
      h_reco_MC_d[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_d_cent"+cent[ibin]));
      h_reco_MC_dbar[ibin] = (TH1D*)closure_histos_pp_MC->Get((TString)("h_reco_corr_dbar_cent"+cent[ibin]));
    }
    else{

      h_eta_MC[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_cent"+cent[ibin-1]));       
      h_eta_MC_q[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_q_cent"+cent[ibin-1]))->Clone((TString)("h_eta_q_cent"+cent[ibin]));
      h_eta_MC_g[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_g_cent"+cent[ibin-1]))->Clone((TString)("h_eta_g_cent"+cent[ibin]));
      //h_eta_data[ibin] = (TH1D*)closure_histos_PbPb_data->Get((TString)("h_eta_full_cent"+cent[ibin-1]));
      h_eta_MC_up[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_up_cent"+cent[ibin-1]));
      h_eta_MC_upbar[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_upbar_cent"+cent[ibin-1]));
      h_eta_MC_d[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_d_cent"+cent[ibin-1]));
      h_eta_MC_dbar[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_eta_full_dbar_cent"+cent[ibin-1]));

      h_reco_MC[ibin] = (TH1D*)closure_histos_PbPb_MC->Get((TString)("h_reco_corr_cent"+cent[ibin-1]));       
      //h_reco_data[ibin] = (TH1D*)closure_histos_PbPb_data->Get((TString)("h_reco_corr_cent"+cent[ibin-1]));
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
    //h_eta_data[ibin]->Scale(1./h_eta_data[ibin]->Integral());
    h_eta_MC_up[ibin]->Scale(1./h_eta_MC_up[ibin]->Integral());
    h_eta_MC_upbar[ibin]->Scale(1./h_eta_MC_upbar[ibin]->Integral());
    h_eta_MC_d[ibin]->Scale(1./h_eta_MC_d[ibin]->Integral());
    h_eta_MC_dbar[ibin]->Scale(1./h_eta_MC_dbar[ibin]->Integral());

    //normalizing jet spectra
    h_reco_MC[ibin]->Scale(1./h_reco_MC[ibin]->Integral());
    //h_reco_data[ibin]->Scale(1./h_reco_data[ibin]->Integral());
    h_reco_MC_q[ibin]->Scale(1./h_reco_MC_q[ibin]->Integral());
    h_reco_MC_g[ibin]->Scale(1./h_reco_MC_g[ibin]->Integral());
    h_reco_MC_up[ibin]->Scale(1./h_reco_MC_up[ibin]->Integral());
    h_reco_MC_upbar[ibin]->Scale(1./h_reco_MC_upbar[ibin]->Integral());
    h_reco_MC_d[ibin]->Scale(1./h_reco_MC_d[ibin]->Integral());
    h_reco_MC_dbar[ibin]->Scale(1./h_reco_MC_dbar[ibin]->Integral());
*/
  }

  TLine *tl1 = new TLine(-1.5,1.,1.5,1.);
  tl1->SetLineStyle(2);

  TLine *tl2 = new TLine(120.,0.,300.,0.);
  tl2->SetLineStyle(2);

  TLatex *tx1;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c_qg_eta = new TCanvas("c_qg_eta","",1500,300);
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

  TCanvas *c_reco = new TCanvas("c_reco","",1500,600);
  c_reco->Divide(5,2,0);
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
    tl2->Draw("same");

  }

  TCanvas *c_reco_gb = new TCanvas("c_reco_gb","",1500,600);
  c_reco_gb->Divide(5,2,0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

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

    if(ibin>0){
      TString tmp1 = cent_tag[ibin-1];
      tx1 = new TLatex(); tx1->SetTextSize(.12);
      tx1->DrawLatexNDC(0.15,0.85, tmp1);
    }

    if(ibin==0) c_reco_gb->cd(ibin+6);
    else c_reco_gb->cd(11-ibin);   
    h_reco_gluon_inc_ratio[ibin] = (TH1D*)h_reco_MC_g[ibin]->Clone("h_reco_gluon_inc_ratio");
    h_reco_gluon_inc_ratio[ibin]->Divide(h_reco_MC_upbar[ibin]);
    h_reco_up_inc_ratio[ibin] = (TH1D*)h_reco_MC_g[ibin]->Clone("h_reco_up_inc_ratio");
    h_reco_up_inc_ratio[ibin]->Divide(h_reco_MC_dbar[ibin]);
    h_reco_gluon_inc_ratio[ibin]->GetYaxis()->SetNdivisions(505);
    h_reco_gluon_inc_ratio[ibin]->GetYaxis()->SetTitleSize(0.05);
    h_reco_gluon_inc_ratio[ibin]->GetYaxis()->SetLabelSize(0.07);
    h_reco_gluon_inc_ratio[ibin]->GetXaxis()->CenterTitle();
    h_reco_gluon_inc_ratio[ibin]->GetXaxis()->SetNdivisions(505);
    h_reco_gluon_inc_ratio[ibin]->GetXaxis()->SetTitleSize(0.05);
    h_reco_gluon_inc_ratio[ibin]->GetXaxis()->SetLabelSize(0.07);
    h_reco_gluon_inc_ratio[ibin]->GetXaxis()->SetRangeUser(120.,400.);
    h_reco_gluon_inc_ratio[ibin]->GetYaxis()->SetRangeUser(10.,30.);
    h_reco_gluon_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_gluon_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_gluon_inc_ratio[ibin]->SetMarkerStyle(22); h_reco_gluon_inc_ratio[ibin]->Draw("e0 same");
    h_reco_up_inc_ratio[ibin]->SetLineColor(kBlue); h_reco_up_inc_ratio[ibin]->SetMarkerColor(kBlue);
    h_reco_up_inc_ratio[ibin]->SetMarkerStyle(23); h_reco_up_inc_ratio[ibin]->Draw("e0 same");
    tl2->Draw("same");

  }


}