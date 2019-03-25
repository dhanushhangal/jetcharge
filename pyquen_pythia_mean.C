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

#define nCBins 1
#define nptBins 3
#define ntrkbins 4
#define njtptbins 10
#define nrmsbins 8

const Double_t pt_low = 120.;
const Double_t pt_high = 600.;

Double_t jt_bin_bounds[9] = {120., 130., 140., 150., 160., 170., 180., 190., 200.};

char saythis[500];

using namespace std;

void pyquen_pythia_mean(){
    
    TFile *closure_histos_pp_pythia_gen = TFile::Open("P+H_jetchg_trkeffconst_no0trkjets_20190220.root");
    TFile *closure_histos_pp_herwig_gen = TFile::Open("PyquenWide_jetchg_trkeffconst_no0trkjets_20190215.root");
    TFile *closure_histos_PbPb_MC_pyquenColl = TFile::Open("PyquenColl_jetchg_trkeffconst_no0trkjets_20190215.root");

    TH1D *h_jetpt_pythia_up_ratio, *h_jetpt_pythia_down_ratio, *h_jetpt_pythia_g_ratio, *h_jetpt_pythia_oth_ratio;
    TH1D *h_jetpt_herwig_up_ratio, *h_jetpt_herwig_down_ratio, *h_jetpt_herwig_g_ratio, *h_jetpt_herwig_oth_ratio;
/*
    TH1D *h_jetpt_pythia = (TH1D*)closure_histos_pp_pythia_gen->Get("h_jetpt");
    TH1D *h_jetpt_pythia_up = (TH1D*)closure_histos_pp_pythia_gen->Get("h_jetpt_up");
    TH1D *h_jetpt_pythia_down = (TH1D*)closure_histos_pp_pythia_gen->Get("h_jetpt_down");
    TH1D *h_jetpt_pythia_g = (TH1D*)closure_histos_pp_pythia_gen->Get("h_jetpt_g");
    TH1D *h_jetpt_pythia_oth = (TH1D*)closure_histos_pp_pythia_gen->Get("h_jetpt_oth");
*/
    TH1D *h_gen_jetpt;
    TH1D *h_gen_jetpt_q;
    TH1D *h_gen_jetpt_g;
    TH1D *h_gen_jetpt_oth;

    TH1D *h_jetchg_pythia_rms[nrmsbins];
    TH1D *h_jetchg_pythia_g_rms[nrmsbins];
    TH1D *h_jetchg_pythia_up_rms[nrmsbins];
    TH1D *h_jetchg_pythia_down_rms[nrmsbins];
    TH1D *h_jetchg_pythia_oth_rms[nrmsbins];

    TH1D *h_jetchg_herwig_rms[nrmsbins];
    TH1D *h_jetchg_herwig_g_rms[nrmsbins];
    TH1D *h_jetchg_herwig_up_rms[nrmsbins];
    TH1D *h_jetchg_herwig_down_rms[nrmsbins];
    TH1D *h_jetchg_herwig_oth_rms[nrmsbins];

    TH1D *h_jetchg_pythia_rms_tot = new TH1D("h_jetchg_pythia_rms_tot","",8,120.,200.); 
    TH1D *h_jetchg_herwig_rms_tot = new TH1D("h_jetchg_herwig_rms_tot","",8,120.,200.);  
    TH1D *h_jetchg_pythia_g_rms_tot = new TH1D("h_jetchg_pythia_g_rms_tot","",8,120.,200.); 
    TH1D *h_jetchg_herwig_g_rms_tot = new TH1D("h_jetchg_herwig_g_rms_tot","",8,120.,200.);  
    TH1D *h_jetchg_pythia_up_rms_tot = new TH1D("h_jetchg_pythia_up_rms_tot","",8,120.,200.); 
    TH1D *h_jetchg_herwig_up_rms_tot = new TH1D("h_jetchg_herwig_up_rms_tot","",8,120.,200.); 
    TH1D *h_jetchg_pythia_down_rms_tot = new TH1D("h_jetchg_pythia_down_rms_tot","",8,120.,200.); 
    TH1D *h_jetchg_herwig_down_rms_tot = new TH1D("h_jetchg_herwig_down_rms_tot","",8,120.,200.); 
/*
    h_jetpt_pythia_up->Scale(1./h_jetpt_pythia->Integral());
    h_jetpt_pythia_down->Scale(1./h_jetpt_pythia->Integral());
    h_jetpt_pythia_g->Scale(1./h_jetpt_pythia->Integral());
    h_jetpt_pythia_oth->Scale(1./h_jetpt_pythia->Integral());
    h_jetpt_pythia->Scale(1./h_jetpt_pythia->Integral());
*/
    //TProfile *h_jetchg_herwig_ratio, *h_jetchg_herwig_up_ratio, *h_jetchg_herwig_down_ratio, *h_jetchg_herwig_g_ratio, *h_jetchg_herwig_oth_ratio;

    TH2D *h_jetchg_recopt_pythia = (TH2D*)closure_histos_pp_pythia_gen->Get("h_chg_refpt_cent0_trk0");
    TH2D *h_jetchg_recopt_pythia_up = (TH2D*)closure_histos_pp_pythia_gen->Get("h_chg_refpt_up_cent0_trk0");
    TH2D *h_jetchg_recopt_pythia_down = (TH2D*)closure_histos_pp_pythia_gen->Get("h_chg_refpt_d_cent0_trk0");
    TH2D *h_jetchg_recopt_pythia_g = (TH2D*)closure_histos_pp_pythia_gen->Get("h_chg_refpt_g_cent0_trk0");

    TProfile *h_jetchg_pythia = h_jetchg_recopt_pythia-> ProfileX("h_jetchg_pythia",h_jetchg_recopt_pythia->GetYaxis()->FindBin(-2.),h_jetchg_recopt_pythia->GetYaxis()->FindBin(2.),"");
    TProfile *h_jetchg_pythia_g = h_jetchg_recopt_pythia_g-> ProfileX("h_jetchg_pythia_g",h_jetchg_recopt_pythia->GetYaxis()->FindBin(-2.),h_jetchg_recopt_pythia->GetYaxis()->FindBin(2.),"");
    TProfile *h_jetchg_pythia_up = h_jetchg_recopt_pythia_up-> ProfileX("h_jetchg_pythia_up",h_jetchg_recopt_pythia->GetYaxis()->FindBin(-2.),h_jetchg_recopt_pythia->GetYaxis()->FindBin(2.),"");
    TProfile *h_jetchg_pythia_down = h_jetchg_recopt_pythia_down-> ProfileX("h_jetchg_pythia_down",h_jetchg_recopt_pythia->GetYaxis()->FindBin(-2.),h_jetchg_recopt_pythia->GetYaxis()->FindBin(2.),"");

    for(int i=0; i<nrmsbins; i++){
        h_jetchg_pythia_rms[i] = h_jetchg_recopt_pythia-> ProjectionY(Form("h_jetchg_pythia_rms_%d",i),h_jetchg_recopt_pythia->GetXaxis()->FindBin(jt_bin_bounds[i]),h_jetchg_recopt_pythia->GetXaxis()->FindBin(jt_bin_bounds[i+1]),"");
        h_jetchg_pythia_g_rms[i] = h_jetchg_recopt_pythia_g-> ProjectionY(Form("h_jetchg_pythia_g_rms_%d",i),h_jetchg_recopt_pythia->GetXaxis()->FindBin(jt_bin_bounds[i]),h_jetchg_recopt_pythia->GetXaxis()->FindBin(jt_bin_bounds[i+1]),"");
        h_jetchg_pythia_up_rms[i] = h_jetchg_recopt_pythia_up-> ProjectionY(Form("h_jetchg_pythia_up_rms_%d",i),h_jetchg_recopt_pythia->GetXaxis()->FindBin(jt_bin_bounds[i]),h_jetchg_recopt_pythia->GetXaxis()->FindBin(jt_bin_bounds[i+1]),"");
        h_jetchg_pythia_down_rms[i] = h_jetchg_recopt_pythia_down-> ProjectionY(Form("h_jetchg_pythia_down_rms_%d",i),h_jetchg_recopt_pythia->GetXaxis()->FindBin(jt_bin_bounds[i]),h_jetchg_recopt_pythia->GetXaxis()->FindBin(jt_bin_bounds[i+1]),"");

        h_jetchg_pythia_rms_tot->SetBinContent(i+1,h_jetchg_pythia_rms[i]->GetRMS());
        h_jetchg_pythia_rms_tot->SetBinError(i+1,h_jetchg_pythia_rms[i]->GetRMSError());
        h_jetchg_pythia_g_rms_tot->SetBinContent(i+1,h_jetchg_pythia_g_rms[i]->GetRMS());
        h_jetchg_pythia_g_rms_tot->SetBinError(i+1,h_jetchg_pythia_g_rms[i]->GetRMSError());
        h_jetchg_pythia_up_rms_tot->SetBinContent(i+1,h_jetchg_pythia_up_rms[i]->GetRMS());
        h_jetchg_pythia_up_rms_tot->SetBinError(i+1,h_jetchg_pythia_up_rms[i]->GetRMSError());
        h_jetchg_pythia_down_rms_tot->SetBinContent(i+1,h_jetchg_pythia_down_rms[i]->GetRMS());
        h_jetchg_pythia_down_rms_tot->SetBinError(i+1,h_jetchg_pythia_down_rms[i]->GetRMSError());

    }
/*
    TH1D *h_jetpt_herwig = (TH1D*)closure_histos_pp_herwig_gen->Get("h_jetpt");
    TH1D *h_jetpt_herwig_up = (TH1D*)closure_histos_pp_herwig_gen->Get("h_jetpt_up");
    TH1D *h_jetpt_herwig_down = (TH1D*)closure_histos_pp_herwig_gen->Get("h_jetpt_down");
    TH1D *h_jetpt_herwig_g = (TH1D*)closure_histos_pp_herwig_gen->Get("h_jetpt_g");
    TH1D *h_jetpt_herwig_oth = (TH1D*)closure_histos_pp_herwig_gen->Get("h_jetpt_oth");

    h_jetpt_herwig_up->Scale(1./h_jetpt_herwig->Integral());
    h_jetpt_herwig_down->Scale(1./h_jetpt_herwig->Integral());
    h_jetpt_herwig_g->Scale(1./h_jetpt_herwig->Integral());
    h_jetpt_herwig_oth->Scale(1./h_jetpt_herwig->Integral());
    h_jetpt_herwig->Scale(1./h_jetpt_herwig->Integral());
*/
    TH2D *h_jetchg_recopt_herwig = (TH2D*)closure_histos_pp_herwig_gen->Get("h_chg_corrpt_cent0_trk0");
    TH2D *h_jetchg_recopt_herwig_up = (TH2D*)closure_histos_pp_herwig_gen->Get("h_chg_corrpt_up_cent0_trk0");
    TH2D *h_jetchg_recopt_herwig_down = (TH2D*)closure_histos_pp_herwig_gen->Get("h_chg_corrpt_d_cent0_trk0");
    TH2D *h_jetchg_recopt_herwig_g = (TH2D*)closure_histos_pp_herwig_gen->Get("h_chg_corrpt_g_cent0_trk0");

    h_gen_jetpt = (TH1D*)closure_histos_pp_pythia_gen->Get("h_gen_full_cent0");
    h_gen_jetpt_q = (TH1D*)closure_histos_pp_pythia_gen->Get("h_gen_full_up_cent0");
    h_gen_jetpt_q-> Add((TH1D*)closure_histos_pp_pythia_gen->Get("h_gen_full_d_cent0"));
    h_gen_jetpt_g = (TH1D*)closure_histos_pp_pythia_gen->Get("h_gen_full_g_cent0");
    h_gen_jetpt_oth = (TH1D*)closure_histos_pp_pythia_gen->Get("h_gen_full_upbar_cent0");
    h_gen_jetpt_oth->Add((TH1D*)closure_histos_pp_pythia_gen->Get("h_gen_full_dbar_cent0")); 
    h_gen_jetpt_oth->Add((TH1D*)closure_histos_pp_pythia_gen->Get("h_gen_full_c_cent0")); 
    h_gen_jetpt_oth->Add((TH1D*)closure_histos_pp_pythia_gen->Get("h_gen_full_s_cent0")); 
    h_gen_jetpt_oth->Add((TH1D*)closure_histos_pp_pythia_gen->Get("h_gen_full_b_cent0")); 

    TProfile *h_jetchg_herwig = h_jetchg_recopt_herwig-> ProfileX("h_jetchg_herwig",h_jetchg_recopt_herwig->GetYaxis()->FindBin(-2.),h_jetchg_recopt_herwig->GetYaxis()->FindBin(2.),"");
    TProfile *h_jetchg_herwig_g = h_jetchg_recopt_herwig_g-> ProfileX("h_jetchg_herwig_g",h_jetchg_recopt_herwig->GetYaxis()->FindBin(-2.),h_jetchg_recopt_herwig->GetYaxis()->FindBin(2.),"");
    TProfile *h_jetchg_herwig_up = h_jetchg_recopt_herwig_up-> ProfileX("h_jetchg_herwig_up",h_jetchg_recopt_herwig->GetYaxis()->FindBin(-2.),h_jetchg_recopt_herwig->GetYaxis()->FindBin(2.),"");
    TProfile *h_jetchg_herwig_down = h_jetchg_recopt_herwig_down-> ProfileX("h_jetchg_herwig_down",h_jetchg_recopt_herwig->GetYaxis()->FindBin(-2.),h_jetchg_recopt_herwig->GetYaxis()->FindBin(2.),"");

    for(int i=0; i<nrmsbins; i++){
        h_jetchg_herwig_rms[i] = h_jetchg_recopt_herwig-> ProjectionY(Form("h_jetchg_herwig_rms_%d",i),h_jetchg_recopt_herwig->GetXaxis()->FindBin(jt_bin_bounds[i]),h_jetchg_recopt_herwig->GetXaxis()->FindBin(jt_bin_bounds[i+1]),"");
        h_jetchg_herwig_g_rms[i] = h_jetchg_recopt_herwig_g-> ProjectionY(Form("h_jetchg_herwig_g_rms_%d",i),h_jetchg_recopt_herwig->GetXaxis()->FindBin(jt_bin_bounds[i]),h_jetchg_recopt_herwig->GetXaxis()->FindBin(jt_bin_bounds[i+1]),"");
        h_jetchg_herwig_up_rms[i] = h_jetchg_recopt_herwig_up-> ProjectionY(Form("h_jetchg_herwig_up_rms_%d",i),h_jetchg_recopt_herwig->GetXaxis()->FindBin(jt_bin_bounds[i]),h_jetchg_recopt_herwig->GetXaxis()->FindBin(jt_bin_bounds[i+1]),"");
        h_jetchg_herwig_down_rms[i] = h_jetchg_recopt_herwig_down-> ProjectionY(Form("h_jetchg_herwig_down_rms_%d",i),h_jetchg_recopt_herwig->GetXaxis()->FindBin(jt_bin_bounds[i]),h_jetchg_recopt_herwig->GetXaxis()->FindBin(jt_bin_bounds[i+1]),"");

        h_jetchg_herwig_rms_tot->SetBinContent(i+1,h_jetchg_herwig_rms[i]->GetRMS());
        h_jetchg_herwig_rms_tot->SetBinError(i+1,h_jetchg_herwig_rms[i]->GetRMSError());
        h_jetchg_herwig_g_rms_tot->SetBinContent(i+1,h_jetchg_herwig_g_rms[i]->GetRMS());
        h_jetchg_herwig_g_rms_tot->SetBinError(i+1,h_jetchg_herwig_g_rms[i]->GetRMSError());
        h_jetchg_herwig_up_rms_tot->SetBinContent(i+1,h_jetchg_herwig_up_rms[i]->GetRMS());
        h_jetchg_herwig_up_rms_tot->SetBinError(i+1,h_jetchg_herwig_up_rms[i]->GetRMSError());
        h_jetchg_herwig_down_rms_tot->SetBinContent(i+1,h_jetchg_herwig_down_rms[i]->GetRMS());
        h_jetchg_herwig_down_rms_tot->SetBinError(i+1,h_jetchg_herwig_down_rms[i]->GetRMSError());
    }

////jet charge

    TCanvas *c_jetchg = new TCanvas("c_jetchg","c_jetchg",600,1200);
    c_jetchg->Divide(1,2,0);
    gStyle->SetOptStat(0);
    //gStyle->SetOptFit(1);

    c_jetchg->cd(1);   

    h_jetchg_pythia_g->GetXaxis()->SetTitle("jet pT");
    h_jetchg_pythia_g->GetXaxis()->SetRangeUser(120.,200.);
    h_jetchg_pythia_g->GetYaxis()->SetRangeUser(-0.25,0.35);
    h_jetchg_pythia_g->GetYaxis()->SetTitle("<jet charge>");
    h_jetchg_pythia_g->GetYaxis()->CenterTitle();
    h_jetchg_pythia->SetLineColor(kBlack); h_jetchg_pythia->SetMarkerColor(kBlack);
    //h_jetchg_pythia->SetMarkerStyle(20); h_jetchg_pythia->Draw("e1 same");
    h_jetchg_pythia_g->SetLineColor(kRed-3); h_jetchg_pythia_g->SetMarkerColor(kRed-3); h_jetchg_pythia_g->SetFillColor(kRed-3);
    h_jetchg_pythia_g->SetMarkerStyle(20); h_jetchg_pythia_g->Draw("e1 same");
    h_jetchg_pythia_up->SetLineColor(kAzure+2); h_jetchg_pythia_up->SetMarkerColor(kAzure+2); h_jetchg_pythia_up->SetFillColor(kAzure+2);
    h_jetchg_pythia_up->SetMarkerStyle(22); h_jetchg_pythia_up->Draw("e1 same");
    h_jetchg_pythia_down->SetLineColor(kGreen-2); h_jetchg_pythia_down->SetMarkerColor(kGreen-2); h_jetchg_pythia_down->SetFillColor(kGreen-2);
    h_jetchg_pythia_down->SetMarkerStyle(23); h_jetchg_pythia_down->Draw("e1 same");
    //h_jetchg_pythia_oth->SetLineColor(kOrange+1); h_jetchg_pythia_oth->SetMarkerColor(kOrange+1); h_jetchg_pythia_oth->SetFillColor(kOrange+1);
    //h_jetchg_pythia_oth->SetMarkerStyle(20); h_jetchg_pythia_oth->Draw("e1 same");

    TLegend *legendq = new TLegend(0.25,0.75,0.75,0.9);
    legendq ->SetLineColor(kWhite);
    //legendq ->AddEntry((TObject*)0, "PYTHIA", "");
    //legendq ->AddEntry(h_jetchg_pythia, "inclusive", "lep");
    legendq ->AddEntry(h_jetchg_pythia_g, "gluon", "lep");
    legendq ->AddEntry(h_jetchg_pythia_up, "up quark", "lep");
    legendq ->AddEntry(h_jetchg_pythia_down, "down quark", "lep");
    legendq ->Draw("same");
/*
    h_jetchg_herwig_g->GetXaxis()->SetTitle("jet pT");
    h_jetchg_herwig_g->GetXaxis()->SetRangeUser(120.,200.);
    h_jetchg_herwig_g->GetYaxis()->SetRangeUser(-0.2,0.3);
    h_jetchg_herwig_g->SetLineColor(kRed-3); h_jetchg_herwig_g->SetMarkerColor(kRed-3); h_jetchg_herwig_g->SetFillColor(kRed-3);
    h_jetchg_herwig_g->SetMarkerStyle(20); h_jetchg_herwig_g->Draw("e1 same");
    h_jetchg_herwig_up->SetLineColor(kAzure+2); h_jetchg_herwig_up->SetMarkerColor(kAzure+2); h_jetchg_herwig_up->SetFillColor(kAzure+2);
    h_jetchg_herwig_up->SetMarkerStyle(22); h_jetchg_herwig_up->Draw("e1 same");
    h_jetchg_herwig_down->SetLineColor(kGreen-2); h_jetchg_herwig_down->SetMarkerColor(kGreen-2); h_jetchg_herwig_down->SetFillColor(kGreen-2);
    h_jetchg_herwig_down->SetMarkerStyle(23); h_jetchg_herwig_down->Draw("e1 same");
    h_jetchg_herwig_oth->SetLineColor(kOrange+1); h_jetchg_herwig_oth->SetMarkerColor(kOrange+1); h_jetchg_herwig_oth->SetFillColor(kOrange+1);
    //h_jetchg_herwig_oth->SetMarkerStyle(4); h_jetchg_herwig_oth->Draw("e1 same");
    h_jetchg_herwig->SetLineColor(kBlack); h_jetchg_herwig->SetMarkerColor(kBlack);
    //h_jetchg_herwig->SetMarkerStyle(4); h_jetchg_herwig->Draw("e1 same");

    TLegend *legendh2x = new TLegend(0.25,0.75,0.75,0.9);
    legendh2x ->SetLineColor(kWhite);
    //legendh2x ->AddEntry((TObject*)0, "HERWIG", "");
    //legendh2x ->AddEntry(h_jetchg_herwig, "inclusive", "lep");
    legendh2x ->AddEntry(h_jetchg_herwig_g, "gluon", "lep");
    legendh2x ->AddEntry(h_jetchg_herwig_up, "up quark", "lep");
    legendh2x ->AddEntry(h_jetchg_herwig_down, "down quark", "lep");
    legendh2x ->Draw("same");
*/
    c_jetchg->cd(2); 
    h_jetchg_pythia_g_rms_tot->GetYaxis()->SetNdivisions(505);
    h_jetchg_pythia_g_rms_tot->GetYaxis()->SetTitleSize(0.05);
    h_jetchg_pythia_g_rms_tot->GetXaxis()->SetRangeUser(120.,200.);
    h_jetchg_pythia_g_rms_tot->GetYaxis()->SetRangeUser(0.2,0.65);
    h_jetchg_pythia_g_rms_tot->GetYaxis()->SetTitle("RMS");
    h_jetchg_pythia_g_rms_tot->GetXaxis()->SetTitle("jet p_{T}");
    h_jetchg_herwig_g_rms_tot->SetLineColor(kBlack); h_jetchg_herwig_rms_tot->SetMarkerColor(kBlack);
    //h_jetchg_herwig_rms_tot->SetMarkerStyle(4); h_jetchg_herwig_rms_tot->Draw("e1 same");
    h_jetchg_pythia_rms_tot->SetLineColor(kBlack); h_jetchg_pythia_rms_tot->SetMarkerColor(kBlack);
    //h_jetchg_pythia_rms_tot->SetMarkerStyle(20); h_jetchg_pythia_rms_tot->Draw("e0 same");

    h_jetchg_herwig_g_rms_tot->SetLineColor(kRed-3); h_jetchg_herwig_g_rms_tot->SetMarkerColor(kRed-3);
    //h_jetchg_herwig_g_rms_tot->SetMarkerStyle(20); h_jetchg_herwig_g_rms_tot->Draw("e1 same");
    h_jetchg_pythia_g_rms_tot->SetLineColor(kRed-3); h_jetchg_pythia_g_rms_tot->SetMarkerColor(kRed-3);
    h_jetchg_pythia_g_rms_tot->SetMarkerStyle(20); h_jetchg_pythia_g_rms_tot->Draw("e1 same");

    h_jetchg_herwig_up_rms_tot->SetLineColor(kAzure+2); h_jetchg_herwig_up_rms_tot->SetMarkerColor(kAzure+2);
    //h_jetchg_herwig_up_rms_tot->SetMarkerStyle(22); h_jetchg_herwig_up_rms_tot->Draw("e1 same");
    h_jetchg_pythia_up_rms_tot->SetLineColor(kAzure+2); h_jetchg_pythia_up_rms_tot->SetMarkerColor(kAzure+2);
    h_jetchg_pythia_up_rms_tot->SetMarkerStyle(22); h_jetchg_pythia_up_rms_tot->Draw("e1 same");

    h_jetchg_herwig_down_rms_tot->SetLineColor(kGreen-2); h_jetchg_herwig_down_rms_tot->SetMarkerColor(kGreen-2);
    //h_jetchg_herwig_down_rms_tot->SetMarkerStyle(23); h_jetchg_herwig_down_rms_tot->Draw("e1 same");
    h_jetchg_pythia_down_rms_tot->SetLineColor(kGreen-2); h_jetchg_pythia_down_rms_tot->SetMarkerColor(kGreen-2);
    h_jetchg_pythia_down_rms_tot->SetMarkerStyle(23); h_jetchg_pythia_down_rms_tot->Draw("e1 same");
/*
    h_jetchg_herwig_ratio->SetLineColor(kBlack); h_jetchg_herwig_ratio->SetMarkerColor(kBlack);
    h_jetchg_herwig_ratio->SetMarkerStyle(25); h_jetchg_herwig_ratio->Draw("e0 same");
    h_jetchg_herwig_up_ratio->SetLineColor(kAzure+2); h_jetchg_herwig_up_ratio->SetMarkerColor(kAzure+2);
    h_jetchg_herwig_up_ratio->SetMarkerStyle(25); h_jetchg_herwig_up_ratio->Draw("e0 same");
    h_jetchg_herwig_down_ratio->SetLineColor(kGreen-2); h_jetchg_herwig_down_ratio->SetMarkerColor(kGreen-2);
    h_jetchg_herwig_down_ratio->SetMarkerStyle(25); h_jetchg_herwig_down_ratio->Draw("e0 same");

    TLegend *legendhs = new TLegend(0.25,0.75,0.75,0.9);
    legendhs ->SetLineColor(kWhite);
    legendhs ->AddEntry((TObject*)0, "HERWIG", "");
    legendhs ->AddEntry(h_jetchg_herwig_g_ratio, "inclusive", "lep");
    legendhs ->AddEntry(h_jetchg_herwig_g_ratio, "gluon", "lep");
    legendhs ->AddEntry(h_jetchg_herwig_up_ratio, "up quark", "lep");
    legendhs ->AddEntry(h_jetchg_herwig_down_ratio, "down quark", "lep");
    legendhs ->Draw("same");
*/
/*
    c_jetchg->cd(2);   
    h_jetchg_herwig_ratio = (TProfile*)h_jetchg_pythia->Clone("h_reco_inc_ratio");
    h_jetchg_herwig_ratio->Divide(h_jetchg_herwig);
    h_jetchg_herwig_g_ratio = (TProfile*)h_jetchg_pythia_g->Clone("h_reco_gluon_inc_ratio");
    h_jetchg_herwig_g_ratio->Divide(h_jetchg_herwig_g);
    h_jetchg_herwig_up_ratio = (TProfile*)h_jetchg_pythia_up->Clone("h_reco_up_inc_ratio");
    h_jetchg_herwig_up_ratio->Divide(h_jetchg_herwig_up);
    h_jetchg_herwig_down_ratio = (TProfile*)h_jetchg_pythia_down->Clone("h_reco_d_inc_ratio");
    h_jetchg_herwig_down_ratio->Divide(h_jetchg_herwig_down);


    h_jetchg_herwig_g_ratio->GetYaxis()->SetNdivisions(505);
    h_jetchg_herwig_g_ratio->GetYaxis()->SetTitleSize(0.05);
    h_jetchg_herwig_g_ratio->GetXaxis()->SetRangeUser(120.,400.);
    h_jetchg_herwig_g_ratio->GetYaxis()->SetRangeUser(0.,0.8);
    h_jetchg_herwig_g_ratio->GetYaxis()->SetTitle("Pythia / Herwig");
    h_jetchg_herwig_g_ratio->GetXaxis()->SetTitle("jet p_{T}");
    h_jetchg_herwig_g_ratio->SetLineColor(kRed-3); h_jetchg_herwig_g_ratio->SetMarkerColor(kRed-3);
    h_jetchg_herwig_g_ratio->SetMarkerStyle(25); h_jetchg_herwig_g_ratio->Draw("e0 same");
    h_jetchg_herwig_ratio->SetLineColor(kBlack); h_jetchg_herwig_ratio->SetMarkerColor(kBlack);
    h_jetchg_herwig_ratio->SetMarkerStyle(25); h_jetchg_herwig_ratio->Draw("e0 same");
    h_jetchg_herwig_up_ratio->SetLineColor(kAzure+2); h_jetchg_herwig_up_ratio->SetMarkerColor(kAzure+2);
    h_jetchg_herwig_up_ratio->SetMarkerStyle(25); h_jetchg_herwig_up_ratio->Draw("e0 same");
    h_jetchg_herwig_down_ratio->SetLineColor(kGreen-2); h_jetchg_herwig_down_ratio->SetMarkerColor(kGreen-2);
    h_jetchg_herwig_down_ratio->SetMarkerStyle(25); h_jetchg_herwig_down_ratio->Draw("e0 same");

    TLegend *legendhs = new TLegend(0.25,0.75,0.75,0.9);
    legendhs ->SetLineColor(kWhite);
    legendhs ->AddEntry((TObject*)0, "HERWIG", "");
    legendhs ->AddEntry(h_jetchg_herwig_g_ratio, "inclusive", "lep");
    legendhs ->AddEntry(h_jetchg_herwig_g_ratio, "gluon", "lep");
    legendhs ->AddEntry(h_jetchg_herwig_up_ratio, "up quark", "lep");
    legendhs ->AddEntry(h_jetchg_herwig_down_ratio, "down quark", "lep");
    legendhs ->Draw("same");
*/

    TCanvas *c_jetpt = new TCanvas("c_jetpt","c_jetpt",600,600);
    gStyle->SetOptStat(0);

    c_jetpt->cd(1);  
    h_gen_jetpt_q->Divide(h_gen_jetpt); 
    h_gen_jetpt_g->Divide(h_gen_jetpt); 
    h_gen_jetpt_oth->Divide(h_gen_jetpt); 

    h_gen_jetpt_q->SetFillColor(kBlue-2); 
    h_gen_jetpt_g->SetFillColor(kGreen-2); 
    h_gen_jetpt_oth->SetFillColor(kRed-3); 

    THStack *h_flav_pt = new THStack("h_flav_pt","");
    h_flav_pt->Add(h_gen_jetpt_q);    
    h_flav_pt->Add(h_gen_jetpt_g);    
    h_flav_pt->Add(h_gen_jetpt_oth);    

    h_flav_pt->Draw("hist");
}