#ifndef fitting_templates_h_
#define fitting_templates_h_

#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TProfile.h"
#include <iostream>

const int nCBins = 5;
const int ntrkbins = 4;

TH1D *h_chg_reco_g[nCBins][ntrkbins];
TH1D *h_chg_reco_up[nCBins][ntrkbins];
TH1D *h_chg_reco_down[nCBins][ntrkbins];
TH1D *h_chg_reco_gubdb[nCBins][ntrkbins];

TH1D *h_chg_reco_c[nCBins][ntrkbins];
TH1D *h_chg_reco_s[nCBins][ntrkbins];
TH1D *h_chg_reco_b[nCBins][ntrkbins];

TH1D *h_chg_ref_q[nCBins][ntrkbins];
TH1D *h_chg_ref_g[nCBins][ntrkbins];
TH1D *h_chg_ref_up[nCBins][ntrkbins];
TH1D *h_chg_ref_down[nCBins][ntrkbins];
TH1D *h_chg_ref_upb[nCBins][ntrkbins];
TH1D *h_chg_ref_downb[nCBins][ntrkbins];
TH1D *h_chg_ref_gubdb[nCBins][ntrkbins];

/*
class fitting_templates{
	
	public:
        double f_udg_reco(Double_t *x, Double_t *par);
	
	private:
        double xx;
        int bin;
        TH1D *h_chg_reco_g[nCBins][ntrkbins];
        TH1D *h_chg_reco_up[nCBins][ntrkbins];
        TH1D *h_chg_reco_down[nCBins][ntrkbins];
        TH1D *h_chg_reco_gubdb[nCBins][ntrkbins];
};
*/
/*
Double_t f_udg_ref_cent0_pt2(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[0][1]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[0][1]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[0][1]->GetBinContent(bin);
   return q+g; 	
}

Double_t f_udg_ref_cent1_pt2(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[1][1]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[1][1]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[1][1]->GetBinContent(bin);
   return q+g; 	
}

Double_t f_udg_ref_cent2_pt2(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[2][1]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[2][1]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[2][1]->GetBinContent(bin);
   return q+g; 	
}

Double_t f_udg_ref_cent3_pt2(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[3][1]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[3][1]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[3][1]->GetBinContent(bin);
   return q+g; 	
}

Double_t f_udg_ref_cent4_pt2(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[4][1]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[4][1]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[4][1]->GetBinContent(bin);
   return q+g; 	
}


Double_t f_udg_ref_cent0_pt4(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[0][2]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[0][2]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[0][2]->GetBinContent(bin);
   return q+g; 	
}

Double_t f_udg_ref_cent1_pt4(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[1][2]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[1][2]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[1][2]->GetBinContent(bin);
   return q+g; 	
}

Double_t f_udg_ref_cent2_pt4(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[2][2]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[2][2]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[2][2]->GetBinContent(bin);
   return q+g; 	
}

Double_t f_udg_ref_cent3_pt4(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[3][2]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[3][2]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[3][2]->GetBinContent(bin);
   return q+g; 	
}

Double_t f_udg_ref_cent4_pt4(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[4][2]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[4][2]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[4][2]->GetBinContent(bin);
   return q+g; 	
}


Double_t f_udg_ref_cent0_pt5(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[0][3]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[0][3]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[0][3]->GetBinContent(bin);
   return q+g; 	
}

Double_t f_udg_ref_cent1_pt5(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[1][3]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[1][3]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[1][3]->GetBinContent(bin);
   return q+g; 	
}

Double_t f_udg_ref_cent2_pt5(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[2][3]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[2][3]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[2][3]->GetBinContent(bin);
   return q+g; 	
}

Double_t f_udg_ref_cent3_pt5(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[3][3]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[3][3]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[3][3]->GetBinContent(bin);
   return q+g; 	
}

Double_t f_udg_ref_cent4_pt5(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[4][3]->FindBin(xx);
   Double_t q = (par[0])*h_chg_ref_q[4][3]->GetBinContent(bin);
   Double_t g = (1-par[0])*h_chg_ref_g[4][3]->GetBinContent(bin);
   return q+g; 	
}
*/

Double_t f_udg_ref_cent0_pt2(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_gubdb[0][1]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[0][1]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[0][1]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[0][1]->GetBinContent(bin);
   return up+dn+g; 	
}

Double_t f_udg_ref_cent1_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[1][1]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[1][1]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[1][1]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[1][1]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_ref_cent2_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[2][1]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[2][1]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[2][1]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[2][1]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_ref_cent3_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[3][1]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[3][1]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[3][1]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[3][1]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_ref_cent4_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[4][1]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[4][1]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[4][1]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[4][1]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_ref_cent0_pt4(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_gubdb[0][2]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[0][2]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[0][2]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[0][2]->GetBinContent(bin);
   return up+dn+g; 	
}

Double_t f_udg_ref_cent1_pt4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[1][2]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[1][2]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[1][2]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[1][2]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_ref_cent2_pt4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[2][2]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[2][2]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[2][2]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[2][2]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_ref_cent3_pt4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[3][2]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[3][2]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[3][2]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[3][2]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_ref_cent4_pt4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[4][2]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[4][2]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[4][2]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[4][2]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_ref_cent0_pt5(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_gubdb[0][3]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[0][3]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[0][3]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[0][3]->GetBinContent(bin);
   return up+dn+g; 	
}

Double_t f_udg_ref_cent1_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[1][3]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[1][3]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[1][3]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[1][3]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_ref_cent2_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[2][3]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[2][3]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[2][3]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[2][3]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_ref_cent3_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[3][3]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[3][3]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[3][3]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[3][3]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_ref_cent4_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_ref_g[4][3]->FindBin(xx);
   Double_t up = (par[0])*h_chg_ref_up[4][3]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_ref_down[4][3]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_ref_gubdb[4][3]->GetBinContent(bin);
   return up+dn+g; 
}







Double_t f_udg_reco_cent0_pt2(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_gubdb[0][1]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[0][1]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[0][1]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[0][1]->GetBinContent(bin);
   return up+dn+g; 	
}

Double_t f_udg_reco_cent1_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g[1][1]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[1][1]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[1][1]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[1][1]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_reco_cent2_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g[2][1]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[2][1]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[2][1]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[2][1]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_reco_cent3_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g[3][1]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[3][1]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[3][1]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[3][1]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_reco_cent4_pt2(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g[4][1]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[4][1]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[4][1]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[4][1]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_reco_cent0_pt4(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_gubdb[0][2]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[0][2]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[0][2]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[0][2]->GetBinContent(bin);
   return up+dn+g; 	
}

Double_t f_udg_reco_cent1_pt4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g[1][2]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[1][2]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[1][2]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[1][2]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_reco_cent2_pt4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g[2][2]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[2][2]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[2][2]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[2][2]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_reco_cent3_pt4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g[3][2]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[3][2]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[3][2]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[3][2]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_reco_cent4_pt4(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g[4][2]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[4][2]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[4][2]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[4][2]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_reco_cent0_pt5(Double_t *x, Double_t *par){
	
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_gubdb[0][3]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[0][3]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[0][3]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[0][3]->GetBinContent(bin);
   return up+dn+g; 	
}

Double_t f_udg_reco_cent1_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g[1][3]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[1][3]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[1][3]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[1][3]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_reco_cent2_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g[2][3]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[2][3]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[2][3]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[2][3]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_reco_cent3_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g[3][3]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[3][3]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[3][3]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[3][3]->GetBinContent(bin);
   return up+dn+g; 
}

Double_t f_udg_reco_cent4_pt5(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = h_chg_reco_g[4][3]->FindBin(xx);
   Double_t up = (par[0])*h_chg_reco_up[4][3]->GetBinContent(bin);
   Double_t dn = (par[1])*h_chg_reco_down[4][3]->GetBinContent(bin);
   Double_t g = (1-(par[0]+par[1]))*h_chg_reco_gubdb[4][3]->GetBinContent(bin);
   return up+dn+g; 
}

#endif
