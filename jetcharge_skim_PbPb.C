#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TF1.h"

#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>
#include "TStopwatch.h"
#include "TEnv.h"

#include "assert.h"
#include <fstream>
#include "TMath.h"
#include <vector>

#include "jffcorr_lib/nCScorr.h"
#include "xiaoTrkCorr/xiaoTrkCorr.h"
#include "TrkCorr_July22_Iterative_pp_eta2p4/getTrkCorr.h"
//#include "TrkCorr_Jun7_Iterative_PbPb_etaLT2p4/getTrkCorr.h"

using namespace std;

float trkPtCut=2.;
const double trketamaxcut = 2.4;
const int nCbins = 4;

int mycbin;

char saythis[500];

float CBins[nCbins+1] = {0, 20, 60, 100, 200};

enum enum_dataset_types {e_Data2015,e_Data_pp,e_HydJet15,e_HydJet30,e_HydJet50, e_HydJet80, e_HydJet120,e_HydJet170,e_HydJet220,e_HydJet280, e_HydJet370,e_Pythia15,e_Pythia30,e_Pythia50, e_Pythia80, e_Pythia120,e_Pythia170,e_Pythia220,e_Pythia280, e_Pythia370, e_n_dataset_types};
TString dataset_type_strs[e_n_dataset_types] = {"Data2015","Data_pp","HydJet15","HydJet30","HydJet50","HydJet80", "HydJet120", "HydJet170","HydJet220","HydJet280","HydJet370","Pythia15","Pythia30","Pythia50","Pythia80", "Pythia120", "Pythia170","Pythia220","Pythia280","Pythia370"};

int dataset_pthats[e_n_dataset_types+1] = {0,0,15,30,50,80,120,170,220,280,370,15,30,50,80,120,170,220,280,370,999};

///pp MC Pythia6 pthat weights
double pthatbins[9] = {50,80,100,120,170,220,280,370,9999};
double xsecs[9] = {3.778E-03, 4.412E-04, 1.511E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
double pthatEntries[8] ={766451, 292063, 91200, 468748, 447938, 259208, 234447, 50972};

///PbPb MC Pythia6+Hydjet pthat weights
double pthatbins_PH[9] = {15,30,50,80,120,170,220,280,9999};
double xsecs_PH[9] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 0};
double pthatEntries_PH[8] ={0, 0, 0, 2571566, 2.85082e+06, 2.68057e+06, 2.89138e+06, 950344};

double kappa[5] = {0.1,0.3,0.5,0.7,0.9};

int dataset_type_code = -999;

bool passTrackCuts(bool ispp, bool useTightCuts, float trkPt, float trkEta, bool highPurity, float trkChi2, int trkNdof, int trkNlayer, int trkNHit, float pfHcal, float pfEcal, float trkDzSigma, float trkDxySigma){
	if(abs(trkEta)>=trketamaxcut) return false;
	if(!highPurity) return false;
    if(useTightCuts && (TMath::Abs(trkDzSigma)>=3.0 ||TMath::Abs(trkDxySigma)>=3.0)) return false;
	if(!ispp && (float)trkChi2/(float)trkNdof/(float)trkNlayer > 0.15) return false;
	if(!ispp && trkNHit<11 && trkPt > 0.7) return false;

	float Et = (pfHcal+pfEcal)/TMath::CosH(trkEta);
	if(!(trkPt<20 || Et > 0.5*trkPt)) return false;
	if(trkPt<=trkPtCut || trkPt > 400) return false;

	return true;
}

bool passGenTrackCuts(float trkPt, float trkEta, int chg){
    if(fabs(trkEta)>=trketamaxcut) return false;
    if(trkPt<=trkPtCut) return false;
    if(chg==0) return false;
    
    return true;
}

//cymbal tune centrality reweighting
double fcent_cymbal(double centrality, TF1* fcent1){ 
  return (centrality < 194) ? fcent1->Eval(centrality) : 1;
}

//arg 1 = which data set, arg 2 = output file number
void jetcharge_skim_PbPb(bool doCrab=0, int jobID=0, int endfile = 3203, int dataset_type_code = 2, int output_file_num = 1)
{
	bool is_pp = false;
	bool do_mixing = false;
    bool do_sube0 = true;

    nCScorr *corrpt = new nCScorr(is_pp,1);

    //TrkCorr* trkCorr = new TrkCorr("TrkCorr_July22_Iterative_pp_eta2p4/");
    //TrkCorr* trkCorr = new TrkCorr("TrkCorr_Jun7_Iterative_PbPb_etaLT2p4/");
    xiaoTrkCorr* trkCorr = new xiaoTrkCorr("xiaoTrkCorr/eta_symmetry_cymbalCorr_FineBin.root");

	bool is_data = false;

	if(dataset_type_code == 0 || dataset_type_code == 1) is_data = true;

	cout << "dataset code: " << dataset_type_code << endl;
	
	bool do_PbPb=1, do_pp_tracking=0;

	int radius = 4;

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

    TF1 *f_res_gauss_lowpt = new TF1("f_res_gauss_lowpt", "exp(-((x-1.)^2)/(2*([0]^2)))", 0., 2.);
    f_res_gauss_lowpt->FixParameter(0,0.0255);
    TF1 *f_res_gauss_highpt = new TF1("f_res_gauss_highpt", "exp(-((x-1.)^2)/(2*([0]^2)))", 0., 2.);
    f_res_gauss_highpt->FixParameter(0,0.015);

    //Pythia8 vz reweighting
    TF1 *fvz_pythia8 = new TF1("fvz_rw","pol4",-15.,15.); 
    fvz_pythia8->SetParameters(9.22005e-01,-1.42259e-02,3.47708e-03,-2.73684e-05,-4.09970e-06);

    //cymbal
    TF1 *fvz_cymbal = new TF1("fvz","pol6",-15,15); 
    fvz_cymbal->SetParameters(1.18472,-0.132675,0.00857998,-0.000326085,-1.48786e-06,4.68665e-07,-7.32942e-09 );

	TEnv *gEnv = new TEnv("HT_Ana_Env");
    gEnv->SetValue("TFile.AsyncPrefetching", 1);

	if(dataset_type_code== 1 || dataset_type_code > 10){do_PbPb = 0;   do_pp_tracking = 1;}

	cout<<"do_PbPb = "<<do_PbPb<<endl;
	cout<<"do_pp_tracking = "<<do_pp_tracking<<endl;

	TTree *inp_tree;
	TTree *inp_tree2;
	TTree *inp_tree3;
	TTree *inp_tree4;
	TTree *inp_tree5;
	TTree *inp_tree6;
	TTree *inp_tree7;
	TTree *inp_tree8;
	TTree *inp_tree_CS;

	string in_file_name;

	if(doCrab){
		in_file_name = Form("job_input_file_list_%d.txt",jobID);
	}
	else if(is_data&&!do_PbPb){
		in_file_name = "ppdata_filelist.txt";
	}else if(is_data&&do_PbPb){
		in_file_name = "PbPb2015_data_Marta_csidfix.txt";
	}else if(dataset_type_code > 10){
		in_file_name = "ppMC_pthat80_filelist.txt";
	}else if(dataset_type_code > 1&&dataset_type_code <11){
		in_file_name = "Pythia6Hydjet_PbPbMC_fullPthat_latestCS.txt";
	}else{
		cerr<<"need to set up to run on that sample..."<<endl;
	}

	cout << "trying a filelist named "<< in_file_name << endl;

	//MC
	TString output_file_base = "./";

	output_file_base +=dataset_type_strs[dataset_type_code];

	TString output_file_extension = "";   
	//output_file_extension += output_file_num;   
	output_file_extension += ".root";
	TFile *output_file = new TFile((TString) (output_file_base+output_file_extension), "RECREATE");
	TTree *ftree = new TTree("unzipMixTree", "");    

    TH1D *h_hibin = new TH1D("h_hibin","",200,0.,200.);
    TH1D *h_pthat = new TH1D("h_pthat","",50,0.,500.);
    TH1D *h_vz = new TH1D("h_vz","",30.,-15.,15.);

    TH1D *h_trk_pt_reco_full[nCbins];
    TH1D *h_trk_pt_reco[nCbins];
    TH1D *h_trk_pt2_reco[nCbins];
    TH1D *h_trk_pt4_reco[nCbins];
    TH1D *h_trk_pt5_reco[nCbins];
    TH1D *h_trk_pt_reco_refl[nCbins];
    TH1D *h_trk_pt2_reco_refl[nCbins];
    TH1D *h_trk_pt4_reco_refl[nCbins];
    TH1D *h_trk_pt5_reco_refl[nCbins];

    TH1D *h_trk_pt_reco_pos[nCbins];
    TH1D *h_trk_pt_reco_neg[nCbins];
    TH1D *h_trk_eta_reco[nCbins];

    TH1D *h_trk_ptk1_reco[nCbins];
    TH1D *h_trk_ptk1p5_reco[nCbins];
    TH1D *h_trk_ptk2_reco[nCbins];
    TH1D *h_trk_ptk2p5_reco[nCbins];

    TH1D *h_trk_ptk1_reco_refl[nCbins];
    TH1D *h_trk_ptk1p5_reco_refl[nCbins];
    TH1D *h_trk_ptk2_reco_refl[nCbins];
    TH1D *h_trk_ptk2p5_reco_refl[nCbins];

    TH1D *h_trk_pt_gen_full[nCbins];
    TH1D *h_trk_pt_gen[nCbins];
    TH1D *h_trk_pt4_gen[nCbins];
    TH1D *h_trk_pt5_gen[nCbins];
    TH1D *h_trk_pt_gen_pos[nCbins];
    TH1D *h_trk_pt_gen_neg[nCbins];
    TH1D *h_trk_eta_gen[nCbins];


    for(int ibin=0; ibin<nCbins; ibin++){
        
        sprintf(saythis,"h_trk_pt_reco_full%d",ibin);
        h_trk_pt_reco_full[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_pt_reco_%d",ibin);
        h_trk_pt_reco[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_pt2_reco_%d",ibin);
        h_trk_pt2_reco[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_pt4_reco_%d",ibin);
        h_trk_pt4_reco[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_pt5_reco_%d",ibin);
        h_trk_pt5_reco[ibin] = new TH1D(saythis,"",200,0.,20.);

        sprintf(saythis,"h_trk_pt_reco_refl_%d",ibin);
        h_trk_pt_reco_refl[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_pt2_reco_refl%d",ibin);
        h_trk_pt2_reco_refl[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_pt4_reco_refl%d",ibin);
        h_trk_pt4_reco_refl[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_pt5_reco_refl%d",ibin);
        h_trk_pt5_reco_refl[ibin] = new TH1D(saythis,"",200,0.,20.);

        sprintf(saythis,"h_trk_pt_reco_pos%d",ibin); 
        h_trk_pt_reco_pos[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_pt_reco_neg%d",ibin);
        h_trk_pt_reco_neg[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_eta_reco%d",ibin);
        h_trk_eta_reco[ibin] = new TH1D(saythis,"",50,-2.5,2.5);

        sprintf(saythis,"h_trk_ptk1_reco_%d",ibin);
        h_trk_ptk1_reco[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_ptk1p5_reco_%d",ibin);
        h_trk_ptk1p5_reco[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_ptk2_reco_%d",ibin);
        h_trk_ptk2_reco[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_ptk2p5_reco_%d",ibin);
        h_trk_ptk2p5_reco[ibin] = new TH1D(saythis,"",200,0.,20.);

        sprintf(saythis,"h_trk_ptk1_reco_refl_%d",ibin);
        h_trk_ptk1_reco_refl[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_ptk1p5_reco_refl_%d",ibin);
        h_trk_ptk1p5_reco_refl[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_ptk2_reco_refl_%d",ibin);
        h_trk_ptk2_reco_refl[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_ptk2p5_reco_refl_%d",ibin);
        h_trk_ptk2p5_reco_refl[ibin] = new TH1D(saythis,"",200,0.,20.);

        sprintf(saythis,"h_trk_pt_gen_full_%d",ibin); 
        h_trk_pt_gen_full[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_pt_gen_%d",ibin);
        h_trk_pt_gen[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_pt4_gen_%d",ibin);
        h_trk_pt4_gen[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_pt5_gen__%d",ibin);
        h_trk_pt5_gen[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_pt_gen_pos_%d",ibin);
        h_trk_pt_gen_pos[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_pt_gen_neg_%d",ibin);
        h_trk_pt_gen_neg[ibin] = new TH1D(saythis,"",200,0.,20.);
        sprintf(saythis,"h_trk_eta_gen_%d",ibin);
        h_trk_eta_gen[ibin] = new TH1D(saythis,"",50,-2.5,2.5);
       
        h_trk_pt_reco_full[ibin]->Sumw2();
        h_trk_pt_reco[ibin]->Sumw2();
        h_trk_pt2_reco[ibin]->Sumw2();
        h_trk_pt4_reco[ibin]->Sumw2();
        h_trk_pt5_reco[ibin]->Sumw2();
        h_trk_pt_reco_refl[ibin]->Sumw2();
        h_trk_pt2_reco_refl[ibin]->Sumw2();
        h_trk_pt4_reco_refl[ibin]->Sumw2();
        h_trk_pt5_reco_refl[ibin]->Sumw2();
        h_trk_pt_reco_pos[ibin]->Sumw2();
        h_trk_pt_reco_neg[ibin]->Sumw2();
        h_trk_eta_reco[ibin]->Sumw2();

        h_trk_ptk1_reco[ibin]->Sumw2();
        h_trk_ptk1p5_reco[ibin]->Sumw2();
        h_trk_ptk2_reco[ibin]->Sumw2();
        h_trk_ptk2p5_reco[ibin]->Sumw2();

        h_trk_ptk1_reco_refl[ibin]->Sumw2();
        h_trk_ptk1p5_reco_refl[ibin]->Sumw2();
        h_trk_ptk2_reco_refl[ibin]->Sumw2();
        h_trk_ptk2p5_reco_refl[ibin]->Sumw2();

        h_trk_pt_gen_full[ibin]->Sumw2();
        h_trk_pt_gen[ibin]->Sumw2();
        h_trk_pt4_gen[ibin]->Sumw2();
        h_trk_pt5_gen[ibin]->Sumw2();
        h_trk_pt_gen_pos[ibin]->Sumw2();
        h_trk_pt_gen_neg[ibin]->Sumw2();
        h_trk_eta_gen[ibin]->Sumw2();
    }
   
    double weight, pthat_weight, weight_vz, weight_cen;    

	float calo_jteta, calo_jtphi, calo_jtpt, calo_corrpt;
	float calo_refpt, calo_refeta, calo_refphi; 
    float gen_totchg, gen_totchg_pt2, gen_totchg_pt3, gen_totchg_pt4, gen_totchg_pt5, trk_totchg ,trk_totchg_pt2, trk_totchg_pt3, trk_totchg_pt4, trk_totchg_pt5;
    float gen_totchg_sube0_pt2, gen_totchg_sube0_pt3, gen_totchg_sube0_pt4, gen_totchg_sube0_pt5, gen_bkgchg_sube0_pt2, gen_bkgchg_sube0_pt3, gen_bkgchg_sube0_pt4, gen_bkgchg_sube0_pt5;
    float trk_mixchg_pt2, trk_mixchg_pt3, trk_mixchg_pt4, trk_mixchg_pt5;
	float gen_bkgchg_pt2, gen_bkgchg_pt3, gen_bkgchg_pt4, gen_bkgchg_pt5, trk_bkgchg_pt2, trk_bkgchg_pt3, trk_bkgchg_pt4, trk_bkgchg_pt5;
	float trk_totchg_pt2_k1, trk_totchg_pt2_k1p5, trk_totchg_pt2_k2, trk_totchg_pt2_k2p5;
    float trk_bkgchg_pt2_k1, trk_bkgchg_pt2_k1p5, trk_bkgchg_pt2_k2, trk_bkgchg_pt2_k2p5;
	int calo_parton_flavor;
	double trk_correction;
    float trk_totchg_pt2_highest, trk_totchg_pt2_highest_2;

    float trk_totchg_pt2_trkuncup, trk_totchg_pt3_trkuncup, trk_totchg_pt4_trkuncup, trk_totchg_pt5_trkuncup, trk_totchg_pt2_trkuncdn, trk_totchg_pt3_trkuncdn, trk_totchg_pt4_trkuncdn, trk_totchg_pt5_trkuncdn;
    float trk_bkgchg_pt2_trkuncup, trk_bkgchg_pt3_trkuncup, trk_bkgchg_pt4_trkuncup, trk_bkgchg_pt5_trkuncup, trk_bkgchg_pt2_trkuncdn, trk_bkgchg_pt3_trkuncdn, trk_bkgchg_pt4_trkuncdn, trk_bkgchg_pt5_trkuncdn;

    float trk_totchg_pt2_trkuncpos, trk_totchg_pt3_trkuncpos, trk_totchg_pt4_trkuncpos, trk_totchg_pt5_trkuncpos;
    float trk_bkgchg_pt2_trkuncpos, trk_bkgchg_pt3_trkuncpos, trk_bkgchg_pt4_trkuncpos, trk_bkgchg_pt5_trkuncpos;

	Int_t HBHENoiseFilter, HBHENoiseFilterResultRun2Loose, eventSelection, pvFilter, phfCoincFilter3, pclusterCompatibilityFilter, pBeamScrapingFilter;
	Int_t HLT_Jet80, HLT_Jet100, HLT_Jet80_Prescl, HLT_Jet100_Prescl, L1_Jet80_Prescl, L1_Jet100_Prescl;
	Int_t pPAcollisionEventSelectionPA = -999;
	Int_t hiBin = -999;
	Float_t pthat = -999;
	Float_t vz = -999, hiHF= -999.;
	Int_t calo_nref, nTrk, nPFpart, nCSpart, nmixTrk;
	Int_t ncs_id1_pt2, ncs_id1_pt3, ncs_id1_pt4, ncs_id1_pt5, ncs_id145_pt2, ncs_id145_pt3, ncs_id145_pt4, ncs_id145_pt5;  
	Int_t npf_id1_pt2, npf_id1_pt3, npf_id1_pt4, npf_id1_pt5, npf_id145_pt2, npf_id145_pt3, npf_id145_pt4, npf_id145_pt5;
    Int_t n_trk, n_bkgtrk, n_gen, n_bkggen, n_mixtrk, n_sube0_gen, n_bkg_sube0_gen;
    Int_t n_trkpt1, n_trkpt2, n_trkpt3, n_trkpt4, n_trkpt5, n_genpt1, n_genpt2, n_genpt3, n_genpt4, n_genpt5;
    Int_t n_bkgtrkpt1, n_bkgtrkpt2, n_bkgtrkpt3, n_bkgtrkpt4, n_bkgtrkpt5, n_bkggenpt1, n_bkggenpt2, n_bkggenpt3, n_bkggenpt4, n_bkggenpt5;
    double highest_idx_pt, highest_idx_2_pt;
    Int_t highest_idx, highest_idx_2;

	ftree->Branch("vz", &vz);
	if(do_PbPb){
	  ftree->Branch("hiBin", &hiBin);
	}

	ftree->Branch("jteta", &calo_jteta);
	ftree->Branch("jtphi", &calo_jtphi);
	ftree->Branch("jtpt", &calo_jtpt);
	ftree->Branch("corrpt", &calo_corrpt);

	if(!is_data){
	  ftree->Branch("refpt", &calo_refpt);
      ftree->Branch("refeta", &calo_refeta);
      ftree->Branch("refphi", &calo_refphi);
	  ftree->Branch("pthat", &pthat);
      ftree->Branch("refparton_flavor",&calo_parton_flavor);

      ftree->Branch("n_gen",&n_gen);
      ftree->Branch("n_genpt1",&n_genpt1);
      ftree->Branch("n_genpt2",&n_genpt2);
      ftree->Branch("n_genpt3",&n_genpt3);
      ftree->Branch("n_genpt4",&n_genpt4);
      ftree->Branch("n_genpt5",&n_genpt5);
      ftree->Branch("n_bkggen",&n_bkggen);    
      ftree->Branch("n_bkggenpt1",&n_bkggenpt1);    
      ftree->Branch("n_bkggenpt2",&n_bkggenpt2);    
      ftree->Branch("n_bkggenpt3",&n_bkggenpt3);    
      ftree->Branch("n_bkggenpt4",&n_bkggenpt4);    
      ftree->Branch("n_bkggenpt5",&n_bkggenpt5);    

      ftree->Branch("gen_totchg_pt2",&gen_totchg_pt2);
      ftree->Branch("gen_totchg_pt3",&gen_totchg_pt3);
      ftree->Branch("gen_totchg_pt4",&gen_totchg_pt4);
      ftree->Branch("gen_totchg_pt5",&gen_totchg_pt5);

      ftree->Branch("gen_bkgchg_pt2",&gen_bkgchg_pt2);
      ftree->Branch("gen_bkgchg_pt3",&gen_bkgchg_pt3);
      ftree->Branch("gen_bkgchg_pt4",&gen_bkgchg_pt4);
      ftree->Branch("gen_bkgchg_pt5",&gen_bkgchg_pt5); 

      if(do_sube0){ 
        ftree->Branch("n_sube0_gen",&n_sube0_gen);
        ftree->Branch("gen_totchg_sube0_pt2",&gen_totchg_sube0_pt2);
        ftree->Branch("gen_totchg_sube0_pt3",&gen_totchg_sube0_pt3);
        ftree->Branch("gen_totchg_sube0_pt4",&gen_totchg_sube0_pt4);
        ftree->Branch("gen_totchg_sube0_pt5",&gen_totchg_sube0_pt5);
        ftree->Branch("n_bkg_sube0_gen",&n_bkg_sube0_gen);
        ftree->Branch("gen_bkgchg_sube0_pt2",&gen_bkgchg_sube0_pt2);
        ftree->Branch("gen_bkgchg_sube0_pt3",&gen_bkgchg_sube0_pt3);
        ftree->Branch("gen_bkgchg_sube0_pt4",&gen_bkgchg_sube0_pt4);
        ftree->Branch("gen_bkgchg_sube0_pt5",&gen_bkgchg_sube0_pt5);
      }
	}

	ftree->Branch("weight", &weight);
  
	ftree->Branch("n_trk",&n_trk);    
    ftree->Branch("n_trkpt1",&n_trkpt1);
    ftree->Branch("n_trkpt2",&n_trkpt2);
    ftree->Branch("n_trkpt3",&n_trkpt3);
    ftree->Branch("n_trkpt4",&n_trkpt4);
    ftree->Branch("n_trkpt5",&n_trkpt5);

    ftree->Branch("n_bkgtrk",&n_bkgtrk);    
    ftree->Branch("n_bkgtrkpt1",&n_bkgtrkpt1);
    ftree->Branch("n_bkgtrkpt2",&n_bkgtrkpt2);
    ftree->Branch("n_bkgtrkpt3",&n_bkgtrkpt3);
    ftree->Branch("n_bkgtrkpt4",&n_bkgtrkpt4);
    ftree->Branch("n_bkgtrkpt5",&n_bkgtrkpt5);

	ftree->Branch("n_mixtrk",&n_mixtrk);    

    ftree->Branch("trk_totchg",&trk_totchg);

    ftree->Branch("trk_totchg_pt2_highest",&trk_totchg_pt2_highest);
    ftree->Branch("trk_totchg_pt2_highest_2",&trk_totchg_pt2_highest_2);

    ftree->Branch("trk_totchg_pt2",&trk_totchg_pt2);
    ftree->Branch("trk_totchg_pt3",&trk_totchg_pt3);
    ftree->Branch("trk_totchg_pt4",&trk_totchg_pt4);
    ftree->Branch("trk_totchg_pt5",&trk_totchg_pt5);

    ftree->Branch("trk_bkgchg_pt2",&trk_bkgchg_pt2);
    ftree->Branch("trk_bkgchg_pt3",&trk_bkgchg_pt3);
    ftree->Branch("trk_bkgchg_pt4",&trk_bkgchg_pt4);
    ftree->Branch("trk_bkgchg_pt5",&trk_bkgchg_pt5);

    ftree->Branch("trk_totchg_pt2_k1",&trk_totchg_pt2_k1);
    ftree->Branch("trk_totchg_pt2_k1p5",&trk_totchg_pt2_k1p5);
    ftree->Branch("trk_totchg_pt2_k2",&trk_totchg_pt2_k2);
    ftree->Branch("trk_totchg_pt2_k2p5",&trk_totchg_pt2_k2p5);

    ftree->Branch("trk_bkgchg_pt2_k1",&trk_bkgchg_pt2_k1);
    ftree->Branch("trk_bkgchg_pt2_k1p5",&trk_bkgchg_pt2_k1p5);
    ftree->Branch("trk_bkgchg_pt2_k2",&trk_bkgchg_pt2_k2);
    ftree->Branch("trk_bkgchg_pt2_k2p5",&trk_bkgchg_pt2_k2p5);
/*
    ftree->Branch("trk_totchg_pt2_trkuncup",&trk_totchg_pt2_trkuncup);
    ftree->Branch("trk_totchg_pt3_trkuncup",&trk_totchg_pt3_trkuncup);
    ftree->Branch("trk_totchg_pt4_trkuncup",&trk_totchg_pt4_trkuncup);
    ftree->Branch("trk_totchg_pt5_trkuncup",&trk_totchg_pt5_trkuncup);

    ftree->Branch("trk_bkgchg_pt2_trkuncup",&trk_bkgchg_pt2_trkuncup);
    ftree->Branch("trk_bkgchg_pt3_trkuncup",&trk_bkgchg_pt3_trkuncup);
    ftree->Branch("trk_bkgchg_pt4_trkuncup",&trk_bkgchg_pt4_trkuncup);
    ftree->Branch("trk_bkgchg_pt5_trkuncup",&trk_bkgchg_pt5_trkuncup);
 
    ftree->Branch("trk_totchg_pt2_trkuncdn",&trk_totchg_pt2_trkuncdn);
    ftree->Branch("trk_totchg_pt3_trkuncdn",&trk_totchg_pt3_trkuncdn);
    ftree->Branch("trk_totchg_pt4_trkuncdn",&trk_totchg_pt4_trkuncdn);
    ftree->Branch("trk_totchg_pt5_trkuncdn",&trk_totchg_pt5_trkuncdn);

    ftree->Branch("trk_bkgchg_pt2_trkuncdn",&trk_bkgchg_pt2_trkuncdn);
    ftree->Branch("trk_bkgchg_pt3_trkuncdn",&trk_bkgchg_pt3_trkuncdn);
    ftree->Branch("trk_bkgchg_pt4_trkuncdn",&trk_bkgchg_pt4_trkuncdn);
    ftree->Branch("trk_bkgchg_pt5_trkuncdn",&trk_bkgchg_pt5_trkuncdn);

    ftree->Branch("trk_totchg_pt2_trkuncpos",&trk_totchg_pt2_trkuncpos);
    ftree->Branch("trk_totchg_pt3_trkuncpos",&trk_totchg_pt3_trkuncpos);
    ftree->Branch("trk_totchg_pt4_trkuncpos",&trk_totchg_pt4_trkuncpos);
    ftree->Branch("trk_totchg_pt5_trkuncpos",&trk_totchg_pt5_trkuncpos);
*/

    if(do_mixing){
      ftree->Branch("trk_mixchg_pt2",&trk_mixchg_pt2);
      ftree->Branch("trk_mixchg_pt3",&trk_mixchg_pt3);
      ftree->Branch("trk_mixchg_pt4",&trk_mixchg_pt4);
      ftree->Branch("trk_mixchg_pt5",&trk_mixchg_pt5);
    }
 
    if(do_PbPb){
    	/*
        ftree->Branch("ncs_id1_pt2",&ncs_id1_pt2);
        ftree->Branch("ncs_id1_pt3",&ncs_id1_pt3);
        ftree->Branch("ncs_id1_pt4",&ncs_id1_pt4);
        ftree->Branch("ncs_id1_pt5",&ncs_id1_pt5);
        */ 
        ftree->Branch("ncs_id145_pt2",&ncs_id145_pt2);
        //ftree->Branch("ncs_id145_pt3",&ncs_id145_pt3);
        //ftree->Branch("ncs_id145_pt4",&ncs_id145_pt4);
        //ftree->Branch("ncs_id145_pt5",&ncs_id145_pt5);
	}
	else{
	    ftree->Branch("npf_id1_pt2",&npf_id1_pt2);
        //ftree->Branch("npf_id1_pt3",&npf_id1_pt3);
        //ftree->Branch("npf_id1_pt4",&npf_id1_pt4);
        //ftree->Branch("npf_id1_pt5",&npf_id1_pt5);
        /*
	    ftree->Branch("npf_id145_pt2",&npf_id145_pt2);
        ftree->Branch("npf_id145_pt3",&npf_id145_pt3);
        ftree->Branch("npf_id145_pt4",&npf_id145_pt4);
        ftree->Branch("npf_id145_pt5",&npf_id145_pt5);
	    */
	}
    
    std::ifstream instr(in_file_name.c_str(), std::ifstream::in);
	if(!instr.is_open()) cout << "filelist not found!! Exiting..." << endl;
	std::string filename;
	int ifile=0;

	const int MAXJETS = 500;
	const int MAXPARTICLES = 100000;
	
	Float_t t_calo_jtpt[MAXJETS], t_calo_jteta[MAXJETS], t_calo_jtphi[MAXJETS], t_calo_trackmax[MAXJETS], t_calo_rawpt[MAXJETS];
    Float_t t_calo_refpt[MAXJETS]; Float_t t_calo_refeta[MAXJETS]; Float_t t_calo_refphi[MAXJETS];
    Int_t t_calo_parton_flavor[MAXJETS];	

	Float_t trkPt[MAXPARTICLES], trkEta[MAXPARTICLES], trkPhi[MAXPARTICLES], trkDxy1[MAXPARTICLES], trkDxyError1[MAXPARTICLES], trkDz1[MAXPARTICLES], trkDzError1[MAXPARTICLES], trkPtError[MAXPARTICLES], pfEcal[MAXPARTICLES], pfHcal[MAXPARTICLES], trkChi2[MAXPARTICLES];
	//Bool_t trkMVALoose[MAXPARTICLES], trkMVATight[MAXPARTICLES];
	UChar_t trkAlgo[MAXPARTICLES], trkNHit[MAXPARTICLES], trkNlayer[MAXPARTICLES], trkNdof[MAXPARTICLES]; //Run2
	Bool_t highPurity[MAXPARTICLES];
    Int_t trk_chg[MAXPARTICLES];

	Float_t me_trkPt[MAXPARTICLES], me_trkEta[MAXPARTICLES], me_trkPhi[MAXPARTICLES], me_trkDxy[MAXPARTICLES], me_trkDxyError[MAXPARTICLES], me_trkDz[MAXPARTICLES], me_trkDzError[MAXPARTICLES], me_trkPtError[MAXPARTICLES], me_pfEcal[MAXPARTICLES], me_pfHcal[MAXPARTICLES], me_trkChi2[MAXPARTICLES];
	//Bool_t trkMVALoose[MAXPARTICLES], trkMVATight[MAXPARTICLES];
	UChar_t me_trkAlgo[MAXPARTICLES], me_trkNHit[MAXPARTICLES], me_trkNlayer[MAXPARTICLES], me_trkNdof[MAXPARTICLES]; //Run2
	Bool_t me_highPurity[MAXPARTICLES];
    Int_t me_trkchg[MAXPARTICLES];

	vector<int> *gen_chg=0, *sube=0;
	vector<float> *genEta=0, *genPhi=0, *genPt=0; 

	vector<float> *pfCandEta=0, *pfCandPhi=0, *pfCandPt=0;
    vector<int> *pfCandId=0;
	vector<float> *csCandEta=0, *csCandPhi=0, *csCandPt=0;
	vector<int> *csCandId=0;

    int entryCalls=0;
	unsigned int meptrig = 1;

	gRandom->SetSeed(0);
	const int nCentMixBins= is_pp ? 1 : 40;
	const int nVzMixBins=30;

	TH1D * centbins = new TH1D("centbins","centbins. JUST A DUMMY REALLY", nCentMixBins, 0.0, 200.0);
	TH1D * vzbins = new TH1D("vzbins","vzbins. JUST A DUMMY REALLY", nVzMixBins, -15., 15.);

	TTree *mixing_event_tree;
	TTree *mixing_track_tree;
	TTree *mixing_hlt_tree;

    if(do_mixing){
      int imix=0;
	  //if(is_pp) me_tree->Add("root://cmsxrootd.fnal.gov///store/user/kjung/ppData_5TeV_inclJet_FinalJFFcorr_skims/Merged/Data_pp_Merged10Files.root");
	  if(is_pp) {
	  	TFile *mixing_file = TFile::Open("root://cmsxrootd.fnal.gov///store/user/tatar/MinimumBias6/Run2015E_PromptReco_v1_Run261553_262328_FOREST/2.root");

        mixing_event_tree = (TTree*) mixing_file->Get("hiEvtAnalyzer/HiTree");
        mixing_track_tree = (TTree*) mixing_file->Get("ppTrack/trackTree");
        mixing_hlt_tree = (TTree*) mixing_file->Get("skimanalysis/HltTree");
      } 
      else{
        TFile *mixing_file = TFile::Open("root://cmsxrootd.fnal.gov///store/user/rbi/merged/HIMinimumBias2-HIRun2015-PromptReco-v1_forest_csjet_v1/0.root");
        
        mixing_event_tree = (TTree*) mixing_file->Get("hiEvtAnalyzer/HiTree");
        mixing_track_tree = (TTree*) mixing_file->Get("anaTrack/trackTree");
        mixing_hlt_tree = (TTree*) mixing_file->Get("skimanalysis/HltTree");
      }
      cout << "total entries: "<<mixing_event_tree->GetEntries()<< endl;
    }
      
	float me_vz, me_pthat;
	int me_hiBin;
	//vector<float> *me_trkPt=0, *me_trkEta=0, *me_trkPhi=0, *me_jtpt=0, *me_jteta=0, *me_jtphi=0, *me_trkChi2=0, *me_pfHcal=0, *me_pfEcal=0;
	//vector<bool> *me_highPurity=0;
	//vector<int> *me_trkNlayer=0, *me_trkNHit=0, *me_trkNdof=0;
    //vector<float> *me_trkDxy=0, *me_trkDxyError=0, *me_trkDz=0, *me_trkDzError=0;
	Int_t me_HBHEFilter, me_collisionEventSelection, me_phfCoincFilter, me_pprimaryVertexFilter;
    std::vector<std::vector<std::vector<unsigned int> > > mixing_lists;
    unsigned int jet_cent, jet_vzbin;

	TStopwatch *mixingWatch = new TStopwatch();

        if(do_mixing){

            mixing_track_tree->SetBranchAddress("nTrk", &nmixTrk);
            mixing_track_tree->SetBranchAddress("trkPt", &me_trkPt);
            mixing_track_tree->SetBranchAddress("trkEta", &me_trkEta);
            mixing_track_tree->SetBranchAddress("trkPhi", &me_trkPhi);
            mixing_track_tree->SetBranchAddress("trkCharge", &me_trkchg);
            mixing_track_tree->SetBranchAddress("highPurity", &me_highPurity);
            mixing_track_tree->SetBranchAddress("trkChi2", &me_trkChi2);
            mixing_track_tree->SetBranchAddress("trkNdof", &me_trkNdof);
            mixing_track_tree->SetBranchAddress("trkNlayer", &me_trkNlayer);
            mixing_track_tree->SetBranchAddress("trkNHit", &me_trkNHit);
            mixing_track_tree->SetBranchAddress("pfHcal", &me_pfHcal);
            mixing_track_tree->SetBranchAddress("pfEcal", &me_pfEcal);
            mixing_track_tree->SetBranchAddress("trkDxy1",&me_trkDxy);
            mixing_track_tree->SetBranchAddress("trkDxyError1",&me_trkDxyError);
            mixing_track_tree->SetBranchAddress("trkDz1",&me_trkDz);
            mixing_track_tree->SetBranchAddress("trkDzError1",&me_trkDzError);

            mixing_event_tree->SetBranchAddress("vz", &me_vz);  
            if(!is_pp) mixing_event_tree->SetBranchAddress("hiBin", &me_hiBin);
            //mixing_event_tree->SetBranchAddress("calo_jteta", &me_jteta);
            //mixing_event_tree->SetBranchAddress("calo_jtphi", &me_jtphi);
            //mixing_event_tree->SetBranchAddress("calo_jtpt", &me_jtpt);
            if(!is_pp){
                mixing_hlt_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&me_HBHEFilter);
                mixing_hlt_tree->SetBranchAddress("pcollisionEventSelection",&me_collisionEventSelection);
                mixing_hlt_tree->SetBranchAddress("pprimaryVertexFilter",&me_pprimaryVertexFilter);
                mixing_hlt_tree->SetBranchAddress("phfCoincFilter3",&me_phfCoincFilter);
            }
            else{
                //mixing_hlt_tree->SetBranchAddress("pcollisionEventSelection",&me_collisionEventSelection);
                mixing_hlt_tree->SetBranchAddress("pPAprimaryVertexFilter",&me_pprimaryVertexFilter);
                mixing_hlt_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&me_HBHEFilter);
            }

            for(int ivz=0; ivz<nVzMixBins; ivz++){
                vector<vector<unsigned int> > dummyVect;
                for(int icent=0; icent<nCentMixBins; icent++){
                    vector<unsigned int> dummyVect2;
                    dummyVect.push_back(dummyVect2);
                }
                mixing_lists.push_back(dummyVect);
            }

            cout << "sorting events into vz and centrality categories..." << endl;
            int totalEntries = mixing_event_tree->GetEntries();
            for(int me_evt=0; me_evt<totalEntries; me_evt++){
            //for(int me_evt=0; me_evt<10000; me_evt++){
                if(me_evt && me_evt%5000==0) cout << "me_evt: "<< me_evt << " of " << mixing_event_tree->GetEntries() << endl;
                mixing_hlt_tree->GetEntry(me_evt);
                mixing_event_tree->GetEntry(me_evt);
                if(!is_pp && (!me_HBHEFilter || !me_collisionEventSelection || !me_phfCoincFilter || !me_pprimaryVertexFilter)) continue;
                //{ cout << "event " << me_evt << " failed" << endl; continue; }
                if(is_pp && !me_pprimaryVertexFilter && !me_HBHEFilter) continue;
                if(abs(me_vz)>=15) { continue; }
                if(!is_data) if(me_pthat<80) { continue; }

                unsigned int me_cent = 0;
                if(!is_pp) me_cent = centbins->FindBin(me_hiBin)-1;
                unsigned int me_vzbin = vzbins->FindBin(me_vz)-1;
                //cout << "pushing back centrality: "<< me_hiBin << " vz: "<< me_vz << " cent bin: "<< me_cent << " vz bin: " << me_vzbin << endl;
                //cout << "cent bin range: "<< centbins->GetBinLowEdge(me_cent+1) << " " << centbins->GetBinLowEdge(me_cent+2) << endl;
                if(mixing_lists.at(me_vzbin).at(me_cent).size()<100){
                    mixing_lists.at(me_vzbin).at(me_cent).push_back(me_evt);
                }
            }
            cout << "sizes..." << endl;
            int totBins = nVzMixBins*nCentMixBins;
            int filledBins=0;
            for(int imix=0; imix<nVzMixBins; imix++){
                for(int icent=0; icent<nCentMixBins; icent++){
                    if(mixing_lists.at(imix).at(icent).size() >= 2) filledBins++;
                }
            }
            
            cout << "lists loaded!" << endl;
            cout << "MIXING DIAGNOSTICS: " << filledBins << " of " << totBins << " bins filled, " << 100.*(double)filledBins/(double)totBins << "%" << endl;
        }

        mixingWatch->Stop();
	    cout << "Mixing load time: "<< mixingWatch->RealTime() << endl;

	while(instr>>filename && ifile<endfile){
		filename.erase(std::remove(filename.begin(), filename.end(), '"'), filename.end());
		filename.erase(std::remove(filename.begin(), filename.end(), ','), filename.end());
		filename.erase(std::remove(filename.begin(), filename.end(), '['), filename.end());
		filename.erase(std::remove(filename.begin(), filename.end(), ']'), filename.end());
		cout<<"File name is "<< filename <<endl;
		ifile++;

		//TFile *my_file = TFile::Open(filename.c_str());

		//if(!my_file){
			int pos = filename.find_first_of('s');
			string reducedfn = filename.substr(pos-1);
			string xrdPrefix = "root://cmsxrootd.fnal.gov//";
			cout << "local file not detected. Trying " << xrdPrefix+reducedfn << endl;
			TFile *my_file = TFile::Open((xrdPrefix+reducedfn).c_str());
			//TFile::Open((xrdPrefix+reducedfn).c_str());
		//}

        //TFile *my_file = TFile::Open("root://cmsxrootd.fnal.gov///store/user/kjung/pp_5TeV_HighPtJet80PD_RecalibJP_Apr2016/HiForestAOD_88.root");
        //TFile *my_file = TFile::Open("root://cmsxrootd.fnal.gov///store/user/kjung/PbPbMC_Pythia6HydjetCymbalTune/Pythia6_Dijet220_pp502_Hydjet_Cymbal_MB/Pythia6_Dijet220_pp502_Hydjet_Cymbal_MB/crab_Pythia6_Dijet220_pp502_Hydjet_Cymbal_MB_reprocess/170629_224823/0000/HiForestAOD_88.root");
        //TFile *my_file = TFile::Open("HiForestAOD_P+H.root");

                if(!my_file){ 
                    cout << "File cannot be found!!" << endl; 
		    continue; 
		}	

		if(my_file->IsZombie()) { 
		    cout << "Is zombie" << endl;
		    continue;
		}

		if(do_PbPb){
			inp_tree = (TTree*)  my_file->Get(Form("akPu%dCaloJetAnalyzer/t",radius));
     		inp_tree_CS = (TTree*) my_file->Get("pfcandAnalyzerCS/pfTree");
	    	if(!inp_tree_CS){ cout << "CSCand Tree not found!! Exiting..." << endl; exit(1); }
		}else{
			inp_tree = (TTree*)  my_file->Get(Form("ak%dCaloJetAnalyzer/t",radius));
		}

		inp_tree2 = (TTree*)  my_file->Get("pfcandAnalyzer/pfTree");
		if(!inp_tree2){ cout << "PFCand Tree not found!! Exiting..." << endl; exit(1); }

		inp_tree3 = (TTree*) my_file->Get("hiEvtAnalyzer/HiTree");
		if(!inp_tree3){ cout << "Evt Tree not found!! Exiting..." << endl; exit(1); }

		inp_tree4 = (TTree*) my_file->Get("skimanalysis/HltTree");
		if(!inp_tree4){ cout << "Skim Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree4);
    
		inp_tree6 = (TTree*) my_file->Get("hltanalysis/HltTree");
		if(is_data && !inp_tree6){ cout << "HLT Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree6);

		inp_tree7 = (TTree*) my_file->Get("HiGenParticleAna/hi");
		if(!is_data && !inp_tree7){ cout << "HiGenParticleAna Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree7);

		if(do_PbPb){
    		inp_tree8 = (TTree*) my_file->Get("anaTrack/trackTree");
		}else{
			inp_tree8 = (TTree*) my_file->Get("ppTrack/trackTree");
		}
		if(!inp_tree8){ cout << "track Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree8);

		cout << "trees loaded" << endl;

		inp_tree3->SetBranchAddress("vz",&vz);
        inp_tree3->SetBranchAddress("hiHF",&hiHF);		

		if(do_PbPb){
			inp_tree3->SetBranchAddress("hiBin",&hiBin);
		}
		inp_tree->SetBranchAddress("nref",&calo_nref);
		inp_tree->SetBranchAddress("jtpt",t_calo_jtpt);
		inp_tree->SetBranchAddress("jteta",t_calo_jteta);
		inp_tree->SetBranchAddress("jtphi",t_calo_jtphi);
		inp_tree->SetBranchAddress("rawpt",t_calo_rawpt);
		inp_tree->SetBranchAddress("trackMax",t_calo_trackmax);

		inp_tree2->SetBranchAddress("nPFpart",&nPFpart);
		inp_tree2->SetBranchAddress("pfPt",&pfCandPt);
        inp_tree2->SetBranchAddress("pfEta",&pfCandEta);
        inp_tree2->SetBranchAddress("pfPhi",&pfCandPhi);
        inp_tree2->SetBranchAddress("pfId",&pfCandId);		
		
		if(!is_data){
    		inp_tree->SetBranchAddress("refpt",t_calo_refpt);
            inp_tree->SetBranchAddress("refeta",t_calo_refeta);
            inp_tree->SetBranchAddress("refphi",t_calo_refphi);            
	     	inp_tree->SetBranchAddress("refparton_flavor",t_calo_parton_flavor);
    	  	inp_tree3->SetBranchAddress("pthat",&pthat);
 
            inp_tree->SetBranchAddress("pt", &genPt);
            inp_tree->SetBranchAddress("eta", &genEta);
            inp_tree->SetBranchAddress("phi", &genPhi);
            inp_tree->SetBranchAddress("chg", &gen_chg);
            inp_tree->SetBranchAddress("sube", &sube);
        }

        inp_tree->SetBranchAddress("trkPt", &trkPt);
        inp_tree->SetBranchAddress("trkEta", &trkEta);
        inp_tree->SetBranchAddress("trkPhi", &trkPhi);
        inp_tree->SetBranchAddress("trkCharge", &trk_chg);
        inp_tree->SetBranchAddress("trkAlgo", &trkAlgo);
        inp_tree->SetBranchAddress("highPurity", &highPurity);

		inp_tree->SetBranchAddress("nTrk",&nTrk);         
        inp_tree->SetBranchAddress("trkDxy1", &trkDxy1);
        inp_tree->SetBranchAddress("trkDxyError1", &trkDxyError1);
        inp_tree->SetBranchAddress("trkDz1", &trkDz1);
        inp_tree->SetBranchAddress("trkDzError1", &trkDzError1);
        inp_tree->SetBranchAddress("trkPtError", &trkPtError);
        inp_tree->SetBranchAddress("trkNHit",&trkNHit);
        inp_tree->SetBranchAddress("trkNlayer",&trkNlayer);
        inp_tree->SetBranchAddress("trkNdof",&trkNdof);
        inp_tree->SetBranchAddress("trkChi2",&trkChi2);
        inp_tree->SetBranchAddress("pfEcal",&pfEcal);
        inp_tree->SetBranchAddress("pfHcal",&pfHcal);
        //inp_tree->SetBranchAddress("trkMVALoose",&trkMVALoose);
        //inp_tree->SetBranchAddress("trkMVATight",&trkMVATight);        
        
        if(do_PbPb){
			inp_tree_CS->SetBranchAddress("nPFpart",&nCSpart);
			inp_tree_CS->SetBranchAddress("pfPt",&csCandPt);
			inp_tree_CS->SetBranchAddress("pfEta",&csCandEta);
			inp_tree_CS->SetBranchAddress("pfPhi",&csCandPhi);
			inp_tree_CS->SetBranchAddress("pfId",&csCandId);
        }

		inp_tree->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilter);
		inp_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHENoiseFilterResultRun2Loose);

		if(!do_PbPb){
			inp_tree->SetBranchAddress("pPAprimaryVertexFilter",&pvFilter);
			inp_tree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
		}
		else{	
			inp_tree->SetBranchAddress("pcollisionEventSelection",&eventSelection);			
			inp_tree->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
                        inp_tree->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3);
 	                inp_tree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);     	
                }

		if(is_data){
			if(!do_PbPb){ 
				inp_tree->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1", &HLT_Jet80);
				inp_tree->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1", &HLT_Jet100);
			}
			else{
				inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1", &HLT_Jet80);
				inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1", &HLT_Jet100);
				inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1_Prescl", &HLT_Jet80_Prescl);
				inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1_Prescl", &HLT_Jet100_Prescl);
				inp_tree->SetBranchAddress("L1_SingleJet44_BptxAND_Prescl", &L1_Jet80_Prescl);
				inp_tree->SetBranchAddress("L1_SingleS1Jet56_BptxAND_Prescl", &L1_Jet100_Prescl);				
			}
		}


		int n_evt = inp_tree->GetEntriesFast();

		cout << "Entries: "<< n_evt << endl;

		TStopwatch evtStopwatch;
        TStopwatch mixStopwatch;
        TStopwatch mixStopwatch_p1;
        TStopwatch mixStopwatch_p2;
        TStopwatch mixStopwatch_p3;

		evtStopwatch.Start(1);

		for(int evi = 0; evi < n_evt; evi++) {
		//for(int evi = 0; evi < 100; evi++) {
			
            inp_tree->GetEntry(evi);
			inp_tree2->GetEntry(evi);
			inp_tree3->GetEntry(evi);
			if(do_PbPb) inp_tree_CS->GetEntry(evi);

            if(fabs(vz) >= 15. || hiHF > 5500.) continue;
			
			if(is_data && !HLT_Jet80 && !HLT_Jet100) continue;
            
            if(!do_PbPb && (!pvFilter || !pBeamScrapingFilter || !HBHENoiseFilterResultRun2Loose)) continue;

            if(do_PbPb && (!pvFilter || !eventSelection || !phfCoincFilter3 || !pclusterCompatibilityFilter || !HBHENoiseFilterResultRun2Loose)) continue;

            weight = 1.;
            
            if(is_data){
    			if(HLT_Jet80 && !HLT_Jet100){
                    weight = L1_Jet80_Prescl*HLT_Jet80_Prescl;
		    	}
                weight = 1.;
            }

     		else {
                int ibin=0;
                while(pthat>pthatbins_PH[ibin+1]) ibin++;
                pthat_weight = (xsecs_PH[ibin]-xsecs_PH[ibin+1])/pthatEntries_PH[ibin];

                //pythia6 
                if(is_pp) weight_vz = 1./fvz_pp->Eval(vz);
                //Pythia6+Hydjet
                else weight_vz = fvz_cymbal->Eval(vz);

                if(do_PbPb) weight_cen = fcent_cymbal(hiBin,f_cent);
                else weight_cen = 1.;

                weight = pthat_weight*weight_vz*weight_cen;
            }
            
			jet_cent = centbins->FindBin(hiBin)-1;
			if(is_pp) jet_cent=0;
            jet_vzbin = vzbins->FindBin(vz)-1;

            if(do_mixing){
				if(mixing_lists.at(jet_vzbin).at(jet_cent).size()<2) continue;
			}

            if( !is_data && pthat<80. ) continue;

            for (int cbin = 0; cbin < nCbins; cbin++){ 
              if (hiBin > CBins[cbin] && hiBin <= CBins[cbin+1]){
                mycbin = cbin; 
              }
            }
            if(is_pp) mycbin = 0;

            h_hibin->Fill(hiBin,weight);
            h_vz->Fill(vz,weight);
            h_pthat->Fill(pthat,weight);

            if(!is_data){
                for(int gen =0; gen < (int) genPt->size(); gen++){
                    if(!passGenTrackCuts(genPt->at(gen), genEta->at(gen), gen_chg->at(gen))) continue;

                    h_trk_pt_gen_full[mycbin]->Fill(genPt->at(gen),weight);
                    if(gen_chg->at(gen) == 1) h_trk_pt_gen_pos[mycbin]->Fill(genPt->at(gen),weight);
                    else if(gen_chg->at(gen) == -1) h_trk_pt_gen_neg[mycbin]->Fill(genPt->at(gen),weight);
                    h_trk_eta_gen[mycbin]->Fill(genEta->at(gen),weight);
                }
            }

            for(int itrk =0; itrk < nTrk; itrk++){
                if(!passTrackCuts(is_pp, 1, trkPt[itrk], trkEta[itrk], highPurity[itrk], trkChi2[itrk], trkNdof[itrk], trkNlayer[itrk], trkNHit[itrk], pfHcal[itrk], pfEcal[itrk], trkDxy1[itrk]/trkDxyError1[itrk], trkDz1[itrk]/trkDzError1[itrk])) continue;
                //float rmin = 999;
                //trk_correction = trkCorr->getTrkCorr(trkPt[itrk],trkEta[itrk],trkPhi[itrk],hiBin,rmin); 
                trk_correction = trkCorr->getTrkCorr(trkPt[itrk],trkEta[itrk],trkPhi[itrk],hiBin); 
                
                h_trk_pt_reco_full[mycbin]->Fill(trkPt[itrk],trk_correction*weight);
                if(trk_chg[itrk] == 1) h_trk_pt_reco_pos[mycbin]->Fill(trkPt[itrk],trk_correction*weight);
                else if(trk_chg[itrk] == -1) h_trk_pt_reco_neg[mycbin]->Fill(trkPt[itrk],trk_correction*weight);
                h_trk_eta_reco[mycbin]->Fill(trkEta[itrk],trk_correction*weight);
            }

			//start calo jet loop
			for(int j4i = 0; j4i < calo_nref ; j4i++) {

				if(fabs(t_calo_jteta[j4i]) > 1.6 || t_calo_jtpt[j4i] < 50.) continue;
				else if(t_calo_trackmax[j4i]/t_calo_rawpt[j4i] > 0.98 ||t_calo_trackmax[j4i]/t_calo_rawpt[j4i] < 0.01) continue;

                calo_jtpt = t_calo_jtpt[j4i];
                calo_jteta = t_calo_jteta[j4i];
                calo_jtphi = t_calo_jtphi[j4i];
                if(!is_data){
                  calo_refpt = t_calo_refpt[j4i];
                  calo_parton_flavor = t_calo_parton_flavor[j4i];
				}
                if (calo_refpt < 0.) continue;
                calo_refeta = t_calo_refeta[j4i];
                calo_refphi = t_calo_refphi[j4i];

                //npf cand for pp and Pythia
				if(!do_PbPb){
    				npf_id1_pt2 = 0; //npf_id1_pt3 = 0; npf_id1_pt4 = 0; npf_id1_pt5 = 0;
    				//npf_id145_pt2 = 0; npf_id145_pt3 = 0; npf_id145_pt4 = 0; npf_id145_pt5 = 0;

	    			for(int icand=0; icand<nPFpart; icand++){
		    			if(abs(pfCandEta->at(icand))>2.4) continue;
			    		double dr = sqrt(pow(pfCandEta->at(icand) - calo_jteta, 2) + pow(acos(cos(calo_jtphi - pfCandPhi->at(icand))),2));
				    	if(dr > 0.4 || pfCandPt->at(icand)<2.) continue;
				    	if(pfCandId->at(icand)==1 || pfCandId->at(icand)==4 || pfCandId->at(icand)==5){
					    	
					    	//npf_id145_pt2++;
					    	if(pfCandId->at(icand)==1) npf_id1_pt2++;
					    	/*
					    	if(pfCandPt->at(icand)>3.){
					    		npf_id145_pt3++;
					    	    if(pfCandId->at(icand)==1) npf_id1_pt3++;
					    	}
					    	if(pfCandPt->at(icand)>4.){
					            npf_id145_pt4++;
					    	    if(pfCandId->at(icand)==1) npf_id1_pt4++;
					    	} 
					    	if(pfCandPt->at(icand)>5.){
					    	    npf_id145_pt5++;
					    	    if(pfCandId->at(icand)==1) npf_id1_pt5++;
					    	}
					    	*/
				        }
			    	}
			    	calo_corrpt = corrpt->getCorrection(is_pp, npf_id1_pt2, 0, calo_jtpt, calo_jteta);
			    }
			    else{
    				//ncs_id1_pt2 = 0; ncs_id1_pt3 = 0; ncs_id1_pt4 = 0; ncs_id1_pt5 = 0;
    				ncs_id145_pt2 = 0; //ncs_id145_pt3 = 0; ncs_id145_pt4 = 0; ncs_id145_pt5 = 0;

	    			for(int icand=0; icand<nCSpart; icand++){
		    			if(abs(csCandEta->at(icand))>2.4) continue;
			    		double dr = sqrt(pow(csCandEta->at(icand) - calo_jteta, 2) + pow(acos(cos(calo_jtphi - csCandPhi->at(icand))),2));
				    	if(dr > 0.4 || csCandPt->at(icand)<2.) continue;
				    	if(csCandId->at(icand)==1 || csCandId->at(icand)==4 || csCandId->at(icand)==5){
					    	
					    	ncs_id145_pt2++;
					    	/*
					    	if(csCandId->at(icand)==1) ncs_id1_pt2++;
					    	
					    	if(csCandPt->at(icand)>3.){
					    		ncs_id145_pt3++;
					    	    if(csCandId->at(icand)==1) ncs_id1_pt3++;
					    	}
					    	if(csCandPt->at(icand)>4.){
					            ncs_id145_pt4++;
					    	    if(csCandId->at(icand)==1) ncs_id1_pt4++;
					    	} 
					    	if(csCandPt->at(icand)>5.){
					    	    ncs_id145_pt5++;
					    	    if(csCandId->at(icand)==1) ncs_id1_pt5++;
					    	}
					    	*/
				        }
			    	}
			    	//calo_corrpt = corrpt->getCorrection(is_pp, ncs_id145_pt2, hiBin, calo_jtpt, calo_jteta);
                    if(nCS_2 > 20 && calo_jtpt <= 120.) calo_corrpt = corrpt->getCorrection(ispp, 20, hiBin, calo_jtpt, calo_jteta);
                    else calo_corrpt = corrpt->getCorrection(ispp, nCS_2, hiBin, calo_jtpt, calo_jteta);
			    }	
                
                if(calo_corrpt <= 0.) continue;
                //if(calo_corrpt<200.) calo_corrpt = calo_corrpt*(f_res_gauss_lowpt->GetRandom());
                //else calo_corrpt = calo_corrpt*(f_res_gauss_highpt->GetRandom());
                //gen particle and gen chg
                if(!is_data){
                  n_gen = 0; n_bkggen = 0; n_sube0_gen = 0; n_bkg_sube0_gen = 0;
                  n_genpt1 = 0; n_genpt2 = 0; n_genpt3 = 0; n_genpt4 = 0; n_genpt5 = 0;
                  n_bkggenpt1 = 0; n_bkggenpt2 = 0; n_bkggenpt3 = 0; n_bkggenpt4 = 0; n_bkggenpt5 = 0;
                  gen_totchg_pt2 = 0.; gen_totchg_pt3 = 0.; gen_totchg_pt4 = 0.; gen_totchg_pt5 = 0.;
                  gen_bkgchg_pt2 = 0.; gen_bkgchg_pt3 = 0.; gen_bkgchg_pt4 = 0.; gen_bkgchg_pt5 = 0.;

                  gen_totchg_sube0_pt2 = 0.; gen_totchg_sube0_pt3 = 0.; gen_totchg_sube0_pt4 = 0.; gen_totchg_sube0_pt5 = 0.;
                  gen_bkgchg_sube0_pt2 = 0.; gen_bkgchg_sube0_pt3 = 0.; gen_bkgchg_sube0_pt4 = 0.; gen_bkgchg_sube0_pt5 = 0.;
 
                  for(int gen =0; gen < (int) genPt->size(); gen++){
                    	if(!passGenTrackCuts(genPt->at(gen), genEta->at(gen), gen_chg->at(gen))) continue;
                        double dr = sqrt(pow(genEta->at(gen) - calo_jteta, 2) + pow(acos(cos(calo_jtphi - genPhi->at(gen))),2));
                        double refl_dr = sqrt(pow(genEta->at(gen) + calo_jteta, 2) + pow(acos(cos(calo_jtphi - genPhi->at(gen))),2));

                        if(dr < 0.4){
                            n_gen++;    
                            h_trk_pt_gen[mycbin]->Fill(genPt->at(gen),weight);

                            if(genPt->at(gen)>1.) n_genpt1++;
                            if(genPt->at(gen)>2.) {
                                n_genpt2++;
                                gen_totchg_pt2 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                            }
                            if(genPt->at(gen)>3.) {
                                n_genpt3++;
                                gen_totchg_pt3 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                            }
                            if(genPt->at(gen)>4.) {
                                n_genpt4++;
                                gen_totchg_pt4 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                            }
                            if(genPt->at(gen)>5.) {
                                n_genpt5++;
                                gen_totchg_pt5 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                            }

                            if(do_sube0 && sube->at(gen)==0){ 
                              if(genPt->at(gen)>2.) {
                                n_sube0_gen++;
                                gen_totchg_sube0_pt2 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                              }
                              if(genPt->at(gen)>3.) gen_totchg_sube0_pt3 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                              if(genPt->at(gen)>4.) gen_totchg_sube0_pt4 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                              if(genPt->at(gen)>5.) gen_totchg_sube0_pt5 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                            }
                        }
                        else if(refl_dr < 0.4){
                            n_bkggen++;
                            if(genPt->at(gen)>1.) n_bkggenpt1++;
                            if(genPt->at(gen)>2.) {
                                n_bkggenpt2++;
                                gen_bkgchg_pt2 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                            }
                            if(genPt->at(gen)>3.) {
                                n_bkggenpt3++;
                                gen_bkgchg_pt3 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                            }
                            if(genPt->at(gen)>4.) {
                                n_bkggenpt4++;
                                gen_bkgchg_pt4 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                            }
                            if(genPt->at(gen)>5.) {
                                n_bkggenpt5++;
                                gen_bkgchg_pt5 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                            }
                            if(do_sube0 && sube->at(gen)==0){ 
                              if(genPt->at(gen)>2.) {
                                gen_bkgchg_sube0_pt2 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                                n_bkg_sube0_gen++;
                              }
                              if(genPt->at(gen)>3.) gen_bkgchg_sube0_pt3 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                              if(genPt->at(gen)>4.) gen_bkgchg_sube0_pt4 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                              if(genPt->at(gen)>5.) gen_bkgchg_sube0_pt5 += pow((genPt->at(gen)/calo_corrpt),0.5)*gen_chg->at(gen);
                            }
 
                        }
                  }
                }
                
                //ntrk and chg
                n_trk = 0; n_bkgtrk = 0; trk_totchg = 0.;
                n_trkpt1 = 0; n_trkpt2 = 0;n_trkpt3 = 0; n_trkpt4 = 0; n_trkpt5 = 0;
                n_bkgtrkpt1 = 0; n_bkgtrkpt2 = 0; n_bkgtrkpt3 = 0; n_bkgtrkpt4 = 0; n_bkgtrkpt5 = 0;
                trk_totchg_pt2 = 0.; trk_totchg_pt3 = 0.; trk_totchg_pt4 = 0.; trk_totchg_pt5 = 0.;
                trk_bkgchg_pt2 = 0.; trk_bkgchg_pt3 = 0.; trk_bkgchg_pt4 = 0.; trk_bkgchg_pt5 = 0.;

                trk_totchg_pt2_k1 = 0.; trk_totchg_pt2_k1p5 = 0.; trk_totchg_pt2_k2 = 0.; trk_totchg_pt2_k2p5 = 0.;
                trk_bkgchg_pt2_k1 = 0.; trk_bkgchg_pt2_k1p5 = 0.; trk_bkgchg_pt2_k2 = 0.; trk_bkgchg_pt2_k2p5 = 0.;
/*
                trk_totchg_pt2_trkuncup = 0.; trk_totchg_pt3_trkuncup = 0.; trk_totchg_pt4_trkuncup = 0.; trk_totchg_pt5_trkuncup = 0.;
                trk_bkgchg_pt2_trkuncup = 0.; trk_bkgchg_pt3_trkuncup = 0.; trk_bkgchg_pt4_trkuncup = 0.; trk_bkgchg_pt5_trkuncup = 0.;

                trk_totchg_pt2_trkuncdn = 0.; trk_totchg_pt3_trkuncdn = 0.; trk_totchg_pt4_trkuncdn = 0.; trk_totchg_pt5_trkuncdn = 0.;
                trk_bkgchg_pt2_trkuncdn = 0.; trk_bkgchg_pt3_trkuncdn = 0.; trk_bkgchg_pt4_trkuncdn = 0.; trk_bkgchg_pt5_trkuncdn = 0.;

                trk_totchg_pt2_trkuncpos = 0.; trk_totchg_pt3_trkuncpos = 0.; trk_totchg_pt4_trkuncpos = 0.; trk_totchg_pt5_trkuncpos = 0.;
                trk_bkgchg_pt2_trkuncpos = 0.; trk_bkgchg_pt3_trkuncpos = 0.; trk_bkgchg_pt4_trkuncpos = 0.; trk_bkgchg_pt5_trkuncpos = 0.;
*/
                highest_idx_pt = 0.; highest_idx_2_pt = 0.;
                highest_idx = 0; highest_idx_2 = 0;

                //search for highest pt track
                for(int itrk =0; itrk < nTrk; itrk++){
                    if(!passTrackCuts(is_pp, 1, trkPt[itrk], trkEta[itrk], highPurity[itrk], trkChi2[itrk], trkNdof[itrk], trkNlayer[itrk], trkNHit[itrk], pfHcal[itrk], pfEcal[itrk], trkDxy1[itrk]/trkDxyError1[itrk], trkDz1[itrk]/trkDzError1[itrk])) continue;
                    double dr = sqrt(pow(trkEta[itrk] - calo_jteta, 2) + pow(acos(cos(calo_jtphi - trkPhi[itrk])),2));
                    if(dr > 0.4) continue;
                    if(trkPt[itrk] > highest_idx_pt){
                      highest_idx_pt = trkPt[itrk];
                      highest_idx = itrk;
                    }
                }
                    
                //search for subleading trk
                for(int itrk =0; itrk < nTrk; itrk++){
                    if(itrk == highest_idx) continue ;
                    if(!passTrackCuts(is_pp, 1, trkPt[itrk], trkEta[itrk], highPurity[itrk], trkChi2[itrk], trkNdof[itrk], trkNlayer[itrk], trkNHit[itrk], pfHcal[itrk], pfEcal[itrk], trkDxy1[itrk]/trkDxyError1[itrk], trkDz1[itrk]/trkDzError1[itrk])) continue;
                    double dr = sqrt(pow(trkEta[itrk] - calo_jteta, 2) + pow(acos(cos(calo_jtphi - trkPhi[itrk])),2));
                    if(dr > 0.4) continue;
                    if(trkPt[itrk] > highest_idx_2_pt){
                      highest_idx_2_pt = trkPt[itrk];
                      highest_idx_2 = itrk;
                    }
                }

                trk_totchg_pt2_highest = -999.; trk_totchg_pt2_highest_2 = -999.;
                if(highest_idx_pt > 0.) trk_totchg_pt2_highest = pow((trkPt[highest_idx]/calo_corrpt),0.5)*trk_chg[highest_idx];
                if(highest_idx_2_pt > 0.) {
                    trk_totchg_pt2_highest_2 = trk_totchg_pt2_highest; 
                    trk_totchg_pt2_highest_2 += pow((trkPt[highest_idx_2]/calo_corrpt),0.5)*trk_chg[highest_idx_2];
                }

                for(int itrk =0; itrk < nTrk; itrk++){
                	if(!passTrackCuts(is_pp, 1, trkPt[itrk], trkEta[itrk], highPurity[itrk], trkChi2[itrk], trkNdof[itrk], trkNlayer[itrk], trkNHit[itrk], pfHcal[itrk], pfEcal[itrk], trkDxy1[itrk]/trkDxyError1[itrk], trkDz1[itrk]/trkDzError1[itrk])) continue;
                	//trkPt[itrk] = trkPt[itrk]*trk_correction;
                    double dr = sqrt(pow(trkEta[itrk] - calo_jteta, 2) + pow(acos(cos(calo_jtphi - trkPhi[itrk])),2));
                    double refl_dr = sqrt(pow(trkEta[itrk] + calo_jteta, 2) + pow(acos(cos(calo_jtphi - trkPhi[itrk])),2));

                    //float rmin = 999;
                    //trk_correction = trkCorr->getTrkCorr(trkPt[itrk],trkEta[itrk],trkPhi[itrk],hiBin,rmin); 
                    trk_correction = trkCorr->getTrkCorr(trkPt[itrk],trkEta[itrk],trkPhi[itrk],hiBin); 
                    //if(trk_chg[itrk]==1 && trkPt[itrk]<4.) trk_correction = 1.01*trk_correction;
                    //else if (trk_chg[itrk]==1 && trkPt[itrk]>4.) trk_correction = 1.005*trk_correction;

                    if(dr < 0.4){
                    
                        if(trkPt[itrk]>highest_idx_pt){ 
                            highest_idx_pt = trkPt[itrk];
                            highest_idx_2_pt = trkPt[itrk];
                        }
                        else if(trkPt[itrk]>highest_idx_2_pt) highest_idx_2_pt = trkPt[itrk];

                        h_trk_pt_reco[mycbin]->Fill(trkPt[itrk],trk_correction*weight);

                    	n_trk++;                        

                        if(trkPt[itrk]>1.) n_trkpt1++;
                        if(trkPt[itrk]>2.) {
                            n_trkpt2++;

                            trk_totchg_pt2 += trk_correction*pow((trkPt[itrk]/calo_corrpt),0.5)*trk_chg[itrk];
                            h_trk_pt2_reco[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),0.5)*trk_correction*weight);

                            h_trk_ptk1_reco[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),1.)*trk_correction*weight);
                            h_trk_ptk1p5_reco[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),1.5)*trk_correction*weight);
                            h_trk_ptk2_reco[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),2.)*trk_correction*weight);
                            h_trk_ptk2p5_reco[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),2.5)*trk_correction*weight);

                            trk_totchg_pt2_k1 += trk_correction*pow((trkPt[itrk]/calo_corrpt),1.)*trk_chg[itrk];
                            trk_totchg_pt2_k1p5 += trk_correction*pow((trkPt[itrk]/calo_corrpt),1.5)*trk_chg[itrk];
                            trk_totchg_pt2_k2 += trk_correction*pow((trkPt[itrk]/calo_corrpt),2.)*trk_chg[itrk];
                            trk_totchg_pt2_k2p5 += trk_correction*pow((trkPt[itrk]/calo_corrpt),2.5)*trk_chg[itrk];
                        }
                        if(trkPt[itrk]>3.) {
                            n_trkpt3++;
                            trk_totchg_pt3 += trk_correction*pow((trkPt[itrk]/calo_corrpt),0.5)*trk_chg[itrk];
                        }
                        if(trkPt[itrk]>4.) {
                            n_trkpt4++;
                            h_trk_pt4_reco[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),0.5)*trk_correction*weight);
                            trk_totchg_pt4 += trk_correction*pow((trkPt[itrk]/calo_corrpt),0.5)*trk_chg[itrk];
                        }
                        if(trkPt[itrk]>5.) {
                            n_trkpt5++;
                            h_trk_pt5_reco[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),0.5)*trk_correction*weight);
                            trk_totchg_pt5 += trk_correction*pow((trkPt[itrk]/calo_corrpt),0.5)*trk_chg[itrk];
                        }
                    }
                    else if(refl_dr < 0.4){
                    	n_bkgtrk++;
                        h_trk_pt_reco_refl[mycbin]->Fill(trkPt[itrk],trk_correction*weight);

                        if(trkPt[itrk]>1.) n_bkgtrkpt1++;
                        if(trkPt[itrk]>2.) {
                            n_bkgtrkpt2++;
                            trk_bkgchg_pt2 += trk_correction*pow((trkPt[itrk]/calo_corrpt),0.5)*trk_chg[itrk];
                            h_trk_pt2_reco_refl[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),0.5)*trk_correction*weight);

                            h_trk_ptk1_reco_refl[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),1.)*trk_correction*weight);
                            h_trk_ptk1p5_reco_refl[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),1.5)*trk_correction*weight);
                            h_trk_ptk2_reco_refl[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),2.)*trk_correction*weight);
                            h_trk_ptk2p5_reco_refl[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),2.5)*trk_correction*weight);

                            trk_bkgchg_pt2_k1 += trk_correction*pow((trkPt[itrk]/calo_corrpt),1.)*trk_chg[itrk];
                            trk_bkgchg_pt2_k1p5 += trk_correction*pow((trkPt[itrk]/calo_corrpt),1.5)*trk_chg[itrk];
                            trk_bkgchg_pt2_k2 += trk_correction*pow((trkPt[itrk]/calo_corrpt),2.)*trk_chg[itrk];
                            trk_bkgchg_pt2_k2p5 += trk_correction*pow((trkPt[itrk]/calo_corrpt),2.5)*trk_chg[itrk];
                        }
                        if(trkPt[itrk]>3.) {
                            n_bkgtrkpt3++;
                            trk_bkgchg_pt3 += trk_correction*pow((trkPt[itrk]/calo_corrpt),0.5)*trk_chg[itrk];
                        }
                        if(trkPt[itrk]>4.) {
                            n_bkgtrkpt4++;
                            h_trk_pt4_reco_refl[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),0.5)*trk_correction*weight);
                            trk_bkgchg_pt4 += trk_correction*pow((trkPt[itrk]/calo_corrpt),0.5)*trk_chg[itrk];
                        }
                        if(trkPt[itrk]>5.) {
                            n_bkgtrkpt5++;
                            h_trk_pt5_reco_refl[mycbin]->Fill(trkPt[itrk],pow((trkPt[itrk]/calo_corrpt),0.5)*trk_correction*weight);
                            trk_bkgchg_pt5 += trk_correction*pow((trkPt[itrk]/calo_corrpt),0.5)*trk_chg[itrk];
                        }
                    }
                }

                //----------------------------------------------------
				//      EVENT MIXING STARTS HERE!  (For DATA ONLY)
				//-----------------------------------------------------

				if(do_mixing){
                        mixStopwatch.Start(0);
						jet_cent = 0;
						if(!is_pp){ jet_cent = centbins->FindBin(hiBin)-1;}
						jet_vzbin = vzbins->FindBin(vz)-1;

						unsigned int pct = gRandom->Rndm()*mixing_lists.at(jet_vzbin).at(jet_cent).size();
			            //cout<<"pct: "<<pct<<"  mixing_lists.at(jet_vzbin).at(jet_cent).size()"<<mixing_lists.at(jet_vzbin).at(jet_cent).size()<<endl;
						bool circleCheck=false;
						
						unsigned int mevi = pct;
						unsigned int dummy=0;																	

                        while(dummy < meptrig){
                        	mixStopwatch_p1.Start(0);
							if((int)mixing_lists[jet_vzbin][jet_cent].size()<2){
								cout << "Warning! vzbin: "<< jet_vzbin-1 << ", cent bin: "<< jet_cent-1 << " is devoid of mixing events!" << endl;
								break;
							}
							if(mevi>mixing_lists.at(jet_vzbin).at(jet_cent).size()){ mevi=1; circleCheck=true; }
					
							if(mevi==0) mevi++;
							
							if(mixing_lists.at(jet_vzbin).at(jet_cent).at(mevi-1) == evi){ 
								mevi++; continue; 
							}
							mixStopwatch_p1.Stop();

                            mixStopwatch_p2.Start(0);
							entryCalls++;							

                            //mixing_event_tree->GetEntry(mixing_lists.at(jet_vzbin).at(jet_cent).at(mevi-1));
                            mixing_track_tree->GetEntry(mixing_lists.at(jet_vzbin).at(jet_cent).at(mevi-1));

                            dummy++; mevi++;
							mixStopwatch_p2.Stop();

                            mixStopwatch_p3.Start(0);
							unsigned int me_cent=0;
							if(!is_pp){ me_cent = centbins->FindBin(me_hiBin)-1;}
							unsigned int me_vzbin = vzbins->FindBin(me_vz)-1;

							//cout << "jet vzbin : "<< jet_vzbin << " me vzbin: "<< me_vzbin << endl;
							//cout << "jet centbin : "<< jet_cent << " me centbin: "<< me_cent << endl;
							//make sure i'm doing this thing correctly...
							assert(me_vzbin==jet_vzbin);
							assert(me_cent==jet_cent);
                            mixStopwatch_p3.Stop();

                            //ntrk and chg
                            n_mixtrk = 0;
                            trk_mixchg_pt2 = 0.; trk_mixchg_pt3 = 0.; trk_mixchg_pt4 = 0.; trk_mixchg_pt5 = 0.;
                            //cout<<"nmixTrk: "<<nmixTrk<<endl;
                            for(int itrk =0; itrk < nmixTrk; itrk++){
                               //cout<<itrk<<endl; 
                               if(!passTrackCuts(is_pp, 1, me_trkPt[itrk], me_trkEta[itrk], me_highPurity[itrk], me_trkChi2[itrk], me_trkNdof[itrk], me_trkNlayer[itrk], me_trkNHit[itrk], me_pfHcal[itrk], me_pfEcal[itrk], me_trkDxy[itrk]/me_trkDxyError[itrk], me_trkDz[itrk]/me_trkDzError[itrk])) continue;
                               //trkPt[itrk] = trkPt[itrk]*trk_correction;
                               double dr = sqrt(pow(me_trkEta[itrk] - calo_jteta, 2) + pow(acos(cos(calo_jtphi - me_trkPhi[itrk])),2));
                               if(dr < 0.4){
                    	         n_mixtrk++;

                             	 //float rmin = 999;
                           	     //trk_correction = trkCorr->getTrkCorr(me_trkPt[itrk],me_trkEta[itrk],me_trkPhi[itrk],hiBin,rmin);
                                 trk_correction = trkCorr->getTrkCorr(trkPt[itrk],trkEta[itrk],trkPhi[itrk],hiBin);

                                 if(me_trkPt[itrk]>2.) trk_mixchg_pt2 += trk_correction*pow((me_trkPt[itrk]/calo_corrpt),0.5)*me_trkchg[itrk];
                                 if(me_trkPt[itrk]>3.) trk_mixchg_pt3 += trk_correction*pow((me_trkPt[itrk]/calo_corrpt),0.5)*me_trkchg[itrk];
                                 if(me_trkPt[itrk]>4.) trk_mixchg_pt4 += trk_correction*pow((me_trkPt[itrk]/calo_corrpt),0.5)*me_trkchg[itrk];
                                 if(me_trkPt[itrk]>5.) trk_mixchg_pt5 += trk_correction*pow((me_trkPt[itrk]/calo_corrpt),0.5)*me_trkchg[itrk];
                               }

                            }//trk loop  
                            //cout<<"mevi: "<<mevi<<endl;       
                        }//meptrig events per trigger
				
                        mixStopwatch.Stop();

				} ///finished mixing 

				///// Fill it
		    	ftree->Fill();
			
			} /// calo jet loop

		}  ///event loop

        evtStopwatch.Stop();
        cout << "Avg event time: " << evtStopwatch.RealTime()/(double)n_evt << " sec/evt" << endl;
        cout << "Avg mix time: "<< mixStopwatch.RealTime()/(double)n_evt << " sec/evt" << endl;
        cout << "mix time, part1: "<< mixStopwatch_p1.RealTime()/(double)n_evt << " sec/evt" << endl;
        cout << "mix time, part2: "<< mixStopwatch_p2.RealTime()/(double)n_evt << " sec/evt" << endl;
        cout << "mix time, part3: "<< mixStopwatch_p3.RealTime()/(double)n_evt << " sec/evt" << endl;
        cout << "total entry calls: " << entryCalls << endl;

		my_file->Close();
	
	} //file loop

	cout<<"writing"<<endl;

	output_file->cd();
	ftree->Write();

    for(int ibin=0; ibin<nCbins; ibin++){
        h_trk_pt_reco[ibin]->Write();
        h_trk_pt_reco_full[ibin]->Write();
        h_trk_pt2_reco[ibin]->Write();
        h_trk_pt4_reco[ibin]->Write();
        h_trk_pt5_reco[ibin]->Write();

        h_trk_pt_reco_refl[ibin]->Write();
        h_trk_pt2_reco_refl[ibin]->Write();
        h_trk_pt4_reco_refl[ibin]->Write();
        h_trk_pt5_reco_refl[ibin]->Write();

        h_trk_pt_reco_pos[ibin]->Write();
        h_trk_pt_reco_neg[ibin]->Write();
        h_trk_eta_reco[ibin]->Write();

        h_trk_ptk1_reco[ibin]->Write();
        h_trk_ptk1p5_reco[ibin]->Write();
        h_trk_ptk2_reco[ibin]->Write();
        h_trk_ptk2p5_reco[ibin]->Write();

        h_trk_ptk1_reco_refl[ibin]->Write();
        h_trk_ptk1p5_reco_refl[ibin]->Write();
        h_trk_ptk2_reco_refl[ibin]->Write();
        h_trk_ptk2p5_reco_refl[ibin]->Write();

        h_trk_pt_gen[ibin]->Write();
        h_trk_pt_gen_full[ibin]->Write();
        h_trk_pt4_gen[ibin]->Write();
        h_trk_pt5_gen[ibin]->Write();
        h_trk_pt_gen_pos[ibin]->Write();
        h_trk_pt_gen_neg[ibin]->Write();
        h_trk_eta_gen[ibin]->Write();
    }

    h_vz->Write();
    h_pthat->Write();
    h_hibin->Write();

	output_file->Close();

	cout<<"done"<<endl;

}
