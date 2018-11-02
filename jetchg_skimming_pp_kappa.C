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
//#include "xiaoTrkCorr/xiaoTrkCorr.h"
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
void jetchg_skimming_pp_kappa(bool doCrab=0, int jobID=0, int endfile = 105, int dataset_type_code = 11, int output_file_num = 1)
{
	bool is_pp = true;
	bool do_mixing = false;
    bool do_sube0 = false;

    nCScorr *corrpt = new nCScorr(is_pp,1);

    TrkCorr* trkCorr = new TrkCorr("TrkCorr_July22_Iterative_pp_eta2p4/");
    //xiaoTrkCorr* trkCorr = new xiaoTrkCorr("xiaoTrkCorr/eta_symmetry_cymbalCorr_FineBin.root");

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
    f_res_gauss_lowpt->FixParameter(0,0.035);
    TF1 *f_res_gauss_highpt = new TF1("f_res_gauss_highpt", "exp(-((x-1.)^2)/(2*([0]^2)))", 0., 2.);
    f_res_gauss_highpt->FixParameter(0,0.03);

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
   
    double weight, pthat_weight, weight_vz, weight_cen;    

	float calo_jteta, calo_jtphi, calo_jtpt, calo_corrpt;
	float calo_refpt, calo_refeta, calo_refphi; 

	float trk_totchg_pt2_k1_pos, trk_totchg_pt2_k1p5_pos, trk_totchg_pt2_k2_pos, trk_totchg_pt2_k1_neg, trk_totchg_pt2_k1p5_neg, trk_totchg_pt2_k2_neg;
    float trk_bkgchg_pt2_k1_pos, trk_bkgchg_pt2_k1p5_pos, trk_bkgchg_pt2_k2_pos, trk_bkgchg_pt2_k1_neg, trk_bkgchg_pt2_k1p5_neg, trk_bkgchg_pt2_k2_neg;

    float trk_totchg_pt2_k0p3_pos, trk_totchg_pt2_k0p5_pos, trk_totchg_pt2_k0p7_pos, trk_totchg_pt2_k0p3_neg, trk_totchg_pt2_k0p5_neg, trk_totchg_pt2_k0p7_neg;
    float trk_bkgchg_pt2_k0p3_pos, trk_bkgchg_pt2_k0p5_pos, trk_bkgchg_pt2_k0p7_pos, trk_bkgchg_pt2_k0p3_neg, trk_bkgchg_pt2_k0p5_neg, trk_bkgchg_pt2_k0p7_neg;

    float trk_totchg_pt4_k0p5_pos, trk_totchg_pt5_k0p5_pos, trk_totchg_pt4_k0p5_neg, trk_totchg_pt5_k0p5_neg;
    float trk_bkgchg_pt4_k0p5_pos, trk_bkgchg_pt5_k0p5_pos, trk_bkgchg_pt4_k0p5_neg, trk_bkgchg_pt5_k0p5_neg;

	int calo_parton_flavor;
	double trk_correction;

	Int_t HBHENoiseFilter, HBHENoiseFilterResultRun2Loose, eventSelection, pvFilter, phfCoincFilter3, pclusterCompatibilityFilter, pBeamScrapingFilter;
	Int_t HLT_Jet80, HLT_Jet100, HLT_Jet80_Prescl, HLT_Jet100_Prescl, L1_Jet80_Prescl, L1_Jet100_Prescl;
	Int_t pPAcollisionEventSelectionPA = -999;
	Int_t hiBin = -999;
	Float_t pthat = -999;
	Float_t vz = -999, hiHF= -999.;
	Int_t calo_nref, nTrk, nPFpart, nCSpart, nmixTrk;
	Int_t ncs_id145_pt2;  
	Int_t npf_id1_pt2;
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

    ftree->Branch("trk_totchg_pt2_k0p3_pos",&trk_totchg_pt2_k0p3_pos);
    ftree->Branch("trk_totchg_pt2_k0p5_pos",&trk_totchg_pt2_k0p5_pos);
    ftree->Branch("trk_totchg_pt2_k0p7_pos",&trk_totchg_pt2_k0p7_pos);
    ftree->Branch("trk_totchg_pt2_k0p3_neg",&trk_totchg_pt2_k0p3_neg);
    ftree->Branch("trk_totchg_pt2_k0p5_neg",&trk_totchg_pt2_k0p5_neg);
    ftree->Branch("trk_totchg_pt2_k0p7_neg",&trk_totchg_pt2_k0p7_neg);

    ftree->Branch("trk_totchg_pt4_k0p5_pos",&trk_totchg_pt4_k0p5_pos);
    ftree->Branch("trk_totchg_pt5_k0p5_pos",&trk_totchg_pt5_k0p5_pos);
    ftree->Branch("trk_totchg_pt4_k0p5_neg",&trk_totchg_pt4_k0p5_neg);
    ftree->Branch("trk_totchg_pt5_k0p5_neg",&trk_totchg_pt5_k0p5_neg);

    ftree->Branch("trk_totchg_pt2_k1_pos",&trk_totchg_pt2_k1_pos);
    ftree->Branch("trk_totchg_pt2_k1p5_pos",&trk_totchg_pt2_k1p5_pos);
    ftree->Branch("trk_totchg_pt2_k2_pos",&trk_totchg_pt2_k2_pos);
    ftree->Branch("trk_totchg_pt2_k1_neg",&trk_totchg_pt2_k1_neg);
    ftree->Branch("trk_totchg_pt2_k1p5_neg",&trk_totchg_pt2_k1p5_neg);
    ftree->Branch("trk_totchg_pt2_k2_neg",&trk_totchg_pt2_k2_neg);

    ftree->Branch("trk_bkgchg_pt2_k0p3_pos",&trk_bkgchg_pt2_k0p3_pos);
    ftree->Branch("trk_bkgchg_pt2_k0p5_pos",&trk_bkgchg_pt2_k0p5_pos);
    ftree->Branch("trk_bkgchg_pt2_k0p7_pos",&trk_bkgchg_pt2_k0p7_pos);
    ftree->Branch("trk_bkgchg_pt2_k0p3_neg",&trk_bkgchg_pt2_k0p3_neg);
    ftree->Branch("trk_bkgchg_pt2_k0p5_neg",&trk_bkgchg_pt2_k0p5_neg);
    ftree->Branch("trk_bkgchg_pt2_k0p7_neg",&trk_bkgchg_pt2_k0p7_neg);

    ftree->Branch("trk_bkgchg_pt2_k1_pos",&trk_bkgchg_pt2_k1_pos);
    ftree->Branch("trk_bkgchg_pt2_k1p5_pos",&trk_bkgchg_pt2_k1p5_pos);
    ftree->Branch("trk_bkgchg_pt2_k2_pos",&trk_bkgchg_pt2_k2_pos);
    ftree->Branch("trk_bkgchg_pt2_k1_neg",&trk_bkgchg_pt2_k1_neg);
    ftree->Branch("trk_bkgchg_pt2_k1p5_neg",&trk_bkgchg_pt2_k1p5_neg);
    ftree->Branch("trk_bkgchg_pt2_k2_neg",&trk_bkgchg_pt2_k2_neg);

    ftree->Branch("trk_bkgchg_pt4_k0p5_pos",&trk_bkgchg_pt4_k0p5_pos);
    ftree->Branch("trk_bkgchg_pt5_k0p5_pos",&trk_bkgchg_pt5_k0p5_pos);
    ftree->Branch("trk_bkgchg_pt4_k0p5_neg",&trk_bkgchg_pt4_k0p5_neg);
    ftree->Branch("trk_bkgchg_pt5_k0p5_neg",&trk_bkgchg_pt5_k0p5_neg);
 
    if(do_PbPb){
        ftree->Branch("ncs_id145_pt2",&ncs_id145_pt2);
	}
	else{
	    ftree->Branch("npf_id1_pt2",&npf_id1_pt2);
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
	UChar_t trkAlgo[MAXPARTICLES], trkNHit[MAXPARTICLES], trkNlayer[MAXPARTICLES], trkNdof[MAXPARTICLES]; //Run2
	Bool_t highPurity[MAXPARTICLES];
    Int_t trk_chg[MAXPARTICLES];

	Float_t me_trkPt[MAXPARTICLES], me_trkEta[MAXPARTICLES], me_trkPhi[MAXPARTICLES], me_trkDxy[MAXPARTICLES], me_trkDxyError[MAXPARTICLES], me_trkDz[MAXPARTICLES], me_trkDzError[MAXPARTICLES], me_trkPtError[MAXPARTICLES], me_pfEcal[MAXPARTICLES], me_pfHcal[MAXPARTICLES], me_trkChi2[MAXPARTICLES];
	UChar_t me_trkAlgo[MAXPARTICLES], me_trkNHit[MAXPARTICLES], me_trkNlayer[MAXPARTICLES], me_trkNdof[MAXPARTICLES]; //Run2
	Bool_t me_highPurity[MAXPARTICLES];
    Int_t me_trkchg[MAXPARTICLES];

	vector<int> *gen_chg=0, *sube=0, *t_sta=0;;
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
        //TFile *my_file = TFile::Open("HiForestAOD_Pythia6.root");

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
            if(do_PbPb)inp_tree->SetBranchAddress("sube", &sube);
            else if(!do_PbPb)inp_tree->SetBranchAddress("sta", &t_sta); 
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

		for(int evi = 0; evi < n_evt; evi++) {
		//for(int evi = 0; evi < 100; evi++) {
			
            inp_tree->GetEntry(evi);
			inp_tree2->GetEntry(evi);
			inp_tree3->GetEntry(evi);
			if(do_PbPb) inp_tree_CS->GetEntry(evi);

            if(fabs(vz) >= 15.) continue;
			
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
                while(pthat>pthatbins[ibin+1]) ibin++;
                pthat_weight = (xsecs[ibin]-xsecs[ibin+1])/pthatEntries[ibin];

                //pythia6 
                if(is_pp) weight_vz = 1./fvz_pp->Eval(vz);
                //Pythia6+Hydjet
                else weight_vz = fvz_cymbal->Eval(vz);

                if(do_PbPb) weight_cen = fcent_cymbal(hiBin,f_cent);
                else weight_cen = 1.;

                weight = pthat_weight*weight_vz*weight_cen;
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
				    	if(pfCandId->at(icand)==1){
					    	npf_id1_pt2++;
				        }
			    	}
			    	calo_corrpt = corrpt->getCorrection(is_pp, npf_id1_pt2, 0, calo_jtpt, calo_jteta);
			    }
			    else{
    				ncs_id145_pt2 = 0;

	    			for(int icand=0; icand<nCSpart; icand++){
		    			if(abs(csCandEta->at(icand))>2.4) continue;
			    		double dr = sqrt(pow(csCandEta->at(icand) - calo_jteta, 2) + pow(acos(cos(calo_jtphi - csCandPhi->at(icand))),2));
				    	if(dr > 0.4 || csCandPt->at(icand)<2.) continue;
				    	if(csCandId->at(icand)==1 || csCandId->at(icand)==4 || csCandId->at(icand)==5){
					    	ncs_id145_pt2++;
				        }
			    	}
			    	//calo_corrpt = corrpt->getCorrection(is_pp, ncs_id145_pt2, hiBin, calo_jtpt, calo_jteta);
                    if(ncs_id145_pt2 > 20 && calo_jtpt <= 120.) calo_corrpt = corrpt->getCorrection(is_pp, 20, hiBin, calo_jtpt, calo_jteta);
                    else calo_corrpt = corrpt->getCorrection(is_pp, ncs_id145_pt2, hiBin, calo_jtpt, calo_jteta);
			    }	
                
                if(calo_corrpt <= 0.) continue;

                //ntrk and chg
                n_trk = 0; n_bkgtrk = 0; 
                n_trkpt1 = 0; n_trkpt2 = 0;n_trkpt3 = 0; n_trkpt4 = 0; n_trkpt5 = 0;
                n_bkgtrkpt1 = 0; n_bkgtrkpt2 = 0; n_bkgtrkpt3 = 0; n_bkgtrkpt4 = 0; n_bkgtrkpt5 = 0;

                trk_totchg_pt2_k1_pos = 0.; trk_totchg_pt2_k1p5_pos = 0.; trk_totchg_pt2_k2_pos = 0.;
                trk_bkgchg_pt2_k1_pos = 0.; trk_bkgchg_pt2_k1p5_pos = 0.; trk_bkgchg_pt2_k2_pos = 0.;
                trk_totchg_pt2_k1_neg = 0.; trk_totchg_pt2_k1p5_neg = 0.; trk_totchg_pt2_k2_neg = 0.;
                trk_bkgchg_pt2_k1_neg = 0.; trk_bkgchg_pt2_k1p5_neg = 0.; trk_bkgchg_pt2_k2_neg = 0.;

                trk_totchg_pt2_k0p3_pos = 0.; trk_totchg_pt2_k0p5_pos = 0.; trk_totchg_pt2_k0p7_pos = 0.;
                trk_bkgchg_pt2_k0p3_pos = 0.; trk_bkgchg_pt2_k0p5_pos = 0.; trk_bkgchg_pt2_k0p7_pos = 0.;
                trk_totchg_pt2_k0p3_neg = 0.; trk_totchg_pt2_k0p5_neg = 0.; trk_totchg_pt2_k0p7_neg = 0.;
                trk_bkgchg_pt2_k0p3_neg = 0.; trk_bkgchg_pt2_k0p5_neg = 0.; trk_bkgchg_pt2_k0p7_neg = 0.;

                trk_totchg_pt4_k0p5_pos = 0.; trk_totchg_pt5_k0p5_pos = 0.;
                trk_bkgchg_pt4_k0p5_pos = 0.; trk_bkgchg_pt5_k0p5_pos = 0.;
                trk_totchg_pt4_k0p5_neg = 0.; trk_totchg_pt5_k0p5_neg = 0.;
                trk_bkgchg_pt4_k0p5_neg = 0.; trk_bkgchg_pt5_k0p5_neg = 0.;

                for(int itrk =0; itrk < nTrk; itrk++){
                	if(!passTrackCuts(is_pp, 1, trkPt[itrk], trkEta[itrk], highPurity[itrk], trkChi2[itrk], trkNdof[itrk], trkNlayer[itrk], trkNHit[itrk], pfHcal[itrk], pfEcal[itrk], trkDxy1[itrk]/trkDxyError1[itrk], trkDz1[itrk]/trkDzError1[itrk])) continue;
                    double dr = sqrt(pow(trkEta[itrk] - calo_jteta, 2) + pow(acos(cos(calo_jtphi - trkPhi[itrk])),2));
                    double refl_dr = sqrt(pow(trkEta[itrk] + calo_jteta, 2) + pow(acos(cos(calo_jtphi - trkPhi[itrk])),2));

                    float rmin = 999.;  
                    trk_correction = trkCorr->getTrkCorr(trkPt[itrk],trkEta[itrk],trkPhi[itrk],1,rmin); 

                    if(dr < 0.4){
                        n_trkpt2++;
                        if(trk_chg[itrk] == 1){
                            trk_totchg_pt2_k0p3_pos += trk_correction*pow(trkPt[itrk],0.3)*trk_chg[itrk];
                            trk_totchg_pt2_k0p5_pos += trk_correction*pow(trkPt[itrk],0.5)*trk_chg[itrk];
                            trk_totchg_pt2_k0p7_pos += trk_correction*pow(trkPt[itrk],0.7)*trk_chg[itrk];

                            trk_totchg_pt2_k1_pos += trk_correction*pow(trkPt[itrk],1.)*trk_chg[itrk];
                            trk_totchg_pt2_k1p5_pos += trk_correction*pow(trkPt[itrk],1.5)*trk_chg[itrk];
                            trk_totchg_pt2_k2_pos += trk_correction*pow(trkPt[itrk],2.)*trk_chg[itrk];
                            if(trkPt[itrk] > 4.){
                            	n_trkpt4++;
                                trk_totchg_pt4_k0p5_pos += trk_correction*pow(trkPt[itrk],0.5)*trk_chg[itrk];
                            }
                            if(trkPt[itrk] > 5.){
                            	n_trkpt5++;
                                trk_totchg_pt5_k0p5_pos += trk_correction*pow(trkPt[itrk],0.5)*trk_chg[itrk];
                            }
                        }                        
                        else if(trk_chg[itrk] == -1){
                            trk_totchg_pt2_k0p3_neg += trk_correction*pow(trkPt[itrk],0.3)*trk_chg[itrk];
                            trk_totchg_pt2_k0p5_neg += trk_correction*pow(trkPt[itrk],0.5)*trk_chg[itrk];
                            trk_totchg_pt2_k0p7_neg += trk_correction*pow(trkPt[itrk],0.7)*trk_chg[itrk];

                            trk_totchg_pt2_k1_neg += trk_correction*pow(trkPt[itrk],1.)*trk_chg[itrk];
                            trk_totchg_pt2_k1p5_neg += trk_correction*pow(trkPt[itrk],1.5)*trk_chg[itrk];
                            trk_totchg_pt2_k2_neg += trk_correction*pow(trkPt[itrk],2.)*trk_chg[itrk];
                            if(trkPt[itrk] > 4.){
                            	n_trkpt4++;
                                trk_totchg_pt4_k0p5_neg += trk_correction*pow(trkPt[itrk],0.5)*trk_chg[itrk];
                            }
                            if(trkPt[itrk] > 5.){
                            	n_trkpt5++;
                                trk_totchg_pt5_k0p5_neg += trk_correction*pow(trkPt[itrk],0.5)*trk_chg[itrk];
                            }
                        }
                    }
                    else if(refl_dr < 0.4){
                    	n_bkgtrk++;
                        if(trk_chg[itrk] == 1){
                            trk_bkgchg_pt2_k0p3_pos += trk_correction*pow(trkPt[itrk],0.3)*trk_chg[itrk];
                            trk_bkgchg_pt2_k0p5_pos += trk_correction*pow(trkPt[itrk],0.5)*trk_chg[itrk];
                            trk_bkgchg_pt2_k0p7_pos += trk_correction*pow(trkPt[itrk],0.7)*trk_chg[itrk];

                            trk_bkgchg_pt2_k1_pos += trk_correction*pow(trkPt[itrk],1.)*trk_chg[itrk];
                            trk_bkgchg_pt2_k1p5_pos += trk_correction*pow(trkPt[itrk],1.5)*trk_chg[itrk];
                            trk_bkgchg_pt2_k2_pos += trk_correction*pow(trkPt[itrk],2.)*trk_chg[itrk];
                            if(trkPt[itrk] > 4.){
                            	n_bkgtrkpt4++;
                                trk_bkgchg_pt4_k0p5_pos += trk_correction*pow(trkPt[itrk],0.5)*trk_chg[itrk];
                            }
                            if(trkPt[itrk] > 5.){
                                trk_bkgchg_pt5_k0p5_pos += trk_correction*pow(trkPt[itrk],0.5)*trk_chg[itrk];
                                n_bkgtrkpt5++;
                            }
                        }
                        else if(trk_chg[itrk] == -1){
                            trk_bkgchg_pt2_k0p3_neg += trk_correction*pow(trkPt[itrk],0.3)*trk_chg[itrk];
                            trk_bkgchg_pt2_k0p5_neg += trk_correction*pow(trkPt[itrk],0.5)*trk_chg[itrk];
                            trk_bkgchg_pt2_k0p7_neg += trk_correction*pow(trkPt[itrk],0.7)*trk_chg[itrk];

                            trk_bkgchg_pt2_k1_neg += trk_correction*pow(trkPt[itrk],1.)*trk_chg[itrk];
                            trk_bkgchg_pt2_k1p5_neg += trk_correction*pow(trkPt[itrk],1.5)*trk_chg[itrk];
                            trk_bkgchg_pt2_k2_neg += trk_correction*pow(trkPt[itrk],2.)*trk_chg[itrk];
                            if(trkPt[itrk] > 4.){
                                trk_bkgchg_pt4_k0p5_neg += trk_correction*pow(trkPt[itrk],0.5)*trk_chg[itrk];
                                n_bkgtrkpt4++;
                            }
                            if(trkPt[itrk] > 5.){
                                trk_bkgchg_pt5_k0p5_neg += trk_correction*pow(trkPt[itrk],0.5)*trk_chg[itrk];
                                n_bkgtrkpt4++;
                            }
                        }
                    }
                }

				///// Fill it
		    	ftree->Fill();
			
			} /// calo jet loop

		}  ///event loop

		my_file->Close();
	
	} //file loop

	cout<<"writing"<<endl;

	output_file->cd();
	ftree->Write();

    h_vz->Write();
    h_pthat->Write();
    h_hibin->Write();

	output_file->Close();

	cout<<"done"<<endl;

}
