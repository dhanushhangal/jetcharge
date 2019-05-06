#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
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

#include "assert.h"
#include <fstream>
#include "TMath.h"
#include <vector>

#include "jffcorr_lib/nCScorr.h"

using namespace std;

float trkPtCut=1.;
const double trketamaxcut = 2.4;

enum enum_dataset_types {e_Data2015, e_Data_pp, e_HydJet, e_Pythia, e_n_dataset_types};
TString dataset_type_strs[e_n_dataset_types] = {"Data2015","Data_pp","HydJet","Pythia"};

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

//arg 1 = which data set, arg 2 = output file number
void make_ntuples2(bool doCrab=0, int jobID=0, int endfile = 3251, int dataset_type_code = 2, int output_file_num = 1)
{

	bool is_data = false;

	if(dataset_type_code == 0 || dataset_type_code == 1) is_data = true;

	cout << "dataset code: " << dataset_type_code << endl;

	//-----------------------------
	// Set JFF-dependent corrections
	//-----------------------------

	float reco_eta, reco_phi, reco_pt, pfPt_temp, pfEta_temp, pfPhi_temp, pfId_temp, pfVsPt_temp;

	double corrected_pt, residual_corrected_pt, r;
	bool do_PbPb=1;

	int radius = 4;

	if(dataset_type_code== 1 || dataset_type_code== 3) do_PbPb = 0;

	cout<<"do_PbPb = "<<do_PbPb<<endl;

	nCScorr *corrpt = new nCScorr(!do_PbPb,1);

	//--------------------------------

	//////////###### centrality Bins ###########///////////////

	TTree *inp_tree;
	TTree *inp_tree2;
	TTree *inp_tree3;
	TTree *inp_tree4;
	TTree *inp_tree5;
	TTree *inp_tree6;
	TTree *inp_tree7=0;
	TTree *pftree;
	TTree *inp_tree_CS;

	string in_file_name;

	if(doCrab){
		in_file_name = Form("job_input_file_list_%d.txt",jobID);
	}
	else if(is_data&&!do_PbPb){
		in_file_name = "pp5TeV_HighPtPD_Apr2016.txt";
	}else if(is_data&&do_PbPb){
		in_file_name = "PbPbData_2015_lxplusskim_MartaAnalysis.txt";
	}else if(dataset_type_code > 10){
		in_file_name = "ppMC_Pythia6_5TeV_partonaxis.txt";
		//in_file_name = "Pythia8_ak4Calo_5TeV.txt";
	}else if(dataset_type_code > 1&&dataset_type_code <11){
		//in_file_name = "pythiaHydjet_2p76TeV_forest.txt";
		in_file_name = "Pythia6_Hydjet_PbPbMix_PurdueList.txt";
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
	TTree *mixing_tree = new TTree("mixing_tree", "");

	const int MAXPARTICLES = 100000;

	Int_t HBHENoiseFilterResultRun2Loose, eventSelection, pvFilter, phfCoincFilter3, pclusterCompatibilityFilter, pBeamScrapingFilter;
	Int_t HLT_Jet80, HLT_Jet100, HLT_Jet80_Prescl, HLT_Jet100_Prescl;
	Int_t pPAcollisionEventSelectionPA = -999;

	vector<float> calo_jteta, calo_jtphi, calo_jtpt, calo_corrpt;
	vector<float> calo_refpt, calo_refeta, calo_refphi;
	vector<int> calo_parton_flavor;
	vector<int> trkChg; //Run2
	vector<float> trkEta, trkPhi, trkPt;
	vector<int> sube, chg, pdg;
	vector<float> pt, phi, eta;

	Int_t hiBin = -999;
	Float_t pthat = -999;
	Float_t vz = -999, hiHF= -999.;//, sumpt[15];

	mixing_tree->Branch("vz", &vz);
	if(do_PbPb){
		mixing_tree->Branch("hiBin", &hiBin);
	}
	if(!is_data){
		mixing_tree->Branch("pthat", &pthat);
    }

	mixing_tree->Branch("calo_jtpt", &calo_jtpt);
	mixing_tree->Branch("calo_corrpt", &calo_corrpt);
	mixing_tree->Branch("calo_jteta", &calo_jteta);
	mixing_tree->Branch("calo_jtphi", &calo_jtphi);

	if(!is_data){
		mixing_tree->Branch("calo_refpt",&calo_refpt);
		mixing_tree->Branch("calo_refeta",&calo_refeta);
		mixing_tree->Branch("calo_refphi",&calo_refphi);
		mixing_tree->Branch("calo_refparton_flavor",&calo_parton_flavor);
	}

	mixing_tree->Branch("trkChg", &trkChg);
	mixing_tree->Branch("trkPt", &trkPt);
	mixing_tree->Branch("trkEta", &trkEta);
	mixing_tree->Branch("trkPhi", &trkPhi);

	if(!is_data){
		mixing_tree->Branch("pt", &pt);
		mixing_tree->Branch("phi", &phi);
		mixing_tree->Branch("eta", &eta);
		mixing_tree->Branch("chg", &chg);
		mixing_tree->Branch("sube", &sube);
		mixing_tree->Branch("pdg",&pdg);
	}


	std::ifstream instr(in_file_name.c_str(), std::ifstream::in);
	if(!instr.is_open()) cout << "filelist not found!! Exiting..." << endl;
	std::string filename;
	int ifile=0;

	const int MAXJETS = 500;

	Float_t t_calo_jtpt[MAXJETS], t_calo_jteta[MAXJETS], t_calo_jtphi[MAXJETS], t_calo_discr_ssvHighEff[MAXJETS], t_calo_discr_ssvHighPur[MAXJETS], t_calo_discr_csvV1[MAXJETS], t_calo_discr_csvV2[MAXJETS], t_calo_discr_prob[MAXJETS], t_calo_svtxm[MAXJETS], t_calo_svtxpt[MAXJETS], t_calo_svtxmcorr[MAXJETS], t_calo_svtxdl[MAXJETS], t_calo_svtxdls[MAXJETS], t_calo_rawpt[MAXJETS], t_calo_trackmax[MAXJETS];
	Float_t t_calo_refpt[MAXJETS], t_calo_refeta[MAXJETS], t_calo_refphi[MAXJETS];

	Int_t t_trkChg[MAXPARTICLES], t_calo_parton_flavor[MAXJETS];

    Float_t t_trkPt[MAXPARTICLES], t_trkEta[MAXPARTICLES], t_trkPhi[MAXPARTICLES], t_trkDxy1[MAXPARTICLES], t_trkDxyError1[MAXPARTICLES], t_trkDz1[MAXPARTICLES], t_trkDzError1[MAXPARTICLES], t_trkPtError[MAXPARTICLES], t_pfEcal[MAXPARTICLES], t_pfHcal[MAXPARTICLES], t_trkChi2[MAXPARTICLES];
	UChar_t t_trkNHit[MAXPARTICLES], t_trkNlayer[MAXPARTICLES], t_trkNdof[MAXPARTICLES]; //Run2
	Bool_t t_highPurity[MAXPARTICLES];

	vector<float> *csCandEta=0, *csCandPhi=0, *csCandPt=0;
	vector<int> *csCandId=0;

	vector<float> *pfEta=0, *pfPhi=0, *pfPt=0;
    vector<int> *pfId=0;

	vector<float> *t_pt=0, *t_phi=0, *t_eta=0, *t_chg=0, *t_sube=0; //Run2
	vector<int> *t_sta=0, *t_pdg=0; //Run2

	Int_t nTrk, mult, calo_nref, nPFpart, nCSpart;
	Int_t ncs_id145_pt2;  
	Int_t npf_id1_pt2;

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
		        //if(!my_file){ cout << "File cannot be found!!" << endl; exit(1); }	
		//}

        //TFile *my_file = TFile::Open("HiForestAOD_herwig.root");	

		if(!my_file){ cout << "File cannot be found!!" << endl; exit(1); }	

		if(my_file->IsZombie()) { 
			cout << "Is zombie" << std::endl;
 		    continue;
        }   

		if(do_PbPb){
			inp_tree = (TTree*)  my_file->Get(Form("akPu%dCaloJetAnalyzer/t",radius));
			inp_tree_CS = (TTree*) my_file->Get("pfcandAnalyzerCS/pfTree");
		}else{
			inp_tree = (TTree*)  my_file->Get(Form("ak%dCaloJetAnalyzer/t",radius));
			inp_tree2 = (TTree*)  my_file->Get("pfcandAnalyzer/pfTree");
			if(!inp_tree2){ cout << "PFCand Tree not found!! Exiting..." << endl; exit(1); }
			else inp_tree->AddFriend(inp_tree2);
		}

		inp_tree3 = (TTree*) my_file->Get("hiEvtAnalyzer/HiTree");
		if(!inp_tree3){ cout << "Evt Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree3);

		inp_tree4 = (TTree*) my_file->Get("skimanalysis/HltTree");
		if(!inp_tree4){ cout << "Skim Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree4);

		if(do_PbPb){
			inp_tree5 = (TTree*) my_file->Get("anaTrack/trackTree");
		}else{
			inp_tree5 = (TTree*) my_file->Get("ppTrack/trackTree");
		}
		if(!inp_tree5){ cout << "track Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree5);

		inp_tree6 = (TTree*) my_file->Get("hltanalysis/HltTree");
		if(is_data && !inp_tree6){ cout << "HLT Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree6);

		if(!is_data){ 
			inp_tree7 = (TTree*) my_file->Get("HiGenParticleAna/hi");
			if(!inp_tree7){ cout << "GenPart Tree not found!! Exiting..." << endl; exit(1); }
			else inp_tree->AddFriend(inp_tree7);
		}

		cout << "trees loaded" << endl;

		inp_tree->SetBranchAddress("trkPt",t_trkPt);
		inp_tree->SetBranchAddress("trkEta",t_trkEta);
		inp_tree->SetBranchAddress("trkPhi",t_trkPhi);
		inp_tree->SetBranchAddress("trkCharge",t_trkChg);

		inp_tree->SetBranchAddress("highPurity",t_highPurity);
		inp_tree->SetBranchAddress("vz",&vz);
		if(do_PbPb){
			inp_tree->SetBranchAddress("hiBin",&hiBin);
            inp_tree->SetBranchAddress("hiHF",&hiHF);		
		}
		inp_tree->SetBranchAddress("nref",&calo_nref);
		inp_tree->SetBranchAddress("jtpt",t_calo_jtpt);
		inp_tree->SetBranchAddress("jteta",t_calo_jteta);
		inp_tree->SetBranchAddress("jtphi",t_calo_jtphi);
		inp_tree->SetBranchAddress("rawpt",t_calo_rawpt);
		inp_tree->SetBranchAddress("trackMax",t_calo_trackmax);
		if(!is_data) {
			inp_tree->SetBranchAddress("refpt",t_calo_refpt);
			inp_tree->SetBranchAddress("refeta",t_calo_refeta);
			inp_tree->SetBranchAddress("refphi",t_calo_refphi);
		}

		inp_tree->SetBranchAddress("nTrk",&nTrk);
		inp_tree->SetBranchAddress("trkDxy1",t_trkDxy1);
		inp_tree->SetBranchAddress("trkDxyError1",t_trkDxyError1);
		inp_tree->SetBranchAddress("trkDz1",t_trkDz1);
		inp_tree->SetBranchAddress("trkDzError1",t_trkDzError1);
		inp_tree->SetBranchAddress("trkPtError",t_trkPtError);
		inp_tree->SetBranchAddress("trkNHit",t_trkNHit);
		inp_tree->SetBranchAddress("trkNlayer",t_trkNlayer);
		inp_tree->SetBranchAddress("trkChi2",t_trkChi2);
		inp_tree->SetBranchAddress("trkNdof",t_trkNdof);
		inp_tree->SetBranchAddress("pfEcal",t_pfEcal);
		inp_tree->SetBranchAddress("pfHcal",t_pfHcal);

		if(do_PbPb){
			inp_tree_CS->SetBranchAddress("nPFpart",&nCSpart);
			inp_tree_CS->SetBranchAddress("pfPt",&csCandPt);
			inp_tree_CS->SetBranchAddress("pfEta",&csCandEta);
			inp_tree_CS->SetBranchAddress("pfPhi",&csCandPhi);
			inp_tree_CS->SetBranchAddress("pfId",&csCandId);
		}
        else{
			inp_tree->SetBranchAddress("nPFpart",&nPFpart);
			inp_tree->SetBranchAddress("pfPt",&pfPt);
	        inp_tree->SetBranchAddress("pfEta",&pfEta);
	        inp_tree->SetBranchAddress("pfPhi",&pfPhi);
	        inp_tree->SetBranchAddress("pfId",&pfId);        	
        } 

		if(!is_data){
			inp_tree->SetBranchAddress("pthat",&pthat);
			inp_tree->SetBranchAddress("refparton_flavor",t_calo_parton_flavor);
		}

		inp_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHENoiseFilterResultRun2Loose);

		if(!do_PbPb){
			inp_tree->SetBranchAddress("pPAprimaryVertexFilter",&pvFilter); //Run2
			inp_tree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter); //Run2
		}
		else{	
			inp_tree->SetBranchAddress("pcollisionEventSelection",&eventSelection);			
			inp_tree->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
			inp_tree->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);
			inp_tree->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3);
		}

		if(is_data){
			if(!do_PbPb){ 
				inp_tree->SetBranchAddress("HLT_AK4PFJet80_Eta5p1_v1", &HLT_Jet80);
				inp_tree->SetBranchAddress("HLT_AK4PFJet100_Eta5p1_v1", &HLT_Jet100);
			}
			else{
				inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1", &HLT_Jet80);
				inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1", &HLT_Jet100);
				//inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1_Prescl", &HLT_Jet80_ps);
				//inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1_Prescl", &HLT_Jet100_ps);
			}
		}

		if(!is_data){
			inp_tree->SetBranchAddress("mult",&mult);
			inp_tree->SetBranchAddress("pt", &t_pt);
			inp_tree->SetBranchAddress("phi", &t_phi);
			inp_tree->SetBranchAddress("eta", &t_eta);
			inp_tree->SetBranchAddress("chg", &t_chg);
			inp_tree->SetBranchAddress("sube", &t_sube);
			inp_tree->SetBranchAddress("pdg", &t_pdg);
			/*if(!do_PbPb)*/ inp_tree->SetBranchAddress("sta", &t_sta);
		}


		int n_evt = inp_tree->GetEntriesFast();

		//n_evt=100;
		cout << "Entries: "<< n_evt << endl;
		for(int evi = 0; evi < n_evt; evi++) {

			inp_tree->GetEntry(evi);
			if(do_PbPb) inp_tree_CS->GetEntry(evi);

            if(fabs(vz) >= 15.) continue;
			if(do_PbPb && hiHF > 5500.) continue;

			if(is_data && !HLT_Jet80 && !HLT_Jet100) continue;

            if(!do_PbPb && (!pvFilter || !pBeamScrapingFilter || !HBHENoiseFilterResultRun2Loose)) continue;

            if(do_PbPb && (!pvFilter || !eventSelection || !phfCoincFilter3 || !pclusterCompatibilityFilter || !HBHENoiseFilterResultRun2Loose)) continue;

            if(!is_data && pthat<80. ) continue;
            
            bool found_inclusive_jet = 0; 

 			//start calo jet loop
			for(int j4i = 0; j4i < calo_nref ; j4i++) {

				if(fabs(t_calo_jteta[j4i]) > 1.6 || t_calo_jtpt[j4i] < 80.) continue;
				else if(t_calo_trackmax[j4i]/t_calo_rawpt[j4i] > 0.98 ||t_calo_trackmax[j4i]/t_calo_rawpt[j4i] < 0.01) continue;
                if (!is_data && t_calo_refpt[j4i] < 0.) continue;

                found_inclusive_jet = 1;

                //npf cand for pp and Pythia
				if(!do_PbPb){
    				npf_id1_pt2 = 0;

	    			for(int icand=0; icand<nPFpart; icand++){
		    			if(abs(pfEta->at(icand))>2.4) continue;
			    		double dr = sqrt(pow(pfEta->at(icand) - t_calo_jteta[j4i], 2) + pow(acos(cos(t_calo_jtphi[j4i] - pfPhi->at(icand))),2));
				    	if(dr > 0.4 || pfPt->at(icand)<2.) continue;
				    	if(pfId->at(icand)==1){
					    	npf_id1_pt2++;
				        }
			    	}
			    	corrected_pt = corrpt->getCorrection(!do_PbPb, npf_id1_pt2, 0, t_calo_jtpt[j4i], t_calo_jteta[j4i]);
			    }
			    else{
    				ncs_id145_pt2 = 0;

	    			for(int icand=0; icand<nCSpart; icand++){
		    			if(abs(csCandEta->at(icand))>2.4) continue;
			    		double dr = sqrt(pow(csCandEta->at(icand) - t_calo_jteta[j4i], 2) + pow(acos(cos(t_calo_jtphi[j4i] - csCandPhi->at(icand))),2));
				    	if(dr > 0.4 || csCandPt->at(icand)<2.) continue;
				    	if(csCandId->at(icand)==1 || csCandId->at(icand)==4 || csCandId->at(icand)==5){
					    	ncs_id145_pt2++;
				        }
			    	}
                    if(ncs_id145_pt2 > 20 && t_calo_jtpt[j4i] <= 120.) corrected_pt = corrpt->getCorrection(!do_PbPb, 20, hiBin, t_calo_jtpt[j4i], t_calo_jteta[j4i]);
                    else corrected_pt = corrpt->getCorrection(!do_PbPb, ncs_id145_pt2, hiBin, t_calo_jtpt[j4i], t_calo_jteta[j4i]);
			    }	

                if(corrected_pt <= 0.) continue;

				calo_jteta.push_back(t_calo_jteta[j4i]);
				calo_jtphi.push_back(t_calo_jtphi[j4i]);
				calo_jtpt.push_back(t_calo_jtpt[j4i]);
				calo_corrpt.push_back(corrected_pt);
				if(!is_data){ 
					calo_refpt.push_back(t_calo_refpt[j4i]);
					calo_refeta.push_back(t_calo_refeta[j4i]);
					calo_refphi.push_back(t_calo_refphi[j4i]);
					calo_parton_flavor.push_back(t_calo_parton_flavor[j4i]);
				}
			} /// calo jet loop

            if(found_inclusive_jet==0) continue; 

			//// reco track loop
			for(int itrk=0;itrk<nTrk;itrk++){
                if(!passTrackCuts(!do_PbPb, 1, t_trkPt[itrk], t_trkEta[itrk], t_highPurity[itrk], t_trkChi2[itrk], t_trkNdof[itrk], t_trkNlayer[itrk], t_trkNHit[itrk], t_pfHcal[itrk], t_pfEcal[itrk], t_trkDxy1[itrk]/t_trkDxyError1[itrk], t_trkDz1[itrk]/t_trkDzError1[itrk])) continue;

                bool in_jetcone = 0;

				for(int j4i = 0; j4i < calo_nref ; j4i++) {

					if(fabs(t_calo_jteta[j4i]) > 1.6 || t_calo_jtpt[j4i] < 80.) continue;
					else if(t_calo_trackmax[j4i]/t_calo_rawpt[j4i] > 0.98 ||t_calo_trackmax[j4i]/t_calo_rawpt[j4i] < 0.01) continue;
	                if (is_data && t_calo_refpt[j4i] < 0.) continue;

	                double dr = sqrt(pow(t_trkEta[itrk] - t_calo_jteta[j4i], 2) + pow(acos(cos(t_calo_jtphi[j4i] - t_trkPhi[itrk])),2));
	                double refl_dr = sqrt(pow(t_trkEta[itrk] + t_calo_jteta[j4i], 2) + pow(acos(cos(t_calo_jtphi[j4i] - t_trkPhi[itrk])),2));
                
                    if(dr <= 0.4 || refl_dr <= 0.4){
                        in_jetcone = 1;
                        break;
                    } 
                } 

                if(in_jetcone == 0) continue;

				trkChg.push_back(t_trkChg[itrk]);
				trkEta.push_back(t_trkEta[itrk]);
				trkPhi.push_back(t_trkPhi[itrk]);
				trkPt.push_back(t_trkPt[itrk]);
			}

			if(!is_data){
				//gen particles loop
				for(int ipart=0;ipart<mult;ipart++){
				  
                    if(/*!do_PbPb && */t_sta->at(ipart)!=1) continue; //remove non final state particles
                    if(!passGenTrackCuts(t_pt->at(ipart), t_eta->at(ipart), t_chg->at(ipart))) continue;

                    bool in_jetcone = 0;

					for(int j4i = 0; j4i < calo_nref ; j4i++) {

						if(fabs(t_calo_jteta[j4i]) > 1.6 || t_calo_jtpt[j4i] < 80.) continue;
						else if(t_calo_trackmax[j4i]/t_calo_rawpt[j4i] > 0.98 ||t_calo_trackmax[j4i]/t_calo_rawpt[j4i] < 0.01) continue;
		                if (t_calo_refpt[j4i] < 0.) continue;

		                double dr = sqrt(pow(t_eta->at(ipart) - t_calo_jteta[j4i], 2) + pow(acos(cos(t_calo_jtphi[j4i] - t_phi->at(ipart))),2));
		                double refl_dr = sqrt(pow(t_eta->at(ipart) + t_calo_jteta[j4i], 2) + pow(acos(cos(t_calo_jtphi[j4i] - t_phi->at(ipart))),2));
		                double gen_dr = sqrt(pow(t_eta->at(ipart) - t_calo_refeta[j4i], 2) + pow(acos(cos(t_calo_refphi[j4i] - t_phi->at(ipart))),2));
	                
	                    if(dr <= 0.4 || refl_dr <= 0.4 || gen_dr <= 0.4){
	                        in_jetcone = 1;
	                        break;
	                    } 
	                } 

	                if(in_jetcone == 0) continue;

					eta.push_back(t_eta->at(ipart));
					phi.push_back(t_phi->at(ipart));
					pt.push_back(t_pt->at(ipart));
					chg.push_back(t_chg->at(ipart));
					sube.push_back(t_sube->at(ipart));
					pdg.push_back(t_pdg->at(ipart));				
				}
			}

			///// Fill it
			mixing_tree->Fill();

			///// Reset
			trkEta.clear();
			trkPhi.clear();
			trkPt.clear();
			trkChg.clear();
			calo_jteta.clear();
			calo_jtphi.clear();
			calo_jtpt.clear();
			calo_corrpt.clear();

			if(!is_data){
     			calo_parton_flavor.clear();
				calo_refpt.clear();
				calo_refeta.clear();
				calo_refphi.clear();

				pt.clear();
				phi.clear();
				eta.clear();
				chg.clear();
				sube.clear();
				pdg.clear();
			}

		}  ///event loop
		my_file->Close();

    }//file loop

	cout<<"writing"<<endl;

	output_file->cd();
	mixing_tree->Write();
	output_file->Close();

	cout<<"done"<<endl;

}
