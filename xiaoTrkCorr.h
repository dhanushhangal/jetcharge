
#ifndef ROOT_TH2F
#include  "TH2F.h"
#endif

#ifndef ROOT_TFile
#include "TFile.h"
#endif

class xiaoTrkCorr{
	public: 
		xiaoTrkCorr(TString f);
		int binarySearch(float key, float* arr, int i_max, int i_min );
		float getTrkCorr(float pt, float eta, float phi, int hibin);
	public:
		TFile* file;
		TH2F* corrTable[22][17];
		int nptbin = 22;
		int ncentbin = 17;
		float ptbin[23] = {0.7,0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2,
			2.5, 3, 3.5, 4,5, 6, 7, 8, 10, 12, 16, 20, 50, 999};
		float centbin[18];
};

xiaoTrkCorr::xiaoTrkCorr(TString f){
	file = TFile::Open(f);
	for (int i=0; i<ncentbin; ++i)centbin[i]= i*10;
	centbin[ncentbin] = 200;
	for(int i=0; i<nptbin; ++i){
		for(int j=0;j<ncentbin; ++j){
			corrTable[i][j]=(TH2F*)file->Get(Form("corr_%d_%d",i,j));
		}
	}
}

int xiaoTrkCorr::binarySearch(float key, float* arr, int i_max, int i_min ){
		if(key> arr[i_max] ) return -1;
		if(key< arr[i_min] ) return -1;
		int mid = floor(float(i_max +i_min)/2);
		//	cout<<mid<<endl;
		if(mid == i_min ) return mid;
		if( arr[mid]> key) return binarySearch(key, arr, mid, i_min);
		else if( arr[mid] < key) return binarySearch(key, arr, i_max, mid);
		else return mid;
}

float xiaoTrkCorr::getTrkCorr(float pt, float eta, float phi, int cent){
	int jpt = binarySearch(pt, ptbin, nptbin,0);
	int jcent = binarySearch(cent, centbin, ncentbin,0);
	if(jpt <0 || jcent <0) {
//		std::cout<<"error!"<<std::endl;
		std::cout<<"jpt="<<jpt<<", jcent="<<jcent<<std::endl;
		std::cout<<"pt="<<pt<<", cent="<<cent<<std::endl;
		return 0;
	}
	return corrTable[jpt][jcent]->GetBinContent(corrTable[jpt][jcent]->FindBin(eta,phi));
}

