#include <iostream>
using namespace std;
#include <memory>
#include <string>
#include <vector>

#include "TTree.h"
#include "TLatex.h"
#include "TMath.h"
#include "TBranch.h"
#include "TFile.h"
#include "TStyle.h"
#include "TPad.h"
#include "TProfile.h"
#include "TPave.h"
#include "TPaveText.h"
#include "TString.h"
#include "TPaveStats.h"

void Best(bool isEB, bool isPi0=true){

    //FILE *f_1;
    FILE *f_2;

    //f_1=fopen("Value_Sorted.txt","w");
    TString input;
    if(isEB){
	  input="Value_notSorted_EB.txt";
	  if(!isPi0) input="Value_notSorted_EB_eta.txt";
    }
    else{
	  input="Value_notSorted_EE.txt";
	  if(!isPi0) input="Value_notSorted_EE_eta.txt";
    }
    f_2=fopen(input.Data(),"r");

    vector< float > values_eff; values_eff.clear();
    vector< float > values_sb; values_sb.clear();
    vector< float > values_mu; values_mu.clear();
    vector< float > values_chi; values_chi.clear();
    vector< int > values_bin; values_bin.clear();
    while(!feof(f_2)){
	  char sb_s, chi_s, eff_s, smu_s, bin_s;
	  float sb, chi, eff, smu;
	  int bin;
	  fscanf(f_2,"%s %i  %s %f  %s %f  %s %f  %s %f", &bin_s, &bin, &sb_s, &sb, &smu_s, &smu, &chi_s, &chi, &eff_s, &eff);
	  //	  cout<<bin_s<<" "<<bin<<"  "<<sb_s<<"  "<<sb<<"  "<<smu_s<<"  "<<smu<<"  "<<chi_s<<"  "<<chi<<"  "<<eff_s<<"  "<<eff<<endl;

	  values_eff.push_back( (float)eff );
	  values_sb.push_back(  (float)sb );
	  values_mu.push_back(  (float)smu );
	  values_chi.push_back( (float)chi );
	  values_bin.push_back( (int)bin );
    }


    float tmp1 = -1.;   int bin1 = -1.;
    float tmp2 = 100.;  int bin2 = -1.;
    for(unsigned int i(0); i<values_eff.size(); i++){

	  if(values_sb[i] > tmp1){
		tmp1 = values_sb[i];
		bin1 = i; cout<<"bin sb  "<<bin1<<" bin "<<values_bin[bin1]<<endl;
	  }
	  if(values_mu[i] < tmp2){
		tmp2 = values_mu[i];
		bin2 = i; cout<<"bin mu  "<<bin2<<" bin "<<values_bin[bin2]<<endl;
	  }
    }

    cout<<"S/B) Bin: "<<values_bin[bin1]<<" sb: "<<values_sb[bin1]<<" smu: "<<values_mu[bin1]<<" eff: "<<values_eff[bin1]<<" chi: "<<values_chi[bin1]<<endl;
    cout<<"M/U) Bin: "<<values_bin[bin2]<<" sb: "<<values_sb[bin2]<<" smu: "<<values_mu[bin2]<<" eff: "<<values_eff[bin2]<<" chi: "<<values_chi[bin2]<<endl;

}
