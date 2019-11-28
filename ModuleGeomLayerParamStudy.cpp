#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THStack.h"
#include "TProfile.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TText.h"
#include "TLegend.h"
#include "TMath.h"
#include "Rtypes.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TRandom.h"


using namespace std;


//FUNCTIONS


float Mean(vector<float> Vect)
{
	float mean = .0;
	float size = Vect.size();
	for(int i=0;i<size;i++)
	{
		mean+=Vect[i];
	}
	return mean/size;
}

float StdDev(vector<float> Vect)
{
	float size = Vect.size();
	float res=.0;
	for(int i=0; i<size; i++)
	{
		res+=pow(Vect[i]-Mean(Vect),2);
	}
	res=pow(res/(float)(size-1),0.5);
	return res;
}

float Chi2(vector<float> Vect,vector<float> VectErr)
{
    float res=0.;
    float mean=Mean(Vect);
    for(int i=0;i<Vect.size();i++)
    {
        cout<<Vect[i]<<" "<<VectErr[i]<<endl;
        res+=pow((Vect[i]-mean)/VectErr[i],2);
    }
    return res;
}





//MAIN 


int main(int argc,char** argv)
{

    TFile* fileLayer = TFile::Open("testMethodLayer.root");
	TTree* treeLayer = (TTree*) fileLayer->Get("tree");
    TFile* fileModulGeom = TFile::Open("testMethodModulGeom.root");
	TTree* treeModulGeom = (TTree*) fileModulGeom->Get("tree");
    TFile* ofile = new TFile("ModulGeomLayerParam_Study.root","RECREATE");

    float p0    = 0.;
    float p1    = 0.;
    float p0err = 0.;
    float p1err = 0.;
    float chi2  = 0.;
    int ndf     = 0;
    float chi2overndf   = 0.;
    int Label           = 0;
    int nstrip          = 0;
    int nstripsat254    = 0;
    int nstripsat255    = 0;




//TREE LAYER

    treeLayer->SetBranchAddress("p0",&p0);
    treeLayer->SetBranchAddress("p1",&p1);
    treeLayer->SetBranchAddress("p0err",&p0err);
    treeLayer->SetBranchAddress("p1err",&p1err);
    treeLayer->SetBranchAddress("chi2",&chi2);
    treeLayer->SetBranchAddress("ndf",&ndf);
    treeLayer->SetBranchAddress("chi2overndf",&chi2overndf);
    treeLayer->SetBranchAddress("Label",&Label);
    treeLayer->SetBranchAddress("nstrip",&nstrip);
    treeLayer->SetBranchAddress("nstripsat254",&nstripsat254);
    treeLayer->SetBranchAddress("nstripsat255",&nstripsat255);

    vector<float> p1__TOB1_4;
    vector<float> p0__TOB1_4;
    vector<float> p1err__TOB1_4;
    vector<float> p0err__TOB1_4;

    for(int counter_treeLayer=0;counter_treeLayer<treeLayer->GetEntries();counter_treeLayer++)
    {
        treeLayer->GetEntry(counter_treeLayer);
        if(Label>=5 && Label<=8)
        {
            if(nstrip==5 && nstripsat254==1 && nstripsat255==0)
            {
                p1__TOB1_4.push_back(p1);
                p0__TOB1_4.push_back(p0);
                p1err__TOB1_4.push_back(p1err);
                p0err__TOB1_4.push_back(p0err);
            }
        }
    }

    float Meanp1 = Mean(p1__TOB1_4);
    float Sigmap1 = StdDev(p1__TOB1_4);
    float Meanp0 = Mean(p0__TOB1_4);
    float Sigmap0 = StdDev(p0__TOB1_4);

    float Chi2_p1 = Chi2(p1__TOB1_4,p1err__TOB1_4);
    cout<<"Chi2_p1 "<<Chi2_p1<<endl;

    float ProbaChi2_p1 = TMath::Prob(Chi2_p1,4);
    cout<<"Prob_p1 "<<ProbaChi2_p1<<endl;

    float Chi2_p0 = Chi2(p0__TOB1_4,p0err__TOB1_4);
    cout<<"Chi2_p0 "<<Chi2_p0<<endl;

    float ProbaChi2_p0 = TMath::Prob(Chi2_p0,4);
    cout<<"Prob_p0 "<<ProbaChi2_p0<<endl;



//TREE MODULGEOM


    treeModulGeom->SetBranchAddress("p0",&p0);
    treeModulGeom->SetBranchAddress("p1",&p1);
    treeModulGeom->SetBranchAddress("p0err",&p0err);
    treeModulGeom->SetBranchAddress("p1err",&p1err);
    treeModulGeom->SetBranchAddress("chi2",&chi2);
    treeModulGeom->SetBranchAddress("ndf",&ndf);
    treeModulGeom->SetBranchAddress("chi2overndf",&chi2overndf);
    treeModulGeom->SetBranchAddress("Label",&Label);
    treeModulGeom->SetBranchAddress("nstrip",&nstrip);
    treeModulGeom->SetBranchAddress("nstripsat254",&nstripsat254);
    treeModulGeom->SetBranchAddress("nstripsat255",&nstripsat255);



    for(int counter_treeModulGeom=0;counter_treeModulGeom<treeModulGeom->GetEntries();counter_treeModulGeom++)
    {
        treeModulGeom->GetEntry(counter_treeModulGeom);
        if(Label==4)
        {
            if(nstrip==5 && nstripsat254==1 && nstripsat255==0)
            {
                cout<<"ModulGeom : "<<p1<<"+-"<<p1err<<"    layer : "<<Meanp1<<"+-"<<Sigmap1<<endl;
                cout<<"ModulGeom : "<<p0<<"+-"<<p0err<<"    layer : "<<Meanp0<<"+-"<<Sigmap0<<endl;
            }
        }
    }


    return 0;

}
