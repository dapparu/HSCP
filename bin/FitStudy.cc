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

using namespace std;

string Label(int i)
{
    if(i==1)   return "TIB1";
    if(i==2)   return "TIB2";
    if(i==3)   return "TIB3";
    if(i==4)   return "TIB4";
    if(i==5)   return "TOB1";
    if(i==6)   return "TOB2";
    if(i==7)   return "TOB3";
    if(i==8)   return "TOB4";
    if(i==9)   return "TOB5";
    if(i==10)  return "TOB6";
    if(i==11)  return "TID1";
    if(i==12)  return "TID2";
    if(i==13)  return "TEC1";
    if(i==14)  return "TEC2";
    if(i==15)  return "TEC3";
    if(i==16)  return "TEC4";
    if(i==17)  return "TEC5";
    if(i==18)  return "TEC6";
    if(i==19)  return "TEC7";
    if(i==20)  return "TEC8";
    if(i==21)  return "TEC9";
}

int main(int argc,char** argv)
{
    string s1 = argv[1];
    string s2 = s1.substr(0,s1.find('.'))+"_Study.root";
    TFile* ifile = TFile::Open(s1.c_str());
    TTree* tree = (TTree*) ifile->Get("tree");
    TFile* ofile = new TFile(s2.c_str(),"RECREATE");


    float p0    = 0.;
    float p1    = 0.;
    float p0err = 0.;
    float p1err = 0.;
    float chi2  = 0.;
    int ndf     = 0;
    float chi2overndf   = 0.;
    int layerLabel      = 0;
    int nstrip          = 0;
    int nstripsat254    = 0;
    int nstripsat255    = 0;

    tree->SetBranchAddress("p0",&p0);
    tree->SetBranchAddress("p1",&p1);
    tree->SetBranchAddress("p0err",&p0err);
    tree->SetBranchAddress("p1err",&p1err);
    tree->SetBranchAddress("chi2",&chi2);
    tree->SetBranchAddress("ndf",&ndf);
    tree->SetBranchAddress("chi2overndf",&chi2overndf);
    tree->SetBranchAddress("layerLabel",&layerLabel);
    tree->SetBranchAddress("nstrip",&nstrip);
    tree->SetBranchAddress("nstripsat254",&nstripsat254);
    tree->SetBranchAddress("nstripsat255",&nstripsat255);

    TH1F* hChi2overNdf = new TH1F("Chi2/Ndf","Chi2/Ndf",200,0,10);
    TH2F* HistoOfInterest;
    TProfile* ProfileOfInterest;
    vector<float> vect_p0;
    vector<float> vect_p1;
    vector<float> vect_p0err;
    vector<float> vect_p1err;
    vector<float> vect_nstripsat;

    /*float tablayer[10]={1,2,3,4,5,6,7,8,9,10};
    float tabp0[10];//={0,0,0,0,0,0,0,0,0,0};
    float tabp0err[10]={0,0,0,0,0,0,0,0,0,0};*/
    float tabx[740];
    float tabp0[740];
    float tabp0err[740];
    float tabp1[740];
    float tabp1err[740];
    for(int i=0;i<740;i++) tabx[i]=i;


int count=0;

    for(int i=0;i<tree->GetEntries();i++)
    {

        tree->GetEntry(i);
        if(layerLabel==5 && nstripsat254>=1 && nstripsat255==0)
        {
            hChi2overNdf->Fill(chi2overndf);
            string str_histo = Label(layerLabel)+" NStrip="+to_string(nstrip)+" NStripSat254="+to_string(nstripsat254)+" NStripSat255="+to_string(nstripsat255);
            HistoOfInterest = (TH2F*) ifile->Get(str_histo.c_str());
            //HistoOfInterest->Write();
            string str_profile = str_histo+"_pfx";
            ProfileOfInterest = (TProfile*) ifile->Get(str_profile.c_str());
            //ProfileOfInterest->Write();
            TCanvas* CanvasOfInterest = new TCanvas(str_histo.c_str(),str_histo.c_str());
            HistoOfInterest->Draw("colz");
            ProfileOfInterest->Draw("same");
            CanvasOfInterest->Write();
        }
        tabp0[i]=p0;
        tabp0err[i]=p0err;
        tabp1[i]=p1;
        tabp1err[i]=p1err;
    }
    //TGraphErrors* gr = new TGraphErrors(vect_nstripsat.size(),&vect_nstripsat[0],&vect_p0[0],0,&vect_p0err[0]);
    TCanvas* cgr = new TCanvas("graph","graph");
    TGraphErrors* grTIB1 = new TGraphErrors(74,tabx,tabp0,0,tabp0err); 
    grTIB1->SetTitle("TIB1 | p0");
    grTIB1->Draw("AP");
    cgr->Write();
    TGraphErrors* gr1 = new TGraphErrors(740,tabx,tabp0,0,tabp0err);
    gr1->SetTitle("p0");
    gr1->GetYaxis()->SetTitle("p0");
    gr1->Draw("AP");
    cgr->Write();
    TGraphErrors* gr2 = new TGraphErrors(740,tabx,tabp1,0,tabp1err);
    gr2->SetTitle("p1");
    gr2->GetYaxis()->SetTitle("p1");
    gr2->Draw("AP");
    cgr->Write();
    ofile->Write();
    ofile->Close();
    delete ofile;


    return 0;
}