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

    for(int i=0;i<tree->GetEntries();i++)
    {
        tree->GetEntry(i);
        hChi2overNdf->Fill(chi2overndf);
        if(chi2overndf>=20)
        {
            string str_histo = Label(layerLabel)+" NStrip="+to_string(nstrip)+" NStripSat254="+to_string(nstripsat254)+" NStripSat255="+to_string(nstripsat255);
            HistoOfInterest = (TH2F*) ifile->Get(str_histo.c_str());
            HistoOfInterest->Write();
            string str_profile = str_histo+"_pfx";
            ProfileOfInterest = (TProfile*) ifile->Get(str_profile.c_str());
            ProfileOfInterest->Write();
        }
    }
    ofile->Write();
    ofile->Close();
    delete ofile;


    return 0;
}