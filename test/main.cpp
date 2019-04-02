#include <iostream>
#include <vector>
#include <stdlib.h>
#include <iterator>
#include <map>
#include "TFile.h"
#include "TTree.h"
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
#include "TExec.h"

#include "../interface/Builder.h"
#include "../interface/Estimator.h"
#include "../interface/PlotterHisto.h"
#include "../interface/Labellizer.h" 

using namespace std;

int main(int argc,char** argv){

    TChain chain;
	chain.SetName("stage/ttree");
	for(int i=1;i<argc;i++) chain.Add(argv[i]);

    string s1 = argv[1];
    string s2 = s1.substr(0,s1.find('.'))+"_results.root";

    Builder* b1 = new Builder(chain);
    b1->SetBranchAdd();
    int nentries = b1->GetEntries();

    TH1F* hDiff_rel_ElossQ_tot = new TH1F("hDiff_rel_ElossQ_tot","hDiff_rel_ElossQ_tot",200,-1,5);
    TH1F* hDiff_rel_ElossQ_NoSat = new TH1F("hDiff_rel_ElossQ_NoSat","hDiff_rel_ElossQ_NoSat",200,-1,5);
    TH1F* hDiff_rel_ElossQ_Sat = new TH1F("hDiff_rel_ElossQ_Sat","hDiff_rel_ElossQ_Sat",200,-1,5);
    TH1F* hDiff_rel_ElossQ_Sat254 = new TH1F("hDiff_rel_ElossQ_Sat254","hDiff_rel_ElossQ_Sat254",200,-1,5);
    TH1F* hDiff_rel_ElossQ_Sat255 = new TH1F("hDiff_rel_ElossQ_Sat255","hDiff_rel_ElossQ_Sat255",200,-1,5);

    TH1F* hRatio_ElossQ_tot = new TH1F("hRatio_ElossQ_tot","hRatio_ElossQ_tot",200,0,8);
    TH1F* hRatio_ElossQ_NoSat = new TH1F("hRatio_ElossQ_NoSat","hRatio_ElossQ_NoSat",200,0,8);
    TH1F* hRatio_ElossQ_Sat = new TH1F("hRatio_ElossQ_Sat","hRatio_ElossQ_Sat",200,0,8);
    TH1F* hRatio_ElossQ_Sat254 = new TH1F("hRatio_ElossQ_Sat254","hRatio_ElossQ_Sat254",200,0,8);
    TH1F* hRatio_ElossQ_Sat255 = new TH1F("hRatio_ElossQ_Sat255","hRatio_ElossQ_Sat255",200,0,8);

    TH2F* h2ElossvQ_tot = new TH2F("h2ElossvQ_tot","h2ElossvQ_tot",300,0,900*pow(10,-6),300,0,900*pow(10,-6));
    TH2F* h2ElossvQ_NoSat = new TH2F("h2ElossvQ_NoSat","h2ElossvQ_NoSat",300,0,900*pow(10,-6),300,0,900*pow(10,-6));
    TH2F* h2ElossvQ_Sat = new TH2F("h2ElossvQ_Sat","h2ElossvQ_Sat",300,0,900*pow(10,-6),300,0,900*pow(10,-6));
    TH2F* h2ElossvQ_Sat254 = new TH2F("h2ElossvQ_Sat254","h2ElossvQ_Sat254",300,0,900*pow(10,-6),300,0,900*pow(10,-6));
    TH2F* h2ElossvQ_Sat255 = new TH2F("h2ElossvQ_Sat255","h2ElossvQ_Sat255",300,0,3000*pow(10,-6),300,0,3000*pow(10,-6));

    TProfile* profElossVsQ_tot = new TProfile("profElossVsQ_tot","profElossVsQ_tot",50,0,900*pow(10,-6),"");
    TProfile* profElossVsQ_NoSat = new TProfile("profElossVsQ_NoSat","profElossVsQ_NoSat",50,0,900*pow(10,-6),"");
    TProfile* profElossVsQ_Sat = new TProfile("profElossVsQ_Sat","profElossVsQ_Sat",50,0,900*pow(10,-6),"");
    TProfile* profElossVsQ_Sat254 = new TProfile("profElossVsQ_Sat254","profElossVsQ_Sat254",50,0,900*pow(10,-6),"");
    TProfile* profElossVsQ_Sat255 = new TProfile("profElossVsQ_Sat255","profElossVsQ_Sat255",50,0,3000*pow(10,-6),"");

    TH1F* hLayer = new TH1F("hLayer","hLayer",80,0,80);
    TH1F* hLayerLabel = new TH1F("hLayerLabel","hLayerLabel",21,0,21);
    TH1F* hLayerLabelSat254 = new TH1F("hLayerLabelSat254","hLayerLabelSat254",21,0,21);
    TH1F* hLayerLabelSat255 = new TH1F("hLayerLabelSat255","hLayerLabelSat255",21,0,21);

    TH1F* hPt_tot = new TH1F("hPt_tot","hPt_tot",100,0,3000);
    TH1F* hPt_NoSat = new TH1F("hPt_NoSat","hPt_NoSat",100,0,3000);
    TH1F* hPt_Sat = new TH1F("hPt_Sat","hPt_Sat",100,0,3000);
    TH1F* hPt_Sat254 = new TH1F("hPt_Sat254","hPt_Sat254",100,0,3000);
    TH1F* hPt_Sat255 = new TH1F("hPt_Sat255","hPt_Sat255",100,0,3000);

    TH1F* hNClusterPerTrack = new TH1F("hNClusterPerTrack","hNClusterPerTrack",25,0,25);
    TH1F* hNClusterSat254PerTrack = new TH1F("hNClusterSat254PerTrack","hNClusterSat254PerTrack",20,0,20);
    TH1F* hNClusterSat255PerTrack = new TH1F("hNClusterSat255PerTrack","hNClusterSat255PerTrack",20,0,20);
    TH1F* hRatio_NClusterSat254 = new TH1F("hRatio_NClusterSat254","hRatio_NClusterSat254",20,0,1.05);
    TH1F* hRatio_NClusterSat255 = new TH1F("hRatio_NClusterSat255","hRatio_NClusterSat255",20,0,1.05);

    TH2F* h2RatioSatPt254 = new TH2F("h2RatioSatPt254","h2RatioSatPt254",100,0,3000,20,0,1.05);
    TProfile* profSatPt254 = new TProfile("profSatPt254","profSatPt254",100,0,3000,"");

    TH2F* h2RatioSatEloss254 = new TH2F("h2RatioSatEloss254","h2RatioSatEloss254",300,0,900*pow(10,-6),20,0,1.05);
    TProfile* profSatEloss254 = new TProfile("profSatEloss254","profSatEloss254",300,0,900*pow(10,-6),"");

    TH2F* h2RatioSatPartID254 = new TH2F("h2RatioSatPartID254","h2RatioSatPartID254",6,0,6,20,0,1.05);
    TProfile* profSatPartID254 = new TProfile("profSatPartID254","profSatPartID254",6,0,6,"");

    TH2F* h2RatioSatPt255 = new TH2F("h2RatioSatPt255","h2RatioSatPt255",100,0,3000,20,0,1.05);
    TProfile* profSatPt255 = new TProfile("profSatPt255","profSatPt255",100,0,3000,"");

    TH2F* h2RatioSatEloss255 = new TH2F("h2RatioSatEloss255","h2RatioSatEloss255",300,0,900*pow(10,-6),20,0,1.05);
    TProfile* profSatEloss255 = new TProfile("profSatEloss255","profSatEloss255",300,0,900*pow(10,-6),"");

    TH2F* h2RatioSatPartID255 = new TH2F("h2RatioSatPartID255","h2RatioSatPartID255",6,0,6,20,0,1.05);
    TProfile* profSatPartID255 = new TProfile("profSatPartID255","profSatPartID255",6,0,6,"");

    vector<TH2F*> VectNStrip_h2EvQ_NoSat;
    vector<TH2F*> VectNStrip_h2EvQ_Sat;
    vector<TH2F*> VectNStrip_h2EvQ_Sat254;
    vector<TH2F*> VectNStrip_h2EvQ_Sat255;
    vector<TH2F*> VectNStripSat_h2EvQ_254;
    vector<TH2F*> VectNStripSat_h2EvQ_255;

    for(int i=1;i<11;i++)
    {   
        string str_i = to_string(i);
        VectNStrip_h2EvQ_NoSat.push_back(new TH2F(("NoSatNstrips_"+str_i).c_str(),("NoSatNstrips_"+str_i).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
        VectNStrip_h2EvQ_Sat.push_back(new TH2F(("SatNstrips_"+str_i).c_str(),("SatNstrips_"+str_i).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
        VectNStrip_h2EvQ_Sat254.push_back(new TH2F(("Sat254Nstrips_"+str_i).c_str(),("Sat254Nstrips_"+str_i).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
        VectNStrip_h2EvQ_Sat255.push_back(new TH2F(("Sat255Nstrips_"+str_i).c_str(),("Sat255Nstrips_"+str_i).c_str(),300,0,3000*pow(10,-6),300,0,3000*pow(10,-6)));
        VectNStripSat_h2EvQ_254.push_back(new TH2F(("Sat254NstripsSat_"+str_i).c_str(),("Sat254NstripsSat_"+str_i).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
        VectNStripSat_h2EvQ_255.push_back(new TH2F(("Sat255NstripsSat_"+str_i).c_str(),("Sat255NstripsSat_"+str_i).c_str(),300,0,3000*pow(10,-6),300,0,3000*pow(10,-6)));
    }

    vector<vector<TH2F*>> VectNStrip_VectNStripSat_h2EvQ;
    vector<vector<TH1F*>> VectNStrip_VectNStripSat_h1EvQ;
    vector<vector<TProfile*>> VectNStrip_VectNStripSat_profEvQ;
    for(int i=1;i<6;i++)
    {
        vector<TH2F*> VectNStripSat_h2EvQ;
        vector<TH1F*> VectNStripSat_h1EvQ;
        vector<TProfile*> VectNStripSat_profEvQ;
        for(int j=0;j<i;j++)
        {
            string str_i = to_string(i);
            string str_j = to_string(j);
            TH2F* h2EvQVectVect = new TH2F(("H2 NStrip="+str_i+" & NStripSat="+str_j+" | TOB1").c_str(),("H2 NStrip="+str_i+" & NStripSat="+str_j+" | TOB1").c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
            TH1F* h1EvQVectVect = new TH1F(("H1 NStrip="+str_i+" & NStripSat="+str_j+" | TOB1").c_str(),("H1 NStrip="+str_i+" & NStripSat="+str_j+" | TOB1").c_str(),200,0,8);
            TProfile* profEvQVectVect = new TProfile(("prof NStrip="+str_i+" & NStripSat="+str_j+" | TOB1").c_str(),("prof NStrip="+str_i+" & NStripSat="+str_j+" | TOB1").c_str(),50,0,1500*pow(10,-6));
            VectNStripSat_h2EvQ.push_back(h2EvQVectVect);
            VectNStripSat_h1EvQ.push_back(h1EvQVectVect);
            VectNStripSat_profEvQ.push_back(profEvQVectVect);
        }
        VectNStrip_VectNStripSat_h2EvQ.push_back(VectNStripSat_h2EvQ);
        VectNStrip_VectNStripSat_h1EvQ.push_back(VectNStripSat_h1EvQ);
        VectNStrip_VectNStripSat_profEvQ.push_back(VectNStripSat_profEvQ);
    }
    
    vector<TH2F*> VectPartID_h2EvQ_NoSat;
    vector<TH2F*> VectPartID_h2EvQ_Sat;
    vector<TH2F*> VectPartID_h2EvQ_Sat254;
    vector<TH2F*> VectPartID_h2EvQ_Sat255;

    for(int i=0;i<4;i++)
    {
        VectPartID_h2EvQ_NoSat.push_back(new TH2F(("NoSat_"+LoopPartID(i+1)).c_str(),("NoSat_"+LoopPartID(i+1)).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
        VectPartID_h2EvQ_Sat.push_back(new TH2F(("Sat_"+LoopPartID(i+1)).c_str(),("Sat_"+LoopPartID(i+1)).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
        VectPartID_h2EvQ_Sat254.push_back(new TH2F(("Sat254_"+LoopPartID(i+1)).c_str(),("Sat254_"+LoopPartID(i+1)).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
        VectPartID_h2EvQ_Sat255.push_back(new TH2F(("Sat255_"+LoopPartID(i+1)).c_str(),("Sat255_"+LoopPartID(i+1)).c_str(),300,0,3000*pow(10,-6),300,0,3000*pow(10,-6)));
    }

    vector<TH2F*> VectLayer_h2EvQ_Tot;
    vector<TH2F*> VectLayer_h2EvQ_NoSat;
    vector<TH2F*> VectLayer_h2EvQ_Sat;
    vector<TH2F*> VectLayer_h2EvQ_Sat254;
    vector<TH2F*> VectLayer_h2EvQ_Sat255;

    vector<TProfile*> VectLayer_profEvQ_Tot;
    vector<TProfile*> VectLayer_profEvQ_NoSat;
    vector<TProfile*> VectLayer_profEvQ_Sat;
    vector<TProfile*> VectLayer_profEvQ_Sat254;
    vector<TProfile*> VectLayer_profEvQ_Sat255;

    vector<TH1F*> VectLayer_h1EvQ_Tot;
    vector<TH1F*> VectLayer_h1EvQ_NoSat;
    vector<TH1F*> VectLayer_h1EvQ_Sat;
    vector<TH1F*> VectLayer_h1EvQ_Sat254;
    vector<TH1F*> VectLayer_h1EvQ_Sat255;

    for(int i=1;i<22;i++)
    {
        VectLayer_h2EvQ_Tot.push_back(new TH2F(("h2_"+Label(i)+"_tot").c_str(),("h2_"+Label(i)+"_tot").c_str(),300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        VectLayer_h1EvQ_Tot.push_back(new TH1F(("h1_"+Label(i)+"_tot").c_str(),("h1_"+Label(i)+"_tot").c_str(),200,0,8));
        VectLayer_profEvQ_Tot.push_back(new TProfile(("pr_"+Label(i)+"_tot").c_str(),("pr_"+Label(i)+"_tot").c_str(),50,0,900*pow(10,-6),""));

        VectLayer_h2EvQ_NoSat.push_back(new TH2F(("h2_"+Label(i)+"_NoSat").c_str(),("h2_"+Label(i)+"_NoSat").c_str(),300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        VectLayer_h1EvQ_NoSat.push_back(new TH1F(("h1_"+Label(i)+"_NoSat").c_str(),("h1_"+Label(i)+"_NoSat").c_str(),200,0,8));
        VectLayer_profEvQ_NoSat.push_back(new TProfile(("pr_"+Label(i)+"_NoSat").c_str(),("pr_"+Label(i)+"_NoSat").c_str(),50,0,900*pow(10,-6),""));

        VectLayer_h2EvQ_Sat.push_back(new TH2F(("h2_"+Label(i)+"_Sat").c_str(),("h2_"+Label(i)+"_Sat").c_str(),300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        VectLayer_h1EvQ_Sat.push_back(new TH1F(("h1_"+Label(i)+"_Sat").c_str(),("h1_"+Label(i)+"_Sat").c_str(),200,0,8));
        VectLayer_profEvQ_Sat.push_back(new TProfile(("pr_"+Label(i)+"_Sat").c_str(),("pr_"+Label(i)+"_Sat").c_str(),50,0,900*pow(10,-6),""));

        VectLayer_h2EvQ_Sat254.push_back(new TH2F(("h2_"+Label(i)+"_Sat254").c_str(),("h2_"+Label(i)+"_Sat254").c_str(),300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        VectLayer_h1EvQ_Sat254.push_back(new TH1F(("h1_"+Label(i)+"_Sat254").c_str(),("h1_"+Label(i)+"_Sat254").c_str(),200,0,8));
        VectLayer_profEvQ_Sat254.push_back(new TProfile(("pr_"+Label(i)+"_Sat254").c_str(),("pr_"+Label(i)+"_Sat254").c_str(),50,0,900*pow(10,-6),""));

        VectLayer_h2EvQ_Sat255.push_back(new TH2F(("h2_"+Label(i)+"_Sat255").c_str(),("h2_"+Label(i)+"_Sat255").c_str(),300,0,3000*pow(10,-6),300,0,3000*pow(10,-6)));
        VectLayer_h1EvQ_Sat255.push_back(new TH1F(("h1_"+Label(i)+"_Sat255").c_str(),("h1_"+Label(i)+"_Sat255").c_str(),200,0,8));
        VectLayer_profEvQ_Sat255.push_back(new TProfile(("pr_"+Label(i)+"_Sat255").c_str(),("pr_"+Label(i)+"_Sat255").c_str(),50,0,3000*pow(10,-6),""));
    }

    vector<TH1F*> VectPartID_Pt;

    VectPartID_Pt.push_back(new TH1F("h1_proton_pt","h1_proton_pt",100,0,3000));
    VectPartID_Pt.push_back(new TH1F("h1_pion_pt","h1_pion_pt",100,0,3000));
    VectPartID_Pt.push_back(new TH1F("h1_gluino-u-dbar_pt","h1_gluino-u-dbar_pt",100,0,3000));
    VectPartID_Pt.push_back(new TH1F("h1_R-hadron_pt","h1_R-hadron_pt",100,0,3000));

    vector<TH1F*> VectPartID_Eloss;

    VectPartID_Eloss.push_back(new TH1F("h_proton_eloss","h_proton_eloss",300,0,3000*pow(10,-6)));
    VectPartID_Eloss.push_back(new TH1F("h_pion_eloss","h_pion_eloss",300,0,3000*pow(10,-6)));
    VectPartID_Eloss.push_back(new TH1F("h_gluino-u-dbar_eloss","h_gluino-u-dbar_eloss",300,0,3000*pow(10,-6)));
    VectPartID_Eloss.push_back(new TH1F("h_R-hadron_eloss","h_R-hadron_eloss",300,0,3000*pow(10,-6)));

    vector<TH1F*> VectPartID_PoverM;

    VectPartID_PoverM.push_back(new TH1F("h_proton_pm","h_proton_pm",1000,0.01,10));
    VectPartID_PoverM.push_back(new TH1F("h_pion_pm","h_pion_pm",1000,0.01,10));
    VectPartID_PoverM.push_back(new TH1F("h_gluino-u-dbar_pm","h_gluino-u-dbar_pm",1000,0.01,10));
    VectPartID_PoverM.push_back(new TH1F("h_R-hadron_pm","h_R-hadron_pm",1000,0.01,10));


    TH1F* hEmpty = new TH1F("hEmpty","hEmpty",21,0,21); //pour tracer les lignes avec la fonction SetHistoLabel

    

    TFile* fileout = new TFile(s2.c_str(),"RECREATE");


    int entries = 1000;


    bool testsat254 = false; 
    bool testsat255 = false;

    int clust=0,clustsat254=0,clustsat255=0;

    for(int i=0;i<entries;i++)
    {
        if(i%1000==0) cout<<"Event "<<i<<endl;
        b1->GetEntry(i);
        for(int track=0;track<b1->GetNtracks();track++)
        {
            vector<int> vect_partID;
            vector<float> vect_eloss;
            testsat254              = false;
            testsat255              = false;
            float pt                = b1->GetVectTrack()[track].GetPt();
            float p                 = b1->GetVectTrack()[track].GetP();
            int NCluster            = b1->GetVectTrack()[track].GetNCluster();
            int NClustSat254        = b1->GetVectTrack()[track].GetNSatCluster(254);
            int NClustSat255        = b1->GetVectTrack()[track].GetNSatCluster(255);
            
            for(int cluster=0;cluster<b1->GetVectTrack()[track].GetNCluster();cluster++)
            {
                float charge        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSclusCharge();
                float eloss         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetEloss();
                int layer           = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetLayer();
                int layerLabel      = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetLayerLabel();
                bool sat254         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSat254();
                bool sat255         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSat255();
                int nsimhits        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSimHits();
                int nstrips         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNStrip();
                int nsatboth        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStripBoth();
                int nsat254         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStrip(254);
                int nsat255         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStrip(255);
                bool edge           = b1->GetVectTrack()[track].GetVectClusters()[cluster].Edge();
                bool cut            = b1->GetVectTrack()[track].GetVectClusters()[cluster].Cut();

                clust++;

                vect_eloss.push_back(eloss);

                h2ElossvQ_tot->Fill(eloss,charge);
                hDiff_rel_ElossQ_tot->Fill((eloss-charge)/charge);
                hRatio_ElossQ_tot->Fill(eloss/charge);
                profElossVsQ_tot->Fill(eloss,charge);
                
                hLayer->Fill(layer);
                hLayerLabel->Fill(layerLabel-1);

                VectLayer_h2EvQ_Tot[layerLabel-1]->Fill(eloss,charge);
                VectLayer_h1EvQ_Tot[layerLabel-1]->Fill(eloss/charge);
                VectLayer_profEvQ_Tot[layerLabel-1]->Fill(eloss,charge);


                for(int simhit=0;simhit<nsimhits;simhit++)
                {
                    int partID      = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectSimHits()[simhit].GetPartId();
                    
                    if(partID==2212 || partID==-2212) //proton
                    {
                        if(layerLabel==5) VectPartID_Eloss[0]->Fill(eloss);
                        if(sat254==false && sat255==false) VectPartID_h2EvQ_NoSat[1]->Fill(eloss,charge);
                        if(sat254==true && GetPoverM(p,2212)<1) VectPartID_h2EvQ_Sat[1]->Fill(eloss,charge);
                        if(sat254==true && sat255==false) VectPartID_h2EvQ_Sat254[1]->Fill(eloss,charge);
                        if(sat255==true) VectPartID_h2EvQ_Sat255[1]->Fill(eloss,charge);
                    }
                    if(partID==211 || partID==-211) //pion
                    {
                        if(layerLabel==5) VectPartID_Eloss[1]->Fill(eloss);
                        if(sat254==false && sat255==false) VectPartID_h2EvQ_NoSat[0]->Fill(eloss,charge);
                        if(sat254==true) VectPartID_h2EvQ_Sat[0]->Fill(eloss,charge);
                        if(sat254==true && sat255==false) VectPartID_h2EvQ_Sat254[0]->Fill(eloss,charge);
                        if(sat255==true) VectPartID_h2EvQ_Sat255[0]->Fill(eloss,charge);
                    }
                    if(partID==1009213 || partID==-1009213) //R-hadron gluino u dbar
                    {
                        if(layerLabel==5) VectPartID_Eloss[2]->Fill(eloss);
                        if(sat254==false && sat255==false) VectPartID_h2EvQ_NoSat[2]->Fill(eloss,charge);
                        if(sat254==true) VectPartID_h2EvQ_Sat[2]->Fill(eloss,charge);
                        if(sat254==true && sat255==false) VectPartID_h2EvQ_Sat254[2]->Fill(eloss,charge);
                        if(sat255==true) VectPartID_h2EvQ_Sat255[2]->Fill(eloss,charge);
                    }
                    if((int)partID/1000==1009 || (int)partID/1000==-1009) //R-hadron
                    {
                        if(layerLabel==5) VectPartID_Eloss[3]->Fill(eloss);
                        if(sat254==false && sat255==false) VectPartID_h2EvQ_NoSat[3]->Fill(eloss,charge);
                        if(sat254==true) VectPartID_h2EvQ_Sat[3]->Fill(eloss,charge);
                        if(sat254==true && sat255==false) VectPartID_h2EvQ_Sat254[3]->Fill(eloss,charge);
                        if(sat255==true) VectPartID_h2EvQ_Sat255[3]->Fill(eloss,charge);
                    }
                    
                }

                for(int nstrip=1;nstrip<6;nstrip++)
                {
                    for(int nstripsat=0;nstripsat<nstrip;nstripsat++)
                    {
                        if(nstrip==nstrips && nstripsat==nsat254 && nsat255==0) //on se place dans TOB1 pour ne pas melanger differents effets + selection sans les bords ou les cuts
                        {
                            if(layerLabel==5 && cut==false && edge==false)
                            {
                                VectNStrip_VectNStripSat_h2EvQ[nstrip-1][nstripsat]->Fill(eloss,charge);
                                VectNStrip_VectNStripSat_h1EvQ[nstrip-1][nstripsat]->Fill(eloss/charge);
                                VectNStrip_VectNStripSat_profEvQ[nstrip-1][nstripsat]->Fill(eloss,charge);
                            }
                        }
                    }
                }

                if(sat254==true && sat255==false) 
                {
                    hRatio_ElossQ_Sat254->Fill(eloss/charge);
                    hDiff_rel_ElossQ_Sat254->Fill((eloss-charge)/charge);
                    h2ElossvQ_Sat254->Fill(eloss,charge);
                    profElossVsQ_Sat254->Fill(eloss,charge);
                    hLayerLabelSat254->Fill(layerLabel-1);
                    VectLayer_h2EvQ_Sat254[layerLabel-1]->Fill(eloss,charge);
                    VectLayer_h1EvQ_Sat254[layerLabel-1]->Fill(eloss/charge);
                    VectLayer_profEvQ_Sat254[layerLabel-1]->Fill(eloss,charge);
                    testsat254=true;
                    clustsat254++;
                    for(int strip=0;strip<10;strip++)
                    {
                        
                        if(nstrips==strip+1)
                        {
                            VectNStrip_h2EvQ_Sat254[strip]->Fill(eloss,charge);
                        }
                    }
                    for(int strip=0;strip<10;strip++)
                    {
                        
                        if(nsat254==strip+1)
                        {
                            VectNStripSat_h2EvQ_254[strip]->Fill(eloss,charge);
                        }
                    }
                }
                if(sat255==true)
                {
                    hRatio_ElossQ_Sat255->Fill(eloss/charge);
                    hDiff_rel_ElossQ_Sat255->Fill((eloss-charge)/charge);
                    h2ElossvQ_Sat255->Fill(eloss,charge);
                    profElossVsQ_Sat255->Fill(eloss,charge);
                    hLayerLabelSat255->Fill(layerLabel-1);
                    VectLayer_h2EvQ_Sat255[layerLabel-1]->Fill(eloss,charge);
                    VectLayer_h1EvQ_Sat255[layerLabel-1]->Fill(eloss/charge);
                    VectLayer_profEvQ_Sat255[layerLabel-1]->Fill(eloss,charge);
                    testsat255=true;
                    clustsat255++;
                    for(int strip=0;strip<10;strip++)
                    {
                        
                        if(nstrips==strip+1)
                        {
                            VectNStrip_h2EvQ_Sat255[strip]->Fill(eloss,charge);
                        }
                    }
                    for(int strip=0;strip<10;strip++)
                    {
                        
                        if(nsat255==strip+1)
                        {
                            VectNStripSat_h2EvQ_255[strip]->Fill(eloss,charge);
                        }
                    }
                }
                if(sat254==true)
                {
                    h2ElossvQ_Sat->Fill(eloss,charge);
                    hDiff_rel_ElossQ_Sat->Fill((eloss-charge)/charge);
                    hRatio_ElossQ_Sat->Fill(eloss/charge);
                    profElossVsQ_Sat->Fill(eloss,charge);
                    VectLayer_h2EvQ_Sat[layerLabel-1]->Fill(eloss,charge);
                    VectLayer_h1EvQ_Sat[layerLabel-1]->Fill(eloss/charge);
                    VectLayer_profEvQ_Sat[layerLabel-1]->Fill(eloss,charge);
                    for(int strip=0;strip<10;strip++)
                    {
                        
                        if(nsatboth==strip+1)
                        {
                            VectNStrip_h2EvQ_Sat[strip]->Fill(eloss,charge);
                        }
                    }
                }
                if(sat254==false && sat255==false)
                {
                    h2ElossvQ_NoSat->Fill(eloss,charge);
                    hDiff_rel_ElossQ_NoSat->Fill((eloss-charge)/charge);
                    hRatio_ElossQ_NoSat->Fill(eloss/charge);
                    profElossVsQ_NoSat->Fill(eloss,charge);
                    VectLayer_h2EvQ_NoSat[layerLabel-1]->Fill(eloss,charge);
                    VectLayer_h1EvQ_NoSat[layerLabel-1]->Fill(eloss/charge);
                    VectLayer_profEvQ_NoSat[layerLabel-1]->Fill(eloss,charge);
                    for(int strip=0;strip<10;strip++)
                    {
                        if(nstrips==strip+1)
                        {
                            VectNStrip_h2EvQ_NoSat[strip]->Fill(eloss,charge);
                        }
                    }
                }
            }
            
            hPt_tot->Fill(pt);
            hNClusterPerTrack->Fill(NCluster);
            hNClusterSat254PerTrack->Fill(NClustSat254);
            hNClusterSat255PerTrack->Fill(NClustSat255);
            hRatio_NClusterSat254->Fill((double)NClustSat254/(double)NCluster);
            hRatio_NClusterSat255->Fill((double)NClustSat255/(double)NCluster);

            float threshold=0;

            int id=GetPartID(b1->GetVectTrack()[track].GetVectClusters(),threshold);
            
            if(id==2212 || id==-2212) 
            {
                VectPartID_Pt[0]->Fill(pt);
                VectPartID_PoverM[0]->Fill(GetPoverM(p,id));
            }
            if(id==211 || id==-211) 
            {
                VectPartID_Pt[1]->Fill(pt);
                VectPartID_PoverM[1]->Fill(GetPoverM(p,id));
            }
            if(id==1009213 || id==-1009213) 
            {
                VectPartID_Pt[2]->Fill(pt);
                VectPartID_PoverM[2]->Fill(GetPoverM(p,id));
            }
            if((int)id/1000==1009 || (int)id/1000==-1009)
            {
                VectPartID_Pt[3]->Fill(pt);
                VectPartID_PoverM[3]->Fill(GetPoverM(p,id));
            }
            
            
            if(testsat254) 
            {
                hPt_Sat254->Fill(pt);
            } 
            if(testsat255)
            {
                hPt_Sat255->Fill(pt);
            }
            if(testsat254==true || testsat255==true)
            {
                hPt_Sat->Fill(pt);
            }
            if(testsat254==false && testsat255==false)
            {
                hPt_NoSat->Fill(pt);
            }

            Estimator estim(vect_eloss);

            h2RatioSatEloss254->Fill(estim.GetMean(),(double)NClustSat254/(double)NCluster);
            profSatEloss254->Fill(estim.GetMean(),(double)NClustSat254/(double)NCluster);
            h2RatioSatPt254->Fill(pt,(double)NClustSat254/(double)NCluster);
            profSatPt254->Fill(pt,(double)NClustSat254/(double)NCluster);
            h2RatioSatPartID254->Fill(ReBinPartID(id),(double)NClustSat254/(double)NCluster);
            profSatPartID254->Fill(ReBinPartID(id),(double)NClustSat254/(double)NCluster);

            h2RatioSatEloss255->Fill(estim.GetMean(),(double)NClustSat255/(double)NCluster);
            profSatEloss255->Fill(estim.GetMean(),(double)NClustSat255/(double)NCluster);
            h2RatioSatPt255->Fill(pt,(double)NClustSat255/(double)NCluster);
            profSatPt255->Fill(pt,(double)NClustSat255/(double)NCluster);
            h2RatioSatPartID255->Fill(ReBinPartID(id),(double)NClustSat255/(double)NCluster);
            profSatPartID255->Fill(ReBinPartID(id),(double)NClustSat255/(double)NCluster);

            /*if(testsat254 && testsat255 )//&& b1->GetVectTrack()[track].GetVectClusters().size() <=3)
            {
                TCanvas* c_profClust = new TCanvas();
                c_profClust = &b1->GetVectTrack()[track].GetProfCluster();
                string str_pt = to_string(pt);
                TText* text_pt = new TText(c_profClust->GetWw()-100,c_profClust->GetWw()-100,str_pt.c_str());
                text_pt->Draw("SAME");
                c_profClust->SaveAs("../data/ProfileClust.pdf");
                getchar();
            }*/
        }
    }

    /*Builder* b2 = new Builder(chain);
    b2->SetBranchAdd();
    for(int i=0;i<entries;i++)
    {
        if(i%1000==0) cout<<"Event "<<i<<endl;
        b2->SetCalibration(factor,i);
        for(int track=0;track<b2->GetNtracks();track++)
        {
            for(int cluster=0;cluster<b2->GetVectTrack()[track].GetVectClusters().size();cluster++)
            {
                float charge = b2->GetVectTrack()[track].GetVectClusters()[cluster].GetSclusCharge();
                float eloss = b2->GetVectTrack()[track].GetVectClusters()[cluster].GetEloss();
                calibProf->Fill(eloss,charge);
                calibh2->Fill(eloss,charge);
                calibh1->Fill((eloss-charge)/charge);
            }
        }
    }

    //TFitResultPtr rcalib = calibProf->Fit("pol1","S");
    TFitResultPtr calibGaus = calibh1->Fit("gaus","S");*/

    TCanvas* cPtPartID = new TCanvas("cPtPartID","cPtPartID",700,400);
    vector<string> VectLegend;
    VectLegend.push_back("p");
    VectLegend.push_back("#pi");
    VectLegend.push_back("R^{+}_{#scale[0.6]{#tilde{g}u#bar{d}}}");
    VectLegend.push_back("R^{+}");
    StackHisto(*cPtPartID,VectPartID_Pt,VectLegend,"Distribution en p_{T} pour differentes particules","p_{T}");
    cPtPartID->Write();

    TCanvas* cElossPartID = new TCanvas("cElossPartID","cElossPartID",700,400);
    StackHisto(*cElossPartID,VectPartID_Eloss,VectLegend,"Distribution E_{loss} pour differentes particules dans TOB1","E_{loss}");
    cElossPartID->Write();

    TCanvas* cPoverMPartID = new TCanvas("cPoverMPartID","cPoverMPartID",700,400);
    StackHisto(*cPoverMPartID,VectPartID_PoverM,VectLegend,"Distribution du P/M pour differentes particules","P/M");
    cPoverMPartID->Write();

    TCanvas* cPoverMPartIDNormal = new TCanvas("cPoverMPartIDNormal","cPoverMPartIDNormal",700,400);
    DrawHistoNormalized(*cPoverMPartIDNormal,VectPartID_PoverM,VectLegend,"Distribution du P/M pour differentes particules","P/M");
    cPoverMPartIDNormal->Write();


    DrawHisto(*fileout,h2RatioSatPt254,"Frac Sat254 = f(Pt)","Pt","#frac{#ClustSat}{#Clust}");
    DrawHisto(*fileout,h2RatioSatEloss254,"Frac Sat254 = f(E_{loss})","E_{loss}","#frac{#ClustSat}{#Clust}");
    DrawHisto(*fileout,h2RatioSatPartID254,"Frac Sat254 = f(PartID)","PartID","#frac{#ClustSat}{#Clust}");
    
    TCanvas* cRatioSatPt254 = new TCanvas("cRatioSatPt254","cRatioSatPt254",700,400);
    h2RatioSatPt254->SetMaximum(200);
    h2RatioSatPt254->Draw();
    h2RatioSatPt254->SetDrawOption("COLZ");
    profSatPt254->SetLineColor(kRed);
    profSatPt254->SetMarkerColor(kRed);
    profSatPt254->SetMarkerStyle(2);
    profSatPt254->Draw("SAME");
    cRatioSatPt254->Write();
    profSatPt254->Write();

    TCanvas* cRatioSatEloss254 = new TCanvas("cRatioSatEloss254","cRatioSatEloss254",700,400);
    h2RatioSatEloss254->Draw();
    h2RatioSatEloss254->SetDrawOption("COLZ");
    profSatEloss254->SetMarkerStyle(2);
    profSatEloss254->SetMarkerColor(kRed);
    profSatEloss254->SetLineColor(kRed);
    profSatEloss254->Draw("SAME");
    cRatioSatEloss254->Write();
    profSatEloss254->Write();

    TCanvas* canvas_RatioSatPartID254 = new TCanvas("canvas_RatioSatPartID254","canvas_RatioSatPartID254",700,400);
    SetHistoLabelPartID(canvas_RatioSatPartID254,h2RatioSatPartID254);
    profSatPartID254->SetMarkerStyle(2);
    profSatPartID254->SetMarkerColor(kRed);
    profSatPartID254->SetLineColor(kRed);
    profSatPartID254->Draw("SAME");
    h2RatioSatPartID254->GetXaxis()->SetTitle("");
    h2RatioSatPartID254->SetDrawOption("COLZ");
    h2RatioSatPartID254->SetMaximum(800);
    canvas_RatioSatPartID254->Write();

    DrawHisto(*fileout,h2RatioSatPt255,"Frac Sat255 = f(Pt)","Pt","#frac{#ClustSat}{#Clust}");
    DrawHisto(*fileout,h2RatioSatEloss255,"Frac Sat255 = f(E_{loss})","E_{loss}","#frac{#ClustSat}{#Clust}");
    DrawHisto(*fileout,h2RatioSatPartID255,"Frac Sat255 = f(PartID)","PartID","#frac{#ClustSat}{#Clust}");
    
    TCanvas* cRatioSatPt255 = new TCanvas("cRatioSatPt255","cRatioSatPt255",700,400);
    h2RatioSatPt255->SetMaximum(200);
    h2RatioSatPt255->Draw();
    h2RatioSatPt255->SetDrawOption("COLZ");
    profSatPt255->SetLineColor(kRed);
    profSatPt255->SetMarkerColor(kRed);
    profSatPt255->SetMarkerStyle(2);
    profSatPt255->Draw("SAME");
    cRatioSatPt255->Write();
    profSatPt255->Write();

    TCanvas* cRatioSatEloss255 = new TCanvas("cRatioSatEloss255","cRatioSatEloss255",700,400);
    h2RatioSatEloss255->Draw();
    h2RatioSatEloss255->SetDrawOption("COLZ");
    profSatEloss255->SetMarkerStyle(2);
    profSatEloss255->SetMarkerColor(kRed);
    profSatEloss255->SetLineColor(kRed);
    profSatEloss255->Draw("SAME");
    cRatioSatEloss255->Write();
    profSatEloss255->Write();

    TCanvas* canvas_RatioSatPartID255 = new TCanvas("canvas_RatioSatPartID255","canvas_RatioSatPartID255",700,400);
    SetHistoLabelPartID(canvas_RatioSatPartID255,h2RatioSatPartID255);
    profSatPartID255->SetMarkerStyle(2);
    profSatPartID255->SetMarkerColor(kRed);
    profSatPartID255->SetLineColor(kRed);
    profSatPartID255->Draw("SAME");
    h2RatioSatPartID255->GetXaxis()->SetTitle("");
    h2RatioSatPartID255->SetDrawOption("COLZ");
    h2RatioSatPartID255->SetMaximum(800);
    canvas_RatioSatPartID255->Write();


    DrawHisto(*fileout,hDiff_rel_ElossQ_tot,"Difference relative","(E-Q)/Q");
    DrawHisto(*fileout,hDiff_rel_ElossQ_NoSat,"Difference relative","(E-Q)/Q");
    DrawHisto(*fileout,hDiff_rel_ElossQ_Sat,"Difference relative","(E-Q)/Q");
    DrawHisto(*fileout,hDiff_rel_ElossQ_Sat254,"Difference relative","(E-Q)/Q");
    DrawHisto(*fileout,hDiff_rel_ElossQ_Sat255,"Difference relative","(E-Q)/Q");

    DrawHisto(*fileout,hRatio_ElossQ_tot,"","E/Q");
    DrawHisto(*fileout,hRatio_ElossQ_NoSat,"","E/Q");
    DrawHisto(*fileout,hRatio_ElossQ_Sat,"","E/Q");
    DrawHisto(*fileout,hRatio_ElossQ_Sat254,"","E/Q");
    DrawHisto(*fileout,hRatio_ElossQ_Sat255,"","E/Q");

    DrawHisto(*fileout,h2ElossvQ_tot,"E_{loss} v. Q","E_{loss}","Q");
    DrawHisto(*fileout,h2ElossvQ_NoSat,"E_{loss} v. Q","E_{loss}","Q");
    DrawHisto(*fileout,h2ElossvQ_Sat,"E_{loss} v. Q","E_{loss}","Q");
    DrawHisto(*fileout,h2ElossvQ_Sat254,"E_{loss} v. Q","E_{loss}","Q");
    DrawHisto(*fileout,h2ElossvQ_Sat255,"E_{loss} v. Q","E_{loss}","Q");

    fileout->Append(profElossVsQ_tot);
    fileout->Append(profElossVsQ_NoSat);
    fileout->Append(profElossVsQ_Sat);
    fileout->Append(profElossVsQ_Sat254);
    fileout->Append(profElossVsQ_Sat255);

    fileout->Append(hLayer);
    fileout->Append(hLayerLabel);
    fileout->Append(hLayerLabelSat254);
    fileout->Append(hLayerLabelSat255);

    DrawHisto(*fileout,hPt_tot,"Distribution du Pt","Pt");
    DrawHisto(*fileout,hPt_NoSat,"Distribution du Pt","Pt");
    DrawHisto(*fileout,hPt_Sat,"Distribution du Pt","Pt");
    DrawHisto(*fileout,hPt_Sat254,"Distribution du Pt","Pt");
    DrawHisto(*fileout,hPt_Sat255,"Distribution du Pt","Pt");

    DrawHisto(*fileout,hNClusterPerTrack,"Nombre de clusters pour une trace","# of clusters");
    DrawHisto(*fileout,hNClusterSat254PerTrack,"Nombre de clusters sat254 pour une trace","# of clusters");
    DrawHisto(*fileout,hNClusterSat255PerTrack,"Nombre de clusters sat255 pour une trace","# of clusters");

    fileout->Append(hRatio_NClusterSat254);
    fileout->Append(hRatio_NClusterSat255);

// ------------------------------------------ ratio saturation fonction de la layer 

    TEfficiency* Ratio_LayerSat254 = new TEfficiency(*hLayerLabelSat254,*hLayerLabel);
    TEfficiency* Ratio_LayerSat255 = new TEfficiency(*hLayerLabelSat255,*hLayerLabel);

    fileout->Append(Ratio_LayerSat254);
    fileout->Append(Ratio_LayerSat255);

    float rat254 = (float)clustsat254/(float)clust;
    float rat255 = (float)clustsat255/(float)clust;
    cout<<"rat254 "<<rat254<<endl;
    cout<<"rat255 "<<rat255<<endl;
    TLine* line254 = new TLine(0,rat254,21,rat254);
    TLine* line255 = new TLine(0,rat255,21,rat255);

    TCanvas* canvas_ratioSatLayer = new TCanvas("canvas_ratioSatLayer","canvas_ratioSatLayer",700,400);
    SetHistoLabel(canvas_ratioSatLayer,hEmpty);
    canvas_ratioSatLayer->cd();
    hEmpty->GetYaxis()->SetRangeUser(pow(10,-4),0.25);
    Ratio_LayerSat254->Draw("SAME");
    Ratio_LayerSat255->SetLineColor(2);
    Ratio_LayerSat255->Draw("SAME");
    line254->SetLineStyle(3);
    line254->Draw("SAME");
    line255->SetLineStyle(3);
    line255->SetLineColor(2);
    line255->Draw("SAME");
    canvas_ratioSatLayer->SetLogy();
    TLegend* legratio = new TLegend(0.6,0.8,0.9,0.9);
    legratio->AddEntry(line254,("Saturation 254 : "+to_string(rat254)).c_str(),"l");
    legratio->AddEntry(line255,("Saturation 255 : "+to_string(rat255)).c_str(),"l");
    legratio->Draw("SAME");
    canvas_ratioSatLayer->Write();

// ------------------------------------------

// ------------------------------------------ TH2F + TProfile 

    TCanvas* canvas_EvQ_h2prof_tot = new TCanvas("canvas_EvQ_h2prof_tot","canvas_EvQ_h2prof_tot",700,400);
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_tot,h2ElossvQ_tot,profElossVsQ_tot,"","E_{loss}","Q");
    canvas_EvQ_h2prof_tot->Write();

    TCanvas* canvas_EvQ_h2prof_NoSat = new TCanvas("canvas_EvQ_h2prof_NoSat","canvas_EvQ_h2prof_NoSat",700,400);
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_NoSat,h2ElossvQ_NoSat,profElossVsQ_NoSat,"","E_{loss}","Q");
    canvas_EvQ_h2prof_NoSat->Write();

    TCanvas* canvas_EvQ_h2prof_Sat = new TCanvas("canvas_EvQ_h2prof_Sat","canvas_EvQ_h2prof_Sat",700,400);
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_Sat,h2ElossvQ_Sat,profElossVsQ_Sat,"","E_{loss}","Q");
    canvas_EvQ_h2prof_Sat->Write();

    TCanvas* canvas_EvQ_h2prof_Sat254 = new TCanvas("canvas_EvQ_h2prof_Sat254","canvas_EvQ_h2prof_Sat254",700,400);
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_Sat254,h2ElossvQ_Sat254,profElossVsQ_Sat254,"","E_{loss}","Q");
    canvas_EvQ_h2prof_Sat254->Write();

    TCanvas* canvas_EvQ_h2prof_Sat255 = new TCanvas("canvas_EvQ_h2prof_Sat255","canvas_EvQ_h2prof_Sat255",700,400);
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_Sat255,h2ElossvQ_Sat255,profElossVsQ_Sat255,"","E_{loss}","Q");
    canvas_EvQ_h2prof_Sat255->Write();

    profElossVsQ_tot->SetMarkerStyle(1);
    profElossVsQ_NoSat->SetMarkerStyle(1);
    profElossVsQ_Sat->SetMarkerStyle(1);
    profElossVsQ_Sat254->SetMarkerStyle(1);
    profElossVsQ_Sat255->SetMarkerStyle(1);

// ------------------------------------------

// ------------------------------------------ TH1F Pt de la trace, avec ou sans saturation

    TCanvas* canvas_Pt = new TCanvas("canvas_Pt","canvas_Pt",700,400);
    vector<TH1F*> VectHistoPt;
    VectHistoPt.push_back(hPt_NoSat);
    VectHistoPt.push_back(hPt_Sat);
    vector<string> VectLegendPt;
    VectLegendPt.push_back("No saturation");
    VectLegendPt.push_back("Saturation");
    DrawHistoNormalized(*canvas_Pt,VectHistoPt,VectLegendPt,"Distribution de P_{T}","P_{T}");
    canvas_Pt->Write();

// ------------------------------------------

// ------------------------------------------ TH1F ratio NClusters qui saturent par trace

    TCanvas* canvas_RatioNClusterSat = new TCanvas("canvas_RatioNClusterSat","canvas_RatioNClusterSat",700,400);
    vector<TH1F*> VectHistoNSat;
    VectHistoNSat.push_back(hRatio_ElossQ_Sat254);
    VectHistoNSat.push_back(hRatio_ElossQ_Sat255);
    vector<string> VectLegendNSat;
    VectLegend.push_back("Saturation 254");
    VectLegend.push_back("Saturation 255");
    DrawHistoNormalized(*canvas_RatioNClusterSat,VectHistoNSat,VectLegendNSat,"Fraction de saturation","#frac{Nombre de cluster sature}{Nombre de cluster dans la trace}");
    canvas_RatioNClusterSat->Write();

// ------------------------------------------

// ------------------------------------------ Ratio E/Q

    TCanvas* canvas_EvQ_Ratio = new TCanvas("canvas_EvQ_Ratio","canvas_EvQ_Ratio",700,400);
    vector<TH1F*> VectHistoEoverQ;
    VectHistoEoverQ.push_back(hRatio_ElossQ_NoSat);
    VectHistoEoverQ.push_back(hRatio_ElossQ_Sat);
    vector<string> VectLegendEoverQ;
    VectLegendEoverQ.push_back("No saturation");
    VectLegendEoverQ.push_back("Saturation");
    DrawHistoNormalized(*canvas_EvQ_Ratio,VectHistoEoverQ,VectLegendEoverQ,"","#frac{E}{Q}");
    canvas_EvQ_Ratio->Write();

// ------------------------------------------

// ------------------------------------------ Fit E/Q

    TF1* f_RatioEvQ = new TF1("f_RatioEvQ","gaus(0)+gaus(3)+gaus(6)+[11]*TMath::Landau(x,[9],[10])",0,8);
    //f_RatioEvQ->SetParameters(1,0.5,0.1,1,1,1,1,2);
    //TF1* f_RatioEvQ = new TF1("f_RatioEvQ","gaus(0)+TMath::Landau(x,[3],[4])",0,8);
    f_RatioEvQ->SetParameters(1,1.03,.09,1,0.45,0.1,1,0.55,1,1);
    TFitResultPtr Fit_RatioEvQ = hRatio_ElossQ_tot->Fit(f_RatioEvQ,"S");
   //hRatio_ElossQ_tot->Write();
    TCanvas* ctest = new TCanvas();
    hRatio_ElossQ_tot->Draw();
    ctest->SetLogy();

// ------------------------------------------

// ------------------------------------------ Fit TProfile

    TFitResultPtr Fit_Profile_EvQ_tot = profElossVsQ_tot->Fit("pol1","S","",0.1*pow(10,-3),0.4*pow(10,-3));
    TFitResultPtr Fit_Profile_EvQ_NoSat = profElossVsQ_NoSat->Fit("pol1","S","",0.2*pow(10,-3),0.5*pow(10,-3));

// ------------------------------------------


/*    for(int i=0;i<21;i++)
    {
        string title = "E_{loss} v. Q  |  "+Label(i+1);
        TCanvas* csuperposed = new TCanvas(title.c_str(),title.c_str(),700,400);
        DrawHisto(*fileout,VectLayer_h2EvQ_Tot[i],title.c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectLayer_h2EvQ_NoSat[i],title.c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectLayer_h2EvQ_Sat[i],title.c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectLayer_h2EvQ_Sat254[i],title.c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectLayer_h2EvQ_Sat255[i],title.c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectLayer_h1EvQ_Tot[i],title.c_str(),"E/Q");
        DrawHisto(*fileout,VectLayer_h1EvQ_NoSat[i],title.c_str(),"E/Q");
        DrawHisto(*fileout,VectLayer_h1EvQ_Sat[i],title.c_str(),"E/Q");
        DrawHisto(*fileout,VectLayer_h1EvQ_Sat254[i],title.c_str(),"E/Q");
        DrawHisto(*fileout,VectLayer_h1EvQ_Sat255[i],title.c_str(),"E/Q");
        SuperposedHisto2DProfile(*csuperposed,VectLayer_h2EvQ_Tot[i],VectLayer_profEvQ_Tot[i],title.c_str(),"E_{loss}","Q");
        csuperposed->Write();
        SuperposedHisto2DProfile(*csuperposed,VectLayer_h2EvQ_NoSat[i],VectLayer_profEvQ_NoSat[i],title.c_str(),"E_{loss}","Q");
        csuperposed->Write();
        SuperposedHisto2DProfile(*csuperposed,VectLayer_h2EvQ_Sat[i],VectLayer_profEvQ_Sat[i],title.c_str(),"E_{loss}","Q");
        csuperposed->Write();
        SuperposedHisto2DProfile(*csuperposed,VectLayer_h2EvQ_Sat254[i],VectLayer_profEvQ_Sat254[i],title.c_str(),"E_{loss}","Q");
        csuperposed->Write();
        SuperposedHisto2DProfile(*csuperposed,VectLayer_h2EvQ_Sat255[i],VectLayer_profEvQ_Sat255[i],title.c_str(),"E_{loss}","Q");
        csuperposed->Write();
    }

    for(int i=0;i<10;i++)
    {
        string title = "NStrip_"+to_string(i+1);
        DrawHisto(*fileout,VectNStrip_h2EvQ_NoSat[i],("NoSat_"+title).c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectNStrip_h2EvQ_Sat[i],("Sat_"+title).c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectNStrip_h2EvQ_Sat254[i],("Sat254_"+title).c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectNStrip_h2EvQ_Sat255[i],("Sat255_"+title).c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectNStripSat_h2EvQ_254[i],("Sat254_"+title+"Sat").c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectNStripSat_h2EvQ_255[i],("Sat255_"+title+"Sat").c_str(),"E_{loss}","Q");
    }

    for(int i=0;i<4;i++)
    {
        DrawHisto(*fileout,VectPartID_h2EvQ_NoSat[i],("NoSat_"+LoopPartID(i+1)).c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectPartID_h2EvQ_Sat[i],("Sat_"+LoopPartID(i+1)).c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectPartID_h2EvQ_Sat254[i],("Sat254_"+LoopPartID(i+1)).c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectPartID_h2EvQ_Sat255[i],("Sat255_"+LoopPartID(i+1)).c_str(),"E_{loss}","Q");
    }

*/

    for(int nstrip=1;nstrip<6;nstrip++)
    {
        TCanvas* cNStripNStripSat_h2prof = new TCanvas(("H2 & profile NStrip="+to_string(nstrip)+" | TOB1").c_str(),("H2 & profile NStrip="+to_string(nstrip)+" | TOB1").c_str());
        TCanvas* cNStripNStripSat_h1 = new TCanvas(("H1 NStrip="+to_string(nstrip)+" | TOB1").c_str(),("H1 NStrip="+to_string(nstrip)+" | TOB1").c_str());
        TLine* linexy = new TLine(0,0,0.0015,0.0015);
        TLine* linesat254 = new TLine(0,254*(3.61*pow(10,-9)*247),0.0015,254*(3.61*pow(10,-9)*247));
        linesat254->SetLineStyle(4);
        linexy->SetLineStyle(3);
        cNStripNStripSat_h2prof->Divide(2,nstrip);
        cNStripNStripSat_h1->Divide(2,nstrip);
        for(int nstripsat=0;nstripsat<nstrip;nstripsat++)
        {
            cNStripNStripSat_h2prof->cd(nstripsat+1);
            SuperposedHisto2DProfile(*cNStripNStripSat_h2prof,VectNStrip_VectNStripSat_h2EvQ[nstrip-1][nstripsat],VectNStrip_VectNStripSat_profEvQ[nstrip-1][nstripsat],("H2 & prof NStrip="+to_string(nstrip)+" & NStripSat="+to_string(nstripsat)+" | TOB1").c_str(),"E_{loss}","Q");
            //VectNStrip_VectNStripSat_h2EvQ[nstrip-1][nstripsat]->Draw("COLZ");
            linexy->Draw("SAME");
            linesat254->Draw("SAME");
            cNStripNStripSat_h1->cd(nstripsat+1);
            VectNStrip_VectNStripSat_h1EvQ[nstrip-1][nstripsat]->Draw();
            //DrawHisto(*fileout,VectNStrip_VectNStripSat_h2EvQ[nstrip-1][nstripsat],("H2 NStrip="+to_string(nstrip)+" & NStripSat="+to_string(nstripsat)).c_str(),"E_{loss}","Q");
            //DrawHisto(*fileout,VectNStrip_VectNStripSat_h1EvQ[nstrip-1][nstripsat],("H1 NStrip="+to_string(nstrip)+" & NStripSat="+to_string(nstripsat)).c_str(),"E_{loss}");

        }
        cNStripNStripSat_h2prof->Write();
        cNStripNStripSat_h1->Write();
    }

    fileout->Write();
    fileout->Close();

    delete fileout;

    system(("mv "+s2+" ./data/.").c_str());
    

    return 0;
}
