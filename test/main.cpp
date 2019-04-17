#include <iostream>
#include <fstream>
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
#include "../interface/CalibCharge.h"
#include "../interface/Calibration.h"

using namespace std;

int main(int argc,char** argv){

    TChain chain;
	chain.SetName("stage/ttree");
	for(int i=1;i<argc;i++) chain.Add(argv[i]);

    string s1 = argv[1];
    string s2 = s1.substr(0,s1.find('.'))+"_results.root";

    vector<vector<vector<vector<TH2F*>>>> VectLayerVectNStripVectNStripSatVectSatTypeProfileFit;
    for(int layer=1;layer<11;layer++)
    {
        vector<vector<vector<TH2F*>>> VectNStripVectNStripSatVectSatTypeProfileFit;
        for(int nstrip=3;nstrip<7;nstrip++)
        {
            vector<vector<TH2F*>> VectNStripSatVectSatTypeProfileFit;
            for(int nstripsat=1;nstripsat<3;nstripsat++)
            {
                vector<TH2F*> VectSatTypeProfileFit;
                for(int sat254=0;sat254<2;sat254++)
                {
                    string SatType="254";
                    if(sat254==1) SatType="255";
                    if(sat254==0) VectSatTypeProfileFit.push_back(new TH2F((Label(layer)+" NStrip="+to_string(nstrip)+" & NStripSat="+to_string(nstripsat)+" Sat"+SatType).c_str(),(Label(layer)+" NStrip="+to_string(nstrip)+" & NStripSat="+to_string(nstripsat)+" Sat"+SatType).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
                    if(sat254==1) VectSatTypeProfileFit.push_back(new TH2F((Label(layer)+" NStrip="+to_string(nstrip)+" & NStripSat="+to_string(nstripsat)+" Sat"+SatType).c_str(),(Label(layer)+" NStrip="+to_string(nstrip)+" & NStripSat="+to_string(nstripsat)+" Sat"+SatType).c_str(),300,0,6000*pow(10,-6),300,0,6000*pow(10,-6)));
                }
                VectNStripSatVectSatTypeProfileFit.push_back(VectSatTypeProfileFit);
            }
            VectNStripVectNStripSatVectSatTypeProfileFit.push_back(VectNStripSatVectSatTypeProfileFit);
        }
        VectLayerVectNStripVectNStripSatVectSatTypeProfileFit.push_back(VectNStripVectNStripSatVectSatTypeProfileFit);
    }
    vector<vector<vector<vector<TH2F*>>>> VectLayerVectNStripVectNStripSat254VectNStripSat255Histo;
    for(int layer=1;layer<11;layer++)
    {
        vector<vector<vector<TH2F*>>> VectNStripVectNStripSat254VectNStripSat255Histo;
        for(int nstrip=3;nstrip<7;nstrip++)
        {
            vector<vector<TH2F*>> VectNStripSat254VectNStripSat255Histo;
            for(int nstripsat254=0;nstripsat254<nstrip+1;nstripsat254++)
            {
                vector<TH2F*> VectNStripSat255Histo;
                for(int nstripsat255=0;nstripsat255<nstrip-nstripsat254+1;nstripsat255++)
                {
                    string title = Label(layer)+" NStrip="+to_string(nstrip)+" NStripSat254="+to_string(nstripsat254)+" NStripSat255="+to_string(nstripsat255);
                    VectNStripSat255Histo.push_back(new TH2F(title.c_str(),title.c_str(),300,0,6000*pow(10,-6),300,0,6000*pow(10,-6)));
                }
                VectNStripSat254VectNStripSat255Histo.push_back(VectNStripSat255Histo);
            }
            VectNStripVectNStripSat254VectNStripSat255Histo.push_back(VectNStripSat254VectNStripSat255Histo);
        }
        VectLayerVectNStripVectNStripSat254VectNStripSat255Histo.push_back(VectNStripVectNStripSat254VectNStripSat255Histo);
    }

    Builder* b1 = new Builder(chain);
    b1->SetBranchAdd();
    int nentries = b1->GetEntries();

    CalibCharge calibration;

    calibration.SetParameters();

    Calibration cal1;
    TFile* file1 = TFile::Open("test.root");
    TTree* tree1 = (TTree*) file1->Get("tree");
    cal1.SetTree(*tree1);

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

    TH2F* h2ElossvQ_tot = new TH2F("h2ElossvQ_tot","h2ElossvQ_tot",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
    TH2F* h2ElossvQ_NoSat = new TH2F("h2ElossvQ_NoSat","h2ElossvQ_NoSat",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
    TH2F* h2ElossvQ_Sat = new TH2F("h2ElossvQ_Sat","h2ElossvQ_Sat",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
    TH2F* h2ElossvQ_Sat254 = new TH2F("h2ElossvQ_Sat254","h2ElossvQ_Sat254",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
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

    TH2F* h2RatioSatPartID254 = new TH2F("h2RatioSatPartID254","h2RatioSatPartID254",7,0,7,20,0,1.05);
    TProfile* profSatPartID254 = new TProfile("profSatPartID254","profSatPartID254",7,0,7,"");

    TH2F* h2RatioSatPoverM254 = new TH2F("h2RatioSatPoverM254","h2RatioSatPoverM254",100,0,10,20,0,1.05);
    TProfile* profSatPoverM254 = new TProfile("profSatPoverM254","profSatPoverM254",100,0,10,"");

    TH2F* h2RatioSatPt255 = new TH2F("h2RatioSatPt255","h2RatioSatPt255",100,0,3000,20,0,1.05);
    TProfile* profSatPt255 = new TProfile("profSatPt255","profSatPt255",100,0,3000,"");

    TH2F* h2RatioSatEloss255 = new TH2F("h2RatioSatEloss255","h2RatioSatEloss255",300,0,900*pow(10,-6),20,0,1.05);
    TProfile* profSatEloss255 = new TProfile("profSatEloss255","profSatEloss255",300,0,900*pow(10,-6),"");

    TH2F* h2RatioSatPartID255 = new TH2F("h2RatioSatPartID255","h2RatioSatPartID255",7,0,7,20,0,1.05);
    TProfile* profSatPartID255 = new TProfile("profSatPartID255","profSatPartID255",7,0,7,"");

    TH2F* h2RatioSatPoverM255 = new TH2F("h2RatioSatPoverM255","h2RatioSatPoverM255",100,0,10,20,0,1.05);
    TProfile* profSatPoverM255 = new TProfile("profSatPoverM255","profSatPoverM255",100,0,10,"");

    vector<vector<TH2F*>> VectNStrip_VectNStripSat254_h2EvQ;
    vector<vector<TH1F*>> VectNStrip_VectNStripSat254_h1EvQ;
    vector<vector<TProfile*>> VectNStrip_VectNStripSat254_profEvQ;
    vector<vector<TH2F*>> VectNStrip_VectNStripSat255_h2EvQ;
    vector<vector<TH1F*>> VectNStrip_VectNStripSat255_h1EvQ;
    vector<vector<TProfile*>> VectNStrip_VectNStripSat255_profEvQ;
    for(int i=1;i<10;i++)
    {
        vector<TH2F*> VectNStripSat254_h2EvQ;
        vector<TH1F*> VectNStripSat254_h1EvQ;
        vector<TProfile*> VectNStripSat254_profEvQ;
        vector<TH2F*> VectNStripSat255_h2EvQ;
        vector<TH1F*> VectNStripSat255_h1EvQ;
        vector<TProfile*> VectNStripSat255_profEvQ;
        for(int j=0;j<i+1;j++)
        {
            string str_i = to_string(i);
            string str_j = to_string(j);
            TH2F* h2EvQVectVectSat254 = new TH2F(("H2 NStrip="+str_i+" & NStripSat254="+str_j+" | TOB1").c_str(),("H2 NStrip="+str_i+" & NStripSat254="+str_j+" | TOB1").c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
            TH1F* h1EvQVectVectSat254 = new TH1F(("H1 NStrip="+str_i+" & NStripSat254="+str_j+" | TOB1").c_str(),("H1 NStrip="+str_i+" & NStripSat254="+str_j+" | TOB1").c_str(),200,0,8);
            TProfile* profEvQVectVectSat254 = new TProfile(("prof NStrip="+str_i+" & NStripSat254="+str_j+" | TOB1").c_str(),("prof NStrip="+str_i+" & NStripSat254="+str_j+" | TOB1").c_str(),50,0,1500*pow(10,-6));
            VectNStripSat254_h2EvQ.push_back(h2EvQVectVectSat254);
            VectNStripSat254_h1EvQ.push_back(h1EvQVectVectSat254);
            VectNStripSat254_profEvQ.push_back(profEvQVectVectSat254);
            TH2F* h2EvQVectVectSat255 = new TH2F(("H2 NStrip="+str_i+" & NStripSat255="+str_j+" | TOB1").c_str(),("H2 NStrip="+str_i+" & NStripSat255="+str_j+" | TOB1").c_str(),300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
            TH1F* h1EvQVectVectSat255 = new TH1F(("H1 NStrip="+str_i+" & NStripSat255="+str_j+" | TOB1").c_str(),("H1 NStrip="+str_i+" & NStripSat255="+str_j+" | TOB1").c_str(),200,0,8);
            TProfile* profEvQVectVectSat255 = new TProfile(("prof NStrip="+str_i+" & NStripSat255="+str_j+" | TOB1").c_str(),("prof NStrip="+str_i+" & NStripSat255="+str_j+" | TOB1").c_str(),50,0,6000*pow(10,-6));
            VectNStripSat255_h2EvQ.push_back(h2EvQVectVectSat255);
            VectNStripSat255_h1EvQ.push_back(h1EvQVectVectSat255);
            VectNStripSat255_profEvQ.push_back(profEvQVectVectSat255);
        }
        VectNStrip_VectNStripSat254_h2EvQ.push_back(VectNStripSat254_h2EvQ);
        VectNStrip_VectNStripSat254_h1EvQ.push_back(VectNStripSat254_h1EvQ);
        VectNStrip_VectNStripSat254_profEvQ.push_back(VectNStripSat254_profEvQ);
        VectNStrip_VectNStripSat255_h2EvQ.push_back(VectNStripSat255_h2EvQ);
        VectNStrip_VectNStripSat255_h1EvQ.push_back(VectNStripSat255_h1EvQ);
        VectNStrip_VectNStripSat255_profEvQ.push_back(VectNStripSat255_profEvQ);
    }
    
    vector<TH2F*> VectPartID_h2EvQ_NoSat;
    vector<TH2F*> VectPartID_h2EvQ_Sat;
    vector<TH2F*> VectPartID_h2EvQ_Sat254;
    vector<TH2F*> VectPartID_h2EvQ_Sat255;
    vector<TH1F*> VectPartID_Pt;
    vector<TH1F*> VectPartID_Eloss;
    vector<TH1F*> VectPartID_PoverM;


    for(int i=0;i<7;i++)
    {
        VectPartID_h2EvQ_NoSat.push_back(new TH2F(("NoSat_"+LoopPartID(i)).c_str(),("NoSat_"+LoopPartID(i)).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
        VectPartID_h2EvQ_Sat.push_back(new TH2F(("Sat_"+LoopPartID(i)).c_str(),("Sat_"+LoopPartID(i)).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
        VectPartID_h2EvQ_Sat254.push_back(new TH2F(("Sat254_"+LoopPartID(i)).c_str(),("Sat254_"+LoopPartID(i)).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
        VectPartID_h2EvQ_Sat255.push_back(new TH2F(("Sat255_"+LoopPartID(i)).c_str(),("Sat255_"+LoopPartID(i)).c_str(),300,0,3000*pow(10,-6),300,0,3000*pow(10,-6)));
        VectPartID_Pt.push_back(new TH1F((LoopPartID(i)+"_p_{T}").c_str(),(LoopPartID(i)+"_p_{T}").c_str(),100,0,3000));
        VectPartID_Eloss.push_back(new TH1F((LoopPartID(i)+"_E_{loss}").c_str(),(LoopPartID(i)+"_E_{loss}").c_str(),300,0,3000*pow(10,-6)));
        VectPartID_PoverM.push_back(new TH1F((LoopPartID(i)+"_P/M").c_str(),(LoopPartID(i)+"_P/M").c_str(),1000,0.01,10));
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

    TH1F* hEmpty = new TH1F("hEmpty","hEmpty",21,0,21); //pour tracer les lignes avec la fonction SetHistoLabel

    TH2F* h2NStrip5Sat1_254_1=new TH2F("h2NStrip5Sat1_254_1","h2NStrip5Sat1_254_1",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
    TH2F* h2NStrip5Sat1_254_2=new TH2F("h2NStrip5Sat1_254_2","h2NStrip5Sat1_254_2",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
    TH2F* h2NStrip5Sat1_254_3=new TH2F("h2NStrip5Sat1_254_3","h2NStrip5Sat1_254_3",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
    TH2F* h2NStrip5Sat1_254_4=new TH2F("h2NStrip5Sat1_254_4","h2NStrip5Sat1_254_4",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
    TH2F* h2NStrip5Sat1_254_5=new TH2F("h2NStrip5Sat1_254_5","h2NStrip5Sat1_254_5",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));

    TProfile* profNstrip5Sat1_254_1=new TProfile("profNstrip5Sat1_254_1","profNstrip5Sat1_254_1",50,0,1500*pow(10,-6),"");
    TProfile* profNstrip5Sat1_254_2=new TProfile("profNstrip5Sat1_254_2","profNstrip5Sat1_254_2",50,0,1500*pow(10,-6),"");
    TProfile* profNstrip5Sat1_254_3=new TProfile("profNstrip5Sat1_254_3","profNstrip5Sat1_254_3",50,0,1500*pow(10,-6),"");
    TProfile* profNstrip5Sat1_254_4=new TProfile("profNstrip5Sat1_254_4","profNstrip5Sat1_254_4",50,0,1500*pow(10,-6),"");
    TProfile* profNstrip5Sat1_254_5=new TProfile("profNstrip5Sat1_254_5","profNstrip5Sat1_254_5",50,0,1500*pow(10,-6),"");

    TH2F* h2TestCutEdge = new TH2F("h2TestCutEdge","h2TestCutEdge",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
    TH1F* h1TestNoCutAndNoEdge = new TH1F("h1TestNoCutAndNoEdge","h1TestNoCutAndNoEdge",200,0,8);
    TH1F* h1TestCutOrEdge = new TH1F("h1TestCutOrEdge","h1TestCutOrEdge",200,0,8);

    TH1F* h1ChargeFauxSimHit = new TH1F("h1FauxSimHit","h1FauxSimHit",255,0,255);

    TH2F* testTOB1 = new TH2F("testTOB1","testTOB1",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
    TH2F* testTOB1_Calib = new TH2F("testTOB1_Calib","testTOB1_Calib",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));

    TH1F* h1DiffRelEvQ = new TH1F("h1DiffRelEvQ","h1DiffRelEvQ",1000,-1.5,5);
    TH1F* h1DiffRelEvQCalib = new TH1F("h1DiffRelEvQCalib","h1DiffRelEvQCalib",1000,-1.5,5);


    int entries = nentries;


    bool testsat254 = false; 
    bool testsat255 = false;

    int clust=0,clustsat254=0,clustsat255=0;

    b1->SetThresholdPartId(0.7); //on veut un minimum de 70% de simhits identiques
    b1->SetThresholdPt(60); //on veut des traces avec un minimum de 60 GeV en pt


    Calibration Charge;
    TFile* file2 = TFile::Open("FitRes.root");
    TTree* tree2 = (TTree*) file2->Get("tree");
    Charge.SetTree(*tree2);

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
            float RatioNClusterSat254 = (double)NClustSat254/(double)NCluster;
            float RatioNClusterSat255 = (double)NClustSat255/(double)NCluster;
            int id                  = b1->GetVectTrack()[track].GetPartId();
            float PoverM            = GetPoverM(p,id);

            if(id==11 || id==-11) 
            {
                VectPartID_Pt[0]->Fill(pt);
                VectPartID_PoverM[0]->Fill(GetPoverM(p,id));
            }
            if(id==211 || id==-211) 
            {
                VectPartID_Pt[1]->Fill(pt);
                VectPartID_PoverM[1]->Fill(GetPoverM(p,id));
            }
            if(id==321 || id==-321) 
            {
                VectPartID_Pt[2]->Fill(pt);
                VectPartID_PoverM[2]->Fill(GetPoverM(p,id));
            }
            if(id==2212 || id==-2212) 
            {
                VectPartID_Pt[3]->Fill(pt);
                VectPartID_PoverM[3]->Fill(GetPoverM(p,id));
            }
            if(id==1009213 || id==-1009213) 
            {
                VectPartID_Pt[5]->Fill(pt);
                VectPartID_PoverM[5]->Fill(GetPoverM(p,id));
            }
            if((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001))
            {
                VectPartID_Pt[6]->Fill(pt);
            }
            else VectPartID_Pt[4]->Fill(pt);
            
            
            hPt_tot->Fill(pt);
            hNClusterPerTrack->Fill(NCluster);
            hNClusterSat254PerTrack->Fill(NClustSat254);
            hNClusterSat255PerTrack->Fill(NClustSat255);
            hRatio_NClusterSat254->Fill(RatioNClusterSat254);
            hRatio_NClusterSat255->Fill(RatioNClusterSat255);

            h2RatioSatPt254->Fill(pt,RatioNClusterSat254);
            profSatPt254->Fill(pt,RatioNClusterSat254);
            h2RatioSatPartID254->Fill(ReBinPartID(id),RatioNClusterSat254);
            profSatPartID254->Fill(ReBinPartID(id),RatioNClusterSat254);
            h2RatioSatPoverM254->Fill(PoverM,RatioNClusterSat254);
            profSatPoverM254->Fill(PoverM,RatioNClusterSat254);

            h2RatioSatPt255->Fill(pt,RatioNClusterSat255);
            profSatPt255->Fill(pt,RatioNClusterSat255);
            h2RatioSatPartID255->Fill(ReBinPartID(id),RatioNClusterSat255);
            profSatPartID255->Fill(ReBinPartID(id),RatioNClusterSat255);
            h2RatioSatPoverM255->Fill(PoverM,RatioNClusterSat255);
            profSatPoverM255->Fill(PoverM,RatioNClusterSat255);
            
            for(int cluster=0;cluster<b1->GetVectTrack()[track].GetNCluster();cluster++)
            {
                float charge        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSclusCharge();
                float eloss         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetEloss();
                int layer           = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetLayer();
                int layerLabel      = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetLayerLabel();
                bool sat254         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSat254();
                bool sat255         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSat255();
                //int nsimhits        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSimHits();
                int nstrips         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNStrip();
                //int nsatboth        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStripBoth();
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

                if(id==11 || id==-11) 
                {
                    if(layerLabel==5) VectPartID_Eloss[0]->Fill(eloss);
                    if(sat254==false && sat255==false) VectPartID_h2EvQ_NoSat[0]->Fill(eloss,charge);
                    if(sat254==true) VectPartID_h2EvQ_Sat[0]->Fill(eloss,charge);
                    if(sat254==true && sat255==false) VectPartID_h2EvQ_Sat254[0]->Fill(eloss,charge);
                    if(sat255==true) VectPartID_h2EvQ_Sat255[0]->Fill(eloss,charge);
                }
                if(id==211 || id==-211) 
                {
                    if(layerLabel==5) VectPartID_Eloss[1]->Fill(eloss);
                    if(sat254==false && sat255==false) VectPartID_h2EvQ_NoSat[1]->Fill(eloss,charge);
                    if(sat254==true) VectPartID_h2EvQ_Sat[1]->Fill(eloss,charge);
                    if(sat254==true && sat255==false) VectPartID_h2EvQ_Sat254[1]->Fill(eloss,charge);
                    if(sat255==true) VectPartID_h2EvQ_Sat255[1]->Fill(eloss,charge);
                }
                if(id==321 || id==-321) 
                {
                    if(layerLabel==5) VectPartID_Eloss[2]->Fill(eloss);
                    if(sat254==false && sat255==false) VectPartID_h2EvQ_NoSat[2]->Fill(eloss,charge);
                    if(sat254==true) VectPartID_h2EvQ_Sat[2]->Fill(eloss,charge);
                    if(sat254==true && sat255==false) VectPartID_h2EvQ_Sat254[2]->Fill(eloss,charge);
                    if(sat255==true) VectPartID_h2EvQ_Sat255[2]->Fill(eloss,charge);
                }
                if(id==2212 || id==-2212) 
                {
                    if(layerLabel==5) VectPartID_Eloss[3]->Fill(eloss);
                    if(sat254==false && sat255==false) VectPartID_h2EvQ_NoSat[3]->Fill(eloss,charge);
                    if(sat254==true) VectPartID_h2EvQ_Sat[3]->Fill(eloss,charge);
                    if(sat254==true && sat255==false) VectPartID_h2EvQ_Sat254[3]->Fill(eloss,charge);
                    if(sat255==true) VectPartID_h2EvQ_Sat255[3]->Fill(eloss,charge);
                }
                if(id==1009213 || id==-1009213) 
                {
                    if(layerLabel==5) VectPartID_Eloss[5]->Fill(eloss);
                    if(sat254==false && sat255==false) VectPartID_h2EvQ_NoSat[5]->Fill(eloss,charge);
                    if(sat254==true) VectPartID_h2EvQ_Sat[5]->Fill(eloss,charge);
                    if(sat254==true && sat255==false) VectPartID_h2EvQ_Sat254[5]->Fill(eloss,charge);
                    if(sat255==true) VectPartID_h2EvQ_Sat255[5]->Fill(eloss,charge);
                }
                if((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001))
                {
                    if(layerLabel==5) VectPartID_Eloss[6]->Fill(eloss);
                    if(sat254==false && sat255==false) VectPartID_h2EvQ_NoSat[6]->Fill(eloss,charge);
                    if(sat254==true) VectPartID_h2EvQ_Sat[6]->Fill(eloss,charge);
                    if(sat254==true && sat255==false) VectPartID_h2EvQ_Sat254[6]->Fill(eloss,charge);
                    if(sat255==true) VectPartID_h2EvQ_Sat255[6]->Fill(eloss,charge);
                }
                else 
                {
                    if(layerLabel==5) VectPartID_Eloss[4]->Fill(eloss);
                    if(sat254==false && sat255==false) VectPartID_h2EvQ_NoSat[4]->Fill(eloss,charge);
                    if(sat254==true) VectPartID_h2EvQ_Sat[4]->Fill(eloss,charge);
                    if(sat254==true && sat255==false) VectPartID_h2EvQ_Sat254[4]->Fill(eloss,charge);
                    if(sat255==true) VectPartID_h2EvQ_Sat255[4]->Fill(eloss,charge);
                }
                if(b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectSimHits().size()==0) h1ChargeFauxSimHit->Fill(charge);

                /*if((layerLabel==1) && (!sat255 && sat254) )
                {
                    bool IsSat254=false; int NStripSat=nsat255;
                    if(sat254 && !sat255) {IsSat254=true; NStripSat=nsat254;}
                    float chargecalib = cal1.CalibCharge(cal1.GetGoodEntry(layerLabel,nstrips,NStripSat,IsSat254),charge);
                    //testTOB1_Calib->Fill(eloss,calibration.Charge(layerLabel,nstrips,NStripSat,IsSat254,charge));
                    testTOB1_Calib->Fill(eloss,chargecalib);
                    testTOB1->Fill(eloss,charge);
                    h1DiffRelEvQ->Fill((charge-eloss)/eloss);
                    //h1DiffRelEvQCalib->Fill((calibration.Charge(layerLabel,nstrips,NStripSat,IsSat254,charge)-eloss)/eloss);
                    h1DiffRelEvQCalib->Fill((chargecalib-eloss)/eloss);
                    //cout<<cal1.CalibCharge(cal1.GetGoodEntry(layerLabel,nstrips,NStripSat,!IsSat254),charge)<<endl;
                    
                }*/
                /*if(layerLabel<=10)
                {
                    float rec_charge = charge;
                    if(RatioNClusterSat255>0. && sat255 && nsat254==0) rec_charge=cal1.RecCharge(charge,layerLabel,nstrips,nsat255,0);
                    else if(RatioNClusterSat254>0. && sat254 && nsat255==0) rec_charge=cal1.RecCharge(charge,layerLabel,nstrips,nsat254,1);
                    testTOB1_Calib->Fill(eloss,rec_charge);
                    testTOB1->Fill(eloss,charge);
                    h1DiffRelEvQ->Fill((charge-eloss)/eloss); 
                    h1DiffRelEvQCalib->Fill((rec_charge-eloss)/eloss);
                }*/

                if(layerLabel<=10)
                {
                    float ChargeCalib=charge;
                    if(RatioNClusterSat254>0.5 || RatioNClusterSat255>0.1)ChargeCalib=Charge.ChargeCalib(charge,layerLabel,nstrips,nsat254,nsat255);
                    testTOB1_Calib->Fill(eloss,ChargeCalib);
                    testTOB1->Fill(eloss,charge);
                    h1DiffRelEvQ->Fill((charge-eloss)/eloss); 
                    h1DiffRelEvQCalib->Fill((ChargeCalib-eloss)/eloss);
                }

                   
                for(int nstrip=1;nstrip<10;nstrip++)
                {
                    for(int nstripsat=0;nstripsat<nstrip+1;nstripsat++)
                    {
                        if(layerLabel==1) //on se place dans TOB1 pour ne pas melanger differents effets + selection sans les bords ou les cuts
                        {
                            if(nstrip==nstrips && nstripsat==nsat254 && nsat255==0)
                            {
                                VectNStrip_VectNStripSat254_h2EvQ[nstrip-1][nstripsat]->Fill(eloss,charge);
                                VectNStrip_VectNStripSat254_h1EvQ[nstrip-1][nstripsat]->Fill(eloss/charge);
                                VectNStrip_VectNStripSat254_profEvQ[nstrip-1][nstripsat]->Fill(eloss,charge);
                                if(nstrip==5 && nstripsat==1)
                                {   
                                    if(0.0<PoverM<=0.2)
                                    {
                                        h2NStrip5Sat1_254_1->Fill(eloss,charge);
                                        profNstrip5Sat1_254_1->Fill(eloss,charge);
                                    }
                                    if(0.2<PoverM<=0.4) 
                                    {
                                        h2NStrip5Sat1_254_2->Fill(eloss,charge);
                                        profNstrip5Sat1_254_2->Fill(eloss,charge);
                                    }
                                    if(0.4<PoverM<=0.6) 
                                    {
                                        h2NStrip5Sat1_254_3->Fill(eloss,charge);
                                        profNstrip5Sat1_254_3->Fill(eloss,charge);
                                    }
                                    if(0.6<PoverM<=0.8) 
                                    {
                                        h2NStrip5Sat1_254_4->Fill(eloss,charge);
                                        profNstrip5Sat1_254_4->Fill(eloss,charge);
                                    }
                                    if(0.8<PoverM<=1.0) 
                                    {
                                        h2NStrip5Sat1_254_5->Fill(eloss,charge);
                                        profNstrip5Sat1_254_5->Fill(eloss,charge);
                                    }
                                }
                            }
                            if(nstrip==nstrips && nstripsat==nsat255 && nsat254==0)
                            {
                                VectNStrip_VectNStripSat255_h2EvQ[nstrip-1][nstripsat]->Fill(eloss,charge);
                                VectNStrip_VectNStripSat255_h1EvQ[nstrip-1][nstripsat]->Fill(eloss/charge);
                                VectNStrip_VectNStripSat255_profEvQ[nstrip-1][nstripsat]->Fill(eloss,charge);
                            }
                        }
                    }
                }
                for(int counter_layer=1;counter_layer<11;counter_layer++)
                {
                    for(int nstrip=3;nstrip<7;nstrip++)
                    {
                        for(int nstripsat=1;nstripsat<3;nstripsat++)
                        {
                            if(counter_layer==layerLabel && nstrip==nstrips && nstripsat==nsat254 && nsat255==0)
                            {
                                VectLayerVectNStripVectNStripSatVectSatTypeProfileFit[counter_layer-1][nstrip-3][nstripsat-1][0]->Fill(eloss,charge);
                            }
                            if(counter_layer==layerLabel && nstrip==nstrips && nstripsat==nsat255 && nsat254==0)
                            {
                                VectLayerVectNStripVectNStripSatVectSatTypeProfileFit[counter_layer-1][nstrip-3][nstripsat-1][1]->Fill(eloss,charge);
                            }
                        }
                    }
                }
                for(int countlayer=1;countlayer<11;countlayer++)
                {
                    for(int countnstrip=3;countnstrip<7;countnstrip++)
                    {
                        for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
                        {
                            for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
                            {
                                if(countlayer==layerLabel && countnstrip==nstrips && countnstripsat254==nsat254 && countnstripsat255==nsat255) VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(eloss,charge);
                            }
                        }
                    }
                }
                if(cut || edge) 
                {
                    h2TestCutEdge->Fill(eloss,charge);
                    h1TestCutOrEdge->Fill(eloss/charge);
                }
                if(cut==false && edge==false)
                {
                    h1TestNoCutAndNoEdge->Fill(eloss/charge);
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
                }
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

            h2RatioSatEloss254->Fill(estim.GetMean(),RatioNClusterSat254);
            profSatEloss254->Fill(estim.GetMean(),RatioNClusterSat254);
            h2RatioSatEloss255->Fill(estim.GetMean(),RatioNClusterSat255);
            profSatEloss255->Fill(estim.GetMean(),RatioNClusterSat255);
            
            

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
    
    TFile* fileout = new TFile(s2.c_str(),"RECREATE");

    ofstream ofile_fit("./data/fit_res.txt",ios::out | ios::trunc);

    /*for(int layer=1;layer<11;layer++)
    {
        for(int counter_sat254=0;counter_sat254<2;counter_sat254++)
        {
        for(int nstrip=3;nstrip<7;nstrip++)
        {
            for(int nstripsat=1;nstripsat<3;nstripsat++)
            {
                
                    string SatType="254";
                    if(counter_sat254==1) SatType="255";
                    if(VectLayerVectNStripVectNStripSatVectSatTypeProfileFit[layer-1][nstrip-3][nstripsat-1][counter_sat254]->GetEntries()>10)
                    {
                        TFitResultPtr fit_res = VectLayerVectNStripVectNStripSatVectSatTypeProfileFit[layer-1][nstrip-3][nstripsat-1][counter_sat254]->Fit("pol1","S");
                        ofile_fit<<(Label(layer)+"_NStrip="+to_string(nstrip)+"_NStripSat="+to_string(nstripsat)+"_SatType="+SatType).c_str()<<"\t"<<fit_res->Parameter(0)<<"\t"<<fit_res->Parameter(1)<<endl;
                    }
                    else if(!(nstrip==3 && nstripsat==2))ofile_fit<<(Label(layer)+"_NStrip="+to_string(nstrip)+"_NStripSat="+to_string(nstripsat)+"_SatType="+SatType).c_str()<<"\t"<<0<<"\t"<<1<<endl;
                
            }
        }
        }
    }*/
    
    ofile_fit.close();

    testTOB1->Write();
    testTOB1_Calib->Write();

    TCanvas* cDiffRelChargeCalib_Charge = new TCanvas("cDiffRelChargeCalib_Charge","cDiffRelChargeCalib_Charge");
    vector<TH1F*> VectDiffRelChargeCalib_Charge;
    VectDiffRelChargeCalib_Charge.push_back(h1DiffRelEvQ);
    VectDiffRelChargeCalib_Charge.push_back(h1DiffRelEvQCalib);
    vector<string> VectLegendDiffRelChargeCalib_Charge;
    VectLegendDiffRelChargeCalib_Charge.push_back("Avant calibration");
    VectLegendDiffRelChargeCalib_Charge.push_back("Apres calibration");
    DrawHistoNormalized(*cDiffRelChargeCalib_Charge,VectDiffRelChargeCalib_Charge,VectLegendDiffRelChargeCalib_Charge,"Calibration","#frac{Q-E_{loss}}{E_{loss}}");
    cDiffRelChargeCalib_Charge->Write();


    TLine* linexy = new TLine(0,0,0.0015,0.0015);

    TCanvas* cTestCutEdge = new TCanvas("cTestCutEdge","cTestCutEdge");
    vector<TH1F*> VectTestCutEdge;
    VectTestCutEdge.push_back(h1TestNoCutAndNoEdge);
    VectTestCutEdge.push_back(h1TestCutOrEdge);
    vector<string> VectLegendTestCutEdge;
    VectLegendTestCutEdge.push_back("No cut & no edge");
    VectLegendTestCutEdge.push_back("Cut or edge");
    DrawHistoNormalized(*cTestCutEdge,VectTestCutEdge,VectLegendTestCutEdge,"Test cut and edge","#frac{E_{loss}}{Q}");
    cTestCutEdge->Write();

    
    TCanvas* cNStrip5Sat1_254_1=new TCanvas("cNStrip5Sat1_254_1","cNStrip5Sat1_254_1");
    SuperposedHisto2DProfile(*cNStrip5Sat1_254_1,h2NStrip5Sat1_254_1,profNstrip5Sat1_254_1,"NStrip=5 & NStripSat=1 | TOB1 | Selected Area | 0.0<P/M<0.2","E_{loss}","Q");
    cNStrip5Sat1_254_1->Write();
    TCanvas* cNStrip5Sat1_254_2=new TCanvas("cNStrip5Sat1_254_2","cNStrip5Sat1_254_2");
    SuperposedHisto2DProfile(*cNStrip5Sat1_254_2,h2NStrip5Sat1_254_2,profNstrip5Sat1_254_2,"NStrip=5 & NStripSat=1 | TOB1 | Selected Area | 0.2<P/M<0.4","E_{loss}","Q");
    cNStrip5Sat1_254_2->Write();
    TCanvas* cNStrip5Sat1_254_3=new TCanvas("cNStrip5Sat1_254_3","cNStrip5Sat1_254_3");
    SuperposedHisto2DProfile(*cNStrip5Sat1_254_3,h2NStrip5Sat1_254_3,profNstrip5Sat1_254_3,"NStrip=5 & NStripSat=1 | TOB1 | Selected Area | 0.4<P/M<0.6","E_{loss}","Q");
    cNStrip5Sat1_254_3->Write();
    TCanvas* cNStrip5Sat1_254_4=new TCanvas("cNStrip5Sat1_254_4","cNStrip5Sat1_254_4");
    SuperposedHisto2DProfile(*cNStrip5Sat1_254_4,h2NStrip5Sat1_254_4,profNstrip5Sat1_254_4,"NStrip=5 & NStripSat=1 | TOB1 | Selected Area | 0.6<P/M<0.8","E_{loss}","Q");
    cNStrip5Sat1_254_4->Write();
    TCanvas* cNStrip5Sat1_254_5=new TCanvas("cNStrip5Sat1_254_5","cNStrip5Sat1_254_5");
    SuperposedHisto2DProfile(*cNStrip5Sat1_254_5,h2NStrip5Sat1_254_5,profNstrip5Sat1_254_5,"NStrip=5 & NStripSat=1 | TOB1 | Selected Area | 0.8<P/M<1.0","E_{loss}","Q");
    cNStrip5Sat1_254_5->Write();

    DrawHisto(*fileout,h1ChargeFauxSimHit,"Distribution de charge des faux SimHits","Q");
    DrawHisto(*fileout,h2TestCutEdge,"Cut or Edge","E_{loss}","Q");
    DrawHisto(*fileout,h1TestCutOrEdge,"Cut or edge","E/Q");

    vector<string> VectLegendPartId;
    VectLegendPartId.push_back("e^{+/-}");
    VectLegendPartId.push_back("#pi^{+/-}");
    VectLegendPartId.push_back("K^{+/-}");
    VectLegendPartId.push_back("p/#bar{p}");
    VectLegendPartId.push_back("SM");
    VectLegendPartId.push_back("R^{+/-}_{#scale[0.6]{#tilde{g}u#bar{d}}}");
    VectLegendPartId.push_back("SUSY");

    TCanvas* cPtPartID = new TCanvas("cPtPartID","cPtPartID");
    StackHisto(*cPtPartID,VectPartID_Pt,VectLegendPartId,"Distribution en p_{T} pour differentes particules","p_{T}");
    cPtPartID->Write();

    TCanvas* cElossPartID = new TCanvas("cElossPartID","cElossPartID");
    StackHisto(*cElossPartID,VectPartID_Eloss,VectLegendPartId,"Distribution E_{loss} pour differentes particules dans TOB1","E_{loss}");
    cElossPartID->Write();

    TCanvas* cPoverMPartID = new TCanvas("cPoverMPartID","cPoverMPartID");
    StackHisto(*cPoverMPartID,VectPartID_PoverM,VectLegendPartId,"Distribution du P/M pour differentes particules","P/M");
    cPoverMPartID->Write();

    TCanvas* cPoverMPartIDNormal = new TCanvas("cPoverMPartIDNormal","cPoverMPartIDNormal");
    vector<TH1F*> VectPartID_PoverM_Normalized;
    for(int vect=0;vect<7;vect++)
    {
        if(vect!=0 && vect!=5 && vect!=7) VectPartID_PoverM_Normalized.push_back(VectPartID_PoverM[vect]);
    }
    vector<string> VectLegendPartIdNorm;
    VectLegendPartIdNorm.push_back("#pi^{+/-}");
    VectLegendPartIdNorm.push_back("K^{+/-}");
    VectLegendPartIdNorm.push_back("p/#bar{p}");
    VectLegendPartIdNorm.push_back("R^{+/-}_{#scale[0.6]{#tilde{g}u#bar{d}}}");
    DrawHistoNormalized(*cPoverMPartIDNormal,VectPartID_PoverM_Normalized,VectLegendPartIdNorm,"Distribution du P/M pour differentes particules","P/M");
    cPoverMPartIDNormal->Write();

    for(int i=0;i<7;i++)
    {
        DrawHisto(*fileout,VectPartID_h2EvQ_NoSat[i],("NoSat_"+LoopPartID(i)).c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectPartID_h2EvQ_Sat[i],("Sat_"+LoopPartID(i)).c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectPartID_h2EvQ_Sat254[i],("Sat254_"+LoopPartID(i)).c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectPartID_h2EvQ_Sat255[i],("Sat255_"+LoopPartID(i)).c_str(),"E_{loss}","Q");
        DrawHisto(*fileout,VectPartID_Pt[i],(LoopPartID(i)+"_p_{T}").c_str(),"p_{T}");
        DrawHisto(*fileout,VectPartID_Eloss[i],(LoopPartID(i)+"_E_{loss}").c_str(),"E_{loss}");
        DrawHisto(*fileout,VectPartID_PoverM[i],(LoopPartID(i)+"_P/M").c_str(),"P/M");
    }
    
    TCanvas* cRatioSatPt254 = new TCanvas("cRatioSatPt254","cRatioSatPt254");
    SuperposedHisto2DProfile(*cRatioSatPt254,h2RatioSatPt254,profSatPt254,"Ratio de clusters qui saturent a 254 en fonction du p_{T}","p_{T}","#frac{#ClusterSat}{#Cluster}");
    cRatioSatPt254->Write();

    TCanvas* cRatioSatEloss254 = new TCanvas("cRatioSatEloss254","cRatioSatEloss254");
    SuperposedHisto2DProfile(*cRatioSatEloss254,h2RatioSatEloss254,profSatEloss254,"Ratio de clusters qui saturent a 254 en fonction de E_{loss}","E_{loss}","#frac{#ClusterSat}{#Cluster}");
    cRatioSatEloss254->Write();

    TCanvas* cRatioSatPartID254 = new TCanvas("cRatioSatPartID254","cRatioSatPartID254");
    SetHistoLabelPartID(cRatioSatPartID254,h2RatioSatPartID254);
    profSatPartID254->SetMarkerStyle(2);
    profSatPartID254->SetMarkerColor(kRed);
    profSatPartID254->SetLineColor(kRed);
    profSatPartID254->Draw("SAME");
    h2RatioSatPartID254->GetXaxis()->SetTitle("");
    h2RatioSatPartID254->SetDrawOption("COLZ");
    h2RatioSatPartID254->SetMaximum(800);
    cRatioSatPartID254->Write();

    TCanvas* cRatioSatPt255 = new TCanvas("cRatioSatPt255","cRatioSatPt255");
    SuperposedHisto2DProfile(*cRatioSatPt255,h2RatioSatPt255,profSatPt255,"Ratio de clusters qui saturent a 255 en fonction du p_{T}","p_{T}","#frac{#ClusterSat}{#Cluster}");
    cRatioSatPt255->Write();

    TCanvas* cRatioSatEloss255 = new TCanvas("cRatioSatEloss255","cRatioSatEloss255");
    SuperposedHisto2DProfile(*cRatioSatEloss255,h2RatioSatEloss255,profSatEloss255,"Ratio de clusters qui saturent a 255 en fonction de E_{loss}","E_{loss}","#frac{#ClusterSat}{#Cluster}");
    cRatioSatEloss255->Write();

    TCanvas* cRatioSatPartID255 = new TCanvas("cRatioSatPartID255","cRatioSatPartID255");
    SetHistoLabelPartID(cRatioSatPartID255,h2RatioSatPartID255);
    profSatPartID255->SetMarkerStyle(2);
    profSatPartID255->SetMarkerColor(kRed);
    profSatPartID255->SetLineColor(kRed);
    profSatPartID255->Draw("SAME");
    h2RatioSatPartID255->GetXaxis()->SetTitle("");
    h2RatioSatPartID255->SetDrawOption("COLZ");
    h2RatioSatPartID255->SetMaximum(800);
    cRatioSatPartID255->Write();

    TCanvas* cRatioSatPoverM_254 = new TCanvas("cRatioSatPoverM_254","cRatioSatPoverM_254");
    SuperposedHisto2DProfile(*cRatioSatPoverM_254,h2RatioSatPoverM254,profSatPoverM254,"Frac Sat254 = f(P/M)","P/M","#frac{#ClustSat}{#Clust}");
    cRatioSatPoverM_254->Write();

    TCanvas* cRatioSatPoverM_255 = new TCanvas("cRatioSatPoverM_255","cRatioSatPoverM_255");
    SuperposedHisto2DProfile(*cRatioSatPoverM_255,h2RatioSatPoverM255,profSatPoverM255,"Frac Sat255 = f(P/M)","P/M","#frac{#ClustSat}{#Clust}");
    cRatioSatPoverM_255->Write();

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
    TLine* line254 = new TLine(0,rat254,21,rat254);
    TLine* line255 = new TLine(0,rat255,21,rat255);

    TCanvas* cratioSatLayer = new TCanvas("cratioSatLayer","cratioSatLayer");
    SetHistoLabel(cratioSatLayer,hEmpty);
    cratioSatLayer->cd();
    //hEmpty->GetYaxis()->SetRangeUser(pow(10,-4),0.25);
    Ratio_LayerSat254->Draw("SAME");
    Ratio_LayerSat255->SetLineColor(2);
    Ratio_LayerSat255->Draw("SAME");
    line254->SetLineStyle(3);
    line254->Draw("SAME");
    line255->SetLineStyle(3);
    line255->SetLineColor(2);
    line255->Draw("SAME");
    //cratioSatLayer->SetLogy();
    TLegend* legratio = new TLegend(0.6,0.8,0.9,0.9);
    legratio->AddEntry(line254,("Saturation 254 : "+to_string(rat254)).c_str(),"l");
    legratio->AddEntry(line255,("Saturation 255 : "+to_string(rat255)).c_str(),"l");
    legratio->Draw("SAME");
    cratioSatLayer->Write();

// ------------------------------------------

// ------------------------------------------ TH2F + TProfile 

    TCanvas* canvas_EvQ_h2prof_tot = new TCanvas("canvas_EvQ_h2prof_tot","canvas_EvQ_h2prof_tot");
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_tot,h2ElossvQ_tot,profElossVsQ_tot,"","E_{loss}","Q");
    canvas_EvQ_h2prof_tot->Write();

    TCanvas* canvas_EvQ_h2prof_NoSat = new TCanvas("canvas_EvQ_h2prof_NoSat","canvas_EvQ_h2prof_NoSat");
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_NoSat,h2ElossvQ_NoSat,profElossVsQ_NoSat,"","E_{loss}","Q");
    canvas_EvQ_h2prof_NoSat->Write();

    TCanvas* canvas_EvQ_h2prof_Sat = new TCanvas("canvas_EvQ_h2prof_Sat","canvas_EvQ_h2prof_Sat");
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_Sat,h2ElossvQ_Sat,profElossVsQ_Sat,"","E_{loss}","Q");
    canvas_EvQ_h2prof_Sat->Write();

    TCanvas* canvas_EvQ_h2prof_Sat254 = new TCanvas("canvas_EvQ_h2prof_Sat254","canvas_EvQ_h2prof_Sat254");
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_Sat254,h2ElossvQ_Sat254,profElossVsQ_Sat254,"","E_{loss}","Q");
    canvas_EvQ_h2prof_Sat254->Write();

    TCanvas* canvas_EvQ_h2prof_Sat255 = new TCanvas("canvas_EvQ_h2prof_Sat255","canvas_EvQ_h2prof_Sat255");
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_Sat255,h2ElossvQ_Sat255,profElossVsQ_Sat255,"","E_{loss}","Q");
    canvas_EvQ_h2prof_Sat255->Write();

    profElossVsQ_tot->SetMarkerStyle(1);
    profElossVsQ_NoSat->SetMarkerStyle(1);
    profElossVsQ_Sat->SetMarkerStyle(1);
    profElossVsQ_Sat254->SetMarkerStyle(1);
    profElossVsQ_Sat255->SetMarkerStyle(1);

// ------------------------------------------

// ------------------------------------------ TH1F Pt de la trace, avec ou sans saturation

    TCanvas* canvas_Pt = new TCanvas("canvas_Pt","canvas_Pt");
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

    TCanvas* canvas_RatioNClusterSat = new TCanvas("canvas_RatioNClusterSat","canvas_RatioNClusterSat");
    vector<TH1F*> VectHistoNSat;
    VectHistoNSat.push_back(hRatio_ElossQ_Sat254);
    VectHistoNSat.push_back(hRatio_ElossQ_Sat255);
    vector<string> VectLegendNSat;
    VectLegendNSat.push_back("Saturation 254");
    VectLegendNSat.push_back("Saturation 255");
    DrawHistoNormalized(*canvas_RatioNClusterSat,VectHistoNSat,VectLegendNSat,"Fraction de saturation","#frac{Nombre de cluster sature}{Nombre de cluster dans la trace}");
    canvas_RatioNClusterSat->Write();

// ------------------------------------------

// ------------------------------------------ Ratio E/Q

    TCanvas* canvas_EvQ_Ratio = new TCanvas("canvas_EvQ_Ratio","canvas_EvQ_Ratio");
    vector<TH1F*> VectHistoEoverQ;
    VectHistoEoverQ.push_back(hRatio_ElossQ_NoSat);
    VectHistoEoverQ.push_back(hRatio_ElossQ_Sat);
    vector<string> VectLegendEoverQ;
    VectLegendEoverQ.push_back("No saturation");
    VectLegendEoverQ.push_back("Saturation");
    DrawHistoNormalized(*canvas_EvQ_Ratio,VectHistoEoverQ,VectLegendEoverQ,"","#frac{E}{Q}");
    canvas_EvQ_Ratio->Write();

// ------------------------------------------
/*
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
*/

/*    for(int i=0;i<21;i++)
    {
        string title = "E_{loss} v. Q  |  "+Label(i+1);
        TCanvas* csuperposed = new TCanvas(title.c_str(),title.c_str());
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
    }*/

    

    
    for(int nstrip=1;nstrip<10;nstrip++)
    {
        TCanvas* cNStripNStripSat254_h2prof = new TCanvas(("H2 & profile NStrip="+to_string(nstrip)+" | Sat254 | TOB1").c_str(),("H2 & profile NStrip="+to_string(nstrip)+" | Sat254 | TOB1").c_str());
        TCanvas* cNStripNStripSat254_h1 = new TCanvas(("H1 NStrip="+to_string(nstrip)+" | Sat254 | TOB1").c_str(),("H1 NStrip="+to_string(nstrip)+" | Sat254 | TOB1").c_str());
        TLine* linesat254 = new TLine(0,254*(3.61*pow(10,-9)*247),0.0015,254*(3.61*pow(10,-9)*247));
        TLine* linesat254times2 = new TLine(0,2*254*(3.61*pow(10,-9)*247),0.0015,2*254*(3.61*pow(10,-9)*247));
        TLine* linesat254times3 = new TLine(0,3*254*(3.61*pow(10,-9)*247),0.0015,3*254*(3.61*pow(10,-9)*247));
        TLine* linesat254times4 = new TLine(0,4*254*(3.61*pow(10,-9)*247),0.0015,4*254*(3.61*pow(10,-9)*247));
        linesat254->SetLineStyle(4);
        linesat254times2->SetLineStyle(4);
        linesat254times3->SetLineStyle(4);
        linesat254times4->SetLineStyle(4);
        linexy->SetLineStyle(3);
        if(nstrip==1)
        {
            cNStripNStripSat254_h2prof->Divide(2,1);
            cNStripNStripSat254_h1->Divide(2,1);
        }
        if(nstrip==2 || nstrip==3)
        {
            cNStripNStripSat254_h2prof->Divide(2,2);
            cNStripNStripSat254_h1->Divide(2,2);
        }
        if(nstrip==4 || nstrip==5)
        {
            cNStripNStripSat254_h2prof->Divide(2,3);
            cNStripNStripSat254_h1->Divide(2,3);
        }
        if(nstrip==6 || nstrip==7)
        {
            cNStripNStripSat254_h2prof->Divide(2,4);
            cNStripNStripSat254_h1->Divide(2,4);
        }
        if(nstrip==8 || nstrip==9)
        {
            cNStripNStripSat254_h2prof->Divide(2,5);
            cNStripNStripSat254_h1->Divide(2,5);
        }
        if(nstrip>=10)
        {
            cNStripNStripSat254_h2prof->Divide(2,6);
            cNStripNStripSat254_h1->Divide(2,6);
        }
        for(int nstripsat=0;nstripsat<nstrip+1;nstripsat++)
        {
            cNStripNStripSat254_h2prof->cd(nstripsat+1);
            SuperposedHisto2DProfile(*cNStripNStripSat254_h2prof,VectNStrip_VectNStripSat254_h2EvQ[nstrip-1][nstripsat],VectNStrip_VectNStripSat254_profEvQ[nstrip-1][nstripsat],("H2 & prof NStrip="+to_string(nstrip)+" & NStripSat254="+to_string(nstripsat)+" | TOB1").c_str(),"E_{loss}","Q");
            //VectNStrip_VectNStripSat254_h2EvQ[nstrip-1][nstripsat]->Draw("COLZ");
            linexy->Draw("SAME");
            if(nstripsat>=0)linesat254->Draw("SAME");
            if(nstripsat>=2)linesat254times2->Draw("SAME");
            if(nstripsat>=3)linesat254times3->Draw("SAME");
            if(nstripsat>=4)linesat254times4->Draw("SAME");
            cNStripNStripSat254_h1->cd(nstripsat+1);
            VectNStrip_VectNStripSat254_h1EvQ[nstrip-1][nstripsat]->Draw();
            DrawHisto(*fileout,VectNStrip_VectNStripSat254_h2EvQ[nstrip-1][nstripsat],("H2 NStrip="+to_string(nstrip)+" & NStripSat254="+to_string(nstripsat)).c_str(),"E_{loss}","Q");
            //DrawHisto(*fileout,VectNStrip_VectNStripSat254_h1EvQ[nstrip-1][nstripsat],("H1 NStrip="+to_string(nstrip)+" & NStripSat="+to_string(nstripsat)).c_str(),"E_{loss}");

        }
        cNStripNStripSat254_h2prof->Write();
        cNStripNStripSat254_h1->Write();
    }

    for(int nstrip=1;nstrip<10;nstrip++)
    {
        TCanvas* cNStripNStripSat255_h2prof = new TCanvas(("H2 & profile NStrip="+to_string(nstrip)+" | Sat255 | TOB1").c_str(),("H2 & profile NStrip="+to_string(nstrip)+" | Sat255 | TOB1").c_str());
        TCanvas* cNStripNStripSat255_h1 = new TCanvas(("H1 NStrip="+to_string(nstrip)+" | Sat255 | TOB1").c_str(),("H1 NStrip="+to_string(nstrip)+" | Sat255 | TOB1").c_str());
        TLine* linesat255 = new TLine(0,255*(3.61*pow(10,-9)*247),0.0015,255*(3.61*pow(10,-9)*247));
        TLine* linesat255times2 = new TLine(0,2*255*(3.61*pow(10,-9)*247),0.0015,2*255*(3.61*pow(10,-9)*247));
        TLine* linesat255times3 = new TLine(0,3*255*(3.61*pow(10,-9)*247),0.0015,3*255*(3.61*pow(10,-9)*247));
        TLine* linesat255times4 = new TLine(0,4*255*(3.61*pow(10,-9)*247),0.0015,4*255*(3.61*pow(10,-9)*247));
        linesat255->SetLineStyle(4);
        linesat255times2->SetLineStyle(4);
        linesat255times3->SetLineStyle(4);
        linesat255times4->SetLineStyle(4);
        linexy->SetLineStyle(3);
        if(nstrip==1)
        {
            cNStripNStripSat255_h2prof->Divide(2,1);
            cNStripNStripSat255_h1->Divide(2,1);
        }
        if(nstrip==2 || nstrip==3)
        {
            cNStripNStripSat255_h2prof->Divide(2,2);
            cNStripNStripSat255_h1->Divide(2,2);
        }
        if(nstrip==4 || nstrip==5)
        {
            cNStripNStripSat255_h2prof->Divide(2,3);
            cNStripNStripSat255_h1->Divide(2,3);
        }
        if(nstrip==6 || nstrip==7)
        {
            cNStripNStripSat255_h2prof->Divide(2,4);
            cNStripNStripSat255_h1->Divide(2,4);
        }
        if(nstrip==8 || nstrip==9)
        {
            cNStripNStripSat255_h2prof->Divide(2,5);
            cNStripNStripSat255_h1->Divide(2,5);
        }
        if(nstrip>=10)
        {
            cNStripNStripSat255_h2prof->Divide(2,6);
            cNStripNStripSat255_h1->Divide(2,6);
        }
        for(int nstripsat=0;nstripsat<nstrip+1;nstripsat++)
        {
            cNStripNStripSat255_h2prof->cd(nstripsat+1);
            SuperposedHisto2DProfile(*cNStripNStripSat255_h2prof,VectNStrip_VectNStripSat255_h2EvQ[nstrip-1][nstripsat],VectNStrip_VectNStripSat255_profEvQ[nstrip-1][nstripsat],("H2 & prof NStrip="+to_string(nstrip)+" & NStripSat255="+to_string(nstripsat)+" | TOB1").c_str(),"E_{loss}","Q");
            //VectNStrip_VectNStripSat255_h2EvQ[nstrip-1][nstripsat]->Draw("COLZ");
            linexy->Draw("SAME");
            if(nstripsat>=0) linesat255->Draw("SAME");
            if(nstripsat>=2) linesat255times2->Draw("SAME");
            if(nstripsat>=3) linesat255times3->Draw("SAME");
            if(nstripsat>=4) linesat255times4->Draw("SAME");
            cNStripNStripSat255_h1->cd(nstripsat+1);
            VectNStrip_VectNStripSat255_h1EvQ[nstrip-1][nstripsat]->Draw();
            DrawHisto(*fileout,VectNStrip_VectNStripSat255_h2EvQ[nstrip-1][nstripsat],("H2 NStrip="+to_string(nstrip)+" & NStripSat255="+to_string(nstripsat)).c_str(),"E_{loss}","Q");
            //DrawHisto(*fileout,VectNStrip_VectNStripSat255_h1EvQ[nstrip-1][nstripsat],("H1 NStrip="+to_string(nstrip)+" & NStripSat="+to_string(nstripsat)).c_str(),"E_{loss}");

        }
        cNStripNStripSat255_h2prof->Write();
        cNStripNStripSat255_h1->Write();
    }

    VectLayerVectNStripVectNStripSatVectSatTypeProfileFit[0][0][0][0]->Write();

    fileout->Write();
    fileout->Close();

    delete fileout;

    Calibration calib;
    calib.SetFileAndTreeName("test.root","tree");
    /*calib.SetHisto(*VectLayerVectNStripVectNStripSatVectSatTypeProfileFit[0][1][0][0]);
    calib.FillHisto(0.1);
    calib.FillProfile();
    calib.FitProfile();
    */

    /*calib.SetBranch();

    for(int layer=1;layer<11;layer++)
    {
        for(int counter_sat254=0;counter_sat254<2;counter_sat254++)
        {
            for(int nstrip=3;nstrip<7;nstrip++)
            {
                for(int nstripsat=1;nstripsat<3;nstripsat++)
                {
                    //cout<<"layer "<<layer<<" countersat "<<counter_sat254<<" nstrip "<<nstrip<<" nstripsat "<<nstripsat<<endl;
                    bool issat255=true;
                    if(counter_sat254!=0) issat255=false; 
                    calib.SetHisto(*VectLayerVectNStripVectNStripSatVectSatTypeProfileFit[layer-1][nstrip-3][nstripsat-1][counter_sat254]);
                    calib.FillHisto(0.1);
                    calib.FillProfile();
                    calib.Write(layer,nstrip,nstripsat,issat255);
                }
            }
        }
    }
    
    calib.WriteFile();*/


    /*Calibration calCharge;
    calCharge.SetFileAndTreeName("FitRes.root","tree");
    calCharge.SetBranch();
    for(int countlayer=1;countlayer<11;countlayer++)
    {
        for(int countnstrip=3;countnstrip<7;countnstrip++)
        {
            for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
            {
                for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
                {
                    calCharge.SetHisto(*VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]);
                    calCharge.FillHisto(0);
                    calCharge.FillProfile();
                    calCharge.Write(countlayer,countnstrip,countnstripsat254,countnstripsat255);
                }
            }
        }
    }
    calCharge.WriteFile();*/


    system(("mv "+s2+" ./data/.").c_str());
    

    return 0;
}
