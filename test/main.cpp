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
    TH1F* hLayerLabel = new TH1F("hLayerLabel","hLayerLabel",25,0,25);
    TH1F* hLayerLabelSat254 = new TH1F("hLayerLabelSat254","hLayerLabelSat254",25,0,25);
    TH1F* hLayerLabelSat255 = new TH1F("hLayerLabelSat255","hLayerLabelSat255",25,0,25);

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

    TH2F* h2RatioSatPt = new TH2F("h2RatioSatPt","h2RatioSatPt",100,0,3000,20,0,1.05);
    TProfile* profSatPt = new TProfile("profSatPt","profSatPt",100,0,3000,"");

    TH2F* h2RatioSatEloss = new TH2F("h2RatioSatEloss","h2RatioSatEloss",300,0,900*pow(10,-6),20,0,1.05);
    TProfile* profSatEloss = new TProfile("profSatEloss","profSatEloss",300,0,900*pow(10,-6),"");

    TH2F* h2RatioSatPartID = new TH2F("h2RatioSatPartID","h2RatioSatPartID",6,0,6,20,0,1.05);
    TProfile* profSatPartID = new TProfile("profSatPartID","profSatPartID",6,0,6,"");

    vector<TH2F*> vect_H2_tot;
    vector<TH2F*> vect_H2_NoSat;
    vector<TH2F*> vect_H2_Sat;
    vector<TH2F*> vect_H2_Sat254;
    vector<TH2F*> vect_H2_Sat255;

    vector<TProfile*> vect_prof_tot;
    vector<TProfile*> vect_prof_NoSat;
    vector<TProfile*> vect_prof_Sat;
    vector<TProfile*> vect_prof_Sat254;
    vector<TProfile*> vect_prof_Sat255;

    vector<TH1F*> vect_H1_tot;
    vector<TH1F*> vect_H1_NoSat;
    vector<TH1F*> vect_H1_Sat;
    vector<TH1F*> vect_H1_Sat254;
    vector<TH1F*> vect_H1_Sat255;

    for(int i=1;i<22;i++)
    {
        string title = "h2_"+Label(i)+"_tot";
        char title_char[title.length()+1];
        strcpy(title_char,title.c_str());
        vect_H2_tot.push_back(new TH2F(title_char,title_char,300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        title = "h1_"+Label(i)+"_tot";
        strcpy(title_char,title.c_str());
        vect_H1_tot.push_back(new TH1F(title_char,title_char,200,0,8));
        title = "pr_"+Label(i)+"_tot";
        strcpy(title_char,title.c_str());
        vect_prof_tot.push_back(new TProfile(title_char,title_char,50,0,900*pow(10,-6),""));

        title = "h2_"+Label(i)+"_NoSat";
        strcpy(title_char,title.c_str());
        vect_H2_NoSat.push_back(new TH2F(title_char,title_char,300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        title = "h1_"+Label(i)+"_NoSat";
        strcpy(title_char,title.c_str());
        vect_H1_NoSat.push_back(new TH1F(title_char,title_char,200,0,8));
        title = "pr_"+Label(i)+"_NoSat";
        strcpy(title_char,title.c_str());
        vect_prof_NoSat.push_back(new TProfile(title_char,title_char,50,0,900*pow(10,-6),""));

        title = "h2_"+Label(i)+"_Sat";
        strcpy(title_char,title.c_str());
        vect_H2_Sat.push_back(new TH2F(title_char,title_char,300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        title = "h1_"+Label(i)+"_Sat";
        strcpy(title_char,title.c_str());
        vect_H1_Sat.push_back(new TH1F(title_char,title_char,200,0,8));
        title = "pr_"+Label(i)+"_Sat";
        strcpy(title_char,title.c_str());
        vect_prof_Sat.push_back(new TProfile(title_char,title_char,50,0,900*pow(10,-6),""));

        title = "h2_"+Label(i)+"_Sat254";
        strcpy(title_char,title.c_str());
        vect_H2_Sat254.push_back(new TH2F(title_char,title_char,300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        title = "h1_"+Label(i)+"_Sat254";
        strcpy(title_char,title.c_str());
        vect_H1_Sat254.push_back(new TH1F(title_char,title_char,200,0,8));
        title = "pr_"+Label(i)+"_Sat254";
        strcpy(title_char,title.c_str());
        vect_prof_Sat254.push_back(new TProfile(title_char,title_char,50,0,900*pow(10,-6),""));

        title = "h2_"+Label(i)+"_Sat255";
        strcpy(title_char,title.c_str());
        vect_H2_Sat255.push_back(new TH2F(title_char,title_char,300,0,3000*pow(10,-6),300,0,3000*pow(10,-6)));
        title = "h1_"+Label(i)+"_Sat255";
        strcpy(title_char,title.c_str());
        vect_H1_Sat255.push_back(new TH1F(title_char,title_char,200,0,8));
        title = "pr_"+Label(i)+"_Sat255";
        strcpy(title_char,title.c_str());
        vect_prof_Sat255.push_back(new TProfile(title_char,title_char,50,0,3000*pow(10,-6),""));
    }

    vector<TH1F*> VectPtPartID;

    VectPtPartID.push_back(new TH1F("h1_proton_pt","h1_proton_pt",100,0,3000));
    VectPtPartID.push_back(new TH1F("h1_pion_pt","h1_pion_pt",100,0,3000));
    VectPtPartID.push_back(new TH1F("h1_gluino-u-dbar_pt","h1_gluino-u-dbar_pt",100,0,3000));
    VectPtPartID.push_back(new TH1F("h1_R-hadron_pt","h1_R-hadron_pt",100,0,3000));

    vector<TH1F*> VectElossPartID;

    VectElossPartID.push_back(new TH1F("h_proton_eloss","h_proton_eloss",300,0,3000*pow(10,-6)));
    VectElossPartID.push_back(new TH1F("h_pion_eloss","h_pion_eloss",300,0,3000*pow(10,-6)));
    VectElossPartID.push_back(new TH1F("h_gluino-u-dbar_eloss","h_gluino-u-dbar_eloss",300,0,3000*pow(10,-6)));
    VectElossPartID.push_back(new TH1F("h_R-hadron_eloss","h_R-hadron_eloss",300,0,3000*pow(10,-6)));

    vector<TH1F*> VectPoverMPartID;

    VectPoverMPartID.push_back(new TH1F("h_proton_pm","h_proton_pm",1000,0.01,10));
    VectPoverMPartID.push_back(new TH1F("h_pion_pm","h_pion_pm",1000,0.01,10));
    VectPoverMPartID.push_back(new TH1F("h_gluino-u-dbar_pm","h_gluino-u-dbar_pm",1000,0.01,10));
    VectPoverMPartID.push_back(new TH1F("h_R-hadron_pm","h_R-hadron_pm",1000,0.01,10));


    TH1F* hEmpty = new TH1F("hEmpty","hEmpty",25,0,25); //pour tracer les lignes avec la fonction SetHistoLabel

    

    TFile* fileout = new TFile(s2.c_str(),"RECREATE");


    int entries = nentries;


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
                clust++;
                
                vect_eloss.push_back(eloss);

                h2ElossvQ_tot->Fill(eloss,charge);
                hDiff_rel_ElossQ_tot->Fill((eloss-charge)/charge);
                hRatio_ElossQ_tot->Fill(eloss/charge);
                profElossVsQ_tot->Fill(eloss,charge);
                
                hLayer->Fill(layer);
                hLayerLabel->Fill(layerLabel);

                vect_H2_tot[layerLabel-1]->Fill(eloss,charge);
                vect_H1_tot[layerLabel-1]->Fill(eloss/charge);
                vect_prof_tot[layerLabel-1]->Fill(eloss,charge);


                for(int simhit=0;simhit<nsimhits;simhit++)
                {
                    int partID      = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectSimHits()[simhit].GetPartId();
                    
                    if(partID==2212 || partID==-2212) //proton
                    {
                        VectElossPartID[0]->Fill(eloss);
                    }
                    if(partID==211 || partID==-211) //pion
                    {
                        VectElossPartID[1]->Fill(eloss);
                    }
                    if(partID==1009213 || partID==-1009213) //R-hadron gluino u dbar
                    {
                        VectElossPartID[2]->Fill(eloss);
                    }
                    if((int)partID/1000==1009 || (int)partID/1000==-1009) //R-hadron
                    {
                        VectElossPartID[3]->Fill(eloss);
                    }
                    
                }

                
                if(sat254==true && sat255==false) 
                {
                    hRatio_ElossQ_Sat254->Fill(eloss/charge);
                    hDiff_rel_ElossQ_Sat254->Fill((eloss-charge)/charge);
                    h2ElossvQ_Sat254->Fill(eloss,charge);
                    profElossVsQ_Sat254->Fill(eloss,charge);
                    hLayerLabelSat254->Fill(layerLabel);
                    vect_H2_Sat254[layerLabel-1]->Fill(eloss,charge);
                    vect_H1_Sat254[layerLabel-1]->Fill(eloss/charge);
                    vect_prof_Sat254[layerLabel-1]->Fill(eloss,charge);
                    testsat254=true;
                    clustsat254++;
                }
                if(sat255==true)
                {
                    hRatio_ElossQ_Sat255->Fill(eloss/charge);
                    hDiff_rel_ElossQ_Sat255->Fill((eloss-charge)/charge);
                    h2ElossvQ_Sat255->Fill(eloss,charge);
                    profElossVsQ_Sat255->Fill(eloss,charge);
                    hLayerLabelSat255->Fill(layerLabel);
                    vect_H2_Sat255[layerLabel-1]->Fill(eloss,charge);
                    vect_H1_Sat255[layerLabel-1]->Fill(eloss/charge);
                    vect_prof_Sat255[layerLabel-1]->Fill(eloss,charge);
                    testsat255=true;
                    clustsat255++;
                }
                if(sat254==true)
                {
                    h2ElossvQ_Sat->Fill(eloss,charge);
                    hDiff_rel_ElossQ_Sat->Fill((eloss-charge)/charge);
                    hRatio_ElossQ_Sat->Fill(eloss/charge);
                    profElossVsQ_Sat->Fill(eloss,charge);
                    vect_H2_Sat[layerLabel-1]->Fill(eloss,charge);
                    vect_H1_Sat[layerLabel-1]->Fill(eloss/charge);
                    vect_prof_Sat[layerLabel-1]->Fill(eloss,charge);
                }
                if(sat254==false && sat255==false)
                {
                    h2ElossvQ_NoSat->Fill(eloss,charge);
                    hDiff_rel_ElossQ_NoSat->Fill((eloss-charge)/charge);
                    hRatio_ElossQ_NoSat->Fill(eloss/charge);
                    profElossVsQ_NoSat->Fill(eloss,charge);
                    vect_H2_NoSat[layerLabel-1]->Fill(eloss,charge);
                    vect_H1_NoSat[layerLabel-1]->Fill(eloss/charge);
                    vect_prof_NoSat[layerLabel-1]->Fill(eloss,charge);
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
                VectPtPartID[0]->Fill(pt);
                VectPoverMPartID[0]->Fill(GetPoverM(p,id));
            }
            if(id==211 || id==-211) 
            {
                VectPtPartID[1]->Fill(pt);
                VectPoverMPartID[1]->Fill(GetPoverM(p,id));
            }
            if(id==1009213 || id==-1009213) 
            {
                VectPtPartID[2]->Fill(pt);
                VectPoverMPartID[2]->Fill(GetPoverM(p,id));
            }
            if((int)id/1000==1009 || (int)id/1000==-1009)
            {
                VectPtPartID[3]->Fill(pt);
                VectPoverMPartID[3]->Fill(GetPoverM(p,id));
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

            h2RatioSatPt->Fill(pt,(double)NClustSat254/(double)NCluster);
            profSatPt->Fill(pt,(double)NClustSat254/(double)NCluster);

            Estimator estim(vect_eloss);

            h2RatioSatEloss->Fill(estim.GetMean(),(double)NClustSat254/(double)NCluster);
            profSatEloss->Fill(estim.GetMean(),(double)NClustSat254/(double)NCluster);

            h2RatioSatPartID->Fill(ReBinPartID(id),(double)NClustSat254/(double)NCluster);
            profSatPartID->Fill(ReBinPartID(id),(double)NClustSat254/(double)NCluster);

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


    //gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    TCanvas* cPtPartID = new TCanvas("cPtPartID","cPtPartID",700,400);
    vector<char*> VectLegend;
    VectLegend.push_back("p");
    VectLegend.push_back("#pi");
    VectLegend.push_back("R^{+}_{#scale[0.6]{#tilde{g}u#bar{d}}}");
    VectLegend.push_back("R^{+}");
    StackHisto(*cPtPartID,VectPtPartID,VectLegend,"Distribution en p_{T} pour differentes particules","p_{T}");
    cPtPartID->Write();

    TCanvas* cElossPartID = new TCanvas("cElossPartID","cElossPartID",700,400);
    StackHisto(*cElossPartID,VectElossPartID,VectLegend,"Distribution E_{loss} pour differentes particules","E_{loss}");
    cElossPartID->Write();

    TCanvas* cPoverMPartID = new TCanvas("cPoverMPartID","cPoverMPartID",700,400);
    StackHisto(*cPoverMPartID,VectPoverMPartID,VectLegend,"Distribution du P/M pour differentes particules","P/M");
    cPoverMPartID->Write();

    TCanvas* cPoverMPartIDNormal = new TCanvas("cPoverMPartIDNormal","cPoverMPartIDNormal",700,400);
    DrawHistoNormalized(*cPoverMPartIDNormal,VectPoverMPartID,VectLegend,"Distribution du P/M pour differentes particules","P/M");
    cPoverMPartIDNormal->Write();


    DrawHisto(fileout,h2RatioSatPt,"Frac Sat254 = f(Pt)",true,"Pt",true,"#frac{#ClustSat}{#Clust}");
    DrawHisto(fileout,h2RatioSatEloss,"Frac Sat254 = f(E_{loss})",true,"E_{loss}",true,"#frac{#ClustSat}{#Clust}");
    DrawHisto(fileout,h2RatioSatPartID,"Frac Sat254 = f(PartID)",true,"PartID",true,"#frac{#ClustSat}{#Clust}");
    

    TCanvas* cRatioSatPt = new TCanvas("cRatioSatPt","cRatioSatPt",700,400);
    h2RatioSatPt->SetMaximum(200);
    h2RatioSatPt->Draw();
    h2RatioSatPt->SetDrawOption("COLZ");
    profSatPt->SetLineColor(kRed);
    profSatPt->SetMarkerColor(kRed);
    profSatPt->SetMarkerStyle(2);
    profSatPt->Draw("SAME");
    cRatioSatPt->Write();
    profSatPt->Write();

    TCanvas* cRatioSatEloss = new TCanvas("cRatioSatEloss","cRatioSatEloss",700,400);
    h2RatioSatEloss->Draw();
    h2RatioSatEloss->SetDrawOption("COLZ");
    profSatEloss->SetMarkerStyle(2);
    profSatEloss->SetMarkerColor(kRed);
    profSatEloss->SetLineColor(kRed);
    profSatEloss->Draw("SAME");
    cRatioSatEloss->Write();
    profSatEloss->Write();

    DrawHisto(fileout,hDiff_rel_ElossQ_tot,"Difference relative",true,"(E-Q)/Q");
    DrawHisto(fileout,hDiff_rel_ElossQ_NoSat,"Difference relative",true,"(E-Q)/Q");
    DrawHisto(fileout,hDiff_rel_ElossQ_Sat,"Difference relative",true,"(E-Q)/Q");
    DrawHisto(fileout,hDiff_rel_ElossQ_Sat254,"Difference relative",true,"(E-Q)/Q");
    DrawHisto(fileout,hDiff_rel_ElossQ_Sat255,"Difference relative",true,"(E-Q)/Q");

    DrawHisto(fileout,hRatio_ElossQ_tot,"",true,"E/Q");
    DrawHisto(fileout,hRatio_ElossQ_NoSat,"",true,"E/Q");
    DrawHisto(fileout,hRatio_ElossQ_Sat,"",true,"E/Q");
    DrawHisto(fileout,hRatio_ElossQ_Sat254,"",true,"E/Q");
    DrawHisto(fileout,hRatio_ElossQ_Sat255,"",true,"E/Q");

    DrawHisto(fileout,h2ElossvQ_tot,"E_{loss} v. Q",true,"E_{loss}",true,"Q");
    DrawHisto(fileout,h2ElossvQ_NoSat,"E_{loss} v. Q",true,"E_{loss}",true,"Q");
    DrawHisto(fileout,h2ElossvQ_Sat,"E_{loss} v. Q",true,"E_{loss}",true,"Q");
    DrawHisto(fileout,h2ElossvQ_Sat254,"E_{loss} v. Q",true,"E_{loss}",true,"Q");
    DrawHisto(fileout,h2ElossvQ_Sat255,"E_{loss} v. Q",true,"E_{loss}",true,"Q");

    fileout->Append(profElossVsQ_tot);
    fileout->Append(profElossVsQ_NoSat);
    fileout->Append(profElossVsQ_Sat);
    fileout->Append(profElossVsQ_Sat254);
    fileout->Append(profElossVsQ_Sat255);

    fileout->Append(hLayer);
    fileout->Append(hLayerLabel);
    fileout->Append(hLayerLabelSat254);
    fileout->Append(hLayerLabelSat255);

    DrawHisto(fileout,hPt_tot,"Distribution du Pt",true,"Pt");
    DrawHisto(fileout,hPt_NoSat,"Distribution du Pt",true,"Pt");
    DrawHisto(fileout,hPt_Sat,"Distribution du Pt",true,"Pt");
    DrawHisto(fileout,hPt_Sat254,"Distribution du Pt",true,"Pt");
    DrawHisto(fileout,hPt_Sat255,"Distribution du Pt",true,"Pt");

    DrawHisto(fileout,hNClusterPerTrack,"Nombre de clusters pour une trace",true,"# of clusters");
    DrawHisto(fileout,hNClusterSat254PerTrack,"Nombre de clusters sat254 pour une trace",true,"# of clusters");
    DrawHisto(fileout,hNClusterSat255PerTrack,"Nombre de clusters sat255 pour une trace",true,"# of clusters");

    fileout->Append(hRatio_NClusterSat254);
    fileout->Append(hRatio_NClusterSat255);

    TCanvas* canvas_RatioSatPartID = new TCanvas("canvas_RatioSatPartID","canvas_RatioSatPartID",700,400);

    SetHistoLabelPartID(canvas_RatioSatPartID,h2RatioSatPartID);
    profSatPartID->SetMarkerStyle(2);
    profSatPartID->SetMarkerColor(kRed);
    profSatPartID->SetLineColor(kRed);
    profSatPartID->Draw("SAME");
    h2RatioSatPartID->GetXaxis()->SetTitle("");
    h2RatioSatPartID->SetDrawOption("COLZ");
    h2RatioSatPartID->SetMaximum(800);
    canvas_RatioSatPartID->Write();


// ------------------------------------------ ratio saturation fonction de la layer 

    TEfficiency* Ratio_LayerSat254 = new TEfficiency(*hLayerLabelSat254,*hLayerLabel);
    TEfficiency* Ratio_LayerSat255 = new TEfficiency(*hLayerLabelSat255,*hLayerLabel);

    fileout->Append(Ratio_LayerSat254);
    fileout->Append(Ratio_LayerSat255);

    float rat254 = (float)clustsat254/(float)clust;
    float rat255 = (float)clustsat255/(float)clust;
    cout<<"rat254 "<<rat254<<endl;
    cout<<"rat255 "<<rat255<<endl;
    TLine* line254 = new TLine(0,rat254,25,rat254);
    TLine* line255 = new TLine(0,rat255,25,rat255);

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
    canvas_ratioSatLayer->Write();

// ------------------------------------------

// ------------------------------------------ TH2F + TProfile 

    TCanvas* canvas_EvQ_h2prof_tot = new TCanvas("canvas_EvQ_h2prof_tot","canvas_EvQ_h2prof_tot",700,400);
    h2ElossvQ_tot->Draw();
    h2ElossvQ_tot->GetXaxis()->SetTitle("E_{loss}");
    h2ElossvQ_tot->GetYaxis()->SetTitle("Q");
    h2ElossvQ_tot->SetDrawOption("COLZ");
    profElossVsQ_tot->SetMarkerStyle(2);
    profElossVsQ_tot->SetMarkerColor(kRed);
    profElossVsQ_tot->SetLineColor(kRed);
    profElossVsQ_tot->Draw("SAME");
    canvas_EvQ_h2prof_tot->Write();

    TCanvas* canvas_EvQ_h2prof_NoSat = new TCanvas("canvas_EvQ_h2prof_NoSat","canvas_EvQ_h2prof_NoSat",700,400);
    h2ElossvQ_NoSat->Draw();
    h2ElossvQ_NoSat->GetXaxis()->SetTitle("E_{loss}");
    h2ElossvQ_NoSat->GetYaxis()->SetTitle("Q");
    h2ElossvQ_NoSat->SetDrawOption("COLZ");
    profElossVsQ_NoSat->SetMarkerStyle(2);
    profElossVsQ_NoSat->SetMarkerColor(kRed);
    profElossVsQ_NoSat->SetLineColor(kRed);
    profElossVsQ_NoSat->Draw("SAME");
    canvas_EvQ_h2prof_NoSat->Write();

    TCanvas* canvas_EvQ_h2prof_Sat = new TCanvas("canvas_EvQ_h2prof_Sat","canvas_EvQ_h2prof_Sat",700,400);
    h2ElossvQ_Sat->Draw();
    h2ElossvQ_Sat->GetXaxis()->SetTitle("E_{loss}");
    h2ElossvQ_Sat->GetYaxis()->SetTitle("Q");
    h2ElossvQ_Sat->SetDrawOption("COLZ");
    profElossVsQ_Sat->SetMarkerStyle(2);
    profElossVsQ_Sat->SetMarkerColor(kRed);
    profElossVsQ_Sat->SetLineColor(kRed);
    profElossVsQ_Sat->Draw("SAME");
    canvas_EvQ_h2prof_Sat->Write();

    TCanvas* canvas_EvQ_h2prof_Sat254 = new TCanvas("canvas_EvQ_h2prof_Sat254","canvas_EvQ_h2prof_Sat254",700,400);
    h2ElossvQ_Sat254->Draw();
    h2ElossvQ_Sat254->GetXaxis()->SetTitle("E_{loss}");
    h2ElossvQ_Sat254->GetYaxis()->SetTitle("Q");
    h2ElossvQ_Sat254->SetDrawOption("COLZ");
    profElossVsQ_Sat254->SetMarkerStyle(2);
    profElossVsQ_Sat254->SetMarkerColor(kRed);
    profElossVsQ_Sat254->SetLineColor(kRed);
    profElossVsQ_Sat254->Draw("SAME");
    canvas_EvQ_h2prof_Sat254->Write();

    TCanvas* canvas_EvQ_h2prof_Sat255 = new TCanvas("canvas_EvQ_h2prof_Sat255","canvas_EvQ_h2prof_Sat255",700,400);
    h2ElossvQ_Sat255->Draw();
    h2ElossvQ_Sat255->GetXaxis()->SetTitle("E_{loss}");
    h2ElossvQ_Sat255->GetYaxis()->SetTitle("Q");
    h2ElossvQ_Sat255->SetDrawOption("COLZ");
    profElossVsQ_Sat255->SetMarkerStyle(2);
    profElossVsQ_Sat255->SetMarkerColor(kRed);
    profElossVsQ_Sat255->SetLineColor(kRed);
    profElossVsQ_Sat255->Draw("SAME");
    canvas_EvQ_h2prof_Sat255->Write();

    profElossVsQ_tot->SetMarkerStyle(1);
    profElossVsQ_NoSat->SetMarkerStyle(1);
    profElossVsQ_Sat->SetMarkerStyle(1);
    profElossVsQ_Sat254->SetMarkerStyle(1);
    profElossVsQ_Sat255->SetMarkerStyle(1);

// ------------------------------------------

// ------------------------------------------ TH1F Pt de la trace, avec ou sans saturation

    TCanvas* canvas_Pt = new TCanvas("canvas_Pt","canvas_Pt",700,400);
    TLegend* leg_pt = new TLegend(0.7,0.7,0.9,0.9);
    leg_pt->AddEntry(hPt_NoSat,"No saturation","l");
    leg_pt->AddEntry(hPt_Sat,"Saturation","l");
    hPt_NoSat->Scale(1./(hPt_NoSat->Integral()));
    hPt_NoSat->Draw();
    hPt_NoSat->GetXaxis()->SetTitle("Pt");
    hPt_Sat->SetLineColor(2);
    hPt_Sat->Scale(1./(hPt_Sat->Integral()));
    hPt_Sat->Draw("SAME");
    leg_pt->Draw("SAME");
    canvas_Pt->SetLogy();
    gStyle->SetOptStat(0);
    canvas_Pt->Write();

// ------------------------------------------

// ------------------------------------------ TH1F ratio NCluster qui saturent par trace

    TCanvas* canvas_RatioNClusterSat = new TCanvas("canvas_RatioNClusterSat","canvas_RatioNClusterSat",700,400);
    TLegend* leg_RatNCluster = new TLegend(0.7,0.7,0.9,0.9);
    leg_RatNCluster->AddEntry(hRatio_NClusterSat255,"Saturation 255","l");
    leg_RatNCluster->AddEntry(hRatio_NClusterSat254,"Saturation 254","l");
    hRatio_NClusterSat255->Scale(1./(hRatio_NClusterSat255->Integral()));
    hRatio_NClusterSat255->Draw();
    hRatio_NClusterSat254->SetLineColor(2);
    hRatio_NClusterSat254->Scale(1./(hRatio_NClusterSat254->Integral()));
    hRatio_NClusterSat254->Draw("SAME");
    leg_RatNCluster->Draw("SAME");
    hRatio_NClusterSat255->GetXaxis()->SetTitle("#frac{Nombre de cluster sature}{Nombre de cluster dans la trace}");
    hRatio_NClusterSat255->SetTitle("Fraction de saturation");
    canvas_RatioNClusterSat->SetLogy();
    gStyle->SetOptStat(0);
    canvas_RatioNClusterSat->Write();

// ------------------------------------------

// ------------------------------------------ Ratio E/Q

    TCanvas* canvas_EvQ_Ratio = new TCanvas("canvas_EvQ_Ratio","canvas_EvQ_Ratio",700,400);
    TLegend* leg_EvQ_Ratio = new TLegend(0.7,0.7,0.9,0.9);
    double scale_RatioEvQ = 1./(hRatio_ElossQ_NoSat->Integral());
    //leg_EvQ_Ratio->AddEntry(hRatio_ElossQ_tot,"Total","l");
    leg_EvQ_Ratio->AddEntry(hRatio_ElossQ_NoSat,"No saturation","l");
    leg_EvQ_Ratio->AddEntry(hRatio_ElossQ_Sat,"Saturation","l");
    /*hRatio_ElossQ_tot->Scale(1./(hRatio_ElossQ_tot->Integral()));
    hRatio_ElossQ_tot->Draw();
    hRatio_ElossQ_tot->GetXaxis()->SetTitle("E/Q");*/
    hRatio_ElossQ_NoSat->SetLineColor(kCyan);
    hRatio_ElossQ_NoSat->Scale(1./(hRatio_ElossQ_NoSat->Integral()));
    hRatio_ElossQ_NoSat->Draw();
    hRatio_ElossQ_Sat->SetLineColor(kViolet);
    hRatio_ElossQ_Sat->Scale(1./(hRatio_ElossQ_Sat->Integral()));
    hRatio_ElossQ_Sat->Draw("SAME");
    leg_EvQ_Ratio->Draw("SAME");
    canvas_EvQ_Ratio->SetLogy();
    gStyle->SetOptStat(0);
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


    for(int i=0;i<21;i++)
    {
        string title = "E_{loss} v. Q  |  "+Label(i+1);
        char title_char[title.length()+1];
        strcpy(title_char,title.c_str());
        TCanvas* csuperposed = new TCanvas(title_char,title_char,700,400);
        DrawHisto(fileout,vect_H2_tot[i],title_char,true,"E_{loss}",true,"Q");
        DrawHisto(fileout,vect_H2_NoSat[i],title_char,true,"E_{loss}",true,"Q");
        DrawHisto(fileout,vect_H2_Sat[i],title_char,true,"E_{loss}",true,"Q");
        DrawHisto(fileout,vect_H2_Sat254[i],title_char,true,"E_{loss}",true,"Q");
        DrawHisto(fileout,vect_H2_Sat255[i],title_char,true,"E_{loss}",true,"Q");
        DrawHisto(fileout,vect_H1_tot[i],title_char,true,"E/Q");
        DrawHisto(fileout,vect_H1_NoSat[i],title_char,true,"E/Q");
        DrawHisto(fileout,vect_H1_Sat[i],title_char,true,"E/Q");
        DrawHisto(fileout,vect_H1_Sat254[i],title_char,true,"E/Q");
        DrawHisto(fileout,vect_H1_Sat255[i],title_char,true,"E/Q");
        SuperposedHisto2DProfile(*csuperposed,vect_H2_tot[i],vect_prof_tot[i],title_char,"E_{loss}","Q");
        DrawCanvas(fileout,*csuperposed);
        SuperposedHisto2DProfile(*csuperposed,vect_H2_NoSat[i],vect_prof_NoSat[i],title_char,"E_{loss}","Q");
        DrawCanvas(fileout,*csuperposed);
        SuperposedHisto2DProfile(*csuperposed,vect_H2_Sat[i],vect_prof_Sat[i],title_char,"E_{loss}","Q");
        DrawCanvas(fileout,*csuperposed);
        SuperposedHisto2DProfile(*csuperposed,vect_H2_Sat254[i],vect_prof_Sat254[i],title_char,"E_{loss}","Q");
        DrawCanvas(fileout,*csuperposed);
        SuperposedHisto2DProfile(*csuperposed,vect_H2_Sat255[i],vect_prof_Sat255[i],title_char,"E_{loss}","Q");
        DrawCanvas(fileout,*csuperposed);
    }
    

    fileout->Write();
    fileout->Close();

    delete fileout;

    string s3 = "mv "+s2+" ./data/.";

    cout<<s3<<endl;

    system(s3.c_str());
    

    return 0;
}
