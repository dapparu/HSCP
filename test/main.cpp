#include <iostream>
#include <vector>
#include <stdlib.h>
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

#include "../interface/Estimator.h"
#include "../interface/Builder.h"
#include "../interface/PlotterHisto.h"

using namespace std;

void SetHistoLabel(TCanvas* canvas,TH1F* histo)
{

    int Wsize = histo->GetMaximum()*1.05;
    histo->GetYaxis()->SetRangeUser(0,histo->GetMaximum()*1.05);

    TLine* line1 = new TLine(5,0,5,0.25);
    TLine* line2 = new TLine(11,0,11,0.25);
    TLine* line3 = new TLine(13,0,13,0.25);

    if(Wsize!=0)
    {
        line1->SetY2(Wsize);
        line2->SetY2(Wsize);
        line3->SetY2(Wsize);
    }


    line1->SetLineStyle(2);
    line2->SetLineStyle(2);
    line3->SetLineStyle(2);
    
    histo->GetXaxis()->SetBinLabel(1+1,"TIB 1");
    histo->GetXaxis()->SetBinLabel(1+2,"TIB 2");
    histo->GetXaxis()->SetBinLabel(1+3,"TIB 3");
    histo->GetXaxis()->SetBinLabel(1+4,"TIB 4");
    histo->GetXaxis()->SetBinLabel(1+5,"TOB 1");
    histo->GetXaxis()->SetBinLabel(1+6,"TOB 2");
    histo->GetXaxis()->SetBinLabel(1+7,"TOB 3");
    histo->GetXaxis()->SetBinLabel(1+8,"TOB 4");
    histo->GetXaxis()->SetBinLabel(1+9,"TOB 5");
    histo->GetXaxis()->SetBinLabel(1+10,"TOB 6");
    histo->GetXaxis()->SetBinLabel(1+11,"TID 1");
    histo->GetXaxis()->SetBinLabel(1+12,"TID 2");
    histo->GetXaxis()->SetBinLabel(1+13,"TEC 1");
    histo->GetXaxis()->SetBinLabel(1+14,"TEC 2");
    histo->GetXaxis()->SetBinLabel(1+15,"TEC 3");
    histo->GetXaxis()->SetBinLabel(1+16,"TEC 4");
    histo->GetXaxis()->SetBinLabel(1+17,"TEC 5");
    histo->GetXaxis()->SetBinLabel(1+18,"TEC 6");
    histo->GetXaxis()->SetBinLabel(1+19,"TEC 7");
    histo->GetXaxis()->SetBinLabel(1+20,"TEC 8");
    histo->GetXaxis()->SetBinLabel(1+21,"TEC 9");
    histo->GetXaxis()->LabelsOption("v");
    canvas->cd();
    histo->Draw();
    line1->Draw();
    line2->Draw();
    line3->Draw();
    
}


string Label(int i)
{
    if(i==1)        return "TIB1";
    else if(i==2)   return "TIB2";
    else if(i==3)   return "TIB3";
    else if(i==4)   return "TIB4";
    else if(i==5)   return "TOB1";
    else if(i==6)   return "TOB2";
    else if(i==7)   return "TOB3";
    else if(i==8)   return "TOB4";
    else if(i==9)   return "TOB5";
    else if(i==10)  return "TOB6";
    else if(i==11)  return "TID1";
    else if(i==12)  return "TID2";
    else if(i==13)  return "TEC1";
    else if(i==14)  return "TEC2";
    else if(i==15)  return "TEC3";
    else if(i==16)  return "TEC4";
    else if(i==17)  return "TEC5";
    else if(i==18)  return "TEC6";
    else if(i==19)  return "TEC7";
    else if(i==20)  return "TEC8";
    else if(i==21)  return "TEC9";
}

string LabelParticle(int i)
{
    if(i==111 || i==211) return "#pi";
    else if(i==2212) return "p";
    else if(i==2112) return "n";
    else if(i==1009213) return "R^{+}_{#tilde{g}u#bar{d}}";
    else if((int)i/1000==1009) return "R^{+}";
}




int GetPartID(const vector<Cluster> &VectClust,float &threshold)
{
    int PartID=0;
    vector<int> VectPartID_Cluster;
    for(int cluster=0;cluster<VectClust.size();cluster++)
    {
        for(int simhit=0;simhit<VectClust[cluster].GetVectSimHits().size();simhit++)
        {
            PartID=VectClust[cluster].GetVectSimHits()[simhit].GetPartId();
            VectPartID_Cluster.push_back(PartID);
        }
    }
    sort(VectPartID_Cluster.begin(),VectPartID_Cluster.end());
    vector<int> VectPartID;
    vector<int> counter;
    VectPartID.push_back(VectPartID_Cluster[0]);
    counter.push_back(1);
    int j=0;
    for(int i=1;i<VectPartID_Cluster.size();i++)
    {
        if(VectPartID_Cluster[i]==VectPartID[j])
        {
            counter[j]++;
        }
        else
        {
            VectPartID.push_back(VectPartID_Cluster[i]);
            counter.push_back(1);
            j++;
        }
    }
    vector<float> ratio;
    for(int i=0;i<counter.size();i++)
    {
        ratio.push_back((float)counter[i]/(float)VectPartID_Cluster.size());
        //cout<<counter[i]<<endl;
        //cout<<VectPartID.size()<<endl;
    }
    float max=0;
    int indice=0;
    for(int i=0;i<ratio.size();i++)
    {
        if(ratio[i]>max)
        {
            max=ratio[i];
            indice=i;
        }
    }
    threshold=max;
    return VectPartID[indice];
}


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
    TH2F* h2ElossvQ_Sat255 = new TH2F("h2ElossvQ_Sat255","h2ElossvQ_Sat255",300,0,900*pow(10,-6),300,0,900*pow(10,-6));

    TProfile* profElossVsQ_tot = new TProfile("profElossVsQ_tot","profElossVsQ_tot",50,0,900*pow(10,-6),"");
    TProfile* profElossVsQ_NoSat = new TProfile("profElossVsQ_NoSat","profElossVsQ_NoSat",50,0,900*pow(10,-6),"");
    TProfile* profElossVsQ_Sat = new TProfile("profElossVsQ_Sat","profElossVsQ_Sat",50,0,900*pow(10,-6),"");
    TProfile* profElossVsQ_Sat254 = new TProfile("profElossVsQ_Sat254","profElossVsQ_Sat254",50,0,900*pow(10,-6),"");
    TProfile* profElossVsQ_Sat255 = new TProfile("profElossVsQ_Sat255","profElossVsQ_Sat255",50,0,900*pow(10,-6),"");

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

    vector<TH2F*> vect_H2_tot;
    vector<TH2F*> vect_H2_NoSat;
    vector<TH2F*> vect_H2_Sat;
    vector<TH2F*> vect_H2_Sat254;
    vector<TH2F*> vect_H2_Sat255;

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

        title = "h2_"+Label(i)+"_NoSat";
        strcpy(title_char,title.c_str());
        vect_H2_NoSat.push_back(new TH2F(title_char,title_char,300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        title = "h1_"+Label(i)+"_NoSat";
        strcpy(title_char,title.c_str());
        vect_H1_NoSat.push_back(new TH1F(title_char,title_char,200,0,8));

        title = "h2_"+Label(i)+"_Sat";
        strcpy(title_char,title.c_str());
        vect_H2_Sat.push_back(new TH2F(title_char,title_char,300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        title = "h1_"+Label(i)+"_Sat";
        strcpy(title_char,title.c_str());
        vect_H1_Sat.push_back(new TH1F(title_char,title_char,200,0,8));

        title = "h2_"+Label(i)+"_Sat254";
        strcpy(title_char,title.c_str());
        vect_H2_Sat254.push_back(new TH2F(title_char,title_char,300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        title = "h1_"+Label(i)+"_Sat254";
        strcpy(title_char,title.c_str());
        vect_H1_Sat254.push_back(new TH1F(title_char,title_char,200,0,8));

        title = "h2_"+Label(i)+"_Sat255";
        strcpy(title_char,title.c_str());
        vect_H2_Sat255.push_back(new TH2F(title_char,title_char,300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        title = "h1_"+Label(i)+"_Sat255";
        strcpy(title_char,title.c_str());
        vect_H1_Sat255.push_back(new TH1F(title_char,title_char,200,0,8));
    }

    vector<TH1F*> VectPtPartID;

    VectPtPartID.push_back(new TH1F("h1_proton_pt","h1_proton_pt",100,0,3000));
    VectPtPartID.push_back(new TH1F("h1_pion_pt","h1_pion_pt",100,0,3000));
    VectPtPartID.push_back(new TH1F("h1_gluino-u-dbar_pt","h1_gluino-u-dbar_pt",100,0,3000));
    VectPtPartID.push_back(new TH1F("h1_R-hadron_pt","h1_R-hadron_pt",100,0,3000));

    vector<TH1F*> VectElossPartID;

    VectElossPartID.push_back(new TH1F("h_proton_eloss","h_proton_eloss",300,0,900*pow(10,-6)));
    VectElossPartID.push_back(new TH1F("h_pion_eloss","h_pion_eloss",300,0,900*pow(10,-6)));
    VectElossPartID.push_back(new TH1F("h_gluino-u-dbar_eloss","h_gluino-u-dbar_eloss",300,0,900*pow(10,-6)));
    VectElossPartID.push_back(new TH1F("h_R-hadron_eloss","h_R-hadron_eloss",300,0,900*pow(10,-6)));

    TH1F* hEmpty = new TH1F("hEmpty","hEmpty",25,0,25); //pour tracer les lignes avec la fonction SetHistoLabel

    

    TFile* fileout = new TFile(s2.c_str(),"RECREATE");


    int entries = 100;


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
            testsat254              = false;
            testsat255              = false;
            float pt                = b1->GetVectTrack()[track].GetPt();
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
                
                h2ElossvQ_tot->Fill(eloss,charge);
                hDiff_rel_ElossQ_tot->Fill((eloss-charge)/charge);
                hRatio_ElossQ_tot->Fill(eloss/charge);
                profElossVsQ_tot->Fill(eloss,charge);
                
                hLayer->Fill(layer);
                hLayerLabel->Fill(layerLabel);

                vect_H2_tot[layerLabel-1]->Fill(eloss,charge);
                vect_H1_tot[layerLabel-1]->Fill(eloss/charge);


                for(int simhit=0;simhit<nsimhits;simhit++)
                {
                    int partID      = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectSimHits()[simhit].GetPartId();
                    
                    if(partID==2212) //proton
                    {
                        VectElossPartID[0]->Fill(eloss);
                    }
                    if(partID==211) //pion
                    {
                        VectElossPartID[1]->Fill(eloss);
                    }
                    if(partID==1009213) //R-hadron gluino u dbar
                    {
                        VectElossPartID[2]->Fill(eloss);
                    }
                    if((int)partID/1000==1009) //R-hadron
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
                }
                else if(sat254==false && sat255==false)
                {
                    h2ElossvQ_NoSat->Fill(eloss,charge);
                    hDiff_rel_ElossQ_NoSat->Fill((eloss-charge)/charge);
                    hRatio_ElossQ_NoSat->Fill(eloss/charge);
                    profElossVsQ_NoSat->Fill(eloss,charge);
                    vect_H2_NoSat[layerLabel-1]->Fill(eloss,charge);
                    vect_H1_NoSat[layerLabel-1]->Fill(eloss/charge);
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
            
            if(id==2212) 
            {
                VectPtPartID[0]->Fill(pt);
            }
            if(id==211) 
            {
                VectPtPartID[1]->Fill(pt);
            }
            if(id==1009213) 
            {
                VectPtPartID[2]->Fill(pt);
            }
            if((int)id/1000==1009)
            {
                VectPtPartID[3]->Fill(pt);
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
            else if(testsat254==false && testsat255==false)
            {
                hPt_NoSat->Fill(pt);
            }

            h2RatioSatPt->Fill(pt,(double)NClustSat254/(double)NCluster);
            profSatPt->Fill(pt,(double)NClustSat254/(double)NCluster);


            /*if(testsat254 && testsat255 )//&& b1->GetVectTrack()[track].GetVectClusters().size() <=3)
            {
                TCanvas* c_profClust = new TCanvas();
                c_profClust = &b1->GetVectTrack()[track].GetProfCluster();
                string str_pt = to_string(pt);
                TText* text_pt = new TText(c_profClust->GetWw()-100,c_profClust->GetWw()-100,str_pt.c_str());
                text_pt->Draw("SAME");
                c_profClust->SaveAs("./Prof_Clust.pdf");
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

    TH2F* h2RatioSatPtClone = (TH2F*) h2RatioSatPt->Clone();

    DrawHisto(fileout,h2RatioSatPt,"Frac Sat254 = f(Pt)",true,"Pt",true,"#frac{#ClustSat}{#Clust}");
    

    TCanvas* cRatioSatPt = new TCanvas("cRatioSatPt","cRatioSatPt",700,400);
    h2RatioSatPtClone->SetMaximum(200);
    h2RatioSatPtClone->Draw();
    h2RatioSatPtClone->SetDrawOption("COLZ");
    profSatPt->SetLineColor(kRed);
    profSatPt->SetMarkerStyle(2);
    profSatPt->Draw("SAME");
    cRatioSatPt->Write();
    profSatPt->Write();

    DrawHisto(fileout,hDiff_rel_ElossQ_tot,"Difference relative",true,"(E-Q)/Q");
    DrawHisto(fileout,hDiff_rel_ElossQ_NoSat,"Difference relative",true,"(E-Q)/Q");
    DrawHisto(fileout,hDiff_rel_ElossQ_Sat,"Difference relative",true,"(E-Q)/Q");
    DrawHisto(fileout,hDiff_rel_ElossQ_Sat254,"Difference relative",true,"(E-Q)/Q");
    DrawHisto(fileout,hDiff_rel_ElossQ_Sat255,"Difference relative",true,"(E-Q)/Q");

    DrawHisto(fileout,hRatio_ElossQ_tot,"Ratio",true,"Q/E");
    DrawHisto(fileout,hRatio_ElossQ_NoSat,"Ratio",true,"Q/E");
    DrawHisto(fileout,hRatio_ElossQ_Sat,"Ratio",true,"Q/E");
    DrawHisto(fileout,hRatio_ElossQ_Sat254,"Ratio",true,"Q/E");
    DrawHisto(fileout,hRatio_ElossQ_Sat255,"Ratio",true,"Q/E");

    DrawHisto(fileout,h2ElossvQ_tot,"E_loss v. Q",true,"E_loss",true,"Q");
    DrawHisto(fileout,h2ElossvQ_NoSat,"E_loss v. Q",true,"E_loss",true,"Q");
    DrawHisto(fileout,h2ElossvQ_Sat,"E_loss v. Q",true,"E_loss",true,"Q");
    DrawHisto(fileout,h2ElossvQ_Sat254,"E_loss v. Q",true,"E_loss",true,"Q");
    DrawHisto(fileout,h2ElossvQ_Sat255,"E_loss v. Q",true,"E_loss",true,"Q");

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
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
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
    h2ElossvQ_tot->GetXaxis()->SetTitle("E_loss");
    h2ElossvQ_tot->GetYaxis()->SetTitle("Q");
    h2ElossvQ_tot->SetDrawOption("COLZ");
    profElossVsQ_tot->SetMarkerStyle(2);
    profElossVsQ_tot->SetMarkerColor(kRed);
    profElossVsQ_tot->Draw("SAME");
    canvas_EvQ_h2prof_tot->Write();

    TCanvas* canvas_EvQ_h2prof_NoSat = new TCanvas("canvas_EvQ_h2prof_NoSat","canvas_EvQ_h2prof_NoSat",700,400);
    h2ElossvQ_NoSat->Draw();
    h2ElossvQ_NoSat->GetXaxis()->SetTitle("E_loss");
    h2ElossvQ_NoSat->GetYaxis()->SetTitle("Q");
    h2ElossvQ_NoSat->SetDrawOption("COLZ");
    profElossVsQ_NoSat->SetMarkerStyle(2);
    profElossVsQ_NoSat->Draw("SAME");
    canvas_EvQ_h2prof_NoSat->Write();

    TCanvas* canvas_EvQ_h2prof_Sat = new TCanvas("canvas_EvQ_h2prof_Sat","canvas_EvQ_h2prof_Sat",700,400);
    h2ElossvQ_Sat->Draw();
    h2ElossvQ_Sat->GetXaxis()->SetTitle("E_loss");
    h2ElossvQ_Sat->GetYaxis()->SetTitle("Q");
    h2ElossvQ_Sat->SetDrawOption("COLZ");
    profElossVsQ_Sat->SetMarkerStyle(2);
    profElossVsQ_Sat->Draw("SAME");
    canvas_EvQ_h2prof_Sat->Write();

    TCanvas* canvas_EvQ_h2prof_Sat254 = new TCanvas("canvas_EvQ_h2prof_Sat254","canvas_EvQ_h2prof_Sat254",700,400);
    h2ElossvQ_Sat254->Draw();
    h2ElossvQ_Sat254->GetXaxis()->SetTitle("E_loss");
    h2ElossvQ_Sat254->GetYaxis()->SetTitle("Q");
    h2ElossvQ_Sat254->SetDrawOption("COLZ");
    profElossVsQ_Sat254->SetMarkerStyle(2);
    profElossVsQ_Sat254->Draw("SAME");
    canvas_EvQ_h2prof_Sat254->Write();

    TCanvas* canvas_EvQ_h2prof_Sat255 = new TCanvas("canvas_EvQ_h2prof_Sat255","canvas_EvQ_h2prof_Sat255",700,400);
    h2ElossvQ_Sat255->Draw();
    h2ElossvQ_Sat255->GetXaxis()->SetTitle("E_loss");
    h2ElossvQ_Sat255->GetYaxis()->SetTitle("Q");
    h2ElossvQ_Sat255->SetDrawOption("COLZ");
    profElossVsQ_Sat255->SetMarkerStyle(2);
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
    double scale_pt = 1./(hPt_NoSat->Integral());
    leg_pt->AddEntry(hPt_NoSat,"No saturation","l");
    leg_pt->AddEntry(hPt_Sat,"Saturation","l");
    hPt_NoSat->Scale(scale_pt);
    hPt_NoSat->Draw();
    hPt_NoSat->GetXaxis()->SetTitle("Pt");
    hPt_Sat->SetLineColor(2);
    hPt_Sat->Scale(scale_pt);
    hPt_Sat->Draw("SAME");
    leg_pt->Draw("SAME");
    canvas_Pt->SetLogy();
    canvas_Pt->Write();

// ------------------------------------------

// ------------------------------------------ TH1F ratio NCluster qui saturent par trace

    TCanvas* canvas_RatioNClusterSat = new TCanvas("canvas_RatioNClusterSat","canvas_RatioNClusterSat",700,400);
    TLegend* leg_RatNCluster = new TLegend(0.7,0.7,0.9,0.9);
    leg_RatNCluster->AddEntry(hRatio_NClusterSat255,"Saturation 255","l");
    leg_RatNCluster->AddEntry(hRatio_NClusterSat254,"Saturation 254","l");
    hRatio_NClusterSat255->Draw();
    hRatio_NClusterSat254->SetLineColor(2);
    hRatio_NClusterSat254->Draw("SAME");
    leg_RatNCluster->Draw("SAME");
    hRatio_NClusterSat255->GetXaxis()->SetTitle("NSat/NTot");
    canvas_RatioNClusterSat->SetLogy();
    canvas_RatioNClusterSat->Write();

// ------------------------------------------

// ------------------------------------------ Ratio E/Q

    TCanvas* canvas_EvQ_Ratio = new TCanvas("canvas_EvQ_Ratio","canvas_EvQ_Ratio",700,400);
    TLegend* leg_EvQ_Ratio = new TLegend(0.7,0.7,0.9,0.9);
    double scale_RatioEvQ = 1./(hRatio_ElossQ_NoSat->Integral());
    leg_EvQ_Ratio->AddEntry(hRatio_ElossQ_tot,"Total","l");
    leg_EvQ_Ratio->AddEntry(hRatio_ElossQ_NoSat,"No saturation","l");
    leg_EvQ_Ratio->AddEntry(hRatio_ElossQ_Sat,"Saturation","l");
    hRatio_ElossQ_tot->Scale(scale_RatioEvQ);
    hRatio_ElossQ_tot->Draw();
    hRatio_ElossQ_tot->GetXaxis()->SetTitle("E/Q");
    hRatio_ElossQ_NoSat->SetLineColor(2);
    hRatio_ElossQ_NoSat->Scale(scale_RatioEvQ);
    hRatio_ElossQ_NoSat->Draw("SAME");
    hRatio_ElossQ_Sat->SetLineColor(3);
    hRatio_ElossQ_Sat->Scale(scale_RatioEvQ);
    hRatio_ElossQ_Sat->Draw("SAME");
    leg_EvQ_Ratio->Draw("SAME");
    canvas_EvQ_Ratio->SetLogy();
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
        string title = "E_loss v. Q  |  "+Label(i+1);
        char title_char[title.length()+1];
        strcpy(title_char,title.c_str());
        DrawHisto(fileout,vect_H2_tot[i],title_char,true,"E_loss",true,"Q");
        DrawHisto(fileout,vect_H2_NoSat[i],title_char,true,"E_loss",true,"Q");
        DrawHisto(fileout,vect_H2_Sat[i],title_char,true,"E_loss",true,"Q");
        DrawHisto(fileout,vect_H2_Sat254[i],title_char,true,"E_loss",true,"Q");
        DrawHisto(fileout,vect_H2_Sat255[i],title_char,true,"E_loss",true,"Q");
        DrawHisto(fileout,vect_H1_tot[i],title_char,true,"E/Q");
        DrawHisto(fileout,vect_H1_NoSat[i],title_char,true,"E/Q");
        DrawHisto(fileout,vect_H1_Sat[i],title_char,true,"E/Q");
        DrawHisto(fileout,vect_H1_Sat254[i],title_char,true,"E/Q");
        DrawHisto(fileout,vect_H1_Sat255[i],title_char,true,"E/Q");
    }
    

    fileout->Write();
    fileout->Close();

    delete fileout;

    string s3 = "mv "+s2+" ./data/.";

    cout<<s3<<endl;

    system(s3.c_str());
    

    return 0;
}
