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
#include "TH1D.h"
#include "TH1D.h"
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
#include "TPaveText.h"
#include "TMatrix.h"

#include "../interface/Builder.h"
#include "../interface/Estimator.h"
#include "../interface/PlotterHisto.h"
#include "../interface/Labellizer.h" 
#include "../interface/Correction.h"

using namespace std;


// Differents facteurs K et C utilises

//const float K = 2.56; //parametres du papier
//const float C = 2.77;

//const float C = 2.64; //obtenus avec fit du proton dqcorr
//const float K = 2.78;

//const float C = 3.20573; //obtenus avec fit eloss susy
//const float K = 1.91189;

const float CsusySim = 4.53; //pour gluino
const float KsusySim = 1.76;

const float KmsSim = 2.28; //pour MS
const float CmsSim = 2.82;

const float CsusyRec = 7.4; //pour gluino
const float KsusyRec = 0.26;

const float KmsRec = 2.49; //pour MS
const float CmsRec = 2.66;

const float CsusyCorr = 1.14; //pour gluino
const float KsusyCorr = 6.7;

const float KmsCorr = 2.5; //pour MS
const float CmsCorr = 2.66;


// constantes utilisees pour la significance 

const float CrossSectionMinBias = 68*pow(10,-3);
const float CrossSectionGluino = 0.128*pow(10,-15);
const float L_data = 140*pow(10,15);
const float factor_signal = L_data*(CrossSectionGluino/15188);
const float factor_bkg = L_data*(CrossSectionMinBias/261000);
const float factor_signal_before = L_data*(CrossSectionGluino/15188);
const float factor_bkg_before = L_data*(CrossSectionMinBias/261000);



std::vector<int> CrossTalkInv(const std::vector<int>&  Q, const float x1=0.10, const float x2=0.04, float threshold=0.25, float thresholdSat=0.20, int label=0, int correctType=3) {
  
    bool DylanCorr=true;
    bool inversion=true;
    if(correctType==1){DylanCorr=false;inversion=false;}
    if(correctType==2){DylanCorr=false;inversion=true;}
    if(correctType==3){DylanCorr=true;inversion=false;}
    if(correctType==4){DylanCorr=true;inversion=true;}

  const unsigned N=Q.size();
  std::vector<int> QII;
  std::vector<float> QI(N,0);
  Double_t a=1-2*x1-2*x2;
  //  bool debugbool=false;
  TMatrix A(N,N);

  //---  que pour 1 max bien net 
 if(Q.size()<2 || Q.size()>8){
	for (unsigned int i=0;i<Q.size();i++){
		QII.push_back((int) Q[i]);
  	}
	return QII;
  }
 if(DylanCorr==false){ 
	  vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end())	;
	  if(*mQ>253){
	 	 if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253 ) return Q ;
	 	 if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40 ){
		     QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2); return QII;}
	  }
  }

  //---

  if(inversion)
  {
  for(unsigned int i=0; i<N; i++) {
          A(i,i) =a;
          if(i<N-1){ A(i+1,i)=x1;A(i,i+1)=x1;}
          else continue; 
          if(i<N-2){ A(i+2,i)=x2;A(i,i+2)=x2;}
    }
  
    if(N==1) A(0,0)=1/a;
    else  A.InvertFast();
  
    for(unsigned int i=0; i<N; i++) {
          for(unsigned int j=0; j<N; j++) {
          QI[i]+=A(i,j)*(float)Q[j];
          }
    }
  
   for (unsigned int i=0;i<QI.size();i++){
  	if(QI[i]<threshold) QI[i]=0; 
  	QII.push_back((int) QI[i]);
    }
  
  return QII;
  }


  if(DylanCorr)
  {
      int NSat254=0;
      int NSat255=0;
      int ClusterCharge=0;
  
      for(int i=0;i<Q.size();i++)
      {
          if(Q[i]==254) NSat254++;
          if(Q[i]==255) NSat255++;
          ClusterCharge+=Q[i];
      }
      TFile* file = TFile::Open("~/Documents/Internship/HSCP/dEdxCorr_Parameters.root");
      TTree* tree = (TTree*) file->Get("tree");
      Correction corr;
      corr.SetTree(*tree);
      
      //ClusterCharge=corr.ChargeCorr(ClusterCharge,label,Q.size(),NSat254,NSat255);
  }

}




int main(int argc,char** argv){

    TChain chain;
	chain.SetName("stage/ttree");
	for(int i=1;i<argc-1;i++) chain.Add(argv[i]); //on utilise les fichiers jusqu'a l'avant dernier 

    string s1 = argv[1];
    string s2 = s1.substr(0,s1.find('.'))+"_results.root"; //notre fichier de sortie prend le nom de l'avant dernier +_results.root
    //string s3 = s1.substr(0,s1.find('.'))+"_FitRes.root";
    string s3 = "~/Documents/Internship/merge_FitRes.root";

    float ratiosat = atof(argv[argc-1]); //le ratio de stauration utilise est le dernier argument donne

    Builder* b1 = new Builder(chain);
    b1->SetBranchAdd();
    int nentries = b1->GetEntries();




    // declaration des differents histogrammes 




    TH1D* hDiff_rel_ElossQ_tot = new TH1D("hDiff_rel_ElossQ_tot","hDiff_rel_ElossQ_tot",200,-1,5);
    TH1D* hDiff_rel_ElossQ_NoSat = new TH1D("hDiff_rel_ElossQ_NoSat","hDiff_rel_ElossQ_NoSat",200,-1,5);
    TH1D* hDiff_rel_ElossQ_Sat = new TH1D("hDiff_rel_ElossQ_Sat","hDiff_rel_ElossQ_Sat",200,-1,5);
    TH1D* hDiff_rel_ElossQ_Sat254 = new TH1D("hDiff_rel_ElossQ_Sat254","hDiff_rel_ElossQ_Sat254",200,-1,5);
    TH1D* hDiff_rel_ElossQ_Sat255 = new TH1D("hDiff_rel_ElossQ_Sat255","hDiff_rel_ElossQ_Sat255",200,-1,5);

    TH1D* hRatio_ElossQ_tot = new TH1D("hRatio_ElossQ_tot","hRatio_ElossQ_tot",200,0,8);
    TH1D* hRatio_ElossQ_NoSat = new TH1D("hRatio_ElossQ_NoSat","hRatio_ElossQ_NoSat",200,0,8);
    TH1D* hRatio_ElossQ_Sat = new TH1D("hRatio_ElossQ_Sat","hRatio_ElossQ_Sat",200,0,8);
    TH1D* hRatio_ElossQ_Sat254 = new TH1D("hRatio_ElossQ_Sat254","hRatio_ElossQ_Sat254",200,0,8);
    TH1D* hRatio_ElossQ_Sat255 = new TH1D("hRatio_ElossQ_Sat255","hRatio_ElossQ_Sat255",200,0,8);

    TH2F* h2ElossvQ_tot = new TH2F("h2ElossvQ_tot","h2ElossvQ_tot",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
    TH2F* h2ElossvQ_NoSat = new TH2F("h2ElossvQ_NoSat","h2ElossvQ_NoSat",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
    TH2F* h2ElossvQ_Sat = new TH2F("h2ElossvQ_Sat","h2ElossvQ_Sat",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
    TH2F* h2ElossvQ_Sat254 = new TH2F("h2ElossvQ_Sat254","h2ElossvQ_Sat254",300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
    TH2F* h2ElossvQ_Sat255 = new TH2F("h2ElossvQ_Sat255","h2ElossvQ_Sat255",300,0,3000*pow(10,-6),300,0,3000*pow(10,-6));

    TProfile* profElossVsQ_tot = new TProfile("profElossVsQ_tot","profElossVsQ_tot",50,0,1500*pow(10,-6),"");
    TProfile* profElossVsQ_NoSat = new TProfile("profElossVsQ_NoSat","profElossVsQ_NoSat",50,0,1500*pow(10,-6),"");
    TProfile* profElossVsQ_Sat = new TProfile("profElossVsQ_Sat","profElossVsQ_Sat",50,0,1500*pow(10,-6),"");
    TProfile* profElossVsQ_Sat254 = new TProfile("profElossVsQ_Sat254","profElossVsQ_Sat254",50,0,1500*pow(10,-6),"");
    TProfile* profElossVsQ_Sat255 = new TProfile("profElossVsQ_Sat255","profElossVsQ_Sat255",50,0,3000*pow(10,-6),"");

    TH1D* hLayer = new TH1D("hLayer","hLayer",80,0,80);
    TH1D* hLayerLabel = new TH1D("hLayerLabel","hLayerLabel",21,0,21);
    TH1D* hLayerLabelSat254 = new TH1D("hLayerLabelSat254","hLayerLabelSat254",21,0,21);
    TH1D* hLayerLabelSat255 = new TH1D("hLayerLabelSat255","hLayerLabelSat255",21,0,21);

    TH1D* hPt_tot = new TH1D("hPt_tot","hPt_tot",100,0,3000);
    TH1D* hPt_NoSat = new TH1D("hPt_NoSat","hPt_NoSat",100,0,3000);
    TH1D* hPt_Sat = new TH1D("hPt_Sat","hPt_Sat",100,0,3000);
    TH1D* hPt_Sat254 = new TH1D("hPt_Sat254","hPt_Sat254",100,0,3000);
    TH1D* hPt_Sat255 = new TH1D("hPt_Sat255","hPt_Sat255",100,0,3000);

    TH1D* hNClusterPerTrack = new TH1D("hNClusterPerTrack","hNClusterPerTrack",25,0,25);
    TH1D* hNClusterSat254PerTrack = new TH1D("hNClusterSat254PerTrack","hNClusterSat254PerTrack",20,0,20);
    TH1D* hNClusterSat255PerTrack = new TH1D("hNClusterSat255PerTrack","hNClusterSat255PerTrack",20,0,20);
    TH1D* hRatio_NClusterSat254 = new TH1D("hRatio_NClusterSat254","hRatio_NClusterSat254",20,0,1.05);
    TH1D* hRatio_NClusterSat255 = new TH1D("hRatio_NClusterSat255","hRatio_NClusterSat255",20,0,1.05);

    TH2F* h2RatioSatPt254 = new TH2F("h2RatioSatPt254","h2RatioSatPt254",100,0,3000,20,0,1.05);
    TProfile* profSatPt254 = new TProfile("profSatPt254","profSatPt254",100,0,3000,"");

    TH2F* h2RatioSatEloss254 = new TH2F("h2RatioSatEloss254","h2RatioSatEloss254",300,0,900*pow(10,-6),20,0,1.05);
    TProfile* profSatEloss254 = new TProfile("profSatEloss254","profSatEloss254",300,0,900*pow(10,-6),"");

    TH2F* h2RatioSatQ254 = new TH2F("h2RatioSatQ254","h2RatioSatQ254",300,0,1200*pow(10,-6),20,0,1.05);
    TProfile* profSatQ254 = new TProfile("profSatQ254","profSatQ254",300,0,1200*pow(10,-6),"");

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
    
    TH2F* h2RatioSatQ255 = new TH2F("h2RatioSatQ255","h2RatioSatQ255",300,0,1200*pow(10,-6),20,0,1.05);
    TProfile* profSatQ255 = new TProfile("profSatQ255","profSatQ255",300,0,1200*pow(10,-6),"");

    /*vector<vector<TH2F*>> VectNStrip_VectNStripSat254_h2EvQ;
    vector<vector<TH1D*>> VectNStrip_VectNStripSat254_h1EvQ;
    vector<vector<TProfile*>> VectNStrip_VectNStripSat254_profEvQ;
    vector<vector<TH2F*>> VectNStrip_VectNStripSat255_h2EvQ;
    vector<vector<TH1D*>> VectNStrip_VectNStripSat255_h1EvQ;
    vector<vector<TProfile*>> VectNStrip_VectNStripSat255_profEvQ;
    for(int i=1;i<10;i++)
    {
        vector<TH2F*> VectNStripSat254_h2EvQ;
        vector<TH1D*> VectNStripSat254_h1EvQ;
        vector<TProfile*> VectNStripSat254_profEvQ;
        vector<TH2F*> VectNStripSat255_h2EvQ;
        vector<TH1D*> VectNStripSat255_h1EvQ;
        vector<TProfile*> VectNStripSat255_profEvQ;
        for(int j=0;j<i+1;j++)
        {
            string str_i = to_string(i);
            string str_j = to_string(j);
            TH2F* h2EvQVectVectSat254 = new TH2F(("H2 NStrip="+str_i+" & NStripSat254="+str_j+" | TOB1").c_str(),("H2 NStrip="+str_i+" & NStripSat254="+str_j+" | TOB1").c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6));
            TH1D* h1EvQVectVectSat254 = new TH1D(("H1 NStrip="+str_i+" & NStripSat254="+str_j+" | TOB1").c_str(),("H1 NStrip="+str_i+" & NStripSat254="+str_j+" | TOB1").c_str(),200,0,8);
            TProfile* profEvQVectVectSat254 = new TProfile(("prof NStrip="+str_i+" & NStripSat254="+str_j+" | TOB1").c_str(),("prof NStrip="+str_i+" & NStripSat254="+str_j+" | TOB1").c_str(),50,0,1500*pow(10,-6));
            VectNStripSat254_h2EvQ.push_back(h2EvQVectVectSat254);
            VectNStripSat254_h1EvQ.push_back(h1EvQVectVectSat254);
            VectNStripSat254_profEvQ.push_back(profEvQVectVectSat254);
            TH2F* h2EvQVectVectSat255 = new TH2F(("H2 NStrip="+str_i+" & NStripSat255="+str_j+" | TOB1").c_str(),("H2 NStrip="+str_i+" & NStripSat255="+str_j+" | TOB1").c_str(),300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
            TH1D* h1EvQVectVectSat255 = new TH1D(("H1 NStrip="+str_i+" & NStripSat255="+str_j+" | TOB1").c_str(),("H1 NStrip="+str_i+" & NStripSat255="+str_j+" | TOB1").c_str(),200,0,8);
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
    }*/
    
    vector<TH2F*> VectPartID_h2EvQ_NoSat;
    vector<TH2F*> VectPartID_h2EvQ_Sat;
    vector<TH2F*> VectPartID_h2EvQ_Sat254;
    vector<TH2F*> VectPartID_h2EvQ_Sat255;
    vector<TH1D*> VectPartID_Pt;
    vector<TH1D*> VectPartID_Eloss;
    vector<TH1D*> VectPartID_PoverM;


    for(int i=0;i<7;i++)
    {
        VectPartID_h2EvQ_NoSat.push_back(new TH2F(("NoSat_"+LoopPartID(i)).c_str(),("NoSat_"+LoopPartID(i)).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
        VectPartID_h2EvQ_Sat.push_back(new TH2F(("Sat_"+LoopPartID(i)).c_str(),("Sat_"+LoopPartID(i)).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
        VectPartID_h2EvQ_Sat254.push_back(new TH2F(("Sat254_"+LoopPartID(i)).c_str(),("Sat254_"+LoopPartID(i)).c_str(),300,0,1500*pow(10,-6),300,0,1500*pow(10,-6)));
        VectPartID_h2EvQ_Sat255.push_back(new TH2F(("Sat255_"+LoopPartID(i)).c_str(),("Sat255_"+LoopPartID(i)).c_str(),300,0,3000*pow(10,-6),300,0,3000*pow(10,-6)));
        VectPartID_Pt.push_back(new TH1D((LoopPartID(i)+"_p_{T}").c_str(),(LoopPartID(i)+"_p_{T}").c_str(),100,0,3000));
        VectPartID_Eloss.push_back(new TH1D((LoopPartID(i)+"_E_{sim}").c_str(),(LoopPartID(i)+"_E_{sim}").c_str(),300,0,3000*pow(10,-6)));
        VectPartID_PoverM.push_back(new TH1D((LoopPartID(i)+"_P/M").c_str(),(LoopPartID(i)+"_P/M").c_str(),1000,0.01,10));
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

    vector<TH1D*> VectLayer_h1EvQ_Tot;
    vector<TH1D*> VectLayer_h1EvQ_NoSat;
    vector<TH1D*> VectLayer_h1EvQ_Sat;
    vector<TH1D*> VectLayer_h1EvQ_Sat254;
    vector<TH1D*> VectLayer_h1EvQ_Sat255;

    for(int i=1;i<22;i++)
    {
        VectLayer_h2EvQ_Tot.push_back(new TH2F(("h2_"+Label(i)+"_tot").c_str(),("h2_"+Label(i)+"_tot").c_str(),300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        VectLayer_h1EvQ_Tot.push_back(new TH1D(("h1_"+Label(i)+"_tot").c_str(),("h1_"+Label(i)+"_tot").c_str(),200,0,8));
        VectLayer_profEvQ_Tot.push_back(new TProfile(("pr_"+Label(i)+"_tot").c_str(),("pr_"+Label(i)+"_tot").c_str(),50,0,900*pow(10,-6),""));

        VectLayer_h2EvQ_NoSat.push_back(new TH2F(("h2_"+Label(i)+"_NoSat").c_str(),("h2_"+Label(i)+"_NoSat").c_str(),300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        VectLayer_h1EvQ_NoSat.push_back(new TH1D(("h1_"+Label(i)+"_NoSat").c_str(),("h1_"+Label(i)+"_NoSat").c_str(),200,0,8));
        VectLayer_profEvQ_NoSat.push_back(new TProfile(("pr_"+Label(i)+"_NoSat").c_str(),("pr_"+Label(i)+"_NoSat").c_str(),50,0,900*pow(10,-6),""));

        VectLayer_h2EvQ_Sat.push_back(new TH2F(("h2_"+Label(i)+"_Sat").c_str(),("h2_"+Label(i)+"_Sat").c_str(),300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        VectLayer_h1EvQ_Sat.push_back(new TH1D(("h1_"+Label(i)+"_Sat").c_str(),("h1_"+Label(i)+"_Sat").c_str(),200,0,8));
        VectLayer_profEvQ_Sat.push_back(new TProfile(("pr_"+Label(i)+"_Sat").c_str(),("pr_"+Label(i)+"_Sat").c_str(),50,0,900*pow(10,-6),""));

        VectLayer_h2EvQ_Sat254.push_back(new TH2F(("h2_"+Label(i)+"_Sat254").c_str(),("h2_"+Label(i)+"_Sat254").c_str(),300,0,900*pow(10,-6),300,0,900*pow(10,-6)));
        VectLayer_h1EvQ_Sat254.push_back(new TH1D(("h1_"+Label(i)+"_Sat254").c_str(),("h1_"+Label(i)+"_Sat254").c_str(),200,0,8));
        VectLayer_profEvQ_Sat254.push_back(new TProfile(("pr_"+Label(i)+"_Sat254").c_str(),("pr_"+Label(i)+"_Sat254").c_str(),50,0,900*pow(10,-6),""));

        VectLayer_h2EvQ_Sat255.push_back(new TH2F(("h2_"+Label(i)+"_Sat255").c_str(),("h2_"+Label(i)+"_Sat255").c_str(),300,0,3000*pow(10,-6),300,0,3000*pow(10,-6)));
        VectLayer_h1EvQ_Sat255.push_back(new TH1D(("h1_"+Label(i)+"_Sat255").c_str(),("h1_"+Label(i)+"_Sat255").c_str(),200,0,8));
        VectLayer_profEvQ_Sat255.push_back(new TProfile(("pr_"+Label(i)+"_Sat255").c_str(),("pr_"+Label(i)+"_Sat255").c_str(),50,0,3000*pow(10,-6),""));
    }

    TH1D* hEmpty = new TH1D("hEmpty","hEmpty",21,0,21); //pour tracer les lignes avec la fonction SetHistoLabel

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
    TH1D* h1TestNoCutAndNoEdge = new TH1D("h1TestNoCutAndNoEdge","h1TestNoCutAndNoEdge",200,0,8);
    TH1D* h1TestCutOrEdge = new TH1D("h1TestCutOrEdge","h1TestCutOrEdge",200,0,8);

    TH1D* h1ChargeFauxSimHit = new TH1D("h1FauxSimHit","h1FauxSimHit",255,0,255);

    TH2F* EvQbarrel = new TH2F("EvQbarrel","EvQbarrel",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
    TH2F* EvQbarrel_corr = new TH2F("EvQbarrel_corr","EvQbarrel_corr",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));

    TH1D* h1DiffRelEvQ = new TH1D("h1DiffRelEvQ","h1DiffRelEvQ",1000,-1,2);
    TH1D* h1DiffRelEvQcorr = new TH1D("h1DiffRelEvQcorr","h1DiffRelEvQcorr",1000,-1,2);

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

    TH1D* hNStripPerCluster = new TH1D("hNStripPerCluster","hNStripPerCluster",10,0,10);
    TH1D* hStripMax = new TH1D("hStripMax","hStripMax",260,0,260);
    TH2F* h2StripMaxEoverQ = new TH2F("h2StripMaxEoverQ","h2StripMaxEoverQ",260,0,260,100,0,5);
    TH1D* hDeDxHarmo = new TH1D("hDeDxHarmo","hDeDxHarmo",1000,0,0.3);
    TH1D* hDeDxTrunc = new TH1D("hDeDxTrunc","hDeDxTrunc",1000,0,0.3);
    
    TH2F* h2PoverMDeDx = new TH2F("h2PoverMDeDx","h2PoverMDeDx",1000,0,2,1000,0,50);
    TH2F* h2PoverMDqDx = new TH2F("h2PoverMDqDx","h2PoverMDqDx",1000,0,2,1000,0,50);
    TH2F* h2PoverMDqcorrDx = new TH2F("h2PoverMDqcorrDx","h2PoverMDqcorrDx",1000,0,2,1000,0,50);
    
    TH2F* h2PoverMDeDxNoEstim = new TH2F("h2PoverMDeDxNoEstim","h2PoverMDeDxNoEstim",1000,0,2,1000,0,50);
    TH2F* h2PoverMDqDxNoEstim = new TH2F("h2PoverMDqDxNoEstim","h2PoverMDqDxNoEstim",1000,0,2,1000,0,50);
    TH2F* h2PoverMDqcorrDxNoEstim = new TH2F("h2PoverMDqcorrDxNoEstim","h2PoverMDqcorrDxNoEstim",1000,0,2,1000,0,50);

    TH2F* h2PoverMDeDxLarge = new TH2F("h2PoverMDeDxLarge","h2PoverMDeDxLarge",1000,0,10,1000,0,50);
    TH2F* h2PoverMDqDxLarge = new TH2F("h2PoverMDqDxLarge","h2PoverMDqDxLarge",1000,0,10,1000,0,50);
    TH2F* h2PoverMDqcorrDxLarge = new TH2F("h2PoverMDqcorrDxLarge","h2PoverMDqcorrDxLarge",1000,0,10,1000,0,50);
    
    TH2F* h2PoverMDeDxLargeNoEstim = new TH2F("h2PoverMDeDxLargeNoEstim","h2PoverMDeDxLargeNoEstim",1000,0,10,1000,0,50);
    TH2F* h2PoverMDqDxLargeNoEstim = new TH2F("h2PoverMDqDxLargeNoEstim","h2PoverMDqDxLargeNoEstim",1000,0,10,1000,0,50);
    TH2F* h2PoverMDqcorrDxLargeNoEstim = new TH2F("h2PoverMDqcorrDxLargeNoEstim","h2PoverMDqcorrDxLargeNoEstim",1000,0,10,1000,0,50);
    
    TH2F* h2PDeDx = new TH2F("h2PDeDx","h2PDeDx",1000,0,5,1000,0,15);
    TH2F* h2PDqDx = new TH2F("h2PDqDx","h2PDqDx",1000,0,5,1000,0,15);
    TH2F* h2PDqcorrDx = new TH2F("h2PDqcorrDx","h2PDqcorrDx",1000,0,5,1000,0,15);

    TH2F* h2PDeDxLarge = new TH2F("h2PDeDxLarge","h2PDeDxLarge",1000,0,3000,1000,0,100);
    TH2F* h2PDqDxLarge = new TH2F("h2PDqDxLarge","h2PDqDxLarge",1000,0,3000,1000,0,100);
    TH2F* h2PDqcorrDxLarge = new TH2F("h2PDqcorrDxLarge","h2PDqcorrDxLarge",1000,0,3000,1000,0,100);
    
    TH2F* h2PDeDxNoEstim = new TH2F("h2PDeDxNoEstim","h2PDeDxNoEstim",1000,0,5,1000,0,15);
    TH2F* h2PDqDxNoEstim = new TH2F("h2PDqDxNoEstim","h2PDqDxNoEstim",1000,0,5,1000,0,15);
    TH2F* h2PDqcorrDxNoEstim = new TH2F("h2PDqcorrDxNoEstim","h2PDqcorrDxNoEstim",1000,0,5,1000,0,15);
    
    TH2F* h2PDeDxLargeNoEstim = new TH2F("h2PDeDxLargeNoEstim","h2PDeDxLargeNoEstim",1000,0,3000,1000,0,100);
    TH2F* h2PDqDxLargeNoEstim = new TH2F("h2PDqDxLargeNoEstim","h2PDqDxLargeNoEstim",1000,0,3000,1000,0,100);
    TH2F* h2PDqcorrDxLargeNoEstim = new TH2F("h2PDqcorrDxLargeNoEstim","h2PDqcorrDxLargeNoEstim",1000,0,3000,1000,0,100);
    

    TH1D* hInvestigationEta = new TH1D("hInvestigationEta","hInvestigationEta",50,-5,5);
    TH1D* hInvestigationPhi = new TH1D("hInvestigationPhi","hInvestigationPhi",50,-5,5);
    TH1D* hInvestigationEloss = new TH1D("hInvestigationEloss","hInvestigationEloss",100,0,0.006);
    TH1D* hInvestigationCharge = new TH1D("hInvestigationCharge","hInvestigationCharge",100,0,0.006);
    TH2F* h2InvestigationEvQ = new TH2F("h2InvestigationEvQ","h2InvestigationEvQ",300,0,0.006,300,0,0.006);
    TH1D* hInvestigationDetId = new TH1D("hInvestigationDetId","hInvestigationDetId",1000,0,500000000);
    TH1D* hInvestigationEntries = new TH1D("hInvestigationEntries","hInvestigationEntries",1000,0,16000);
    TH1D* hInvestigationShape = new TH1D("hInvestigationShape","hInvestigationShape",2,0,2);
    TH1D* hInvestigationSimHits = new TH1D("hInvestigationSimHits","hInvestigationSimHits",5,0,5);
    TH1D* hInvestigationPt = new TH1D("hInvestigationPt","hInvestigationPt",100,0,3000);
    TH2F* h2Shape = new TH2F("h2Shape","h2Shape",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
    TH1D* hShapePartID = new TH1D("ShapePartID","ShapePartID",7,0,7);

    TH1D* hDistribMassSMeloss = new TH1D("hDistribMassSMeloss","hDistribMassSMeloss",50,0,2);
    TH1D* hDistribMassSUSYeloss = new TH1D("hDistribMassSUSYeloss","hDistribMassSUSYeloss",50,0,5000);
    TH1D* hDistribMassSMq = new TH1D("hDistribMassSMq","hDistribMassSMq",50,0,5000);
    TH1D* hDistribMassSUSYq = new TH1D("hDistribMassSUSYq","hDistribMassSUSYq",50,0,5000);
    TH1D* hDistribMassSMqcorr = new TH1D("hDistribMassSMqcorr","hDistribMassSMqcorr",50,0,5000);
    TH1D* hDistribMassSUSYqcorr = new TH1D("hDistribMassSUSYqcorr","hDistribMassSUSYqcorr",50,0,5000);
    TH1D* hDistribMassttbarqcorr = new TH1D("hDistribMassttbarqcorr","hDistribMassttbarqcorr",100,0,200);
    TH1D* hDistribMassSUSYqcorrinf = new TH1D("hDistribMassSUSYqcorrinf","hDistribMassSUSYqcorrinf",50,0,5000);
    TH1D* hDistribMassSUSYqcorrsup = new TH1D("hDistribMassSUSYqcorrsup","hDistribMassSUSYqcorrsup",50,0,5000);

    TH2F* h2PoverMDeDxProton = new TH2F("h2PoverMDeDxProton","h2PoverMDeDxProton",50,0,3,50,0,50);
    //TH2F* h2PoverMDeDxPion = new TH2F("h2PoverMDeDxPion","h2PoverMDeDxPion",1000,0,1000,50,0,7);
    TH2F* h2PoverMDeDxKaon = new TH2F("h2PoverMDeDxKaon","h2PoverMDeDxKaon",50,0,3,50,0,50);
    TH2F* h2PoverMDeDxRhad = new TH2F("h2PoverMDeDxRhad","h2PoverMDeDxRhad",50,0,3,50,0,50);
    TH2F* h2PoverMDeDxSUSY = new TH2F("h2PoverMDeDxSUSY","h2PoverMDeDxSUSY",50,0,3,50,0,50);

    TH1D* h1TestDeDxProtonFixMom = new TH1D("h1TestDeDxProtonFixMom","h1TestDeDxProtonFixMom",100,0,15);
    TH1D* h1TestDeDxProtonFixMomNoEstim = new TH1D("h1TestDeDxProtonFixMomNoEstim","h1TestDeDxProtonFixMomNoEstim",100,0,15);

    TH1D* hDistribMassDeuteroneloss = new TH1D("hDistribMassDeuteroneloss","hDistribMassDeuteroneloss",200,0,2);
    TH1D* hDistribMassProtoneloss = new TH1D("hDistribMassProtoneloss","hDistribMassProtoneloss",200,0,2);
    TH1D* hDistribMassKaoneloss = new TH1D("hDistribMassKaoneloss","hDistribMassKaoneloss",200,0,2);
    TH1D* hDistribMassPioneloss = new TH1D("hDistribMassPioneloss","hDistribMassPioneloss",200,0,2);
    TH1D* hDistribMassMSeloss = new TH1D("hDistribMassMSeloss","hDistribMassMSeloss",200,0,2);

    TH1D* elosspathTOB1 = new TH1D("elosspathTOB1","elosspathTOB1",100,0,15);
    TH1D* elosspathTIB1 = new TH1D("elosspathTIB1","elosspathTIB1",100,0,15);
    TH1D* elossTOB1 = new TH1D("elossTOB1","elossTOB1",100,0,0.5);
    TH1D* elossTIB1 = new TH1D("elossTIB1","elossTIB1",100,0,0.5);

    TH1D* ratsatpion = new TH1D("ratsatpion","ratsatpion",15,0,1.2);
    TH1D* ratsatkaon = new TH1D("ratsatkaon","ratsatkaon",15,0,1.2);
    TH1D* ratsatproton = new TH1D("ratsatproton","ratsatproton",15,0,1.2);
    TH1D* ratsatMS = new TH1D("ratsatMS","ratsatMS",15,0,1.2);
    TH1D* ratsatrhad = new TH1D("ratsatrhad","ratsatrhad",15,0,1.2);
    TH1D* ratsatSUSY = new TH1D("ratsatSUSY","ratsatSUSY",15,0,1.2);
    
    TH1D* ratsatpion255 = new TH1D("ratsatpion255","ratsatpion255",15,0,1.2);
    TH1D* ratsatkaon255 = new TH1D("ratsatkaon255","ratsatkaon255",15,0,1.2);
    TH1D* ratsatproton255 = new TH1D("ratsatproton255","ratsatproton255",15,0,1.2);
    TH1D* ratsatMS255 = new TH1D("ratsatMS255","ratsatMS255",15,0,1.2);
    TH1D* ratsatrhad255 = new TH1D("ratsatrhad255","ratsatrhad255",15,0,1.2);
    TH1D* ratsatSUSY255 = new TH1D("ratsatSUSY255","ratsatSUSY255",15,0,1.2);

    TH2D* MassVsIas = new TH2D("MassVsIas","MassVsIas",50,0,5000,10,0,1);
    TH1D* MassPt5 = new TH1D("MassPt5","MassPt5",50,0,5000);
    TH1D* MassPt20 = new TH1D("MassPt20","MassPt20",50,0,5000);
    TH1D* MassPt35 = new TH1D("MassPt35","MassPt35",50,0,5000);
    TH1D* MassPt50 = new TH1D("MassPt50","MassPt50",50,0,5000);
    TH1D* MassPtcorr5 = new TH1D("MassPtcorr5","MassPtcorr5",50,0,5000);
    TH1D* MassPtcorr20 = new TH1D("MassPtcorr20","MassPtcorr20",50,0,5000);
    TH1D* MassPtcorr35 = new TH1D("MassPtcorr35","MassPtcorr35",50,0,5000);
    TH1D* MassPtcorr50 = new TH1D("MassPtcorr50","MassPtcorr50",50,0,5000);
    TH1D* MassBeforeFactorWithCorrection = new TH1D("MassBeforeFactorWithCorrection","MassBeforeFactorWithCorrection",50,0,5000);
    TH1D* MassBeforeFactorWithoutCorrection = new TH1D("MassBeforeFactorWithoutCorrection","MassBeforeFactorWithoutCorrection",50,0,5000);


     
    int nbins = 72;
	double xbinsPion[nbins];
	/*for(int i=1;i<nbins+1;i++)
	{

        
		xbinsPion[i]=pow(10,-2);
	}*/

    int counterbin=0;
    for(int i=0;i<4;i++)
    {
        for(int j=2;j<20;j++)
        {
            cout<<counterbin<<"   "<<pow(10,-1+i)*(float)j/2<<endl;
            xbinsPion[counterbin]=pow(10,-1+i)*(float)j/2;
            counterbin++;
        }
    }

    //double xbinsPion[31]={0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,220,240,260,280,300,350,400,450,500,550,600,700,800,900,1000};

    

    
    TH2F* h2PvsPt = new TH2F("h2PvsPt","h2PvsPt",100,0,100,100,0,100);
    TProfile* profPvsPt = new TProfile("profPvsPt","profPvsPt",50,0,100,"");

    TH2F* h2PoverMDeDxPionHarmonic2 = new TH2F("h2PoverMDeDxPionHarmonic2","h2PoverMDeDxPionHarmonic2",nbins-1,xbinsPion,50,0,10);
    TProfile* profPoverMDeDxPionHarmonic2 = new TProfile("profPoverMDeDxPionHarmonic2","profPoverMDeDxPionHarmonic2",nbins-1,xbinsPion,"S");
    TH2F* h2PoverMDeDxPionHarmonic2Clone = (TH2F*) h2PoverMDeDxPionHarmonic2->Clone("pionHarmonic2clone");
    h2PoverMDeDxPionHarmonic2Clone->Reset();

    TH2F* h2PoverMDeDxPionTrunc40 = new TH2F("h2PoverMDeDxPionTrunc40","h2PoverMDeDxPionTrunc40",nbins-1,xbinsPion,50,0,10);
    TProfile* profPoverMDeDxPionTrunc40 = new TProfile("profPoverMDeDxPionTrunc40","profPoverMDeDxPionTrunc40",nbins-1,xbinsPion,"S");
    TH2F* h2PoverMDeDxPionTrunc40Clone = (TH2F*) h2PoverMDeDxPionTrunc40->Clone("pionTrunc40clone");
    h2PoverMDeDxPionTrunc40Clone->Reset();

    TH2F* h2PoverMDeDxPionRMS = new TH2F("h2PoverMDeDxPionRMS","h2PoverMDeDxPionRMS",nbins-1,xbinsPion,50,0,10);
    TProfile* profPoverMDeDxPionRMS = new TProfile("profPoverMDeDxPionRMS","profPoverMDeDxPionRMS",nbins-1,xbinsPion,"S");
    TH2F* h2PoverMDeDxPionRMSClone = (TH2F*) h2PoverMDeDxPionRMS->Clone("pionRMSclone");
    h2PoverMDeDxPionRMSClone->Reset();

    TH2F* h2PoverMvsEtaPion = new TH2F("h2PoverMvsEtaPion","h2PoverMvsEtaPion",nbins-1,xbinsPion,50,0,5);
    TProfile* profPoverMvsEtaPion = new TProfile("profPoverMvsEtaPion","profPoverMvsEtaPion",nbins-1,xbinsPion,"");
    TH2F* h2PoverMvsPathlengthPion = new TH2F("h2PoverMvsPathlengthPion","h2PoverMvsPathlengthPion",nbins-1,xbinsPion,50,0,0.4);
    TProfile* profPoverMvsPathlengthPion = new TProfile("profPoverMvsPathlengthPion","profPoverMvsPathlengthPion",nbins-1,xbinsPion,"");
    TH2F* h2PoverMDeDxPionCluster = new TH2F("h2PoverMDeDxPionCluster","h2PoverMDeDxPionCluster",nbins-1,xbinsPion,50,0,10);
    TProfile* profPoverMDeDxPionCluster = new TProfile("profPoverMDeDxPionCluster","profPoverMDeDxPionCluster",nbins-1,xbinsPion,"S");
    TH2F* h2PoverMDeDxPionClusterClone = (TH2F*) h2PoverMDeDxPionCluster->Clone("pionclusterclone");
    h2PoverMDeDxPionClusterClone->Reset();

    TH3F* h3plengthEtaPt = new TH3F("h3plengthEtaP","h3plengthEtaP",50,0,100,50,0,3,50,0,0.4);


    TH2F* h2testmodulgeom = new TH2F("h2testmodulgeom","h2testmodulgeom",20,0,20,20,0,20);

   

    TF1* fKC = new TF1("fKC","[0]+[1]*(1/pow(x,2))");


    vector<TH1D*> Vectelosspath;
    vector<string> VectelosspathLegend;
    for(int i=5;i<10;i++)
    {
        Vectelosspath.push_back(new TH1D("","",100,0,15));
        VectelosspathLegend.push_back(("Taille : "+to_string(i)).c_str());
    }


    int entries = nentries;


    float K=0.;
    float C=0.;

    bool testsat254 = false; 
    bool testsat255 = false;
    bool testshape  = true;

    int clust=0,clustsat254=0,clustsat255=0;

    b1->SetThresholdPartId(0.7); //on veut un minimum de 70% de simhits identiques
    b1->SetThresholdPt(0); //on veut des traces avec un minimum de 60 GeV en pt
    b1->SetThresholdP(0);
    b1->SetThresholdEta(5);

    Correction Charge;
    TFile* file2 = TFile::Open(s3.c_str());
    TTree* tree2 = (TTree*) file2->Get("tree");
    Charge.SetTree(*tree2);
    int count=0;
    bool testpt=false;

    string comment = "p>="+to_string(b1->GetThresholdP())+"GeV pt>="+to_string(b1->GetThresholdPt())+"GeV";

    for(int i=0;i<entries;i++)
    {
        if(i%1000==0) cout<<"Event "<<i<<endl;
        b1->GetEntry(i);
        for(int track=0;track<b1->GetNtracks();track++)
        {
            testpt=true;
            vector<int> vect_partID;
            vector<float> vect_eloss;
            vector<float> vect_charge;
            vector<float> vect_dedx;
            vector<float> vect_dqdx;
            vector<float> vect_dqcorrdx;
            vector<float> vect_dqcorrdx2;
            vector<float> vect_plength;
            testsat254              = false;
            testsat255              = false;
            testshape               = true;
            float pt                = b1->GetVectTrack()[track].GetPt();
            float p                 = b1->GetVectTrack()[track].GetP();
            float eta               = b1->GetVectTrack()[track].GetEta();
            float phi               = b1->GetVectTrack()[track].GetPhi();
            float ias_ampl          = b1->GetVectTrack()[track].GetIasAmpl();
            int NCluster            = b1->GetVectTrack()[track].GetNCluster();
            int NClustSat254        = b1->GetVectTrack()[track].GetNSatCluster(254);//+b1->GetVectTrack()[track].GetNSatCluster(255);
            int NClustSat255        = b1->GetVectTrack()[track].GetNSatCluster(255);
            float RatioNClusterSat254 = (double)NClustSat254/(double)NCluster;
            float RatioNClusterSat255 = (double)NClustSat255/(double)NCluster;
            int id                  = b1->GetVectTrack()[track].GetPartId();
            float PoverM            = GetPoverM(p,id);

            h2PvsPt->Fill(pt,p);
            profPvsPt->Fill(pt,p);

            //if(((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001))) {C=CsusyRec; K=KsusyRec;}
            //if(!((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001))) {C=CmsRec; K=KmsRec;}
            //C=CmsRec;K=KmsRec;


//if(p>=1)
//if(/*((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001)) && */(RatioNClusterSat254>0 || RatioNClusterSat255>0))
//if(id==211 || id==-211)
//if(((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001)) || id==211 || id==-211)

//if(((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001)))
{
    count++;

            if(id==11 || id==-11) 
            {
                VectPartID_Pt[0]->Fill(pt);
                VectPartID_PoverM[0]->Fill(GetPoverM(p,id));
            }
            if(id==211 || id==-211) 
            {
                VectPartID_Pt[1]->Fill(pt);
                VectPartID_PoverM[1]->Fill(GetPoverM(p,id));
                ratsatpion->Fill(RatioNClusterSat254);
                ratsatpion255->Fill(RatioNClusterSat255);
            }
            if(id==321 || id==-321) 
            {
                VectPartID_Pt[2]->Fill(pt);
                VectPartID_PoverM[2]->Fill(GetPoverM(p,id));
                ratsatkaon->Fill(RatioNClusterSat254);
                ratsatkaon255->Fill(RatioNClusterSat255);
            }
            if(id==2212 || id==-2212) 
            {
                VectPartID_Pt[3]->Fill(pt);
                VectPartID_PoverM[3]->Fill(GetPoverM(p,id));
                ratsatproton->Fill(RatioNClusterSat254);
                ratsatproton255->Fill(RatioNClusterSat255);
            }    
            if(id==1009213 || id==-1009213) 
            {
                VectPartID_Pt[5]->Fill(pt);
                VectPartID_PoverM[5]->Fill(GetPoverM(p,id));
                ratsatrhad->Fill(RatioNClusterSat254);
                ratsatrhad255->Fill(RatioNClusterSat255);
            }
            else if((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001))
            {
                VectPartID_Pt[6]->Fill(pt);
            }
            else 
            {
                VectPartID_Pt[4]->Fill(pt);
            }
            if((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001)) 
            {
                ratsatSUSY->Fill(RatioNClusterSat254);
                ratsatSUSY255->Fill(RatioNClusterSat255);
            }
            if(!((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001))) 
            {
                ratsatMS->Fill(RatioNClusterSat254);
                ratsatMS255->Fill(RatioNClusterSat255);
            }
            if(id==1000010020 || id==-1000010020)
            {

            }
            
            
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
            {//cout<<"Entry "<<i<<" track "<<track<<" cluster "<<cluster<<endl;
                float charge        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSclusCharge();
                float eloss         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetEloss();
                float pathlength    = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetPathLength();
                int layer           = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetLayer();
                int layerLabel      = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetLayerLabel();
                bool sat254         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSat254();
                bool sat255         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSat255();
                bool shape          = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetShape();
                int nsimhits        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSimHits();
                int nstrips         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNStrip();
                //int nsatboth        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStripBoth();
                int nsat254         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStrip(254);
                int nsat255         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStrip(255);
                bool edge           = b1->GetVectTrack()[track].GetVectClusters()[cluster].Edge();
                bool cut            = b1->GetVectTrack()[track].GetVectClusters()[cluster].Cut();
                int detid           = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetDetId();
                int subdetid        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSubDetId();
                int modulgeom      = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetModulGeom();
                int maxstrip        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetMaxStrip();
                //bool test_ilot      = SelectedArea()
                float Chargecorr   = charge;








                h2testmodulgeom->Fill(modulgeom,detid);

                

                vect_plength.push_back(pathlength);

                clust++;

                /*vect_eloss.push_back(eloss);
                vect_charge.push_back(charge);
                vect_dedx.push_back(eloss/pathlength);
                vect_dqdx.push_back(charge/pathlength);

                hNStripPerCluster->Fill(nstrips);
                hStripMax->Fill(maxstrip);
                h2StripMaxEoverQ->Fill(maxstrip,eloss/charge);*/

                if(layerLabel==1) {elosspathTIB1->Fill(eloss*pow(10,3)/pathlength);elossTIB1->Fill(eloss*pow(10,3));}
                if(layerLabel==5) {elosspathTOB1->Fill(eloss*pow(10,3)/pathlength);elossTOB1->Fill(eloss*pow(10,3));}

                for(int taille=5;taille<10;taille++)
                {
                    if(taille==NCluster) Vectelosspath[taille-5]->Fill(eloss*pow(10,3));
                }

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
                else if((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001))
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

                if(/*(nsat254==0 && nsat255==0) &&*/ layerLabel<=10)// && ((RatioNClusterSat254>=ratiosat || RatioNClusterSat255>=ratiosat) || (RatioNClusterSat254<ratiosat && RatioNClusterSat255<ratiosat && sat254==false && sat255==false)))//(nsat254==0 && nsat255==0)) )/*&& RatioNClusterSat254>=ratiosat*//* && (1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001)*/)
                {
                    if((RatioNClusterSat254>=ratiosat || RatioNClusterSat255>=ratiosat) )//|| ((RatioNClusterSat254<ratiosat && RatioNClusterSat255<ratiosat && sat254==false && sat255==false)))
                    {
                        //Chargecorr=Charge.ChargeCorr(charge,layerLabel,nstrips,nsat254,nsat255);
                    }
                    /* else
                    {
                        if(!(sat254==false && sat255==false)) Chargecorr=0.;
                        else Chargecorr=charge;
                    }*/

                    //if(sat254==false && sat255==false) 
                    vect_dqdx.push_back(charge/pathlength);

                    if(shape==0)
                    {
                        h2Shape->Fill(eloss,charge);
                        testshape=false;
                    } 
                    EvQbarrel_corr->Fill(eloss,Chargecorr);
                    EvQbarrel->Fill(eloss,charge);
                    h1DiffRelEvQ->Fill((charge-eloss)/eloss);  
                    h1DiffRelEvQcorr->Fill((Chargecorr-eloss)/eloss);
                    vect_eloss.push_back(eloss);
                    vect_charge.push_back(charge);
                    
                    
                    vect_dqcorrdx.push_back(Chargecorr/pathlength);
    
                    hNStripPerCluster->Fill(nstrips);
                    hStripMax->Fill(maxstrip);
                    h2StripMaxEoverQ->Fill(maxstrip,eloss/charge);


                    /*if(charge<=0.00085)
                    {
                        hInvestigationCharge->Fill(charge);
                        hInvestigationEloss->Fill(eloss);
                        hInvestigationEta->Fill(eta);
                        hInvestigationPhi->Fill(phi);
                        hInvestigationDetId->Fill(detid);
                        h2InvestigationEvQ->Fill(eloss,charge);
                        hInvestigationEntries->Fill(i);
                        hInvestigationShape->Fill(shape);
                        hInvestigationSimHits->Fill(nsimhits);
                        cout<<"particle : "<<LabelParticle(id)<<endl;
                        cout<<"detid="<<detid<<" subdetid="<<subdetid<<endl;
                        cout<<"eta="<<eta<<" phi="<<phi<<endl;
                        cout<<"eloss="<<eloss/(3.61*pow(10,-9)*247)<<" q="<<charge<<endl;
                        TCanvas* c_profClust = new TCanvas();
                        TProfile* profDistribStrip = &b1->GetVectTrack()[track].GetVectClusters()[cluster].GetDistribStrip();
                        profDistribStrip->SetLineColor(kBlue);
                        profDistribStrip->Draw();
                        c_profClust->SaveAs("./data/ProfileClust.pdf");
                        getchar();
                    }*/
                }
                if(layerLabel<=10)
                {
                    
                    vect_dedx.push_back(eloss/pathlength);
                } 
                
                //vect_dqcorrdx.push_back(Chargecorr/pathlength);
                //if(nsimhits>=2 && (id==11 || id==-11)) DrawClusterProfile(&b1->GetVectTrack()[track].GetVectClusters()[cluster]);

                   
                /*for(int nstrip=1;nstrip<10;nstrip++)
                {
                    for(int nstripsat=0;nstripsat<nstrip+1;nstripsat++)
                    {
                        if(layerLabel==5) //on se place dans TOB1 pour ne pas melanger differents effets + selection sans les bords ou les cuts
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
                }*/



                for(int countlayer=1;countlayer<11;countlayer++)
                {
                    for(int countnstrip=3;countnstrip<7;countnstrip++)
                    {
                        for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
                        {
                            for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
                            {
                                if(countlayer==layerLabel && countnstrip==nstrips && countnstripsat254==nsat254 && countnstripsat255==nsat255)
                                {
                                    if(layerLabel<=4) VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(eloss,charge);
                                    if(layerLabel>=5 && layerLabel<=10 && shape ) VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(eloss,charge);
                                } 
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

            estim.SetVect(vect_charge);
            h2RatioSatQ254->Fill(estim.GetMean(),RatioNClusterSat254);
            profSatQ254->Fill(estim.GetMean(),RatioNClusterSat254);
            h2RatioSatQ255->Fill(estim.GetMean(),RatioNClusterSat255);
            profSatQ255->Fill(estim.GetMean(),RatioNClusterSat255);
            
            estim.SetVect(vect_dedx);
            hDeDxHarmo->Fill(estim.GetHarmonic2());
            hDeDxTrunc->Fill(estim.GetTrunc40());

            h2PoverMDeDx->Fill(PoverM,estim.GetHarmonic2()*pow(10,3));
            h2PoverMDeDxLarge->Fill(PoverM,estim.GetHarmonic2()*pow(10,3));
            h2PDeDx->Fill(p,estim.GetHarmonic2()*pow(10,3));
            h2PDeDxLarge->Fill(p,estim.GetHarmonic2()*pow(10,3));
            //if(estim.GetHarmonic2()>0)
            {
                //estim.SetVect(vect_dqcorrdx);
                //if(id==211 || id==-211) 
                {
                    if(estim.GetHarmonic2()>0)
                    {
                        h2PoverMDeDxPionHarmonic2->Fill(PoverM,estim.GetHarmonic2()*pow(10,3));
                        profPoverMDeDxPionHarmonic2->Fill(PoverM,estim.GetHarmonic2()*pow(10,3));
                    }

                    if(estim.GetTrunc40()>0)
                    {
                        h2PoverMDeDxPionTrunc40->Fill(PoverM,estim.GetTrunc40()*pow(10,3));
                        profPoverMDeDxPionTrunc40->Fill(PoverM,estim.GetTrunc40()*pow(10,3));
                    }

                    if(estim.GetRMS()>0)
                    {
                        h2PoverMDeDxPionRMS->Fill(PoverM,estim.GetRMS()*pow(10,3));
                        profPoverMDeDxPionRMS->Fill(PoverM,estim.GetRMS()*pow(10,3));
                    }
                }
                if(id==321 || id==-321) h2PoverMDeDxKaon->Fill(PoverM,estim.GetHarmonic2()*pow(10,3));
                if(id==2212 || id==-2212) h2PoverMDeDxProton->Fill(PoverM,estim.GetHarmonic2()*pow(10,3));
                if(id==1009213 || id==-1009213) h2PoverMDeDxRhad->Fill(PoverM,estim.GetHarmonic2()*pow(10,3));
                if((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001)) h2PoverMDeDxSUSY->Fill(PoverM,estim.GetHarmonic2()*pow(10,3));

            }



            
            h2PoverMvsEtaPion->Fill(PoverM,abs(eta));
            profPoverMvsEtaPion->Fill(PoverM,abs(eta));

            

            for(int plength=0;plength<vect_plength.size();plength++)
            {
                h2PoverMvsPathlengthPion->Fill(PoverM,vect_plength[plength]);
                profPoverMvsPathlengthPion->Fill(PoverM,vect_plength[plength]);
                h3plengthEtaPt->Fill(pt,abs(eta),vect_plength[plength]);
            } 
            
            for(int nclust=0;nclust<vect_dedx.size();nclust++) 
            {
                h2PoverMDeDxPionCluster->Fill(PoverM,vect_dedx[nclust]*pow(10,3));
                profPoverMDeDxPionCluster->Fill(PoverM,vect_dedx[nclust]*pow(10,3));
            }









            
            if(p>=0.75 && p<=0.85)
            {
                h1TestDeDxProtonFixMom->Fill(estim.GetHarmonic2()*pow(10,3));
                for(int nclust=0;nclust<vect_dedx.size();nclust++) h1TestDeDxProtonFixMomNoEstim->Fill(vect_dedx[nclust]*pow(10,3));
            }

            for(int nclust=0;nclust<vect_dedx.size();nclust++) 
            {
                h2PoverMDeDxNoEstim->Fill(PoverM,vect_dedx[nclust]*pow(10,3));
                h2PoverMDeDxLargeNoEstim->Fill(PoverM,vect_dedx[nclust]*pow(10,3));
                h2PDeDxNoEstim->Fill(p,vect_dedx[nclust]*pow(10,3));
                h2PDeDxLargeNoEstim->Fill(p,vect_dedx[nclust]*pow(10,3));
            }

            estim.SetVect(vect_dqdx);
            h2PoverMDqDx->Fill(PoverM,estim.GetHarmonic2()*pow(10,3));
            h2PoverMDqDxLarge->Fill(PoverM,estim.GetHarmonic2()*pow(10,3));
            h2PDqDx->Fill(p,estim.GetHarmonic2()*pow(10,3));
            h2PDqDxLarge->Fill(p,estim.GetHarmonic2()*pow(10,3));

            for(int nclust=0;nclust<vect_dqdx.size();nclust++) 
            {
                h2PoverMDqDxNoEstim->Fill(PoverM,vect_dqdx[nclust]*pow(10,3));
                h2PoverMDqDxLargeNoEstim->Fill(PoverM,vect_dqdx[nclust]*pow(10,3));
                h2PDqDxNoEstim->Fill(p,vect_dqdx[nclust]*pow(10,3));
                h2PDqDxLargeNoEstim->Fill(p,vect_dqdx[nclust]*pow(10,3));
            }

            estim.SetVect(vect_dqcorrdx);
            h2PoverMDqcorrDx->Fill(PoverM,estim.GetHarmonic2()*pow(10,3));
            h2PoverMDqcorrDxLarge->Fill(PoverM,estim.GetHarmonic2()*pow(10,3));
            h2PDqcorrDx->Fill(p,estim.GetHarmonic2()*pow(10,3));
            h2PDqcorrDxLarge->Fill(p,estim.GetHarmonic2()*pow(10,3));
            for(int nclust=0;nclust<vect_dqcorrdx.size();nclust++) 
            {
                h2PoverMDqcorrDxNoEstim->Fill(PoverM,vect_dqcorrdx[nclust]*pow(10,3));
                h2PoverMDqcorrDxLargeNoEstim->Fill(PoverM,vect_dqcorrdx[nclust]*pow(10,3));
                h2PDqcorrDxNoEstim->Fill(p,vect_dqcorrdx[nclust]*pow(10,3));
                h2PDqcorrDxLargeNoEstim->Fill(p,vect_dqcorrdx[nclust]*pow(10,3));
            }
            estim.SetVect(vect_dedx);
            float dedx = estim.GetHarmonic2()*pow(10,3);
            if(((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001))) {K=KsusySim;C=CsusySim;}
            if(!((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001))) {K=KmsSim;C=CmsSim;}
            float mass_dedx = sqrt((dedx-C)*pow(p,2)/K);
            if(dedx<=C && dedx>=0) mass_dedx=0;
            estim.SetVect(vect_dqdx);
            float dqdx = estim.GetHarmonic2()*pow(10,3);
            K=KmsSim;C=CmsSim;
            float mass_dqdx = sqrt((dqdx-C)*pow(p,2)/K);
            estim.SetVect(vect_dqcorrdx);
            float dqcorrdx = estim.GetHarmonic2()*pow(10,3);
            float mass_dqcorrdx = sqrt((dqcorrdx-C)*pow(p,2)/K);

            
            if(dqdx<=C && dqdx>=0) mass_dqdx=0;
            if(dqcorrdx<=C && dqcorrdx>=0) mass_dqcorrdx=0;
           

            //if(mass_dqcorrdx>0)
            {

            
                if(!((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001)))
                {
                    hDistribMassSMeloss->Fill(mass_dedx);
                    hDistribMassSMq->Fill(mass_dqdx);
                    hDistribMassSMqcorr->Fill(mass_dqcorrdx);
                    hDistribMassttbarqcorr->Fill(mass_dqcorrdx);
                    
                    MassVsIas->Fill(mass_dqcorrdx,ias_ampl);
                    if(mass_dqdx>=1150 && mass_dqdx<=2300)
                    {
                        MassBeforeFactorWithoutCorrection->Fill(mass_dqdx);
                    }
                    if(mass_dqcorrdx>=1550 && mass_dqcorrdx<=2650)
                    {
                        MassBeforeFactorWithCorrection->Fill(mass_dqcorrdx);
                    }
                    if(pt>=5)
                    {
                        MassPt5->Fill(mass_dqdx,factor_bkg_before);
                        MassPtcorr5->Fill(mass_dqcorrdx,factor_bkg);
                    }
                    if(pt>=20)
                    {
                        MassPt20->Fill(mass_dqdx,factor_bkg_before);
                        MassPtcorr20->Fill(mass_dqcorrdx,factor_bkg);
                    }
                    if(pt>=35)
                    {
                        MassPt35->Fill(mass_dqdx,factor_bkg_before);
                        MassPtcorr35->Fill(mass_dqcorrdx,factor_bkg);
                    }
                    if(pt>=50)
                    {
                        MassPt50->Fill(mass_dqdx,factor_bkg_before);
                        MassPtcorr50->Fill(mass_dqcorrdx,factor_bkg);
                    }
                }
                
                if(((1000001<=id && id<=2000015) || (-2000015<=id && id<=-1000001)))
                {
                    hDistribMassSUSYeloss->Fill(mass_dedx);
                    hDistribMassSUSYq->Fill(mass_dqdx);
                    if((RatioNClusterSat254>=ratiosat || RatioNClusterSat255>=ratiosat))// || ((RatioNClusterSat254==0 && RatioNClusterSat255==0)))
                    {
                        hDistribMassSUSYqcorr->Fill(mass_dqcorrdx);//,factor_signal);
                    }
                    else
                    {
                        hDistribMassSUSYqcorr->Fill(mass_dqdx);
                    }
                    if((RatioNClusterSat254>=ratiosat || RatioNClusterSat255>=ratiosat))// || ((RatioNClusterSat254==0 && RatioNClusterSat255==0)))
                    {
                        hDistribMassSUSYqcorrsup->Fill(mass_dqcorrdx);
                    }
                    else if(RatioNClusterSat254==0 && RatioNClusterSat255==0)
                    {
                        hDistribMassSUSYqcorrsup->Fill(mass_dqdx);
                    }
                    else
                    {
                        hDistribMassSUSYqcorrinf->Fill(mass_dqcorrdx);
                    }
                }
                /*if(!(dqcorrdx>C))
                {
                    cerr<<"size clust "<<b1->GetVectTrack()[track].GetNCluster()<<endl;
                    cerr<<"id "<<id<<endl;
                    cerr<<"size qcorr "<<vect_dqcorrdx.size()<<endl;
                    cerr<<"size e "<<vect_dedx.size()<<endl;
                    cerr<<"dqdx "<<dqcorrdx<<endl;
                    cerr<<"p2/K "<<pow(p,2)/K<<endl;
                }*/
            }
            if(mass_dedx>0)
            {

            
            if(id==2212 || id==-2212) hDistribMassProtoneloss->Fill(mass_dedx);
            if(id==211 || id==-211) hDistribMassPioneloss->Fill(mass_dedx);
            if(id==321 || id==-321) hDistribMassKaoneloss->Fill(mass_dedx);
            if(id==1000010020 || id==-1000010020) hDistribMassDeuteroneloss->Fill(mass_dedx);
            else hDistribMassMSeloss->Fill(mass_dedx);
            }
            
            


            /*if(testsat254 && testsat255 )//&& b1->GetVectTrack()[track].GetVectClusters().size() <=3)
            {
                TCanvas* c_profClust = new TCanvas();
                c_profClust = &b1->GetVectTrack()[track].GetProfCluster();
                string str_pt = to_string(pt);
                TText* text_pt = new TText(c_profClust->GetWw()-100,c_profClust->GetWw()-100,str_pt.c_str());
                text_pt->Draw("SAME");
                c_profClust->SaveAs("ProfileClust.pdf");
                getchar();
            }*/
}
        }
        //if(testpt) count++;
    }
    cout<<entries<<endl;
    cout<<count<<endl;
    
    ofstream ofile("parameters.txt");
    ofile<<"Avant corr"<<endl;
    ofile<<"Mean : "<<h1DiffRelEvQ->GetMean()<<"\tError : "<<h1DiffRelEvQ->GetMeanError()<<endl;
    ofile<<"StdDev : "<<h1DiffRelEvQ->GetStdDev()<<"\tError : "<<h1DiffRelEvQ->GetStdDevError()<<endl;
    ofile<<"Apres corr"<<endl;
    ofile<<"Mean : "<<h1DiffRelEvQcorr->GetMean()<<"\tError : "<<h1DiffRelEvQcorr->GetMeanError()<<endl;
    ofile<<"StdDev : "<<h1DiffRelEvQcorr->GetStdDev()<<"\tError : "<<h1DiffRelEvQcorr->GetStdDevError()<<endl;


    TFile* ofilecorr = new TFile(("results_ratiosat"+to_string(ratiosat)+".root").c_str(),"RECREATE");




    DrawHisto(*ofilecorr,hDistribMassSMeloss,"","m [GeV/c^2]");
    DrawHisto(*ofilecorr,hDistribMassSUSYeloss,"","m [GeV/c^2]");
    DrawHisto(*ofilecorr,hDistribMassSMq,"","m [GeV/c^2]");
    DrawHisto(*ofilecorr,hDistribMassSUSYq,"","m [GeV/c^2]");
    DrawHisto(*ofilecorr,hDistribMassSMqcorr,"","m [GeV/c^2]");
    DrawHisto(*ofilecorr,hDistribMassSUSYqcorr,"","m [GeV/c^2]");
    DrawHisto(*ofilecorr,hDistribMassttbarqcorr,"","m [GeV/c^2]");
    DrawHisto(*ofilecorr,hDistribMassSUSYqcorrsup,"","m [GeV/c^2]");
    DrawHisto(*ofilecorr,hDistribMassSUSYqcorrinf,"","m [GeV/c^2]");
    //cDiffRelChargecorr_Charge->Write();
    h1DiffRelEvQ->Write();
    h1DiffRelEvQcorr->Write();
    ofilecorr->Write();
    ofilecorr->Close();
    delete ofilecorr;






    TFile* fileout = new TFile(s2.c_str(),"RECREATE");

    TF1* Log = new TF1("Log","[0]*TMath::Log(x)+[1]",1,1000);
    TF1* fitPvsPt = new TF1("fitPvsPt","pow(x,[0])+[1]");


    h3plengthEtaPt->Write();

    h2testmodulgeom->Write();


    //TFitResultPtr FitlogHarmonic2 = profPoverMDeDxPionHarmonic2->Fit("Log","S");
    //TFitResultPtr FitlogTrunc40 = profPoverMDeDxPionTrunc40->Fit("Log","S");
    //TFitResultPtr FitlogRMS = profPoverMDeDxPionRMS->Fit("Log","S");



    TFitResultPtr FitPvsPt = profPvsPt->Fit("fitPvsPt","S");
    TFitResultPtr FitlogDeDxCluster = profPoverMDeDxPionCluster->Fit("Log","S");


    TCanvas* csuperposed = new TCanvas();

    h2PoverMvsEtaPion->Write();
    profPoverMvsEtaPion->Write();
    h2PoverMvsPathlengthPion->Write();
    profPoverMvsPathlengthPion->Write();
    h2PoverMDeDxPionCluster->Write();
    profPoverMDeDxPionCluster->Write();

    SuperposedHisto2DProfile(*csuperposed,h2PoverMvsEtaPion,profPoverMvsEtaPion,"","p/M","#eta","[a.u.]");
    csuperposed->Write();
    SuperposedHisto2DProfile(*csuperposed,h2PoverMvsPathlengthPion,profPoverMvsPathlengthPion,"","p/M","Pathlength [cm]","[a.u.]");
    csuperposed->Write();
    SuperposedHisto2DProfile(*csuperposed,h2PoverMDeDxPionCluster,profPoverMDeDxPionCluster,"","p/M","dE/dx_{cluster} [MeV/cm]","[a.u.]");
    csuperposed->Write();

    for(int index_xbin=0;index_xbin<h2PoverMDeDxPionCluster->GetNbinsX()+1;index_xbin++)
    {
        for(int index_ybin=0;index_ybin<h2PoverMDeDxPionCluster->GetNbinsY()+1;index_ybin++)
        {
            if(h2PoverMDeDxPionCluster->ProjectionY()->GetBinLowEdge(index_ybin)>=profPoverMDeDxPionCluster->GetBinContent(index_xbin)+0.7) 
            {
                h2PoverMDeDxPionClusterClone->SetBinContent(index_xbin,index_ybin,h2PoverMDeDxPionCluster->GetBinContent(index_xbin,index_ybin));
            }
        }
    }
    DrawHisto(*fileout,h2PoverMDeDxPionClusterClone,"","p/M","dE/dx_{cluster} [MeV/cm]");

    h2PoverMDeDxPionClusterClone->ProjectionX()->Write();
    h2PoverMDeDxPionCluster->ProjectionX()->Write();

    TEfficiency* eff_pionCluster = new TEfficiency(*h2PoverMDeDxPionClusterClone->ProjectionX(),*h2PoverMDeDxPionCluster->ProjectionX());
 
    eff_pionCluster->Write();



    profPvsPt->Write();

    //profPoverMDeDxPion->Write();

    

    DrawHisto(*fileout,h2PvsPt,"","p_{T} [GeV/c]","p [GeV/c]");
    DrawHisto(*fileout,h2PoverMDeDxPionHarmonic2,"","p/M","dE/dx_{track} [MeV/cm]");

    
    SuperposedHisto2DProfile(*csuperposed,h2PvsPt,profPvsPt,"","p_{T} [GeV/c]","p [GeV/c]","[u.a.]");
    csuperposed->Write();
    csuperposed->Clear();
    SuperposedHisto2DProfile(*csuperposed,h2PoverMDeDxPionHarmonic2,profPoverMDeDxPionHarmonic2,"","p/M","dE/dx_{track} [MeV/cm]","[u.a.]");
    csuperposed->Write();


    for(int index_xbin=0;index_xbin<h2PoverMDeDxPionHarmonic2->GetNbinsX()+1;index_xbin++)
    {
        for(int index_ybin=0;index_ybin<h2PoverMDeDxPionHarmonic2->GetNbinsY()+1;index_ybin++)
        {
            if(h2PoverMDeDxPionHarmonic2->ProjectionY()->GetBinLowEdge(index_ybin)>=profPoverMDeDxPionHarmonic2->GetBinContent(index_xbin)+0.7) 
            {
                h2PoverMDeDxPionHarmonic2Clone->SetBinContent(index_xbin,index_ybin,h2PoverMDeDxPionHarmonic2->GetBinContent(index_xbin,index_ybin));
            }
        }
    }
    DrawHisto(*fileout,h2PoverMDeDxPionHarmonic2Clone,"","p/M","dE/dx_{track} [MeV/cm]");

    h2PoverMDeDxPionHarmonic2Clone->ProjectionX()->Write();
    h2PoverMDeDxPionHarmonic2->ProjectionX()->Write();

    TEfficiency* eff_pionHarmonic2 = new TEfficiency(*h2PoverMDeDxPionHarmonic2Clone->ProjectionX(),*h2PoverMDeDxPionHarmonic2->ProjectionX());
 
    eff_pionHarmonic2->Write();



    profPoverMDeDxPionTrunc40->Write();

    SuperposedHisto2DProfile(*csuperposed,h2PoverMDeDxPionTrunc40,profPoverMDeDxPionTrunc40,"","p/M","dE/dx_{track} [MeV/cm]","[u.a.]");
    csuperposed->Write();


    for(int index_xbin=0;index_xbin<h2PoverMDeDxPionTrunc40->GetNbinsX()+1;index_xbin++)
    {
        for(int index_ybin=0;index_ybin<h2PoverMDeDxPionTrunc40->GetNbinsY()+1;index_ybin++)
        {
            if(h2PoverMDeDxPionTrunc40->ProjectionY()->GetBinLowEdge(index_ybin)>=profPoverMDeDxPionTrunc40->GetBinContent(index_xbin)+0.7) 
            {
                h2PoverMDeDxPionTrunc40Clone->SetBinContent(index_xbin,index_ybin,h2PoverMDeDxPionTrunc40->GetBinContent(index_xbin,index_ybin));
            }
        }
    }
    DrawHisto(*fileout,h2PoverMDeDxPionTrunc40Clone,"","p/M","dE/dx_{track} [MeV/cm]");

    h2PoverMDeDxPionTrunc40Clone->ProjectionX()->Write();
    h2PoverMDeDxPionTrunc40->ProjectionX()->Write();

    TEfficiency* eff_pionTrunc40 = new TEfficiency(*h2PoverMDeDxPionTrunc40Clone->ProjectionX(),*h2PoverMDeDxPionTrunc40->ProjectionX());
 
    eff_pionTrunc40->Write();





    SuperposedHisto2DProfile(*csuperposed,h2PoverMDeDxPionRMS,profPoverMDeDxPionRMS,"","p/M","dE/dx_{track} [MeV/cm]","[u.a.]");
    csuperposed->Write();


    for(int index_xbin=0;index_xbin<h2PoverMDeDxPionRMS->GetNbinsX()+1;index_xbin++)
    {
        for(int index_ybin=0;index_ybin<h2PoverMDeDxPionRMS->GetNbinsY()+1;index_ybin++)
        {
            if(h2PoverMDeDxPionRMS->ProjectionY()->GetBinLowEdge(index_ybin)>=profPoverMDeDxPionRMS->GetBinContent(index_xbin)+0.7) 
            {
                h2PoverMDeDxPionRMSClone->SetBinContent(index_xbin,index_ybin,h2PoverMDeDxPionRMS->GetBinContent(index_xbin,index_ybin));
            }
        }
    }
    DrawHisto(*fileout,h2PoverMDeDxPionRMSClone,"","p/M","dE/dx_{track} [MeV/cm]");

    h2PoverMDeDxPionRMSClone->ProjectionX()->Write();
    h2PoverMDeDxPionRMS->ProjectionX()->Write();

    TEfficiency* eff_pionRMS = new TEfficiency(*h2PoverMDeDxPionRMSClone->ProjectionX(),*h2PoverMDeDxPionRMS->ProjectionX());
 
    eff_pionRMS->Write();


    



    TLegend* LegendEff = new TLegend(0.7,0.7,0.9,0.9);
    LegendEff->AddEntry(eff_pionCluster,"Cluster","l");
    LegendEff->AddEntry(eff_pionRMS,"RMS","l");
    LegendEff->AddEntry(eff_pionTrunc40,"Trunc40","l");
    LegendEff->AddEntry(eff_pionHarmonic2,"Harmonic2","l");
    
    TCanvas* csuperposedEfficiency = new TCanvas();
    
    eff_pionCluster->Draw();
    eff_pionRMS->Draw("same");
    eff_pionRMS->SetLineColor(2);
    eff_pionTrunc40->Draw("same");
    eff_pionTrunc40->SetLineColor(4);
    eff_pionHarmonic2->Draw("same");
    eff_pionHarmonic2->SetLineColor(6);
    LegendEff->Draw("same");
    csuperposedEfficiency->Write();
    





    MassVsIas->Write();
    MassBeforeFactorWithCorrection->Write();
    MassBeforeFactorWithoutCorrection->Write();
    MassPt5->Write();
    MassPt20->Write();
    MassPt35->Write();
    MassPt50->Write();
    MassPtcorr5->Write();
    MassPtcorr20->Write();
    MassPtcorr35->Write();
    MassPtcorr50->Write();







    TCanvas* cSUSY3 = new TCanvas();
    TLegend* legSUSY3 = new TLegend(0.7,0.7,0.9,0.9);
    TLine* lineSUSY3 = new TLine(2400,0,2400,400);
    legSUSY3->AddEntry(hDistribMassSUSYeloss,"M(E_{sim})","l");
    legSUSY3->AddEntry(hDistribMassSUSYqcorr,"M(E_{corr})","l");
    legSUSY3->AddEntry(hDistribMassSUSYqcorrinf,("#Xi < "+to_string((int)(ratiosat*100))+" %").c_str(),"l");
    legSUSY3->AddEntry(hDistribMassSUSYqcorrsup,("#Xi > "+to_string((int)(ratiosat*100))+" %").c_str(),"l");
    //hDistribMassSUSYqcorrinf->Scale(1./(hDistribMassSUSYqcorr->Integral()-hDistribMassSUSYqcorrsup->Integral()));
    //hDistribMassSUSYqcorrsup->Scale(1./(hDistribMassSUSYqcorr->Integral()-hDistribMassSUSYqcorrinf->Integral()));
    hDistribMassSUSYeloss->Draw();
    hDistribMassSUSYqcorr->Draw("same");
    hDistribMassSUSYqcorrinf->Draw("same");
    hDistribMassSUSYqcorrsup->Draw("same");
    legSUSY3->Draw("same");
    lineSUSY3->Draw("same");
    hDistribMassSUSYeloss->GetXaxis()->SetTitle("m [GeV/c^{2}]");
    hDistribMassSUSYeloss->GetYaxis()->SetTitle("[u.a.]");    
    cSUSY3->Write();

    DrawHistoNormalized(*fileout,h2ElossvQ_Sat,"","E_r","E_s","u.a.");

    TCanvas *celosspath = new TCanvas();
    DrawHistoNormalized(*celosspath,Vectelosspath,VectelosspathLegend,"","");
    celosspath->Write();

    elosspathTIB1->Write();
    elosspathTOB1->Write();
    elossTIB1->Write();
    elossTOB1->Write();

    vector<TH1D*> VectMass;
    VectMass.push_back(hDistribMassPioneloss);
    VectMass.push_back(hDistribMassKaoneloss);
    VectMass.push_back(hDistribMassProtoneloss);
    
    //VectMass.push_back(hDistribMassMSeloss);
    vector<string> VectLegendMass;
    VectLegendMass.push_back("Pion");
    VectLegendMass.push_back("Kaon");
    VectLegendMass.push_back("Proton");
    
    //VectLegendMass.push_back("MS");
    TCanvas* cMass = new TCanvas();
    StackHisto(*cMass,VectMass,VectLegendMass,"","m [GeV/c^{2}]");
    cMass->Write();


    DrawHisto(*fileout,h1TestDeDxProtonFixMom,"Estimateur harmonique-2","#frac{dE_{sim}}{dx}");
    DrawHisto(*fileout,h1TestDeDxProtonFixMomNoEstim,"Sans estimateur","#frac{dE_{sim}}{dx}");

    /*h2PoverMDeDxPion->SetMarkerColor(36);
    h2PoverMDeDxPion->SetMarkerStyle(25);
    h2PoverMDeDxKaon->SetMarkerColor(38);
    h2PoverMDeDxKaon->SetMarkerStyle(25);
    h2PoverMDeDxProton->SetMarkerColor(40);
    h2PoverMDeDxProton->SetMarkerStyle(25);
    h2PoverMDeDxRhad->SetMarkerColor(42);
    h2PoverMDeDxRhad->SetMarkerStyle(25);
    h2PoverMDeDxSUSY->SetMarkerColor(44);
    h2PoverMDeDxSUSY->SetMarkerStyle(25);*/

    vector<TH2F*> VectPoverMpart;
    //VectPoverMpart.push_back(h2PoverMDeDxPion);
    //VectPoverMpart.push_back(h2PoverMDeDxKaon);
    VectPoverMpart.push_back(h2PoverMDeDxProton);
    //VectPoverMpart.push_back(h2PoverMDeDxRhad);
    VectPoverMpart.push_back(h2PoverMDeDxSUSY);
    vector<string> VectLegendPoverMpart;
    //VectLegendPoverMpart.push_back("Pion");
    //VectLegendPoverMpart.push_back("Kaon");
    VectLegendPoverMpart.push_back("Proton");
    //VectLegendPoverMpart.push_back("Rhad");
    VectLegendPoverMpart.push_back("SUSY");
    DrawHisto(*fileout,VectPoverMpart,VectLegendPoverMpart,"","p/M","dE/dx [MeV/cm]");

    TCanvas* cpoverm = new TCanvas();
    h2PoverMDeDxPionHarmonic2->Draw("same");
    h2PoverMDeDxKaon->Draw("same");
    h2PoverMDeDxProton->Draw("same");
    //h2PoverMDeDxRhad->Draw("same");
    h2PoverMDeDxSUSY->Draw("same");
    cpoverm->Write();


    DrawHisto(*fileout,hDistribMassSMeloss,"","m [GeV/c^{2}]");
    DrawHisto(*fileout,hDistribMassSUSYeloss,"","m [GeV/c^{2}]");
    DrawHisto(*fileout,hDistribMassSMq,"","m [GeV/c^{2}]");
    DrawHisto(*fileout,hDistribMassSUSYq,"","m [GeV/c^{2}]");
    DrawHisto(*fileout,hDistribMassSMqcorr,"","m [GeV/c^{2}]");
    DrawHisto(*fileout,hDistribMassSUSYqcorr,"","m [GeV/c^{2}]");
    DrawHisto(*fileout,hDistribMassttbarqcorr,"","m [GeV/c^{2}]");

    
    
    vector<TH1D*> VectDistribMassSM;
    VectDistribMassSM.push_back(hDistribMassSMeloss);
    VectDistribMassSM.push_back(hDistribMassSMq);
    VectDistribMassSM.push_back(hDistribMassSMqcorr);
    
    vector<string> VectLegendDistribMass;
    VectLegendDistribMass.push_back("M(E_{sim})");
    VectLegendDistribMass.push_back("M(E_{rec})");
    vector<TH1D*> VectDistribMassSUSY;
    VectDistribMassSUSY.push_back(hDistribMassSUSYeloss);
    VectDistribMassSUSY.push_back(hDistribMassSUSYq);

    DrawHisto(*fileout,VectDistribMassSM,VectLegendDistribMass,"","m [GeV/c^{2}]");
    TCanvas* cSusy = new TCanvas();
    DrawHisto(*cSusy,VectDistribMassSUSY,VectLegendDistribMass,"","m [GeV/c^{2}]");
    TLine* line2400 = new TLine(2400,0,2400,600);
    line2400->Draw("same");
    cSusy->Write();

    vector<string> VectLegendDistribMass2;
    VectLegendDistribMass2.push_back("M(E_{sim})");
    //VectLegendDistribMass2.push_back(("M(E_{corr}, Seuil > "+to_string((int)(ratiosat*100))+" %)").c_str());
    //VectLegendDistribMass2.push_back(("M(E_{corr}, Seuil < "+to_string((int)(ratiosat*100))+" %)").c_str());
    VectLegendDistribMass2.push_back("M(E_{full})");
    vector<TH1D*> VectDistribMassSUSY2;
    VectDistribMassSUSY2.push_back(hDistribMassSUSYeloss);
    //VectDistribMassSUSY2.push_back(hDistribMassSUSYqcorrsup);
    //VectDistribMassSUSY2.push_back(hDistribMassSUSYqcorrinf);
    VectDistribMassSUSY2.push_back(hDistribMassSUSYqcorr);
    TPaveText* pave2 = new TPaveText(0.4,0.4,0.6,0.6,"NDC");
    pave2->AddText(("Seuil : "+to_string((int)(ratiosat*100))+" %").c_str());

    DrawHisto(*fileout,VectDistribMassSM,VectLegendDistribMass2,"","m [GeV/c^{2}]");
    TCanvas* cSusy2 = new TCanvas();
    DrawHisto(*cSusy2,VectDistribMassSUSY2,VectLegendDistribMass2,"","m [GeV/c^{2}]");
    line2400->Draw("same");
    pave2->Draw("same");
    cSusy2->Write();



    hInvestigationCharge->Write();
    hInvestigationEloss->Write();
    h2InvestigationEvQ->Write();
    hInvestigationEta->Write();
    hInvestigationPhi->Write();
    hInvestigationDetId->Write();
    hInvestigationEntries->Write();
    hInvestigationShape->Write();
    hInvestigationSimHits->Write();

    TLine* linexy = new TLine(0,0,0.006,0.006);
    linexy->SetLineStyle(4);
    TCanvas* cEvQbarrel = new TCanvas("EvQbarrel","EvQbarrel");
    EvQbarrel->Draw("colz");
    EvQbarrel->GetXaxis()->SetTitle("E_{sim} [MeV]");
    EvQbarrel->GetYaxis()->SetTitle("E_{rec} [MeV]");
    linexy->Draw("SAME");
    cEvQbarrel->Write();
    TCanvas* cEvQbarrel_corr = new TCanvas("EvQbarrel_corr","EvQbarrel_corr");
    EvQbarrel_corr->Draw("colz");
    EvQbarrel_corr->GetXaxis()->SetTitle("E_{sim} [MeV]");
    EvQbarrel_corr->GetYaxis()->SetTitle("E_{rec} [MeV]");
    linexy->Draw("SAME");
    cEvQbarrel_corr->Write();

    DrawHisto(*fileout,h2Shape,"Shape=false","E_{sim}","E_{rec}");

    TPaveText *pave1 = new TPaveText(0.6,0.5,0.8,0.7,"NDC");
    pave1->SetFillColor(0);
    pave1->AddText(("Ratio de saturation >= "+to_string((int)(ratiosat*100))+"%").c_str());
    TCanvas* cDiffRelChargecorr_Charge = new TCanvas("cDiffRelChargecorr_Charge","cDiffRelChargecorr_Charge");
    vector<TH1D*> VectDiffRelChargecorr_Charge;
    VectDiffRelChargecorr_Charge.push_back(h1DiffRelEvQ);
    VectDiffRelChargecorr_Charge.push_back(h1DiffRelEvQcorr);
    vector<string> VectLegendDiffRelChargecorr_Charge;
    VectLegendDiffRelChargecorr_Charge.push_back("Avant correction");
    VectLegendDiffRelChargecorr_Charge.push_back("Apres correction");
    DrawHistoNormalized(*cDiffRelChargecorr_Charge,VectDiffRelChargecorr_Charge,VectLegendDiffRelChargecorr_Charge,"Correction","#frac{E_{corr}-E_{sim}}{E_{sim}}");
    pave1->Draw("same");
    VectDiffRelChargecorr_Charge[0]->GetXaxis()->SetRangeUser(-1,2);
    cDiffRelChargecorr_Charge->SetGrid();
    cDiffRelChargecorr_Charge->Write();


    

    TCanvas* cTestCutEdge = new TCanvas("cTestCutEdge","cTestCutEdge");
    vector<TH1D*> VectTestCutEdge;
    VectTestCutEdge.push_back(h1TestNoCutAndNoEdge);
    VectTestCutEdge.push_back(h1TestCutOrEdge);
    vector<string> VectLegendTestCutEdge;
    VectLegendTestCutEdge.push_back("No cut & no edge");
    VectLegendTestCutEdge.push_back("Cut or edge");
    DrawHistoNormalized(*cTestCutEdge,VectTestCutEdge,VectLegendTestCutEdge,"Test cut and edge","#frac{E_{sim}}{E_{rec}}");
    cTestCutEdge->Write();

    
    TCanvas* cNStrip5Sat1_254_1=new TCanvas("cNStrip5Sat1_254_1","cNStrip5Sat1_254_1");
    SuperposedHisto2DProfile(*cNStrip5Sat1_254_1,h2NStrip5Sat1_254_1,profNstrip5Sat1_254_1,"NStrip=5 & NStripSat=1 | TOB1 | Selected Area | 0.0<P/M<0.2","E_{sim}","E_{rec}");
    cNStrip5Sat1_254_1->Write();
    TCanvas* cNStrip5Sat1_254_2=new TCanvas("cNStrip5Sat1_254_2","cNStrip5Sat1_254_2");
    SuperposedHisto2DProfile(*cNStrip5Sat1_254_2,h2NStrip5Sat1_254_2,profNstrip5Sat1_254_2,"NStrip=5 & NStripSat=1 | TOB1 | Selected Area | 0.2<P/M<0.4","E_{sim}","E_{rec}");
    cNStrip5Sat1_254_2->Write();
    TCanvas* cNStrip5Sat1_254_3=new TCanvas("cNStrip5Sat1_254_3","cNStrip5Sat1_254_3");
    SuperposedHisto2DProfile(*cNStrip5Sat1_254_3,h2NStrip5Sat1_254_3,profNstrip5Sat1_254_3,"NStrip=5 & NStripSat=1 | TOB1 | Selected Area | 0.4<P/M<0.6","E_{sim}","E_{rec}");
    cNStrip5Sat1_254_3->Write();
    TCanvas* cNStrip5Sat1_254_4=new TCanvas("cNStrip5Sat1_254_4","cNStrip5Sat1_254_4");
    SuperposedHisto2DProfile(*cNStrip5Sat1_254_4,h2NStrip5Sat1_254_4,profNstrip5Sat1_254_4,"NStrip=5 & NStripSat=1 | TOB1 | Selected Area | 0.6<P/M<0.8","E_{sim}","E_{rec}");
    cNStrip5Sat1_254_4->Write();
    TCanvas* cNStrip5Sat1_254_5=new TCanvas("cNStrip5Sat1_254_5","cNStrip5Sat1_254_5");
    SuperposedHisto2DProfile(*cNStrip5Sat1_254_5,h2NStrip5Sat1_254_5,profNstrip5Sat1_254_5,"NStrip=5 & NStripSat=1 | TOB1 | Selected Area | 0.8<P/M<1.0","E_{sim}","E_{rec}");
    cNStrip5Sat1_254_5->Write();

    DrawHisto(*fileout,h1ChargeFauxSimHit,"Distribution de charge des faux SimHits","E_{rec}");
    DrawHisto(*fileout,h2TestCutEdge,"Cut or Edge","E_{sim}","E_{rec}");
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
    TCanvas* cPtPartIDNormalized = new TCanvas("cPtPartIDNormalized","cPtPartIDNormalized");
    DrawHistoNormalized(*cPtPartIDNormalized,VectPartID_Pt,VectLegendPartId,"Distribution en p_{T} pour differentes particules","p_{T}");
    cPtPartIDNormalized->Write();

    TCanvas* cElossPartID = new TCanvas("cElossPartID","cElossPartID");
    StackHisto(*cElossPartID,VectPartID_Eloss,VectLegendPartId,"Distribution E_{sim} pour differentes particules dans TOB1","E_{sim}");
    cElossPartID->Write();
    TCanvas* cElossPartIDNormalized = new TCanvas("cElossPartIDNormalized","cElossPartIDNormalized");
    DrawHistoNormalized(*cElossPartIDNormalized,VectPartID_Eloss,VectLegendPartId,"Distribution E_{sim} pour differentes particules dans TOB1","E_{sim}");
    cElossPartIDNormalized->Write();

    TCanvas* cPoverMPartID = new TCanvas("cPoverMPartID","cPoverMPartID");
    StackHisto(*cPoverMPartID,VectPartID_PoverM,VectLegendPartId,"Distribution du P/M pour differentes particules","P/M");
    cPoverMPartID->Write();
    TCanvas* cPoverMPartIDNormalized = new TCanvas("cPoverMPartIDNormalized","cPoverMPartIDNormalized");
    DrawHistoNormalized(*cPoverMPartIDNormalized,VectPartID_PoverM,VectLegendPartId,"Distribution du P/M pour differentes particules","P/M");
    cPoverMPartIDNormalized->Write();

    DrawHisto(*fileout,hNStripPerCluster,"Distribution du nombre de strips par cluster","# Strip");
    DrawHisto(*fileout,hStripMax,"Distribution strip max","Q_{max}");
    DrawHisto(*fileout,h2StripMaxEoverQ,"E/Q=f(Q_{max})","Q_{max}","E/Q");
    vector<TH1D*> VectDeDx;
    VectDeDx.push_back(hDeDxHarmo);
    VectDeDx.push_back(hDeDxTrunc);
    vector<string> VectLegendDeDx;
    VectLegendDeDx.push_back("Harmonique");
    VectLegendDeDx.push_back("Troncature");
    TCanvas* cDeDx = new TCanvas("cDeDx","cDeDx");
    DrawHistoNormalized(*cDeDx,VectDeDx,VectLegendDeDx,"Distribution de #frac{dE_{sim}}{dx} pour differents estimateurs","#frac{dE_{sim}}{dx}");
    cDeDx->Write();

    DrawHisto(*fileout,h2PoverMDeDx,"Bethe & Bloch","P/M","#frac{dE_{sim}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PoverMDqDx,"Bethe & Bloch","P/M","#frac{dE_{rec}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PoverMDqcorrDx,"Bethe & Bloch","P/M","#frac{dE_{corr}}{dx} [MeV/cm]");

    DrawHisto(*fileout,h2PoverMDeDxNoEstim,"Bethe & Bloch","P/M","#frac{dE_{sim}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PoverMDqDxNoEstim,"Bethe & Bloch","P/M","#frac{dE_{rec}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PoverMDqcorrDxNoEstim,"Bethe & Bloch","P/M","#frac{dE_{corr}}{dx} [MeV/cm]");

    DrawHisto(*fileout,h2PDeDx,"","impulsion [GeV/c]","#frac{dE_{sim}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PDqDx,"","impulsion [GeV/c]","#frac{dE_{rec}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PDqcorrDx,"","impulsion [GeV/c]","#frac{dE_{corr}}{dx} [MeV/cm]");

    DrawHisto(*fileout,h2PDeDxNoEstim,"","impulsion [GeV/c]","#frac{dE_{sim}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PDqDxNoEstim,"","impulsion [GeV/c]","#frac{dE_{rec}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PDqcorrDxNoEstim,"","impulsion [GeV/c]","#frac{dE_{corr}}{dx} [MeV/cm]");

    DrawHisto(*fileout,h2PoverMDeDxLarge,"Bethe & Bloch","P/M","#frac{dE_{sim}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PoverMDqDxLarge,"Bethe & Bloch","P/M","#frac{dE_{rec}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PoverMDqcorrDxLarge,"Bethe & Bloch","P/M","#frac{dE_{corr}}{dx} [MeV/cm]");

    DrawHisto(*fileout,h2PoverMDeDxLargeNoEstim,"Bethe & Bloch","P/M","#frac{dE_{sim}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PoverMDqDxLargeNoEstim,"Bethe & Bloch","P/M","#frac{dE_{rec}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PoverMDqcorrDxLargeNoEstim,"Bethe & Bloch","P/M","#frac{dE_{corr}}{dx} [MeV/cm]");

    DrawHisto(*fileout,h2PDeDxLarge,"","impulsion [GeV/c]","#frac{dE_{sim}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PDqDxLarge,"","impulsion [GeV/c]","#frac{dE_{rec}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PDqcorrDxLarge,"","impulsion [GeV/c]","#frac{dE_{corr}}{dx} [MeV/cm]");

    DrawHisto(*fileout,h2PDeDxLargeNoEstim,"","impulsion [GeV/c]","#frac{dE_{sim}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PDqDxLargeNoEstim,"","impulsion [GeV/c]","#frac{dE_{rec}}{dx} [MeV/cm]");
    DrawHisto(*fileout,h2PDqcorrDxLargeNoEstim,"","impulsion [GeV/c]","#frac{dE_{corr}}{dx} [MeV/cm]");
    

    TCanvas* cPoverMPartIDNormal = new TCanvas("cPoverMPartIDNormal","cPoverMPartIDNormal");
    vector<TH1D*> VectPartID_PoverM_Normalized;
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
        DrawHisto(*fileout,VectPartID_h2EvQ_NoSat[i],("NoSat_"+LoopPartID(i)).c_str(),"E_{sim}","E_{rec}");
        DrawHisto(*fileout,VectPartID_h2EvQ_Sat[i],("Sat_"+LoopPartID(i)).c_str(),"E_{sim}","E_{rec}");
        DrawHisto(*fileout,VectPartID_h2EvQ_Sat254[i],("Sat254_"+LoopPartID(i)).c_str(),"E_{sim}","E_{rec}");
        DrawHisto(*fileout,VectPartID_h2EvQ_Sat255[i],("Sat255_"+LoopPartID(i)).c_str(),"E_{sim}","E_{rec}");
        DrawHisto(*fileout,VectPartID_Pt[i],(LoopPartID(i)+"_p_{T}").c_str(),"p_{T}");
        DrawHisto(*fileout,VectPartID_Eloss[i],(LoopPartID(i)+"_E_{sim}").c_str(),"E_{sim}");
        DrawHisto(*fileout,VectPartID_PoverM[i],(LoopPartID(i)+"_P/M").c_str(),"P/M");
    }
    
    TCanvas* cRatioSatPt254 = new TCanvas("cRatioSatPt254","cRatioSatPt254");
    SuperposedHisto2DProfile(*cRatioSatPt254,h2RatioSatPt254,profSatPt254,"Ratio de clusters qui saturent a 254 en fonction du p_{T}","p_{T}","#frac{Nombre de clusters qui saturent}{Nombre de clusters de la trace}");
    cRatioSatPt254->Write();

    TCanvas* cRatioSatEloss254 = new TCanvas("cRatioSatEloss254","cRatioSatEloss254");
    SuperposedHisto2DProfile(*cRatioSatEloss254,h2RatioSatEloss254,profSatEloss254,"Ratio de clusters qui saturent a 254 en fonction de E_{sim}","E_{sim}","#frac{Nombre de clusters qui saturent}{Nombre de clusters de la trace}");
    cRatioSatEloss254->Write();

    TCanvas* cRatioSatQ254 = new TCanvas("cRatioSatQ254","cRatioSatQ254");
    SuperposedHisto2DProfile(*cRatioSatQ254,h2RatioSatQ254,profSatQ254,"Ratio de clusters qui saturent a 254 en fonction de E_{rec}","E_{rec}","#frac{Nombre de clusters qui saturent}{Nombre de clusters de la trace}");
    cRatioSatQ254->Write();

    vector<TH1D*> VectSatPartID;

    for(int bin1=1;bin1<h2RatioSatPartID254->GetXaxis()->GetNbins()+1;bin1++)
    {
        int nbre_entries=0;
        TH1D* h1D = new TH1D("","",15,0,1.2);
        for(int bin2=1;bin2<h2RatioSatPartID254->GetYaxis()->GetNbins()+1;bin2++) 
        {
            nbre_entries+=h2RatioSatPartID254->GetBinContent(bin1,bin2);
        }
        for(int bin2=1;bin2<h2RatioSatPartID254->GetYaxis()->GetNbins()+1;bin2++)
        {  
            h2RatioSatPartID254->SetBinContent(bin1,bin2,h2RatioSatPartID254->GetBinContent(bin1,bin2)/nbre_entries);
        }

        for(int bin2=1;bin2<h2RatioSatPartID254->GetYaxis()->GetNbins()+1;bin2++)
        {
            h1D->SetBinContent(bin2,(h2RatioSatPartID254->GetBinContent(bin1,bin2)));
        }
        h1D->Write();
        VectSatPartID.push_back(h1D);
    }

    DrawHisto(*fileout,VectSatPartID,VectLegendPartId,"","#frac{Nombre de clusters qui saturent}{Nombre de clusters de la trace}");


    /*ratsatpion->Scale(1./ratsatpion->GetEntries());
    ratsatkaon->Scale(1./ratsatkaon->GetEntries());
    ratsatproton->Scale(1./ratsatproton->GetEntries());
    ratsatMS->Scale(1./ratsatMS->GetEntries());
    ratsatrhad->Scale(1./ratsatrhad->GetEntries());
    ratsatSUSY->Scale(1./ratsatSUSY->GetEntries());*/


    for(int bin1=1;bin1<ratsatpion->GetNbinsX()+1;bin1++) ratsatpion->SetBinContent(bin1,ratsatpion->GetBinContent(bin1)/ratsatpion->GetEntries());
    for(int bin1=1;bin1<ratsatkaon->GetNbinsX()+1;bin1++) ratsatkaon->SetBinContent(bin1,ratsatkaon->GetBinContent(bin1)/ratsatkaon->GetEntries());
    for(int bin1=1;bin1<ratsatproton->GetNbinsX()+1;bin1++) ratsatproton->SetBinContent(bin1,ratsatproton->GetBinContent(bin1)/ratsatproton->GetEntries());
    for(int bin1=1;bin1<ratsatMS->GetNbinsX()+1;bin1++) ratsatMS->SetBinContent(bin1,ratsatMS->GetBinContent(bin1)/ratsatMS->GetEntries());
    for(int bin1=1;bin1<ratsatrhad->GetNbinsX()+1;bin1++) ratsatrhad->SetBinContent(bin1,ratsatrhad->GetBinContent(bin1)/ratsatrhad->GetEntries());
    for(int bin1=1;bin1<ratsatSUSY->GetNbinsX()+1;bin1++) ratsatSUSY->SetBinContent(bin1,ratsatSUSY->GetBinContent(bin1)/ratsatSUSY->GetEntries());

    for(int bin1=1;bin1<ratsatpion255->GetNbinsX()+1;bin1++) ratsatpion255->SetBinContent(bin1,ratsatpion255->GetBinContent(bin1)/ratsatpion255->GetEntries());
    for(int bin1=1;bin1<ratsatkaon255->GetNbinsX()+1;bin1++) ratsatkaon255->SetBinContent(bin1,ratsatkaon255->GetBinContent(bin1)/ratsatkaon255->GetEntries());
    for(int bin1=1;bin1<ratsatproton255->GetNbinsX()+1;bin1++) ratsatproton255->SetBinContent(bin1,ratsatproton255->GetBinContent(bin1)/ratsatproton255->GetEntries());
    for(int bin1=1;bin1<ratsatMS255->GetNbinsX()+1;bin1++) ratsatMS255->SetBinContent(bin1,ratsatMS255->GetBinContent(bin1)/ratsatMS255->GetEntries());
    for(int bin1=1;bin1<ratsatrhad255->GetNbinsX()+1;bin1++) ratsatrhad255->SetBinContent(bin1,ratsatrhad255->GetBinContent(bin1)/ratsatrhad255->GetEntries());
    for(int bin1=1;bin1<ratsatSUSY255->GetNbinsX()+1;bin1++) ratsatSUSY255->SetBinContent(bin1,ratsatSUSY255->GetBinContent(bin1)/ratsatSUSY255->GetEntries());


    vector<TH1D*> VectSatPartId2;
    VectSatPartId2.push_back(ratsatpion);
    VectSatPartId2.push_back(ratsatkaon);
    VectSatPartId2.push_back(ratsatproton);
    VectSatPartId2.push_back(ratsatMS);
    VectSatPartId2.push_back(ratsatrhad);
    VectSatPartId2.push_back(ratsatSUSY);
    vector<TH1D*> VectSatPartId2255;
    VectSatPartId2255.push_back(ratsatpion255);
    VectSatPartId2255.push_back(ratsatkaon255);
    VectSatPartId2255.push_back(ratsatproton255);
    VectSatPartId2255.push_back(ratsatMS255);
    VectSatPartId2255.push_back(ratsatrhad255);
    VectSatPartId2255.push_back(ratsatSUSY255);
    vector<string> VectSatPartId2Legend;
    VectSatPartId2Legend.push_back("Pion");
    VectSatPartId2Legend.push_back("Kaon");
    VectSatPartId2Legend.push_back("Proton");
    VectSatPartId2Legend.push_back("MS");
    VectSatPartId2Legend.push_back("R^{+}_{#tilde{g}u#bar{d}}");
    VectSatPartId2Legend.push_back("SUSY");
    DrawHisto(*fileout,VectSatPartId2,VectSatPartId2Legend,"Saturation 254","Rat.Sat.");
    DrawHisto(*fileout,VectSatPartId2255,VectSatPartId2Legend,"Saturation 255","Rat.Sat.");


    TCanvas* cRatioSatPartID254 = new TCanvas("cRatioSatPartID254","cRatioSatPartID254");
    SetHistoLabelPartID(cRatioSatPartID254,h2RatioSatPartID254);
    profSatPartID254->SetMarkerStyle(2);
    profSatPartID254->SetMarkerColor(kRed);
    profSatPartID254->SetLineColor(kRed);
    profSatPartID254->Draw("SAME");
    h2RatioSatPartID254->GetXaxis()->SetTitle("");
    h2RatioSatPartID254->SetDrawOption("COLZ");
    h2RatioSatPartID254->GetYaxis()->SetTitle("#frac{Nombre de clusters qui saturent}{Nombre de clusters de la trace}");
    h2RatioSatPartID254->GetZaxis()->SetTitle("[u.a.]");
    //h2RatioSatPartID254->SetMaximum(800);
    cRatioSatPartID254->Write();

    TCanvas* cRatioSatPt255 = new TCanvas("cRatioSatPt255","cRatioSatPt255");
    SuperposedHisto2DProfile(*cRatioSatPt255,h2RatioSatPt255,profSatPt255,"Ratio de clusters qui saturent a 255 en fonction du p_{T}","p_{T}","#frac{#ClusterSat}{#Cluster}");
    cRatioSatPt255->Write();

    TCanvas* cRatioSatEloss255 = new TCanvas("cRatioSatEloss255","cRatioSatEloss255");
    SuperposedHisto2DProfile(*cRatioSatEloss255,h2RatioSatEloss255,profSatEloss255,"Ratio de clusters qui saturent a 255 en fonction de E_{sim}","E_{sim}","#frac{#ClusterSat}{#Cluster}");
    cRatioSatEloss255->Write();

    TCanvas* cRatioSatQ255 = new TCanvas("cRatioSatQ255","cRatioSatQ255");
    SuperposedHisto2DProfile(*cRatioSatQ255,h2RatioSatQ255,profSatQ255,"Ratio de clusters qui saturent a 255 en fonction de E_{rec}","E_{rec}","#frac{#ClusterSat}{#Cluster}");
    cRatioSatQ255->Write();

    for(int bin1=0;bin1<h2RatioSatPartID255->GetXaxis()->GetNbins()+1;bin1++)
    {
        for(int bin2=0;bin2<h2RatioSatPartID255->GetYaxis()->GetNbins()+1;bin2++)
        {
            h2RatioSatPartID255->SetBinContent(bin1,bin2,h2RatioSatPartID255->GetBinContent(bin1,bin2)/h2RatioSatPartID255->GetEntries());
        }
    }

    TCanvas* cRatioSatPartID255 = new TCanvas("cRatioSatPartID255","cRatioSatPartID255");
    SetHistoLabelPartID(cRatioSatPartID255,h2RatioSatPartID255);
    profSatPartID255->SetMarkerStyle(2);
    profSatPartID255->SetMarkerColor(kRed);
    profSatPartID255->SetLineColor(kRed);
    profSatPartID255->Draw("SAME");
    h2RatioSatPartID255->GetXaxis()->SetTitle("");
    h2RatioSatPartID255->SetDrawOption("COLZ");
    h2RatioSatPartID255->GetZaxis()->SetTitle("[u.a.]");
    //h2RatioSatPartID255->SetMaximum(800);
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

    DrawHisto(*fileout,h2ElossvQ_tot,"E_{sim} v. Q","E_{sim} [MeV]","E_{rec} [MeV]");
    DrawHisto(*fileout,h2ElossvQ_NoSat,"E_{sim} v. Q","E_{sim} [MeV]","E_{rec} [MeV]");
    DrawHisto(*fileout,h2ElossvQ_Sat,"E_{sim} v. Q","E_{sim} [MeV]","E_{rec} [MeV]");
    DrawHisto(*fileout,h2ElossvQ_Sat254,"E_{sim} v. Q","E_{sim} [MeV]","E_{rec} [MeV]");
    DrawHisto(*fileout,h2ElossvQ_Sat255,"E_{sim} v. Q","E_{sim} [MeV]","E_{rec} [MeV]");

    DrawProfile(*fileout,profElossVsQ_tot,"E_{sim} v. Q","E_{sim} [MeV]","E_{rec} [MeV]");
    DrawProfile(*fileout,profElossVsQ_NoSat,"E_{sim} v. Q","E_{sim} [MeV]","E_{rec} [MeV]");
    DrawProfile(*fileout,profElossVsQ_Sat,"E_{sim} v. Q","E_{sim} [MeV]","E_{rec} [MeV]");
    DrawProfile(*fileout,profElossVsQ_Sat254,"E_{sim} v. Q","E_{sim} [MeV]","E_{rec} [MeV]");
    DrawProfile(*fileout,profElossVsQ_Sat255,"E_{sim} v. Q","E_{sim} [MeV]","E_{rec} [MeV]");

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
    hEmpty->GetYaxis()->SetTitle("#frac{Nombre total de clusters qui saturent}{Nombre total de clusters}");
    Ratio_LayerSat254->Draw("SAME");
    Ratio_LayerSat255->SetLineColor(2);
    Ratio_LayerSat255->Draw("SAME");
    line254->SetLineStyle(3);
    line254->Draw("SAME");
    line255->SetLineStyle(3);
    line255->SetLineColor(2);
    line255->Draw("SAME");
    //cratioSatLayer->SetLogy();
    TLegend* legratio = new TLegend(0.6,0.7,0.9,0.9);
    legratio->AddEntry(line254,("Saturation 254 : "+to_string(rat254)).c_str(),"l");
    legratio->AddEntry(line255,("Saturation 255 : "+to_string(rat255)).c_str(),"l");
    legratio->Draw("SAME");
    cratioSatLayer->Write();

// ------------------------------------------

// ------------------------------------------ TH2F + TProfile 

    TCanvas* canvas_EvQ_h2prof_tot = new TCanvas("canvas_EvQ_h2prof_tot","canvas_EvQ_h2prof_tot");
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_tot,h2ElossvQ_tot,profElossVsQ_tot,"E_{sim} v. Q","E_{sim} [GeV]","E_{rec} [GeV]","[u.a.]");
    canvas_EvQ_h2prof_tot->Write();

    TCanvas* canvas_EvQ_h2prof_NoSat = new TCanvas("canvas_EvQ_h2prof_NoSat","canvas_EvQ_h2prof_NoSat");
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_NoSat,h2ElossvQ_NoSat,profElossVsQ_NoSat,"E_{sim} v. Q","E_{sim}","E_{rec}");
    canvas_EvQ_h2prof_NoSat->Write();

    TCanvas* canvas_EvQ_h2prof_Sat = new TCanvas("canvas_EvQ_h2prof_Sat","canvas_EvQ_h2prof_Sat");
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_Sat,h2ElossvQ_Sat,profElossVsQ_Sat,"E_{sim} v. Q","E_{sim}","E_{rec}");
    canvas_EvQ_h2prof_Sat->Write();

    TCanvas* canvas_EvQ_h2prof_Sat254 = new TCanvas("canvas_EvQ_h2prof_Sat254","canvas_EvQ_h2prof_Sat254");
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_Sat254,h2ElossvQ_Sat254,profElossVsQ_Sat254,"E_{sim} v. Q","E_{sim}","E_{rec}");
    canvas_EvQ_h2prof_Sat254->Write();

    TCanvas* canvas_EvQ_h2prof_Sat255 = new TCanvas("canvas_EvQ_h2prof_Sat255","canvas_EvQ_h2prof_Sat255");
    SuperposedHisto2DProfile(*canvas_EvQ_h2prof_Sat255,h2ElossvQ_Sat255,profElossVsQ_Sat255,"E_{sim} v. Q","E_{sim}","E_{rec}");
    canvas_EvQ_h2prof_Sat255->Write();

    profElossVsQ_tot->SetMarkerStyle(1);
    profElossVsQ_NoSat->SetMarkerStyle(1);
    profElossVsQ_Sat->SetMarkerStyle(1);
    profElossVsQ_Sat254->SetMarkerStyle(1);
    profElossVsQ_Sat255->SetMarkerStyle(1);

// ------------------------------------------

// ------------------------------------------ TH1D Pt de la trace, avec ou sans saturation

    TCanvas* canvas_Pt = new TCanvas("canvas_Pt","canvas_Pt");
    vector<TH1D*> VectHistoPt;
    VectHistoPt.push_back(hPt_NoSat);
    VectHistoPt.push_back(hPt_Sat);
    vector<string> VectLegendPt;
    VectLegendPt.push_back("No saturation");
    VectLegendPt.push_back("Saturation");
    DrawHistoNormalized(*canvas_Pt,VectHistoPt,VectLegendPt,"Distribution de P_{T}","P_{T}");
    canvas_Pt->Write();

// ------------------------------------------

// ------------------------------------------ TH1D ratio NClusters qui saturent par trace

    TCanvas* canvas_RatioNClusterSat = new TCanvas("canvas_RatioNClusterSat","canvas_RatioNClusterSat");
    vector<TH1D*> VectHistoNSat;
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
    vector<TH1D*> VectHistoEoverQ;
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
        string title = "E_{sim} v. Q  |  "+Label(i+1);
        TCanvas* csuperposed = new TCanvas(title.c_str(),title.c_str());
        DrawHisto(*fileout,VectLayer_h2EvQ_Tot[i],title.c_str(),"E_{sim}","Q");
        DrawHisto(*fileout,VectLayer_h2EvQ_NoSat[i],title.c_str(),"E_{sim}","Q");
        DrawHisto(*fileout,VectLayer_h2EvQ_Sat[i],title.c_str(),"E_{sim}","Q");
        DrawHisto(*fileout,VectLayer_h2EvQ_Sat254[i],title.c_str(),"E_{sim}","Q");
        DrawHisto(*fileout,VectLayer_h2EvQ_Sat255[i],title.c_str(),"E_{sim}","Q");
        DrawHisto(*fileout,VectLayer_h1EvQ_Tot[i],title.c_str(),"E/Q");
        DrawHisto(*fileout,VectLayer_h1EvQ_NoSat[i],title.c_str(),"E/Q");
        DrawHisto(*fileout,VectLayer_h1EvQ_Sat[i],title.c_str(),"E/Q");
        DrawHisto(*fileout,VectLayer_h1EvQ_Sat254[i],title.c_str(),"E/Q");
        DrawHisto(*fileout,VectLayer_h1EvQ_Sat255[i],title.c_str(),"E/Q");
        SuperposedHisto2DProfile(*csuperposed,VectLayer_h2EvQ_Tot[i],VectLayer_profEvQ_Tot[i],title.c_str(),"E_{sim}","Q");
        csuperposed->Write();
        SuperposedHisto2DProfile(*csuperposed,VectLayer_h2EvQ_NoSat[i],VectLayer_profEvQ_NoSat[i],title.c_str(),"E_{sim}","Q");
        csuperposed->Write();
        SuperposedHisto2DProfile(*csuperposed,VectLayer_h2EvQ_Sat[i],VectLayer_profEvQ_Sat[i],title.c_str(),"E_{sim}","Q");
        csuperposed->Write();
        SuperposedHisto2DProfile(*csuperposed,VectLayer_h2EvQ_Sat254[i],VectLayer_profEvQ_Sat254[i],title.c_str(),"E_{sim}","Q");
        csuperposed->Write();
        SuperposedHisto2DProfile(*csuperposed,VectLayer_h2EvQ_Sat255[i],VectLayer_profEvQ_Sat255[i],title.c_str(),"E_{sim}","Q");
        csuperposed->Write();
    }*/

    

    
    /*for(int nstrip=1;nstrip<10;nstrip++)
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
            SuperposedHisto2DProfile(*cNStripNStripSat254_h2prof,VectNStrip_VectNStripSat254_h2EvQ[nstrip-1][nstripsat],VectNStrip_VectNStripSat254_profEvQ[nstrip-1][nstripsat],("H2 & prof NStrip="+to_string(nstrip)+" & NStripSat254="+to_string(nstripsat)+" | TOB1").c_str(),"E_{sim}","Q");
            //VectNStrip_VectNStripSat254_h2EvQ[nstrip-1][nstripsat]->Draw("COLZ");
            linexy->Draw("SAME");
            if(nstripsat>=0)linesat254->Draw("SAME");
            if(nstripsat>=2)linesat254times2->Draw("SAME");
            if(nstripsat>=3)linesat254times3->Draw("SAME");
            if(nstripsat>=4)linesat254times4->Draw("SAME");
            cNStripNStripSat254_h1->cd(nstripsat+1);
            VectNStrip_VectNStripSat254_h1EvQ[nstrip-1][nstripsat]->Draw();
            DrawHisto(*fileout,VectNStrip_VectNStripSat254_h2EvQ[nstrip-1][nstripsat],("H2 NStrip="+to_string(nstrip)+" & NStripSat254="+to_string(nstripsat)).c_str(),"E_{sim}","Q");
            //DrawHisto(*fileout,VectNStrip_VectNStripSat254_h1EvQ[nstrip-1][nstripsat],("H1 NStrip="+to_string(nstrip)+" & NStripSat="+to_string(nstripsat)).c_str(),"E_{sim}");

        }
        cNStripNStripSat254_h2prof->Write();
        cNStripNStripSat254_h1->Write();
    }*/

    /*for(int nstrip=1;nstrip<10;nstrip++)
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
            SuperposedHisto2DProfile(*cNStripNStripSat255_h2prof,VectNStrip_VectNStripSat255_h2EvQ[nstrip-1][nstripsat],VectNStrip_VectNStripSat255_profEvQ[nstrip-1][nstripsat],("H2 & prof NStrip="+to_string(nstrip)+" & NStripSat255="+to_string(nstripsat)+" | TOB1").c_str(),"E_{sim}","Q");
            //VectNStrip_VectNStripSat255_h2EvQ[nstrip-1][nstripsat]->Draw("COLZ");
            linexy->Draw("SAME");
            if(nstripsat>=0) linesat255->Draw("SAME");
            if(nstripsat>=2) linesat255times2->Draw("SAME");
            if(nstripsat>=3) linesat255times3->Draw("SAME");
            if(nstripsat>=4) linesat255times4->Draw("SAME");
            cNStripNStripSat255_h1->cd(nstripsat+1);
            VectNStrip_VectNStripSat255_h1EvQ[nstrip-1][nstripsat]->Draw();
            DrawHisto(*fileout,VectNStrip_VectNStripSat255_h2EvQ[nstrip-1][nstripsat],("H2 NStrip="+to_string(nstrip)+" & NStripSat255="+to_string(nstripsat)).c_str(),"E_{sim}","Q");
            //DrawHisto(*fileout,VectNStrip_VectNStripSat255_h1EvQ[nstrip-1][nstripsat],("H1 NStrip="+to_string(nstrip)+" & NStripSat="+to_string(nstripsat)).c_str(),"E_{sim}");

        }
        cNStripNStripSat255_h2prof->Write();
        cNStripNStripSat255_h1->Write();
    }*/

    fileout->Write();
    fileout->Close();

    delete fileout;

/*    Correction calCharge;
    calCharge.SetFileAndTreeName(s3.c_str(),"tree");
    calCharge.SetBranch();
    for(int countlayer=1;countlayer<11;countlayer++)
    {
        for(int countnstrip=3;countnstrip<7;countnstrip++)
        {
            for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
            //for(int countnstripsat254=0;countnstripsat254<3;countnstripsat254++)
            {
                for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
                //for(int countnstripsat255=0;countnstripsat255<3;countnstripsat255++)
                {
                    calCharge.SetHisto(*VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]);
                    calCharge.FillHisto(countlayer,countnstrip,countnstripsat254,countnstripsat255);
                    calCharge.FillProfile();
                    calCharge.Write(countlayer,countnstrip,countnstripsat254,countnstripsat255);
                }
            }
        }
    }
    calCharge.WriteFile();*/
    
    

    system(("cp "+s2+" ./data/.").c_str());
    

    return 0;
}
