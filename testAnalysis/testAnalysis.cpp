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

using namespace std;


#include "../interface/Builder.h"
#include "../interface/Correction.h"
#include "../interface/Estimator.h"






//CONSTANTES K ET C

const float CsusyCorr = 1.14; //pour gluino
const float KsusyCorr = 6.7;



const float CprotonEloss=3.27;
const float KprotonEloss=2.49;

const float CprotonEcorrNew=6.73;
const float KprotonEcorrNew=1.18;

const float CprotonEcorr=6.15;
const float KprotonEcorr=2.91; 

const float MGluino = 2400.;







//FUNCTIONS -- LABEL HISTO

void StudyGraph(ofstream &ofile, TH1F* histo)
{
    float DistribMean = histo->GetMean();
    float DistribStdDev = histo->GetStdDev();
    TFitResultPtr fit = histo->Fit("gaus","QRS","",DistribMean-DistribStdDev,DistribMean+DistribStdDev);
    float FitMean = fit->Parameter(1);
    float FitSigma = fit->Parameter(2);
    int CountPass=0;
    for(int i=0;i<histo->GetNbinsX();i++)
    {
        //cout<<i<<"   "<<histo->GetBinCenter(i)<<"    "<<histo->GetBinContent(i)<<endl;
        if(histo->GetBinCenter(i)>=FitMean-FitSigma && histo->GetBinCenter(i)<=FitMean+FitSigma) 
        {
            CountPass+=histo->GetBinContent(i);
        }
    }
    float FracInterval = (float)CountPass/(float)histo->GetEntries();
    ofile<<histo->GetName()<<"\t"<<DistribMean<<"\t"<<DistribStdDev<<"\t"<<FitMean<<"\t"<<FitSigma<<"\t"<<FracInterval<<endl;
}

void SetHistoLabel(TH2I* histo)
{
	histo->GetXaxis()->SetBinLabel(1,"TIB 1");
    histo->GetXaxis()->SetBinLabel(2,"TIB 2");
    histo->GetXaxis()->SetBinLabel(3,"TIB 3");
    histo->GetXaxis()->SetBinLabel(4,"TIB 4");
    histo->GetXaxis()->SetBinLabel(5,"TOB 1");
    histo->GetXaxis()->SetBinLabel(6,"TOB 2");
    histo->GetXaxis()->SetBinLabel(7,"TOB 3");
    histo->GetXaxis()->SetBinLabel(8,"TOB 4");
    histo->GetXaxis()->SetBinLabel(9,"TOB 5");
    histo->GetXaxis()->SetBinLabel(10,"TOB 6");
    histo->GetXaxis()->SetBinLabel(11,"TID 1");
    histo->GetXaxis()->SetBinLabel(12,"TID 2");
    histo->GetXaxis()->SetBinLabel(13,"TEC 1");
    histo->GetXaxis()->SetBinLabel(14,"TEC 2");
    histo->GetXaxis()->SetBinLabel(15,"TEC 3");
    histo->GetXaxis()->SetBinLabel(16,"TEC 4");
    histo->GetXaxis()->SetBinLabel(17,"TEC 5");
    histo->GetXaxis()->SetBinLabel(18,"TEC 6");
    histo->GetXaxis()->SetBinLabel(19,"TEC 7");
    histo->GetXaxis()->SetBinLabel(20,"TEC 8");
    histo->GetXaxis()->SetBinLabel(21,"TEC 9");
    histo->GetXaxis()->LabelsOption("v");

	histo->GetYaxis()->SetBinLabel(1,"IB1");
    histo->GetYaxis()->SetBinLabel(2,"IB2");
    histo->GetYaxis()->SetBinLabel(3,"OB1");
    histo->GetYaxis()->SetBinLabel(4,"OB2");
    histo->GetYaxis()->SetBinLabel(5,"W1A");
    histo->GetYaxis()->SetBinLabel(6,"W2A");
    histo->GetYaxis()->SetBinLabel(7,"W3A");
    histo->GetYaxis()->SetBinLabel(8,"W1B");
    histo->GetYaxis()->SetBinLabel(9,"W2B");
    histo->GetYaxis()->SetBinLabel(10,"W3B");
    histo->GetYaxis()->SetBinLabel(11,"W4");
    histo->GetYaxis()->SetBinLabel(12,"W5");
    histo->GetYaxis()->SetBinLabel(13,"W6");
    histo->GetYaxis()->SetBinLabel(14,"W7");
     
}

void SetHistoLabel(TCanvas* canvas,TH1D* histo)
{

    int Wsize = histo->GetMaximum()*1.05;
    histo->GetYaxis()->SetRangeUser(0,histo->GetMaximum()*1.05);

    TLine* line1 = new TLine(4,0,4,0.25);
    TLine* line2 = new TLine(10,0,10,0.25);
    TLine* line3 = new TLine(12,0,12,0.25);

    if(Wsize!=0)
    {
        line1->SetY2(Wsize);
        line2->SetY2(Wsize);
        line3->SetY2(Wsize);
    }


    line1->SetLineStyle(2);
    line2->SetLineStyle(2);
    line3->SetLineStyle(2);
    
    histo->GetXaxis()->SetBinLabel(1,"TIB 1");
    histo->GetXaxis()->SetBinLabel(2,"TIB 2");
    histo->GetXaxis()->SetBinLabel(3,"TIB 3");
    histo->GetXaxis()->SetBinLabel(4,"TIB 4");
    histo->GetXaxis()->SetBinLabel(5,"TOB 1");
    histo->GetXaxis()->SetBinLabel(6,"TOB 2");
    histo->GetXaxis()->SetBinLabel(7,"TOB 3");
    histo->GetXaxis()->SetBinLabel(8,"TOB 4");
    histo->GetXaxis()->SetBinLabel(9,"TOB 5");
    histo->GetXaxis()->SetBinLabel(10,"TOB 6");
    histo->GetXaxis()->SetBinLabel(11,"TID 1");
    histo->GetXaxis()->SetBinLabel(12,"TID 2");
    histo->GetXaxis()->SetBinLabel(13,"TEC 1");
    histo->GetXaxis()->SetBinLabel(14,"TEC 2");
    histo->GetXaxis()->SetBinLabel(15,"TEC 3");
    histo->GetXaxis()->SetBinLabel(16,"TEC 4");
    histo->GetXaxis()->SetBinLabel(17,"TEC 5");
    histo->GetXaxis()->SetBinLabel(18,"TEC 6");
    histo->GetXaxis()->SetBinLabel(19,"TEC 7");
    histo->GetXaxis()->SetBinLabel(20,"TEC 8");
    histo->GetXaxis()->SetBinLabel(21,"TEC 9");
    histo->GetXaxis()->LabelsOption("v");
    canvas->cd();
    histo->Draw();
    line1->Draw();
    line2->Draw();
    line3->Draw();
    
}


string LabelModulGeom(int modulgeom)
{
	if(modulgeom==1) return "IB1";
	if(modulgeom==2) return "IB2";
	if(modulgeom==3) return "OB1";
	if(modulgeom==4) return "OB2";
	if(modulgeom==5) return "W1A";
	if(modulgeom==6) return "W2A";
	if(modulgeom==7) return "W3A";
	if(modulgeom==8) return "W1B";
	if(modulgeom==9) return "W2B";
	if(modulgeom==10) return "W3B";
	if(modulgeom==11) return "W4";
	if(modulgeom==12) return "W5";
	if(modulgeom==13) return "W6";
	if(modulgeom==14) return "W7";
}

string LabelLayer(int i)
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

vector<float> FactorKC(TH2F* histo,float mass) //histo pre-filtre sur le particleID --> il faut donner la masse correspondante 
{

    TProfile profile;       
    TH2F* histoclone = (TH2F*) histo->Clone();
    histoclone->Reset();
    TFitResultPtr fit_res;
    TH1D* projY;
    TH1F* historec = new TH1F("","",1000,0,5);
    int divider=50;
    for(int i=1;i<=histo->GetNbinsX();i++)
    {
        for(int j=2;j<=histo->GetNbinsY();j++)
        {
            if(histo->GetBinContent(i,j)>0) histoclone->SetBinContent(i,j,histo->GetBinContent(i,j));
        }
        if(i%divider==0)
        {
            projY = histoclone->ProjectionY();
            projY->Rebin(5);
            if(projY->GetEntries()>0)
            {
                fit_res = projY->Fit("gaus","QS");
                TH1D* projX = histoclone->ProjectionX();
                //projX->Rebin(5);
                //TFitResultPtr fit_resX = projX->Fit("gaus","QS");
                //historec.SetBinContent(projX->GetMean(),fit_res->Parameter(1));
                //historec.SetBinError(projX->GetMean(),fit_res->Error(1));
                historec->SetBinContent(i-(float)divider/2,fit_res->Parameter(1));
                historec->SetBinError(i-(float)divider/2,fit_res->Error(1));
                historec->SetMarkerStyle(3);
                historec->SetMarkerColor(4);

            }
            histoclone->Reset();
        }
    }
    TF1* InvSquare = new TF1("InvSquare","[0]+[1]*(1/pow(x,2))");
    TFitResultPtr fit = historec->Fit("InvSquare","QRS","",0,3000);
    vector<float> vectres;
    vectres.push_back(fit->Parameter(1)/(mass*mass)); //K
    vectres.push_back(fit->Parameter(0)); //C
    vectres.push_back(fit->Error(1)/(mass*mass));
    vectres.push_back(fit->Error(0));
    return vectres;
}


vector<int> CrossTalkInv(Correction corr,const std::vector<int>&  Q, const float x1, const float x2, float threshold, float thresholdSat, int label, int correctType, float RatioSat254, float RatioSat255, float ThresholdRatioSat) //Methode de correction + XTalkInv
{
  
    bool DylanCorr=true;
    bool inversion=true;
    if(correctType==1){DylanCorr=false;inversion=false;}
    if(correctType==2){DylanCorr=false;}
    if(correctType==3){inversion=false;}

    float charge=0.;

    const unsigned N=Q.size();
    std::vector<int> QII;
    std::vector<float> QI(N,0);
    Double_t a=1-2*x1-2*x2;
    //  bool debugbool=false;
    TMatrix A(N,N);

    float RatioSatTotal = RatioSat254+RatioSat255;

    vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());

    //----------------------
  
    if(!DylanCorr)
    { 
        //---  que pour 1 max bien net 
        if(Q.size()<2 || Q.size()>8)
        {
    	    for (unsigned int i=0;i<Q.size();i++)
            {
    	    	QII.push_back(Q[i]);
      	    }
            return QII;
        }
        
    	if(*mQ>253){
            
    		if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253) return Q;

    		if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40 )
            {
                QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2);
                return QII;
            }
    	}
    }

    //----------------------

    if(DylanCorr && *mQ>253 && RatioSatTotal>=ThresholdRatioSat)
    {
        int NSat254=0;
        int NSat255=0;
        float ClusterCharge=0;

        for(int i=0;i<Q.size();i++)
        {     
            if(Q[i]==254) NSat254++;
            if(Q[i]==255) NSat255++;
            ClusterCharge+=(Q[i]);
        }

        float ClusterChargeCorr=corr.ChargeCorr(ClusterCharge*(3.61*pow(10,-9)*247),label,Q.size(),NSat254,NSat255)/(3.61*pow(10,-9)*247);
        float DiffClusterCharge=ClusterChargeCorr-ClusterCharge;
        for (unsigned int i=0;i<Q.size();i++)
        {
    	    QII.push_back(Q[i]);
      	}
        //Incomplet pour le moment --> a voir au moment de ClusterCleaning
        //Le surplus de charge ne se fait pas que sur une piste 

        vector<int>::iterator maxQ = max_element(QII.begin(), QII.end());
        if(DiffClusterCharge>=0) QII.at(std::distance(QII.begin(),maxQ))+=DiffClusterCharge;
        return QII;
    }

    //----------------------

    if(inversion)
    {
        for(unsigned int i=0; i<N; i++) 
        {
            A(i,i)=a;
            if(i<N-1){ A(i+1,i)=x1;A(i,i+1)=x1;}
            else continue; 
            if(i<N-2){ A(i+2,i)=x2;A(i,i+2)=x2;}
        }
    
        if(N==1) A(0,0)=1/a;
        else  A.InvertFast();

        for(unsigned int i=0; i<N; i++) 
        {
            for(unsigned int j=0; j<N; j++) 
            {
                QI[i]+=A(i,j)*(float)Q[j];
            }
        }
    
        for (unsigned int i=0;i<QI.size();i++)
        {
            if(QI[i]<threshold) QI[i]=0; 
            QII.push_back((int) QI[i]);
        }
        return QII;
    
    }

    if(!inversion) return Q;

}







//MAIN



int main(int argc,char** argv)
{

	TChain chain;
	chain.SetName("stage/ttree");
	for(int i=1;i<argc-1;i++) chain.Add(argv[i]); 
    string s1 = argv[1];
	string s2 = s1.substr(0,s1.find('.'))+"_results.root";

    Builder* b1 = new Builder(chain);
    b1->SetBranchAdd();


	int nentries = b1->GetEntries();


//OUVERTURE FICHIERS DE CONFIGS POUR LA METHODE

	Correction EcorrLayer;
	TFile* fileLayer = TFile::Open("testMethodLayer.root");
	TTree* treeLayer = (TTree*) fileLayer->Get("tree");
	EcorrLayer.SetTree(*treeLayer);
    EcorrLayer.Read();
	Correction EcorrModulGeom;
	TFile* fileModulGeom = TFile::Open("testMethodModulGeom.root");
	TTree* treeModulGeom = (TTree*) fileModulGeom->Get("tree");
	EcorrModulGeom.SetTree(*treeModulGeom);
	Correction EcorrModulGeomNoSat;
	TFile* fileModulGeomNoSat = TFile::Open("testMethodModulGeomNoSat.root");
	TTree* treeModulGeomNoSat = (TTree*) fileModulGeomNoSat->Get("tree");
	EcorrModulGeomNoSat.SetTree(*treeModulGeomNoSat);

// OUVERTURE FICHIER CONFIG FACTEURS K ET C

    ifstream inputfileFactorKC ("factKC.txt");
    float K_Esim, C_Esim;
    float K_Erec, C_Erec;
    float K_Ecorr, C_Ecorr;
    float K_EfullLayer, C_EfullLayer;
    float K_EfullModulGeom, C_EfullModulGeom;

    inputfileFactorKC>>K_Esim; inputfileFactorKC>>C_Esim;
    inputfileFactorKC>>K_Erec; inputfileFactorKC>>C_Erec;
    inputfileFactorKC>>K_Ecorr; inputfileFactorKC>>C_Ecorr;
    inputfileFactorKC>>K_EfullLayer; inputfileFactorKC>>C_EfullLayer;
    inputfileFactorKC>>K_EfullModulGeom; inputfileFactorKC>>C_EfullModulGeom;
    



	bool method=false;
	bool corr=false;
    bool FactKC=false;
    bool XTalkInvStudy=false;

    bool gluino=false;

    float idlow=2212,idup=2212;

    float ratsat=0.;

//CHOIX OPTION



	if(atof(argv[argc-1])==1) method=true;
	if(atof(argv[argc-1])==2) corr=true;
    if(atof(argv[argc-1])==3) {corr=true;FactKC=true;}
    if(atof(argv[argc-1])==4) XTalkInvStudy=true;


    if(atof(argv[argc-2])==1) {gluino=true;idlow=1000001;idup=2000015;}

    ratsat=atof(argv[argc-3]);
    

//------------------------------------------------------
//
//
//          Le premier 1 sert a prendre les gluinos
//
//
//
//          Option 1 : creer les fichiers avec les parametres de correction 
//          Option 2 : correction de l'energie
//          Option 3 : enregistre les facteurs K et C (necessite en plus l'option 2)
//          Option 4 : CrossTalkInvAlgo study
//
//
//
//
//
//------------------------------------------------------




//DECLARATION TABLEAUX POUR LES METHODES 

vector<vector<vector<vector<TH2F*>>>> VectLayerVectNStripVectNStripSat254VectNStripSat255Histo;
vector<vector<vector<vector<TH2F*>>>> VectModulGeomVectNStripVectNStripSat254VectNStripSat255Histo;
vector<vector<vector<vector<TH2F*>>>> VectModulGeomVectNStripVectNStripSat254VectNStripSat255HistoNoSat;
if(method)
{
//METHODE 1 -- CATEGORIES LAYER
    for(int layer=1;layer<22;layer++)
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
                    string title = LabelLayer(layer)+" NStrip="+to_string(nstrip)+" NStripSat254="+to_string(nstripsat254)+" NStripSat255="+to_string(nstripsat255);
                    VectNStripSat255Histo.push_back(new TH2F(title.c_str(),title.c_str(),300,0,6000*pow(10,-6),300,0,6000*pow(10,-6)));
                }
                VectNStripSat254VectNStripSat255Histo.push_back(VectNStripSat255Histo);
            }
            VectNStripVectNStripSat254VectNStripSat255Histo.push_back(VectNStripSat254VectNStripSat255Histo);
        }
        VectLayerVectNStripVectNStripSat254VectNStripSat255Histo.push_back(VectNStripVectNStripSat254VectNStripSat255Histo);
    }

//METHODE 2 -- CATEGORIES XTALK

	for(int ModulGeom=1;ModulGeom<15;ModulGeom++)
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
                    string title = LabelModulGeom(ModulGeom)+" NStrip="+to_string(nstrip)+" NStripSat254="+to_string(nstripsat254)+" NStripSat255="+to_string(nstripsat255);
                    VectNStripSat255Histo.push_back(new TH2F(title.c_str(),title.c_str(),300,0,6000*pow(10,-6),300,0,6000*pow(10,-6)));
                }
                VectNStripSat254VectNStripSat255Histo.push_back(VectNStripSat255Histo);
            }
            VectNStripVectNStripSat254VectNStripSat255Histo.push_back(VectNStripSat254VectNStripSat255Histo);
        }
        VectModulGeomVectNStripVectNStripSat254VectNStripSat255Histo.push_back(VectNStripVectNStripSat254VectNStripSat255Histo);
    }



}


//DECLARATION HISTOGRAMMES POUR ETUDES

TH1I* hmodulgeom = new TH1I("modulgeom","",25,0,25);
TH1I* hlayerlabel = new TH1I("layerlabel","",25,0,25);

TH2I* hmodulgeomvslayer = new TH2I("hmodulgeomvslayer","",21,0,21,14,0,14);

TH1F* hEcorrLayer = new TH1F("hEcorrLayer","",100,-2,2);
TH1F* hEcorrModulGeom = new TH1F("hEcorrModulGeom","",100,-2,2);
TH1F* hEcorrModulGeomNoSat = new TH1F("hEcorrModulGeomNoSat","",100,-2,2);

TH1F* hErecLayer = new TH1F("hErecLayer","",100,-2,2);
TH1F* hErecModulGeom = new TH1F("hErecModulGeom","",100,-2,2);
TH1F* hErecModulGeomNoSat = new TH1F("hErecModulGeomNoSat","",100,-2,2);

TH2F* h2ElossEcorrSatLayer = new TH2F("h2ElossEcorrSatLayer","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH2F* h2ElossEcorrModulGeom = new TH2F("h2ElossEcorrModulGeom","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH2F* h2ElossEcorrModulGeomNoSat = new TH2F("h2ElossEcorrModulGeomNoSat","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));

TH2F* h2ElossErec = new TH2F("h2ElossErec","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));

TH2F* h2ElossErecLayer = new TH2F("h2ElossErecLayer","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH2F* h2ElossErecModulGeom = new TH2F("h2ElossErecModulGeom","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH2F* h2ElossErecModulGeomNoSat = new TH2F("h2ElossErecModulGeomNoSat","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));

TH2F* h2EcorrLayerVsEcorrModulGeom = new TH2F("h2EcorrLayerVsEcorrModulGeom","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH1F* h1EcorrLayerVsEcorrModulGeom = new TH1F("h1EcorrLayerVsEcorrModulGeom","",100,-0.1,0.1);

TH2F* h2PoverMDecorrDxLayer = new TH2F("h2PoverMDecorrDxLayer","",100,0,10,50,0,50);
TH2F* h2PoverMDecorrDxModulGeom = new TH2F("h2PoverMDecorrDxModulGeom","",100,0,10,50,0,50);

TH2F* h2PvsDecorrDxLayer = new TH2F("h2PvsDecorrDxLayer","",100,0,3000,50,0,50);
TH2F* h2PvsDecorrDxModulGeom = new TH2F("h2PvsDecorrDxModulGeom","",100,0,3000,50,0,50);

TH2F* h2PvsDelossDxProtonKC = new TH2F("h2PvsDelossDxProtonKC","",1000,0,5,1000,0,15);
TH2F* h2PvsDecorrDxProtonKC = new TH2F("h2PvsDecorrDxProtonKC","",1000,0,5,1000,0,15);
TH2F* h2PvsDecorrDxProtonKCnewSatLayer = new TH2F("h2PvsDecorrDxProtonKCnewSatLayer","",1000,0,5,1000,0,15);

TH1F* historeceloss = new TH1F("","",1000,0,5);
TH1F* historececorr = new TH1F("","",1000,0,5);
TH1F* historececorrnew = new TH1F("","",1000,0,5);


TH1F* distribmassprotonEloss = new TH1F("distribmassprotonEloss","",100,0,5000);
TH1F* distribmassprotonEcorr = new TH1F("distribmassprotonEcorr","",100,0,5000);
TH1F* distribmassprotonEcorrNew = new TH1F("distribmassprotonEcorrNew","",100,0,5000);


TH2F* pVsdEsimdx = new TH2F("pVsdEsimdx","",1000,0,5,1000,0,15);
TH2F* pVsdErecdx = new TH2F("pVsdErecdx","",1000,0,5,1000,0,15);
TH2F* pVsdEcorrdx = new TH2F("pVsdEcorrdx","",1000,0,5,1000,0,15);
TH2F* pVsdEfulldx_Layer = new TH2F("pVsdEfulldx_Layer","",1000,0,5,1000,0,15);
TH2F* pVsdEfulldx_ModulGeom = new TH2F("pVsdEfulldx_ModulGeom","",1000,0,5,1000,0,15);


float rangemass=2;
if(gluino) rangemass=5000;

TH1F* massEsim = new TH1F("massEsim","",100,0,rangemass);
TH1F* massErec = new TH1F("massErec","",100,0,rangemass);
TH1F* hmassEcorr = new TH1F("massEcorr","",100,0,rangemass);
TH1F* massEfullLayer = new TH1F("massEfullLayer","",100,0,rangemass);
TH1F* massEfullModulGeom = new TH1F("massEfullModulGeom","",100,0,rangemass);

//Distrib d'interet (saturation)

TH1F* h_ClustSat_DistribP = new TH1F("h_ClustSat_DistribP","",100,0,10000);
TH1F* h_ClustSat_DistribEta = new TH1F("h_ClustSat_DistribEta","",100,0,3);
TH1I* h_ClustSat_DistribSize = new TH1I("h_ClustSat_DistribSize","",8,0,8);
TH1I* h_ClustSat_DistribNSat = new TH1I("h_ClustSat_DistribNsat","",5,0,5);

TH1F* h_DistribP = new TH1F("h_DistribP","",100,0,10000);
TH1F* h_DistribEta = new TH1F("h_DistribEta","",100,0,3);
TH1I* h_DistribSize = new TH1I("h_DistribSize","",8,0,8);
TH1I* h_DistribNSat = new TH1I("h_DistribNsat","",5,0,5);



//Etude de la methode actuelle + XTalkInv + comparaison avec DylanCorr

TH1F* hWithoutCorr = new TH1F("hWithoutCorr","",100,-2,2);
TH1F* hCorrHSCP = new TH1F("hCorrHSCP","",100,-2,2);
TH1F* hCorrHSCP_WithoutInversion = new TH1F("hCorrHSCP_WithoutInversion","",100,-2,2);
TH1F* hMyCorr = new TH1F("hMyCorr","",100,-2,2);
TH1F* hMyCorr_WithoutInversion = new TH1F("MyCorr_WithoutInversion","",100,-2,2);
TH1F* hMyCorr_RatioSat = new TH1F("hMyCorr_RatioSat","",100,-2,2);
TH1F* hMyCorr_WithoutInversion_RatioSat = new TH1F("hMyCorr_WithoutInversion_RatioSat","",100,-2,2);

TH1F* h_ClustSat_WithoutCorr = new TH1F("h_ClustSat_WithoutCorr","",100,-2,2);
TH1F* h_ClustSat_CorrHSCP = new TH1F("h_ClustSat_CorrHSCP","",100,-2,2);
TH1F* h_ClustSat_CorrHSCP_WithoutInversion = new TH1F("h_ClustSat_CorrHSCP_WithoutInversion","",100,-2,2);
TH1F* h_ClustSat_MyCorr = new TH1F("h_ClustSat_MyCorr","",100,-2,2);
TH1F* h_ClustSat_MyCorr_WithoutInversion = new TH1F("h_ClustSat_MyCorr_WithoutInversion","",100,-2,2);
TH1F* h_ClustSat_MyCorr_RatioSat = new TH1F("h_ClustSat_MyCorr_RatioSat","",100,-2,2);
TH1F* h_ClustSat_MyCorr_WithoutInversion_RatioSat = new TH1F("h_ClustSat_MyCorr_WithoutInversion_RatioSat","",100,-2,2);

TH1F* h_ClustNoSat_WithoutCorr = new TH1F("h_ClustNoSat_WithoutCorr","",100,-2,2);
TH1F* h_ClustNoSat_CorrHSCP = new TH1F("h_ClustNoSat_CorrHSCP","",100,-2,2);
TH1F* h_ClustNoSat_CorrHSCP_WithoutInversion = new TH1F("h_ClustNoSat_CorrHSCP_WithoutInversion","",100,-2,2);
TH1F* h_ClustNoSat_MyCorr = new TH1F("h_ClustNoSat_MyCorr","",100,-2,2);
TH1F* h_ClustNoSat_MyCorr_WithoutInversion = new TH1F("h_ClustNoSat_MyCorr_WithoutInversion","",100,-2,2);
TH1F* h_ClustNoSat_MyCorr_RatioSat = new TH1F("h_ClustNoSat_MyCorr_RatioSat","",100,-2,2);
TH1F* h_ClustNoSat_MyCorr_WithoutInversion_RatioSat = new TH1F("h_ClustNoSat_MyCorr_WithoutInversion_RatioSat","",100,-2,2);

TH2F* h2_ClustSat_MyCorrWOinv_p = new TH2F("h2_ClustSat_MyCorrWOinv_p","",100,-2,2,100,0,4000);
TH2F* h2_ClustSat_MyCorrWOinv_eta = new TH2F("h2_ClustSat_MyCorrWOinv_eta","",100,-2,2,100,0,2.5);

TH1F* h_ClustSat_MyCorrWOinv_0p500 = new TH1F("h_ClustSat_MyCorrWOinv_0p500","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_500p1000 = new TH1F("h_ClustSat_MyCorrWOinv_500p1000","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_1000p1500 = new TH1F("h_ClustSat_MyCorrWOinv_1000p1500","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_1500p2000 = new TH1F("h_ClustSat_MyCorrWOinv_1500p2000","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_2000p4000 = new TH1F("h_ClustSat_MyCorrWOinv_2000p4000","",50,-2,2);

TH1F* h_ClustSat_MyCorrWOinv_0eta05 = new TH1F("h_ClustSat_MyCorrWOinv_0eta05","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_05eta10 = new TH1F("h_ClustSat_MyCorrWOinv_05eta10","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_10eta15 = new TH1F("h_ClustSat_MyCorrWOinv_10eta15","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_15eta25 = new TH1F("h_ClustSat_MyCorrWOinv_15eta25","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_20eta25 = new TH1F("h_ClustSat_MyCorrWOinv_20eta25","",50,-2,2);

TH1F* h_ClustSat_HscpCorrWOinv_0p500 = new TH1F("h_ClustSat_HscpCorrWOinv_0p500","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_500p1000 = new TH1F("h_ClustSat_HscpCorrWOinv_500p1000","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_1000p1500 = new TH1F("h_ClustSat_HscpCorrWOinv_1000p1500","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_1500p2000 = new TH1F("h_ClustSat_HscpCorrWOinv_1500p2000","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_2000p4000 = new TH1F("h_ClustSat_HscpCorrWOinv_2000p4000","",50,-2,2);

TH1F* h_ClustSat_HscpCorrWOinv_0eta05 = new TH1F("h_ClustSat_HscpCorrWOinv_0eta05","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_05eta10 = new TH1F("h_ClustSat_HscpCorrWOinv_05eta10","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_10eta15 = new TH1F("h_ClustSat_HscpCorrWOinv_10eta15","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_15eta25 = new TH1F("h_ClustSat_HscpCorrWOinv_15eta25","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_20eta25 = new TH1F("h_ClustSat_HscpCorrWOinv_20eta25","",50,-2,2);

//Selon les differentes couches 

TH1F* h_ClustSat_MyCorrWOinv_TIB = new TH1F("h_ClustSat_MyCorrWOinv_TIB","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_TID = new TH1F("h_ClustSat_MyCorrWOinv_TID","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_TOB = new TH1F("h_ClustSat_MyCorrWOinv_TOB","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_TEC = new TH1F("h_ClustSat_MyCorrWOinv_TEC","",50,-2,2);

TH1F* h_ClustSat_HscpCorrWOinv_TIB = new TH1F("h_ClustSat_HscpCorrWOinv_TIB","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_TID = new TH1F("h_ClustSat_HscpCorrWOinv_TID","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_TOB = new TH1F("h_ClustSat_HscpCorrWOinv_TOB","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_TEC = new TH1F("h_ClustSat_HscpCorrWOinv_TEC","",50,-2,2);

TH1F* h_ClustSat_MyCorrWOinv_Layer1 = new TH1F("h_ClustSat_MyCorrWOinv_Layer1","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer2 = new TH1F("h_ClustSat_MyCorrWOinv_Layer2","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer3 = new TH1F("h_ClustSat_MyCorrWOinv_Layer3","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer4 = new TH1F("h_ClustSat_MyCorrWOinv_Layer4","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer5 = new TH1F("h_ClustSat_MyCorrWOinv_Layer5","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer6 = new TH1F("h_ClustSat_MyCorrWOinv_Layer6","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer7 = new TH1F("h_ClustSat_MyCorrWOinv_Layer7","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer8 = new TH1F("h_ClustSat_MyCorrWOinv_Layer8","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer9 = new TH1F("h_ClustSat_MyCorrWOinv_Layer9","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer10 = new TH1F("h_ClustSat_MyCorrWOinv_Layer10","",50,-2,2);

TH1F* h_ClustSat_HscpCorrWOinv_Layer1 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer1","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer2 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer2","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer3 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer3","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer4 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer4","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer5 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer5","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer6 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer6","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer7 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer7","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer8 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer8","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer9 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer9","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer10 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer10","",50,-2,2);

//Selon taille clusters 

TH1F* h_ClustSat_MyCorrWOinv_size1 = new TH1F("h_ClustSat_MyCorrWOinv_size1","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size2 = new TH1F("h_ClustSat_MyCorrWOinv_size2","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size3 = new TH1F("h_ClustSat_MyCorrWOinv_size3","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size4 = new TH1F("h_ClustSat_MyCorrWOinv_size4","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size5 = new TH1F("h_ClustSat_MyCorrWOinv_size5","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size6 = new TH1F("h_ClustSat_MyCorrWOinv_size6","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size7 = new TH1F("h_ClustSat_MyCorrWOinv_size7","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size8 = new TH1F("h_ClustSat_MyCorrWOinv_size8","",50,-2,2);

TH1F* h_ClustSat_HscpCorrWOinv_size1 = new TH1F("h_ClustSat_HscpCorrWOinv_size1","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size2 = new TH1F("h_ClustSat_HscpCorrWOinv_size2","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size3 = new TH1F("h_ClustSat_HscpCorrWOinv_size3","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size4 = new TH1F("h_ClustSat_HscpCorrWOinv_size4","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size5 = new TH1F("h_ClustSat_HscpCorrWOinv_size5","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size6 = new TH1F("h_ClustSat_HscpCorrWOinv_size6","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size7 = new TH1F("h_ClustSat_HscpCorrWOinv_size7","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size8 = new TH1F("h_ClustSat_HscpCorrWOinv_size8","",50,-2,2);

//Selon nombre de strips saturees 

TH1F* h_ClustSat_MyCorrWOinv_sat1 = new TH1F("h_ClustSat_MyCorrWOinv_sat1","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_sat2 = new TH1F("h_ClustSat_MyCorrWOinv_sat2","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_sat3 = new TH1F("h_ClustSat_MyCorrWOinv_sat3","",50,-2,2);

TH1F* h_ClustSat_HscpCorrWOinv_sat1 = new TH1F("h_ClustSat_HscpCorrWOinv_sat1","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_sat2 = new TH1F("h_ClustSat_HscpCorrWOinv_sat2","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_sat3 = new TH1F("h_ClustSat_HscpCorrWOinv_sat3","",50,-2,2);

//Etude sous correction a bas p 

TH1I* h_ClustSat_LowP_ClusterSize = new TH1I("h_ClustSat_LowP_ClusterSize","",8,0,8);
TH1I* h_ClustSat_LowP_Sat254 = new TH1I("h_ClustSat_LowP_Sat254","",5,0,5);
TH1I* h_ClustSat_LowP_Sat255 = new TH1I("h_ClustSat_LowP_Sat255","",5,0,5);
TH1I* h_ClustSat_LowP_NSat = new TH1I("h_ClustSat_LowP_NSat","",5,0,5);
TH1D* h_ClustSat_LowP_layerlabel = new TH1D("h_ClustSat_LowP_layerlabel","",21,0,21);
TH1F* h_ClustSat_LowP_Eta = new TH1F("h_ClustSat_LowP_Eta","",10,0,2.5);
TH1F* h_ClustSat_LowP_Esim = new TH1F("h_ClustSat_LowP_Esim","",300,0,6000*pow(10,-6));
TH1F* h_ClustSat_LowP_EcorrHscp = new TH1F("h_ClustSat_LowP_EcorrHscp","",300,0,6000*pow(10,-6));
TH1F* h_ClustSat_LowP_Emycorr = new TH1F("h_ClustSat_LowP_Emycorr","",300,0,6000*pow(10,-6));

//Etude seuil 1 

TH1F* h_ClustSat_HscpCorrWinv_firstThreshold20 = new TH1F("h_ClustSat_HscpCorrWinv_firstThreshold20","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstThreshold50 = new TH1F("h_ClustSat_HscpCorrWinv_firstThreshold50","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstThreshold100 = new TH1F("h_ClustSat_HscpCorrWinv_firstThreshold100","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstThreshold150 = new TH1F("h_ClustSat_HscpCorrWinv_firstThreshold150","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstThreshold200 = new TH1F("h_ClustSat_HscpCorrWinv_firstThreshold200","",50,-2,2);

//Etude seuil 2

TH1F* h_ClustSat_HscpCorrWinv_secondThreshold25 = new TH1F("h_ClustSat_HscpCorrWinv_secondThreshold25","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondThreshold50 = new TH1F("h_ClustSat_HscpCorrWinv_secondThreshold50","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondThreshold100 = new TH1F("h_ClustSat_HscpCorrWinv_secondThreshold100","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondThreshold150 = new TH1F("h_ClustSat_HscpCorrWinv_secondThreshold150","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondThreshold200 = new TH1F("h_ClustSat_HscpCorrWinv_secondThreshold200","",50,-2,2);

//Etude x-talk1

TH1F* h_ClustSat_HscpCorrWinv_firstXTalk10 = new TH1F("h_ClustSat_HscpCorrWinv_firstXTalk10","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstXTalk11 = new TH1F("h_ClustSat_HscpCorrWinv_firstXTalk11","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstXTalk12 = new TH1F("h_ClustSat_HscpCorrWinv_firstXTalk12","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstXTalk9 = new TH1F("h_ClustSat_HscpCorrWinv_firstXTalk9","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstXTalk8 = new TH1F("h_ClustSat_HscpCorrWinv_firstXTalk8","",50,-2,2);

//Etude x-talk2

TH1F* h_ClustSat_HscpCorrWinv_secondXTalk4 = new TH1F("h_ClustSat_HscpCorrWinv_secondXTalk4","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondXTalk5 = new TH1F("h_ClustSat_HscpCorrWinv_secondXTalk5","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondXTalk6 = new TH1F("h_ClustSat_HscpCorrWinv_secondXTalk6","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondXTalk3 = new TH1F("h_ClustSat_HscpCorrWinv_secondXTalk3","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondXTalk2 = new TH1F("h_ClustSat_HscpCorrWinv_secondXTalk2","",50,-2,2);

//Etude Ih

TH2F* h2_IhHarmonic_sim = new TH2F("h2_IhHarmonic_sim","",100,0,3,50,0,50);
TH2F* h2_IhHarmonic_MyCorr = new TH2F("h2_IhHarmonic_MyCorr","",100,0,3,50,0,50);
TH2F* h2_IhHarmonic_HscpCorrWOinv = new TH2F("h2_IhHarmonic_HscpCorrWOinv","",100,0,3,50,0,50);
TH2F* h2_IhHarmonic_HscpCorr = new TH2F("h2_IhHarmonic_HscpCorr","",100,0,3,50,0,50);





//COUPURES SUR LES EVENTS


b1->SetThresholdPartId(0.7); 
b1->SetThresholdPt(65); 
b1->SetThresholdP(0);
b1->SetThresholdEta(5);




//BOUCLE SUR LES EVENTS

int entries = nentries;

for(int i=0;i<entries;i++)
{
    if(i%1000==0) cout<<"Event "<<i<<endl;
    b1->GetEntry(i);

	//BOUCLE SUR LES TRACES

    for(int track=0;track<b1->GetNtracks();track++)
    {
        vector<int> vect_partID;
        vector<float> vect_Eloss;
        vector<float> vect_charge;
        vector<float> vect_dedx;
        vector<float> vect_dqdx;
        vector<float> vect_dqcorrdx;
        vector<float> vect_dqcorrdx2;
        vector<float> vect_plength;

        vector<float> vect_decorrdxnewSatLayer;
        vector<float> vect_decorrdx;

        vector<float> dEsimdx;
        vector<float> dErecdx;
        vector<float> dEcorrdx;
        vector<float> dEfulldx_Layer;
        vector<float> dEfulldx_ModulGeom;



        vector<float> VectDecorrDxLayer;
        vector<float> VectDecorrDxModulGeom;
        float pt                = b1->GetVectTrack()[track].GetPt();
        float p                 = b1->GetVectTrack()[track].GetP();
        float eta               = b1->GetVectTrack()[track].GetEta();
        float phi               = b1->GetVectTrack()[track].GetPhi();
        float ias_ampl          = b1->GetVectTrack()[track].GetIasAmpl();
        int NCluster            = b1->GetVectTrack()[track].GetNCluster();
        int NClustSat254        = b1->GetVectTrack()[track].GetNSatCluster(254);
        int NClustSat255        = b1->GetVectTrack()[track].GetNSatCluster(255);
        float RatioNClusterSat254 = (double)NClustSat254/(double)NCluster;
        float RatioNClusterSat255 = (double)NClustSat255/(double)NCluster;
        int id                  = b1->GetVectTrack()[track].GetPartId();
        float PoverM            = GetPoverM(p,id);


        std::vector<float> dEdxsim;
        std::vector<float> dEdxMyCorr;
        std::vector<float> dEdxHscpCorrWOinv;
        std::vector<float> dEdxHscpCorr;


    if(idlow<=abs(id) && abs(id)<=idup)

    {
		//BOUCLE SUR LES CLUSTERS

        for(int cluster=0;cluster<b1->GetVectTrack()[track].GetNCluster();cluster++)
        {
            float charge        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSclusCharge();
            float charge_corr   = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSclusChargeCorr();
            float EcorrBef      = charge_corr*(3.61*pow(10,-9)*247);
			float Erec			= charge*(3.61*pow(10,-9)*247);
			float chargeNoSat	= b1->GetVectTrack()[track].GetVectClusters()[cluster].GetChargeWithoutSaturation();
			float ErecNoSat		= chargeNoSat*(3.61*pow(10,-9)*247);
            float Eloss         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetEloss();
            float pathlength    = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetPathLength();
            int layer           = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetLayer();
            int layerLabel      = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetLayerLabel();
            bool sat254         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSat254();
            bool sat255         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSat255();
            bool shape          = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetShape();
            int nsimhits        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSimHits();
            int nstrips         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNStrip();
            int nsat254         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStrip(254);
            int nsat255         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStrip(255);
            bool edge           = b1->GetVectTrack()[track].GetVectClusters()[cluster].Edge();
            bool cut            = b1->GetVectTrack()[track].GetVectClusters()[cluster].Cut();
            int detid           = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetDetId();
            int subdetid        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSubDetId();
            int modulgeom		= b1->GetVectTrack()[track].GetVectClusters()[cluster].GetModulGeom();
            int maxstrip        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetMaxStrip();
			float EcorrSatLayer	= Erec;
            float EcorrSatModul = Erec;
			float EcorrNoSat	= ErecNoSat;

			h2ElossErec->Fill(Eloss,Erec);


			hmodulgeom->Fill(modulgeom);
			hlayerlabel->Fill(layerLabel);
			hmodulgeomvslayer->Fill(layerLabel-1,modulgeom-1);

            vect_dedx.push_back(Eloss/pathlength);

            //BOUCLE SUR LES STRIPS D'UN CLUSTER

            int ChargeCluster=0;
            vector<int> VectChargeStrip;
            for(int strip=0;strip<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNStrip();strip++)
            {
                VectChargeStrip.push_back(b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip].GetAmpl());
                ChargeCluster+=b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip].GetAmpl();
            }
            
            //--------------------------------

            if(XTalkInvStudy)
            {
            
                

                std::vector<int> VectChargeCorrHSCP_WithoutInversion = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,1,RatioNClusterSat254,RatioNClusterSat255,0);
                std::vector<int> VectChargeCorrHSCP = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                std::vector<int> VectChargeMyCorr_WithoutInversion = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,0);
                std::vector<int> VectChargeMyCorr = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,4,RatioNClusterSat254,RatioNClusterSat255,0);
                std::vector<int> VectChargeMyCorr_WithoutInversion_RatioSat = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,ratsat);
                std::vector<int> VectChargeMyCorr_RatioSat = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,4,RatioNClusterSat254,RatioNClusterSat255,ratsat);
                

                float ChargeCorrHSCP_WithoutInversion=0.;
                float ChargeCorrHSCP=0.;
                float ChargeMyCorr_WithoutInversion=0.;
                float ChargeMyCorr=0.;
                float ChargeMyCorr_WithoutInversion_RatioSat=0.;
                float ChargeMyCorr_RatioSat=0.;

                for(int i=0;i<VectChargeCorrHSCP_WithoutInversion.size();i++) ChargeCorrHSCP_WithoutInversion+=VectChargeCorrHSCP_WithoutInversion.at(i);
                for(int i=0;i<VectChargeCorrHSCP.size();i++) ChargeCorrHSCP+=VectChargeCorrHSCP.at(i);
                for(int i=0;i<VectChargeMyCorr_WithoutInversion.size();i++) ChargeMyCorr_WithoutInversion+=VectChargeMyCorr_WithoutInversion.at(i);
                for(int i=0;i<VectChargeMyCorr.size();i++) ChargeMyCorr+=VectChargeMyCorr.at(i);
                for(int i=0;i<VectChargeMyCorr_WithoutInversion_RatioSat.size();i++) ChargeMyCorr_WithoutInversion_RatioSat+=VectChargeMyCorr_WithoutInversion_RatioSat.at(i);
                for(int i=0;i<VectChargeMyCorr_RatioSat.size();i++) ChargeMyCorr_RatioSat+=VectChargeMyCorr_RatioSat.at(i);

                float ECorrHSCP_WithoutInversion=(3.61*pow(10,-9)*247)*ChargeCorrHSCP_WithoutInversion;
                float ECorrHSCP=(3.61*pow(10,-9)*247)*ChargeCorrHSCP;
                //float EMyCorr_WithoutInversion=(3.61*pow(10,-9)*247)*ChargeMyCorr_WithoutInversion;
                float EMyCorr=(3.61*pow(10,-9)*247)*ChargeMyCorr;
                
                float EMyCorr_WithoutInversion=(3.61*pow(10,-9)*247)*ChargeMyCorr_WithoutInversion_RatioSat;
                
                float EMyCorr_WithoutInversion_RatioSat=(3.61*pow(10,-9)*247)*ChargeMyCorr_WithoutInversion_RatioSat;
                float EMyCorr_RatioSat=(3.61*pow(10,-9)*247)*ChargeMyCorr_RatioSat;


                hWithoutCorr->Fill((Erec-Eloss)/Eloss);
                hCorrHSCP_WithoutInversion->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                hCorrHSCP->Fill((ECorrHSCP-Eloss)/Eloss);
                hMyCorr_WithoutInversion->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                hMyCorr->Fill((EMyCorr-Eloss)/Eloss);
                hMyCorr_WithoutInversion_RatioSat->Fill((EMyCorr_WithoutInversion_RatioSat-Eloss)/Eloss);
                hMyCorr_RatioSat->Fill((EMyCorr_RatioSat-Eloss)/Eloss);     

                h_DistribP->Fill(p);
                h_DistribEta->Fill(abs(eta));
                h_DistribSize->Fill(nstrips);
                h_DistribNSat->Fill(nsat254+nsat255);

                //First Threshold study

                std::vector<int> VectChargeCorrHSCPfirstThreshold20 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstThreshold20=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstThreshold20.size();i++) ChargeCorrHSCPfirstThreshold20+=VectChargeCorrHSCPfirstThreshold20.at(i);
                float ECorrHSCPfirstThreshold20=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstThreshold20;

                std::vector<int> VectChargeCorrHSCPfirstThreshold50 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,50,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstThreshold50=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstThreshold50.size();i++) ChargeCorrHSCPfirstThreshold50+=VectChargeCorrHSCPfirstThreshold50.at(i);
                float ECorrHSCPfirstThreshold50=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstThreshold50;

                std::vector<int> VectChargeCorrHSCPfirstThreshold100 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,100,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstThreshold100=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstThreshold100.size();i++) ChargeCorrHSCPfirstThreshold100+=VectChargeCorrHSCPfirstThreshold100.at(i);
                float ECorrHSCPfirstThreshold100=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstThreshold100;

                std::vector<int> VectChargeCorrHSCPfirstThreshold150 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,150,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstThreshold150=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstThreshold150.size();i++) ChargeCorrHSCPfirstThreshold150+=VectChargeCorrHSCPfirstThreshold150.at(i);
                float ECorrHSCPfirstThreshold150=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstThreshold150;

                std::vector<int> VectChargeCorrHSCPfirstThreshold200 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,200,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstThreshold200=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstThreshold200.size();i++) ChargeCorrHSCPfirstThreshold200+=VectChargeCorrHSCPfirstThreshold200.at(i);
                float ECorrHSCPfirstThreshold200=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstThreshold200;

                //Second Threshold Study

                std::vector<int> VectChargeCorrHSCPsecondThreshold25 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondThreshold25=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondThreshold25.size();i++) ChargeCorrHSCPsecondThreshold25+=VectChargeCorrHSCPsecondThreshold25.at(i);
                float ECorrHSCPsecondThreshold25=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondThreshold25;

                std::vector<int> VectChargeCorrHSCPsecondThreshold50 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,50,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondThreshold50=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondThreshold50.size();i++) ChargeCorrHSCPsecondThreshold50+=VectChargeCorrHSCPsecondThreshold50.at(i);
                float ECorrHSCPsecondThreshold50=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondThreshold50;

                std::vector<int> VectChargeCorrHSCPsecondThreshold100 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,100,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondThreshold100=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondThreshold100.size();i++) ChargeCorrHSCPsecondThreshold100+=VectChargeCorrHSCPsecondThreshold100.at(i);
                float ECorrHSCPsecondThreshold100=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondThreshold100;

                std::vector<int> VectChargeCorrHSCPsecondThreshold150 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,150,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondThreshold150=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondThreshold150.size();i++) ChargeCorrHSCPsecondThreshold150+=VectChargeCorrHSCPsecondThreshold150.at(i);
                float ECorrHSCPsecondThreshold150=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondThreshold150;

                std::vector<int> VectChargeCorrHSCPsecondThreshold200 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,200,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondThreshold200=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondThreshold200.size();i++) ChargeCorrHSCPsecondThreshold200+=VectChargeCorrHSCPsecondThreshold200.at(i);
                float ECorrHSCPsecondThreshold200=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondThreshold200;

                //First XTalk study 

                std::vector<int> VectChargeCorrHSCPfirstXTalk10 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstXTalk10=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstXTalk10.size();i++) ChargeCorrHSCPfirstXTalk10+=VectChargeCorrHSCPfirstXTalk10.at(i);
                float ECorrHSCPfirstXTalk10=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstXTalk10;

                std::vector<int> VectChargeCorrHSCPfirstXTalk11 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.11,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstXTalk11=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstXTalk11.size();i++) ChargeCorrHSCPfirstXTalk11+=VectChargeCorrHSCPfirstXTalk11.at(i);
                float ECorrHSCPfirstXTalk11=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstXTalk11;

                std::vector<int> VectChargeCorrHSCPfirstXTalk12 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.12,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstXTalk12=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstXTalk12.size();i++) ChargeCorrHSCPfirstXTalk12+=VectChargeCorrHSCPfirstXTalk12.at(i);
                float ECorrHSCPfirstXTalk12=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstXTalk12;

                std::vector<int> VectChargeCorrHSCPfirstXTalk9 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.09,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstXTalk9=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstXTalk9.size();i++) ChargeCorrHSCPfirstXTalk9+=VectChargeCorrHSCPfirstXTalk9.at(i);
                float ECorrHSCPfirstXTalk9=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstXTalk9;

                std::vector<int> VectChargeCorrHSCPfirstXTalk8 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.08,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstXTalk8=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstXTalk8.size();i++) ChargeCorrHSCPfirstXTalk8+=VectChargeCorrHSCPfirstXTalk8.at(i);
                float ECorrHSCPfirstXTalk8=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstXTalk8;


                //second XTalk study 

                std::vector<int> VectChargeCorrHSCPsecondXTalk4 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondXTalk4=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondXTalk4.size();i++) ChargeCorrHSCPsecondXTalk4+=VectChargeCorrHSCPsecondXTalk4.at(i);
                float ECorrHSCPsecondXTalk4=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondXTalk4;

                std::vector<int> VectChargeCorrHSCPsecondXTalk5 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.05,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondXTalk5=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondXTalk5.size();i++) ChargeCorrHSCPsecondXTalk5+=VectChargeCorrHSCPsecondXTalk5.at(i);
                float ECorrHSCPsecondXTalk5=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondXTalk5;

                std::vector<int> VectChargeCorrHSCPsecondXTalk6 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.06,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondXTalk6=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondXTalk6.size();i++) ChargeCorrHSCPsecondXTalk6+=VectChargeCorrHSCPsecondXTalk6.at(i);
                float ECorrHSCPsecondXTalk6=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondXTalk6;

                std::vector<int> VectChargeCorrHSCPsecondXTalk3 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.03,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondXTalk3=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondXTalk3.size();i++) ChargeCorrHSCPsecondXTalk3+=VectChargeCorrHSCPsecondXTalk3.at(i);
                float ECorrHSCPsecondXTalk3=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondXTalk3;

                std::vector<int> VectChargeCorrHSCPsecondXTalk2 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.02,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondXTalk2=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondXTalk2.size();i++) ChargeCorrHSCPsecondXTalk2+=VectChargeCorrHSCPsecondXTalk2.at(i);
                float ECorrHSCPsecondXTalk2=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondXTalk2;


                if(sat254 || sat255)
                {
                    h_ClustSat_HscpCorrWinv_firstThreshold20->Fill((ECorrHSCPfirstThreshold20-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstThreshold50->Fill((ECorrHSCPfirstThreshold50-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstThreshold100->Fill((ECorrHSCPfirstThreshold100-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstThreshold150->Fill((ECorrHSCPfirstThreshold150-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstThreshold200->Fill((ECorrHSCPfirstThreshold200-Eloss)/Eloss);

                    h_ClustSat_HscpCorrWinv_secondThreshold25->Fill((ECorrHSCPsecondThreshold25-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondThreshold50->Fill((ECorrHSCPsecondThreshold50-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondThreshold100->Fill((ECorrHSCPsecondThreshold100-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondThreshold150->Fill((ECorrHSCPsecondThreshold150-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondThreshold200->Fill((ECorrHSCPsecondThreshold200-Eloss)/Eloss);

                    h_ClustSat_HscpCorrWinv_firstXTalk10->Fill((ECorrHSCPfirstXTalk10-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstXTalk11->Fill((ECorrHSCPfirstXTalk11-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstXTalk12->Fill((ECorrHSCPfirstXTalk12-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstXTalk9->Fill((ECorrHSCPfirstXTalk9-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstXTalk8->Fill((ECorrHSCPfirstXTalk8-Eloss)/Eloss);

                    h_ClustSat_HscpCorrWinv_secondXTalk4->Fill((ECorrHSCPsecondXTalk4-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondXTalk5->Fill((ECorrHSCPsecondXTalk5-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondXTalk6->Fill((ECorrHSCPsecondXTalk6-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondXTalk3->Fill((ECorrHSCPsecondXTalk3-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondXTalk2->Fill((ECorrHSCPsecondXTalk2-Eloss)/Eloss);


                    h_ClustSat_WithoutCorr->Fill((Erec-Eloss)/Eloss);
                    h_ClustSat_CorrHSCP_WithoutInversion->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    h_ClustSat_CorrHSCP->Fill((ECorrHSCP-Eloss)/Eloss);            
                    h_ClustSat_MyCorr_WithoutInversion->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    h_ClustSat_MyCorr->Fill((EMyCorr-Eloss)/Eloss);
                    h_ClustSat_MyCorr_WithoutInversion_RatioSat->Fill((EMyCorr_WithoutInversion_RatioSat-Eloss)/Eloss);
                    h_ClustSat_MyCorr_RatioSat->Fill((EMyCorr_RatioSat-Eloss)/Eloss);

                    h2_ClustSat_MyCorrWOinv_p->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss,p);
                    h2_ClustSat_MyCorrWOinv_eta->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss,abs(eta));

                    h_ClustSat_DistribP->Fill(p);
                    h_ClustSat_DistribEta->Fill(abs(eta));
                    h_ClustSat_DistribSize->Fill(nstrips);
                    h_ClustSat_DistribNSat->Fill(nsat254+nsat255);

                    if(p<=500) h_ClustSat_MyCorrWOinv_0p500->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(p<=1000 && p>500) h_ClustSat_MyCorrWOinv_500p1000->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(p<=1500 && p>1000) h_ClustSat_MyCorrWOinv_1000p1500->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(p<=2000 && p>1500) h_ClustSat_MyCorrWOinv_1500p2000->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(p<=4000 && p>2000) h_ClustSat_MyCorrWOinv_2000p4000->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);

                    if(abs(eta)<=0.5) h_ClustSat_MyCorrWOinv_0eta05->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=1 && abs(eta)>0.5) h_ClustSat_MyCorrWOinv_05eta10->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=1.5 && abs(eta)>1) h_ClustSat_MyCorrWOinv_10eta15->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=2.5 && abs(eta)>1.5) h_ClustSat_MyCorrWOinv_15eta25->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=2.5 && abs(eta)>2) h_ClustSat_MyCorrWOinv_20eta25->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);

                    if(p<=500) h_ClustSat_HscpCorrWOinv_0p500->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(p<=1000 && p>500) h_ClustSat_HscpCorrWOinv_500p1000->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(p<=1500 && p>1000) h_ClustSat_HscpCorrWOinv_1000p1500->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(p<=2000 && p>1500) h_ClustSat_HscpCorrWOinv_1500p2000->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(p<=4000 && p>2000) h_ClustSat_HscpCorrWOinv_2000p4000->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);

                    if(abs(eta)<=0.5) h_ClustSat_HscpCorrWOinv_0eta05->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=1 && abs(eta)>0.5) h_ClustSat_HscpCorrWOinv_05eta10->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=1.5 && abs(eta)>1) h_ClustSat_HscpCorrWOinv_10eta15->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=2.5 && abs(eta)>1.5) h_ClustSat_HscpCorrWOinv_15eta25->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=2.5 && abs(eta)>2) h_ClustSat_HscpCorrWOinv_20eta25->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);

                    if(p<=500 && (EMyCorr_WithoutInversion-Eloss)/Eloss<-0.25)
                    {
                        h_ClustSat_LowP_ClusterSize->Fill(nstrips);
                        h_ClustSat_LowP_Sat254->Fill(nsat254);
                        h_ClustSat_LowP_Sat255->Fill(nsat255);
                        h_ClustSat_LowP_NSat->Fill(nsat254+nsat255);
                        h_ClustSat_LowP_layerlabel->Fill(layerLabel);
                        h_ClustSat_LowP_Eta->Fill(abs(eta));
                        h_ClustSat_LowP_Esim->Fill(Eloss);
                        h_ClustSat_LowP_Emycorr->Fill(EMyCorr_WithoutInversion);
                        h_ClustSat_LowP_EcorrHscp->Fill(ECorrHSCP_WithoutInversion);
                    }

                    if(subdetid==3)
                    {
                        h_ClustSat_MyCorrWOinv_TIB->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_TIB->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(subdetid==4)
                    {
                        h_ClustSat_MyCorrWOinv_TID->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_TID->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(subdetid==5)
                    {
                        h_ClustSat_MyCorrWOinv_TOB->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_TOB->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(subdetid==6)
                    {
                        h_ClustSat_MyCorrWOinv_TEC->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_TEC->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==1) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer1->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer1->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==2) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer2->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer2->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==3) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer3->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer3->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==4) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer4->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer4->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==5) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer5->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer5->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==6) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer6->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer6->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==7) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer7->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer7->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==8) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer8->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer8->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==9) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer9->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer9->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==10) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer10->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer10->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nstrips==1)
                    {
                        h_ClustSat_MyCorrWOinv_size1->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_size1->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nstrips==2)
                    {
                        h_ClustSat_MyCorrWOinv_size2->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_size2->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nstrips==3)
                    {
                        h_ClustSat_MyCorrWOinv_size3->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_size3->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nstrips==4)
                    {
                        h_ClustSat_MyCorrWOinv_size4->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_size4->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nstrips==5)
                    {
                        h_ClustSat_MyCorrWOinv_size5->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_size5->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nstrips==6)
                    {
                        h_ClustSat_MyCorrWOinv_size6->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_size6->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nsat254+nsat255==1)
                    {
                        h_ClustSat_MyCorrWOinv_sat1->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_sat1->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nsat254+nsat255==2)
                    {
                        h_ClustSat_MyCorrWOinv_sat2->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_sat2->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nsat254+nsat255==3)
                    {
                        h_ClustSat_MyCorrWOinv_sat3->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_sat3->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }

                }
                if(!(sat254 || sat255))
                {
                    h_ClustNoSat_WithoutCorr->Fill((Erec-Eloss)/Eloss);
                    h_ClustNoSat_CorrHSCP_WithoutInversion->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    h_ClustNoSat_CorrHSCP->Fill((ECorrHSCP-Eloss)/Eloss);            
                    h_ClustNoSat_MyCorr_WithoutInversion->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    h_ClustNoSat_MyCorr->Fill((EMyCorr-Eloss)/Eloss);
                    h_ClustNoSat_MyCorr_WithoutInversion_RatioSat->Fill((EMyCorr_WithoutInversion_RatioSat-Eloss)/Eloss);
                    h_ClustNoSat_MyCorr_RatioSat->Fill((EMyCorr_RatioSat-Eloss)/Eloss);
                }




                //Ih study

                dEdxsim.push_back(Eloss/pathlength);
                dEdxMyCorr.push_back(EMyCorr_WithoutInversion/pathlength);
                dEdxHscpCorrWOinv.push_back(ECorrHSCP_WithoutInversion/pathlength);
                dEdxHscpCorr.push_back(ECorrHSCP/pathlength);

            }

            

            if(corr)
            {
                if(layerLabel<=22) dEsimdx.push_back(Eloss/pathlength);
                dErecdx.push_back(Erec/pathlength);
                dEcorrdx.push_back(EcorrBef/pathlength);
                if(nstrips>=3)
                {
                    if((RatioNClusterSat254+RatioNClusterSat255>=ratsat)) 
                    {
                        EcorrSatLayer = EcorrLayer.ChargeCorr(Erec,layerLabel,nstrips,nsat254,nsat255);
                        dEfulldx_Layer.push_back(EcorrSatLayer/pathlength);
                        EcorrSatModul = EcorrModulGeom.ChargeCorr(Erec,modulgeom,nstrips,nsat254,nsat255);
                        dEfulldx_ModulGeom.push_back(EcorrSatModul/pathlength);



                    }
                    else
                    {
                        if(!(sat254 || sat255))
                        {
                            dEfulldx_Layer.push_back(Erec/pathlength);
                            dEfulldx_ModulGeom.push_back(Erec/pathlength);
                        }
                    }
                }
            }
            

            if((sat254 || sat255))
            {
                hErecLayer->Fill((Erec-Eloss)/Eloss);
                //EcorrSatLayer = EcorrLayer.ChargeCorr(Erec,layerLabel,nstrips,nsat254,nsat255);
                hEcorrLayer->Fill((EcorrSatLayer-Eloss)/Eloss);
                h2ElossEcorrSatLayer->Fill(Eloss,EcorrSatLayer);
                hEcorrModulGeom->Fill((EcorrSatModul-Eloss)/Eloss);
                h2ElossEcorrModulGeom->Fill(Eloss,EcorrSatModul);

            }




			if(method)
			{
				//METHODE LAYER

				for(int countlayer=1;countlayer<22;countlayer++)
            	{
            	    for(int countnstrip=3;countnstrip<7;countnstrip++)
            	    {
            	        for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
            	        {
            	            for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
            	            {
            	                if(countlayer==layerLabel && countnstrip==nstrips && countnstripsat254==nsat254 && countnstripsat255==nsat255)
            	                {
            	                    VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(Eloss,Erec);
            	                } 
            	            }
            	        }
            	    }
				}

				//METHODE MODULGEOM

				for(int countlayer=1;countlayer<15;countlayer++)
            	{
            	    for(int countnstrip=3;countnstrip<7;countnstrip++)
            	    {
            	        for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
            	        {
            	            for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
            	            {
            	                if(countlayer==modulgeom && countnstrip==nstrips && countnstripsat254==nsat254 && countnstripsat255==nsat255)
            	                {
            	                    VectModulGeomVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(Eloss,Erec);
            	                } 
            	            }
            	        }
            	    }
				}

				

			}
            




		}

        Estimator estimLayer(VectDecorrDxLayer);
        Estimator estimModulGeom(VectDecorrDxModulGeom);
        Estimator estimDelossDx(vect_dedx);

        h2PoverMDecorrDxLayer->Fill(p/2400,estimLayer.GetHarmonic2()*pow(10,3));
        h2PoverMDecorrDxModulGeom->Fill(p/2400,estimModulGeom.GetHarmonic2()*pow(10,3));

        h2PvsDecorrDxLayer->Fill(p,estimLayer.GetHarmonic2()*pow(10,3));
        h2PvsDecorrDxModulGeom->Fill(p,estimModulGeom.GetHarmonic2()*pow(10,3));

        float K = KsusyCorr;
        float C = CsusyCorr;
        //float dedx = estimLayer.GetHarmonic2()*pow(10,3);
        float dedx = estimDelossDx.GetHarmonic2()*pow(10,3);

        float mass_dedx = sqrt((dedx-C)*pow(p,2)/K);

        if(dedx>0) h2PvsDelossDxProtonKC->Fill(p,dedx);




        Estimator estimDecorrDx(vect_decorrdx);
        Estimator estimDecorrDxnewSatLayer(vect_decorrdxnewSatLayer);
        float decorrdx = estimDecorrDx.GetHarmonic2()*pow(10,3);
        float decorrdxnewsatlayer = estimDecorrDxnewSatLayer.GetHarmonic2()*pow(10,3);

        if(decorrdx>0) h2PvsDecorrDxProtonKC->Fill(p,decorrdx);
        if(decorrdxnewsatlayer>0) h2PvsDecorrDxProtonKCnewSatLayer->Fill(p,decorrdxnewsatlayer);





        float massEloss = sqrt((dedx-CprotonEloss)*pow(p,2)/KprotonEloss);
        float massEcorr = sqrt((decorrdx-CprotonEcorr)*pow(p,2)/KprotonEcorr);
        float massEcorrNew = sqrt((decorrdxnewsatlayer-CprotonEcorrNew)*pow(p,2)/KprotonEcorrNew);
        //float massEloss = sqrt((dedx-CprotonEloss)*pow(p,2)/KprotonEloss);
        //float massEcorr = sqrt((decorrdx-CprotonEcorr)*pow(p,2)/KprotonEcorr);
        //float massEcorrNew = sqrt((decorrdxnewsatlayer-CprotonEcorrNew)*pow(p,2)/KprotonEcorrNew);
        
        //if(dedx<=CprotonEloss && dedx>=0) massEloss=0.;
        //if(decorrdx<=CprotonEcorr && decorrdx>=0) massEcorr=0.;
        //if(decorrdxnewsatlayer<=CprotonEcorrNew && decorrdxnewsatlayer>=0) massEcorrNew=0.;
        distribmassprotonEloss->Fill(massEloss);
        distribmassprotonEcorr->Fill(massEcorr);
        distribmassprotonEcorrNew->Fill(massEcorrNew);



        float CsusySim = 4.53; //pour gluino
        float KsusySim = 1.76;




        Estimator estimDeDx(dEsimdx);
        pVsdEsimdx->Fill(p,estimDeDx.GetHarmonic2()*pow(10,3));
        //massEsim->Fill(sqrt((estimDeDx.GetHarmonic2()*pow(10,3)-C_Esim)*pow(p,2)/K_Esim));
        massEsim->Fill(sqrt((estimDeDx.GetHarmonic2()*pow(10,3)-CsusySim)*pow(p,2)/KsusySim));

        estimDeDx.SetVect(dErecdx);
        pVsdErecdx->Fill(p,estimDeDx.GetHarmonic2()*pow(10,3));
        massErec->Fill(sqrt((estimDeDx.GetHarmonic2()*pow(10,3)-C_Erec)*pow(p,2)/K_Erec));

        estimDeDx.SetVect(dEcorrdx);
        pVsdEcorrdx->Fill(p,estimDeDx.GetHarmonic2()*pow(10,3));
        hmassEcorr->Fill(sqrt((estimDeDx.GetHarmonic2()*pow(10,3)-C_Ecorr)*pow(p,2)/K_Ecorr));

        estimDeDx.SetVect(dEfulldx_Layer);
        pVsdEfulldx_Layer->Fill(p,estimDeDx.GetHarmonic2()*pow(10,3));
        massEfullLayer->Fill(sqrt((estimDeDx.GetHarmonic2()*pow(10,3)-C_EfullLayer)*pow(p,2)/K_EfullLayer));

        estimDeDx.SetVect(dEfulldx_ModulGeom);
        pVsdEfulldx_ModulGeom->Fill(p,estimDeDx.GetHarmonic2()*pow(10,3));
        massEfullModulGeom->Fill(sqrt((estimDeDx.GetHarmonic2()*pow(10,3)-C_EfullModulGeom)*pow(p,2)/K_EfullModulGeom));




        //Study Ih

        Estimator estimdEdxsim(dEdxsim);
        Estimator estimdEdxMyCorr(dEdxMyCorr);
        Estimator estimdEdxHscpCorrWOinv(dEdxHscpCorrWOinv);
        Estimator estimdEdxHscpCorr(dEdxHscpCorr);


        h2_IhHarmonic_sim->Fill(p/MGluino,estimdEdxsim.GetHarmonic2()*pow(10,3));
        h2_IhHarmonic_MyCorr->Fill(p/MGluino,estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
        h2_IhHarmonic_HscpCorrWOinv->Fill(p/MGluino,estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));
        h2_IhHarmonic_HscpCorr->Fill(p/MGluino,estimdEdxHscpCorr.GetHarmonic2()*pow(10,3));





    }

	}
} //END LOOP ENTRIES

if(FactKC)
{
    vector<float> VectKCEsim = FactorKC(pVsdEsimdx,0.938);
    cout<<VectKCEsim[0]<<" "<<VectKCEsim[1]<<endl;

    vector<float> VectKCErec = FactorKC(pVsdErecdx,0.938);
    cout<<VectKCErec[0]<<" "<<VectKCErec[1]<<endl;

    vector<float> VectKCEcorr = FactorKC(pVsdEcorrdx,0.938);
    cout<<VectKCEcorr[0]<<" "<<VectKCEcorr[1]<<endl;

    vector<float> VectKCEfull_Layer = FactorKC(pVsdEfulldx_Layer,0.938);
    cout<<VectKCEfull_Layer[0]<<" "<<VectKCEfull_Layer[1]<<endl;

    vector<float> VectKCEfull_ModulGeom = FactorKC(pVsdEfulldx_ModulGeom,0.938);
    cout<<VectKCEfull_ModulGeom[0]<<" "<<VectKCEfull_ModulGeom[1]<<endl;

    ofstream outputfileFactorKC ("factKC.txt");
    outputfileFactorKC<<VectKCEsim[0]<<"\t"<<VectKCEsim[1]<<endl;
    outputfileFactorKC<<VectKCErec[0]<<"\t"<<VectKCErec[1]<<endl;
    outputfileFactorKC<<VectKCEcorr[0]<<"\t"<<VectKCEcorr[1]<<endl;
    outputfileFactorKC<<VectKCEfull_Layer[0]<<"\t"<<VectKCEfull_Layer[1]<<endl;
    outputfileFactorKC<<VectKCEfull_ModulGeom[0]<<"\t"<<VectKCEfull_ModulGeom[1]<<endl;
}




//ECRITURE PARAMETRES METHODE


if(method)
{
//METHODE LAYER 

	Correction MethodLayer; 
    MethodLayer.SetRebin(10);
    //MethodLayer.SetRange(0.7);
	MethodLayer.SetFileAndTreeName("testMethodLayer.root","tree");
	MethodLayer.SetBranch();
	for(int countlayer=1;countlayer<22;countlayer++)
    {
        for(int countnstrip=3;countnstrip<7;countnstrip++)
        {
        	for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
            {
                for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
				{
                    MethodLayer.SetHisto(*VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]);
                    MethodLayer.FillHisto(countlayer,countnstrip,countnstripsat254,countnstripsat255);
                    MethodLayer.FillProfile();
                    MethodLayer.Write(countlayer,countnstrip,countnstripsat254,countnstripsat255);
                    cout<<LabelLayer(countlayer)+" NStrip="+to_string(countnstrip)+" NStripSat254="+to_string(countnstripsat254)+" NStripSat255="+to_string(countnstripsat255)<<endl;
                }
            }
        }
    }
	MethodLayer.WriteFile();

//METHODE MODULGEOM 

	Correction MethodModulGeom; 
	MethodModulGeom.SetFileAndTreeName("testMethodModulGeom.root","tree");
	MethodModulGeom.SetBranch();
	for(int countlayer=1;countlayer<15;countlayer++)
    {
        for(int countnstrip=3;countnstrip<7;countnstrip++)
        {
            for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
            {
                for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
                {
                    MethodModulGeom.SetHisto(*VectModulGeomVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]);
                    MethodModulGeom.FillHisto(countlayer,countnstrip,countnstripsat254,countnstripsat255);
                    MethodModulGeom.FillProfile();
                    MethodModulGeom.Write(countlayer,countnstrip,countnstripsat254,countnstripsat255);
                }
            }
        }
    }
	MethodModulGeom.WriteFile();


}



//
if(XTalkInvStudy)
{
    for(int j=0;j<h2_ClustSat_MyCorrWOinv_p->GetYaxis()->GetNbins()+1;j++)
    {   
        int Xentries=h2_ClustSat_MyCorrWOinv_p->ProjectionX()->GetEntries();
        for(int i=0;i<h2_ClustSat_MyCorrWOinv_p->GetXaxis()->GetNbins()+1;i++)
        {
            h2_ClustSat_MyCorrWOinv_p->SetBinContent(i,j,(float)h2_ClustSat_MyCorrWOinv_p->GetBinContent(i,j)/Xentries);
        }
    }

    for(int j=0;j<h2_ClustSat_MyCorrWOinv_eta->GetYaxis()->GetNbins()+1;j++)
    {   
        int Xentries=h2_ClustSat_MyCorrWOinv_eta->ProjectionX()->GetEntries();
        for(int i=0;i<h2_ClustSat_MyCorrWOinv_eta->GetXaxis()->GetNbins()+1;i++)
        {
            h2_ClustSat_MyCorrWOinv_eta->SetBinContent(i,j,(float)h2_ClustSat_MyCorrWOinv_eta->GetBinContent(i,j)/Xentries);
        }
    }
}

ofstream ofile_Layer;
ofile_Layer.open("comparison_Layer.txt");
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer1);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer1);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer2);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer2);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer3);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer3);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer4);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer4);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer5);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer5);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer6);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer6);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer7);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer7);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer8);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer8);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer9);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer9);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer10);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer10);
ofile_Layer.close();

ofstream ofile_ModulGeom;
ofile_ModulGeom.open("comparison_ModulGeom.txt");
StudyGraph(ofile_ModulGeom,h_ClustSat_MyCorrWOinv_TIB);
StudyGraph(ofile_ModulGeom,h_ClustSat_HscpCorrWOinv_TIB);
StudyGraph(ofile_ModulGeom,h_ClustSat_MyCorrWOinv_TID);
StudyGraph(ofile_ModulGeom,h_ClustSat_HscpCorrWOinv_TID);
StudyGraph(ofile_ModulGeom,h_ClustSat_MyCorrWOinv_TOB);
StudyGraph(ofile_ModulGeom,h_ClustSat_HscpCorrWOinv_TOB);
StudyGraph(ofile_ModulGeom,h_ClustSat_MyCorrWOinv_TEC);
StudyGraph(ofile_ModulGeom,h_ClustSat_HscpCorrWOinv_TEC);
ofile_ModulGeom.close();

ofstream ofile_Size;
ofile_Size.open("comparison_Size.txt");
StudyGraph(ofile_Size,h_ClustSat_MyCorrWOinv_size1);
StudyGraph(ofile_Size,h_ClustSat_HscpCorrWOinv_size1);
StudyGraph(ofile_Size,h_ClustSat_MyCorrWOinv_size2);
StudyGraph(ofile_Size,h_ClustSat_HscpCorrWOinv_size2);
StudyGraph(ofile_Size,h_ClustSat_MyCorrWOinv_size3);
StudyGraph(ofile_Size,h_ClustSat_HscpCorrWOinv_size3);
StudyGraph(ofile_Size,h_ClustSat_MyCorrWOinv_size4);
StudyGraph(ofile_Size,h_ClustSat_HscpCorrWOinv_size4);
StudyGraph(ofile_Size,h_ClustSat_MyCorrWOinv_size5);
StudyGraph(ofile_Size,h_ClustSat_HscpCorrWOinv_size5);
StudyGraph(ofile_Size,h_ClustSat_MyCorrWOinv_size6);
StudyGraph(ofile_Size,h_ClustSat_HscpCorrWOinv_size6);
ofile_Size.close();

ofstream ofile_NSat;
ofile_NSat.open("comparison_NSat.txt");
StudyGraph(ofile_NSat,h_ClustSat_MyCorrWOinv_sat1);
StudyGraph(ofile_NSat,h_ClustSat_HscpCorrWOinv_sat1);
StudyGraph(ofile_NSat,h_ClustSat_MyCorrWOinv_sat2);
StudyGraph(ofile_NSat,h_ClustSat_HscpCorrWOinv_sat2);
StudyGraph(ofile_NSat,h_ClustSat_MyCorrWOinv_sat3);
StudyGraph(ofile_NSat,h_ClustSat_HscpCorrWOinv_sat3);
ofile_NSat.close();

ofstream ofile_p;
ofile_p.open("comparison_p.txt");
StudyGraph(ofile_p,h_ClustSat_MyCorrWOinv_0p500);
StudyGraph(ofile_p,h_ClustSat_HscpCorrWOinv_0p500);
StudyGraph(ofile_p,h_ClustSat_MyCorrWOinv_500p1000);
StudyGraph(ofile_p,h_ClustSat_HscpCorrWOinv_500p1000);
StudyGraph(ofile_p,h_ClustSat_MyCorrWOinv_1000p1500);
StudyGraph(ofile_p,h_ClustSat_HscpCorrWOinv_1000p1500);
StudyGraph(ofile_p,h_ClustSat_MyCorrWOinv_1500p2000);
StudyGraph(ofile_p,h_ClustSat_HscpCorrWOinv_1500p2000);
StudyGraph(ofile_p,h_ClustSat_MyCorrWOinv_2000p4000);
StudyGraph(ofile_p,h_ClustSat_HscpCorrWOinv_2000p4000);
ofile_p.close();

ofstream ofile_eta;
ofile_eta.open("comparison_eta.txt");
StudyGraph(ofile_eta,h_ClustSat_MyCorrWOinv_0eta05);
StudyGraph(ofile_eta,h_ClustSat_HscpCorrWOinv_0eta05);
StudyGraph(ofile_eta,h_ClustSat_MyCorrWOinv_05eta10);
StudyGraph(ofile_eta,h_ClustSat_HscpCorrWOinv_05eta10);
StudyGraph(ofile_eta,h_ClustSat_MyCorrWOinv_10eta15);
StudyGraph(ofile_eta,h_ClustSat_HscpCorrWOinv_10eta15);
StudyGraph(ofile_eta,h_ClustSat_MyCorrWOinv_15eta25);
StudyGraph(ofile_eta,h_ClustSat_HscpCorrWOinv_15eta25);
ofile_eta.close();

ofstream ofile_firstThreshold;
ofile_firstThreshold.open("comparison_firstThreshold.txt");
StudyGraph(ofile_firstThreshold,h_ClustSat_HscpCorrWinv_firstThreshold20);
StudyGraph(ofile_firstThreshold,h_ClustSat_HscpCorrWinv_firstThreshold50);
StudyGraph(ofile_firstThreshold,h_ClustSat_HscpCorrWinv_firstThreshold100);
StudyGraph(ofile_firstThreshold,h_ClustSat_HscpCorrWinv_firstThreshold150);
StudyGraph(ofile_firstThreshold,h_ClustSat_HscpCorrWinv_firstThreshold200);
ofile_firstThreshold.close();

ofstream ofile_secondThreshold;
ofile_secondThreshold.open("comparison_secondThreshold.txt");
StudyGraph(ofile_secondThreshold,h_ClustSat_HscpCorrWinv_secondThreshold25);
StudyGraph(ofile_secondThreshold,h_ClustSat_HscpCorrWinv_secondThreshold50);
StudyGraph(ofile_secondThreshold,h_ClustSat_HscpCorrWinv_secondThreshold100);
StudyGraph(ofile_secondThreshold,h_ClustSat_HscpCorrWinv_secondThreshold150);
StudyGraph(ofile_secondThreshold,h_ClustSat_HscpCorrWinv_secondThreshold200);
ofile_secondThreshold.close();

ofstream ofile_firstXTalk;
ofile_firstXTalk.open("comparison_firstXTalk.txt");
StudyGraph(ofile_firstXTalk,h_ClustSat_HscpCorrWinv_firstXTalk10);
StudyGraph(ofile_firstXTalk,h_ClustSat_HscpCorrWinv_firstXTalk11);
StudyGraph(ofile_firstXTalk,h_ClustSat_HscpCorrWinv_firstXTalk12);
StudyGraph(ofile_firstXTalk,h_ClustSat_HscpCorrWinv_firstXTalk9);
StudyGraph(ofile_firstXTalk,h_ClustSat_HscpCorrWinv_firstXTalk8);
ofile_firstXTalk.close();

ofstream ofile_secondXTalk;
ofile_secondXTalk.open("comparison_secondXTalk.txt");
StudyGraph(ofile_secondXTalk,h_ClustSat_HscpCorrWinv_secondXTalk4);
StudyGraph(ofile_secondXTalk,h_ClustSat_HscpCorrWinv_secondXTalk5);
StudyGraph(ofile_secondXTalk,h_ClustSat_HscpCorrWinv_secondXTalk6);
StudyGraph(ofile_secondXTalk,h_ClustSat_HscpCorrWinv_secondXTalk3);
StudyGraph(ofile_secondXTalk,h_ClustSat_HscpCorrWinv_secondXTalk2);
ofile_secondXTalk.close();

//ECRITURE RESULTATS ETUDES

TFile* outputfile = new TFile("testEtudeScdeMethod.root","RECREATE");

h2_IhHarmonic_sim->Write();
h2_IhHarmonic_MyCorr->Write();
h2_IhHarmonic_HscpCorrWOinv->Write();
h2_IhHarmonic_HscpCorr->Write();

h_ClustSat_HscpCorrWinv_firstThreshold20->Write();
h_ClustSat_HscpCorrWinv_firstThreshold50->Write();
h_ClustSat_HscpCorrWinv_firstThreshold100->Write();
h_ClustSat_HscpCorrWinv_firstThreshold150->Write();
h_ClustSat_HscpCorrWinv_firstThreshold200->Write();

h_ClustSat_HscpCorrWinv_secondThreshold25->Write();
h_ClustSat_HscpCorrWinv_secondThreshold50->Write();
h_ClustSat_HscpCorrWinv_secondThreshold100->Write();
h_ClustSat_HscpCorrWinv_secondThreshold150->Write();
h_ClustSat_HscpCorrWinv_secondThreshold200->Write();

h_ClustSat_HscpCorrWinv_firstXTalk10->Write();
h_ClustSat_HscpCorrWinv_firstXTalk11->Write();
h_ClustSat_HscpCorrWinv_firstXTalk12->Write();
h_ClustSat_HscpCorrWinv_firstXTalk9->Write();
h_ClustSat_HscpCorrWinv_firstXTalk8->Write();

h_ClustSat_HscpCorrWinv_secondXTalk4->Write();
h_ClustSat_HscpCorrWinv_secondXTalk5->Write();
h_ClustSat_HscpCorrWinv_secondXTalk6->Write();
h_ClustSat_HscpCorrWinv_secondXTalk3->Write();
h_ClustSat_HscpCorrWinv_secondXTalk2->Write();


h_ClustSat_LowP_ClusterSize->Write();
h_ClustSat_LowP_Sat254->Write();
h_ClustSat_LowP_Sat255->Write();
h_ClustSat_LowP_NSat->Write();
h_ClustSat_LowP_layerlabel->Write();
h_ClustSat_LowP_Eta->Write();
h_ClustSat_LowP_Esim->Write();
h_ClustSat_LowP_Emycorr->Write();
h_ClustSat_LowP_EcorrHscp->Write();

h_ClustSat_MyCorrWOinv_TIB->Write();
h_ClustSat_MyCorrWOinv_TID->Write();
h_ClustSat_MyCorrWOinv_TOB->Write();
h_ClustSat_MyCorrWOinv_TEC->Write();

h_ClustSat_HscpCorrWOinv_TIB->Write();
h_ClustSat_HscpCorrWOinv_TID->Write();
h_ClustSat_HscpCorrWOinv_TOB->Write();
h_ClustSat_HscpCorrWOinv_TEC->Write();

h_ClustSat_MyCorrWOinv_Layer1->Write();
h_ClustSat_MyCorrWOinv_Layer2->Write();
h_ClustSat_MyCorrWOinv_Layer3->Write();
h_ClustSat_MyCorrWOinv_Layer4->Write();
h_ClustSat_MyCorrWOinv_Layer5->Write();
h_ClustSat_MyCorrWOinv_Layer6->Write();
h_ClustSat_MyCorrWOinv_Layer7->Write();
h_ClustSat_MyCorrWOinv_Layer8->Write();
h_ClustSat_MyCorrWOinv_Layer9->Write();
h_ClustSat_MyCorrWOinv_Layer10->Write();

h_ClustSat_HscpCorrWOinv_Layer1->Write();
h_ClustSat_HscpCorrWOinv_Layer2->Write();
h_ClustSat_HscpCorrWOinv_Layer3->Write();
h_ClustSat_HscpCorrWOinv_Layer4->Write();
h_ClustSat_HscpCorrWOinv_Layer5->Write();
h_ClustSat_HscpCorrWOinv_Layer6->Write();
h_ClustSat_HscpCorrWOinv_Layer7->Write();
h_ClustSat_HscpCorrWOinv_Layer8->Write();
h_ClustSat_HscpCorrWOinv_Layer9->Write();
h_ClustSat_HscpCorrWOinv_Layer10->Write();

h_ClustSat_MyCorrWOinv_size1->Write();
h_ClustSat_MyCorrWOinv_size2->Write();
h_ClustSat_MyCorrWOinv_size3->Write();
h_ClustSat_MyCorrWOinv_size4->Write();
h_ClustSat_MyCorrWOinv_size5->Write();
h_ClustSat_MyCorrWOinv_size6->Write();
h_ClustSat_MyCorrWOinv_size7->Write();
h_ClustSat_MyCorrWOinv_size8->Write();

h_ClustSat_HscpCorrWOinv_size1->Write();
h_ClustSat_HscpCorrWOinv_size2->Write();
h_ClustSat_HscpCorrWOinv_size3->Write();
h_ClustSat_HscpCorrWOinv_size4->Write();
h_ClustSat_HscpCorrWOinv_size5->Write();
h_ClustSat_HscpCorrWOinv_size6->Write();
h_ClustSat_HscpCorrWOinv_size7->Write();
h_ClustSat_HscpCorrWOinv_size8->Write();

h_ClustSat_MyCorrWOinv_sat1->Write();
h_ClustSat_MyCorrWOinv_sat2->Write();
h_ClustSat_MyCorrWOinv_sat3->Write();

h_ClustSat_HscpCorrWOinv_sat1->Write();
h_ClustSat_HscpCorrWOinv_sat2->Write();
h_ClustSat_HscpCorrWOinv_sat3->Write();


h_ClustSat_MyCorrWOinv_0p500->Write();
h_ClustSat_MyCorrWOinv_500p1000->Write();
h_ClustSat_MyCorrWOinv_1000p1500->Write();
h_ClustSat_MyCorrWOinv_1500p2000->Write();
h_ClustSat_MyCorrWOinv_2000p4000->Write();
h_ClustSat_MyCorrWOinv_0eta05->Write();
h_ClustSat_MyCorrWOinv_05eta10->Write();
h_ClustSat_MyCorrWOinv_10eta15->Write();
h_ClustSat_MyCorrWOinv_15eta25->Write();
h_ClustSat_MyCorrWOinv_20eta25->Write();

h_ClustSat_HscpCorrWOinv_0p500->Write();
h_ClustSat_HscpCorrWOinv_500p1000->Write();
h_ClustSat_HscpCorrWOinv_1000p1500->Write();
h_ClustSat_HscpCorrWOinv_1500p2000->Write();
h_ClustSat_HscpCorrWOinv_2000p4000->Write();
h_ClustSat_HscpCorrWOinv_0eta05->Write();
h_ClustSat_HscpCorrWOinv_05eta10->Write();
h_ClustSat_HscpCorrWOinv_10eta15->Write();
h_ClustSat_HscpCorrWOinv_15eta25->Write();
h_ClustSat_HscpCorrWOinv_20eta25->Write();


TCanvas* cMyCorrDiffP = new TCanvas();
TLegend* legendMyCorrP = new TLegend(0.7,0.7,0.9,0.9);
legendMyCorrP->AddEntry(h_ClustSat_MyCorrWOinv_0p500,"0<p<500","l");
legendMyCorrP->AddEntry(h_ClustSat_MyCorrWOinv_500p1000,"500<p<1000","l");
legendMyCorrP->AddEntry(h_ClustSat_MyCorrWOinv_1000p1500,"1000<p<1500","l");
legendMyCorrP->AddEntry(h_ClustSat_MyCorrWOinv_1500p2000,"1500<p<2000","l");
legendMyCorrP->AddEntry(h_ClustSat_MyCorrWOinv_2000p4000,"2000<p<4000","l");
h_ClustSat_MyCorrWOinv_0p500->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_500p1000->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_1000p1500->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_1500p2000->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_2000p4000->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_0p500->SetLineColor(30);
h_ClustSat_MyCorrWOinv_500p1000->SetLineColor(38);
h_ClustSat_MyCorrWOinv_1000p1500->SetLineColor(42);
h_ClustSat_MyCorrWOinv_1500p2000->SetLineColor(46);
h_ClustSat_MyCorrWOinv_2000p4000->SetLineColor(1);
h_ClustSat_MyCorrWOinv_0p500->Scale(1./h_ClustSat_MyCorrWOinv_0p500->GetEntries());
h_ClustSat_MyCorrWOinv_500p1000->Scale(1./h_ClustSat_MyCorrWOinv_500p1000->GetEntries());
h_ClustSat_MyCorrWOinv_1000p1500->Scale(1./h_ClustSat_MyCorrWOinv_1000p1500->GetEntries());
h_ClustSat_MyCorrWOinv_1500p2000->Scale(1./h_ClustSat_MyCorrWOinv_1500p2000->GetEntries());
h_ClustSat_MyCorrWOinv_2000p4000->Scale(1./h_ClustSat_MyCorrWOinv_2000p4000->GetEntries());
h_ClustSat_MyCorrWOinv_0p500->Draw();
h_ClustSat_MyCorrWOinv_500p1000->Draw("same");
h_ClustSat_MyCorrWOinv_1000p1500->Draw("same");
h_ClustSat_MyCorrWOinv_1500p2000->Draw("same");
h_ClustSat_MyCorrWOinv_2000p4000->Draw("same");
legendMyCorrP->Draw("same");
cMyCorrDiffP->Write();

TCanvas* cMyCorrDiffeta = new TCanvas();
TLegend* legendMyCorreta = new TLegend(0.7,0.7,0.9,0.9);
legendMyCorreta->AddEntry(h_ClustSat_MyCorrWOinv_0eta05,"|#eta|<0.5","l");
legendMyCorreta->AddEntry(h_ClustSat_MyCorrWOinv_05eta10,"0.5<|#eta|<1.0","l");
legendMyCorreta->AddEntry(h_ClustSat_MyCorrWOinv_10eta15,"1.0<|#eta|<1.5","l");
legendMyCorreta->AddEntry(h_ClustSat_MyCorrWOinv_15eta25,"1.5<|#eta|<2.5","l");
//legendMyCorreta->AddEntry(h_ClustSat_MyCorrWOinv_20eta25,"2.0<|#eta|<2.5","l");
h_ClustSat_MyCorrWOinv_0eta05->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_05eta10->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_10eta15->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_15eta25->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_20eta25->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_0eta05->SetLineColor(30);
h_ClustSat_MyCorrWOinv_05eta10->SetLineColor(38);
h_ClustSat_MyCorrWOinv_10eta15->SetLineColor(42);
h_ClustSat_MyCorrWOinv_15eta25->SetLineColor(46);
h_ClustSat_MyCorrWOinv_20eta25->SetLineColor(1);
h_ClustSat_MyCorrWOinv_0eta05->Scale(1./h_ClustSat_MyCorrWOinv_0eta05->GetEntries());
h_ClustSat_MyCorrWOinv_05eta10->Scale(1./h_ClustSat_MyCorrWOinv_05eta10->GetEntries());
h_ClustSat_MyCorrWOinv_10eta15->Scale(1./h_ClustSat_MyCorrWOinv_10eta15->GetEntries());
h_ClustSat_MyCorrWOinv_15eta25->Scale(1./h_ClustSat_MyCorrWOinv_15eta25->GetEntries());
h_ClustSat_MyCorrWOinv_20eta25->Scale(1./h_ClustSat_MyCorrWOinv_20eta25->GetEntries());
h_ClustSat_MyCorrWOinv_0eta05->Draw();
h_ClustSat_MyCorrWOinv_05eta10->Draw("same");
h_ClustSat_MyCorrWOinv_10eta15->Draw("same");
h_ClustSat_MyCorrWOinv_15eta25->Draw("same");
//h_ClustSat_MyCorrWOinv_20eta25->Draw("same");
legendMyCorreta->Draw("same");
cMyCorrDiffeta->Write();

TCanvas* cHscpCorrDiffP = new TCanvas();
TLegend* legendHscpCorrP = new TLegend(0.7,0.7,0.9,0.9);
legendHscpCorrP->AddEntry(h_ClustSat_HscpCorrWOinv_0p500,"0<p<500","l");
legendHscpCorrP->AddEntry(h_ClustSat_HscpCorrWOinv_500p1000,"500<p<1000","l");
legendHscpCorrP->AddEntry(h_ClustSat_HscpCorrWOinv_1000p1500,"1000<p<1500","l");
legendHscpCorrP->AddEntry(h_ClustSat_HscpCorrWOinv_1500p2000,"1500<p<2000","l");
legendHscpCorrP->AddEntry(h_ClustSat_HscpCorrWOinv_2000p4000,"2000<p<4000","l");
h_ClustSat_HscpCorrWOinv_0p500->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_500p1000->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_1000p1500->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_1500p2000->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_2000p4000->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_0p500->SetLineColor(30);
h_ClustSat_HscpCorrWOinv_500p1000->SetLineColor(38);
h_ClustSat_HscpCorrWOinv_1000p1500->SetLineColor(42);
h_ClustSat_HscpCorrWOinv_1500p2000->SetLineColor(46);
h_ClustSat_HscpCorrWOinv_2000p4000->SetLineColor(1);
h_ClustSat_HscpCorrWOinv_0p500->Scale(1./h_ClustSat_HscpCorrWOinv_0p500->GetEntries());
h_ClustSat_HscpCorrWOinv_500p1000->Scale(1./h_ClustSat_HscpCorrWOinv_500p1000->GetEntries());
h_ClustSat_HscpCorrWOinv_1000p1500->Scale(1./h_ClustSat_HscpCorrWOinv_1000p1500->GetEntries());
h_ClustSat_HscpCorrWOinv_1500p2000->Scale(1./h_ClustSat_HscpCorrWOinv_1500p2000->GetEntries());
h_ClustSat_HscpCorrWOinv_2000p4000->Scale(1./h_ClustSat_HscpCorrWOinv_2000p4000->GetEntries());
h_ClustSat_HscpCorrWOinv_0p500->Draw();
h_ClustSat_HscpCorrWOinv_500p1000->Draw("same");
h_ClustSat_HscpCorrWOinv_1000p1500->Draw("same");
h_ClustSat_HscpCorrWOinv_1500p2000->Draw("same");
h_ClustSat_HscpCorrWOinv_2000p4000->Draw("same");
legendHscpCorrP->Draw("same");
cHscpCorrDiffP->Write();

TCanvas* cHscpCorrDiffeta = new TCanvas();
TLegend* legendHscpCorreta = new TLegend(0.7,0.7,0.9,0.9);
legendHscpCorreta->AddEntry(h_ClustSat_HscpCorrWOinv_0eta05,"|#eta|<0.5","l");
legendHscpCorreta->AddEntry(h_ClustSat_HscpCorrWOinv_05eta10,"0.5<|#eta|<1.0","l");
legendHscpCorreta->AddEntry(h_ClustSat_HscpCorrWOinv_10eta15,"1.0<|#eta|<1.5","l");
legendHscpCorreta->AddEntry(h_ClustSat_HscpCorrWOinv_15eta25,"1.5<|#eta|<2.5","l");
//legendHscpCorreta->AddEntry(h_ClustSat_HscpCorrWOinv_20eta25,"2.0<|#eta|<2.5","l");
h_ClustSat_HscpCorrWOinv_0eta05->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_05eta10->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_10eta15->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_15eta25->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_20eta25->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_0eta05->SetLineColor(30);
h_ClustSat_HscpCorrWOinv_05eta10->SetLineColor(38);
h_ClustSat_HscpCorrWOinv_10eta15->SetLineColor(42);
h_ClustSat_HscpCorrWOinv_15eta25->SetLineColor(46);
h_ClustSat_HscpCorrWOinv_20eta25->SetLineColor(1);
h_ClustSat_HscpCorrWOinv_0eta05->Scale(1./h_ClustSat_HscpCorrWOinv_0eta05->GetEntries());
h_ClustSat_HscpCorrWOinv_05eta10->Scale(1./h_ClustSat_HscpCorrWOinv_05eta10->GetEntries());
h_ClustSat_HscpCorrWOinv_10eta15->Scale(1./h_ClustSat_HscpCorrWOinv_10eta15->GetEntries());
h_ClustSat_HscpCorrWOinv_15eta25->Scale(1./h_ClustSat_HscpCorrWOinv_15eta25->GetEntries());
h_ClustSat_HscpCorrWOinv_20eta25->Scale(1./h_ClustSat_HscpCorrWOinv_20eta25->GetEntries());
h_ClustSat_HscpCorrWOinv_0eta05->Draw();
h_ClustSat_HscpCorrWOinv_05eta10->Draw("same");
h_ClustSat_HscpCorrWOinv_10eta15->Draw("same");
h_ClustSat_HscpCorrWOinv_15eta25->Draw("same");
//h_ClustSat_HscpCorrWOinv_20eta25->Draw("same");
legendHscpCorreta->Draw("same");
cHscpCorrDiffeta->Write();

h_ClustSat_DistribP->Write();
h_ClustSat_DistribEta->Write();
h_ClustSat_DistribSize->Write();
h_ClustSat_DistribNSat->Write();

h_DistribP->Write();
h_DistribEta->Write();
h_DistribSize->Write();
h_DistribNSat->Write();

h2_ClustSat_MyCorrWOinv_p->Write();
h2_ClustSat_MyCorrWOinv_eta->Write();

hWithoutCorr->Write();
hCorrHSCP_WithoutInversion->Write();
hCorrHSCP->Write();
hMyCorr_WithoutInversion->Write();
hMyCorr->Write();
hMyCorr_WithoutInversion_RatioSat->Write();
hMyCorr_RatioSat->Write();

h_ClustSat_WithoutCorr->Write();
h_ClustSat_CorrHSCP_WithoutInversion->Write();
h_ClustSat_CorrHSCP->Write();
h_ClustSat_MyCorr_WithoutInversion->Write();
h_ClustSat_MyCorr->Write();
h_ClustSat_MyCorr_WithoutInversion_RatioSat->Write();
h_ClustSat_MyCorr_RatioSat->Write();

h_ClustNoSat_WithoutCorr->Write();
h_ClustNoSat_CorrHSCP_WithoutInversion->Write();
h_ClustNoSat_CorrHSCP->Write();
h_ClustNoSat_MyCorr_WithoutInversion->Write();
h_ClustNoSat_MyCorr->Write();
h_ClustNoSat_MyCorr_WithoutInversion_RatioSat->Write();
h_ClustNoSat_MyCorr_RatioSat->Write();

TCanvas* cStudyCorr = new TCanvas();
TLegend* legendStudyCorr = new TLegend(0.7,0.7,0.9,0.9);
hWithoutCorr->SetLineColor(1);
hWithoutCorr->SetLineWidth(3);
hWithoutCorr->SetLineStyle(3);
legendStudyCorr->AddEntry(hWithoutCorr,"W/O Correction","l");
hCorrHSCP_WithoutInversion->SetLineColor(40);
hCorrHSCP_WithoutInversion->SetLineWidth(3);
legendStudyCorr->AddEntry(hCorrHSCP_WithoutInversion,"HSCP-Corr W/O Inversion","l");
hCorrHSCP->SetLineColor(42);
hCorrHSCP->SetLineWidth(3);
legendStudyCorr->AddEntry(hCorrHSCP,"HSCP-Corr","l");
hMyCorr_WithoutInversion->SetLineColor(38);
hMyCorr_WithoutInversion->SetLineWidth(3);
hMyCorr_WithoutInversion->SetLineStyle(9);
legendStudyCorr->AddEntry(hMyCorr_WithoutInversion,"My correction W/O Inversion","l");
hMyCorr->SetLineColor(46);
hMyCorr->SetLineWidth(3);
hMyCorr->SetLineStyle(9);
legendStudyCorr->AddEntry(hMyCorr,"My correction","l");
hMyCorr_RatioSat->SetLineColor(30);
hMyCorr_RatioSat->SetLineWidth(3);
hMyCorr_RatioSat->SetLineStyle(9);
legendStudyCorr->AddEntry(hMyCorr_RatioSat,"My correction W/ RatioSat","l");
hWithoutCorr->GetYaxis()->SetRangeUser(0,10000);
hWithoutCorr->Draw();
hCorrHSCP_WithoutInversion->Draw("same");
hCorrHSCP->Draw("same");
hMyCorr_WithoutInversion->Draw("same");
hMyCorr->Draw("same");
hMyCorr_RatioSat->Draw("same");
legendStudyCorr->Draw("same");
cStudyCorr->Write();

TCanvas* cmass = new TCanvas();
TLegend* legendmass = new TLegend(0.7,0.7,0.9,0.9);
TLine* line2400 = new TLine(2400,0,2400,100);
legendmass->AddEntry(massEsim,"M(E_{sim})","l");
legendmass->AddEntry(massEfullLayer,"M(E_{full}^{Layer})","l");
massEsim->GetXaxis()->SetTitle("m [GeV/c^{2}]");
massEsim->GetYaxis()->SetTitle("[u.a.]");
massEsim->Draw();
massEfullLayer->Draw("SAME");
legendmass->Draw("same");
line2400->Draw("same");
cmass->Write();

TCanvas* cmassModul = new TCanvas();
TLegend* legendmassModul = new TLegend(0.7,0.7,0.9,0.9);
legendmassModul->AddEntry(massEsim,"M(E_{sim})","l");
legendmassModul->AddEntry(massEfullLayer,"M(E_{full}^{ModulGeom})","l");
massEsim->GetXaxis()->SetTitle("m [GeV/c^{2}]");
massEsim->GetYaxis()->SetTitle("[u.a.]");
massEsim->Draw();
massEfullModulGeom->Draw("SAME");
legendmassModul->Draw("same");
line2400->Draw("same");
cmassModul->Write();


cout<<massEfullLayer->GetMean()<<endl;
cout<<massEfullLayer->GetStdDev()<<endl;
cout<<massEfullModulGeom->GetMean()<<endl;
cout<<massEfullModulGeom->GetStdDev()<<endl;





massEsim->Write();
massErec->Write();
hmassEcorr->Write();
massEfullLayer->Write();
massEfullModulGeom->Write();


pVsdEsimdx->Write();
pVsdErecdx->Write();
pVsdEcorrdx->Write();
pVsdEfulldx_Layer->Write();
pVsdEfulldx_ModulGeom->Write();




distribmassprotonEloss->Write();
distribmassprotonEcorr->Write();
distribmassprotonEcorrNew->Write();


TCanvas* canvasKCecorrnew = new TCanvas();
h2PvsDecorrDxProtonKCnewSatLayer->Draw();
historececorrnew->Draw("same");
canvasKCecorrnew->Write();

TCanvas* canvasKCecorr = new TCanvas();
h2PvsDecorrDxProtonKC->Draw();
historececorr->Draw("same");
canvasKCecorr->Write();

TCanvas* canvasKC = new TCanvas();
h2PvsDelossDxProtonKC->Draw();
historeceloss->Draw("same");
canvasKC->Write();


h2PoverMDecorrDxLayer->Write();
h2PoverMDecorrDxModulGeom->Write();

h2PvsDecorrDxLayer->Write();
h2PvsDecorrDxModulGeom->Write();

hmodulgeom->Write();
hlayerlabel->Write();
SetHistoLabel(hmodulgeomvslayer);
hmodulgeomvslayer->Write();

TCanvas* cLayer = new TCanvas("layer","");
hEcorrLayer->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
hEcorrLayer->Draw();
hErecLayer->SetLineColor(2);
hErecLayer->Draw("same");
cLayer->Write();

TCanvas* cModulGeom = new TCanvas("ModulGeom","");
hEcorrModulGeom->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
hEcorrModulGeom->Draw();
hErecModulGeom->SetLineColor(2);
hErecModulGeom->Draw("same");
cModulGeom->Write();

TCanvas* cModulGeomNoSat = new TCanvas("ModulGeomNoSat","");
hEcorrModulGeomNoSat->Draw();
hErecModulGeomNoSat->SetLineColor(2);
hErecModulGeomNoSat->Draw("same");
cModulGeomNoSat->Write();

TCanvas* cLayerModulGeom = new TCanvas("LayerModulGeom","");
TLegend* legEcorrLayerModulGeom = new TLegend(0.7,0.7,0.9,0.9);
legEcorrLayerModulGeom->AddEntry(hEcorrLayer,"Layer","l");
legEcorrLayerModulGeom->AddEntry(hEcorrModulGeom,"ModulGeom","l");
hEcorrLayer->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
hEcorrLayer->Draw();
hEcorrModulGeom->SetLineColor(2);
hEcorrModulGeom->Draw("same");
legEcorrLayerModulGeom->Draw("same");
cLayerModulGeom->Write();

h2ElossErec->Write();
h2ElossErecLayer->Write();
h2PvsDecorrDxProtonKC->Write();
h2ElossErecModulGeom->Write();
h2ElossEcorrModulGeom->Write();
h2ElossErecModulGeomNoSat->Write();
h2ElossEcorrModulGeomNoSat->Write();


h2EcorrLayerVsEcorrModulGeom->Write();
h1EcorrLayerVsEcorrModulGeom->Write();


h2ElossEcorrSatLayer->Write();







outputfile->Write();
outputfile->Close();

delete outputfile;




return 0;
}
