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

using namespace std;


#include "../interface/Builder.h"
#include "../interface/Correction.h"


//FUNCTIONS -- LABEL HISTO

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

	b1->SetThresholdPartId(0.7); 
    b1->SetThresholdPt(0); 
    b1->SetThresholdP(0);
    b1->SetThresholdEta(5);

	int nentries = b1->GetEntries();

	Correction EcorrLayer;
	TFile* fileLayer = TFile::Open("testMethodLayer.root");
	TTree* treeLayer = (TTree*) fileLayer->Get("tree");
	EcorrLayer.SetTree(*treeLayer);
	Correction EcorrModulGeom;
	TFile* fileModulGeom = TFile::Open("testMethodModulGeom.root");
	TTree* treeModulGeom = (TTree*) fileModulGeom->Get("tree");
	EcorrModulGeom.SetTree(*treeModulGeom);
	Correction EcorrModulGeomNoSat;
	TFile* fileModulGeomNoSat = TFile::Open("testMethodModulGeomNoSat.root");
	TTree* treeModulGeomNoSat = (TTree*) fileModulGeomNoSat->Get("tree");
	EcorrModulGeomNoSat.SetTree(*treeModulGeomNoSat);


	bool method=false;
	bool corr=false;



	if(atof(argv[argc-1])==1) method=true;
	if(atof(argv[argc-1])==2) corr=true;


//DECLARATION TABLEAUX POUR LES METHODES 

vector<vector<vector<vector<TH2F*>>>> VectLayerVectNStripVectNStripSat254VectNStripSat255Histo;
vector<vector<vector<vector<TH2F*>>>> VectModulGeomVectNStripVectNStripSat254VectNStripSat255Histo;
vector<vector<vector<vector<TH2F*>>>> VectModulGeomVectNStripVectNStripSat254VectNStripSat255HistoNoSat;
if(method)
{
//METHODE 1 -- CATEGORIES LAYER
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

	for(int ModulGeom=1;ModulGeom<5;ModulGeom++)
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

//METHODE 3 -- CATEGORIES XTALK ET PISTES SANS SATURATION

    for(int ModulGeom=1;ModulGeom<5;ModulGeom++)
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
        VectModulGeomVectNStripVectNStripSat254VectNStripSat255HistoNoSat.push_back(VectNStripVectNStripSat254VectNStripSat255Histo);
    }


}

//DECLARATION HISTOGRAMMES POUR ETUDES

TH1I* hmodulgeom = new TH1I("modulgeom","",25,0,25);
TH1I* hlayerlabel = new TH1I("layerlabel","",25,0,25);

TH2I* hmodulgeomvslayer = new TH2I("hmodulgeomvslayer","",21,0,21,14,0,14);

TH1F* hEcorrLayer = new TH1F("hEcorrLayer","",100,-2,2);
TH1F* hEcorrModulGeom = new TH1F("hEcorrModulGeom","",100,-2,2);
TH1F* hEcorrModulGeomNoSat = new TH1F("hEcorrModulGeomNoSat","",100,-2,2);

TH2F* h2ElossEcorrLayer = new TH2F("h2ElossEcorrLayer","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH2F* h2ElossEcorrModulGeom = new TH2F("h2ElossEcorrModulGeom","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH2F* h2ElossEcorrModulGeomNoSat = new TH2F("h2ElossEcorrModulGeomNoSat","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));

TH2F* h2ElossErec = new TH2F("h2ElossErec","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));

TH2F* h2ElossErecLayer = new TH2F("h2ElossErecLayer","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH2F* h2ElossErecModulGeom = new TH2F("h2ElossErecModulGeom","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH2F* h2ElossErecModulGeomNoSat = new TH2F("h2ElossErecModulGeomNoSat","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));


//BOUCLE SUR LES EVENTS


int entries = 100;

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


		//BOUCLE SUR LES CLUSTERS

        for(int cluster=0;cluster<b1->GetVectTrack()[track].GetNCluster();cluster++)
        {
            float charge        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSclusCharge();
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
			float EcorrSat		= Erec;
			float EcorrNoSat	= ErecNoSat;

			h2ElossErec->Fill(Eloss,Erec);


			hmodulgeom->Fill(modulgeom);
			hlayerlabel->Fill(layerLabel);
			hmodulgeomvslayer->Fill(layerLabel-1,modulgeom-1);


			if(corr)
			{
				if(layerLabel<=2)
				{
					EcorrSat = EcorrLayer.ChargeCorr(Erec,layerLabel,nstrips,nsat254,nsat255);
					hEcorrLayer->Fill((EcorrSat-Erec)/Erec);
					h2ElossEcorrLayer->Fill(Eloss,EcorrSat);
					h2ElossErecLayer->Fill(Eloss,Erec);
				}
				if(modulgeom<=1)
				{
					EcorrSat = EcorrModulGeom.ChargeCorr(Erec,modulgeom,nstrips,nsat254,nsat255);
					hEcorrModulGeom->Fill((EcorrSat-Erec)/Erec);
					h2ElossEcorrModulGeom->Fill(Eloss,EcorrSat);
					h2ElossErecModulGeom->Fill(Eloss,Erec);
					EcorrNoSat = EcorrModulGeomNoSat.ChargeCorr(Erec,modulgeom,nstrips,nsat254,nsat255);
					hEcorrModulGeomNoSat->Fill((EcorrNoSat-Erec)/Erec);
					h2ElossEcorrModulGeomNoSat->Fill(Eloss,EcorrNoSat);
					h2ElossErecModulGeomNoSat->Fill(Eloss,Erec);
				}
			}

			if(method)
			{
				//METHODE LAYER

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
            	                    VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(Eloss,Erec);
            	                } 
            	            }
            	        }
            	    }
				}

				//METHODE MODULGEOM

				for(int countlayer=1;countlayer<5;countlayer++)
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

				//METHODE MODULGEOM SANS PISTES SATUREES

				for(int countlayer=1;countlayer<5;countlayer++)
            	{
            	    for(int countnstrip=3;countnstrip<7;countnstrip++)
            	    {
            	        for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
            	        {
            	            for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
            	            {
            	                if(countlayer==modulgeom && countnstrip==nstrips && countnstripsat254==nsat254 && countnstripsat255==nsat255)
            	                {
            	                    VectModulGeomVectNStripVectNStripSat254VectNStripSat255HistoNoSat[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(Eloss,ErecNoSat);
            	                } 
            	            }
            	        }
            	    }
				}

			}

		}
	}
}







//ECRITURE PARAMETRES METHODE


if(method)
{
//METHODE LAYER 

	Correction MethodLayer; 
	MethodLayer.SetFileAndTreeName("testMethodLayer.root","tree");
	MethodLayer.SetBranch();
	for(int countlayer=1;countlayer<11;countlayer++)
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
                }
            }
        }
    }
	MethodLayer.WriteFile();

//METHODE MODULGEOM 

	Correction MethodModulGeom; 
	MethodModulGeom.SetFileAndTreeName("testMethodModulGeom.root","tree");
	MethodModulGeom.SetBranch();
	for(int countlayer=1;countlayer<5;countlayer++)
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

//METHODE MODULGEOM 

	Correction MethodModulGeomNoSat; 
	MethodModulGeomNoSat.SetFileAndTreeName("testMethodModulGeomNoSat.root","tree");
	MethodModulGeomNoSat.SetBranch();
	for(int countlayer=1;countlayer<5;countlayer++)
    {
        for(int countnstrip=3;countnstrip<7;countnstrip++)
        {
            for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
            {
                for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
                {
                    MethodModulGeomNoSat.SetHisto(*VectModulGeomVectNStripVectNStripSat254VectNStripSat255HistoNoSat[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]);
                    MethodModulGeomNoSat.FillHisto(countlayer,countnstrip,countnstripsat254,countnstripsat255);
                    MethodModulGeomNoSat.FillProfile();
                    MethodModulGeomNoSat.Write(countlayer,countnstrip,countnstripsat254,countnstripsat255);
                }
            }
        }
    }
	MethodModulGeomNoSat.WriteFile();
}










//ECRITURE RESULTATS ETUDES

TFile* outputfile = new TFile("testEtudeScdeMethod.root","RECREATE");

hmodulgeom->Write();
hlayerlabel->Write();
SetHistoLabel(hmodulgeomvslayer);
hmodulgeomvslayer->Write();


hEcorrLayer->Write();
hEcorrModulGeom->Write();
hEcorrModulGeomNoSat->Write();

h2ElossErec->Write();
h2ElossErecLayer->Write();
h2ElossEcorrLayer->Write();
h2ElossErecModulGeom->Write();
h2ElossEcorrModulGeom->Write();
h2ElossErecModulGeomNoSat->Write();
h2ElossEcorrModulGeomNoSat->Write();





outputfile->Write();
outputfile->Close();

delete outputfile;




	return 0;
}
