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
#include "TGraphErrors.h"

using namespace std;

int main(int argc,char** argv)
{
    TProfile* profdedx;
	TFile* _file0 = TFile::Open(argv[1]);
	TH2F* _hdEdx = (TH2F*) _file0->Get("h2PDeDx");
    //TH2F* _hdEdxRec = (TH2F*) _file0->Get("h2PDqcalibDx");
    TH2F* _hdEdxRec = (TH2F*) _file0->Get("h2PDqcorrDx");
    TH2F* _hdedxclone = (TH2F*) _hdEdxRec->Clone();
    _hdedxclone->Reset();
    TH1F* _hId = new TH1F("id","id",100,0,_hdEdxRec->GetNbinsX());
    vector<TH1D*> VectSubHisto;
    TH1D* _hdEdxProjY;
    vector<float> VectFit;
    vector<float> VectFitErr;
    vector<float> vectx;
    int nbre=20;
    TH1F* hrec = new TH1F("","",1000,0,2);
    
    TFile* ofile = new TFile("MassReco.root","RECREATE");

    _hdEdx->Write();
    _hdEdxRec->Write();
    int start=100;
    TFitResultPtr fit_res;

    for(int i=1;i<=_hdEdxRec->GetNbinsX();i++)
    {
        for(int j=1;j<=_hdEdxRec->GetNbinsY();j++)
        {
            _hdedxclone->SetBinContent(i,j,_hdEdxRec->GetBinContent(i,j));
        }
        if(i%10==0)
        {
            _hdEdxProjY = _hdedxclone->ProjectionY();
            _hdEdxProjY->Rebin(2);
            if(_hdEdxProjY->GetEntries()>0)
            {
                fit_res = _hdEdxProjY->Fit("landau","QS");
                hrec->SetBinContent(i,fit_res->Parameter(1));
                hrec->SetBinError(i,fit_res->Error(1));
                hrec->SetMarkerStyle(3);
            }
            _hdedxclone->Write();
            _hdEdxProjY->Write();
            _hdedxclone->Reset();
            
        }
    }
    //TF1* InvSquare = new TF1("InvSquare","[0]/(x*x)+[1]");
    //TF1* InvSquare = new TF1("InvSquare","[0]*(1/x)*(-pow(x,3)/16+pow(x,2)/3-x+1)+[1]/(x*x)");
    //TF1* InvSquare = new TF1("InvSquare","[0]/pow(x,4)+[1]/pow(x,2)+[2]");
    //TF1* InvSquare = new TF1("InvSquare","[0]/pow(x,6)+[1]/pow(x,4)+[2]/pow(x,2)+[3]");
    TF1* InvSquare = new TF1("InvSquare","[0]+[1]*(1/pow(x,2))");
    TFitResultPtr fit = hrec->Fit("InvSquare","S");
    hrec->GetYaxis()->SetRangeUser(0,100);

    TH1F* h_diff = new TH1F("hdiff","hdiff",hrec->GetNbinsX(),0,hrec->GetBinLowEdge(hrec->GetNbinsX()+1));
    TH1F* h_rat = new TH1F("hrat","hrat",hrec->GetNbinsX(),0,hrec->GetBinLowEdge(hrec->GetNbinsX()+1));
    

    for(int i=1;i<_hdEdxRec->GetNbinsX();i++)
    {
        if(i%10==0)
        {
            
            h_diff->SetBinContent(i,hrec->GetBinContent(i)-(fit->Parameter(0)+fit->Parameter(1)*(1/pow(hrec->GetBinLowEdge(i+1),2))));
            h_rat->SetBinContent(i,hrec->GetBinContent(i)/(fit->Parameter(0)+fit->Parameter(1)*(1/pow(hrec->GetBinLowEdge(i+1),2))));
            //cout<<fit_res->Parameter(1)-(fit->Parameter(0)+fit->Parameter(1)*(1/pow(hrec->GetBinLowEdge(i),2)))<<endl;
            cout<<hrec->GetBinLowEdge(i+1)<<endl;
        }
    }

    //TF1* FuncFitDiff = new TF1("FuncFitDiff","")
    //TFItResultPtr = fit_diff = h_diff->Fit("FuncFitDiff","S");


    /*for(int i=0;i<nbre && start<_hdEdx->GetNbinsX();i++)
    {   
        _hdEdxProjY = _hdEdx->ProjectionY("_py",start,start+_hdEdx->GetNbinsX()/nbre);
        _hdEdxProjY->Rebin(10);     
        VectSubHisto.push_back(_hdEdxProjY);
        start+=_hdEdx->GetNbinsX()/nbre;
        if(_hdEdxProjY->GetEntries()>0) fit_res = _hdEdxProjY->Fit("landau","QSR","",_hdEdxProjY->GetBinCenter(_hdEdxProjY->GetMaximumBin()-_hdEdxProjY->GetStdDev()),_hdEdxProjY->GetNbinsX());
        _hdEdxProjY->Write();
        hproj->SetBinContent(i,fit_res->Parameter(0));
        hproj->SetBinError(i,fit_res->Error(0));
        
    }
    TF1* InvExp = new TF1("InvExp","[0]/(x*x)+[1]",0,nbre);
    TFitResultPtr fit = hproj->Fit("InvExp","SR");
    hproj->Write();*/
    hrec->Write();
    h_diff->Write();
    h_rat->Write();

    TCanvas* c1 = new TCanvas();
    _hdEdxRec->Draw();
    hrec->Draw("SAME");
    c1->Write();

    ofile->Write();
    ofile->Close();
    delete ofile;

    float mass=0.938; //en GeV

    cout<<"par0 "<<fit->Parameter(0)<<" par1 "<<fit->Parameter(1)<<" par2 "<<fit->Parameter(2)<<endl;
    cout<<"par0 "<<fit->Parameter(0)<<" par1/M^2 "<<fit->Parameter(1)/pow(mass,2)<<" par2/M^4 "<<fit->Parameter(2)/pow(mass,4)<<endl;

    return 0;
}