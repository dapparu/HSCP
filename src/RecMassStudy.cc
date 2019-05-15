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
    TFile* _file0 = TFile::Open(argv[1]);
    TFile* _file1 = TFile::Open(argv[2]);
    //TFile* _file2 = TFile::Open(argv[3]);
    TH1F* histoMassSignal = (TH1F*) _file0->Get("hDistribMassSUSYqcalib");
    TH1F* histoMassBkg = (TH1F*) _file1->Get("hDistribMassSUSYqcalib"); //s'appelle SUSY car SUSY correspond au grand range --> on donne le fichier ttbar comme entree
    TH1F* histoMassSignalBefore = (TH1F*) _file0->Get("hDistribMassSUSYq");
    TH1F* histoMassBkgBefore = (TH1F*) _file1->Get("hDistribMassSUSYq");
    //TH1F* histoMassBkgttbar = (TH1F*) _file2->Get("hDistribMassttbarqcalib");

    TFile* fileout = new TFile("RecMassStudy.root","RECREATE");
    histoMassBkg->SetTitle("Background");
    histoMassSignal->SetTitle("Signal");
    histoMassSignal->Rebin(5);
    //histoMassBkg->Rebin(5);
    histoMassSignalBefore->Rebin(5);
    TF1* puislaw = new TF1("puislaw","[0]*pow(x,[1])",0,400);
    //TFitResultPtr fit_bkg = histoMassBkg->Fit("puislaw","SR","",5,1000);
    //TFitResultPtr fit_bkg = histoMassBkg->Fit("expo","SR0","",100,4000);
    histoMassBkg->Scale(1./histoMassBkg->Integral());
    TFitResultPtr fit_bkg = histoMassBkg->Fit("expo","SR","",20,140);
    TFitResultPtr fit_signal = histoMassSignal->Fit("gaus","SR","",1500,3000);
    float constant_bkg = fit_bkg->Parameter(0);
    float puis_bkg = fit_bkg->Parameter(1);
    float mean_signal = fit_signal->Parameter(1);
    float sigma_signal = fit_signal->Parameter(2)/2;
    TF1* expo_bkg = new TF1("expo_bkg","exp([1]*x+[0])",0,4000);
    expo_bkg->SetParameter(0,constant_bkg);
    expo_bkg->SetParameter(1,puis_bkg);

    histoMassBkgBefore->Scale(1./histoMassBkgBefore->Integral());
    TFitResultPtr fit_bkg_before = histoMassBkgBefore->Fit("expo","SR","",20,140);
    float constant_bkg_before = fit_bkg_before->Parameter(0);
    float puis_bkg_before = fit_bkg_before->Parameter(1);
    TF1* expo_bkg_before = new TF1("expo_bkg_before","exp([1]*x+[0])",0,4000);
    expo_bkg_before->SetParameter(0,constant_bkg_before);
    expo_bkg_before->SetParameter(1,puis_bkg_before);

    //histoMassBkg->GetYaxis()->SetRangeUser(0,histoMassBkg->GetBinContent(histoMassBkg->FindBin(10))*1.05);
    //histoMassBkg->GetXaxis()->SetRangeUser(0,100);
    histoMassSignal->GetYaxis()->SetRangeUser(0,(histoMassSignal->GetBinContent(histoMassSignal->FindBin(mean_signal)))*1.10);
    

    float binning = histoMassSignal->GetBinContent(histoMassSignal->GetNbinsX()+1)/histoMassSignal->GetNbinsX();
    float integral_signal=0.;
    float integral_signal_before=0.;
    float NbinsInterval = histoMassSignal->FindBin(mean_signal+sigma_signal)-histoMassSignal->FindBin(mean_signal-sigma_signal)+1;
    float Nbkg=0.;
    float Nbkg_before=0.;
    for(int i=histoMassSignal->FindBin(mean_signal-sigma_signal);i<histoMassSignal->FindBin(mean_signal+sigma_signal)+1;i++)
    //for(int i=histoMassBkg->FindBin(350);i<histoMassBkg->FindBin(400);i++)
    {
        integral_signal+=histoMassSignal->GetBinContent(i);
        integral_signal_before+=histoMassSignalBefore->GetBinContent(i);
        //Nbkg+=(constant_bkg*pow(histoMassBkg->GetBinLowEdge(i),puis_bkg)/5);
        //Nbkg+=(constant_bkg*TMath::Exp(histoMassBkg->GetBinLowEdge(i)*puis_bkg)*50);
        Nbkg+=(TMath::Exp(histoMassBkg->GetBinLowEdge(i)*puis_bkg+constant_bkg));
        Nbkg_before+=(TMath::Exp(histoMassBkgBefore->GetBinLowEdge(i)*puis_bkg_before+constant_bkg_before));

    }
    cout<<expo_bkg->Integral(300,400)<<endl;
    cout<<expo_bkg_before->Integral(300,400)<<endl;

    float CrossSectionMinBias = 68*pow(10,-3);
    float CrossSectionGluino = 0.128*pow(10,-15);
    float L_data = 140*pow(10,15);
    float factor_signal = L_data*(CrossSectionGluino/15188);
    float factor_bkg = L_data*(CrossSectionMinBias/43163);
    float factor_signal_before = L_data*(CrossSectionGluino/15188);
    float factor_bkg_before = L_data*(CrossSectionMinBias/43163);

    cout<<"interval nbins "<<NbinsInterval<<endl;
    cout<<"integral signal "<<integral_signal<<endl;
    cout<<"integral signal before "<<integral_signal_before<<endl;
    cout<<"nbkg "<<Nbkg<<endl;
    cout<<"nbkg_before "<<Nbkg_before<<endl;

    float significance=(float)integral_signal/sqrt(Nbkg);
    cout<<"significance "<<significance<<endl;

    cout<<"factor_signal "<<factor_signal<<endl;
    cout<<"factor_bkg "<<factor_bkg<<endl;

    float SignalWithFactor = factor_signal*integral_signal;
    float SignalWithFactorBefore = factor_signal_before*integral_signal_before;
    float BkgWithFactor = factor_bkg*Nbkg;
    float BkgWithFactorBefore = factor_bkg_before*Nbkg_before;

    cout<<"factor_signal*signal "<<SignalWithFactor<<endl;
    cout<<"factor_signal*signal_before "<<SignalWithFactorBefore<<endl;
    cout<<"factor_bkg*bkg "<<BkgWithFactor<<endl;
    cout<<"factor_bkg*bkg_before "<<BkgWithFactorBefore<<endl;

    float significanceWithFactors = SignalWithFactor/sqrt(BkgWithFactor);
    float significanceWithFactorsBefore = SignalWithFactorBefore/sqrt(BkgWithFactorBefore);

    cout<<"Significance methode "<<significanceWithFactors<<endl;
    cout<<"Significance avant "<<significanceWithFactorsBefore<<endl;

    cout<<"evolution signifiance : significance/significance_before "<<endl;
    cout<<significanceWithFactors/significanceWithFactorsBefore<<endl;

    histoMassBkg->Write();
    histoMassBkgBefore->Write();
    histoMassSignal->Write();
    histoMassSignalBefore->Write();
    TCanvas* c1 = new TCanvas();
    //expo_bkg->Draw();
    histoMassBkg->Scale(1./histoMassBkg->Integral());
    histoMassBkg->Draw("same");
    c1->Write();
 
    fileout->Write();
    fileout->Close();
    delete fileout;


    return 0;
}