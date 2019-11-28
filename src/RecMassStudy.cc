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
    TH1F* histoMassSignal = (TH1F*) _file0->Get("hDistribMassSUSYqcorr");
    TH1F* histoMassBkg = (TH1F*) _file1->Get("hDistribMassSMqcorr"); //s'appelle SUSY car SUSY correspond au grand range --> on donne le fichier ttbar comme entree
    TH1F* histoMassSignalBefore = (TH1F*) _file0->Get("hDistribMassSUSYq");
    TH1F* histoMassBkgBefore = (TH1F*) _file1->Get("hDistribMassSMq");
    //TH1F* histoMassBkgttbar = (TH1F*) _file2->Get("hDistribMassttbarqcorr");
    

    TFile* fileout = new TFile("RecMassStudy.root","RECREATE");
    histoMassBkg->SetTitle("Background");
    histoMassSignal->SetTitle("Signal");
    //histoMassSignal->Rebin(5);
    //histoMassBkg->Rebin(5);
    //histoMassSignalBefore->Rebin(5);
    //TF1* puislaw = new TF1("puislaw","[0]*pow(x,[1])",0,400);
    //TFitResultPtr fit_bkg = histoMassBkg->Fit("puislaw","SR","",5,1000);
    //TFitResultPtr fit_bkg = histoMassBkg->Fit("expo","SR0","",100,4000);
    //histoMassBkg->Scale(1./histoMassBkg->Integral());
    TFitResultPtr fit_bkg = histoMassBkg->Fit("expo","SR","",70,1000);
    TFitResultPtr fit_signal = histoMassSignal->Fit("gaus","SRQ","",1000,3000);
    TFitResultPtr fit_signal_before = histoMassSignalBefore->Fit("gaus","SRQ","",1000,3000);
    float constant_bkg = fit_bkg->Parameter(0);
    float puis_bkg = fit_bkg->Parameter(1);
    float mean_signal = fit_signal->Parameter(1);
    float sigma_signal = fit_signal->Parameter(2)/2;
    TF1* expo_bkg = new TF1("expo_bkg","exp([1]*x+[0])",0,4000);
    expo_bkg->SetParameter(0,constant_bkg);
    expo_bkg->SetParameter(1,puis_bkg);


    int integ_mass = 0;
    for(int i=histoMassSignal->FindBin(550);i<histoMassSignal->GetNbinsX();i++)
    {
        //cout<<histoMassSignal->GetBinLowEdge(i)<<endl;
        integ_mass+=histoMassSignal->GetBinContent(i);
    } 

    int integ_mass_bkg = 0;
    for(int i=400;i<3000;i++)
    {
        integ_mass_bkg+=exp(puis_bkg*i+constant_bkg);
    } 

    float ecart_signal = 2400-fit_signal->Parameter(1);
    float ecart_signal_before = 2400-fit_signal_before->Parameter(1);

    int ecart_signal_bin = histoMassSignal->FindBin(2400)-histoMassSignal->FindBin(fit_signal->Parameter(1));
    int ecart_signal_before_bin = histoMassSignalBefore->FindBin(2400)-histoMassSignalBefore->FindBin(fit_signal_before->Parameter(1));

cout<<"signal "<<ecart_signal<<endl;
cout<<ecart_signal_bin<<endl;
cout<<"signal before "<<ecart_signal_before<<endl;
cout<<ecart_signal_before_bin<<endl;


    //histoMassBkgBefore->Scale(1./histoMassBkgBefore->Integral());
    TFitResultPtr fit_bkg_before = histoMassBkgBefore->Fit("expo","SRQ","",70,1000);
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
    double Nbkg=0.;
    double Nbkg_before=0.;
    float IntOfInterest=0.68;
    float intSignal=0.;
    int BorneGaucheSignal=0;
    int BorneDroiteSignal=0;
    bool test1=true;
    bool test2=true;
    for(int i=1;i<histoMassSignal->GetNbinsX()+1;i++)
    {
        intSignal+=histoMassSignal->GetBinContent(i)/histoMassSignal->Integral();
        if(intSignal>=(1-IntOfInterest)/2 && test1) {BorneGaucheSignal=i;test1=false;}
        if(intSignal>=1-(1-IntOfInterest)/2 && test2) {BorneDroiteSignal=i;test2=false;}
    }

    cout<<"intsignal "<<intSignal<<endl;
    cout<<BorneGaucheSignal<<" "<<BorneDroiteSignal<<endl;
    cout<<histoMassSignal->Integral(BorneGaucheSignal,BorneDroiteSignal)/histoMassSignal->Integral()<<endl;


    //for(int i=histoMassSignal->FindBin(mean_signal-1*sigma_signal);i<histoMassSignal->FindBin(mean_signal+1*sigma_signal)+1;i++) cout<<i<<endl;
    //for(int i=histoMassBkg->FindBin(350);i<histoMassBkg->FindBin(400);i++)
    //for(int i=BorneGaucheSignal;i<BorneDroiteSignal+1;i++)
    for(int i=BorneGaucheSignal;i<BorneDroiteSignal+1;i++)
    {
        integral_signal+=histoMassSignal->GetBinContent(i);
        //integral_signal_before+=histoMassSignalBefore->GetBinContent(i);
        //Nbkg+=(constant_bkg*pow(histoMassBkg->GetBinLowEdge(i),puis_bkg)/5);
        //Nbkg+=(constant_bkg*TMath::Exp(histoMassBkg->GetBinLowEdge(i)*puis_bkg)*50);
        //cout<<histoMassSignal->GetBinLowEdge(i)<<endl;
        Nbkg+=(TMath::Exp(histoMassSignal->GetBinLowEdge(i)*puis_bkg+constant_bkg));
        //Nbkg_before+=(TMath::Exp(histoMassBkgBefore->GetBinLowEdge(i)*puis_bkg_before+constant_bkg_before));

    }
    cout<<"expo integral "<<expo_bkg->Integral(histoMassSignal->GetBinLowEdge(BorneGaucheSignal),histoMassSignal->GetBinLowEdge(BorneDroiteSignal))<<endl;
    cout<<"expo integral2 "<<expo_bkg->Integral(200,400)/10<<endl;
    cout<<"expo integral3 "<<expo_bkg->Integral(100,600)/10<<endl;

    /*Nbkg=0.;
    for(int i=BorneGaucheSignal+ecart_signal_bin;i<BorneDroiteSignal+ecart_signal_bin+1;i++) 
    {
        Nbkg+=(TMath::Exp(histoMassSignal->GetBinLowEdge(i)*puis_bkg+constant_bkg));
    }*/
cout<<"1     gauchebefore      "<<histoMassSignal->GetBinLowEdge(BorneGaucheSignal)<<"     droite     "<<histoMassSignal->GetBinLowEdge(BorneDroiteSignal)<<endl;
    Nbkg=expo_bkg->Integral(histoMassSignal->GetBinLowEdge(BorneGaucheSignal),histoMassSignal->GetBinLowEdge(BorneDroiteSignal))/10;

    cout<<"Nbkg    "<<Nbkg<<endl;

    intSignal=0.;
    BorneGaucheSignal=0;
    BorneDroiteSignal=0;
    test1=true;
    test2=true;
    for(int i=1;i<histoMassSignalBefore->GetNbinsX()+1;i++)
    {
        intSignal+=histoMassSignalBefore->GetBinContent(i)/histoMassSignalBefore->Integral();
        if(intSignal>=(1-IntOfInterest)/2 && test1) {BorneGaucheSignal=i;test1=false;}
        if(intSignal>=1-(1-IntOfInterest)/2 && test2) {BorneDroiteSignal=i;test2=false;}
    }

    //cout<<"intsignal "<<intSignal<<endl;
    cout<<BorneGaucheSignal<<" "<<BorneDroiteSignal<<endl;
    cout<<histoMassSignalBefore->Integral(BorneGaucheSignal,BorneDroiteSignal)/histoMassSignalBefore->Integral()<<endl;
   
    for(int i=BorneGaucheSignal;i<BorneDroiteSignal+1;i++)
    {
        integral_signal_before+=histoMassSignalBefore->GetBinContent(i);
        //integral_signal_before+=histoMassSignalBefore->GetBinContent(i);
        //Nbkg+=(constant_bkg*pow(histoMassBkg->GetBinLowEdge(i),puis_bkg)/5);
        //Nbkg+=(constant_bkg*TMath::Exp(histoMassBkg->GetBinLowEdge(i)*puis_bkg)*50);
        //Nbkg_before+=(TMath::Exp(histoMassBkg->GetBinLowEdge(i)*puis_bkg+constant_bkg));
        //cout<<histoMassSignalBefore->GetBinLowEdge(i)<<endl;
        Nbkg_before+=(TMath::Exp(histoMassSignalBefore->GetBinLowEdge(i)*puis_bkg_before+constant_bkg_before));

    }

    //Nbkg_before=0.;

    //cout<<"              here"<<endl;
    
    /*for(int i=BorneGaucheSignal+ecart_signal_before_bin;i<BorneDroiteSignal+ecart_signal_before_bin+1;i++) 
    {
        cout<<histoMassSignalBefore->GetBinLowEdge(i)<<endl;
        Nbkg_before+=(TMath::Exp(histoMassSignalBefore->GetBinLowEdge(i)*puis_bkg_before+constant_bkg_before));
    }*/
cout<<"gauchebefore      "<<histoMassSignal->GetBinLowEdge(BorneGaucheSignal)<<"     droite     "<<histoMassSignal->GetBinLowEdge(BorneDroiteSignal)<<endl;
    Nbkg_before=expo_bkg_before->Integral(histoMassSignal->GetBinLowEdge(BorneGaucheSignal),histoMassSignal->GetBinLowEdge(BorneDroiteSignal))/10;

    cout<<"Nbkg_before    "<<Nbkg_before<<endl;

    float CrossSectionMinBias = 68*pow(10,-3);
    float CrossSectionGluino = 0.128*pow(10,-15);
    float L_data = 140*pow(10,15);
    float factor_signal = L_data*(CrossSectionGluino/15188);
    //float factor_bkg = L_data*(CrossSectionMinBias/69763);
    float factor_bkg = L_data*(CrossSectionMinBias/261000);
    float factor_signal_before = L_data*(CrossSectionGluino/15188);
    //float factor_bkg_before = L_data*(CrossSectionMinBias/69763);
    float factor_bkg_before = L_data*(CrossSectionMinBias/261000);

    cout<<"integral signal "<<integral_signal<<endl;
    cout<<"integral signal before "<<integral_signal_before<<endl;
    cout<<"nbkg "<<Nbkg<<endl;
    cout<<"nbkg_before "<<Nbkg_before<<endl;

    float significance=(float)integral_signal/sqrt(Nbkg);
    cout<<"significance "<<significance<<endl;

    cout<<"factor_signal "<<factor_signal<<endl;
    cout<<"factor_bkg "<<factor_bkg<<endl;
    cout<<"factor_signal_before "<<factor_signal_before<<endl;
    cout<<"factor_bkg_before "<<factor_bkg_before<<endl;

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

    float signi_papier = sqrt(2*((integral_signal+Nbkg)*log(1+integral_signal/Nbkg)-integral_signal));
    //float signi_papier = (2*((SignalWithFactor+BkgWithFactor)*log(1+SignalWithFactor/BkgWithFactor)-SignalWithFactor));
    //cout<<"signi_papier "<<signi_papier<<endl;

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