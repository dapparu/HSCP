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
	vector<TH1F*> VectHisto;
	vector<TFitResultPtr> VectFit;
	vector<float> VectChi2Gaus;
	vector<float> VectConstGaus;
	vector<float> VectConstErrGaus;
	vector<float> VectMeanGaus;
	vector<float> VectMeanErrGaus;
	vector<float> VectSigmaGaus;
	vector<float> VectSigmaErrGaus;
	vector<float> VectRightIntGaus;
	vector<float> VectLeftIntGaus;
	TH1F* hChi2Gaus = new TH1F();
	for(int i=1;i<argc;i++)
	{
		TFile* _file0 = TFile::Open(argv[i]);
		TH1F* _histo0 = (TH1F*) _file0->Get("h1DiffRelEvQCalib");
		TFitResultPtr _fit0 = _histo0->Fit("gaus","QRS","",_histo0->GetBinCenter(_histo0->GetMaximumBin())-_histo0->GetStdDev()/2,_histo0->GetBinCenter(_histo0->GetMaximumBin())+_histo0->GetStdDev()/2);
		VectHisto.push_back(_histo0);
		VectFit.push_back(_fit0);
		VectChi2Gaus.push_back(_fit0->Chi2());
		hChi2Gaus->Fill(_fit0->Chi2());
		VectConstGaus.push_back(_fit0->Parameter(0));
		VectConstErrGaus.push_back(_fit0->ParError(0));
		VectMeanGaus.push_back(_fit0->Parameter(1));
		VectMeanErrGaus.push_back(_fit0->ParError(1));
		VectSigmaGaus.push_back(_fit0->Parameter(2));
		VectSigmaErrGaus.push_back(_fit0->ParError(2));
		VectRightIntGaus.push_back(_histo0->Integral(_histo0->FindBin(_fit0->Parameter(1)+2*_fit0->Parameter(2)),_histo0->GetNbinsX()));
		VectLeftIntGaus.push_back(_histo0->Integral(1,_histo0->FindBin(_fit0->Parameter(1)-2*_fit0->Parameter(2))));
	}
	TH1F* hMean = new TH1F("mean","mean",200,-1,1);
	vector<float> xaxis;
	vector<float> vectmean;
	vector<float> vectstddev;
	for(int i=0;i<VectHisto.size();i++)
	{
		xaxis.push_back((float)i/10);
		hMean->Fill(VectHisto[i]->GetMean());
		vectmean.push_back(VectHisto[i]->GetMean());
		vectstddev.push_back(VectHisto[i]->GetStdDev());
	}
	TFile* fileout = new TFile("study.root","RECREATE");
	for(int i=0;i<VectHisto.size();i++) VectHisto[i]->Write();
	TGraphErrors* grMean = new TGraphErrors(argc-1,&xaxis[0],&vectmean[0]);
	TCanvas* cmean = new TCanvas("mean","mean");
	grMean->Draw("ALP");
	cmean->Write();
	TGraphErrors* grStdDev = new TGraphErrors(argc-1,&xaxis[0],&vectstddev[0]);
	TCanvas* cstddev = new TCanvas("stddev","stddev");
	grStdDev->Draw("ALP");
	cstddev->Write();
	TGraphErrors* grMeanGaus = new TGraphErrors(argc-1,&xaxis[0],&VectMeanGaus[0],0,&VectMeanErrGaus[0]);
	grMeanGaus->SetTitle("Mean");
	grMeanGaus->GetXaxis()->SetTitle("Seuil minimum ratio de saturation");
	TCanvas* cmeangaus = new TCanvas("MeanGaus","MeanGaus");
	grMeanGaus->Draw("ALP");
	cmeangaus->Write();
	TGraphErrors* grSigmaGaus = new TGraphErrors(argc-1,&xaxis[0],&VectSigmaGaus[0],0,&VectSigmaErrGaus[0]);
	grSigmaGaus->SetTitle("Sigma");
	grSigmaGaus->GetXaxis()->SetTitle("Seuil minimum ratio de saturation");
	TCanvas* cSigmagaus = new TCanvas("SigmaGaus","SigmaGaus");
	grSigmaGaus->Draw("ALP");
	cSigmagaus->Write();
	TGraphErrors* grRightIntGaus = new TGraphErrors(argc-1,&xaxis[0],&VectRightIntGaus[0],0,0);
	grRightIntGaus->SetTitle("Integrale droite");
	grRightIntGaus->GetXaxis()->SetTitle("Seuil minimum ratio de saturation");
	TCanvas* cRightIntGaus = new TCanvas("RightIntGaus","RightIntGaus");
	grRightIntGaus->Draw("ALP");
	cRightIntGaus->Write();
	TGraphErrors* grLeftIntGaus = new TGraphErrors(argc-1,&xaxis[0],&VectLeftIntGaus[0],0,0);
	grLeftIntGaus->SetTitle("Integrale gauche");
	grLeftIntGaus->GetXaxis()->SetTitle("Seuil minimum ratio de saturation");
	TCanvas* cLeftIntGaus = new TCanvas("LeftIntGaus","LeftIntGaus");
	grLeftIntGaus->Draw("ALP");
	cLeftIntGaus->Write();
	TGraphErrors* gr1v2 = new TGraphErrors(argc-1,&VectLeftIntGaus[0],&VectRightIntGaus[0],0,0);
	gr1v2->SetTitle("LeftVsRight");
	gr1v2->GetXaxis()->SetTitle("Int left");
	gr1v2->GetXaxis()->SetRangeUser(0,1);
	gr1v2->GetYaxis()->SetTitle("Int right");
	gr1v2->GetYaxis()->SetRangeUser(0,1);
	TCanvas* c1v2 = new TCanvas("1v2","1v2");
	gr1v2->Draw("ALP");
	c1v2->Write();
	fileout->Write();
	fileout->Close();
	delete fileout;


	return 0;
}
