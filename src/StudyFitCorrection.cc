#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
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
#include "TGraphErrors.h"
#include "TPaveText.h"

using namespace std;

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

string Label(int modulgeom)
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

int main(int argc,char** argv)
{
    string s1 = argv[1];
    string s2 = s1.substr(0,s1.find('.'))+"_Study.root";
    TFile* ifile = TFile::Open(s1.c_str());
    TTree* tree = (TTree*) ifile->Get("tree");
    TFile* ofile = new TFile(s2.c_str(),"RECREATE");


    float p0    = 0.;
    float p1    = 0.;
    float p0err = 0.;
    float p1err = 0.;
    float chi2  = 0.;
    int ndf     = 0;
    float chi2overndf   = 0.;
    int label      = 0;
    int nstrip          = 0;
    int nstripsat254    = 0;
    int nstripsat255    = 0;

    tree->SetBranchAddress("p0",&p0);
    tree->SetBranchAddress("p1",&p1);
    tree->SetBranchAddress("p0err",&p0err);
    tree->SetBranchAddress("p1err",&p1err);
    tree->SetBranchAddress("chi2",&chi2);
    tree->SetBranchAddress("ndf",&ndf);
    tree->SetBranchAddress("chi2overndf",&chi2overndf);
    tree->SetBranchAddress("label",&label);
    tree->SetBranchAddress("nstrip",&nstrip);
    tree->SetBranchAddress("nstripsat254",&nstripsat254);
    tree->SetBranchAddress("nstripsat255",&nstripsat255);

    TH1F* hChi2overNdf = new TH1F("Chi2/Ndf","Chi2/Ndf",200,0,10);
    TH2F* HistoOfInterest;
    TH1F* ProfileOfInterest;
    vector<float> vect_p0;
    vector<float> vect_p1;
    vector<float> vect_p0err;
    vector<float> vect_p1err;
    vector<float> vect_nstripsat;

    /*float tablayer[10]={1,2,3,4,5,6,7,8,9,10};
    float tabp0[10];//={0,0,0,0,0,0,0,0,0,0};
    float tabp0err[10]={0,0,0,0,0,0,0,0,0,0};*/


    float tabx[10000];
    float tabp0[10000];
    float tabp0err[10000];
    float tabp1[10000];
    float tabp1err[10000];
    for(int i=0;i<10000;i++) tabx[i]=i;

    float tabxTIB1[24];
    float tapp0TIB1[24];
    float tabp0errTIB2[24];
    float tabp1TIB1[24];
    float tabp1errTIB2[24];

    TLine *line = new TLine(0,0,0.0015,0.0015);
    line->SetLineStyle(2);

    float p0_interest=0;
    float p1_interest=0;


int count=0;

    for(int i=0;i<tree->GetEntries();i++)
    {

        tree->GetEntry(i);
        bool test=false;


        if(label==3 && nstrip==5 && nstripsat254==1 && nstripsat255==0) 
        {
            p0_interest=p0;
            p1_interest=p1;
        }
        
        if((chi2overndf<5 && chi2overndf>0) && p1>0.02) // on choisit ici nos criteres d'etude !
        {

            TPaveText *pave = new TPaveText(0.15,0.6,0.3,0.7,"NDC");
            pave->SetFillColor(0);
            pave->AddText(("Chi2 = "+to_string(chi2)).c_str());
            pave->AddText(("Ndf = "+to_string(ndf)).c_str());
            pave->AddText(("Chi2/Ndf = "+to_string(chi2overndf)).c_str());
            pave->AddText(("p1 = "+to_string(p1)).c_str());
            pave->AddText(("p0 = "+to_string(p0)).c_str());
        
            TPaveText *pavetitle = new TPaveText(0.15,0.75,0.6,0.85,"NDC");
            pavetitle->AddText((Label(label)+" Taille du cluster : "+to_string(nstrip)).c_str());
            pavetitle->AddText(("Nombre de strips saturees a 254 : "+to_string(nstripsat254)).c_str());
            pavetitle->AddText(("Nombre de strips saturees a 255 : "+to_string(nstripsat255)).c_str());

            hChi2overNdf->Fill(chi2overndf);
            string str_histo = Label(label)+" NStrip="+to_string(nstrip)+" NStripSat254="+to_string(nstripsat254)+" NStripSat255="+to_string(nstripsat255);
            HistoOfInterest = (TH2F*) ifile->Get(str_histo.c_str());
            cout<<str_histo<<endl;
            HistoOfInterest->Write();
            string str_profile = str_histo+"_NoDiag_rec";
            ProfileOfInterest = (TH1F*) ifile->Get(str_profile.c_str());
            cout<<str_profile<<endl;
            ProfileOfInterest->Write();
            TCanvas* CanvasOfInterest = new TCanvas(str_histo.c_str(),str_histo.c_str());

            /*for(int bin=0;bin<HistoOfInterest->GetXaxis()->GetNbins()+1;bin++)
            {
                for(int bin2=0;bin2<HistoOfInterest->GetYaxis()->GetNbins()+1;bin2++)
                {
                    HistoOfInterest->SetBinContent(bin,bin2,HistoOfInterest->GetBinContent(bin,bin2)/HistoOfInterest->GetEntries());
                }
            }*/
            ProfileOfInterest->SetLineColor(kRed);
            ProfileOfInterest->SetMarkerColor(kRed);
            HistoOfInterest->Draw("colz");
            TF1* myfunc = new TF1("myfunc","x*[0]+[1]",0,5000);
            myfunc->SetParameters(p1,p0);
            myfunc->Draw("same");
            //ProfileOfInterest->Draw("same");
            line->Draw("same");
            pave->Draw("same");
            pavetitle->Draw("same");
            HistoOfInterest->GetXaxis()->SetTitle("Q_{sim} [GeV]");
            HistoOfInterest->GetYaxis()->SetTitle("Q_{rec} [GeV]");
            HistoOfInterest->GetZaxis()->SetTitle("[u.a.]");
            ProfileOfInterest->GetXaxis()->SetTitle("Q_{sim} [GeV]");
            ProfileOfInterest->GetYaxis()->SetTitle("Q_{rec} [GeV]");
            CanvasOfInterest->Write();

            test=true;

        }
        //tabp0[i]=p0;
        //tabp0err[i]=p0err;
        //tabp1[i]=p1;
        //tabp1err[i]=p1err;

/*if(test)
{       
        TCanvas* CanvasMultiple = new TCanvas();
        CanvasMultiple->Divide(4,4);
        for(int canv=1;canv<17;canv++)
        {
            for(int countlayer=1;countlayer<22;countlayer++)
            {
                for(int countnstrip=3;countnstrip<7;countnstrip++)
                {
                    for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
                    {
                        for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
                        {
                            if(countlayer==label && countnstrip==nstrip && countnstripsat254==nstripsat254 && countnstripsat255==nstripsat255)
                            {
                                if(countnstripsat254<3 && countnstripsat255<3)
                                {
                                    CanvasMultiple->cd(canv);
                                    string str_histo2 = Label(label)+" NStrip="+to_string(nstrip)+" NStripSat254="+to_string(nstripsat254)+" NStripSat255="+to_string(nstripsat255);
                                    TH2F* histo = (TH2F*) ifile->Get(str_histo2.c_str());
                                    histo->Draw("colz");
                                    TF1* myfunc2 = new TF1("myfunc2","x*[0]+[1]",0,5000);
                                    myfunc2->SetParameters(p1,p0);
                                    myfunc2->Draw("same");
                                }
                            } 
                        }
                    }
                }
			}
            
        }
        CanvasMultiple->Write();
        delete CanvasMultiple;
}*/

    }
    TH1F* reco = new TH1F("reco","",100,-2,2);
    TH2F* h2_ClosureTest = (TH2F*) ifile->Get("OB1 NStrip=5 NStripSat254=1 NStripSat255=0_NoDiag");
    for(int binx=0;binx<h2_ClosureTest->GetNbinsX();binx++)
    {
        for(int biny=0;biny<h2_ClosureTest->GetNbinsY();biny++)
        {   
            if(h2_ClosureTest->GetBinContent(binx,biny)>0)
            {
                float esim = h2_ClosureTest->GetXaxis()->GetBinCenter(binx);
                float ecorr = (h2_ClosureTest->GetYaxis()->GetBinCenter(biny)-p0_interest)/p1_interest;
                cout<<esim<<"    "<<h2_ClosureTest->GetYaxis()->GetBinCenter(biny)<<"    "<<ecorr<<"    "<<p0_interest<<"   "<<p1_interest<<endl;
                for(int i=1;i<h2_ClosureTest->GetBinContent(binx,biny);i++)reco->Fill((ecorr-esim)/esim);
            }
            
        }
    }
    h2_ClosureTest->Write();
    reco->GetXaxis()->SetRangeUser(-0.2,0.2);
    reco->Write();




    TCanvas* cgr = new TCanvas("graph","graph");
    TGraphErrors* gr1 = new TGraphErrors(10000,tabx,tabp0,0,tabp0err);
    gr1->SetTitle("p0");
    gr1->GetYaxis()->SetTitle("p0");
    gr1->Draw("AP");
    cgr->Write();
    TGraphErrors* gr2 = new TGraphErrors(10000,tabx,tabp1,0,tabp1err);
    gr2->SetTitle("p1");
    gr2->GetYaxis()->SetTitle("p1");
    gr2->Draw("AP");
    cgr->Write();
    ofile->Write();
    ofile->Close();
    delete ofile;


    return 0;
}