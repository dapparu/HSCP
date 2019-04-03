#include "../interface/PlotterHisto.h"

using namespace std;

void DrawHisto(TFile &fileout,TH1F* histo,string title,string x_title)
{
    histo->SetTitle(title.c_str());
    histo->GetXaxis()->SetTitle(x_title.c_str());
    fileout.Append(histo);
}

void DrawProfile(TFile &fileout,TProfile* prof,string title,string x_title)
{
    prof->SetTitle(title.c_str());
    prof->GetXaxis()->SetTitle(x_title.c_str());
    fileout.Append(prof);
}

void DrawHisto(TFile &fileout,TH2F* histo,string title,string x_title,string y_title)
{
    histo->SetTitle(title.c_str());
    histo->GetXaxis()->SetTitle(x_title.c_str());
    histo->GetYaxis()->SetTitle(y_title.c_str());
    fileout.Append(histo);
}

void SuperposedHisto2DProfile(TCanvas &canvas,TH2F* histo,TProfile* prof,string title,string xtitle, string ytitle)
{
    prof->SetMarkerStyle(2);
    prof->SetMarkerColor(kRed);
    prof->SetLineColor(kRed);
    histo->Draw("COLZ");
    histo->SetTitle(title.c_str());
    histo->GetXaxis()->SetTitle(xtitle.c_str());
    histo->GetYaxis()->SetTitle(ytitle.c_str());
    prof->Draw("SAME");
}

void StackHisto(TCanvas &canvas,vector<TH1F*> VectHisto,vector<string> VectLegend,string title,string x_title)
{
    THStack* Stack = new THStack(title.c_str(),title.c_str());
    for(int i=0;i<VectHisto.size();i++)
    {
        VectHisto[i]->SetFillColor(36+2*i);
        VectHisto[i]->SetLineColor(36+2*i);
        VectHisto[i]->GetXaxis()->SetTitle(title.c_str());
        Stack->Add(VectHisto[i]);
    }
    canvas.cd();
    Stack->Draw();
    Stack->GetXaxis()->SetTitle(x_title.c_str());
    TLegend* leg_stack = new TLegend(0.7,0.7,0.9,0.9);
    for(int i=0;i<VectLegend.size();i++)
    {
        leg_stack->AddEntry(VectHisto[i],VectLegend[i].c_str(),"f");
    }
    leg_stack->Draw("SAME");
}

void DrawHistoNormalized(TCanvas &canvas,vector<TH1F*> VectHisto,vector<string> VectLegend,string title,string x_title)
{
    canvas.cd();
    for(int i=0;i<VectHisto.size();i++)
    {
        VectHisto[i]->Scale(1./(VectHisto[i]->Integral()));
        VectHisto[i]->SetLineColor(2*(i+1));
        VectHisto[i]->SetFillColor(0);
        VectHisto[i]->GetXaxis()->SetTitle(x_title.c_str());
        VectHisto[i]->Draw("SAME");
    }
    VectHisto[0]->SetTitle(title.c_str());
    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    for(int i=0;i<VectLegend.size();i++)
    {
        leg->AddEntry(VectHisto[i],VectLegend[i].c_str(),"l");
    }
    leg->Draw("SAME");
}

bool SelectedArea(float x1,float y1,float x2,float y2,float x3,float y3,float x4,float y4,float x,float y)
{
    bool res=false;
    float coef1 = (y3-y1)/(x3-x1);
    float coef2 = (y4-y2)/(x4-x2);
    float coef3 = (y2-y1)/(x2-x1);
    float coef4 = (y4-y3)/(x4-x3);
    float orign1 = y1-coef1*x1;
    float orign2 = y4-coef2*x4;
    float orign3 = y2-coef3*x2;
    float orign4 = y3-coef4*x3;
    //if(x>=(y-orign3)/coef3 && x<=(y-orign4)/coef4 && y>=coef1*x+orign1 && y<=coef2*x+orign2) res=true;
    if(x>=x1 && x>=x2 && x<=x3 && x<=x4 && y>=coef1*x+orign1 && y<=coef2*x+orign2) res=true;
    return res;    
}