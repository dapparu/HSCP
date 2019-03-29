#include "../interface/PlotterHisto.h"

using namespace std;

void DrawCanvas(TFile* _fileout,TCanvas &canvas)
{
    TFile* fileout = _fileout;
    canvas.Write();
}

void DrawHisto(TFile* _fileout,TH1F* histo,string title,bool x_bool=false,string x_title="")
{
    TFile* fileout = _fileout;
    char title_char[title.length()+1];
    char x_title_char[x_title.length()+1];
    strcpy(title_char,title.c_str());
    strcpy(x_title_char,x_title.c_str());
    histo->SetTitle(title_char);
    if(x_bool) histo->GetXaxis()->SetTitle(x_title_char);
    fileout->Append(histo);
}

void DrawProfile(TFile* _fileout,TProfile* prof,string title,bool x_bool=false,string x_title="")
{
    TFile* fileout = _fileout;
    char title_char[title.length()+1];
    char x_title_char[x_title.length()+1];
    strcpy(title_char,title.c_str());
    strcpy(x_title_char,x_title.c_str());
    prof->SetTitle(title_char);
    if(x_bool) prof->GetXaxis()->SetTitle(x_title_char);
    fileout->Append(prof);
}

void DrawHisto(TFile* _fileout,TH2F* histo,string title,bool x_bool=false,string x_title="",bool y_bool=false,string y_title="")
{
    TFile* fileout = _fileout;
    char title_char[title.length()+1];
    char x_title_char[x_title.length()+1];
    char y_title_char[y_title.length()+1];
    strcpy(title_char,title.c_str());
    strcpy(x_title_char,x_title.c_str());
    strcpy(y_title_char,y_title.c_str());
    histo->SetTitle(title_char);
    if(x_bool) histo->GetXaxis()->SetTitle(x_title_char);
    if(y_bool) histo->GetYaxis()->SetTitle(y_title_char);
    fileout->Append(histo);
}

void SuperposedHisto2DProfile(TCanvas &canvas,TH2F* histo,TProfile* prof,string title,string xtitle, string ytitle)
{
    char title_c[title.length()+1];
    strcpy(title_c,title.c_str());
    char xtitle_c[xtitle.length()+1];
    strcpy(xtitle_c,xtitle.c_str());
    char ytitle_c[ytitle.length()+1];
    strcpy(ytitle_c,ytitle.c_str());
    
    prof->SetMarkerStyle(2);
    prof->SetMarkerColor(kRed);
    prof->SetLineColor(kRed);

    canvas.cd();
    histo->Draw();
    histo->SetDrawOption("COLZ");
    histo->SetTitle(title_c);
    histo->GetXaxis()->SetTitle(xtitle_c);
    histo->GetYaxis()->SetTitle(ytitle_c);
    prof->Draw("SAME");
}

void StackHisto(TCanvas &canvas,vector<TH1F*> VectHisto,vector<char*> VectLegend,string title,string x_title)
{
    char title_char[title.length()+1];
    strcpy(title_char,title.c_str());
    THStack* Stack = new THStack(title_char,title_char);
    title_char[x_title.length()+1];
    strcpy(title_char,x_title.c_str());
    for(int i=0;i<VectHisto.size();i++)
    {
        VectHisto[i]->SetFillColor(36+2*i);
        VectHisto[i]->SetLineColor(36+2*i);
        VectHisto[i]->GetXaxis()->SetTitle(title_char);
        Stack->Add(VectHisto[i]);
    }
    canvas.cd();
    Stack->Draw();
    Stack->GetXaxis()->SetTitle(title_char);
    TLegend* leg_stack = new TLegend(0.7,0.7,0.9,0.9);
    for(int i=0;i<VectLegend.size();i++)
    {
        leg_stack->AddEntry(VectHisto[i],VectLegend[i],"f");
    }
    leg_stack->Draw("SAME");
}