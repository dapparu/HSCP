#include "../interface/PlotterHisto.h"

using namespace std;


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