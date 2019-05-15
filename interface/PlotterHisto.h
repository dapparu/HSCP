#include <iostream>
#include <vector>
#include <stdlib.h>
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

#include "Labellizer.h"

using namespace std;

void DrawHisto(TFile &fileout,TH1F* histo,string title,string x_title="");
void DrawProfile(TFile &fileout,TProfile* prof,string title,string x_title="");
void DrawHisto(TFile &fileout,TH2F* histo,string title,string x_title="",string y_title="");
void DrawHisto(TFile &fileout,vector<TH2F*> VectHisto,vector<string> VectLegend,string title,string x_title="",string y_title="");
void SuperposedHisto2DProfile(TCanvas &canvas,TH2F* histo,TProfile* prof,string title="",string x_title="", string y_title="");
void StackHisto(TCanvas &canvas,vector<TH1F*> VectHisto,vector<string> VectLegend,string title,string x_title);
void DrawHistoNormalized(TCanvas &canvas,vector<TH1F*> VectHisto,vector<string> VectLegend,string title,string x_title);
void DrawHisto(TCanvas &canvas,vector<TH1F*> VectHisto,vector<string> VectLegend,string title,string x_title);
void DrawHisto(TFile &fileout,vector<TH1F*> VectHisto,vector<string> VectLegend,string title,string x_title);
void DrawClusterProfile(const Cluster* clus);
bool SelectedArea(float x1,float y1,float x2,float y2,float x3,float y3,float x4,float y4,float x,float y);