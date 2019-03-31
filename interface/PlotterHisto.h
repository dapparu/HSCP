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

using namespace std;

void DrawCanvas(TFile* _fileout,TCanvas &canvas);
void DrawHisto(TFile* _fileout,TH1F* histo,string title,bool x_bool,string x_title);
void DrawProfile(TFile* _fileout,TProfile* prof,string title,bool x_bool,string x_title);
void DrawHisto(TFile* _fileout,TH2F* histo,string title,bool x_bool,string x_title,bool y_bool,string y_title);
void SuperposedHisto2DProfile(TCanvas &canvas,TH2F* histo,TProfile* prof,string title="",string x_title="", string y_title="");
void StackHisto(TCanvas &canvas,vector<TH1F*> VectHisto,vector<char*> VectLegend,string title,string x_title);
void DrawHistoNormalized(TCanvas &canvas,vector<TH1F*> VectHisto,vector<char*> VectLegend,string title,string x_title);