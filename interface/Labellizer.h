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

#include "../interface/Builder.h"


using namespace std;


const char* StringToChar(string str);
void SetHistoLabel(TCanvas* canvas,TH1F* histo);
void SetHistoLabelPartID(TCanvas* canvas,TH2F* histo);
string Label(int i);
string LabelParticle(int i);
int GetPartID(const vector<Cluster> &VectClust,float &threshold);
int ReBinPartID(int i);
float GetPoverM(float p,int i);
string LoopPartID(int i);
int GetThickness(int layerLabel);
int GetPitch(int layerLabel);
