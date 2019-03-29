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

const float m_proton = 938.27; //masse du proton en MeV
const float m_pion = 139.57; //masse du pion en MeV
const float m_Rhadrons = 2400*pow(10,3); //masse des R-hadrons (proche masse gluino) en MeV 

using namespace std;


void SetHistoLabel(TCanvas* canvas,TH1F* histo);
void SetHistoLabelPartID(TCanvas* canvas,TH2F* histo);
string Label(int i);
string LabelParticle(int i);
int GetPartID(const vector<Cluster> &VectClust,float &threshold);
int ReBinPartID(int i);
float GetPoverM(float p,int i);
