#include <iostream>
#include <fstream>
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

class CalibCharge
{

    private:
        string TOB1NStrip3NStripSat1;           float p0_TOB1NStrip3NStripSat1, p1_TOB1NStrip3NStripSat1;
        string TOB1NStrip4NStripSat1;           float p0_TOB1NStrip4NStripSat1, p1_TOB1NStrip4NStripSat1;
        string TOB1NStrip5NStripSat1;           float p0_TOB1NStrip5NStripSat1, p1_TOB1NStrip5NStripSat1;
        string TOB1NStrip6orMoreNStripSat1;     float p0_TOB1NStrip6orMoreNStripSat1, p1_TOB1NStrip6orMoreNStripSat1;

    public: 
        CalibCharge();
        ~CalibCharge();
        bool SelectedArea(float x1,float y1,float x2,float y2,float x3,float y3,float x4,float y4,float x,float y);
        bool Area(int layerLabel,int NStrip,float eloss,float charge);
        void SetParameters();
        float Charge(int layerLabel,int NStrip,float NStripSat,float charge);

};



