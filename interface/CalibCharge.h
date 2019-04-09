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

        string TIB1_254_NStrip3_NStripSat1;           float p0_TIB1_254_NStrip3_NStripSat1, p1_TIB1_254_NStrip3_NStripSat1;
        string TIB1_254_NStrip4_NStripSat1;           float p0_TIB1_254_NStrip4_NStripSat1, p1_TIB1_254_NStrip4_NStripSat1;
        string TIB1_254_NStrip4_NStripSat2;           float p0_TIB1_254_NStrip4_NStripSat2, p1_TIB1_254_NStrip4_NStripSat2;
        string TIB1_254_NStrip5_NStripSat1;           float p0_TIB1_254_NStrip5_NStripSat1, p1_TIB1_254_NStrip5_NStripSat1;
        string TIB1_254_NStrip5_NStripSat2;           float p0_TIB1_254_NStrip5_NStripSat2, p1_TIB1_254_NStrip5_NStripSat2;
        string TIB1_254_NStrip6_NStripSat1;           float p0_TIB1_254_NStrip6_NStripSat1, p1_TIB1_254_NStrip6_NStripSat1;
        string TIB1_254_NStrip6_NStripSat2;           float p0_TIB1_254_NStrip6_NStripSat2, p1_TIB1_254_NStrip6_NStripSat2;

        string TIB1_255_NStrip3_NStripSat1;           float p0_TIB1_255_NStrip3_NStripSat1, p1_TIB1_255_NStrip3_NStripSat1;
        string TIB1_255_NStrip4_NStripSat1;           float p0_TIB1_255_NStrip4_NStripSat1, p1_TIB1_255_NStrip4_NStripSat1;
        string TIB1_255_NStrip4_NStripSat2;           float p0_TIB1_255_NStrip4_NStripSat2, p1_TIB1_255_NStrip4_NStripSat2;
        string TIB1_255_NStrip5_NStripSat1;           float p0_TIB1_255_NStrip5_NStripSat1, p1_TIB1_255_NStrip5_NStripSat1;
        string TIB1_255_NStrip5_NStripSat2;           float p0_TIB1_255_NStrip5_NStripSat2, p1_TIB1_255_NStrip5_NStripSat2;
        string TIB1_255_NStrip6_NStripSat1;           float p0_TIB1_255_NStrip6_NStripSat1, p1_TIB1_255_NStrip6_NStripSat1;
        string TIB1_255_NStrip6_NStripSat2;           float p0_TIB1_255_NStrip6_NStripSat2, p1_TIB1_255_NStrip6_NStripSat2;

        string TIB2_254_NStrip3_NStripSat1;           float p0_TIB2_254_NStrip3_NStripSat1, p1_TIB2_254_NStrip3_NStripSat1;
        string TIB2_254_NStrip4_NStripSat1;           float p0_TIB2_254_NStrip4_NStripSat1, p1_TIB2_254_NStrip4_NStripSat1;
        string TIB2_254_NStrip4_NStripSat2;           float p0_TIB2_254_NStrip4_NStripSat2, p1_TIB2_254_NStrip4_NStripSat2;
        string TIB2_254_NStrip5_NStripSat1;           float p0_TIB2_254_NStrip5_NStripSat1, p1_TIB2_254_NStrip5_NStripSat1;
        string TIB2_254_NStrip5_NStripSat2;           float p0_TIB2_254_NStrip5_NStripSat2, p1_TIB2_254_NStrip5_NStripSat2;
        string TIB2_254_NStrip6_NStripSat1;           float p0_TIB2_254_NStrip6_NStripSat1, p1_TIB2_254_NStrip6_NStripSat1;
        string TIB2_254_NStrip6_NStripSat2;           float p0_TIB2_254_NStrip6_NStripSat2, p1_TIB2_254_NStrip6_NStripSat2;

        string TIB2_255_NStrip3_NStripSat1;           float p0_TIB2_255_NStrip3_NStripSat1, p1_TIB2_255_NStrip3_NStripSat1;
        string TIB2_255_NStrip4_NStripSat1;           float p0_TIB2_255_NStrip4_NStripSat1, p1_TIB2_255_NStrip4_NStripSat1;
        string TIB2_255_NStrip4_NStripSat2;           float p0_TIB2_255_NStrip4_NStripSat2, p1_TIB2_255_NStrip4_NStripSat2;
        string TIB2_255_NStrip5_NStripSat1;           float p0_TIB2_255_NStrip5_NStripSat1, p1_TIB2_255_NStrip5_NStripSat1;
        string TIB2_255_NStrip5_NStripSat2;           float p0_TIB2_255_NStrip5_NStripSat2, p1_TIB2_255_NStrip5_NStripSat2;
        string TIB2_255_NStrip6_NStripSat1;           float p0_TIB2_255_NStrip6_NStripSat1, p1_TIB2_255_NStrip6_NStripSat1;
        string TIB2_255_NStrip6_NStripSat2;           float p0_TIB2_255_NStrip6_NStripSat2, p1_TIB2_255_NStrip6_NStripSat2;

        string TIB3_254_NStrip3_NStripSat1;           float p0_TIB3_254_NStrip3_NStripSat1, p1_TIB3_254_NStrip3_NStripSat1;
        string TIB3_254_NStrip4_NStripSat1;           float p0_TIB3_254_NStrip4_NStripSat1, p1_TIB3_254_NStrip4_NStripSat1;
        string TIB3_254_NStrip4_NStripSat2;           float p0_TIB3_254_NStrip4_NStripSat2, p1_TIB3_254_NStrip4_NStripSat2;
        string TIB3_254_NStrip5_NStripSat1;           float p0_TIB3_254_NStrip5_NStripSat1, p1_TIB3_254_NStrip5_NStripSat1;
        string TIB3_254_NStrip5_NStripSat2;           float p0_TIB3_254_NStrip5_NStripSat2, p1_TIB3_254_NStrip5_NStripSat2;
        string TIB3_254_NStrip6_NStripSat1;           float p0_TIB3_254_NStrip6_NStripSat1, p1_TIB3_254_NStrip6_NStripSat1;
        string TIB3_254_NStrip6_NStripSat2;           float p0_TIB3_254_NStrip6_NStripSat2, p1_TIB3_254_NStrip6_NStripSat2;

        string TIB3_255_NStrip3_NStripSat1;           float p0_TIB3_255_NStrip3_NStripSat1, p1_TIB3_255_NStrip3_NStripSat1;
        string TIB3_255_NStrip4_NStripSat1;           float p0_TIB3_255_NStrip4_NStripSat1, p1_TIB3_255_NStrip4_NStripSat1;
        string TIB3_255_NStrip4_NStripSat2;           float p0_TIB3_255_NStrip4_NStripSat2, p1_TIB3_255_NStrip4_NStripSat2;
        string TIB3_255_NStrip5_NStripSat1;           float p0_TIB3_255_NStrip5_NStripSat1, p1_TIB3_255_NStrip5_NStripSat1;
        string TIB3_255_NStrip5_NStripSat2;           float p0_TIB3_255_NStrip5_NStripSat2, p1_TIB3_255_NStrip5_NStripSat2;
        string TIB3_255_NStrip6_NStripSat1;           float p0_TIB3_255_NStrip6_NStripSat1, p1_TIB3_255_NStrip6_NStripSat1;
        string TIB3_255_NStrip6_NStripSat2;           float p0_TIB3_255_NStrip6_NStripSat2, p1_TIB3_255_NStrip6_NStripSat2;

        string TIB4_254_NStrip3_NStripSat1;           float p0_TIB4_254_NStrip3_NStripSat1, p1_TIB4_254_NStrip3_NStripSat1;
        string TIB4_254_NStrip4_NStripSat1;           float p0_TIB4_254_NStrip4_NStripSat1, p1_TIB4_254_NStrip4_NStripSat1;
        string TIB4_254_NStrip4_NStripSat2;           float p0_TIB4_254_NStrip4_NStripSat2, p1_TIB4_254_NStrip4_NStripSat2;
        string TIB4_254_NStrip5_NStripSat1;           float p0_TIB4_254_NStrip5_NStripSat1, p1_TIB4_254_NStrip5_NStripSat1;
        string TIB4_254_NStrip5_NStripSat2;           float p0_TIB4_254_NStrip5_NStripSat2, p1_TIB4_254_NStrip5_NStripSat2;
        string TIB4_254_NStrip6_NStripSat1;           float p0_TIB4_254_NStrip6_NStripSat1, p1_TIB4_254_NStrip6_NStripSat1;
        string TIB4_254_NStrip6_NStripSat2;           float p0_TIB4_254_NStrip6_NStripSat2, p1_TIB4_254_NStrip6_NStripSat2;

        string TIB4_255_NStrip3_NStripSat1;           float p0_TIB4_255_NStrip3_NStripSat1, p1_TIB4_255_NStrip3_NStripSat1;
        string TIB4_255_NStrip4_NStripSat1;           float p0_TIB4_255_NStrip4_NStripSat1, p1_TIB4_255_NStrip4_NStripSat1;
        string TIB4_255_NStrip4_NStripSat2;           float p0_TIB4_255_NStrip4_NStripSat2, p1_TIB4_255_NStrip4_NStripSat2;
        string TIB4_255_NStrip5_NStripSat1;           float p0_TIB4_255_NStrip5_NStripSat1, p1_TIB4_255_NStrip5_NStripSat1;
        string TIB4_255_NStrip5_NStripSat2;           float p0_TIB4_255_NStrip5_NStripSat2, p1_TIB4_255_NStrip5_NStripSat2;
        string TIB4_255_NStrip6_NStripSat1;           float p0_TIB4_255_NStrip6_NStripSat1, p1_TIB4_255_NStrip6_NStripSat1;
        string TIB4_255_NStrip6_NStripSat2;           float p0_TIB4_255_NStrip6_NStripSat2, p1_TIB4_255_NStrip6_NStripSat2;

        string TOB1_254_NStrip3_NStripSat1;           float p0_TOB1_254_NStrip3_NStripSat1, p1_TOB1_254_NStrip3_NStripSat1;
        string TOB1_254_NStrip4_NStripSat1;           float p0_TOB1_254_NStrip4_NStripSat1, p1_TOB1_254_NStrip4_NStripSat1;
        string TOB1_254_NStrip4_NStripSat2;           float p0_TOB1_254_NStrip4_NStripSat2, p1_TOB1_254_NStrip4_NStripSat2;
        string TOB1_254_NStrip5_NStripSat1;           float p0_TOB1_254_NStrip5_NStripSat1, p1_TOB1_254_NStrip5_NStripSat1;
        string TOB1_254_NStrip5_NStripSat2;           float p0_TOB1_254_NStrip5_NStripSat2, p1_TOB1_254_NStrip5_NStripSat2;
        string TOB1_254_NStrip6_NStripSat1;           float p0_TOB1_254_NStrip6_NStripSat1, p1_TOB1_254_NStrip6_NStripSat1;
        string TOB1_254_NStrip6_NStripSat2;           float p0_TOB1_254_NStrip6_NStripSat2, p1_TOB1_254_NStrip6_NStripSat2;

        string TOB1_255_NStrip3_NStripSat1;           float p0_TOB1_255_NStrip3_NStripSat1, p1_TOB1_255_NStrip3_NStripSat1;
        string TOB1_255_NStrip4_NStripSat1;           float p0_TOB1_255_NStrip4_NStripSat1, p1_TOB1_255_NStrip4_NStripSat1;
        string TOB1_255_NStrip4_NStripSat2;           float p0_TOB1_255_NStrip4_NStripSat2, p1_TOB1_255_NStrip4_NStripSat2;
        string TOB1_255_NStrip5_NStripSat1;           float p0_TOB1_255_NStrip5_NStripSat1, p1_TOB1_255_NStrip5_NStripSat1;
        string TOB1_255_NStrip5_NStripSat2;           float p0_TOB1_255_NStrip5_NStripSat2, p1_TOB1_255_NStrip5_NStripSat2;
        string TOB1_255_NStrip6_NStripSat1;           float p0_TOB1_255_NStrip6_NStripSat1, p1_TOB1_255_NStrip6_NStripSat1;
        string TOB1_255_NStrip6_NStripSat2;           float p0_TOB1_255_NStrip6_NStripSat2, p1_TOB1_255_NStrip6_NStripSat2;

        string TOB2_254_NStrip3_NStripSat1;           float p0_TOB2_254_NStrip3_NStripSat1, p1_TOB2_254_NStrip3_NStripSat1;
        string TOB2_254_NStrip4_NStripSat1;           float p0_TOB2_254_NStrip4_NStripSat1, p1_TOB2_254_NStrip4_NStripSat1;
        string TOB2_254_NStrip4_NStripSat2;           float p0_TOB2_254_NStrip4_NStripSat2, p1_TOB2_254_NStrip4_NStripSat2;
        string TOB2_254_NStrip5_NStripSat1;           float p0_TOB2_254_NStrip5_NStripSat1, p1_TOB2_254_NStrip5_NStripSat1;
        string TOB2_254_NStrip5_NStripSat2;           float p0_TOB2_254_NStrip5_NStripSat2, p1_TOB2_254_NStrip5_NStripSat2;
        string TOB2_254_NStrip6_NStripSat1;           float p0_TOB2_254_NStrip6_NStripSat1, p1_TOB2_254_NStrip6_NStripSat1;
        string TOB2_254_NStrip6_NStripSat2;           float p0_TOB2_254_NStrip6_NStripSat2, p1_TOB2_254_NStrip6_NStripSat2;

        string TOB2_255_NStrip3_NStripSat1;           float p0_TOB2_255_NStrip3_NStripSat1, p1_TOB2_255_NStrip3_NStripSat1;
        string TOB2_255_NStrip4_NStripSat1;           float p0_TOB2_255_NStrip4_NStripSat1, p1_TOB2_255_NStrip4_NStripSat1;
        string TOB2_255_NStrip4_NStripSat2;           float p0_TOB2_255_NStrip4_NStripSat2, p1_TOB2_255_NStrip4_NStripSat2;
        string TOB2_255_NStrip5_NStripSat1;           float p0_TOB2_255_NStrip5_NStripSat1, p1_TOB2_255_NStrip5_NStripSat1;
        string TOB2_255_NStrip5_NStripSat2;           float p0_TOB2_255_NStrip5_NStripSat2, p1_TOB2_255_NStrip5_NStripSat2;
        string TOB2_255_NStrip6_NStripSat1;           float p0_TOB2_255_NStrip6_NStripSat1, p1_TOB2_255_NStrip6_NStripSat1;
        string TOB2_255_NStrip6_NStripSat2;           float p0_TOB2_255_NStrip6_NStripSat2, p1_TOB2_255_NStrip6_NStripSat2;

        string TOB3_254_NStrip3_NStripSat1;           float p0_TOB3_254_NStrip3_NStripSat1, p1_TOB3_254_NStrip3_NStripSat1;
        string TOB3_254_NStrip4_NStripSat1;           float p0_TOB3_254_NStrip4_NStripSat1, p1_TOB3_254_NStrip4_NStripSat1;
        string TOB3_254_NStrip4_NStripSat2;           float p0_TOB3_254_NStrip4_NStripSat2, p1_TOB3_254_NStrip4_NStripSat2;
        string TOB3_254_NStrip5_NStripSat1;           float p0_TOB3_254_NStrip5_NStripSat1, p1_TOB3_254_NStrip5_NStripSat1;
        string TOB3_254_NStrip5_NStripSat2;           float p0_TOB3_254_NStrip5_NStripSat2, p1_TOB3_254_NStrip5_NStripSat2;
        string TOB3_254_NStrip6_NStripSat1;           float p0_TOB3_254_NStrip6_NStripSat1, p1_TOB3_254_NStrip6_NStripSat1;
        string TOB3_254_NStrip6_NStripSat2;           float p0_TOB3_254_NStrip6_NStripSat2, p1_TOB3_254_NStrip6_NStripSat2;

        string TOB3_255_NStrip3_NStripSat1;           float p0_TOB3_255_NStrip3_NStripSat1, p1_TOB3_255_NStrip3_NStripSat1;
        string TOB3_255_NStrip4_NStripSat1;           float p0_TOB3_255_NStrip4_NStripSat1, p1_TOB3_255_NStrip4_NStripSat1;
        string TOB3_255_NStrip4_NStripSat2;           float p0_TOB3_255_NStrip4_NStripSat2, p1_TOB3_255_NStrip4_NStripSat2;
        string TOB3_255_NStrip5_NStripSat1;           float p0_TOB3_255_NStrip5_NStripSat1, p1_TOB3_255_NStrip5_NStripSat1;
        string TOB3_255_NStrip5_NStripSat2;           float p0_TOB3_255_NStrip5_NStripSat2, p1_TOB3_255_NStrip5_NStripSat2;
        string TOB3_255_NStrip6_NStripSat1;           float p0_TOB3_255_NStrip6_NStripSat1, p1_TOB3_255_NStrip6_NStripSat1;
        string TOB3_255_NStrip6_NStripSat2;           float p0_TOB3_255_NStrip6_NStripSat2, p1_TOB3_255_NStrip6_NStripSat2;

        string TOB4_254_NStrip3_NStripSat1;           float p0_TOB4_254_NStrip3_NStripSat1, p1_TOB4_254_NStrip3_NStripSat1;
        string TOB4_254_NStrip4_NStripSat1;           float p0_TOB4_254_NStrip4_NStripSat1, p1_TOB4_254_NStrip4_NStripSat1;
        string TOB4_254_NStrip4_NStripSat2;           float p0_TOB4_254_NStrip4_NStripSat2, p1_TOB4_254_NStrip4_NStripSat2;
        string TOB4_254_NStrip5_NStripSat1;           float p0_TOB4_254_NStrip5_NStripSat1, p1_TOB4_254_NStrip5_NStripSat1;
        string TOB4_254_NStrip5_NStripSat2;           float p0_TOB4_254_NStrip5_NStripSat2, p1_TOB4_254_NStrip5_NStripSat2;
        string TOB4_254_NStrip6_NStripSat1;           float p0_TOB4_254_NStrip6_NStripSat1, p1_TOB4_254_NStrip6_NStripSat1;
        string TOB4_254_NStrip6_NStripSat2;           float p0_TOB4_254_NStrip6_NStripSat2, p1_TOB4_254_NStrip6_NStripSat2;

        string TOB4_255_NStrip3_NStripSat1;           float p0_TOB4_255_NStrip3_NStripSat1, p1_TOB4_255_NStrip3_NStripSat1;
        string TOB4_255_NStrip4_NStripSat1;           float p0_TOB4_255_NStrip4_NStripSat1, p1_TOB4_255_NStrip4_NStripSat1;
        string TOB4_255_NStrip4_NStripSat2;           float p0_TOB4_255_NStrip4_NStripSat2, p1_TOB4_255_NStrip4_NStripSat2;
        string TOB4_255_NStrip5_NStripSat1;           float p0_TOB4_255_NStrip5_NStripSat1, p1_TOB4_255_NStrip5_NStripSat1;
        string TOB4_255_NStrip5_NStripSat2;           float p0_TOB4_255_NStrip5_NStripSat2, p1_TOB4_255_NStrip5_NStripSat2;
        string TOB4_255_NStrip6_NStripSat1;           float p0_TOB4_255_NStrip6_NStripSat1, p1_TOB4_255_NStrip6_NStripSat1;
        string TOB4_255_NStrip6_NStripSat2;           float p0_TOB4_255_NStrip6_NStripSat2, p1_TOB4_255_NStrip6_NStripSat2;

        string TOB5_254_NStrip3_NStripSat1;           float p0_TOB5_254_NStrip3_NStripSat1, p1_TOB5_254_NStrip3_NStripSat1;
        string TOB5_254_NStrip4_NStripSat1;           float p0_TOB5_254_NStrip4_NStripSat1, p1_TOB5_254_NStrip4_NStripSat1;
        string TOB5_254_NStrip4_NStripSat2;           float p0_TOB5_254_NStrip4_NStripSat2, p1_TOB5_254_NStrip4_NStripSat2;
        string TOB5_254_NStrip5_NStripSat1;           float p0_TOB5_254_NStrip5_NStripSat1, p1_TOB5_254_NStrip5_NStripSat1;
        string TOB5_254_NStrip5_NStripSat2;           float p0_TOB5_254_NStrip5_NStripSat2, p1_TOB5_254_NStrip5_NStripSat2;
        string TOB5_254_NStrip6_NStripSat1;           float p0_TOB5_254_NStrip6_NStripSat1, p1_TOB5_254_NStrip6_NStripSat1;
        string TOB5_254_NStrip6_NStripSat2;           float p0_TOB5_254_NStrip6_NStripSat2, p1_TOB5_254_NStrip6_NStripSat2;

        string TOB5_255_NStrip3_NStripSat1;           float p0_TOB5_255_NStrip3_NStripSat1, p1_TOB5_255_NStrip3_NStripSat1;
        string TOB5_255_NStrip4_NStripSat1;           float p0_TOB5_255_NStrip4_NStripSat1, p1_TOB5_255_NStrip4_NStripSat1;
        string TOB5_255_NStrip4_NStripSat2;           float p0_TOB5_255_NStrip4_NStripSat2, p1_TOB5_255_NStrip4_NStripSat2;
        string TOB5_255_NStrip5_NStripSat1;           float p0_TOB5_255_NStrip5_NStripSat1, p1_TOB5_255_NStrip5_NStripSat1;
        string TOB5_255_NStrip5_NStripSat2;           float p0_TOB5_255_NStrip5_NStripSat2, p1_TOB5_255_NStrip5_NStripSat2;
        string TOB5_255_NStrip6_NStripSat1;           float p0_TOB5_255_NStrip6_NStripSat1, p1_TOB5_255_NStrip6_NStripSat1;
        string TOB5_255_NStrip6_NStripSat2;           float p0_TOB5_255_NStrip6_NStripSat2, p1_TOB5_255_NStrip6_NStripSat2;

        string TOB6_254_NStrip3_NStripSat1;           float p0_TOB6_254_NStrip3_NStripSat1, p1_TOB6_254_NStrip3_NStripSat1;
        string TOB6_254_NStrip4_NStripSat1;           float p0_TOB6_254_NStrip4_NStripSat1, p1_TOB6_254_NStrip4_NStripSat1;
        string TOB6_254_NStrip4_NStripSat2;           float p0_TOB6_254_NStrip4_NStripSat2, p1_TOB6_254_NStrip4_NStripSat2;
        string TOB6_254_NStrip5_NStripSat1;           float p0_TOB6_254_NStrip5_NStripSat1, p1_TOB6_254_NStrip5_NStripSat1;
        string TOB6_254_NStrip5_NStripSat2;           float p0_TOB6_254_NStrip5_NStripSat2, p1_TOB6_254_NStrip5_NStripSat2;
        string TOB6_254_NStrip6_NStripSat1;           float p0_TOB6_254_NStrip6_NStripSat1, p1_TOB6_254_NStrip6_NStripSat1;
        string TOB6_254_NStrip6_NStripSat2;           float p0_TOB6_254_NStrip6_NStripSat2, p1_TOB6_254_NStrip6_NStripSat2;

        string TOB6_255_NStrip3_NStripSat1;           float p0_TOB6_255_NStrip3_NStripSat1, p1_TOB6_255_NStrip3_NStripSat1;
        string TOB6_255_NStrip4_NStripSat1;           float p0_TOB6_255_NStrip4_NStripSat1, p1_TOB6_255_NStrip4_NStripSat1;
        string TOB6_255_NStrip4_NStripSat2;           float p0_TOB6_255_NStrip4_NStripSat2, p1_TOB6_255_NStrip4_NStripSat2;
        string TOB6_255_NStrip5_NStripSat1;           float p0_TOB6_255_NStrip5_NStripSat1, p1_TOB6_255_NStrip5_NStripSat1;
        string TOB6_255_NStrip5_NStripSat2;           float p0_TOB6_255_NStrip5_NStripSat2, p1_TOB6_255_NStrip5_NStripSat2;
        string TOB6_255_NStrip6_NStripSat1;           float p0_TOB6_255_NStrip6_NStripSat1, p1_TOB6_255_NStrip6_NStripSat1;
        string TOB6_255_NStrip6_NStripSat2;           float p0_TOB6_255_NStrip6_NStripSat2, p1_TOB6_255_NStrip6_NStripSat2;

        



    public: 
        CalibCharge();
        ~CalibCharge();
        bool SelectedArea(float x1,float y1,float x2,float y2,float x3,float y3,float x4,float y4,float x,float y);
        bool Area(int layerLabel,int NStrip,int NStripSat,bool sat254,float eloss,float charge);
        void SetParameters();
        float Charge(int layerLabel,int NStrip,float NStripSat,bool sat254,float charge);

};



