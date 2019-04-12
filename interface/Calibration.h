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

#include "../interface/Builder.h"

using namespace std;

class Calibration
{

    private: 
        TFile* file_;
        TTree* tree_;
        float p0_;
        float p1_;
        float p0err_;
        float p1err_;
        float chi2_;
        int layerLabel_;
        int nstrip_;
        int nstripsat_;
        bool IsSat255_;
        TProfile* profile_;

    public:
        Calibration();
        Calibration(TFile &file);
        ~Calibration();
        float CalibCharge(int entry,float charge);
        void SetFileAndTree(string file_name,string tree_name);
        void SetBranch();
        void Write(float p0,float p1,float chi2);
        void WriteFile();
        void Read(TFile *file);
        int GetGoodEntry(int layerLabel,int nstrip,int nstripsat,bool IsSat255);
        void FillProfile(float E,float Q);

};