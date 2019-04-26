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
        int Ndf_;
        float Chi2overNdf_;
        int layerLabel_;
        int nstrip_;
        int nstripsat_;
        bool IsSat255_;
        int nstripsat254_;
        int nstripsat255_;
        TProfile* profile_;
        TH2F* histo_;
        //TFitResultPtr FitRes_;

    public:
        Calibration();
        Calibration(TH2F &histo);
        ~Calibration();
        float CalibCharge(int entry,float charge);
        void FillHisto(int layer,int nstrip,int nstripsat254,int nstripsat255);
        void SetFileAndTreeName(string file_name,string tree_name);
        void SetTree(TTree &tree);
        void SetBranch();
        void Write(int layerLabel,int NStrip,int NStripSat,bool IsSat255);
        void Write(int layerLabel,int nstrip,int nstripsat254,int nstripsat255);
        void WriteFile();
        void Read(TFile *file);
        int GetGoodEntry(int layerLabel,int nstrip,int nstripsat,bool IsSat255);
        void FillProfile();
        void FitProfile();
        void SetHisto(TH2F &histo);
        float RecCharge(float charge,int layerlabel,int nstrip,int nstripsat,bool issat255);
        float ChargeCalib(float charge,int layerLabel,int nstrip,int nstripsat254,int nstripsat255);
};