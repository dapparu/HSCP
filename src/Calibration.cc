#include "../interface/Calibration.h"

using namespace std;

Calibration::Calibration()
{
    p0_     = 0.;
    p1_     = 0.;
    p0err_  = 0.;
    p1err_  = 0.;
    chi2_   = 0.;
    layerLabel_ = 0;
    nstrip_ = 0;
    nstripsat_ = 0;
    IsSat255_ = false;
}

Calibration::Calibration(TH2F &histo)
{
    histo_ = &histo;
    p0_     = 0.;
    p1_     = 0.;
    p0err_  = 0.;
    p1err_  = 0.;
    chi2_   = 0.;
    layerLabel_ = 0;
    nstrip_ = 0;
    nstripsat_ = 0;
    IsSat255_ = false;
}

Calibration::~Calibration()
{

}

void Calibration::SetHisto(TH2F &histo)
{
    histo_ = &histo;
}

float Calibration::CalibCharge(int entry,float charge)
{
    TTree* tree = (TTree*) file_->Get("tree");   
    tree->SetBranchAddress("p0",&p0_);
    tree->SetBranchAddress("p1",&p1_);
    tree->GetEntry(entry);
    return -(p0_-charge)/p1_;
}

void Calibration::FillHisto(float threshold)
{
    TH2F* histo_clone = (TH2F*) histo_->Clone();
    histo_clone->Reset();
    for(int i=0;i<histo_clone->GetNbinsX()+2;i++)
    {
        float LowEdgeX = histo_->GetXaxis()->GetBinLowEdge(i);
        float WidthX = histo_->GetXaxis()->GetBinWidth(i);
        for(int j=0;j<histo_clone->GetNbinsY()+2;j++)
        {
            float LowEdgeY = histo_->GetXaxis()->GetBinLowEdge(j);
            float WidthY = histo_->GetXaxis()->GetBinWidth(j);
            if(LowEdgeX>LowEdgeY+threshold)
            {
                histo_clone->SetBinContent(i,j,histo_->GetBinContent(i,j));
                histo_clone->SetBinError(i,j,histo_->GetBinError(i,j));
            }
        }
    }
    histo_ = histo_clone;
}

void Calibration::FillProfile()
{
    profile_ = histo_->ProfileX();
}

void Calibration::FitProfile()
{
    FitRes_ = profile_->Fit("pol1","S");
}

void Calibration::SetFileAndTree(string file_name,string tree_name)
{
    file_ = new TFile(file_name.c_str(),"RECREATE");
    tree_ = new TTree(tree_name.c_str(),tree_name.c_str());
}

void Calibration::SetBranch()
{
    tree_->Branch("p0",&p0_,"p0/F");
    tree_->Branch("p1",&p1_,"p1/F");
    tree_->Branch("p0err",&p0err_,"p0err/F");
    tree_->Branch("p1err",&p1err_,"p1err/F");
    tree_->Branch("chi2",&chi2_,"chi2/F");
    tree_->Branch("layerLabel",&layerLabel_,"layerLabel/I");
    tree_->Branch("nstrip",&nstrip_,"nstrip/I");
    tree_->Branch("nstripsat",&nstripsat_,"nstripsat/I");
    tree_->Branch("issat255",&IsSat255_,"issat255/O");
    tree_->Branch("histo","TH2F",&histo_,32000,0);
    tree_->Branch("profile","TProfile",&profile_,32000,0);
}

void Calibration::Write(int layerLabel,int NStrip,int NStripSat,bool IsSat255)
{
    p0_     = FitRes_->Parameter(0);
    p1_     = FitRes_->Parameter(1);
    p0err_  = FitRes_->Error(0);
    p1err_  = FitRes_->Error(1);
    chi2_   = FitRes_->Chi2();
    layerLabel_ = layerLabel;
    nstrip_ = NStrip;
    nstripsat_  = NStripSat;
    IsSat255_   = IsSat255;
    tree_->Fill();
}

void Calibration::WriteFile()
{
    file_->Write();
    file_->Close();
}

void Calibration::Read(TFile *file)
{
    TTree* tree = (TTree*) file->Get("tree");   
    tree_ = tree;
    tree->SetBranchAddress("p0",&p0_);
    tree->SetBranchAddress("p1",&p1_);
    for(int i=0;i<tree->GetEntries();i++)
    {
        tree->GetEntry(i);

    }
    delete tree;
}

int Calibration::GetGoodEntry(int layerLabel,int nstrip,int nstripsat,bool IsSat255)
{
    int indice=0;
    tree_->SetBranchAddress("layerLabel",&layerLabel_);
    tree_->SetBranchAddress("nstrip",&nstrip_);
    tree_->SetBranchAddress("nstripsat",&nstripsat_);
    tree_->SetBranchAddress("issat255",&IsSat255_);
    for(int i=0;i<tree_->GetEntries();i++)
    {
        tree_->GetEntry(i);
        if(layerLabel_==layerLabel && nstrip_==nstrip && nstripsat_==nstripsat && IsSat255_==IsSat255) indice=i;
    }
    return indice;
}
