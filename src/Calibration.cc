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
    profile_ = new TProfile("","",50,0,4500*pow(10,-6),"");
}

Calibration::~Calibration()
{

}

float Calibration::CalibCharge(int entry,float charge)
{
    TTree* tree = (TTree*) file_->Get("tree");   
    tree->SetBranchAddress("p0",&p0_);
    tree->SetBranchAddress("p1",&p1_);
    tree->GetEntry(entry);
    return -(p0_-charge)/p1_;
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
    tree_->Branch("profile","TProfile",&profile_,32000,0);
}

void Calibration::Write(float p0,float p1,float chi2)
{
    p0_     = p0;
    p1_     = p1;
    chi2_   = chi2;
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

void Calibration::FillProfile(float E,float Q)
{
    if(E>Q*(1+0.1)) profile_->Fill(E,Q);
}
