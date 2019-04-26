#include "../interface/Calibration.h"

using namespace std;

Calibration::Calibration()
{
    p0_     = 0.;
    p1_     = 1.;
    p0err_  = 0.;
    p1err_  = 0.;
    chi2_   = 0.;
    Ndf_    = 0;
    Chi2overNdf_ = 1.;
    layerLabel_ = 0;
    nstrip_ = 0;
    nstripsat_ = 0;
    IsSat255_ = false;
    nstripsat254_   = 0;
    nstripsat255_   = 0;
}

Calibration::Calibration(TH2F &histo)
{
    histo_ = &histo;
    p0_     = 0.;
    p1_     = 1.;
    p0err_  = 0.;
    p1err_  = 0.;
    chi2_   = 0.;
    Ndf_    = 0;
    Chi2overNdf_ = 1.;
    layerLabel_ = 0;
    nstrip_ = 0;
    nstripsat_ = 0;
    IsSat255_ = false;
    nstripsat254_   = 0;
    nstripsat255_   = 0;
}

Calibration::~Calibration()
{

}

void Calibration::SetHisto(TH2F &histo)
{

    histo_ = &histo;
    histo_->Write();
}

void Calibration::FillHisto(int layer,int nstrip,int nstripsat254,int nstripsat255)
{
    float lowedgexu = 1.;
    float lowedgexd = 0.;
    float lowedgeyu = 1.;
    /*if(nstrip>=5)
    {
        if(nstripsat255==0) 
        {
            lowedgexu=0.0015;lowedgexd=0.0005;
            if(nstripsat254==2) {lowedgexu=0.001;lowedgexd=0.0006;}
        }
    }*/
    if(layer<=10 && nstrip==4 && nstripsat254==1 && nstripsat255==0) lowedgeyu=0.0004;
    if(layer<=10 && nstrip==5 && nstripsat254==1 && nstripsat255==0) lowedgeyu=0.00048;
    if(layer>=5 && layer<=10 && nstrip>=3 && nstripsat254==1 && nstripsat255==0) lowedgexu=0.004;
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
            if(LowEdgeX>LowEdgeY /*&& LowEdgeX>lowedgexd*/ && LowEdgeX<lowedgexu && LowEdgeY<lowedgeyu)
            {
                histo_clone->SetBinContent(i,j,histo_->GetBinContent(i,j));
                histo_clone->SetBinError(i,j,histo_->GetBinError(i,j));
            }
        }
    }
    histo_ = histo_clone;
    histo_->Write();
}

void Calibration::FillProfile()
{
    profile_ = histo_->ProfileX("_pfx",1,-1,"");
    profile_->Rebin(5);
    delete histo_;
}

/*void Calibration::FitProfile()
{
    if(profile_->GetEntries()>30) FitRes_ = profile_->Fit("pol1","S");
    else FitRes_=0;
}*/

void Calibration::SetFileAndTreeName(string file_name,string tree_name)
{
    file_ = new TFile(file_name.c_str(),"RECREATE");
    tree_ = new TTree(tree_name.c_str(),tree_name.c_str());
}

void Calibration::SetTree(TTree &tree)
{
    tree_ = &tree;
}

void Calibration::SetBranch()
{
    tree_->Branch("p0",&p0_,"p0/F");
    tree_->Branch("p1",&p1_,"p1/F");
    tree_->Branch("p0err",&p0err_,"p0err/F");
    tree_->Branch("p1err",&p1err_,"p1err/F");
    tree_->Branch("chi2",&chi2_,"chi2/F");
    tree_->Branch("ndf",&Ndf_,"ndf/I");
    tree_->Branch("chi2overndf",&Chi2overNdf_,"chi2overndf/F");
    tree_->Branch("layerLabel",&layerLabel_,"layerLabel/I");
    tree_->Branch("nstrip",&nstrip_,"nstrip/I");
    tree_->Branch("nstripsat",&nstripsat_,"nstripsat/I");
    tree_->Branch("issat255",&IsSat255_,"issat255/O");
//    tree_->Branch("histo","TH2F",&histo_,0,0);
//    tree_->Branch("profile","TProfile",&profile_,0,0);
    tree_->Branch("nstripsat254",&nstripsat254_,"nstripsat254/I");
    tree_->Branch("nstripsat255",&nstripsat255_,"nstripsat255/I");
}

void Calibration::Write(int layerLabel,int NStrip,int NStripSat,bool IsSat255)
{
    if(profile_->GetEntries()>5)
    {    
        TFitResultPtr FitRes_ = profile_->Fit("pol1","SQ");
        p0_     = FitRes_->Parameter(0);
        p1_     = FitRes_->Parameter(1);
        p0err_  = FitRes_->Error(0);
        p1err_  = FitRes_->Error(1);
        chi2_   = FitRes_->Chi2();
        Ndf_    = FitRes_->Ndf();
        Chi2overNdf_ = (float)chi2_/(float)Ndf_;
    //    tree_->Branch("fitresult","TFitResultPtr",&FitRes_,0,0);
    }
    else
    {
        p0_     = 0.;
        p1_     = 1.;
        p0err_  = 0.;
        p1err_  = 0.;
        chi2_   = 0.;
        Ndf_    = 0;
        Chi2overNdf_ = 1.;
    }
    layerLabel_ = layerLabel;
    nstrip_ = NStrip;
    nstripsat_  = NStripSat;
    IsSat255_   = IsSat255;
    tree_->Fill();
}

void Calibration::Write(int layerLabel,int nstrip,int nstripsat254,int nstripsat255)
{
    if(profile_->GetEntries()>10 && (nstripsat254!=0 || nstripsat255!=0))
    {    
        TFitResultPtr FitRes_ = profile_->Fit("pol1","SQ");
        p0_     = FitRes_->Parameter(0);
        p1_     = FitRes_->Parameter(1);
        p0err_  = FitRes_->Error(0);
        p1err_  = FitRes_->Error(1);
        chi2_   = FitRes_->Chi2();
        Ndf_    = FitRes_->Ndf();
        Chi2overNdf_ = (float)chi2_/(float)Ndf_;
    //    tree_->Branch("fitresult","TFitResultPtr",&FitRes_,0,0);
    }
    else
    {
        p0_     = 0.;
        p1_     = 1.;
        p0err_  = 0.;
        p1err_  = 0.;
        chi2_   = 0.;
        Ndf_    = 0;
        Chi2overNdf_ = 0.;
    }
    layerLabel_ = layerLabel;
    nstrip_ = nstrip;
    nstripsat254_ = nstripsat254;
    nstripsat255_ = nstripsat255;
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

float Calibration::CalibCharge(int entry,float charge)
{
    if(entry==0) return charge;
    else
    { 
    tree_->SetBranchAddress("p0",&p0_);
    tree_->SetBranchAddress("p1",&p1_);
    tree_->GetEntry(entry);
    return -(p0_-charge)/p1_;
    }
    
}

float Calibration::RecCharge(float charge,int layerlabel,int nstrip,int nstripsat,bool issat255)
{
    float res=charge;
    tree_->SetBranchAddress("p0",&p0_);
    tree_->SetBranchAddress("p1",&p1_);
    tree_->SetBranchAddress("layerLabel",&layerLabel_);
    tree_->SetBranchAddress("nstrip",&nstrip_);
    tree_->SetBranchAddress("nstripsat",&nstripsat_);
    tree_->SetBranchAddress("issat255",&IsSat255_);
    for(int i=0;i<tree_->GetEntries();i++)
    {
        tree_->GetEntry(i);
        if(nstrip>=6) nstrip=6;
        if(nstripsat>=2) nstripsat=2;
        //if(nstrip==4 && nstripsat==1) {p0_=0.000294361;p1_=0.148644;}
        if(layerlabel==layerLabel_ && nstrip==nstrip_ && nstripsat==nstripsat_ && issat255==IsSat255_) {res=-(p0_-charge)/p1_;}
    }
    return res;
}

float Calibration::ChargeCalib(float charge,int layerLabel,int nstrip,int nstripsat254,int nstripsat255)
{
    float res=charge;
    tree_->SetBranchAddress("p0",&p0_);
    tree_->SetBranchAddress("p1",&p1_);
    tree_->SetBranchAddress("layerLabel",&layerLabel_);
    tree_->SetBranchAddress("nstrip",&nstrip_);
    tree_->SetBranchAddress("nstripsat254",&nstripsat254_);
    tree_->SetBranchAddress("nstripsat255",&nstripsat255_);
    tree_->SetBranchAddress("chi2overndf",&Chi2overNdf_);
    for(int i=0;i<tree_->GetEntries();i++)
    {
        tree_->GetEntry(i);
        if(nstrip>=6) nstrip=6;
        if(nstripsat254>=2) nstripsat254=2;
        if(nstripsat255>=2) nstripsat255=2;
        //if(layerLabel==2) layerLabel=1;
        //if(layerLabel==4) layerLabel=3;
        //if(layerLabel>=6 && layerLabel<=10) layerLabel=5;
        if(layerLabel==layerLabel_ && nstrip==nstrip_ && nstripsat254==nstripsat254_ && nstripsat255==nstripsat255_) {res=-(p0_-charge)/p1_;}
    }
    return res;
}