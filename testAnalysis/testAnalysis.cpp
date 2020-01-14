#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <iterator>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH1D.h"
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
#include "TExec.h"
#include "TPaveText.h"
#include "TMatrix.h"

using namespace std;


#include "../interface/Builder.h"
#include "../interface/Correction.h"
#include "../interface/Estimator.h"


#include "../SaturationCorrection.h"



//CONSTANTES K ET C

const float CsusyCorr = 1.14; //pour gluino
const float KsusyCorr = 6.7;



const float CprotonEloss=3.27;
const float KprotonEloss=2.49;

const float CprotonEcorrNew=6.73;
const float KprotonEcorrNew=1.18;

const float CprotonEcorr=6.15;
const float KprotonEcorr=2.91; 

const float MGluino2400 = 2400.;
const float MGluino1400 = 1400.;
const float MStau651 = 651.;







//FUNCTIONS -- LABEL HISTO

void StudyGraph(ofstream &ofile, TH1F* histo)
{
    float DistribMean = histo->GetMean();
    float DistribStdDev = histo->GetStdDev();
    float DistribMeanErr = histo->GetMeanError();
    float DistribStdDevErr = histo->GetStdDevError();
    int FitStatus = histo->Fit("gaus","QRN","",-0.1,0.1);
    //TFitResultPtr fit = histo->Fit("gaus","QRS","",DistribMean-DistribStdDev,DistribMean+DistribStdDev);
    TFitResultPtr fit = histo->Fit("gaus","QRSN","",-0.1,0.1);
    float FitMean = 0;
    float FitSigma = 0;
    float FitMeanErr = 0;
    float FitSigmaErr = 0;
    if(FitStatus>-1) 
    {
        FitMean = fit->Parameter(1);
        FitSigma = fit->Parameter(2);
        FitMeanErr = fit->Error(1);
        FitSigmaErr = fit->Error(2);
    }
    int CountPass=0;
    int CountPassGlissantSigma=0;
    int CountPassGlissant=0;
    int CountPassRight=0;
    int CountPassLeft=0;
    int bornegche=0;
    int bornedrte=0;
    for(int i=0;i<histo->GetNbinsX();i++)
    {
        if(histo->GetBinCenter(i)>=FitMean-FitSigma && histo->GetBinCenter(i)<=FitMean+FitSigma) 
        {
            CountPassGlissantSigma+=histo->GetBinContent(i);
        }
        if(histo->GetBinCenter(i)>=FitMean-0.34 && histo->GetBinCenter(i)<=FitMean+0.34) 
        {
            CountPassGlissant+=histo->GetBinContent(i);
        }
        if(histo->GetBinCenter(i)>=-0.25 && histo->GetBinCenter(i)<=0.25)
        {
            CountPass+=histo->GetBinContent(i);
        }
        if(histo->GetBinCenter(i)>=0.4) bornedrte=i;
        if(histo->GetBinCenter(i)<=-0.4) bornegche=i;
    }
    float FracInterval = (float)CountPass/(float)histo->GetEntries();
    float FracIntervalGlissantSigma = (float)CountPassGlissantSigma/(float)histo->GetEntries();
    float FracIntervalGlissant = (float)CountPassGlissant/(float)histo->GetEntries();
    //float FracRight = (float)CountPassRight/(float)histo->GetEntries();
    //float FracLeft = (float)CountPassLeft/(float)histo->GetEntries();
    float FracRight = (float)histo->Integral(bornedrte,histo->GetNbinsX()+1)/(float)histo->Integral();
    float FracLeft = (float)histo->Integral(0,bornegche)/(float)histo->Integral();
    //ofile<<"Name"<<"\t\t"<<"DistribMean"<<"\t"<<"DistribMeanErr"<<"\t"<<"DistribStdDev"<<"\t"<<"DistribStdDevErr"<<"\t"<<"FitMean"<<"\t"<<"FitMeanErr"<<"\t"<<"FitSigma"<<"\t"<<"FitSigmaErr"<<"\t"<<"FracIntervalGlissantSigma"<<"\t"<<"FracInterval"<<endl;
    ofile<<endl;
    ofile<<histo->GetName()<<"\t\t"<<histo->GetMaximum()<<"\t"<<DistribMean<<"\t"<<DistribMeanErr<<"\t"<<DistribStdDev<<"\t"<<DistribStdDevErr<<"\t FitMean : "<<FitMean<<"\t"<<FitMeanErr<<"\t FitSigma : "<<FitSigma<<"\t"<<FitSigmaErr<<"\t"<<FracIntervalGlissantSigma<<"\t"<<FracIntervalGlissant<<"\t"<<FracInterval<<"\t FracRight : "<<FracRight<<"\t FracLeft : "<<FracLeft<<endl;
    
}

void SetHistoLabel(TH2I* histo)
{
	histo->GetXaxis()->SetBinLabel(1,"TIB 1");
    histo->GetXaxis()->SetBinLabel(2,"TIB 2");
    histo->GetXaxis()->SetBinLabel(3,"TIB 3");
    histo->GetXaxis()->SetBinLabel(4,"TIB 4");
    histo->GetXaxis()->SetBinLabel(5,"TOB 1");
    histo->GetXaxis()->SetBinLabel(6,"TOB 2");
    histo->GetXaxis()->SetBinLabel(7,"TOB 3");
    histo->GetXaxis()->SetBinLabel(8,"TOB 4");
    histo->GetXaxis()->SetBinLabel(9,"TOB 5");
    histo->GetXaxis()->SetBinLabel(10,"TOB 6");
    histo->GetXaxis()->SetBinLabel(11,"TID 1");
    histo->GetXaxis()->SetBinLabel(12,"TID 2");
    histo->GetXaxis()->SetBinLabel(13,"TEC 1");
    histo->GetXaxis()->SetBinLabel(14,"TEC 2");
    histo->GetXaxis()->SetBinLabel(15,"TEC 3");
    histo->GetXaxis()->SetBinLabel(16,"TEC 4");
    histo->GetXaxis()->SetBinLabel(17,"TEC 5");
    histo->GetXaxis()->SetBinLabel(18,"TEC 6");
    histo->GetXaxis()->SetBinLabel(19,"TEC 7");
    histo->GetXaxis()->SetBinLabel(20,"TEC 8");
    histo->GetXaxis()->SetBinLabel(21,"TEC 9");
    histo->GetXaxis()->LabelsOption("v");

	histo->GetYaxis()->SetBinLabel(1,"IB1");
    histo->GetYaxis()->SetBinLabel(2,"IB2");
    histo->GetYaxis()->SetBinLabel(3,"OB1");
    histo->GetYaxis()->SetBinLabel(4,"OB2");
    histo->GetYaxis()->SetBinLabel(5,"W1A");
    histo->GetYaxis()->SetBinLabel(6,"W2A");
    histo->GetYaxis()->SetBinLabel(7,"W3A");
    histo->GetYaxis()->SetBinLabel(8,"W1B");
    histo->GetYaxis()->SetBinLabel(9,"W2B");
    histo->GetYaxis()->SetBinLabel(10,"W3B");
    histo->GetYaxis()->SetBinLabel(11,"W4");
    histo->GetYaxis()->SetBinLabel(12,"W5");
    histo->GetYaxis()->SetBinLabel(13,"W6");
    histo->GetYaxis()->SetBinLabel(14,"W7");
     
}

void SetHistoLabel(TCanvas* canvas,TH1D* histo)
{

    int Wsize = histo->GetMaximum()*1.05;
    histo->GetYaxis()->SetRangeUser(0,histo->GetMaximum()*1.05);

    TLine* line1 = new TLine(4,0,4,0.25);
    TLine* line2 = new TLine(10,0,10,0.25);
    TLine* line3 = new TLine(12,0,12,0.25);

    if(Wsize!=0)
    {
        line1->SetY2(Wsize);
        line2->SetY2(Wsize);
        line3->SetY2(Wsize);
    }


    line1->SetLineStyle(2);
    line2->SetLineStyle(2);
    line3->SetLineStyle(2);
    
    histo->GetXaxis()->SetBinLabel(1,"TIB 1");
    histo->GetXaxis()->SetBinLabel(2,"TIB 2");
    histo->GetXaxis()->SetBinLabel(3,"TIB 3");
    histo->GetXaxis()->SetBinLabel(4,"TIB 4");
    histo->GetXaxis()->SetBinLabel(5,"TOB 1");
    histo->GetXaxis()->SetBinLabel(6,"TOB 2");
    histo->GetXaxis()->SetBinLabel(7,"TOB 3");
    histo->GetXaxis()->SetBinLabel(8,"TOB 4");
    histo->GetXaxis()->SetBinLabel(9,"TOB 5");
    histo->GetXaxis()->SetBinLabel(10,"TOB 6");
    histo->GetXaxis()->SetBinLabel(11,"TID 1");
    histo->GetXaxis()->SetBinLabel(12,"TID 2");
    histo->GetXaxis()->SetBinLabel(13,"TEC 1");
    histo->GetXaxis()->SetBinLabel(14,"TEC 2");
    histo->GetXaxis()->SetBinLabel(15,"TEC 3");
    histo->GetXaxis()->SetBinLabel(16,"TEC 4");
    histo->GetXaxis()->SetBinLabel(17,"TEC 5");
    histo->GetXaxis()->SetBinLabel(18,"TEC 6");
    histo->GetXaxis()->SetBinLabel(19,"TEC 7");
    histo->GetXaxis()->SetBinLabel(20,"TEC 8");
    histo->GetXaxis()->SetBinLabel(21,"TEC 9");
    histo->GetXaxis()->LabelsOption("v");
    canvas->cd();
    histo->Draw();
    line1->Draw();
    line2->Draw();
    line3->Draw();
    
}


string LabelModulGeom(int modulgeom)
{
	if(modulgeom==1) return "IB1";
	if(modulgeom==2) return "IB2";
	if(modulgeom==3) return "OB1";
	if(modulgeom==4) return "OB2";
	if(modulgeom==5) return "W1A";
	if(modulgeom==6) return "W2A";
	if(modulgeom==7) return "W3A";
	if(modulgeom==8) return "W1B";
	if(modulgeom==9) return "W2B";
	if(modulgeom==10) return "W3B";
	if(modulgeom==11) return "W4";
	if(modulgeom==12) return "W5";
	if(modulgeom==13) return "W6";
	if(modulgeom==14) return "W7";
}

string LabelLayer(int i)
{
    if(i==1)   return "TIB1";
    if(i==2)   return "TIB2";
    if(i==3)   return "TIB3";
    if(i==4)   return "TIB4";
    if(i==5)   return "TOB1";
    if(i==6)   return "TOB2";
    if(i==7)   return "TOB3";
    if(i==8)   return "TOB4";
    if(i==9)   return "TOB5";
    if(i==10)  return "TOB6";
    if(i==11)  return "TID1";
    if(i==12)  return "TID2";
    if(i==13)  return "TEC1";
    if(i==14)  return "TEC2";
    if(i==15)  return "TEC3";
    if(i==16)  return "TEC4";
    if(i==17)  return "TEC5";
    if(i==18)  return "TEC6";
    if(i==19)  return "TEC7";
    if(i==20)  return "TEC8";
    if(i==21)  return "TEC9";
}




vector<int> CrossTalkInv(Correction corr,const std::vector<int>&  Q, const float x1, const float x2, float threshold, float thresholdSat, int label, int correctType, float RatioSat254, float RatioSat255, float ThresholdRatioSat) //Methode de correction + XTalkInv
{
  
    bool DylanCorr=true;
    bool inversion=true;
    if(correctType==1){DylanCorr=false;inversion=false;}
    if(correctType==2){DylanCorr=false;}
    if(correctType==3){inversion=false;}

    float charge=0.;

    const unsigned N=Q.size();
    std::vector<int> QII;
    std::vector<float> QI(N,0);
    Double_t a=1-2*x1-2*x2;
    //  bool debugbool=false;
    TMatrix A(N,N);

    float RatioSatTotal = RatioSat254+RatioSat255;

    vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());

    //----------------------
  
    if(!DylanCorr)
    { 
        //---  que pour 1 max bien net 
        if(Q.size()<2 || Q.size()>8)
        {
    	    for (unsigned int i=0;i<Q.size();i++)
            {
    	    	QII.push_back(Q[i]);
      	    }
            return QII;
        }
        
    	if(*mQ>253){
            
    		if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253) return Q;

    		if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40 )
            //if(*(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40)
            //if(RatioSatTotal>=ThresholdRatioSat)
            {
                QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2);
                return QII;
            }
    	}
    }

    //----------------------

    if(DylanCorr && RatioSatTotal>=ThresholdRatioSat && *mQ>253 && (*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat))
    {
        int NSat254=0;
        int NSat255=0;
        float ClusterCharge=0;

        for(int i=0;i<Q.size();i++)
        {     
            if(Q[i]==254) NSat254++;
            if(Q[i]==255) NSat255++;
            ClusterCharge+=(Q[i]);
        }

        float ClusterChargeCorr=corr.ChargeCorr(ClusterCharge,label,Q.size(),NSat254,NSat255)*1.038;
        float DiffClusterCharge=ClusterChargeCorr-ClusterCharge;
        for (unsigned int i=0;i<Q.size();i++)
        {
    	    QII.push_back(Q[i]);
      	}
        //Incomplet pour le moment --> a voir au moment de ClusterCleaning
        //Le surplus de charge ne se fait pas que sur une piste 

        vector<int>::iterator maxQ = max_element(QII.begin(), QII.end());
        if(DiffClusterCharge>=0) QII.at(std::distance(QII.begin(),maxQ))+=DiffClusterCharge;
        bool test=false;
        if(NSat254+NSat255==1 && abs(*(mQ-1) - *(mQ+1)) < 40 ) test=true;
        if(NSat254+NSat255>1 && abs(*(mQ-1) - *(mQ+1)) < 40 ) test=true;
        if(test) return QII;
        if(!test) return Q;
    }

    //----------------------

    if(inversion)
    {
        for(unsigned int i=0; i<N; i++) 
        {
            A(i,i)=a;
            if(i<N-1){ A(i+1,i)=x1;A(i,i+1)=x1;}
            else continue; 
            if(i<N-2){ A(i+2,i)=x2;A(i,i+2)=x2;}
        }
    
        if(N==1) A(0,0)=1/a;
        else  A.InvertFast();

        for(unsigned int i=0; i<N; i++) 
        {
            for(unsigned int j=0; j<N; j++) 
            {
                QI[i]+=A(i,j)*(float)Q[j];
            }
        }
    
        for (unsigned int i=0;i<QI.size();i++)
        {
            if(QI[i]<threshold) QI[i]=0; 
            QII.push_back((int) QI[i]);
        }
        return QII;
    
    }

    if(!inversion) return Q;

}

vector<int> CrossTalkInvStudy(Correction corr,const std::vector<int>&  Q, const float x1, const float x2, float threshold, float thresholdSat, float thresholdDiff, int label, int correctType, float RatioSat254, float RatioSat255, float ThresholdRatioSat, bool& testThresholdSat, bool& testThresholdDiff, bool& testBothThreshold) //Methode de correction + XTalkInv
{
    testThresholdSat=false;
    testThresholdDiff=false;
    testBothThreshold=false;
  
    bool DylanCorr=true;
    bool inversion=true;
    if(correctType==1){DylanCorr=false;inversion=false;}
    if(correctType==2){DylanCorr=false;}
    if(correctType==3){inversion=false;}

    float charge=0.;

    const unsigned N=Q.size();
    std::vector<int> QII;
    std::vector<float> QI(N,0);
    Double_t a=1-2*x1-2*x2;
    //  bool debugbool=false;
    TMatrix A(N,N);

    float RatioSatTotal = RatioSat254+RatioSat255;

    vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());

    //----------------------
  
    if(!DylanCorr)
    { 
        //---  que pour 1 max bien net 
        if(Q.size()<2 || Q.size()>8)
        {
    	    for (unsigned int i=0;i<Q.size();i++)
            {
    	    	QII.push_back(Q[i]);
      	    }
            return QII;
        }
        
    	if(*mQ>253){
            
    		if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253) return Q;

            if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 )
            {
                testThresholdSat=true;
            } 

            if(abs(*(mQ-1) - *(mQ+1)) < thresholdDiff && *(mQ-1)<254 && *(mQ+1)<254 )
            {
                testThresholdDiff=true;
            } 

    		if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < thresholdDiff )
            {
                testBothThreshold=true;
            }
            if(testThresholdSat )
            {
                QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2);
                return QII;
            }
            if(testThresholdDiff )
            {
                QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2);
                return QII;
            }
            if(testBothThreshold)
            {
                QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2);
                return QII;
            }

    	}
    }

    //----------------------

    if(DylanCorr && RatioSatTotal>=ThresholdRatioSat && *mQ>253 && (*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < thresholdDiff ))
    {
        int NSat254=0;
        int NSat255=0;
        float ClusterCharge=0;

        for(int i=0;i<Q.size();i++)
        {     
            if(Q[i]==254) NSat254++;
            if(Q[i]==255) NSat255++;
            ClusterCharge+=(Q[i]);
        }

        //float ClusterChargeCorr=corr.ChargeCorr(ClusterCharge*(3.61*pow(10,-9)*247),label,Q.size(),NSat254,NSat255)/(3.61*pow(10,-9)*247);
        float ClusterChargeCorr=corr.ChargeCorr(ClusterCharge,label,Q.size(),NSat254,NSat255);
        float DiffClusterCharge=ClusterChargeCorr-ClusterCharge;
        for (unsigned int i=0;i<Q.size();i++)
        {
    	    QII.push_back(Q[i]);
      	}
        //Incomplet pour le moment --> a voir au moment de ClusterCleaning
        //Le surplus de charge ne se fait pas que sur une piste 

        vector<int>::iterator maxQ = max_element(QII.begin(), QII.end());
        if(DiffClusterCharge>=0) QII.at(std::distance(QII.begin(),maxQ))+=DiffClusterCharge;
        return QII;
    }


    if(!inversion) return Q;

    //----------------------

    if(inversion)
    {
        for(unsigned int i=0; i<N; i++) 
        {
            A(i,i)=a;
            if(i<N-1){ A(i+1,i)=x1;A(i,i+1)=x1;}
            else continue; 
            if(i<N-2){ A(i+2,i)=x2;A(i,i+2)=x2;}
        }
    
        if(N==1) A(0,0)=1/a;
        else  A.InvertFast();

        for(unsigned int i=0; i<N; i++) 
        {
            for(unsigned int j=0; j<N; j++) 
            {
                QI[i]+=A(i,j)*(float)Q[j];
            }
        }
    
        for (unsigned int i=0;i<QI.size();i++)
        {
            if(QI[i]<threshold) QI[i]=0; 
            QII.push_back((int) QI[i]);
        }
        return QII;
    
    }


}

TH1F* ProfileMPV(TH2F* histo_,float &p0,float &p1,float &Chi2)
{
    int firstdiv = histo_->GetXaxis()->GetBinCenter(histo_->GetNbinsX())/2;
    int divider = histo_->GetXaxis()->GetBinCenter(histo_->GetNbinsX())/firstdiv;
    TH1F* historec_ = new TH1F(histo_->GetTitle(),histo_->GetTitle(),histo_->GetNbinsX()/divider,0,histo_->GetXaxis()->GetBinCenter(histo_->GetNbinsX()));
    TH2F* histoclone = (TH2F*) histo_->Clone();
    histoclone->Reset();
    int binrec=0;

    for(int x=1;x<=histo_->GetNbinsX();x++)
    {
        for(int y=1;y<=histo_->GetNbinsY();y++)
        {
            if(histo_->GetBinContent(x,y)>0) histoclone->SetBinContent(x,y,histo_->GetBinContent(x,y));
        }
        if(x%2==0)
        {
            binrec++;
            float mpv=1.;
            float mpverror=0.;
            TH1F* projY = (TH1F*) histoclone->ProjectionY();
            if(projY->GetEntries()>0) 
            {
                int fitStatus = projY->Fit("gaus","QL","");
                TFitResultPtr func_fit = projY->Fit("gaus","QLS","");
                if(fitStatus>-1) 
                {
                    mpv=func_fit->Parameter(1);
                    mpverror=func_fit->Error(1);
                }
                historec_->SetBinContent(binrec,mpv);
                historec_->SetBinError(binrec,mpverror);
                //historec_->SetBinError(binrec,sigma);
            }
            histoclone->Reset();
            delete projY;
        }
    }
    delete histoclone;
    TFitResultPtr FitRes_ = historec_->Fit("pol1","SQR","",500,1500);
    int fitStat = historec_->Fit("pol1","QR","",500,1500);
    if(fitStat>-1)
    {
        float p0_     = FitRes_->Parameter(0);
        float p1_     = FitRes_->Parameter(1);
        float p0err_  = FitRes_->Error(0);
        float p1err_  = FitRes_->Error(1);
        float chi2_   = FitRes_->Chi2();
        float Ndf_    = FitRes_->Ndf();
        float Chi2overNdf_ = (float)chi2_/(float)Ndf_;
        float chicalc=0.;
        int nentries=0;
        for(int bin=0;bin<historec_->GetNbinsX();bin++)
        {
            if(historec_->GetBinContent(bin)>0)
            {
                nentries++;
                float valueCalc=(historec_->GetXaxis()->GetBinCenter(bin)*p1_+p0_);
                chicalc+=pow(valueCalc-historec_->GetBinContent(bin),2);
            }
            
        }
        chicalc=(sqrt(chicalc)/nentries);
        Chi2overNdf_=chicalc;
        p0=p0_;
        p1=p1_;
        Chi2=Chi2overNdf_;
    }
    return historec_;
}

TH1F* ProfilerMPV(TH2F* histo_)
{
    int firstdiv = histo_->GetXaxis()->GetBinCenter(histo_->GetNbinsX())/2;
    int divider = histo_->GetXaxis()->GetBinCenter(histo_->GetNbinsX())/firstdiv;
    TH1F* historec_ = new TH1F(histo_->GetTitle(),histo_->GetTitle(),histo_->GetNbinsX()/divider,0,histo_->GetXaxis()->GetBinCenter(histo_->GetNbinsX()));
    TH2F* histoclone = (TH2F*) histo_->Clone();
    histoclone->Reset();
    int binrec=0;

    for(int x=1;x<=histo_->GetNbinsX();x++)
    {
        for(int y=1;y<=histo_->GetNbinsY();y++)
        {
            if(histo_->GetBinContent(x,y)>0) histoclone->SetBinContent(x,y,histo_->GetBinContent(x,y));
        }
        if(x%2==0)
        {
            binrec++;
            float mpv=1.;
            float sigma=0.;
            float mpverror=0.;
            float sigmaerror=0.;
            TH1F* projY = (TH1F*) histoclone->ProjectionY();
            if(projY->GetEntries()>0) 
            {
                int fitStatus = projY->Fit("gaus","QR","",projY->GetMean()-3*projY->GetStdDev(),projY->GetMean()+3*projY->GetStdDev());
                TFitResultPtr func_fit = projY->Fit("gaus","QRS","",projY->GetMean()-3*projY->GetStdDev(),projY->GetMean()+3*projY->GetStdDev());
                if(fitStatus==0) 
                {
                    mpv=func_fit->Parameter(1);
                    sigma=func_fit->Parameter(2);
                    mpverror=func_fit->Error(1);
                    sigmaerror=func_fit->Error(2);
                }
                historec_->SetBinContent(binrec,mpv);
                historec_->SetBinError(binrec,mpverror);
            }
            histoclone->Reset();
            delete projY;
        }
    }
    delete histoclone;
    return historec_;
}

TH1F* ProfilerResolution(TH2F* histo_)
{
    int firstdiv = histo_->GetXaxis()->GetBinCenter(histo_->GetNbinsX())/2;
    int divider = histo_->GetXaxis()->GetBinCenter(histo_->GetNbinsX())/firstdiv;
    TH1F* historec_ = new TH1F(histo_->GetTitle(),histo_->GetTitle(),histo_->GetNbinsX()/divider,0,histo_->GetXaxis()->GetBinCenter(histo_->GetNbinsX()));
    TH2F* histoclone = (TH2F*) histo_->Clone();
    histoclone->Reset();
    int binrec=0;

    for(int x=1;x<=histo_->GetNbinsX();x++)
    {
        for(int y=1;y<=histo_->GetNbinsY();y++)
        {
            if(histo_->GetBinContent(x,y)>0) histoclone->SetBinContent(x,y,histo_->GetBinContent(x,y));
        }
        if(x%2==0)
        {
            binrec++;
            float mpv=1.;
            float sigma=0.;
            float mpverror=0.;
            float sigmaerror=0.;
            TH1F* projY = (TH1F*) histoclone->ProjectionY();
            if(projY->GetEntries()>0) 
            {
                int fitStatus = projY->Fit("gaus","QR","",projY->GetMean()-3*projY->GetStdDev(),projY->GetMean()+3*projY->GetStdDev());
                TFitResultPtr func_fit = projY->Fit("gaus","QRS","",projY->GetMean()-3*projY->GetStdDev(),projY->GetMean()+3*projY->GetStdDev());
                if(fitStatus==0) 
                {
                    mpv=func_fit->Parameter(1);
                    sigma=func_fit->Parameter(2);
                    mpverror=func_fit->Error(1);
                    sigmaerror=func_fit->Error(2);
                }
                // sigma=projY->GetStdDev();
                // sigmaerror=projY->GetStdDevError();
                historec_->SetBinContent(binrec,sigma);
                historec_->SetBinError(binrec,sigmaerror);
            }
            histoclone->Reset();
            delete projY;
        }
    }
    delete histoclone;
    return historec_;
}

TH1F* FactorKC(TH2F* histo,float& K,float& C) //histo pre-filtre sur le particleID --> il faut donner la masse correspondante 
{
    TH1F* historec = (TH1F*) ProfilerMPV(histo);
    TF1* InvSquare = new TF1("InvSquare","[0]+[1]*(1/pow(x,2))");
    TFitResultPtr fit = historec->Fit("InvSquare","QRS","",0,3000);
    int fitstatus = historec->Fit("InvSquare","QR","",0,3000);
    vector<float> vectres;
    if(fitstatus==0)
    {
        vectres.push_back(fit->Parameter(1)); //K
        vectres.push_back(fit->Parameter(0)); //C
        vectres.push_back(fit->Error(1));
        vectres.push_back(fit->Error(0));   
    }
    else{vectres.push_back(0);vectres.push_back(0);vectres.push_back(0);vectres.push_back(0);}
    K=vectres[0];C=vectres[2];
    return historec;
}

void FitProfile(TProfile* profile,float &p0,float &p1)
{
    TFitResultPtr FitRes_ = profile->Fit("pol1","SQR","",400,2000);
    int fitStat = profile->Fit("pol1","QR","",400,2000);
    if(fitStat>-1)
    {
        float p0_     = FitRes_->Parameter(0);
        float p1_     = FitRes_->Parameter(1);
        float p0err_  = FitRes_->Error(0);
        float p1err_  = FitRes_->Error(1);
        float chi2_   = FitRes_->Chi2();
        float Ndf_    = FitRes_->Ndf();
        float Chi2overNdf_ = (float)chi2_/(float)Ndf_;
        float chicalc=0.;
        int nentries=0;
        for(int bin=0;bin<profile->GetNbinsX();bin++)
        {
            if(profile->GetBinContent(bin)>0)
            {
                nentries++;
                float valueCalc=(profile->GetXaxis()->GetBinCenter(bin)*p1_+p0_);
                chicalc+=pow(valueCalc-profile->GetBinContent(bin),2);
            }
            
        }
        chicalc=(sqrt(chicalc)/nentries);
        Chi2overNdf_=chicalc;
        p0=p1_;
        p1=p1err_;
    }
}

float DeDxMoyen(float p,float K,float C)
{
    float PoverM=p/2400;
    return K*pow(PoverM,-2)+C;
}

float QtoE(float Q) 
{
    return Q*(3.61*pow(10,-9)*247);
}

TCanvas* DrawTransparency(vector<TH1F*> histo1,TH1F* histo2,string title,string xaxis,string yaxis)
{
    TCanvas* c1 = new TCanvas(title.c_str(),"");
    histo1.at(0)->GetXaxis()->SetTitle(xaxis.c_str());
    histo1.at(0)->GetYaxis()->SetTitle(yaxis.c_str());
    c1->cd();
    for(int i=0;i<histo1.size();i++)
    {
        histo1.at(i)->SetLineColor(2*(i+1));
        histo1.at(i)->Draw("same");
    }
    histo2->Scale(0.05);
    histo2->SetFillStyle(3003);
    histo2->SetFillColor(1);
    histo2->SetLineColor(1);
    histo2->Draw("same,bar");
    return c1;
}





//MAIN



int main(int argc,char** argv)
{

	TChain chain;
	chain.SetName("stage/ttree");
	for(int i=1;i<argc-1;i++) chain.Add(argv[i]); 
    string s1 = argv[1];
	string s2 = s1.substr(0,s1.find('.'))+"_results.root";

    Builder* b1 = new Builder(chain);
    b1->SetBranchAdd();


	int nentries = b1->GetEntries();


//OUVERTURE FICHIERS DE CONFIGS POUR LA METHODE

	Correction EcorrLayer;
	TFile* fileLayer = TFile::Open("testMethodLayer.root");
    //TFile* fileLayer = TFile::Open("testMethodModulGeom.root");
	TTree* treeLayer = (TTree*) fileLayer->Get("tree");
	EcorrLayer.SetTree(*treeLayer);
    EcorrLayer.ReadLayer();
    //EcorrLayer.ReadModulGeom();

    Correction EcorrModulGeom;
	TFile* fileModulGeom = TFile::Open("testMethodModulGeom.root");
    //TFile* fileModulGeom = TFile::Open("testMethodModulGeom.root");
	TTree* treeModulGeom = (TTree*) fileModulGeom->Get("tree");
	EcorrModulGeom.SetTree(*treeModulGeom);
    EcorrModulGeom.ReadModulGeom();


    SaturationCorrection sc;
    TFile* fileparameters = TFile::Open("FileParameters.root");
    TTree* treeparameters = (TTree*) fileparameters->Get("tree");

    SaturationCorrection sc2;
    sc2.SetTree(*treeparameters);
    sc2.ReadParameters();


	/*Correction EcorrModulGeom;
	TFile* fileModulGeom = TFile::Open("testMethodModulGeom.root");
	TTree* treeModulGeom = (TTree*) fileModulGeom->Get("tree");
	EcorrModulGeom.SetTree(*treeModulGeom);
	Correction EcorrModulGeomNoSat;
	TFile* fileModulGeomNoSat = TFile::Open("testMethodModulGeomNoSat.root");
	TTree* treeModulGeomNoSat = (TTree*) fileModulGeomNoSat->Get("tree");
	EcorrModulGeomNoSat.SetTree(*treeModulGeomNoSat);*/

// OUVERTURE FICHIER CONFIG FACTEURS K ET C

    ifstream inputfileFactorKC ("factKC.txt");
    float K_Esim, C_Esim;
    float K_Erec, C_Erec;
    float K_Ecorr, C_Ecorr;
    float K_EfullLayer, C_EfullLayer;
    float K_EfullModulGeom, C_EfullModulGeom;

    inputfileFactorKC>>K_Esim; inputfileFactorKC>>C_Esim;
    inputfileFactorKC>>K_Erec; inputfileFactorKC>>C_Erec;
    inputfileFactorKC>>K_Ecorr; inputfileFactorKC>>C_Ecorr;
    inputfileFactorKC>>K_EfullLayer; inputfileFactorKC>>C_EfullLayer;
    inputfileFactorKC>>K_EfullModulGeom; inputfileFactorKC>>C_EfullModulGeom;
    

    float Calibp0=0., Calibp1=0.;
    ifstream inputfileCalib("calib.txt");
    inputfileCalib>>Calibp0;
    inputfileCalib>>Calibp1;


	bool method=false;
	bool corr=false;
    bool FactKC=false;
    bool XTalkInvStudy=false;

    bool gluino=false;

    float idlow=2212,idup=2212;

    float ratsat=0.;

//CHOIX OPTION



	if(atof(argv[argc-1])==1) method=true;
	if(atof(argv[argc-1])==2) corr=true;
    if(atof(argv[argc-1])==3) {corr=true;FactKC=true;}
    if(atof(argv[argc-1])==4) XTalkInvStudy=true;


    if(atof(argv[argc-2])==1) {gluino=true;idlow=1000001;idup=2000015;}

    ratsat=atof(argv[argc-3]);
    

//------------------------------------------------------
//
//
//          Le premier 1 sert a prendre les gluinos
//
//
//
//          Option 1 : creer les fichiers avec les parametres de correction 
//          Option 2 : correction de l'energie
//          Option 3 : enregistre les facteurs K et C (necessite en plus l'option 2)
//          Option 4 : CrossTalkInvAlgo study
//
//
//
//
//
//------------------------------------------------------




//DECLARATION TABLEAUX POUR LES METHODES 

vector<vector<vector<vector<TH2F*>>>> VectLayerVectNStripVectNStripSat254VectNStripSat255Histo;
vector<vector<vector<vector<TH2F*>>>> VectModulGeomVectNStripVectNStripSat254VectNStripSat255Histo;
vector<vector<vector<vector<TH2F*>>>> VectModulGeomVectNStripVectNStripSat254VectNStripSat255HistoNoSat;
if(method)
{

    sc.InitHisto();




//METHODE 1 -- CATEGORIES LAYER
    for(int layer=1;layer<22;layer++)
    {
        vector<vector<vector<TH2F*>>> VectNStripVectNStripSat254VectNStripSat255Histo;
        for(int nstrip=3;nstrip<7;nstrip++)
        {
            vector<vector<TH2F*>> VectNStripSat254VectNStripSat255Histo;
            for(int nstripsat254=0;nstripsat254<nstrip+1;nstripsat254++)
            {
                vector<TH2F*> VectNStripSat255Histo;
                for(int nstripsat255=0;nstripsat255<nstrip-nstripsat254+1;nstripsat255++)
                {
                    string title = LabelLayer(layer)+" NStrip="+to_string(nstrip)+" NStripSat254="+to_string(nstripsat254)+" NStripSat255="+to_string(nstripsat255);
                    //VectNStripSat255Histo.push_back(new TH2F(title.c_str(),title.c_str(),300,0,6000*pow(10,-6),300,0,6000*pow(10,-6)));
                    VectNStripSat255Histo.push_back(new TH2F(title.c_str(),title.c_str(),200,0,5000,200,0,5000));
                }
                VectNStripSat254VectNStripSat255Histo.push_back(VectNStripSat255Histo);
            }
            VectNStripVectNStripSat254VectNStripSat255Histo.push_back(VectNStripSat254VectNStripSat255Histo);
        }
        VectLayerVectNStripVectNStripSat254VectNStripSat255Histo.push_back(VectNStripVectNStripSat254VectNStripSat255Histo);
    }

//METHODE 2 -- CATEGORIES XTALK

	for(int ModulGeom=1;ModulGeom<15;ModulGeom++)
    {
        vector<vector<vector<TH2F*>>> VectNStripVectNStripSat254VectNStripSat255Histo;
        for(int nstrip=3;nstrip<7;nstrip++)
        {
            vector<vector<TH2F*>> VectNStripSat254VectNStripSat255Histo;
            for(int nstripsat254=0;nstripsat254<nstrip+1;nstripsat254++)
            {
                vector<TH2F*> VectNStripSat255Histo;
                for(int nstripsat255=0;nstripsat255<nstrip-nstripsat254+1;nstripsat255++)
                {
                    string title = LabelModulGeom(ModulGeom)+" NStrip="+to_string(nstrip)+" NStripSat254="+to_string(nstripsat254)+" NStripSat255="+to_string(nstripsat255);
                    //VectNStripSat255Histo.push_back(new TH2F(title.c_str(),title.c_str(),300,0,6000*pow(10,-6),300,0,6000*pow(10,-6)));
                    VectNStripSat255Histo.push_back(new TH2F(title.c_str(),title.c_str(),200,0,5000,200,0,5000));
                }
                VectNStripSat254VectNStripSat255Histo.push_back(VectNStripSat255Histo);
            }
            VectNStripVectNStripSat254VectNStripSat255Histo.push_back(VectNStripSat254VectNStripSat255Histo);
        }
        VectModulGeomVectNStripVectNStripSat254VectNStripSat255Histo.push_back(VectNStripVectNStripSat254VectNStripSat255Histo);
    }
    vector<vector<vector<TH2F*>>> vModulGeomvNStripvNSat;
    for(int ModulGeom=1;ModulGeom<15;ModulGeom++)
    {
        vector<vector<TH2F*>> vNStripvNSat;
        for(int nstrip=3;nstrip<7;nstrip++)
        {
            vector<TH2F*> vNSat;
            for(int nsat=1;nsat<3;nsat++)
            {
                string title = LabelModulGeom(ModulGeom)+" NStrip="+to_string(nstrip)+" NSat="+to_string(nsat);
                vNSat.push_back(new TH2F(title.c_str(),title.c_str(),200,0,5000,200,0,5000));
            }
            vNStripvNSat.push_back(vNSat);
        }
        vModulGeomvNStripvNSat.push_back(vNStripvNSat);
    }

}


//DECLARATION HISTOGRAMMES POUR ETUDES

TH1I* hmodulgeom = new TH1I("modulgeom","",25,0,25);
TH1I* hlayerlabel = new TH1I("layerlabel","",25,0,25);

TH2I* hmodulgeomvslayer = new TH2I("hmodulgeomvslayer","",21,0,21,14,0,14);

TH1F* hEcorrLayer = new TH1F("hEcorrLayer","",100,-2,2);
TH1F* hEcorrModulGeom = new TH1F("hEcorrModulGeom","",100,-2,2);
TH1F* hEcorrModulGeomNoSat = new TH1F("hEcorrModulGeomNoSat","",100,-2,2);

TH1F* hErecLayer = new TH1F("hErecLayer","",100,-2,2);
TH1F* hErecModulGeom = new TH1F("hErecModulGeom","",100,-2,2);
TH1F* hErecModulGeomNoSat = new TH1F("hErecModulGeomNoSat","",100,-2,2);

TH2F* h2ElossEcorrSatLayer = new TH2F("h2ElossEcorrSatLayer","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH2F* h2ElossEcorrModulGeom = new TH2F("h2ElossEcorrModulGeom","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH2F* h2ElossEcorrModulGeomNoSat = new TH2F("h2ElossEcorrModulGeomNoSat","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));

TH2F* h2ElossErec = new TH2F("h2ElossErec","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));

TH2F* h2ElossErecLayer = new TH2F("h2ElossErecLayer","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH2F* h2ElossErecModulGeom = new TH2F("h2ElossErecModulGeom","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH2F* h2ElossErecModulGeomNoSat = new TH2F("h2ElossErecModulGeomNoSat","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));

TH2F* h2EcorrLayerVsEcorrModulGeom = new TH2F("h2EcorrLayerVsEcorrModulGeom","",300,0,6000*pow(10,-6),300,0,6000*pow(10,-6));
TH1F* h1EcorrLayerVsEcorrModulGeom = new TH1F("h1EcorrLayerVsEcorrModulGeom","",100,-0.1,0.1);

TH2F* h2PoverMDecorrDxLayer = new TH2F("h2PoverMDecorrDxLayer","",100,0,10,50,0,50);
TH2F* h2PoverMDecorrDxModulGeom = new TH2F("h2PoverMDecorrDxModulGeom","",100,0,10,50,0,50);

TH2F* h2PvsDecorrDxLayer = new TH2F("h2PvsDecorrDxLayer","",100,0,3000,50,0,50);
TH2F* h2PvsDecorrDxModulGeom = new TH2F("h2PvsDecorrDxModulGeom","",100,0,3000,50,0,50);

TH2F* h2PvsDelossDxProtonKC = new TH2F("h2PvsDelossDxProtonKC","",1000,0,5,1000,0,15);
TH2F* h2PvsDecorrDxProtonKC = new TH2F("h2PvsDecorrDxProtonKC","",1000,0,5,1000,0,15);
TH2F* h2PvsDecorrDxProtonKCnewSatLayer = new TH2F("h2PvsDecorrDxProtonKCnewSatLayer","",1000,0,5,1000,0,15);

TH1F* historeceloss = new TH1F("","",1000,0,5);
TH1F* historececorr = new TH1F("","",1000,0,5);
TH1F* historececorrnew = new TH1F("","",1000,0,5);


TH1F* distribmassprotonEloss = new TH1F("distribmassprotonEloss","",100,0,5000);
TH1F* distribmassprotonEcorr = new TH1F("distribmassprotonEcorr","",100,0,5000);
TH1F* distribmassprotonEcorrNew = new TH1F("distribmassprotonEcorrNew","",100,0,5000);


TH2F* pVsdEsimdx = new TH2F("pVsdEsimdx","",1000,0,5,1000,0,15);
TH2F* pVsdErecdx = new TH2F("pVsdErecdx","",1000,0,5,1000,0,15);
TH2F* pVsdEcorrdx = new TH2F("pVsdEcorrdx","",1000,0,5,1000,0,15);
TH2F* pVsdEfulldx_Layer = new TH2F("pVsdEfulldx_Layer","",1000,0,5,1000,0,15);
TH2F* pVsdEfulldx_ModulGeom = new TH2F("pVsdEfulldx_ModulGeom","",1000,0,5,1000,0,15);


float rangemass=2;
if(gluino) rangemass=5000;

TH1F* massEsim = new TH1F("massEsim","",100,0,rangemass);
TH1F* massErec = new TH1F("massErec","",100,0,rangemass);
TH1F* hmassEcorr = new TH1F("massEcorr","",100,0,rangemass);
TH1F* massEfullLayer = new TH1F("massEfullLayer","",100,0,rangemass);
TH1F* massEfullModulGeom = new TH1F("massEfullModulGeom","",100,0,rangemass);

//Distrib d'interet (saturation)

TH1F* h_ClustSat_DistribP = new TH1F("h_ClustSat_DistribP","",100,0,10000);
TH1F* h_ClustSat_DistribEta = new TH1F("h_ClustSat_DistribEta","",100,0,3);
TH1I* h_ClustSat_DistribSize = new TH1I("h_ClustSat_DistribSize","",8,0,8);
TH1I* h_ClustSat_DistribNSat = new TH1I("h_ClustSat_DistribNsat","",5,0,5);

TH1F* h_DistribP = new TH1F("h_DistribP","",100,0,10000);
TH1F* h_DistribEta = new TH1F("h_DistribEta","",100,0,3);
TH1I* h_DistribSize = new TH1I("h_DistribSize","",8,0,8);
TH1I* h_DistribNSat = new TH1I("h_DistribNsat","",5,0,5);



//Etude de la methode actuelle + XTalkInv + comparaison avec DylanCorr

TH1F* hWithoutCorr = new TH1F("hWithoutCorr","",200,-2,2);
TH1F* hCorrHSCP = new TH1F("hCorrHSCP","",200,-2,2);
TH1F* hCorrHSCP_WithoutInversion = new TH1F("hCorrHSCP_WithoutInversion","",200,-2,2);
TH1F* hMyCorr = new TH1F("hMyCorr","",200,-2,2);
TH1F* hMyCorr_WithoutInversion = new TH1F("MyCorr_WithoutInversion","",200,-2,2);
TH1F* hMyCorr_RatioSat = new TH1F("hMyCorr_RatioSat","",200,-2,2);
TH1F* hMyCorr_WithoutInversion_RatioSat = new TH1F("hMyCorr_WithoutInversion_RatioSat","",200,-2,2);



TH1F* hMyCorrModulGoem = new TH1F("hMyCorrModulGoem","",200,-2,2);



TH1F* h_ClustSat_WithoutCorr = new TH1F("h_ClustSat_WithoutCorr","",200,-2,2);
TH1F* h_ClustSat_MyFullCorr = new TH1F("h_ClustSat_MyFullCorr","",200,-2,2);
TH1F* h_ClustSat_HscpCorrNoInv = new TH1F("h_ClustSat_HscpCorrNoInv","",200,-2,2);

TH1F* h_ClustSatAndHscpTest_WithoutCorr = new TH1F("h_ClustSatAndHscpTest_WithoutCorr","",200,-2,2);
TH1F* h_ClustSatAndHscpTest_MyFullCorr = new TH1F("h_ClustSatAndHscpTest_MyFullCorr","",200,-2,2);
TH1F* h_ClustSatAndHscpTest_HscpCorrNoInv = new TH1F("h_ClustSatAndHscpTest_HscpCorrNoInv","",200,-2,2);
TH1F* h_ClustSatAndHscpTest_MyFullCorrModulGeom = new TH1F("h_ClustSatAndHscpTest_MyFullCorrModulGeom","",200,-2,2);


TH1F* h_ClustSat_CorrHSCP = new TH1F("h_ClustSat_CorrHSCP","",200,-2,2);
TH1F* h_ClustSat_CorrHSCP_WithoutInversion = new TH1F("h_ClustSat_CorrHSCP_WithoutInversion","",200,-2,2);
TH1F* h_ClustSat_MyCorr = new TH1F("h_ClustSat_MyCorr","",200,-2,2);
TH1F* h_ClustSat_MyCorr_WithoutInversion = new TH1F("h_ClustSat_MyCorr_WithoutInversion","",200,-2,2);
TH1F* h_ClustSat_MyCorr_RatioSat = new TH1F("h_ClustSat_MyCorr_RatioSat","",200,-2,2);
TH1F* h_ClustSat_MyCorr_WithoutInversion_RatioSat = new TH1F("h_ClustSat_MyCorr_WithoutInversion_RatioSat","",200,-2,2);

TH1F* h_ClustNoSat_WithoutCorr = new TH1F("h_ClustNoSat_WithoutCorr","",200,-2,2);
TH1F* h_ClustNoSat_CorrHSCP = new TH1F("h_ClustNoSat_CorrHSCP","",200,-2,2);
TH1F* h_ClustNoSat_CorrHSCP_WithoutInversion = new TH1F("h_ClustNoSat_CorrHSCP_WithoutInversion","",200,-2,2);
TH1F* h_ClustNoSat_MyCorr = new TH1F("h_ClustNoSat_MyCorr","",200,-2,2);
TH1F* h_ClustNoSat_MyCorr_WithoutInversion = new TH1F("h_ClustNoSat_MyCorr_WithoutInversion","",200,-2,2);
TH1F* h_ClustNoSat_MyCorr_RatioSat = new TH1F("h_ClustNoSat_MyCorr_RatioSat","",200,-2,2);
TH1F* h_ClustNoSat_MyCorr_WithoutInversion_RatioSat = new TH1F("h_ClustNoSat_MyCorr_WithoutInversion_RatioSat","",200,-2,2);

TH2F* h2_ClustSat_MyCorrWOinv_p = new TH2F("h2_ClustSat_MyCorrWOinv_p","",50,-2,2,50,0,4000);
TH2F* h2_ClustSat_MyCorrWOinv_eta = new TH2F("h2_ClustSat_MyCorrWOinv_eta","",50,-2,2,50,0,2.5);

TH1F* h_ClustSat_MyCorrWOinv_0p500 = new TH1F("h_ClustSat_MyCorrWOinv_0p500","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_500p1000 = new TH1F("h_ClustSat_MyCorrWOinv_500p1000","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_1000p1500 = new TH1F("h_ClustSat_MyCorrWOinv_1000p1500","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_1500p2000 = new TH1F("h_ClustSat_MyCorrWOinv_1500p2000","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_2000p4000 = new TH1F("h_ClustSat_MyCorrWOinv_2000p4000","",50,-2,2);

TH1F* h_ClustSat_MyCorrWOinv_0eta05 = new TH1F("h_ClustSat_MyCorrWOinv_0eta05","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_05eta10 = new TH1F("h_ClustSat_MyCorrWOinv_05eta10","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_10eta15 = new TH1F("h_ClustSat_MyCorrWOinv_10eta15","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_15eta25 = new TH1F("h_ClustSat_MyCorrWOinv_15eta25","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_20eta25 = new TH1F("h_ClustSat_MyCorrWOinv_20eta25","",50,-2,2);

TH1F* h_ClustSat_HscpCorrWOinv_0p500 = new TH1F("h_ClustSat_HscpCorrWOinv_0p500","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_500p1000 = new TH1F("h_ClustSat_HscpCorrWOinv_500p1000","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_1000p1500 = new TH1F("h_ClustSat_HscpCorrWOinv_1000p1500","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_1500p2000 = new TH1F("h_ClustSat_HscpCorrWOinv_1500p2000","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_2000p4000 = new TH1F("h_ClustSat_HscpCorrWOinv_2000p4000","",50,-2,2);

TH1F* h_ClustSat_HscpCorrWOinv_0eta05 = new TH1F("h_ClustSat_HscpCorrWOinv_0eta05","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_05eta10 = new TH1F("h_ClustSat_HscpCorrWOinv_05eta10","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_10eta15 = new TH1F("h_ClustSat_HscpCorrWOinv_10eta15","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_15eta25 = new TH1F("h_ClustSat_HscpCorrWOinv_15eta25","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_20eta25 = new TH1F("h_ClustSat_HscpCorrWOinv_20eta25","",50,-2,2);

//Selon les differentes couches 

TH1F* h_ClustSat_MyCorrWOinv_TIB = new TH1F("h_ClustSat_MyCorrWOinv_TIB","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_TID = new TH1F("h_ClustSat_MyCorrWOinv_TID","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_TOB = new TH1F("h_ClustSat_MyCorrWOinv_TOB","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_TEC = new TH1F("h_ClustSat_MyCorrWOinv_TEC","",50,-2,2);

TH1F* h_ClustSat_HscpCorrWOinv_TIB = new TH1F("h_ClustSat_HscpCorrWOinv_TIB","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_TID = new TH1F("h_ClustSat_HscpCorrWOinv_TID","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_TOB = new TH1F("h_ClustSat_HscpCorrWOinv_TOB","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_TEC = new TH1F("h_ClustSat_HscpCorrWOinv_TEC","",50,-2,2);

TH1F* h_ClustSat_MyCorrWOinv_Layer1 = new TH1F("h_ClustSat_MyCorrWOinv_Layer1","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer2 = new TH1F("h_ClustSat_MyCorrWOinv_Layer2","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer3 = new TH1F("h_ClustSat_MyCorrWOinv_Layer3","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer4 = new TH1F("h_ClustSat_MyCorrWOinv_Layer4","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer5 = new TH1F("h_ClustSat_MyCorrWOinv_Layer5","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer6 = new TH1F("h_ClustSat_MyCorrWOinv_Layer6","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer7 = new TH1F("h_ClustSat_MyCorrWOinv_Layer7","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer8 = new TH1F("h_ClustSat_MyCorrWOinv_Layer8","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer9 = new TH1F("h_ClustSat_MyCorrWOinv_Layer9","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_Layer10 = new TH1F("h_ClustSat_MyCorrWOinv_Layer10","",50,-2,2);

TH1F* h_ClustSat_HscpCorrWOinv_Layer1 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer1","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer2 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer2","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer3 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer3","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer4 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer4","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer5 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer5","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer6 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer6","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer7 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer7","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer8 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer8","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer9 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer9","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_Layer10 = new TH1F("h_ClustSat_HscpCorrWOinv_Layer10","",50,-2,2);

//Selon taille clusters 

TH1F* h_ClustSat_MyCorrWOinv_size1 = new TH1F("h_ClustSat_MyCorrWOinv_size1","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size2 = new TH1F("h_ClustSat_MyCorrWOinv_size2","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size3 = new TH1F("h_ClustSat_MyCorrWOinv_size3","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size4 = new TH1F("h_ClustSat_MyCorrWOinv_size4","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size5 = new TH1F("h_ClustSat_MyCorrWOinv_size5","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size6 = new TH1F("h_ClustSat_MyCorrWOinv_size6","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size7 = new TH1F("h_ClustSat_MyCorrWOinv_size7","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_size8 = new TH1F("h_ClustSat_MyCorrWOinv_size8","",50,-2,2);

TH1F* h_ClustSat_HscpCorrWOinv_size1 = new TH1F("h_ClustSat_HscpCorrWOinv_size1","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size2 = new TH1F("h_ClustSat_HscpCorrWOinv_size2","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size3 = new TH1F("h_ClustSat_HscpCorrWOinv_size3","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size4 = new TH1F("h_ClustSat_HscpCorrWOinv_size4","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size5 = new TH1F("h_ClustSat_HscpCorrWOinv_size5","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size6 = new TH1F("h_ClustSat_HscpCorrWOinv_size6","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size7 = new TH1F("h_ClustSat_HscpCorrWOinv_size7","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_size8 = new TH1F("h_ClustSat_HscpCorrWOinv_size8","",50,-2,2);

//Selon nombre de strips saturees 

TH1F* h_ClustSat_MyCorrWOinv_sat1 = new TH1F("h_ClustSat_MyCorrWOinv_sat1","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_sat2 = new TH1F("h_ClustSat_MyCorrWOinv_sat2","",50,-2,2);
TH1F* h_ClustSat_MyCorrWOinv_sat3 = new TH1F("h_ClustSat_MyCorrWOinv_sat3","",50,-2,2);

TH1F* h_ClustSat_HscpCorrWOinv_sat1 = new TH1F("h_ClustSat_HscpCorrWOinv_sat1","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_sat2 = new TH1F("h_ClustSat_HscpCorrWOinv_sat2","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWOinv_sat3 = new TH1F("h_ClustSat_HscpCorrWOinv_sat3","",50,-2,2);

//Etude sous correction a bas p 

TH1I* h_ClustSat_LowP_ClusterSize = new TH1I("h_ClustSat_LowP_ClusterSize","",8,0,8);
TH1I* h_ClustSat_LowP_Sat254 = new TH1I("h_ClustSat_LowP_Sat254","",5,0,5);
TH1I* h_ClustSat_LowP_Sat255 = new TH1I("h_ClustSat_LowP_Sat255","",5,0,5);
TH1I* h_ClustSat_LowP_NSat = new TH1I("h_ClustSat_LowP_NSat","",5,0,5);
TH1D* h_ClustSat_LowP_layerlabel = new TH1D("h_ClustSat_LowP_layerlabel","",21,0,21);
TH1F* h_ClustSat_LowP_Eta = new TH1F("h_ClustSat_LowP_Eta","",10,0,2.5);
TH1F* h_ClustSat_LowP_Esim = new TH1F("h_ClustSat_LowP_Esim","",300,0,6000*pow(10,-6));
TH1F* h_ClustSat_LowP_EcorrHscp = new TH1F("h_ClustSat_LowP_EcorrHscp","",300,0,6000*pow(10,-6));
TH1F* h_ClustSat_LowP_Emycorr = new TH1F("h_ClustSat_LowP_Emycorr","",300,0,6000*pow(10,-6));

//Etude seuil 1 

TH1F* h_ClustSat_HscpCorrWinv_firstThreshold20 = new TH1F("h_ClustSat_HscpCorrWinv_firstThreshold20","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstThreshold50 = new TH1F("h_ClustSat_HscpCorrWinv_firstThreshold50","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstThreshold100 = new TH1F("h_ClustSat_HscpCorrWinv_firstThreshold100","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstThreshold150 = new TH1F("h_ClustSat_HscpCorrWinv_firstThreshold150","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstThreshold200 = new TH1F("h_ClustSat_HscpCorrWinv_firstThreshold200","",50,-2,2);

//Etude seuil 2

TH1F* h_ClustSat_HscpCorrWinv_secondThreshold25 = new TH1F("h_ClustSat_HscpCorrWinv_secondThreshold25","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondThreshold50 = new TH1F("h_ClustSat_HscpCorrWinv_secondThreshold50","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondThreshold100 = new TH1F("h_ClustSat_HscpCorrWinv_secondThreshold100","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondThreshold150 = new TH1F("h_ClustSat_HscpCorrWinv_secondThreshold150","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondThreshold200 = new TH1F("h_ClustSat_HscpCorrWinv_secondThreshold200","",50,-2,2);

//Etude x-talk1

TH1F* h_ClustSat_HscpCorrWinv_firstXTalk10 = new TH1F("h_ClustSat_HscpCorrWinv_firstXTalk10","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstXTalk11 = new TH1F("h_ClustSat_HscpCorrWinv_firstXTalk11","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstXTalk12 = new TH1F("h_ClustSat_HscpCorrWinv_firstXTalk12","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstXTalk9 = new TH1F("h_ClustSat_HscpCorrWinv_firstXTalk9","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_firstXTalk8 = new TH1F("h_ClustSat_HscpCorrWinv_firstXTalk8","",50,-2,2);

//Etude x-talk2

TH1F* h_ClustSat_HscpCorrWinv_secondXTalk4 = new TH1F("h_ClustSat_HscpCorrWinv_secondXTalk4","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondXTalk5 = new TH1F("h_ClustSat_HscpCorrWinv_secondXTalk5","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondXTalk6 = new TH1F("h_ClustSat_HscpCorrWinv_secondXTalk6","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondXTalk3 = new TH1F("h_ClustSat_HscpCorrWinv_secondXTalk3","",50,-2,2);
TH1F* h_ClustSat_HscpCorrWinv_secondXTalk2 = new TH1F("h_ClustSat_HscpCorrWinv_secondXTalk2","",50,-2,2);

//Etude Ih

TH2F* h2_IhHarmonic_sim = new TH2F("h2_IhHarmonic_sim","",100,0,3,50,0,50);
TH2F* h2_IhHarmonic_MyCorr = new TH2F("h2_IhHarmonic_MyCorr","",100,0,3,50,0,50);
TH2F* h2_IhHarmonic_HscpCorrWOinv = new TH2F("h2_IhHarmonic_HscpCorrWOinv","",100,0,3,50,0,50);
TH2F* h2_IhHarmonic_HscpCorr = new TH2F("h2_IhHarmonic_HscpCorr","",100,0,3,50,0,50);

//Etude ratio sat 

TH1F* h_ClustSat_MyCorrRatioSat0 = new TH1F("h_ClustSat_MyCorrRatioSat0","",50,-2,2);
TH1F* h_ClustSat_MyCorrRatioSat10 = new TH1F("h_ClustSat_MyCorrRatioSat10","",50,-2,2);
TH1F* h_ClustSat_MyCorrRatioSat20 = new TH1F("h_ClustSat_MyCorrRatioSat20","",50,-2,2);
TH1F* h_ClustSat_MyCorrRatioSat30 = new TH1F("h_ClustSat_MyCorrRatioSat30","",50,-2,2);
TH1F* h_ClustSat_MyCorrRatioSat40 = new TH1F("h_ClustSat_MyCorrRatioSat40","",50,-2,2);
TH1F* h_ClustSat_MyCorrRatioSat50 = new TH1F("h_ClustSat_MyCorrRatioSat50","",50,-2,2);
TH1F* h_ClustSat_MyCorrRatioSat60 = new TH1F("h_ClustSat_MyCorrRatioSat60","",50,-2,2);
TH1F* h_ClustSat_MyCorrRatioSat70 = new TH1F("h_ClustSat_MyCorrRatioSat70","",50,-2,2);
TH1F* h_ClustSat_MyCorrRatioSat80 = new TH1F("h_ClustSat_MyCorrRatioSat80","",50,-2,2);
TH1F* h_ClustSat_MyCorrRatioSat90 = new TH1F("h_ClustSat_MyCorrRatioSat90","",50,-2,2);
TH1F* h_ClustSat_MyCorrRatioSat100 = new TH1F("h_ClustSat_MyCorrRatioSat100","",50,-2,2);

//Etude correction+Ih

TH2F* h2_ClustSat_EsimEMyCorrRatioSat = new TH2F("h2_ClustSat_EsimEMyCorrRatioSat","",50,0,6000*pow(10,-6),50,0,6000*pow(10,-6));
TH2F* h2_ClustSat_EsimECorrHscp = new TH2F("h2_ClustSat_EsimECorrHscp","",50,0,6000*pow(10,-6),50,0,6000*pow(10,-6));
TProfile* prof_ClustSat_EsimEMyCorrRatioSat = new TProfile("prof_ClustSat_EsimEMyCorrRatioSat","",50,0,6000*pow(10,-6),"");
TProfile* prof_ClustSat_EsimECorrHscp = new TProfile("prof_ClustSat_EsimECorrHscp","",50,0,6000*pow(10,-6),"");

TH2F* h2_ClustSat_IhSimIhMyCorrRatioSat = new TH2F("h2_ClustSat_IhSimIhMyCorrRatioSat","",50,0,50,50,0,50);
TH2F* h2_ClustSat_IhSimIhHscpCorr = new TH2F("h2_ClustSat_IhSimIhHscpCorr","",50,0,50,50,0,50);
TProfile* prof_ClustSat_IhSimIhMyCorrRatioSat = new TProfile("prof_ClustSat_IhSimIhMyCorrRatioSat","",50,0,50,"");
TProfile* prof_ClustSat_IhSimIhHscpCorr = new TProfile("prof_ClustSat_IhSimIhHscpCorr","",50,0,50,"");

//Etude Sur-corr

//TH1F* h_ClustSat_Esim = new TH1F("h_ClustSat_Esim","",300,0,6000*pow(10,-6));
TH1F* h_ClustSat_Esim = new TH1F("h_ClustSat_Esim","",50,0,2500);
TH1F* h_ClustSat_Erec = new TH1F("h_ClustSat_Erec","",50,0,2500);
TH1I* h_ClustSat_ClusterSize = new TH1I("h_ClustSat_ClusterSize","",8,0,8);
TH1I* h_ClustSat_Sat254 = new TH1I("h_ClustSat_Sat254","",5,0,5);
TH1I* h_ClustSat_Sat255 = new TH1I("h_ClustSat_Sat255","",5,0,5);
TH1I* h_ClustSat_NSat = new TH1I("h_ClustSat_NSat","",5,0,5);
TH1D* h_ClustSat_layerlabel = new TH1D("h_ClustSat_layerlabel","",21,0,21);
TH1F* h_ClustSat_Eta = new TH1F("h_ClustSat_Eta","",10,0,2.5);
//TH1F* h_ClustSat_SurCorr_Esim = new TH1F("h_ClustSat_SurCorr_Esim","",300,0,6000*pow(10,-6));
TH1F* h_ClustSat_SurCorr_Esim = new TH1F("h_ClustSat_SurCorr_Esim","",50,0,2500);
TH1F* h_ClustSat_SurCorr_Erec = new TH1F("h_ClustSat_SurCorr_Erec","",50,0,2500);
TH1I* h_ClustSat_SurCorr_ClusterSize = new TH1I("h_ClustSat_SurCorr_ClusterSize","",8,0,8);
TH1I* h_ClustSat_SurCorr_Sat254 = new TH1I("h_ClustSat_SurCorr_Sat254","",5,0,5);
TH1I* h_ClustSat_SurCorr_Sat255 = new TH1I("h_ClustSat_SurCorr_Sat255","",5,0,5);
TH1I* h_ClustSat_SurCorr_NSat = new TH1I("h_ClustSat_SurCorr_NSat","",5,0,5);
TH1D* h_ClustSat_SurCorr_layerlabel = new TH1D("h_ClustSat_SurCorr_layerlabel","",21,0,21);
TH1F* h_ClustSat_SurCorr_Eta = new TH1F("h_ClustSat_SurCorr_Eta","",10,0,2.5);


TH2F* h2_ClustSatSurCorr_EsimEMyCorrRatioSat = new TH2F("h2_ClustSatSurCorr_EsimEMyCorrRatioSat","",50,0,6000*pow(10,-6),50,0,6000*pow(10,-6));
TH2F* h2_ClustSatSurCorr_EsimECorrHscp = new TH2F("h2_ClustSatSurCorr_EsimECorrHscp","",50,0,6000*pow(10,-6),50,0,6000*pow(10,-6));
TProfile* prof_ClustSatSurCorr_EsimEMyCorrRatioSat = new TProfile("prof_ClustSatSurCorr_EsimEMyCorrRatioSat","",50,0,6000*pow(10,-6),"");
TProfile* prof_ClustSatSurCorr_EsimECorrHscp = new TProfile("prof_ClustSatSurCorr_EsimECorrHscp","",50,0,6000*pow(10,-6),"");

TH2F* h2_ClustSatSurCorr_IhSimIhMyCorrRatioSat = new TH2F("h2_ClustSatSurCorr_IhSimIhMyCorrRatioSat","",50,0,50,50,0,50);
TH2F* h2_ClustSatSurCorr_IhSimIhHscpCorr = new TH2F("h2_ClustSatSurCorr_IhSimIhHscpCorr","",50,0,50,50,0,50);
TProfile* prof_ClustSatSurCorr_IhSimIhMyCorrRatioSat = new TProfile("prof_ClustSatSurCorr_IhSimIhMyCorrRatioSat","",50,0,50,"");
TProfile* prof_ClustSatSurCorr_IhSimIhHscpCorr = new TProfile("prof_ClustSatSurCorr_IhSimIhHscpCorr","",50,0,50,"");

TH2F* h2_ClustSatNoSurCorr_EsimEMyCorrRatioSat = new TH2F("h2_ClustSatNoSurCorr_EsimEMyCorrRatioSat","",50,0,6000*pow(10,-6),50,0,6000*pow(10,-6));
TH2F* h2_ClustSatNoSurCorr_EsimECorrHscp = new TH2F("h2_ClustSatNoSurCorr_EsimECorrHscp","",50,0,6000*pow(10,-6),50,0,6000*pow(10,-6));
TProfile* prof_ClustSatNoSurCorr_EsimEMyCorrRatioSat = new TProfile("prof_ClustSatNoSurCorr_EsimEMyCorrRatioSat","",50,0,6000*pow(10,-6),"");
TProfile* prof_ClustSatNoSurCorr_EsimECorrHscp = new TProfile("prof_ClustSatNoSurCorr_EsimECorrHscp","",50,0,6000*pow(10,-6),"");

TH2F* h2_ClustSatNoSurCorr_IhSimIhMyCorrRatioSat = new TH2F("h2_ClustSatNoSurCorr_IhSimIhMyCorrRatioSat","",50,0,50,50,0,50);
TH2F* h2_ClustSatNoSurCorr_IhSimIhHscpCorr = new TH2F("h2_ClustSatNoSurCorr_IhSimIhHscpCorr","",50,0,50,50,0,50);
TProfile* prof_ClustSatNoSurCorr_IhSimIhMyCorrRatioSat = new TProfile("prof_ClustSatNoSurCorr_IhSimIhMyCorrRatioSat","",50,0,50,"");
TProfile* prof_ClustSatNoSurCorr_IhSimIhHscpCorr = new TProfile("prof_ClustSatNoSurCorr_IhSimIhHscpCorr","",50,0,50,"");

TH2F* h2_ClustSatMoreThanOneSat_EsimEMyCorrRatioSat = new TH2F("h2_ClustSatMoreThanOneSat_EsimEMyCorrRatioSat","",50,0,6000*pow(10,-6),50,0,6000*pow(10,-6));
TH2F* h2_ClustSatMoreThanOneSat_EsimECorrHscp = new TH2F("h2_ClustSatMoreThanOneSat_EsimECorrHscp","",50,0,6000*pow(10,-6),50,0,6000*pow(10,-6));
TProfile* prof_ClustSatMoreThanOneSat_EsimEMyCorrRatioSat = new TProfile("prof_ClustSatMoreThanOneSat_EsimEMyCorrRatioSat","",50,0,6000*pow(10,-6),"");
TProfile* prof_ClustSatMoreThanOneSat_EsimECorrHscp = new TProfile("prof_ClustSatMoreThanOneSat_EsimECorrHscp","",50,0,6000*pow(10,-6),"");

TH2F* h2_ClustSatMoreThanOneSat_IhSimIhMyCorrRatioSat = new TH2F("h2_ClustSatMoreThanOneSat_IhSimIhMyCorrRatioSat","",50,0,50,50,0,50);
TH2F* h2_ClustSatMoreThanOneSat_IhSimIhHscpCorr = new TH2F("h2_ClustSatMoreThanOneSat_IhSimIhHscpCorr","",50,0,50,50,0,50);
TProfile* prof_ClustSatMoreThanOneSat_IhSimIhMyCorrRatioSat = new TProfile("prof_ClustSatMoreThanOneSat_IhSimIhMyCorrRatioSat","",50,0,50,"");
TProfile* prof_ClustSatMoreThanOneSat_IhSimIhHscpCorr = new TProfile("prof_ClustSatMoreThanOneSat_IhSimIhHscpCorr","",50,0,50,"");

TH2F* h2_ClustSatOneSat_EsimEMyCorrRatioSat = new TH2F("h2_ClustSatOneSat_EsimEMyCorrRatioSat","",50,0,6000*pow(10,-6),50,0,6000*pow(10,-6));
TH2F* h2_ClustSatOneSat_EsimECorrHscp = new TH2F("h2_ClustSatOneSat_EsimECorrHscp","",50,0,6000*pow(10,-6),50,0,6000*pow(10,-6));
TProfile* prof_ClustSatOneSat_EsimEMyCorrRatioSat = new TProfile("prof_ClustSatOneSat_EsimEMyCorrRatioSat","",50,0,6000*pow(10,-6),"");
TProfile* prof_ClustSatOneSat_EsimECorrHscp = new TProfile("prof_ClustSatOneSat_EsimECorrHscp","",50,0,6000*pow(10,-6),"");

TH2F* h2_ClustSatOneSat_IhSimIhMyCorrRatioSat = new TH2F("h2_ClustSatOneSat_IhSimIhMyCorrRatioSat","",50,0,50,50,0,50);
TH2F* h2_ClustSatOneSat_IhSimIhHscpCorr = new TH2F("h2_ClustSatOneSat_IhSimIhHscpCorr","",50,0,50,50,0,50);
TProfile* prof_ClustSatOneSat_IhSimIhMyCorrRatioSat = new TProfile("prof_ClustSatOneSat_IhSimIhMyCorrRatioSat","",50,0,50,"");
TProfile* prof_ClustSatOneSat_IhSimIhHscpCorr = new TProfile("prof_ClustSatOneSat_IhSimIhHscpCorr","",50,0,50,"");

TH2F* h2_ClustSatPassingHscpTest_EsimEMyCorrRatioSat = new TH2F("h2_ClustSatPassingHscpTest_EsimEMyCorrRatioSat","",50,0,6000*pow(10,-6),50,0,6000*pow(10,-6));
TH2F* h2_ClustSatPassingHscpTest_EsimECorrHscp = new TH2F("h2_ClustSatPassingHscpTest_EsimECorrHscp","",50,0,6000*pow(10,-6),50,0,6000*pow(10,-6));
TProfile* prof_ClustSatPassingHscpTest_EsimEMyCorrRatioSat = new TProfile("prof_ClustSatPassingHscpTest_EsimEMyCorrRatioSat","",50,0,6000*pow(10,-6),"");
TProfile* prof_ClustSatPassingHscpTest_EsimECorrHscp = new TProfile("prof_ClustSatPassingHscpTest_EsimECorrHscp","",50,0,6000*pow(10,-6),"");

TH2F* h2_ClustSatPassingHscpTest_IhSimIhMyCorrRatioSat = new TH2F("h2_ClustSatPassingHscpTest_IhSimIhMyCorrRatioSat","",50,0,50,50,0,50);
TH2F* h2_ClustSatPassingHscpTest_IhSimIhHscpCorr = new TH2F("h2_ClustSatPassingHscpTest_IhSimIhHscpCorr","",50,0,50,50,0,50);
TProfile* prof_ClustSatPassingHscpTest_IhSimIhMyCorrRatioSat = new TProfile("prof_ClustSatPassingHscpTest_IhSimIhMyCorrRatioSat","",50,0,50,"");
TProfile* prof_ClustSatPassingHscpTest_IhSimIhHscpCorr = new TProfile("prof_ClustSatPassingHscpTest_IhSimIhHscpCorr","",50,0,50,"");

//Criteria studies

TH1F* h1_ClustSat_testThresholdSat_HscpCorr = new TH1F("h1_ClustSat_testThresholdSat_HscpCorr","",200,-2,2);
TH1F* h1_ClustSat_testThresholdDiff_HscpCorr = new TH1F("h1_ClustSat_testThresholdDiff_HscpCorr","",200,-2,2);
TH1F* h1_ClustSat_testBothThreshold_HscpCorr = new TH1F("h1_ClustSat_testBothThreshold_HscpCorr","",200,-2,2);


//Fake clust

TH1F* h1_FakeClusters = new TH1F("h1_FakeClusters","",50,0,1);
TH1F* h1_FakeClusters_testThresholdSat = new TH1F("h1_FakeClusters_testThresholdSat","",50,0,1);
TH1F* h1_FakeClusters_testThresholdDiff = new TH1F("h1_FakeClusters_testThresholdDiff","",50,0,1);

TH1F* h1_DeltaQ_True = new TH1F("h1_DeltaQ_True","",100,0,100);
TH1F* h1_DeltaQ_False = new TH1F("h1_DeltaQ_False","",100,0,100);

//DeltaQ Criterion study

TH1F* h1_ClustSat_DeltaQ = new TH1F("h1_ClustSat_DeltaQ","",200,-2,2);
TH1F* h1_ClustSat_DeltaQ35 = new TH1F("h1_ClustSat_DeltaQ35","",200,-2,2);
TH1F* h1_ClustSat_DeltaQ30 = new TH1F("h1_ClustSat_DeltaQ30","",200,-2,2);
TH1F* h1_ClustSat_DeltaQ25 = new TH1F("h1_ClustSat_DeltaQ25","",200,-2,2);
TH1F* h1_ClustSat_DeltaQ20 = new TH1F("h1_ClustSat_DeltaQ20","",200,-2,2);
TH1F* h1_ClustSat_DeltaQ15 = new TH1F("h1_ClustSat_DeltaQ15","",200,-2,2);

TH1F* h1_ClustSatDeltaQ_nsat2_TrueClusters = new TH1F("h1_ClustSatDeltaQ_nsat2_TrueClusters","",100,0,100);
TH1F* h1_ClustSatDeltaQ_nsat2_FakeClusters = new TH1F("h1_ClustSatDeltaQ_nsat2_FakeClusters","",100,0,100);


TH2F* h2_ClustSat_DeltaQ = new TH2F("h2_ClustSat_DeltaQ","",40,-2,2,250,0,250);

//ThresholdSat crit. study

TH2F* h2_ClustSat_ThresholdSatCrit = new TH2F("h2_ClustSat_ThresholdSatCrit","",40,-2,2,150,-10,140);


//Study crit. MyCorr

TH2F* h2_ClustSat_QsimBetaGamma_Gluino2400 = new TH2F("h2_ClustSat_QsimBetaGamma_Gluino2400","",100,0,2,100,0,4000);
TProfile* prof_ClustSat_QsimBetaGamma_Gluino2400 = new TProfile("prof_ClustSat_QsimBetaGamma_Gluino2400","",100,0,2,0,4000,"");

TH2F* h2_ClustSat_QsimBetaGamma_Gluino1400 = new TH2F("h2_ClustSat_QsimBetaGamma_Gluino1400","",100,0,2,100,0,4000);
TProfile* prof_ClustSat_QsimBetaGamma_Gluino1400 = new TProfile("prof_ClustSat_QsimBetaGamma_Gluino1400","",100,0,2,0,4000,"");

TH2F* h2_ClustSat_QsimBetaGamma_Stau651 = new TH2F("h2_ClustSat_QsimBetaGamma_Stau651","",100,0,2,100,0,4000);
TProfile* prof_ClustSat_QsimBetaGamma_Stau651 = new TProfile("prof_ClustSat_QsimBetaGamma_Stau651","",100,0,2,0,4000,"");



TH1F* h1_ClustSat_Qsim = new TH1F("h1_ClustSat_Qsim","",100,0,4000);

TProfile* profile_ClustSat_HscpNoInv = new TProfile("profile_ClustSat_HscpNoInv","",100,0,4000);

TH1F* h1_ClustSat_MyCorrModulGeom_WihtoutCrit = new TH1F("h1_ClustSat_MyCorrModulGeom_WihtoutCrit","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_WihtoutCrit = new TProfile("profile_ClustSat_MyCorrModulGeom_WihtoutCrit","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_DeltaQ = new TH1F("h1_ClustSat_MyCorrModulGeom_DeltaQ","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_DeltaQ = new TProfile("profile_ClustSat_MyCorrModulGeom_DeltaQ","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_DeltaQ40 = new TH1F("h1_ClustSat_MyCorrModulGeom_DeltaQ40","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_DeltaQ40 = new TProfile("profile_ClustSat_MyCorrModulGeom_DeltaQ40","",100,0,4000,"");
TH2F* h2_ClustSat_MyCorrModulGeom_DeltaQ40 = new TH2F("h2_ClustSat_MyCorrModulGeom_DeltaQ40","",200,0,4000,400,0,4000);

TH1F* h1_ClustSat_MyCorrModulGeom_DeltaQ30 = new TH1F("h1_ClustSat_MyCorrModulGeom_DeltaQ30","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_DeltaQ30 = new TProfile("profile_ClustSat_MyCorrModulGeom_DeltaQ30","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_DeltaQ20 = new TH1F("h1_ClustSat_MyCorrModulGeom_DeltaQ20","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_DeltaQ20 = new TProfile("profile_ClustSat_MyCorrModulGeom_DeltaQ20","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_DeltaQ10 = new TH1F("h1_ClustSat_MyCorrModulGeom_DeltaQ10","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_DeltaQ10 = new TProfile("profile_ClustSat_MyCorrModulGeom_DeltaQ10","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_DeltaQ15 = new TH1F("h1_ClustSat_MyCorrModulGeom_DeltaQ15","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_DeltaQ15 = new TProfile("profile_ClustSat_MyCorrModulGeom_DeltaQ15","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_DeltaQ5 = new TH1F("h1_ClustSat_MyCorrModulGeom_DeltaQ5","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_DeltaQ5 = new TProfile("profile_ClustSat_MyCorrModulGeom_DeltaQ5","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_DeltaQ50 = new TH1F("h1_ClustSat_MyCorrModulGeom_DeltaQ50","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_DeltaQ50 = new TProfile("profile_ClustSat_MyCorrModulGeom_DeltaQ50","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_DeltaQ60 = new TH1F("h1_ClustSat_MyCorrModulGeom_DeltaQ60","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_DeltaQ60 = new TProfile("profile_ClustSat_MyCorrModulGeom_DeltaQ60","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ThresholdSat25 = new TH1F("h1_ClustSat_MyCorrModulGeom_ThresholdSat25","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ThresholdSat25 = new TProfile("profile_ClustSat_MyCorrModulGeom_ThresholdSat25","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_RatioSat90 = new TH1F("h1_ClustSat_MyCorrModulGeom_RatioSat90","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_RatioSat90 = new TProfile("profile_ClustSat_MyCorrModulGeom_RatioSat90","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_RatioSat60 = new TH1F("h1_ClustSat_MyCorrModulGeom_RatioSat60","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_RatioSat60 = new TProfile("profile_ClustSat_MyCorrModulGeom_RatioSat60","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat90 = new TH1F("h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat90","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_DeltaQRatioSat90 = new TProfile("profile_ClustSat_MyCorrModulGeom_DeltaQRatioSat90","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat90 = new TH1F("h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat90","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat90 = new TProfile("profile_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat90","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat60 = new TH1F("h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat60","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_DeltaQRatioSat60 = new TProfile("profile_ClustSat_MyCorrModulGeom_DeltaQRatioSat60","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat60 = new TH1F("h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat60","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat60 = new TProfile("profile_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat60","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQ = new TH1F("h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQ","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQ = new TProfile("profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQ","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60 = new TH1F("h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60 = new TProfile("profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat90 = new TH1F("h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat90","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat90 = new TProfile("profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat90","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ChargeMinEstimated = new TH1F("h1_ClustSat_MyCorrModulGeom_ChargeMinEstimated","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ChargeMinEstimated = new TProfile("profile_ClustSat_MyCorrModulGeom_ChargeMinEstimated","",100,0,4000,"");

TH2F* h2_ClustSat_Overlap = new TH2F("h2_ClustSat_Overlap","",100,0,100,100,0,1);
TProfile* profile_ClustSat_Overlap = new TProfile("profile_ClustSat_Overlap","",100,0,100,0,1,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEta = new TH1F("h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEta","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEta = new TProfile("profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEta","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TIB1 = new TH1F("h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TIB1","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TIB1 = new TProfile("profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TIB1","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_OB1 = new TH1F("h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_OB1","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_OB1 = new TProfile("profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_OB1","",100,0,4000,"");

TH1F* h1_ClustSat_HscpNoInv_CutEta = new TH1F("h1_ClustSat_HscpNoInv_CutEta","",200,-2,2);
TProfile* profile_ClustSat_HscpNoInv_CutEta = new TProfile("profile_ClustSat_HscpNoInv_CutEta","",100,0,4000,"");

TH1F* h1_ClustSat_HscpNoInv_TIB1 = new TH1F("h1_ClustSat_HscpNoInv_TIB1","",200,-2,2);
TProfile* profile_ClustSat_HscpNoInv_TIB1 = new TProfile("profile_ClustSat_HscpNoInv_TIB1","",100,0,4000,"");

TH1F* h1_ClustSat_HscpNoInv_OB1 = new TH1F("h1_ClustSat_HscpNoInv_OB1","",200,-2,2);
TProfile* profile_ClustSat_HscpNoInv_OB1 = new TProfile("profile_ClustSat_HscpNoInv_OB1","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEtaSup = new TH1F("h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEtaSup","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEtaSup = new TProfile("profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEtaSup","",100,0,4000,"");

TH1F* h1_ClustSat_HscpNoInv_CutEtaSup = new TH1F("h1_ClustSat_HscpNoInv_CutEtaSup","",200,-2,2);
TProfile* profile_ClustSat_HscpNoInv_CutEtaSup = new TProfile("profile_ClustSat_HscpNoInv_CutEtaSup","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TEC = new TH1F("h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TEC","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TEC = new TProfile("profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TEC","",100,0,4000,"");

TH1F* h1_ClustSat_HscpNoInv_TEC = new TH1F("h1_ClustSat_HscpNoInv_TEC","",200,-2,2);
TProfile* profile_ClustSat_HscpNoInv_TEC = new TProfile("profile_ClustSat_HscpNoInv_TEC","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TID = new TH1F("h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TID","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TID = new TProfile("profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TID","",100,0,4000,"");

TH1F* h1_ClustSat_HscpNoInv_TID = new TH1F("h1_ClustSat_HscpNoInv_TID","",200,-2,2);
TProfile* profile_ClustSat_HscpNoInv_TID = new TProfile("profile_ClustSat_HscpNoInv_TID","",100,0,4000,"");

TH1F* h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_Barrel = new TH1F("h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_Barrel","",200,-2,2);
TProfile* profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_Barrel = new TProfile("profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_Barrel","",100,0,4000,"");

TH1F* h1_ClustSat_HscpNoInv_Barrel = new TH1F("h1_ClustSat_HscpNoInv_Barrel","",200,-2,2);
TProfile* profile_ClustSat_HscpNoInv_Barrel = new TProfile("profile_ClustSat_HscpNoInv_Barrel","",100,0,4000,"");

TH1F* h1_Proposal = new TH1F("h1_Proposal","",200,-2,2);
TProfile* profile_Proposal = new TProfile("profile_Proposal","",100,0,4000,"");
TH1F* h1_ratioProfile = new TH1F("h1_ratioProfile","",100,0,4000);


TH1F* h1_IhWTEC = new TH1F("h1_IhWTEC","",100,0,15);
TH1F* h1_IhWOTEC = new TH1F("h1_IhWOTEC","",100,0,15);

TH2F* distribdeltaQ = new TH2F("distribdeltaQ","",100,0,100,100,-2,2);
TH2F* distribminQ = new TH2F("distribminQ","",100,0,100,100,-2,2); 

TH1F* distribdedx = new TH1F("distribdedx","",100,0,40);
TH1F* distribdedx2 = new TH1F("distribdedx2","",100,0,40);

TH1F* h1_NewProcedure = new TH1F("h1_NewProcedure","",200,-2,2);
TProfile* profile_NewProcedure = new TProfile("profile_NewProcedure","",100,0,4000,"");

TH2F* IhPoverM = new TH2F("IhPoverM","",100,0,5,100,0,50);

TH1F* h1_IhSim = new TH1F("h1_IhSim","",200,-2,2);
TH1F* h1_IhHscpCorr = new TH1F("h1_IhHscpCorr","",200,-2,2);
TH1F* h1_IhProposal = new TH1F("h1_IhProposal","",200,-2,2);
TH1F* h1_IhScdProposal = new TH1F("h1_IhScdProposal","",200,-2,2);
TH1F* h1_IhThirdProposal = new TH1F("h1_IhThirdProposal","",200,-2,2);
TH1F* h1_IhFourthProposal = new TH1F("h1_IhFourthProposal","",200,-2,2);
TH1F* h1_IhRecWOSat = new TH1F("h1_IhRecWOSat","",200,-2,2);

TH2F* h2_IhSim = new TH2F("h2_IhSim","",100,0,20,100,0,20);
TH2F* h2_IhHscpCorr = new TH2F("h2_IhHscpCorr","",100,0,20,100,0,20);
TH2F* h2_IhProposal = new TH2F("h2_IhProposal","",100,0,20,100,0,20);
TH2F* h2_IhScdProposal = new TH2F("h2_IhScdProposal","",100,0,20,100,0,20);
TH2F* h2_IhThirdProposal = new TH2F("h2_IhThirdProposal","",100,0,20,100,0,20);
TH2F* h2_IhFourthProposal = new TH2F("h2_IhFourthProposal","",100,0,20,100,0,20);
TH2F* h2_IhRecWOSat = new TH2F("h2_IhRecWOSat","",100,0,20,100,0,20);

TH2F* h2_Resolution_IhSim = new TH2F("h2_Resolution_IhSim","",100,0,20,100,-2,2);
TH2F* h2_Resolution_IhHscpCorr = new TH2F("h2_Resolution_IhHscpCorr","",100,0,20,100,-2,2);
TH2F* h2_Resolution_IhProposal = new TH2F("h2_Resolution_IhProposal","",100,0,20,100,-2,2);
TH2F* h2_Resolution_IhScdProposal = new TH2F("h2_Resolution_IhScdProposal","",100,0,20,100,-2,2);
TH2F* h2_Resolution_IhThirdProposal = new TH2F("h2_Resolution_IhThirdProposal","",100,0,20,100,-2,2);
TH2F* h2_Resolution_IhFourthProposal = new TH2F("h2_Resolution_IhFourthProposal","",100,0,20,100,-2,2);
TH2F* h2_Resolution_IhRecWOSat = new TH2F("h2_Resolution_IhRecWOSat","",100,0,20,100,-2,2);


TH1F* distribIhSim = new TH1F("distribIhSim","",50,0,20);

//COUPURES SUR LES EVENTS


b1->SetThresholdPartId(0.7); 
b1->SetThresholdPt(65); 
b1->SetThresholdP(0);
b1->SetThresholdEta(5);



//Ratio nombre de pistes satures 

int ncluster=0;
int sat1=0;
int sat2=0;
int sat3=0;
int sat1HscpTest=0;
int nSat=0;
int nHscpTest=0;
int nRatSatTest=0;
int nHscpAndRatSatTest=0;
int sat2RatSatTest=0;
int surcorr=0;
int surcorrHSCP=0;


int nParticle=0;
int nPatriclePassingMyCrit=0;


int counterpass=0;




int countThresholdSat=0;
int countThresholdDiff=0;
int countBothThreshold=0;
int countSatCluster=0;



//Efficiency criteria

float nPassingNoCrit=0.;
float nPassingHSCPtest=0.;
float nPassingDeltaQ=0.;
float nPassingDeltaQ40=0.;
float nPassingDeltaQ30=0.;
float nPassingDeltaQ20=0.;
float nPassingDeltaQ10=0.;
float nPassingDeltaQ50=0.;
float nPassingDeltaQ60=0.;
float nPassingDeltaQ15=0.;
float nPassingDeltaQ5=0.;
float nPassingThresholdSat25=0.;
float nPassingRatioSat60=0.;
float nPassingRatioSat90=0.;
float nPassingDeltaQRatioSat60=0.;
float nPassingDeltaQRatioSat90=0.;
float nPassingThresholdSat25RatioSat60=0.;
float nPassingThresholdSat25RatioSat90=0.;
float nPassingThresholdSat25DeltaQ=0.;
float nPassingThresholdSat25DeltaQRatioSat60=0.;
float nPassingThresholdSat25DeltaQRatioSat90=0.;
float nPassingChargeMinEstimated=0.;
float nPassingNewProcedure=0.;

vector<string> nopass;

//BOUCLE SUR LES EVENTS

int entries = nentries;

for(int i=0;i<entries;i++)
{
    if(i%1000==0) cout<<"Event "<<i<<endl;
    b1->GetEntry(i);

	//BOUCLE SUR LES TRACES

    for(int track=0;track<b1->GetNtracks();track++)
    {
        vector<int> vect_partID;
        vector<float> vect_Eloss;
        vector<float> vect_charge;
        vector<float> vect_dedx;
        vector<float> vect_dqdx;
        vector<float> vect_dqcorrdx;
        vector<float> vect_dqcorrdx2;
        vector<float> vect_plength;

        vector<float> vect_decorrdxnewSatLayer;
        vector<float> vect_decorrdx;

        vector<float> dEsimdx;
        vector<float> dErecdx;
        vector<float> dEcorrdx;
        vector<float> dEfulldx_Layer;
        vector<float> dEfulldx_ModulGeom;

        vector<float> dEdxWTEC;
        vector<float> dEdxWOTEC;



        vector<float> VectDecorrDxLayer;
        vector<float> VectDecorrDxModulGeom;
        float pt                = b1->GetVectTrack()[track].GetPt();
        float p                 = b1->GetVectTrack()[track].GetP();
        float eta               = b1->GetVectTrack()[track].GetEta();
        float phi               = b1->GetVectTrack()[track].GetPhi();
        float ias_ampl          = b1->GetVectTrack()[track].GetIasAmpl();
        int NCluster            = b1->GetVectTrack()[track].GetNCluster();
        int NClustSat254        = b1->GetVectTrack()[track].GetNSatCluster(254);
        int NClustSat255        = b1->GetVectTrack()[track].GetNSatCluster(255);
        float RatioNClusterSat254 = (double)NClustSat254/(double)NCluster;
        float RatioNClusterSat255 = (double)NClustSat255/(double)NCluster;
        int id                  = b1->GetVectTrack()[track].GetPartId();
        float PoverM            = GetPoverM(p,id);


        nParticle++;
        if(RatioNClusterSat254+RatioNClusterSat255>=ratsat) nPatriclePassingMyCrit++;


        bool testSurCorr=false;
        bool testMoreThanOne=false;
        bool testHscptrack=false;


        std::vector<float> dEdxsim;
        std::vector<float> dEdxMyCorr;
        std::vector<float> dEdxHscpCorrWOinv;
        std::vector<float> dEdxHscpCorr;


        std::vector<float> DeDxSim;
        std::vector<float> DeDxHscpCorr;
        std::vector<float> DeDxProposal;
        std::vector<float> DeDxScdProposal;
        std::vector<float> DeDxThirdProposal;
        std::vector<float> DeDxFourthProposal;
        std::vector<float> DeDxRecWOSat;


    if(idlow<=abs(id) && abs(id)<=idup)

    {
		//BOUCLE SUR LES CLUSTERS

        for(int cluster=0;cluster<b1->GetVectTrack()[track].GetNCluster();cluster++)
        {
            int fakecluster=0;

            float charge        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSclusCharge();
            float charge_corr   = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSclusChargeCorr();
            float EcorrBef      = charge_corr*(3.61*pow(10,-9)*247);
			float Erec			= charge*(3.61*pow(10,-9)*247);
			float chargeNoSat	= b1->GetVectTrack()[track].GetVectClusters()[cluster].GetChargeWithoutSaturation();
			float ErecNoSat		= chargeNoSat*(3.61*pow(10,-9)*247);
            float Eloss         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetEloss();
            float pathlength    = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetPathLength();
            int layer           = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetLayer();
            int layerLabel      = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetLayerLabel();
            bool sat254         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSat254();
            bool sat255         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSat255();
            bool shape          = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetShape();
            int nsimhits        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSimHits();
            int nstrips         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNStrip();
            int nsat254         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStrip(254);
            int nsat255         = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSatStrip(255);
            bool edge           = b1->GetVectTrack()[track].GetVectClusters()[cluster].Edge();
            bool cut            = b1->GetVectTrack()[track].GetVectClusters()[cluster].Cut();
            int detid           = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetDetId();
            int subdetid        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetSubDetId();
            int modulgeom		= b1->GetVectTrack()[track].GetVectClusters()[cluster].GetModulGeom();
            int maxstrip        = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetMaxStrip();
            int cluster_id      = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetPartId();
			float EcorrSatLayer	= Erec;
            float EcorrSatModul = Erec;
			float EcorrNoSat	= ErecNoSat;
            int charge_sim      = Eloss/(3.61*pow(10,-9)*247);
            DeDxSim.push_back(Eloss/pathlength);

            bool fakeclust = true;

            for(int simhit=0;simhit<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSimHits();simhit++)
            {
                int simhit_id = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectSimHits()[simhit].GetPartId();
                //if(simhit_id!=cluster_id) cout<<cluster<<"    "<<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNSimHits()<<"    "<<simhit<<"    "<<simhit_id<<"    "<<cluster_id<<endl;
                if(simhit_id==id) fakeclust=false;
            }

            float Ratio_fakeclust = (float)fakecluster/(float)NCluster;
            h1_FakeClusters->Fill(Ratio_fakeclust);

            bool testThresholdSat=false;
            bool testThresholdDiff=false;
            bool testBothThreshold=false;


			h2ElossErec->Fill(Eloss,Erec);

			hmodulgeom->Fill(modulgeom);
			hlayerlabel->Fill(layerLabel);
			hmodulgeomvslayer->Fill(layerLabel-1,modulgeom-1);

            vect_dedx.push_back(Eloss/pathlength);

            //BOUCLE SUR LES STRIPS D'UN CLUSTER

            int ChargeCluster=0;
            vector<int> VectChargeStrip;

            int CentralChargeAndNeighbours=0;



            float Qminus=0.;
            float Qplus=0.;
            bool testQminus=true;

            float StripNoSat=0;
            float StripSat=0;
            

            for(int strip=0;strip<b1->GetVectTrack()[track].GetVectClusters()[cluster].GetNStrip();strip++)
            {
                float ChargeStrip = b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip].GetAmpl();
                VectChargeStrip.push_back(b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip].GetAmpl());
                ChargeCluster+=b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip].GetAmpl();

                if(ChargeStrip>253 && testQminus)
                {
                    testQminus=false;
                    Qminus=b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip-1].GetAmpl();
                }
                if(ChargeStrip>253)
                {
                    Qplus=b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip+1].GetAmpl();
                }
                if(ChargeStrip<254) StripNoSat+=ChargeStrip;
                if(ChargeStrip>253) StripSat+=ChargeStrip;

                if(ChargeStrip>253)
                {
                    CentralChargeAndNeighbours+=ChargeStrip;
                    CentralChargeAndNeighbours+=b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip-1].GetAmpl();
                    if(b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip+1].GetAmpl()>253)
                    {
                        CentralChargeAndNeighbours+=b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip+1].GetAmpl();
                        if(b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip+2].GetAmpl()>253)
                        {
                            CentralChargeAndNeighbours+=b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip+2].GetAmpl();
                            CentralChargeAndNeighbours+=b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip+3].GetAmpl();
                        }
                        else CentralChargeAndNeighbours+=b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip+2].GetAmpl();
                    }
                    else CentralChargeAndNeighbours+=b1->GetVectTrack()[track].GetVectClusters()[cluster].GetVectStrips()[strip+1].GetAmpl();
                }

            }

            bool testCorrHscp = false;
            vector<int>::const_iterator mQ = max_element(VectChargeStrip.begin(), VectChargeStrip.end());
            if(*mQ>253 && *(mQ-1)>25 && *(mQ+1)>25 && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40 ) testCorrHscp=true;
        
            float DeltaQ = abs(*(mQ-1) - *(mQ+1));
            float NewTestEnergy = CentralChargeAndNeighbours*(3.61*pow(10,-9)*247);
            float DeltaQCalc = abs(Qplus-Qminus);

            int minQ=9999;

            if(*(mQ-1)<minQ) minQ=(*(mQ-1));
            if(*(mQ+1)<minQ) minQ=(*(mQ+1));

            bool testProposal=false;
            if(Qminus>=25 && Qplus>=25 && RatioNClusterSat254+RatioNClusterSat255>=0.6 && DeltaQCalc<=40) testProposal=true;

            // if(!testProposal) continue;

            h2_ClustSat_Overlap->Fill(minQ,RatioNClusterSat254+RatioNClusterSat255);
            profile_ClustSat_Overlap->Fill(minQ,RatioNClusterSat254+RatioNClusterSat255);



            float ChargeMinEstimated = StripNoSat+1.2*StripSat;


            //--------------------------------

            if(XTalkInvStudy)
            {
                

                std::vector<int> VectChargeCorrHSCP_WithoutInversion = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,1,RatioNClusterSat254,RatioNClusterSat255,0);
                std::vector<int> VectChargeCorrHSCP = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                std::vector<int> VectChargeMyCorr_WithoutInversion = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,0);
                std::vector<int> VectChargeMyCorr = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,4,RatioNClusterSat254,RatioNClusterSat255,0);
                std::vector<int> VectChargeMyCorr_WithoutInversion_RatioSat = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,ratsat);
                std::vector<int> VectChargeMyCorr_RatioSat = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,4,RatioNClusterSat254,RatioNClusterSat255,ratsat);
                
                std::vector<int> VectChargeMyCorrModulGeom_WithoutInversion_RatioSat = CrossTalkInv(EcorrModulGeom,VectChargeStrip,0.10,0.04,20,25,modulgeom,3,RatioNClusterSat254,RatioNClusterSat255,ratsat);

                float ChargeCorrHSCP_WithoutInversion=0.;
                float ChargeCorrHSCP=0.;
                float ChargeMyCorr_WithoutInversion=0.;
                float ChargeMyCorr=0.;
                float ChargeMyCorr_WithoutInversion_RatioSat=0.;
                float ChargeMyCorr_RatioSat=0.;


                float ChargeMyCorrModulGeom_WithoutInversion_RatioSat=0.;

                for(int i=0;i<VectChargeCorrHSCP_WithoutInversion.size();i++) ChargeCorrHSCP_WithoutInversion+=VectChargeCorrHSCP_WithoutInversion.at(i);
                for(int i=0;i<VectChargeCorrHSCP.size();i++) ChargeCorrHSCP+=VectChargeCorrHSCP.at(i);
                for(int i=0;i<VectChargeMyCorr_WithoutInversion.size();i++) ChargeMyCorr_WithoutInversion+=VectChargeMyCorr_WithoutInversion.at(i);
                for(int i=0;i<VectChargeMyCorr.size();i++) ChargeMyCorr+=VectChargeMyCorr.at(i);
                for(int i=0;i<VectChargeMyCorr_WithoutInversion_RatioSat.size();i++) ChargeMyCorr_WithoutInversion_RatioSat+=VectChargeMyCorr_WithoutInversion_RatioSat.at(i);
                for(int i=0;i<VectChargeMyCorr_RatioSat.size();i++) ChargeMyCorr_RatioSat+=VectChargeMyCorr_RatioSat.at(i);

                for(int i=0;i<VectChargeMyCorrModulGeom_WithoutInversion_RatioSat.size();i++) ChargeMyCorrModulGeom_WithoutInversion_RatioSat+=VectChargeMyCorrModulGeom_WithoutInversion_RatioSat.at(i);

                float ECorrHSCP_WithoutInversion=(3.61*pow(10,-9)*247)*ChargeCorrHSCP_WithoutInversion;
                float ECorrHSCP=(3.61*pow(10,-9)*247)*ChargeCorrHSCP;
                //float EMyCorr_WithoutInversion=(3.61*pow(10,-9)*247)*ChargeMyCorr_WithoutInversion;
                float EMyCorr=(3.61*pow(10,-9)*247)*ChargeMyCorr;
                
                float EMyCorr_WithoutInversion=(3.61*pow(10,-9)*247)*ChargeMyCorr_WithoutInversion_RatioSat;
                
                float EMyCorr_WithoutInversion_RatioSat=(3.61*pow(10,-9)*247)*ChargeMyCorr_WithoutInversion_RatioSat;
                float EMyCorr_RatioSat=(3.61*pow(10,-9)*247)*ChargeMyCorr_RatioSat;

                float EMyCorrModulGeom_WithoutInversion_RatioSat=(3.61*pow(10,-9)*247)*ChargeMyCorrModulGeom_WithoutInversion_RatioSat;


                hWithoutCorr->Fill((Erec-Eloss)/Eloss);
                hCorrHSCP_WithoutInversion->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                hMyCorr_WithoutInversion->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                hMyCorr_WithoutInversion_RatioSat->Fill((EMyCorr_WithoutInversion_RatioSat-Eloss)/Eloss);
 
                h_DistribP->Fill(p);
                h_DistribEta->Fill(abs(eta));
                h_DistribSize->Fill(nstrips);
                h_DistribNSat->Fill(nsat254+nsat255);


            //First Threshold study

                std::vector<int> VectChargeCorrHSCPfirstThreshold20 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstThreshold20=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstThreshold20.size();i++) ChargeCorrHSCPfirstThreshold20+=VectChargeCorrHSCPfirstThreshold20.at(i);
                float ECorrHSCPfirstThreshold20=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstThreshold20;

                std::vector<int> VectChargeCorrHSCPfirstThreshold50 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,50,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstThreshold50=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstThreshold50.size();i++) ChargeCorrHSCPfirstThreshold50+=VectChargeCorrHSCPfirstThreshold50.at(i);
                float ECorrHSCPfirstThreshold50=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstThreshold50;

                std::vector<int> VectChargeCorrHSCPfirstThreshold100 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,100,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstThreshold100=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstThreshold100.size();i++) ChargeCorrHSCPfirstThreshold100+=VectChargeCorrHSCPfirstThreshold100.at(i);
                float ECorrHSCPfirstThreshold100=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstThreshold100;

                std::vector<int> VectChargeCorrHSCPfirstThreshold150 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,150,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstThreshold150=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstThreshold150.size();i++) ChargeCorrHSCPfirstThreshold150+=VectChargeCorrHSCPfirstThreshold150.at(i);
                float ECorrHSCPfirstThreshold150=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstThreshold150;

                std::vector<int> VectChargeCorrHSCPfirstThreshold200 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,200,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstThreshold200=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstThreshold200.size();i++) ChargeCorrHSCPfirstThreshold200+=VectChargeCorrHSCPfirstThreshold200.at(i);
                float ECorrHSCPfirstThreshold200=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstThreshold200;

            //Second Threshold Study

                std::vector<int> VectChargeCorrHSCPsecondThreshold25 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondThreshold25=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondThreshold25.size();i++) ChargeCorrHSCPsecondThreshold25+=VectChargeCorrHSCPsecondThreshold25.at(i);
                float ECorrHSCPsecondThreshold25=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondThreshold25;

                std::vector<int> VectChargeCorrHSCPsecondThreshold50 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,50,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondThreshold50=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondThreshold50.size();i++) ChargeCorrHSCPsecondThreshold50+=VectChargeCorrHSCPsecondThreshold50.at(i);
                float ECorrHSCPsecondThreshold50=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondThreshold50;

                std::vector<int> VectChargeCorrHSCPsecondThreshold100 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,100,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondThreshold100=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondThreshold100.size();i++) ChargeCorrHSCPsecondThreshold100+=VectChargeCorrHSCPsecondThreshold100.at(i);
                float ECorrHSCPsecondThreshold100=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondThreshold100;

                std::vector<int> VectChargeCorrHSCPsecondThreshold150 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,150,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondThreshold150=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondThreshold150.size();i++) ChargeCorrHSCPsecondThreshold150+=VectChargeCorrHSCPsecondThreshold150.at(i);
                float ECorrHSCPsecondThreshold150=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondThreshold150;

                std::vector<int> VectChargeCorrHSCPsecondThreshold200 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,200,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondThreshold200=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondThreshold200.size();i++) ChargeCorrHSCPsecondThreshold200+=VectChargeCorrHSCPsecondThreshold200.at(i);
                float ECorrHSCPsecondThreshold200=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondThreshold200;

            //First XTalk study 

                std::vector<int> VectChargeCorrHSCPfirstXTalk10 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstXTalk10=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstXTalk10.size();i++) ChargeCorrHSCPfirstXTalk10+=VectChargeCorrHSCPfirstXTalk10.at(i);
                float ECorrHSCPfirstXTalk10=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstXTalk10;

                std::vector<int> VectChargeCorrHSCPfirstXTalk11 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.11,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstXTalk11=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstXTalk11.size();i++) ChargeCorrHSCPfirstXTalk11+=VectChargeCorrHSCPfirstXTalk11.at(i);
                float ECorrHSCPfirstXTalk11=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstXTalk11;

                std::vector<int> VectChargeCorrHSCPfirstXTalk12 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.12,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstXTalk12=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstXTalk12.size();i++) ChargeCorrHSCPfirstXTalk12+=VectChargeCorrHSCPfirstXTalk12.at(i);
                float ECorrHSCPfirstXTalk12=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstXTalk12;

                std::vector<int> VectChargeCorrHSCPfirstXTalk9 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.09,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstXTalk9=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstXTalk9.size();i++) ChargeCorrHSCPfirstXTalk9+=VectChargeCorrHSCPfirstXTalk9.at(i);
                float ECorrHSCPfirstXTalk9=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstXTalk9;

                std::vector<int> VectChargeCorrHSCPfirstXTalk8 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.08,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPfirstXTalk8=0.;
                for(int i=0;i<VectChargeCorrHSCPfirstXTalk8.size();i++) ChargeCorrHSCPfirstXTalk8+=VectChargeCorrHSCPfirstXTalk8.at(i);
                float ECorrHSCPfirstXTalk8=(3.61*pow(10,-9)*247)*ChargeCorrHSCPfirstXTalk8;


            //second XTalk study 

                std::vector<int> VectChargeCorrHSCPsecondXTalk4 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondXTalk4=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondXTalk4.size();i++) ChargeCorrHSCPsecondXTalk4+=VectChargeCorrHSCPsecondXTalk4.at(i);
                float ECorrHSCPsecondXTalk4=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondXTalk4;

                std::vector<int> VectChargeCorrHSCPsecondXTalk5 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.05,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondXTalk5=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondXTalk5.size();i++) ChargeCorrHSCPsecondXTalk5+=VectChargeCorrHSCPsecondXTalk5.at(i);
                float ECorrHSCPsecondXTalk5=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondXTalk5;

                std::vector<int> VectChargeCorrHSCPsecondXTalk6 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.06,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondXTalk6=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondXTalk6.size();i++) ChargeCorrHSCPsecondXTalk6+=VectChargeCorrHSCPsecondXTalk6.at(i);
                float ECorrHSCPsecondXTalk6=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondXTalk6;

                std::vector<int> VectChargeCorrHSCPsecondXTalk3 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.03,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondXTalk3=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondXTalk3.size();i++) ChargeCorrHSCPsecondXTalk3+=VectChargeCorrHSCPsecondXTalk3.at(i);
                float ECorrHSCPsecondXTalk3=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondXTalk3;

                std::vector<int> VectChargeCorrHSCPsecondXTalk2 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.02,20,25,layerLabel,2,RatioNClusterSat254,RatioNClusterSat255,0);
                float ChargeCorrHSCPsecondXTalk2=0.;
                for(int i=0;i<VectChargeCorrHSCPsecondXTalk2.size();i++) ChargeCorrHSCPsecondXTalk2+=VectChargeCorrHSCPsecondXTalk2.at(i);
                float ECorrHSCPsecondXTalk2=(3.61*pow(10,-9)*247)*ChargeCorrHSCPsecondXTalk2;

            //Etude ratio sat

                std::vector<int> VectChargeMyCorrRatioSat0 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,0.);
                float ChargeMyCorrRatioSat0=0.;
                for(int i=0;i<VectChargeMyCorrRatioSat0.size();i++) ChargeMyCorrRatioSat0+=VectChargeMyCorrRatioSat0.at(i);
                float EMyCorrRatioSat0=(3.61*pow(10,-9)*247)*ChargeMyCorrRatioSat0;

                std::vector<int> VectChargeMyCorrRatioSat10 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,0.1);
                float ChargeMyCorrRatioSat10=0.;
                for(int i=0;i<VectChargeMyCorrRatioSat10.size();i++) ChargeMyCorrRatioSat10+=VectChargeMyCorrRatioSat10.at(i);
                float EMyCorrRatioSat10=(3.61*pow(10,-9)*247)*ChargeMyCorrRatioSat10;

                std::vector<int> VectChargeMyCorrRatioSat20 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,0.2);
                float ChargeMyCorrRatioSat20=0.;
                for(int i=0;i<VectChargeMyCorrRatioSat20.size();i++) ChargeMyCorrRatioSat20+=VectChargeMyCorrRatioSat20.at(i);
                float EMyCorrRatioSat20=(3.61*pow(10,-9)*247)*ChargeMyCorrRatioSat20;

                std::vector<int> VectChargeMyCorrRatioSat30 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,0.3);
                float ChargeMyCorrRatioSat30=0.;
                for(int i=0;i<VectChargeMyCorrRatioSat30.size();i++) ChargeMyCorrRatioSat30+=VectChargeMyCorrRatioSat30.at(i);
                float EMyCorrRatioSat30=(3.61*pow(10,-9)*247)*ChargeMyCorrRatioSat30;

                std::vector<int> VectChargeMyCorrRatioSat40 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,0.4);
                float ChargeMyCorrRatioSat40=0.;
                for(int i=0;i<VectChargeMyCorrRatioSat40.size();i++) ChargeMyCorrRatioSat40+=VectChargeMyCorrRatioSat40.at(i);
                float EMyCorrRatioSat40=(3.61*pow(10,-9)*247)*ChargeMyCorrRatioSat40;

                std::vector<int> VectChargeMyCorrRatioSat50 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,0.5);
                float ChargeMyCorrRatioSat50=0.;
                for(int i=0;i<VectChargeMyCorrRatioSat50.size();i++) ChargeMyCorrRatioSat50+=VectChargeMyCorrRatioSat50.at(i);
                float EMyCorrRatioSat50=(3.61*pow(10,-9)*247)*ChargeMyCorrRatioSat50;

                std::vector<int> VectChargeMyCorrRatioSat60 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,0.6);
                float ChargeMyCorrRatioSat60=0.;
                for(int i=0;i<VectChargeMyCorrRatioSat60.size();i++) ChargeMyCorrRatioSat60+=VectChargeMyCorrRatioSat60.at(i);
                float EMyCorrRatioSat60=(3.61*pow(10,-9)*247)*ChargeMyCorrRatioSat60;

                std::vector<int> VectChargeMyCorrRatioSat70 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,0.7);
                float ChargeMyCorrRatioSat70=0.;
                for(int i=0;i<VectChargeMyCorrRatioSat70.size();i++) ChargeMyCorrRatioSat70+=VectChargeMyCorrRatioSat70.at(i);
                float EMyCorrRatioSat70=(3.61*pow(10,-9)*247)*ChargeMyCorrRatioSat70;

                std::vector<int> VectChargeMyCorrRatioSat80 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,0.8);
                float ChargeMyCorrRatioSat80=0.;
                for(int i=0;i<VectChargeMyCorrRatioSat80.size();i++) ChargeMyCorrRatioSat80+=VectChargeMyCorrRatioSat80.at(i);
                float EMyCorrRatioSat80=(3.61*pow(10,-9)*247)*ChargeMyCorrRatioSat80;

                std::vector<int> VectChargeMyCorrRatioSat90 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,0.9);
                float ChargeMyCorrRatioSat90=0.;
                for(int i=0;i<VectChargeMyCorrRatioSat90.size();i++) ChargeMyCorrRatioSat90+=VectChargeMyCorrRatioSat90.at(i);
                float EMyCorrRatioSat90=(3.61*pow(10,-9)*247)*ChargeMyCorrRatioSat90;

                std::vector<int> VectChargeMyCorrRatioSat100 = CrossTalkInv(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,layerLabel,3,RatioNClusterSat254,RatioNClusterSat255,1.);
                float ChargeMyCorrRatioSat100=0.;
                for(int i=0;i<VectChargeMyCorrRatioSat100.size();i++) ChargeMyCorrRatioSat100+=VectChargeMyCorrRatioSat100.at(i);
                float EMyCorrRatioSat100=(3.61*pow(10,-9)*247)*ChargeMyCorrRatioSat100;



                float QMyCorr = 1.038*EcorrModulGeom.ChargeCorr(charge,modulgeom,nstrips,nsat254,nsat255);

                dEdxWTEC.push_back((QMyCorr*(3.61*pow(10,-9)*247))/pathlength);
                if(subdetid!=6) dEdxWOTEC.push_back((QMyCorr*(3.61*pow(10,-9)*247))/pathlength);

                
                
                
            //--------
                    
                ncluster++;

                if(!(sat254 || sat255))
                {
                    DeDxHscpCorr.push_back(Erec/pathlength);
                    DeDxProposal.push_back(Erec/pathlength);
                    DeDxScdProposal.push_back(Erec/pathlength);
                    DeDxThirdProposal.push_back(Erec/pathlength);
                    DeDxFourthProposal.push_back(Erec/pathlength);
                    DeDxRecWOSat.push_back(Erec/pathlength);
                }

                if((sat254 || sat255))
                {

                distribdedx->Fill((QMyCorr*(3.61*pow(10,-9)*247)*pow(10,3))/pathlength);
                if(subdetid!=6)distribdedx2->Fill((QMyCorr*(3.61*pow(10,-9)*247)*pow(10,3))/pathlength);
                //if(subdetid==6 && NClustSat254+NClustSat255<=1)distribdedx2->Fill((QMyCorr*(3.61*pow(10,-9)*247)*pow(10,3))/pathlength);

                //if(nsat254+nsat255>=2)
                {

                    if(fakeclust) h1_ClustSatDeltaQ_nsat2_FakeClusters->Fill(DeltaQCalc);
                    if(!fakeclust) h1_ClustSatDeltaQ_nsat2_TrueClusters->Fill(DeltaQCalc);

                    //Study crit. MyCorrModulGeom

                    h2_ClustSat_QsimBetaGamma_Gluino2400->Fill(p/MGluino2400,charge_sim);
                    prof_ClustSat_QsimBetaGamma_Gluino2400->Fill(p/MGluino2400,charge_sim);

                    h2_ClustSat_QsimBetaGamma_Gluino1400->Fill(p/MGluino1400,charge_sim);
                    prof_ClustSat_QsimBetaGamma_Gluino1400->Fill(p/MGluino1400,charge_sim);

                    h2_ClustSat_QsimBetaGamma_Stau651->Fill(p/MStau651,charge_sim);
                    prof_ClustSat_QsimBetaGamma_Stau651->Fill(p/MStau651,charge_sim);

                    bool testRatioSat60 = false;
                    bool testRatioSat90 = false;
                    bool testThresholdSat25 = false;
                    bool testThresholdDeltaQ = false;
                    bool testThresholdDeltaQ30 = false;
                    bool testThresholdDeltaQ40 = false;
                    bool testThresholdDeltaQ50 = false;
                    bool testThresholdDeltaQ60 = false;
                    bool testThresholdDeltaQ20 = false;
                    bool testThresholdDeltaQ10 = false;
                    bool testThresholdDeltaQ15 = false;
                    bool testThresholdDeltaQ5 = false;
                    bool testChargeMinEstimated = false;
                    
                    
                    

                    if(RatioNClusterSat254+RatioNClusterSat255>=0.6) testRatioSat60=true;
                    if(RatioNClusterSat254+RatioNClusterSat255>=0.9) testRatioSat90=true;
                    if(Qminus>=25 && Qplus>=25) testThresholdSat25=true;
                    if(DeltaQCalc<=40 && nsat254+nsat255>=1) testThresholdDeltaQ=true;
                    //if(DeltaQCalc<=150 && nsat254+nsat255>1) testThresholdDeltaQ=true;
                    if(DeltaQCalc<=30) testThresholdDeltaQ30=true;
                    if(DeltaQCalc<=40) testThresholdDeltaQ40=true;
                    if(DeltaQCalc<=50) testThresholdDeltaQ50=true;
                    if(DeltaQCalc<=60) testThresholdDeltaQ60=true;
                    if(DeltaQCalc<=20) testThresholdDeltaQ20=true;
                    if(DeltaQCalc<=10) testThresholdDeltaQ10=true;
                    if(DeltaQCalc<=15) testThresholdDeltaQ15=true;
                    if(DeltaQCalc<=5) testThresholdDeltaQ5=true;

                    
                    

                    float QMyCorrModulGeom = 1.038*EcorrModulGeom.ChargeCorr(charge,modulgeom,nstrips,nsat254,nsat255);
                    //if(QMyCorrModulGeom>400) QMyCorrModulGeom = QMyCorrModulGeom*Calibp0+Calibp1;
                    bool testCorrModulGeom = EcorrModulGeom.TestCorr(modulgeom,nstrips,nsat254,nsat255);
                    string nopassing = (LabelModulGeom(modulgeom)+","+to_string(nstrips)+","+to_string(nsat254)+","+to_string(nsat255)).c_str();
                    if(!testCorrModulGeom && nstrips<=6 && nsat254<=2 && nsat255<=2) nopass.push_back(nopassing);
                    float QDeltaQ = charge;
                    float QDeltaQ30 = charge;
                    float QDeltaQ40 = charge;
                    float QDeltaQ50 = charge;
                    float QDeltaQ60 = charge;
                    float QDeltaQ20 = charge;
                    float QDeltaQ10 = charge;
                    float QDeltaQ15 = charge;
                    float QDeltaQ5 = charge;
                    float QThresholdSat25 = charge;
                    float QRatioSat60 = charge;
                    float QRatioSat90 = charge;
                    float QDeltaQRatioSat60 = charge;
                    float QThresholdSat25RatioSat60 = charge;
                    float QDeltaQRatioSat90 = charge;
                    float QThresholdSat25RatioSat90 = charge;
                    float QThresholdSat25DeltaQ = charge;
                    float QThresholdSat25DeltaQRatioSat60 = charge;
                    float QThresholdSat25DeltaQRatioSat90 = charge;
                    float QChargeMinEstimated = charge;

                    h1_ClustSat_Qsim->Fill(charge_sim);
                    
                    
                    
                    h1_ClustSat_MyCorrModulGeom_WihtoutCrit->Fill((QMyCorrModulGeom-charge_sim)/charge_sim);
                    profile_ClustSat_MyCorrModulGeom_WihtoutCrit->Fill(charge_sim,QMyCorrModulGeom);
                    nPassingNoCrit++;

                    float Qnewproc = 1.03*sc2.ChargeCorrected(charge,modulgeom,nstrips,nsat254+nsat255);

                    DeDxFourthProposal.push_back(QtoE(Qnewproc)/pathlength);

                    if(Qminus>=25 && Qplus>=25 && DeltaQCalc<=40 && Qnewproc-charge!=0) DeDxThirdProposal.push_back(QtoE(Qnewproc)/pathlength);

                    if(Qminus>=25 && Qplus>=25 && testRatioSat60 && DeltaQCalc<=40 && Qnewproc-charge!=0)
                    // if(Qminus>=25 && Qplus>=25 && DeltaQCalc<=40 && Qnewproc-charge!=0)
                    // if(testCorrHscp)
                    {
                        //if(nsat254+nsat255==2 && modulgeom>=5)
                        {
                        h1_NewProcedure->Fill((Qnewproc-charge_sim)/charge_sim);
                        profile_NewProcedure->Fill(charge_sim,Qnewproc);
                        nPassingNewProcedure++;
                        }
                        DeDxProposal.push_back(QtoE(Qnewproc)/pathlength);
                        DeDxScdProposal.push_back(QtoE(Qnewproc)/pathlength);
                    }
                    else DeDxScdProposal.push_back(Erec/pathlength);

                    if(Qminus>=25 && Qplus>=25 && testRatioSat60 && DeltaQCalc<=40 && QMyCorrModulGeom-charge!=0)
                    {
                        /*if(modulgeom>=5 && nsat254+nsat255==1 && nstrips>3)
                        {
                            h1_Proposal->Fill((QMyCorrModulGeom-charge_sim)/charge_sim);
                            profile_Proposal->Fill(charge_sim,QMyCorrModulGeom);
                        }*/
                        /*if(modulgeom<=4 && nsat254+nsat255==1 && nsat255!=1)
                        {
                            h1_Proposal->Fill((QMyCorrModulGeom-charge_sim)/charge_sim);
                            profile_Proposal->Fill(charge_sim,QMyCorrModulGeom);
                        }
                        if(modulgeom<=4 && nsat254+nsat255==2)
                        {
                            h1_Proposal->Fill((QMyCorrModulGeom-charge_sim)/charge_sim);
                            profile_Proposal->Fill(charge_sim,QMyCorrModulGeom);
                        }*/
                        //if(modulgeom>4 && modulgeom<8 && nstrips==5 && nsat254+nsat255==1 && nsat255!=1)
                        {
                            h1_Proposal->Fill((QMyCorrModulGeom-charge_sim)/charge_sim);
                            profile_Proposal->Fill(charge_sim,QMyCorrModulGeom);
                        }
                    }




                    if(testCorrHscp)
                    {
                        h_ClustSat_HscpCorrNoInv->Fill((ChargeCorrHSCP_WithoutInversion-charge_sim)/charge_sim);
                        profile_ClustSat_HscpNoInv->Fill(charge_sim,ChargeCorrHSCP_WithoutInversion);
                        DeDxHscpCorr.push_back(QtoE(ChargeCorrHSCP_WithoutInversion)/pathlength);
                        if(abs(eta)<=1.4)
                        {
                            h1_ClustSat_HscpNoInv_CutEta->Fill((ChargeCorrHSCP_WithoutInversion-charge_sim)/charge_sim);
                            profile_ClustSat_HscpNoInv_CutEta->Fill(charge_sim,ChargeCorrHSCP_WithoutInversion);
                        }
                        if(layerLabel==1)
                        {
                            h1_ClustSat_HscpNoInv_TIB1->Fill((ChargeCorrHSCP_WithoutInversion-charge_sim)/charge_sim);
                            profile_ClustSat_HscpNoInv_TIB1->Fill(charge_sim,ChargeCorrHSCP_WithoutInversion);
                        }
                        if(modulgeom==3)
                        {
                            h1_ClustSat_HscpNoInv_OB1->Fill((ChargeCorrHSCP_WithoutInversion-charge_sim)/charge_sim);
                            profile_ClustSat_HscpNoInv_OB1->Fill(charge_sim,ChargeCorrHSCP_WithoutInversion);
                        }
                        if(abs(eta)>1.4)
                        {
                            h1_ClustSat_HscpNoInv_CutEtaSup->Fill((ChargeCorrHSCP_WithoutInversion-charge_sim)/charge_sim);
                            profile_ClustSat_HscpNoInv_CutEtaSup->Fill(charge_sim,ChargeCorrHSCP_WithoutInversion);
                        }
                        if(subdetid==6)
                        {
                            h1_ClustSat_HscpNoInv_TEC->Fill((ChargeCorrHSCP_WithoutInversion-charge_sim)/charge_sim);
                            profile_ClustSat_HscpNoInv_TEC->Fill(charge_sim,ChargeCorrHSCP_WithoutInversion);
                        }
                        if(subdetid==4)
                        {
                            h1_ClustSat_HscpNoInv_TID->Fill((ChargeCorrHSCP_WithoutInversion-charge_sim)/charge_sim);
                            profile_ClustSat_HscpNoInv_TID->Fill(charge_sim,ChargeCorrHSCP_WithoutInversion);
                        }
                        if(subdetid==3 || subdetid==5)
                        {
                            h1_ClustSat_HscpNoInv_Barrel->Fill((ChargeCorrHSCP_WithoutInversion-charge_sim)/charge_sim);
                            profile_ClustSat_HscpNoInv_Barrel->Fill(charge_sim,ChargeCorrHSCP_WithoutInversion);
                        }
                        nPassingHSCPtest++;
                    }
                    if(testThresholdDeltaQ)
                    {
                        QDeltaQ=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_DeltaQ->Fill((QDeltaQ-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_DeltaQ->Fill(charge_sim,QDeltaQ);
                        nPassingDeltaQ++;
                    }
                    if(testThresholdDeltaQ40)
                    {
                        QDeltaQ40=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_DeltaQ40->Fill((QDeltaQ40-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_DeltaQ40->Fill(charge_sim,QDeltaQ40);
                        h2_ClustSat_MyCorrModulGeom_DeltaQ40->Fill(charge_sim,QDeltaQ40);
                        nPassingDeltaQ40++;
                    }
                    if(testThresholdDeltaQ30)
                    {
                        QDeltaQ30=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_DeltaQ30->Fill((QDeltaQ30-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_DeltaQ30->Fill(charge_sim,QDeltaQ30);
                        nPassingDeltaQ30++;
                    }
                    if(testThresholdDeltaQ50)
                    {
                        QDeltaQ50=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_DeltaQ50->Fill((QDeltaQ50-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_DeltaQ50->Fill(charge_sim,QDeltaQ50);
                        nPassingDeltaQ50++;
                    }
                    if(testThresholdDeltaQ60)
                    {
                        QDeltaQ60=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_DeltaQ60->Fill((QDeltaQ60-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_DeltaQ60->Fill(charge_sim,QDeltaQ60);
                        nPassingDeltaQ60++;
                    }
                    if(testThresholdDeltaQ20)
                    {
                        QDeltaQ20=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_DeltaQ20->Fill((QDeltaQ20-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_DeltaQ20->Fill(charge_sim,QDeltaQ20);
                        nPassingDeltaQ20++;
                    }
                    if(testThresholdDeltaQ10)
                    {
                        QDeltaQ10=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_DeltaQ10->Fill((QDeltaQ10-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_DeltaQ10->Fill(charge_sim,QDeltaQ10);
                        nPassingDeltaQ10++;
                    }
                    if(testThresholdDeltaQ15)
                    {
                        QDeltaQ15=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_DeltaQ15->Fill((QDeltaQ15-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_DeltaQ15->Fill(charge_sim,QDeltaQ15);
                        nPassingDeltaQ15++;
                    }
                    if(nstrips>=6)
                    {
                        QDeltaQ5=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_DeltaQ5->Fill((QDeltaQ5-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_DeltaQ5->Fill(charge_sim,QDeltaQ5);
                        nPassingDeltaQ5++;
                    }
                    if(nstrips==5)
                    {
                        QThresholdSat25=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_ThresholdSat25->Fill((QThresholdSat25-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_ThresholdSat25->Fill(charge_sim,QThresholdSat25);
                        nPassingThresholdSat25++;
                    }
                    if(nstrips==3) 
                    {
                        QRatioSat60=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_RatioSat60->Fill((QRatioSat60-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_RatioSat60->Fill(charge_sim,QRatioSat60);
                        nPassingRatioSat60++;
                    }
                    if(nstrips==4) 
                    {
                        QRatioSat90=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_RatioSat90->Fill((QRatioSat90-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_RatioSat90->Fill(charge_sim,QRatioSat90);
                        nPassingRatioSat90++;
                    }
                    if(testThresholdDeltaQ && testRatioSat60) 
                    {
                        QDeltaQRatioSat60=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat60->Fill((QDeltaQRatioSat60-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_DeltaQRatioSat60->Fill(charge_sim,QDeltaQRatioSat60);
                        nPassingDeltaQRatioSat60++;
                    }
                    if(testThresholdDeltaQ && testRatioSat90) 
                    {
                        QDeltaQRatioSat90=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat90->Fill((QDeltaQRatioSat90-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_DeltaQRatioSat90->Fill(charge_sim,QDeltaQRatioSat90);
                        nPassingDeltaQRatioSat90++;
                    }
                    if(testThresholdSat25 && testRatioSat60) 
                    {
                        QThresholdSat25RatioSat60=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat60->Fill((QThresholdSat25RatioSat60-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat60->Fill(charge_sim,QThresholdSat25RatioSat60);
                        nPassingThresholdSat25RatioSat60++;
                    }
                    if(testThresholdSat25 && testRatioSat90) 
                    {
                        QThresholdSat25RatioSat90=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat90->Fill((QThresholdSat25RatioSat90-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat90->Fill(charge_sim,QThresholdSat25RatioSat90);
                        nPassingThresholdSat25RatioSat90++;
                    }
                    if(testThresholdSat25 && testThresholdDeltaQ) 
                    {
                        QThresholdSat25DeltaQ=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQ->Fill((QThresholdSat25DeltaQ-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQ->Fill(charge_sim,QThresholdSat25DeltaQ);
                        nPassingThresholdSat25DeltaQ++;
                    }
                    if(testThresholdSat25 && testThresholdDeltaQ && testRatioSat60) 
                    {
                        
                        QThresholdSat25DeltaQRatioSat60=QMyCorrModulGeom;
                        if(subdetid!=6)
                        {
                            h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60->Fill((QThresholdSat25DeltaQRatioSat60-charge_sim)/charge_sim);
                            profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60->Fill(charge_sim,QThresholdSat25DeltaQRatioSat60);
                        }
                        if(modulgeom>=5 && nsat254+nsat255==1)
                        {
                            h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60->Fill((QThresholdSat25DeltaQRatioSat60-charge_sim)/charge_sim);
                            profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60->Fill(charge_sim,QThresholdSat25DeltaQRatioSat60);
                        }
                        
                        if(abs(eta)<=1.4) 
                        {
                            h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEta->Fill((QThresholdSat25DeltaQRatioSat60-charge_sim)/charge_sim);
                            profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEta->Fill(charge_sim,QThresholdSat25DeltaQRatioSat60);
                        }
                        if(layerLabel==1)
                        {
                            h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TIB1->Fill((QThresholdSat25DeltaQRatioSat60-charge_sim)/charge_sim);
                            profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TIB1->Fill(charge_sim,QThresholdSat25DeltaQRatioSat60);
                        }
                        if(modulgeom==3)
                        {
                            h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_OB1->Fill((QThresholdSat25DeltaQRatioSat60-charge_sim)/charge_sim);
                            profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_OB1->Fill(charge_sim,QThresholdSat25DeltaQRatioSat60);
                        }
                        if(abs(eta)>1.4)
                        {
                            h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEtaSup->Fill((QThresholdSat25DeltaQRatioSat60-charge_sim)/charge_sim);
                            profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEtaSup->Fill(charge_sim,QThresholdSat25DeltaQRatioSat60);
                        }
                        if(subdetid==6 && nsat254+nsat255==1 && QThresholdSat25DeltaQRatioSat60>0 && DeltaQCalc<=25)
                        {
                            h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TEC->Fill((QThresholdSat25DeltaQRatioSat60-charge_sim)/charge_sim);
                            profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TEC->Fill(charge_sim,QThresholdSat25DeltaQRatioSat60);
                            //if((QThresholdSat25DeltaQRatioSat60-charge_sim)/charge_sim<-0.3) 
                            {
                                distribdeltaQ->Fill(DeltaQCalc,(QThresholdSat25DeltaQRatioSat60-charge_sim)/charge_sim);
                                distribminQ->Fill(minQ,(QThresholdSat25DeltaQRatioSat60-charge_sim)/charge_sim);
                                //cout<<charge_sim<<"    "<<charge<<"      "<<QThresholdSat25DeltaQRatioSat60<<"     "<<(QThresholdSat25DeltaQRatioSat60-charge_sim)/charge_sim<<endl;
                                //cout<<"modul "<<LabelModulGeom(modulgeom)<<"   Layer  "<<LabelLayer(layerLabel)<<"   NStrip    "<<nstrips<<"   nsat254   "<<nsat254<<"    nsat255    "<<nsat255<<endl;
                            }
                            
                        }
                        if(subdetid==4 && nsat254+nsat255==1 && QThresholdSat25DeltaQRatioSat60>0 && DeltaQCalc<=25)
                        {
                            h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TID->Fill((QThresholdSat25DeltaQRatioSat60-charge_sim)/charge_sim);
                            profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TID->Fill(charge_sim,QThresholdSat25DeltaQRatioSat60);
                        }
                        if(subdetid==3 || subdetid==5)
                        {
                            h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_Barrel->Fill((QThresholdSat25DeltaQRatioSat60-charge_sim)/charge_sim);
                            profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_Barrel->Fill(charge_sim,QThresholdSat25DeltaQRatioSat60);
                        }
                        nPassingThresholdSat25DeltaQRatioSat60++;
                    }
                    if(testThresholdSat25 && testThresholdDeltaQ && testRatioSat90) 
                    {
                        QThresholdSat25DeltaQRatioSat90=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat90->Fill((QThresholdSat25DeltaQRatioSat90-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat90->Fill(charge_sim,QThresholdSat25DeltaQRatioSat90);
                        nPassingThresholdSat25DeltaQRatioSat90++;
                    }
                    if(QMyCorrModulGeom>ChargeMinEstimated) 
                    {
                        QChargeMinEstimated=QMyCorrModulGeom;
                        h1_ClustSat_MyCorrModulGeom_ChargeMinEstimated->Fill((QChargeMinEstimated-charge_sim)/charge_sim);
                        profile_ClustSat_MyCorrModulGeom_ChargeMinEstimated->Fill(charge_sim,QChargeMinEstimated);
                        nPassingChargeMinEstimated++;
                    }

                }
                    //-------------------------------------



                    if(!fakeclust) h1_DeltaQ_False->Fill(DeltaQ);
                    if(fakeclust) h1_DeltaQ_True->Fill(DeltaQ);
countSatCluster++;


                    //Criteria studies

std::vector<int> VectChargeCorrHSCP_WithoutInversion_StudyCrit = CrossTalkInvStudy(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,40,layerLabel,1,RatioNClusterSat254,RatioNClusterSat255,0,testThresholdSat,testThresholdDiff,testBothThreshold);
float ChargeCorrHSCP_WithoutInversion_StudyCrit=0.;
for(int i=0;i<VectChargeCorrHSCP_WithoutInversion_StudyCrit.size();i++) ChargeCorrHSCP_WithoutInversion_StudyCrit+=VectChargeCorrHSCP_WithoutInversion_StudyCrit.at(i);
float ECorrHSCP_WithoutInversion_StudyCrit=(3.61*pow(10,-9)*247)*ChargeCorrHSCP_WithoutInversion_StudyCrit;
if(testThresholdSat ) 
{
    h1_ClustSat_testThresholdSat_HscpCorr->Fill((ECorrHSCP_WithoutInversion_StudyCrit-Eloss)/Eloss);
    h1_FakeClusters_testThresholdSat->Fill(Ratio_fakeclust);
    countThresholdSat++;
}
if(testThresholdDiff ) 
{
    h1_ClustSat_testThresholdDiff_HscpCorr->Fill((ECorrHSCP_WithoutInversion_StudyCrit-Eloss)/Eloss);
    h1_FakeClusters_testThresholdDiff->Fill(Ratio_fakeclust);
    countThresholdDiff++;
}
if(testBothThreshold) 
{
    h1_ClustSat_testBothThreshold_HscpCorr->Fill((ECorrHSCP_WithoutInversion_StudyCrit-Eloss)/Eloss);
    countBothThreshold++;
}  

//DeltaQ criterion study

std::vector<int> VectChargeCorrHSCP_WithoutInversion_DeltaQ = CrossTalkInvStudy(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,40,layerLabel,1,RatioNClusterSat254,RatioNClusterSat255,0,testThresholdSat,testThresholdDiff,testBothThreshold);
float ChargeCorrHSCP_WithoutInversion_DeltaQ=0.;
for(int i=0;i<VectChargeCorrHSCP_WithoutInversion_DeltaQ.size();i++) ChargeCorrHSCP_WithoutInversion_DeltaQ+=VectChargeCorrHSCP_WithoutInversion_DeltaQ.at(i);
float ECorrHSCP_WithoutInversion_DeltaQ=(3.61*pow(10,-9)*247)*ChargeCorrHSCP_WithoutInversion_DeltaQ;
if(testThresholdDiff ) h1_ClustSat_DeltaQ->Fill((ECorrHSCP_WithoutInversion_DeltaQ-Eloss)/Eloss);

std::vector<int> VectChargeCorrHSCP_WithoutInversion_DeltaQ35 = CrossTalkInvStudy(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,35,layerLabel,1,RatioNClusterSat254,RatioNClusterSat255,0,testThresholdSat,testThresholdDiff,testBothThreshold);
float ChargeCorrHSCP_WithoutInversion_DeltaQ35=0.;
for(int i=0;i<VectChargeCorrHSCP_WithoutInversion_DeltaQ35.size();i++) ChargeCorrHSCP_WithoutInversion_DeltaQ35+=VectChargeCorrHSCP_WithoutInversion_DeltaQ35.at(i);
float ECorrHSCP_WithoutInversion_DeltaQ35=(3.61*pow(10,-9)*247)*ChargeCorrHSCP_WithoutInversion_DeltaQ35;
if(testThresholdDiff ) h1_ClustSat_DeltaQ35->Fill((ECorrHSCP_WithoutInversion_DeltaQ35-Eloss)/Eloss);

std::vector<int> VectChargeCorrHSCP_WithoutInversion_DeltaQ30 = CrossTalkInvStudy(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,30,layerLabel,1,RatioNClusterSat254,RatioNClusterSat255,0,testThresholdSat,testThresholdDiff,testBothThreshold);
float ChargeCorrHSCP_WithoutInversion_DeltaQ30=0.;
for(int i=0;i<VectChargeCorrHSCP_WithoutInversion_DeltaQ30.size();i++) ChargeCorrHSCP_WithoutInversion_DeltaQ30+=VectChargeCorrHSCP_WithoutInversion_DeltaQ30.at(i);
float ECorrHSCP_WithoutInversion_DeltaQ30=(3.61*pow(10,-9)*247)*ChargeCorrHSCP_WithoutInversion_DeltaQ30;
if(testThresholdDiff ) h1_ClustSat_DeltaQ30->Fill((ECorrHSCP_WithoutInversion_DeltaQ30-Eloss)/Eloss);

std::vector<int> VectChargeCorrHSCP_WithoutInversion_DeltaQ25 = CrossTalkInvStudy(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,25,layerLabel,1,RatioNClusterSat254,RatioNClusterSat255,0,testThresholdSat,testThresholdDiff,testBothThreshold);
float ChargeCorrHSCP_WithoutInversion_DeltaQ25=0.;
for(int i=0;i<VectChargeCorrHSCP_WithoutInversion_DeltaQ25.size();i++) ChargeCorrHSCP_WithoutInversion_DeltaQ25+=VectChargeCorrHSCP_WithoutInversion_DeltaQ25.at(i);
float ECorrHSCP_WithoutInversion_DeltaQ25=(3.61*pow(10,-9)*247)*ChargeCorrHSCP_WithoutInversion_DeltaQ25;
if(testThresholdDiff ) h1_ClustSat_DeltaQ25->Fill((ECorrHSCP_WithoutInversion_DeltaQ25-Eloss)/Eloss);

std::vector<int> VectChargeCorrHSCP_WithoutInversion_DeltaQ20 = CrossTalkInvStudy(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,20,layerLabel,1,RatioNClusterSat254,RatioNClusterSat255,0,testThresholdSat,testThresholdDiff,testBothThreshold);
float ChargeCorrHSCP_WithoutInversion_DeltaQ20=0.;
for(int i=0;i<VectChargeCorrHSCP_WithoutInversion_DeltaQ20.size();i++) ChargeCorrHSCP_WithoutInversion_DeltaQ20+=VectChargeCorrHSCP_WithoutInversion_DeltaQ20.at(i);
float ECorrHSCP_WithoutInversion_DeltaQ20=(3.61*pow(10,-9)*247)*ChargeCorrHSCP_WithoutInversion_DeltaQ20;
if(testThresholdDiff ) h1_ClustSat_DeltaQ20->Fill((ECorrHSCP_WithoutInversion_DeltaQ20-Eloss)/Eloss);

std::vector<int> VectChargeCorrHSCP_WithoutInversion_DeltaQ15 = CrossTalkInvStudy(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,15,layerLabel,1,RatioNClusterSat254,RatioNClusterSat255,0,testThresholdSat,testThresholdDiff,testBothThreshold);
float ChargeCorrHSCP_WithoutInversion_DeltaQ15=0.;
for(int i=0;i<VectChargeCorrHSCP_WithoutInversion_DeltaQ15.size();i++) ChargeCorrHSCP_WithoutInversion_DeltaQ15+=VectChargeCorrHSCP_WithoutInversion_DeltaQ15.at(i);
float ECorrHSCP_WithoutInversion_DeltaQ15=(3.61*pow(10,-9)*247)*ChargeCorrHSCP_WithoutInversion_DeltaQ15;
if(testThresholdDiff ) h1_ClustSat_DeltaQ15->Fill((ECorrHSCP_WithoutInversion_DeltaQ15-Eloss)/Eloss);


/*for(int DQ=0;DQ<250;DQ++)
{
    std::vector<int> VectChargeCorrHSCP_WithoutInversion_DeltaQStudy = CrossTalkInvStudy(EcorrLayer,VectChargeStrip,0.10,0.04,20,25,DQ,layerLabel,1,RatioNClusterSat254,RatioNClusterSat255,0,testThresholdSat,testThresholdDiff,testBothThreshold);
    float ChargeCorrHSCP_WithoutInversion_DeltaQStudy=0.;
    for(int i=0;i<VectChargeCorrHSCP_WithoutInversion_DeltaQStudy.size();i++) ChargeCorrHSCP_WithoutInversion_DeltaQStudy+=VectChargeCorrHSCP_WithoutInversion_DeltaQStudy.at(i);
    float ECorrHSCP_WithoutInversion_DeltaQStudy=(3.61*pow(10,-9)*247)*ChargeCorrHSCP_WithoutInversion_DeltaQStudy;
    if(testBothThreshold ) h2_ClustSat_DeltaQ->Fill((ECorrHSCP_WithoutInversion_DeltaQStudy-Eloss)/Eloss,DQ);
}

for(int TSat=-10;TSat<140;TSat++)
{
    std::vector<int> VectChargeCorrHSCP_WithoutInversion_ThresholdSatStudy = CrossTalkInvStudy(EcorrLayer,VectChargeStrip,0.10,0.04,20,TSat,40,layerLabel,1,RatioNClusterSat254,RatioNClusterSat255,0,testThresholdSat,testThresholdDiff,testBothThreshold);
    float ChargeCorrHSCP_WithoutInversion_ThresholdSatStudy=0.;
    for(int i=0;i<VectChargeCorrHSCP_WithoutInversion_ThresholdSatStudy.size();i++) ChargeCorrHSCP_WithoutInversion_ThresholdSatStudy+=VectChargeCorrHSCP_WithoutInversion_ThresholdSatStudy.at(i);
    float ECorrHSCP_WithoutInversion_ThresholdSatStudy=(3.61*pow(10,-9)*247)*ChargeCorrHSCP_WithoutInversion_ThresholdSatStudy;
    if(testBothThreshold ) h2_ClustSat_ThresholdSatCrit->Fill((ECorrHSCP_WithoutInversion_ThresholdSatStudy-Eloss)/Eloss,TSat);
}*/


                    if(EMyCorrRatioSat0==Erec) counterpass++;




                    nSat++;
                    if(nsat254+nsat255==1) sat1++;
                    if(nsat254+nsat255==2) sat2++;
                    if(nsat254+nsat255>=3) sat3++;
                    if(RatioNClusterSat254+RatioNClusterSat255>=ratsat) nRatSatTest++;
                    if(nsat254+nsat255>=2 && RatioNClusterSat254+RatioNClusterSat255>=ratsat) sat2RatSatTest++;


                    //Ih study

                    dEdxsim.push_back(Eloss/pathlength);
                    dEdxMyCorr.push_back(EMyCorr_WithoutInversion/pathlength);
                    dEdxHscpCorrWOinv.push_back(ECorrHSCP_WithoutInversion/pathlength);
                    dEdxHscpCorr.push_back(ECorrHSCP/pathlength);


                    h2_ClustSat_EsimEMyCorrRatioSat->Fill(Eloss,EMyCorr_WithoutInversion_RatioSat);
                    h2_ClustSat_EsimECorrHscp->Fill(Eloss,ECorrHSCP_WithoutInversion);
                    prof_ClustSat_EsimEMyCorrRatioSat->Fill(Eloss,EMyCorr_WithoutInversion_RatioSat);
                    prof_ClustSat_EsimECorrHscp->Fill(Eloss,ECorrHSCP_WithoutInversion);

                    h_ClustSat_Esim->Fill(Eloss/(3.61*pow(10,-9)*247));
                    h_ClustSat_Erec->Fill(Erec/(3.61*pow(10,-9)*247));
                    h_ClustSat_ClusterSize->Fill(nstrips);
                    h_ClustSat_Sat254->Fill(nsat254);
                    h_ClustSat_Sat255->Fill(nsat255);
                    h_ClustSat_NSat->Fill(nsat254+nsat255);
                    h_ClustSat_layerlabel->Fill(layerLabel);
                    h_ClustSat_Eta->Fill(abs(eta));

                    if(((EMyCorr_WithoutInversion_RatioSat-Eloss)/Eloss)>0.25)
                    {

                        

                        h_ClustSat_SurCorr_Esim->Fill(Eloss/(3.61*pow(10,-9)*247));
                        h_ClustSat_SurCorr_Erec->Fill(Erec/(3.61*pow(10,-9)*247));
                        h_ClustSat_SurCorr_ClusterSize->Fill(nstrips);
                        h_ClustSat_SurCorr_Sat254->Fill(nsat254);
                        h_ClustSat_SurCorr_Sat255->Fill(nsat255);
                        h_ClustSat_SurCorr_NSat->Fill(nsat254+nsat255);
                        h_ClustSat_SurCorr_layerlabel->Fill(layerLabel);
                        h_ClustSat_SurCorr_Eta->Fill(abs(eta));

                        


                    }
                    if(((ECorrHSCP_WithoutInversion-Eloss)/Eloss)>0.25) surcorrHSCP++;

                    if(((EMyCorr_WithoutInversion_RatioSat-Eloss)/Eloss)>0.25)
                    {
                        surcorr++;
                        testSurCorr=true;
                        h2_ClustSatSurCorr_EsimECorrHscp->Fill(Eloss,ECorrHSCP_WithoutInversion);
                        h2_ClustSatSurCorr_EsimEMyCorrRatioSat->Fill(Eloss,EMyCorr_WithoutInversion_RatioSat);
                        prof_ClustSatSurCorr_EsimECorrHscp->Fill(Eloss,ECorrHSCP_WithoutInversion);
                        prof_ClustSatSurCorr_EsimEMyCorrRatioSat->Fill(Eloss,EMyCorr_WithoutInversion_RatioSat);
                    }

                    if(((EMyCorr_WithoutInversion_RatioSat-Eloss)/Eloss)<0.25)
                    {
                        h2_ClustSatNoSurCorr_EsimECorrHscp->Fill(Eloss,ECorrHSCP_WithoutInversion);
                        h2_ClustSatNoSurCorr_EsimEMyCorrRatioSat->Fill(Eloss,EMyCorr_WithoutInversion_RatioSat);
                        prof_ClustSatNoSurCorr_EsimECorrHscp->Fill(Eloss,ECorrHSCP_WithoutInversion);
                        prof_ClustSatNoSurCorr_EsimEMyCorrRatioSat->Fill(Eloss,EMyCorr_WithoutInversion_RatioSat);
                    }

                    if(nsat254+nsat255>1)
                    {

                        testMoreThanOne=true;

                        h2_ClustSatMoreThanOneSat_EsimECorrHscp->Fill(Eloss,ECorrHSCP_WithoutInversion);
                        h2_ClustSatMoreThanOneSat_EsimEMyCorrRatioSat->Fill(Eloss,EMyCorr_WithoutInversion_RatioSat);
                        prof_ClustSatMoreThanOneSat_EsimECorrHscp->Fill(Eloss,ECorrHSCP_WithoutInversion);
                        prof_ClustSatMoreThanOneSat_EsimEMyCorrRatioSat->Fill(Eloss,EMyCorr_WithoutInversion_RatioSat);
                    }

                    if(nsat254+nsat255==1)
                    {

                        h2_ClustSatOneSat_EsimECorrHscp->Fill(Eloss,ECorrHSCP_WithoutInversion);
                        h2_ClustSatOneSat_EsimEMyCorrRatioSat->Fill(Eloss,EMyCorr_WithoutInversion_RatioSat);
                        prof_ClustSatOneSat_EsimECorrHscp->Fill(Eloss,ECorrHSCP_WithoutInversion);
                        prof_ClustSatOneSat_EsimEMyCorrRatioSat->Fill(Eloss,EMyCorr_WithoutInversion_RatioSat);
                    }

                    if(testCorrHscp)
                    {
                        testHscptrack=true;
                        h2_ClustSatPassingHscpTest_EsimECorrHscp->Fill(Eloss,ECorrHSCP_WithoutInversion);
                        h2_ClustSatPassingHscpTest_EsimEMyCorrRatioSat->Fill(Eloss,EMyCorr_WithoutInversion_RatioSat);
                        prof_ClustSatPassingHscpTest_EsimECorrHscp->Fill(Eloss,ECorrHSCP_WithoutInversion);
                        prof_ClustSatPassingHscpTest_EsimEMyCorrRatioSat->Fill(Eloss,EMyCorr_WithoutInversion_RatioSat);
                    }

                    h_ClustSat_MyCorrRatioSat0->Fill((EMyCorrRatioSat0-Eloss)/Eloss);
                    h_ClustSat_MyCorrRatioSat10->Fill((EMyCorrRatioSat10-Eloss)/Eloss);
                    h_ClustSat_MyCorrRatioSat20->Fill((EMyCorrRatioSat20-Eloss)/Eloss);
                    h_ClustSat_MyCorrRatioSat30->Fill((EMyCorrRatioSat30-Eloss)/Eloss);
                    h_ClustSat_MyCorrRatioSat40->Fill((EMyCorrRatioSat40-Eloss)/Eloss);
                    h_ClustSat_MyCorrRatioSat50->Fill((EMyCorrRatioSat50-Eloss)/Eloss);
                    h_ClustSat_MyCorrRatioSat60->Fill((EMyCorrRatioSat60-Eloss)/Eloss);
                    h_ClustSat_MyCorrRatioSat70->Fill((EMyCorrRatioSat70-Eloss)/Eloss);
                    h_ClustSat_MyCorrRatioSat80->Fill((EMyCorrRatioSat80-Eloss)/Eloss);
                    h_ClustSat_MyCorrRatioSat90->Fill((EMyCorrRatioSat90-Eloss)/Eloss);
                    h_ClustSat_MyCorrRatioSat100->Fill((EMyCorrRatioSat100-Eloss)/Eloss);

                    if(testCorrHscp)
                    {
                        nHscpTest++;
                        if(nsat254+nsat255==1) sat1HscpTest++;
                        if(RatioNClusterSat254+RatioNClusterSat255>=ratsat) nHscpAndRatSatTest++;

                        h_ClustSatAndHscpTest_MyFullCorr->Fill((EMyCorr_WithoutInversion_RatioSat-Eloss)/Eloss);
                        h_ClustSatAndHscpTest_HscpCorrNoInv->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                        h_ClustSatAndHscpTest_WithoutCorr->Fill((Erec-Eloss)/Eloss);
                        h_ClustSatAndHscpTest_MyFullCorrModulGeom->Fill((EMyCorrModulGeom_WithoutInversion_RatioSat-Eloss)/Eloss);
                    }
//if(EcorrLayer.TestCorr(layerLabel,nstrips,nsat254,nsat255))
{
                    h_ClustSat_MyFullCorr->Fill((EMyCorr_WithoutInversion_RatioSat-Eloss)/Eloss);
                    
                    h_ClustSat_WithoutCorr->Fill((Erec-Eloss)/Eloss);

                    
                    hMyCorrModulGoem->Fill((EMyCorrModulGeom_WithoutInversion_RatioSat-Eloss)/Eloss);
                    h_ClustSat_CorrHSCP_WithoutInversion->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    h_ClustSat_CorrHSCP->Fill((ECorrHSCP-Eloss)/Eloss);            
                    h_ClustSat_MyCorr_WithoutInversion->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    h_ClustSat_MyCorr->Fill((EMyCorr-Eloss)/Eloss);
                    h_ClustSat_MyCorr_WithoutInversion_RatioSat->Fill((EMyCorr_WithoutInversion_RatioSat-Eloss)/Eloss);
                    h_ClustSat_MyCorr_RatioSat->Fill((EMyCorr_RatioSat-Eloss)/Eloss);
}
                    h2_ClustSat_MyCorrWOinv_p->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss,p);
                    h2_ClustSat_MyCorrWOinv_eta->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss,abs(eta));

                    h_ClustSat_DistribP->Fill(p);
                    h_ClustSat_DistribEta->Fill(abs(eta));
                    h_ClustSat_DistribSize->Fill(nstrips);
                    h_ClustSat_DistribNSat->Fill(nsat254+nsat255);

                    if(p<=500) h_ClustSat_MyCorrWOinv_0p500->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(p<=1000 && p>500) h_ClustSat_MyCorrWOinv_500p1000->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(p<=1500 && p>1000) h_ClustSat_MyCorrWOinv_1000p1500->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(p<=2000 && p>1500) h_ClustSat_MyCorrWOinv_1500p2000->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(p<=4000 && p>2000) h_ClustSat_MyCorrWOinv_2000p4000->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);

                    if(abs(eta)<=0.5) h_ClustSat_MyCorrWOinv_0eta05->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=1 && abs(eta)>0.5) h_ClustSat_MyCorrWOinv_05eta10->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=1.5 && abs(eta)>1) h_ClustSat_MyCorrWOinv_10eta15->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=2.5 && abs(eta)>1.5) h_ClustSat_MyCorrWOinv_15eta25->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=2.5 && abs(eta)>2) h_ClustSat_MyCorrWOinv_20eta25->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);

                    if(p<=500) h_ClustSat_HscpCorrWOinv_0p500->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(p<=1000 && p>500) h_ClustSat_HscpCorrWOinv_500p1000->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(p<=1500 && p>1000) h_ClustSat_HscpCorrWOinv_1000p1500->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(p<=2000 && p>1500) h_ClustSat_HscpCorrWOinv_1500p2000->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(p<=4000 && p>2000) h_ClustSat_HscpCorrWOinv_2000p4000->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);

                    if(abs(eta)<=0.5) h_ClustSat_HscpCorrWOinv_0eta05->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=1 && abs(eta)>0.5) h_ClustSat_HscpCorrWOinv_05eta10->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=1.5 && abs(eta)>1) h_ClustSat_HscpCorrWOinv_10eta15->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=2.5 && abs(eta)>1.5) h_ClustSat_HscpCorrWOinv_15eta25->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    if(abs(eta)<=2.5 && abs(eta)>2) h_ClustSat_HscpCorrWOinv_20eta25->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);

                    if(p<=500 && (EMyCorr_WithoutInversion-Eloss)/Eloss<-0.25)
                    {
                        h_ClustSat_LowP_ClusterSize->Fill(nstrips);
                        h_ClustSat_LowP_Sat254->Fill(nsat254);
                        h_ClustSat_LowP_Sat255->Fill(nsat255);
                        h_ClustSat_LowP_NSat->Fill(nsat254+nsat255);
                        h_ClustSat_LowP_layerlabel->Fill(layerLabel);
                        h_ClustSat_LowP_Eta->Fill(abs(eta));
                        h_ClustSat_LowP_Esim->Fill(Eloss);
                        h_ClustSat_LowP_Emycorr->Fill(EMyCorr_WithoutInversion);
                        h_ClustSat_LowP_EcorrHscp->Fill(ECorrHSCP_WithoutInversion);
                    }

                    if(subdetid==3)
                    {
                        h_ClustSat_MyCorrWOinv_TIB->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_TIB->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(subdetid==4)
                    {
                        h_ClustSat_MyCorrWOinv_TID->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_TID->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(subdetid==5)
                    {
                        h_ClustSat_MyCorrWOinv_TOB->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_TOB->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(subdetid==6)
                    {
                        h_ClustSat_MyCorrWOinv_TEC->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_TEC->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==1) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer1->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer1->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==2) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer2->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer2->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==3) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer3->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer3->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==4) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer4->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer4->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==5) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer5->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer5->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==6) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer6->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer6->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==7) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer7->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer7->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==8) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer8->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer8->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==9) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer9->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer9->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(layerLabel==10) 
                    {
                        h_ClustSat_MyCorrWOinv_Layer10->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_Layer10->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nstrips==1)
                    {
                        h_ClustSat_MyCorrWOinv_size1->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_size1->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nstrips==2)
                    {
                        h_ClustSat_MyCorrWOinv_size2->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_size2->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nstrips==3)
                    {
                        h_ClustSat_MyCorrWOinv_size3->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_size3->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nstrips==4)
                    {
                        h_ClustSat_MyCorrWOinv_size4->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_size4->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nstrips==5)
                    {
                        h_ClustSat_MyCorrWOinv_size5->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_size5->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nstrips==6)
                    {
                        h_ClustSat_MyCorrWOinv_size6->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_size6->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nsat254+nsat255==1)
                    {
                        h_ClustSat_MyCorrWOinv_sat1->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_sat1->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nsat254+nsat255==2)
                    {
                        h_ClustSat_MyCorrWOinv_sat2->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_sat2->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }
                    if(nsat254+nsat255==3)
                    {
                        h_ClustSat_MyCorrWOinv_sat3->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                        h_ClustSat_HscpCorrWOinv_sat3->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    }

                }


                if(!(sat254 || sat255))
                {

                    h_ClustSat_HscpCorrWinv_firstThreshold20->Fill((ECorrHSCPfirstThreshold20-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstThreshold50->Fill((ECorrHSCPfirstThreshold50-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstThreshold100->Fill((ECorrHSCPfirstThreshold100-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstThreshold150->Fill((ECorrHSCPfirstThreshold150-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstThreshold200->Fill((ECorrHSCPfirstThreshold200-Eloss)/Eloss);

                    h_ClustSat_HscpCorrWinv_secondThreshold25->Fill((ECorrHSCPsecondThreshold25-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondThreshold50->Fill((ECorrHSCPsecondThreshold50-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondThreshold100->Fill((ECorrHSCPsecondThreshold100-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondThreshold150->Fill((ECorrHSCPsecondThreshold150-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondThreshold200->Fill((ECorrHSCPsecondThreshold200-Eloss)/Eloss);

                    h_ClustSat_HscpCorrWinv_firstXTalk10->Fill((ECorrHSCPfirstXTalk10-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstXTalk11->Fill((ECorrHSCPfirstXTalk11-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstXTalk12->Fill((ECorrHSCPfirstXTalk12-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstXTalk9->Fill((ECorrHSCPfirstXTalk9-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstXTalk8->Fill((ECorrHSCPfirstXTalk8-Eloss)/Eloss);

                    h_ClustSat_HscpCorrWinv_secondXTalk4->Fill((ECorrHSCPsecondXTalk4-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondXTalk5->Fill((ECorrHSCPsecondXTalk5-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondXTalk6->Fill((ECorrHSCPsecondXTalk6-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondXTalk3->Fill((ECorrHSCPsecondXTalk3-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondXTalk2->Fill((ECorrHSCPsecondXTalk2-Eloss)/Eloss);


                    hCorrHSCP->Fill((ECorrHSCP-Eloss)/Eloss);
                    hMyCorr->Fill((EMyCorr-Eloss)/Eloss);
                    hMyCorr_RatioSat->Fill((EMyCorr_RatioSat-Eloss)/Eloss);

                    h_ClustNoSat_WithoutCorr->Fill((Erec-Eloss)/Eloss);
                    h_ClustNoSat_CorrHSCP_WithoutInversion->Fill((ECorrHSCP_WithoutInversion-Eloss)/Eloss);
                    h_ClustNoSat_CorrHSCP->Fill((ECorrHSCP-Eloss)/Eloss);            
                    h_ClustNoSat_MyCorr_WithoutInversion->Fill((EMyCorr_WithoutInversion-Eloss)/Eloss);
                    h_ClustNoSat_MyCorr->Fill((EMyCorr-Eloss)/Eloss);
                    h_ClustNoSat_MyCorr_WithoutInversion_RatioSat->Fill((EMyCorr_WithoutInversion_RatioSat-Eloss)/Eloss);
                    h_ClustNoSat_MyCorr_RatioSat->Fill((EMyCorr_RatioSat-Eloss)/Eloss);

                    h_ClustSat_HscpCorrWinv_firstThreshold20->Fill((ECorrHSCPfirstThreshold20-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstThreshold50->Fill((ECorrHSCPfirstThreshold50-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstThreshold100->Fill((ECorrHSCPfirstThreshold100-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstThreshold150->Fill((ECorrHSCPfirstThreshold150-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstThreshold200->Fill((ECorrHSCPfirstThreshold200-Eloss)/Eloss);

                    h_ClustSat_HscpCorrWinv_secondThreshold25->Fill((ECorrHSCPsecondThreshold25-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondThreshold50->Fill((ECorrHSCPsecondThreshold50-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondThreshold100->Fill((ECorrHSCPsecondThreshold100-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondThreshold150->Fill((ECorrHSCPsecondThreshold150-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondThreshold200->Fill((ECorrHSCPsecondThreshold200-Eloss)/Eloss);

                    h_ClustSat_HscpCorrWinv_firstXTalk10->Fill((ECorrHSCPfirstXTalk10-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstXTalk11->Fill((ECorrHSCPfirstXTalk11-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstXTalk12->Fill((ECorrHSCPfirstXTalk12-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstXTalk9->Fill((ECorrHSCPfirstXTalk9-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_firstXTalk8->Fill((ECorrHSCPfirstXTalk8-Eloss)/Eloss);

                    h_ClustSat_HscpCorrWinv_secondXTalk4->Fill((ECorrHSCPsecondXTalk4-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondXTalk5->Fill((ECorrHSCPsecondXTalk5-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondXTalk6->Fill((ECorrHSCPsecondXTalk6-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondXTalk3->Fill((ECorrHSCPsecondXTalk3-Eloss)/Eloss);
                    h_ClustSat_HscpCorrWinv_secondXTalk2->Fill((ECorrHSCPsecondXTalk2-Eloss)/Eloss);
                }



            }

            if(corr)
            {
                if(layerLabel<=22) dEsimdx.push_back(Eloss/pathlength);
                dErecdx.push_back(Erec/pathlength);
                dEcorrdx.push_back(EcorrBef/pathlength);
                if(nstrips>=3)
                {
                    if((RatioNClusterSat254+RatioNClusterSat255>=ratsat)) 
                    {
                        EcorrSatLayer = EcorrLayer.ChargeCorr(Erec,layerLabel,nstrips,nsat254,nsat255);
                        dEfulldx_Layer.push_back(EcorrSatLayer/pathlength);
                        EcorrSatModul = EcorrModulGeom.ChargeCorr(Erec,modulgeom,nstrips,nsat254,nsat255);
                        dEfulldx_ModulGeom.push_back(EcorrSatModul/pathlength);



                    }
                    else
                    {
                        if(!(sat254 || sat255))
                        {
                            dEfulldx_Layer.push_back(Erec/pathlength);
                            dEfulldx_ModulGeom.push_back(Erec/pathlength);
                        }
                    }
                }
            }
            

            if((sat254 || sat255))
            {
                hErecLayer->Fill((Erec-Eloss)/Eloss);
                //EcorrSatLayer = EcorrLayer.ChargeCorr(Erec,layerLabel,nstrips,nsat254,nsat255);
                hEcorrLayer->Fill((EcorrSatLayer-Eloss)/Eloss);
                h2ElossEcorrSatLayer->Fill(Eloss,EcorrSatLayer);
                hEcorrModulGeom->Fill((EcorrSatModul-Eloss)/Eloss);
                h2ElossEcorrModulGeom->Fill(Eloss,EcorrSatModul);

            }




			if(method)
			{
				//METHODE LAYER

				for(int countlayer=1;countlayer<22;countlayer++)
            	{
            	    for(int countnstrip=3;countnstrip<7;countnstrip++)
            	    {
            	        for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
            	        {
            	            for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
            	            {
            	                if(countlayer==layerLabel && countnstrip==nstrips && countnstripsat254==nsat254 && countnstripsat255==nsat255)
            	                {
            	                    if(countnstripsat254+countnstripsat255==1 && testCorrHscp) VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(charge_sim,charge);
                                    if(countnstripsat254+countnstripsat255>1) VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(charge_sim,charge);
                                    //if(countnstripsat254+countnstripsat255==1 && testCorrHscp) VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(Eloss,Erec);
                                    //if(countnstripsat254+countnstripsat255>1) VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(Eloss,Erec);
                                    //VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(Eloss,Erec);
            	                } 
            	            }
            	        }
            	    }
				}

				//METHODE MODULGEOM

				for(int countlayer=1;countlayer<15;countlayer++)
            	{
            	    for(int countnstrip=3;countnstrip<7;countnstrip++)
            	    {
            	        for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
            	        {
            	            for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
            	            {
            	                if(countlayer==modulgeom && countnstrip==nstrips && countnstripsat254==nsat254 && countnstripsat255==nsat255)
            	                {
            	                    //VectModulGeomVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(Eloss,Erec);
                                    if(countnstripsat254+countnstripsat255==1 && testCorrHscp) VectModulGeomVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(charge_sim,charge);
                                    if(countnstripsat254+countnstripsat255>1) VectModulGeomVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]->Fill(charge_sim,charge);
            	                } 
            	            }
            	        }
            	    }
				}

                if(Qminus>=25 && Qplus>=25 && RatioNClusterSat254+RatioNClusterSat255>=0.6 && DeltaQCalc<=40) sc.FillHisto(charge_sim,charge,modulgeom,nstrips,nsat254+nsat255);
				

			}
            




		}


        Estimator estimLayer(VectDecorrDxLayer);
        Estimator estimModulGeom(VectDecorrDxModulGeom);
        Estimator estimDelossDx(vect_dedx);

        h2PoverMDecorrDxLayer->Fill(p/2400,estimLayer.GetHarmonic2()*pow(10,3));
        h2PoverMDecorrDxModulGeom->Fill(p/2400,estimModulGeom.GetHarmonic2()*pow(10,3));

        h2PvsDecorrDxLayer->Fill(p,estimLayer.GetHarmonic2()*pow(10,3));
        h2PvsDecorrDxModulGeom->Fill(p,estimModulGeom.GetHarmonic2()*pow(10,3));

        float K = KsusyCorr;
        float C = CsusyCorr;
        //float dedx = estimLayer.GetHarmonic2()*pow(10,3);
        float dedx = estimDelossDx.GetHarmonic2()*pow(10,3);

        float mass_dedx = sqrt((dedx-C)*pow(p,2)/K);

        if(dedx>0) h2PvsDelossDxProtonKC->Fill(p,dedx);




        Estimator estimDecorrDx(vect_decorrdx);
        Estimator estimDecorrDxnewSatLayer(vect_decorrdxnewSatLayer);
        float decorrdx = estimDecorrDx.GetHarmonic2()*pow(10,3);
        float decorrdxnewsatlayer = estimDecorrDxnewSatLayer.GetHarmonic2()*pow(10,3);

        if(decorrdx>0) h2PvsDecorrDxProtonKC->Fill(p,decorrdx);
        if(decorrdxnewsatlayer>0) h2PvsDecorrDxProtonKCnewSatLayer->Fill(p,decorrdxnewsatlayer);





        float massEloss = sqrt((dedx-CprotonEloss)*pow(p,2)/KprotonEloss);
        float massEcorr = sqrt((decorrdx-CprotonEcorr)*pow(p,2)/KprotonEcorr);
        float massEcorrNew = sqrt((decorrdxnewsatlayer-CprotonEcorrNew)*pow(p,2)/KprotonEcorrNew);
        //float massEloss = sqrt((dedx-CprotonEloss)*pow(p,2)/KprotonEloss);
        //float massEcorr = sqrt((decorrdx-CprotonEcorr)*pow(p,2)/KprotonEcorr);
        //float massEcorrNew = sqrt((decorrdxnewsatlayer-CprotonEcorrNew)*pow(p,2)/KprotonEcorrNew);
        
        //if(dedx<=CprotonEloss && dedx>=0) massEloss=0.;
        //if(decorrdx<=CprotonEcorr && decorrdx>=0) massEcorr=0.;
        //if(decorrdxnewsatlayer<=CprotonEcorrNew && decorrdxnewsatlayer>=0) massEcorrNew=0.;
        distribmassprotonEloss->Fill(massEloss);
        distribmassprotonEcorr->Fill(massEcorr);
        distribmassprotonEcorrNew->Fill(massEcorrNew);



        float CsusySim = 4.53; //pour gluino
        float KsusySim = 1.76;




        Estimator estimDeDx(dEsimdx);
        pVsdEsimdx->Fill(p,estimDeDx.GetHarmonic2()*pow(10,3));
        //massEsim->Fill(sqrt((estimDeDx.GetHarmonic2()*pow(10,3)-C_Esim)*pow(p,2)/K_Esim));
        massEsim->Fill(sqrt((estimDeDx.GetHarmonic2()*pow(10,3)-CsusySim)*pow(p,2)/KsusySim));

        estimDeDx.SetVect(dErecdx);
        pVsdErecdx->Fill(p,estimDeDx.GetHarmonic2()*pow(10,3));
        massErec->Fill(sqrt((estimDeDx.GetHarmonic2()*pow(10,3)-C_Erec)*pow(p,2)/K_Erec));

        estimDeDx.SetVect(dEcorrdx);
        pVsdEcorrdx->Fill(p,estimDeDx.GetHarmonic2()*pow(10,3));
        hmassEcorr->Fill(sqrt((estimDeDx.GetHarmonic2()*pow(10,3)-C_Ecorr)*pow(p,2)/K_Ecorr));

        estimDeDx.SetVect(dEfulldx_Layer);
        pVsdEfulldx_Layer->Fill(p,estimDeDx.GetHarmonic2()*pow(10,3));
        massEfullLayer->Fill(sqrt((estimDeDx.GetHarmonic2()*pow(10,3)-C_EfullLayer)*pow(p,2)/K_EfullLayer));

        estimDeDx.SetVect(dEfulldx_ModulGeom);
        pVsdEfulldx_ModulGeom->Fill(p,estimDeDx.GetHarmonic2()*pow(10,3));
        massEfullModulGeom->Fill(sqrt((estimDeDx.GetHarmonic2()*pow(10,3)-C_EfullModulGeom)*pow(p,2)/K_EfullModulGeom));


        Estimator estimdEdxWTEC(dEdxWTEC);
        Estimator estimdEdxWOTEC(dEdxWOTEC);
        h1_IhWTEC->Fill(estimdEdxWTEC.GetHarmonic2()*pow(10,3));
        h1_IhWOTEC->Fill(estimdEdxWOTEC.GetHarmonic2()*pow(10,3));



        //Study Ih

        Estimator estimdEdxsim(dEdxsim);
        Estimator estimdEdxMyCorr(dEdxMyCorr);
        Estimator estimdEdxHscpCorrWOinv(dEdxHscpCorrWOinv);
        Estimator estimdEdxHscpCorr(dEdxHscpCorr);


        h2_IhHarmonic_sim->Fill(p/MGluino2400,estimdEdxsim.GetHarmonic2()*pow(10,3));
        h2_IhHarmonic_MyCorr->Fill(p/MGluino2400,estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
        h2_IhHarmonic_HscpCorrWOinv->Fill(p/MGluino2400,estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));
        h2_IhHarmonic_HscpCorr->Fill(p/MGluino2400,estimdEdxHscpCorr.GetHarmonic2()*pow(10,3));


        h2_ClustSat_IhSimIhMyCorrRatioSat->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
        h2_ClustSat_IhSimIhHscpCorr->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));
        prof_ClustSat_IhSimIhMyCorrRatioSat->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
        prof_ClustSat_IhSimIhHscpCorr->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));

        if(testSurCorr)
        {
            h2_ClustSatSurCorr_IhSimIhMyCorrRatioSat->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
            h2_ClustSatSurCorr_IhSimIhHscpCorr->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));
            prof_ClustSatSurCorr_IhSimIhMyCorrRatioSat->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
            prof_ClustSatSurCorr_IhSimIhHscpCorr->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));
        }

        if(!testSurCorr)
        {
            h2_ClustSatNoSurCorr_IhSimIhMyCorrRatioSat->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
            h2_ClustSatNoSurCorr_IhSimIhHscpCorr->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));
            prof_ClustSatNoSurCorr_IhSimIhMyCorrRatioSat->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
            prof_ClustSatNoSurCorr_IhSimIhHscpCorr->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));
        }

        if(testMoreThanOne)
        {
            h2_ClustSatMoreThanOneSat_IhSimIhMyCorrRatioSat->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
            h2_ClustSatMoreThanOneSat_IhSimIhHscpCorr->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));
            prof_ClustSatMoreThanOneSat_IhSimIhMyCorrRatioSat->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
            prof_ClustSatMoreThanOneSat_IhSimIhHscpCorr->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));
        }

        if(!testMoreThanOne)
        {
            h2_ClustSatOneSat_IhSimIhMyCorrRatioSat->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
            h2_ClustSatOneSat_IhSimIhHscpCorr->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));
            prof_ClustSatOneSat_IhSimIhMyCorrRatioSat->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
            prof_ClustSatOneSat_IhSimIhHscpCorr->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));
        }

        if(testHscptrack)
        {
            h2_ClustSatPassingHscpTest_IhSimIhMyCorrRatioSat->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
            h2_ClustSatPassingHscpTest_IhSimIhHscpCorr->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));
            prof_ClustSatPassingHscpTest_IhSimIhMyCorrRatioSat->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxMyCorr.GetHarmonic2()*pow(10,3));
            prof_ClustSatPassingHscpTest_IhSimIhHscpCorr->Fill(estimdEdxsim.GetHarmonic2()*pow(10,3),estimdEdxHscpCorrWOinv.GetHarmonic2()*pow(10,3));
        }


        //Ih
    
    {
        Estimator estimDeDxSim(DeDxSim);
        Estimator estimDeDxHscpCorr(DeDxHscpCorr);
        Estimator estimDeDxProposal(DeDxProposal);
        Estimator estimDeDxScdProposal(DeDxScdProposal);
        Estimator estimDeDxThirdProposal(DeDxThirdProposal);
        Estimator estimDeDxFourthProposal(DeDxFourthProposal);
        Estimator estimDeDxRecWOSat(DeDxRecWOSat);

        if(p/2400>0.2)IhPoverM->Fill(p/2400,estimDeDxSim.GetHarmonic2()*pow(10,3));
        float Csim=4.47,Ksim=1.77;
        float comparatif = DeDxMoyen(p,Ksim,Csim);
        comparatif=estimDeDxSim.GetHarmonic2()*pow(10,3);

        // if(DeDxHscpCorr.size()>=8)
        {
            h1_IhHscpCorr->Fill((estimDeDxHscpCorr.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);
            h2_IhHscpCorr->Fill(comparatif,estimDeDxHscpCorr.GetHarmonic2()*pow(10,3));
            h2_Resolution_IhHscpCorr->Fill(comparatif,(estimDeDxHscpCorr.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);
        }
        // if(DeDxProposal.size()>=8)
        {
            h1_IhProposal->Fill((estimDeDxProposal.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);
            h2_IhProposal->Fill(comparatif,estimDeDxProposal.GetHarmonic2()*pow(10,3));
            h2_Resolution_IhProposal->Fill(comparatif,(estimDeDxProposal.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);
        }

        h1_IhSim->Fill((estimDeDxSim.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);
        
        distribIhSim->Fill(estimDeDxSim.GetHarmonic2()*pow(10,3));


        h1_IhScdProposal->Fill((estimDeDxScdProposal.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);
        h1_IhThirdProposal->Fill((estimDeDxThirdProposal.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);
        h1_IhFourthProposal->Fill((estimDeDxFourthProposal.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);
        h1_IhRecWOSat->Fill((estimDeDxRecWOSat.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);


        h2_IhSim->Fill(comparatif,estimDeDxSim.GetHarmonic2()*pow(10,3));
        
        
        h2_IhScdProposal->Fill(comparatif,estimDeDxScdProposal.GetHarmonic2()*pow(10,3));
        h2_IhThirdProposal->Fill(comparatif,estimDeDxThirdProposal.GetHarmonic2()*pow(10,3));
        h2_IhFourthProposal->Fill(comparatif,estimDeDxFourthProposal.GetHarmonic2()*pow(10,3));
        h2_IhRecWOSat->Fill(comparatif,estimDeDxRecWOSat.GetHarmonic2()*pow(10,3));

        h2_Resolution_IhSim->Fill(comparatif,(estimDeDxSim.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);
        
        
        h2_Resolution_IhScdProposal->Fill(comparatif,(estimDeDxScdProposal.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);
        h2_Resolution_IhThirdProposal->Fill(comparatif,(estimDeDxThirdProposal.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);
        h2_Resolution_IhFourthProposal->Fill(comparatif,(estimDeDxFourthProposal.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);
        h2_Resolution_IhRecWOSat->Fill(comparatif,(estimDeDxRecWOSat.GetHarmonic2()*pow(10,3)-comparatif)/comparatif);


    }

    }

	}
} //END LOOP ENTRIES

//Ratio number of saturated strip


float RatioNSat = (float)((float)nSat/(float)ncluster);
float RatioSat1 = (float)((float)sat1/(float)nSat);
float RatioSat2 = (float)((float)sat2/(float)nSat);
float RatioSat3 = (float)((float)sat3/(float)nSat);
float RatioSat1HscpTest = (float)((float)sat1HscpTest/(float)nHscpTest);
float RatioNHscpTest = (float)((float)nHscpTest/(float)nSat);
float RatioNRatSatTest = (float)((float)nRatSatTest/(float)nSat);
float RatioNHscpAndRatSatTest = (float)((float)nHscpAndRatSatTest/(float)nHscpTest);
float RatioSat2RatSatTest = (float)((float)sat2RatSatTest/(float)(sat2+sat3));
float RatioSurCorr = (float)((float)surcorr/(float)nSat);
float RatioSurCorrHSCP = (float)((float)surcorrHSCP/(float)nSat);
float RatioNparticlePassingMyCrit = (float)((float)nPatriclePassingMyCrit/(float)nParticle);

cout<<"nsat   "<<RatioNSat<<endl;
cout<<"sat1   "<<RatioSat1<<endl;
cout<<"sat2   "<<RatioSat2<<endl;
cout<<"sat3   "<<RatioSat3<<endl;
cout<<"sat1HscpTest   "<<RatioSat1HscpTest<<endl;
cout<<"nHscpTest   "<<RatioNHscpTest<<endl;
cout<<"nRatSatTest   "<<RatioNRatSatTest<<endl;
cout<<"nHscpRatSatTest   "<<RatioNHscpAndRatSatTest<<endl;
cout<<"nsat2RatSatTest   "<<RatioSat2RatSatTest<<endl;
cout<<"nSurcorr   "<<RatioSurCorr<<endl;
cout<<"nSurcorrHscp   "<<RatioSurCorrHSCP<<endl;
cout<<"nParticlePassing   "<<RatioNparticlePassingMyCrit<<endl;

cout<<"sat    "<<nSat<<endl;
cout<<"counterpass    "<<counterpass<<endl;


float EffThresholdSat = (float)countThresholdSat/(float)countSatCluster;
float EffThresholdDiff = (float)countThresholdDiff/(float)countSatCluster;
float EffThresholdBoth = (float)countBothThreshold/(float)countSatCluster;

cout<<"EffSat    "<<EffThresholdSat<<endl;
cout<<"EffDiff    "<<EffThresholdDiff<<endl;
cout<<"EffBoth    "<<EffThresholdBoth<<endl;






// if(FactKC)
// {
//     vector<float> VectKCEsim = FactorKC(pVsdEsimdx,0.938);
//     cout<<VectKCEsim[0]<<" "<<VectKCEsim[1]<<endl;

//     vector<float> VectKCErec = FactorKC(pVsdErecdx,0.938);
//     cout<<VectKCErec[0]<<" "<<VectKCErec[1]<<endl;

//     vector<float> VectKCEcorr = FactorKC(pVsdEcorrdx,0.938);
//     cout<<VectKCEcorr[0]<<" "<<VectKCEcorr[1]<<endl;

//     vector<float> VectKCEfull_Layer = FactorKC(pVsdEfulldx_Layer,0.938);
//     cout<<VectKCEfull_Layer[0]<<" "<<VectKCEfull_Layer[1]<<endl;

//     vector<float> VectKCEfull_ModulGeom = FactorKC(pVsdEfulldx_ModulGeom,0.938);
//     cout<<VectKCEfull_ModulGeom[0]<<" "<<VectKCEfull_ModulGeom[1]<<endl;

//     ofstream outputfileFactorKC ("factKC.txt");
//     outputfileFactorKC<<VectKCEsim[0]<<"\t"<<VectKCEsim[1]<<endl;
//     outputfileFactorKC<<VectKCErec[0]<<"\t"<<VectKCErec[1]<<endl;
//     outputfileFactorKC<<VectKCEcorr[0]<<"\t"<<VectKCEcorr[1]<<endl;
//     outputfileFactorKC<<VectKCEfull_Layer[0]<<"\t"<<VectKCEfull_Layer[1]<<endl;
//     outputfileFactorKC<<VectKCEfull_ModulGeom[0]<<"\t"<<VectKCEfull_ModulGeom[1]<<endl;
// }

// vector<float> vKC = FactorKC(IhPoverM,2400);
// ofstream fileKC("KC.txt");
// fileKC<<"K  "<<vKC[0]<<"  +/-   "<<vKC[2]<<endl;
// fileKC<<"C  "<<vKC[1]<<"  +/-   "<<vKC[3]<<endl;
// fileKC.close();


//ECRITURE PARAMETRES METHODE


if(method)
{
    sc.SetFileAndTreeName("FileParameters.root","tree");
    sc.SetBranch();
    sc.WriteParameters();
    sc.WriteFile();
    

//METHODE LAYER 
/*
	Correction MethodLayer; 
    MethodLayer.SetRebin(10);
    //MethodLayer.SetRange(0.7);
	MethodLayer.SetFileAndTreeName("testMethodLayer.root","tree");
	MethodLayer.SetBranch();
	for(int countlayer=1;countlayer<22;countlayer++)
    {
        for(int countnstrip=3;countnstrip<7;countnstrip++)
        {
        	for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
            {
                for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
				{
                    //cout<<LabelLayer(countlayer)+" NStrip="+to_string(countnstrip)+" NStripSat254="+to_string(countnstripsat254)+" NStripSat255="+to_string(countnstripsat255)<<endl;
                    MethodLayer.SetHisto(*VectLayerVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]);
                    MethodLayer.FillHisto(countlayer,countnstrip,countnstripsat254,countnstripsat255);
                    MethodLayer.FillProfile();
                    MethodLayer.Write(countlayer,countnstrip,countnstripsat254,countnstripsat255);
                }
            }
        }
    }
	MethodLayer.WriteFile();
*/
//METHODE MODULGEOM 

	/*Correction MethodModulGeom; 
	MethodModulGeom.SetFileAndTreeName("testMethodModulGeom.root","tree");
	MethodModulGeom.SetBranch();
	for(int countlayer=1;countlayer<15;countlayer++)
    {
        for(int countnstrip=3;countnstrip<7;countnstrip++)
        {
            for(int countnstripsat254=0;countnstripsat254<countnstrip+1;countnstripsat254++)
            {
                for(int countnstripsat255=0;countnstripsat255<countnstrip-countnstripsat254+1;countnstripsat255++)
                {
                    if()
                    {
                        MethodModulGeom.SetHisto(*VectModulGeomVectNStripVectNStripSat254VectNStripSat255Histo[countlayer-1][countnstrip-3][countnstripsat254][countnstripsat255]);
                        MethodModulGeom.FillHisto(countlayer,countnstrip,countnstripsat254,countnstripsat255);
                        MethodModulGeom.FillProfile();
                        MethodModulGeom.Write(countlayer,countnstrip,countnstripsat254,countnstripsat255);
                    }
                }
            }
        }
    }
	MethodModulGeom.WriteFile();*/

}



//
/*if(XTalkInvStudy)
{
    for(int j=0;j<h2_ClustSat_MyCorrWOinv_p->GetYaxis()->GetNbins()+1;j++)
    {   
        int Xentries=h2_ClustSat_MyCorrWOinv_p->ProjectionX()->GetEntries();
        for(int i=0;i<h2_ClustSat_MyCorrWOinv_p->GetXaxis()->GetNbins()+1;i++)
        {
            h2_ClustSat_MyCorrWOinv_p->SetBinContent(i,j,(float)h2_ClustSat_MyCorrWOinv_p->GetBinContent(i,j)/Xentries);
        }
    }

    for(int j=0;j<h2_ClustSat_MyCorrWOinv_eta->GetYaxis()->GetNbins()+1;j++)
    {   
        int Xentries=h2_ClustSat_MyCorrWOinv_eta->ProjectionX()->GetEntries();
        for(int i=0;i<h2_ClustSat_MyCorrWOinv_eta->GetXaxis()->GetNbins()+1;i++)
        {
            h2_ClustSat_MyCorrWOinv_eta->SetBinContent(i,j,(float)h2_ClustSat_MyCorrWOinv_eta->GetBinContent(i,j)/Xentries);
        }
    }


ofstream ofile_Layer;
ofile_Layer.open("comparison_Layer.txt");
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer1);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer1);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer2);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer2);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer3);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer3);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer4);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer4);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer5);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer5);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer6);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer6);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer7);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer7);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer8);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer8);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer9);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer9);
StudyGraph(ofile_Layer,h_ClustSat_MyCorrWOinv_Layer10);
StudyGraph(ofile_Layer,h_ClustSat_HscpCorrWOinv_Layer10);
ofile_Layer.close();

ofstream ofile_ModulGeom;
ofile_ModulGeom.open("comparison_ModulGeom.txt");
ofile_ModulGeom<<"Name"<<"\t\t"<<"MPV"<<"\t"<<"DistribMean"<<"\t"<<"DistribMeanErr"<<"\t"<<"DistribStdDev"<<"\t"<<"DistribStdDevErr"<<"\t"<<"FitMean"<<"\t"<<"FitMeanErr"<<"\t"<<"FitSigma"<<"\t"<<"FitSigmaErr"<<"\t"<<"FracIntervalGlissantSigma"<<"\t"<<"FracIntervalGlissant68"<<"\t"<<"FracInterval"<<endl;
StudyGraph(ofile_ModulGeom,h_ClustSat_MyCorrWOinv_TIB);
StudyGraph(ofile_ModulGeom,h_ClustSat_HscpCorrWOinv_TIB);
StudyGraph(ofile_ModulGeom,h_ClustSat_MyCorrWOinv_TID);
StudyGraph(ofile_ModulGeom,h_ClustSat_HscpCorrWOinv_TID);
StudyGraph(ofile_ModulGeom,h_ClustSat_MyCorrWOinv_TOB);
StudyGraph(ofile_ModulGeom,h_ClustSat_HscpCorrWOinv_TOB);
StudyGraph(ofile_ModulGeom,h_ClustSat_MyCorrWOinv_TEC);
StudyGraph(ofile_ModulGeom,h_ClustSat_HscpCorrWOinv_TEC);
ofile_ModulGeom.close();

ofstream ofile_Size;
ofile_Size.open("comparison_Size.txt");
ofile_Size<<"Name"<<"\t\t"<<"MPV"<<"\t"<<"DistribMean"<<"\t"<<"DistribMeanErr"<<"\t"<<"DistribStdDev"<<"\t"<<"DistribStdDevErr"<<"\t"<<"FitMean"<<"\t"<<"FitMeanErr"<<"\t"<<"FitSigma"<<"\t"<<"FitSigmaErr"<<"\t"<<"FracIntervalGlissantSigma"<<"\t"<<"FracIntervalGlissant68"<<"\t"<<"FracInterval"<<endl;
StudyGraph(ofile_Size,h_ClustSat_MyCorrWOinv_size1);
StudyGraph(ofile_Size,h_ClustSat_HscpCorrWOinv_size1);
StudyGraph(ofile_Size,h_ClustSat_MyCorrWOinv_size2);
StudyGraph(ofile_Size,h_ClustSat_HscpCorrWOinv_size2);
StudyGraph(ofile_Size,h_ClustSat_MyCorrWOinv_size3);
StudyGraph(ofile_Size,h_ClustSat_HscpCorrWOinv_size3);
StudyGraph(ofile_Size,h_ClustSat_MyCorrWOinv_size4);
StudyGraph(ofile_Size,h_ClustSat_HscpCorrWOinv_size4);
StudyGraph(ofile_Size,h_ClustSat_MyCorrWOinv_size5);
StudyGraph(ofile_Size,h_ClustSat_HscpCorrWOinv_size5);
StudyGraph(ofile_Size,h_ClustSat_MyCorrWOinv_size6);
StudyGraph(ofile_Size,h_ClustSat_HscpCorrWOinv_size6);
ofile_Size.close();

ofstream ofile_NSat;
ofile_NSat.open("comparison_NSat.txt");
StudyGraph(ofile_NSat,h_ClustSat_MyCorrWOinv_sat1);
StudyGraph(ofile_NSat,h_ClustSat_HscpCorrWOinv_sat1);
StudyGraph(ofile_NSat,h_ClustSat_MyCorrWOinv_sat2);
StudyGraph(ofile_NSat,h_ClustSat_HscpCorrWOinv_sat2);
StudyGraph(ofile_NSat,h_ClustSat_MyCorrWOinv_sat3);
StudyGraph(ofile_NSat,h_ClustSat_HscpCorrWOinv_sat3);
ofile_NSat.close();

ofstream ofile_p;
ofile_p.open("comparison_p.txt");
ofile_p<<"Name"<<"\t\t"<<"MPV"<<"\t"<<"DistribMean"<<"\t"<<"DistribMeanErr"<<"\t"<<"DistribStdDev"<<"\t"<<"DistribStdDevErr"<<"\t"<<"FitMean"<<"\t"<<"FitMeanErr"<<"\t"<<"FitSigma"<<"\t"<<"FitSigmaErr"<<"\t"<<"FracIntervalGlissantSigma"<<"\t"<<"FracIntervalGlissant68"<<"\t"<<"FracInterval"<<endl;
StudyGraph(ofile_p,h_ClustSat_MyCorrWOinv_0p500);
StudyGraph(ofile_p,h_ClustSat_HscpCorrWOinv_0p500);
StudyGraph(ofile_p,h_ClustSat_MyCorrWOinv_500p1000);
StudyGraph(ofile_p,h_ClustSat_HscpCorrWOinv_500p1000);
StudyGraph(ofile_p,h_ClustSat_MyCorrWOinv_1000p1500);
StudyGraph(ofile_p,h_ClustSat_HscpCorrWOinv_1000p1500);
StudyGraph(ofile_p,h_ClustSat_MyCorrWOinv_1500p2000);
StudyGraph(ofile_p,h_ClustSat_HscpCorrWOinv_1500p2000);
StudyGraph(ofile_p,h_ClustSat_MyCorrWOinv_2000p4000);
StudyGraph(ofile_p,h_ClustSat_HscpCorrWOinv_2000p4000);
ofile_p.close();

ofstream ofile_eta;
ofile_eta.open("comparison_eta.txt");
StudyGraph(ofile_eta,h_ClustSat_MyCorrWOinv_0eta05);
StudyGraph(ofile_eta,h_ClustSat_HscpCorrWOinv_0eta05);
StudyGraph(ofile_eta,h_ClustSat_MyCorrWOinv_05eta10);
StudyGraph(ofile_eta,h_ClustSat_HscpCorrWOinv_05eta10);
StudyGraph(ofile_eta,h_ClustSat_MyCorrWOinv_10eta15);
StudyGraph(ofile_eta,h_ClustSat_HscpCorrWOinv_10eta15);
StudyGraph(ofile_eta,h_ClustSat_MyCorrWOinv_15eta25);
StudyGraph(ofile_eta,h_ClustSat_HscpCorrWOinv_15eta25);
ofile_eta.close();

ofstream ofile_firstThreshold;
ofile_firstThreshold.open("comparison_firstThreshold.txt");
StudyGraph(ofile_firstThreshold,h_ClustSat_HscpCorrWinv_firstThreshold20);
StudyGraph(ofile_firstThreshold,h_ClustSat_HscpCorrWinv_firstThreshold50);
StudyGraph(ofile_firstThreshold,h_ClustSat_HscpCorrWinv_firstThreshold100);
StudyGraph(ofile_firstThreshold,h_ClustSat_HscpCorrWinv_firstThreshold150);
StudyGraph(ofile_firstThreshold,h_ClustSat_HscpCorrWinv_firstThreshold200);
ofile_firstThreshold.close();

ofstream ofile_secondThreshold;
ofile_secondThreshold.open("comparison_secondThreshold.txt");
StudyGraph(ofile_secondThreshold,h_ClustSat_HscpCorrWinv_secondThreshold25);
StudyGraph(ofile_secondThreshold,h_ClustSat_HscpCorrWinv_secondThreshold50);
StudyGraph(ofile_secondThreshold,h_ClustSat_HscpCorrWinv_secondThreshold100);
StudyGraph(ofile_secondThreshold,h_ClustSat_HscpCorrWinv_secondThreshold150);
StudyGraph(ofile_secondThreshold,h_ClustSat_HscpCorrWinv_secondThreshold200);
ofile_secondThreshold.close();

ofstream ofile_firstXTalk;
ofile_firstXTalk.open("comparison_firstXTalk.txt");
StudyGraph(ofile_firstXTalk,h_ClustSat_HscpCorrWinv_firstXTalk10);
StudyGraph(ofile_firstXTalk,h_ClustSat_HscpCorrWinv_firstXTalk11);
StudyGraph(ofile_firstXTalk,h_ClustSat_HscpCorrWinv_firstXTalk12);
StudyGraph(ofile_firstXTalk,h_ClustSat_HscpCorrWinv_firstXTalk9);
StudyGraph(ofile_firstXTalk,h_ClustSat_HscpCorrWinv_firstXTalk8);
ofile_firstXTalk.close();

ofstream ofile_secondXTalk;
ofile_secondXTalk.open("comparison_secondXTalk.txt");
StudyGraph(ofile_secondXTalk,h_ClustSat_HscpCorrWinv_secondXTalk4);
StudyGraph(ofile_secondXTalk,h_ClustSat_HscpCorrWinv_secondXTalk5);
StudyGraph(ofile_secondXTalk,h_ClustSat_HscpCorrWinv_secondXTalk6);
StudyGraph(ofile_secondXTalk,h_ClustSat_HscpCorrWinv_secondXTalk3);
StudyGraph(ofile_secondXTalk,h_ClustSat_HscpCorrWinv_secondXTalk2);
ofile_secondXTalk.close();

ofstream ofile_Methods;
ofile_Methods.open("comparison_Methods.txt");
ofile_Methods<<"Name"<<"\t\t"<<"MPV"<<"\t"<<"DistribMean"<<"\t"<<"DistribMeanErr"<<"\t"<<"DistribStdDev"<<"\t"<<"DistribStdDevErr"<<"\t"<<"FitMean"<<"\t"<<"FitMeanErr"<<"\t"<<"FitSigma"<<"\t"<<"FitSigmaErr"<<"\t"<<"FracIntervalGlissantSigma"<<"\t"<<"FracIntervalGlissant68"<<"\t"<<"FracInterval"<<endl;
StudyGraph(ofile_Methods,h_ClustSat_HscpCorrNoInv);
StudyGraph(ofile_Methods,h_ClustSat_MyFullCorr);
ofile_Methods.close();

ofstream ofile_Methods_HscpTest;
ofile_Methods_HscpTest.open("comparison_Methods_HscpTest.txt");
ofile_Methods_HscpTest<<"Name"<<"\t\t"<<"MPV"<<"\t"<<"DistribMean"<<"\t"<<"DistribMeanErr"<<"\t"<<"DistribStdDev"<<"\t"<<"DistribStdDevErr"<<"\t"<<"FitMean"<<"\t"<<"FitMeanErr"<<"\t"<<"FitSigma"<<"\t"<<"FitSigmaErr"<<"\t"<<"FracIntervalGlissantSigma"<<"\t"<<"FracIntervalGlissant68"<<"\t"<<"FracInterval"<<endl;
StudyGraph(ofile_Methods_HscpTest,h_ClustSatAndHscpTest_HscpCorrNoInv);
StudyGraph(ofile_Methods_HscpTest,h_ClustSatAndHscpTest_MyFullCorr);
ofile_Methods_HscpTest.close();



}
*/ 

ofstream ofile_Criteria;
ofile_Criteria.open("comparison_criteria.txt");
ofile_Criteria<<"Name"<<"\t\t"<<"MPV"<<"\t"<<"DistribMean"<<"\t"<<"DistribMeanErr"<<"\t"<<"DistribStdDev"<<"\t"<<"DistribStdDevErr"<<"\t"<<"FitMean"<<"\t"<<"FitMeanErr"<<"\t"<<"FitSigma"<<"\t"<<"FitSigmaErr"<<"\t"<<"FracIntervalGlissantSigma"<<"\t"<<"FracIntervalGlissant68"<<"\t"<<"FracInterval"<<"\t"<<"FracRight"<<"\t"<<"FracLeft"<<endl;
StudyGraph(ofile_Criteria,h_ClustSat_HscpCorrNoInv);
StudyGraph(ofile_Criteria,h1_ClustSat_MyCorrModulGeom_WihtoutCrit);
StudyGraph(ofile_Criteria,h1_ClustSat_MyCorrModulGeom_DeltaQ);
StudyGraph(ofile_Criteria,h1_ClustSat_MyCorrModulGeom_ThresholdSat25);
StudyGraph(ofile_Criteria,h1_ClustSat_MyCorrModulGeom_RatioSat60);
StudyGraph(ofile_Criteria,h1_ClustSat_MyCorrModulGeom_RatioSat90);
StudyGraph(ofile_Criteria,h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat60);
StudyGraph(ofile_Criteria,h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat90);
StudyGraph(ofile_Criteria,h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat60);
StudyGraph(ofile_Criteria,h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat90);
StudyGraph(ofile_Criteria,h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQ);
StudyGraph(ofile_Criteria,h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60);
StudyGraph(ofile_Criteria,h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat90);
StudyGraph(ofile_Criteria,h1_ClustSat_MyCorrModulGeom_ChargeMinEstimated);
StudyGraph(ofile_Criteria,h1_NewProcedure);
ofile_Criteria.close();

ofstream ofile_EfficiencyCriteria;
ofile_EfficiencyCriteria.open("efficiency_criteria.txt");
ofile_EfficiencyCriteria<<"Only saturated clusters"<<endl;
ofile_EfficiencyCriteria<<"Passing no criteria \t\t"<<nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing HSCP criteria \t\t"<<nPassingHSCPtest/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing DeltaQ \t\t"<<nPassingDeltaQ/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing DeltaQ40 \t\t"<<nPassingDeltaQ40/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing DeltaQ30 \t\t"<<nPassingDeltaQ30/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing DeltaQ20 \t\t"<<nPassingDeltaQ20/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing DeltaQ10 \t\t"<<nPassingDeltaQ10/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing DeltaQ50 \t\t"<<nPassingDeltaQ50/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing DeltaQ60 \t\t"<<nPassingDeltaQ60/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing DeltaQ15 \t\t"<<nPassingDeltaQ15/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing DeltaQ5 \t\t"<<nPassingDeltaQ5/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing ThresholdSat25 \t\t"<<nPassingThresholdSat25/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing RatioSat60 \t\t"<<nPassingRatioSat60/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing RatioSat90 \t\t"<<nPassingRatioSat90/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing DeltaQRatioSat60 \t\t"<<nPassingDeltaQRatioSat60/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing DeltaQRatioSat90 \t\t"<<nPassingDeltaQRatioSat90/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing ThresholdSat25RatioSat60 \t\t"<<nPassingThresholdSat25RatioSat60/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing ThresholdSat25RatioSat90 \t\t"<<nPassingThresholdSat25RatioSat90/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing ThresholdSat25DeltaQRatioSat60 \t\t"<<nPassingThresholdSat25DeltaQRatioSat60/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing ThresholdSat25DeltaQRatioSat90 \t\t"<<nPassingThresholdSat25DeltaQRatioSat90/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing ChargeMinEstimated \t\t"<<nPassingChargeMinEstimated/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria<<"Passing NewProcedure \t\t"<<nPassingNewProcedure/nPassingNoCrit<<endl;
ofile_EfficiencyCriteria.close();

ofstream ofile_FitProfileCriteria;
ofile_FitProfileCriteria.open("FitProfileCrit.txt");
float p00=-1,p01=-1;

FitProfile(profile_ClustSat_HscpNoInv,p00,p01);
ofile_FitProfileCriteria<<"Passing HSCP criteria \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_HscpNoInv_CutEta,p00,p01);
ofile_FitProfileCriteria<<"Passing HSCP criteria CutEta \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_HscpNoInv_TIB1,p00,p01);
ofile_FitProfileCriteria<<"Passing HSCP criteria TIB1 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_HscpNoInv_OB1,p00,p01);
ofile_FitProfileCriteria<<"Passing HSCP criteria OB1 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_WihtoutCrit,p00,p01);
ofile_FitProfileCriteria<<"Passing no criteria \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_DeltaQ,p00,p01);
ofile_FitProfileCriteria<<"Passing DeltaQ \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_DeltaQ40,p00,p01);
ofile_FitProfileCriteria<<"Passing DeltaQ40 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_DeltaQ30,p00,p01);
ofile_FitProfileCriteria<<"Passing DeltaQ30 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_ThresholdSat25,p00,p01);
ofile_FitProfileCriteria<<"Passing ThresholdSat25 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_RatioSat60,p00,p01);
ofile_FitProfileCriteria<<"Passing RatioSat60 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_RatioSat90,p00,p01);
ofile_FitProfileCriteria<<"Passing RatioSat90 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_DeltaQRatioSat60,p00,p01);
ofile_FitProfileCriteria<<"Passing DeltaQRatioSat60 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_DeltaQRatioSat90,p00,p01);
ofile_FitProfileCriteria<<"Passing DeltaQRatioSat90 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat60,p00,p01);
ofile_FitProfileCriteria<<"Passing ThresholdSat25RatioSat60 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat90,p00,p01);
ofile_FitProfileCriteria<<"Passing ThresholdSat25RatioSat90 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60,p00,p01);
ofile_FitProfileCriteria<<"Passing ThresholdSat25DeltaQRatioSat60 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEta,p00,p01);
ofile_FitProfileCriteria<<"Passing ThresholdSat25DeltaQRatioSat60_CutEta \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TIB1,p00,p01);
ofile_FitProfileCriteria<<"Passing ThresholdSat25DeltaQRatioSat60_TIB1 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_OB1,p00,p01);
ofile_FitProfileCriteria<<"Passing ThresholdSat25DeltaQRatioSat60_OB1 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat90,p00,p01);
ofile_FitProfileCriteria<<"Passing ThresholdSat25DeltaQRatioSat90 \t\t"<<p00<<"\t"<<p01<<endl;

FitProfile(profile_NewProcedure,p00,p01);
ofile_FitProfileCriteria<<"New procedure \t\t"<<p00<<"\t"<<p01<<endl;

ofile_FitProfileCriteria<<"CorrelationFactor ThresholdSat25RatioSat60 \t\t"<<h2_ClustSat_Overlap->GetCorrelationFactor()<<endl;

ofile_FitProfileCriteria.close();


/*ofstream ofile_Nopassing;
ofile_Nopassing.open("NoPass.txt");
for(int i=0;i<nopass.size();i++)
{
    ofile_Nopassing<<nopass[i]<<endl;
}


for(int xbin=1;xbin<profile_Proposal->GetNbinsX()+1;xbin++)
{
    h1_ratioProfile->SetBinContent(xbin,profile_Proposal->GetBinContent(xbin)/profile_Proposal->GetBinCenter(xbin));
}
TFitResultPtr fitRatio = h1_ratioProfile->Fit("pol1","QRS","",400,4000);
ofstream ofile_calib;
ofile_calib.open("calib.txt");
ofile_calib<<fitRatio->Parameter(0)<<" "<<fitRatio->Parameter(1)<<endl;
ofile_calib.close();
*/ 

//ECRITURE RESULTATS ETUDES

TFile* outputfile = new TFile("testEtudeScdeMethod.root","RECREATE");

float Kcalc,Ccalc;
FactorKC(IhPoverM,Kcalc,Ccalc)->Write();
IhPoverM->Write();
cout<<"K  "<<Kcalc<<"  C   "<<Ccalc<<endl;

h1_IhSim->Write();
h1_IhHscpCorr->Write();
h1_IhProposal->Write();
h1_IhScdProposal->Write();
h1_IhThirdProposal->Write();
h1_IhFourthProposal->Write();
h1_IhRecWOSat->Write();

h2_IhSim->Write();
h2_IhHscpCorr->Write();
h2_IhProposal->Write();
h2_IhScdProposal->Write();
h2_IhThirdProposal->Write();
h2_IhFourthProposal->Write();
h2_IhRecWOSat->Write();

TH1F* profilerMPVsim = ProfilerMPV(h2_IhSim);
TH1F* profilerMPVhscp = ProfilerMPV(h2_IhHscpCorr);
TH1F* profilerMPVprop = ProfilerMPV(h2_IhProposal);
TH1F* profilerMPVscd = ProfilerMPV(h2_IhScdProposal);
TH1F* profilerMPVthird = ProfilerMPV(h2_IhThirdProposal);
TH1F* profilerMPVfourth = ProfilerMPV(h2_IhFourthProposal);
TH1F* profilerMPVrecwosat = ProfilerMPV(h2_IhRecWOSat);

h2_Resolution_IhSim->Write();
h2_Resolution_IhHscpCorr->Write();
h2_Resolution_IhProposal->Write();

TH1F* profilerResolutionMPVsim = ProfilerResolution(h2_Resolution_IhSim);
TH1F* profilerResolutionMPVhscp = ProfilerResolution(h2_Resolution_IhHscpCorr);
TH1F* profilerResolutionMPVprop = ProfilerResolution(h2_Resolution_IhProposal);
TH1F* profilerResolutionMPVscd = ProfilerResolution(h2_Resolution_IhScdProposal);
TH1F* profilerResolutionMPVthird = ProfilerResolution(h2_Resolution_IhThirdProposal);
TH1F* profilerResolutionMPVfourth = ProfilerResolution(h2_Resolution_IhFourthProposal);
TH1F* profilerResolutionMPVrecwosat = ProfilerResolution(h2_Resolution_IhRecWOSat);

vector<TH1F*> vProfileMPV;
vector<TH1F*> vProfileResolution;

vProfileMPV.push_back(profilerMPVhscp);
vProfileMPV.push_back(profilerMPVprop);
vProfileResolution.push_back(profilerResolutionMPVprop);
vProfileResolution.push_back(profilerResolutionMPVhscp);


DrawTransparency(vProfileMPV,distribIhSim,"linearity","I_{h}^{sim}","linearity")->Write();
DrawTransparency(vProfileResolution,distribIhSim,"resolution","I_{h}^{sim}","resolution")->Write();

profilerMPVsim->SetName("mpvsim");
profilerMPVhscp->SetName("mpvhscp");
profilerMPVprop->SetName("mpvprop");
profilerMPVscd->SetName("mpvscd");
profilerMPVthird->SetName("mpvthird");
profilerMPVfourth->SetName("mpvfourth");
profilerMPVrecwosat->SetName("mpvrecwosat");


profilerResolutionMPVsim->SetName("ressim");
profilerResolutionMPVhscp->SetName("reshscp");
profilerResolutionMPVprop->SetName("resprop");
profilerResolutionMPVscd->SetName("resscd");
profilerResolutionMPVthird->SetName("resthird");
profilerResolutionMPVfourth->SetName("resfourth");
profilerResolutionMPVrecwosat->SetName("resrecwosat");

profilerMPVsim->Write();
profilerMPVhscp->Write();
profilerMPVprop->Write();
profilerMPVscd->Write();
profilerMPVthird->Write();
profilerMPVfourth->Write();
profilerMPVrecwosat->Write();

profilerResolutionMPVsim->Write();
profilerResolutionMPVhscp->Write();
profilerResolutionMPVprop->Write();
profilerResolutionMPVscd->Write();
profilerResolutionMPVthird->Write();
profilerResolutionMPVfourth->Write();
profilerResolutionMPVrecwosat->Write();








TLine* line4000 = new TLine(0,0,4000,4000);
line4000->Write();

h1_NewProcedure->Scale(1./h1_NewProcedure->GetEntries());
h1_NewProcedure->Write();
profile_NewProcedure->Write();

h1_Proposal->Scale(1./h1_Proposal->GetEntries());
h1_Proposal->Write();
profile_Proposal->Write();
h1_ratioProfile->Write();

distribdedx->Write();
distribdedx2->Write();

distribdeltaQ->Write();
distribminQ->Write();

h1_IhWTEC->Write();
h1_IhWOTEC->Write();

h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEtaSup->Write();
h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TEC->Write();
h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TID->Write();
h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_Barrel->Write();

h1_ClustSat_HscpNoInv_CutEtaSup->Write();
h1_ClustSat_HscpNoInv_TEC->Write();
h1_ClustSat_HscpNoInv_TID->Write();
h1_ClustSat_HscpNoInv_Barrel->Write();

profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEtaSup->Write();
profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TEC->Write();
profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TID->Write();
profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_Barrel->Write();

profile_ClustSat_HscpNoInv_CutEtaSup->Write();
profile_ClustSat_HscpNoInv_TEC->Write();
profile_ClustSat_HscpNoInv_TID->Write();
profile_ClustSat_HscpNoInv_Barrel->Write();

h2_ClustSat_QsimBetaGamma_Gluino2400->Write();
h2_ClustSat_QsimBetaGamma_Gluino1400->Write();
h2_ClustSat_QsimBetaGamma_Stau651->Write();
prof_ClustSat_QsimBetaGamma_Gluino2400->Write();
prof_ClustSat_QsimBetaGamma_Gluino1400->Write();
prof_ClustSat_QsimBetaGamma_Stau651->Write();

h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEta->Write();
profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_CutEta->Write();
h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TIB1->Write();
profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_TIB1->Write();
h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_OB1->Write();
profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60_OB1->Write();

h1_ClustSat_HscpNoInv_CutEta->Write();
profile_ClustSat_HscpNoInv_CutEta->Write();
h1_ClustSat_HscpNoInv_TIB1->Write();
profile_ClustSat_HscpNoInv_TIB1->Write();
h1_ClustSat_HscpNoInv_OB1->Write();
profile_ClustSat_HscpNoInv_OB1->Write();





h1_ClustSat_Qsim->Write();
h2_ClustSat_Overlap->Write();
profile_ClustSat_Overlap->Write();

h1_ClustSatDeltaQ_nsat2_FakeClusters->Scale(1./h1_ClustSatDeltaQ_nsat2_FakeClusters->GetEntries());
h1_ClustSatDeltaQ_nsat2_FakeClusters->Write();
h1_ClustSatDeltaQ_nsat2_TrueClusters->Scale(1./h1_ClustSatDeltaQ_nsat2_TrueClusters->GetEntries());
h1_ClustSatDeltaQ_nsat2_TrueClusters->Write();

h_ClustSat_HscpCorrNoInv->Scale(1./h_ClustSat_HscpCorrNoInv->GetEntries());
h_ClustSat_HscpCorrNoInv->Write();
profile_ClustSat_HscpNoInv->Write();

h1_ClustSat_MyCorrModulGeom_WihtoutCrit->Scale(1./h1_ClustSat_MyCorrModulGeom_WihtoutCrit->GetEntries());
h1_ClustSat_MyCorrModulGeom_WihtoutCrit->Write();
profile_ClustSat_MyCorrModulGeom_WihtoutCrit->Write();

h1_ClustSat_MyCorrModulGeom_DeltaQ->Scale(1./h1_ClustSat_MyCorrModulGeom_DeltaQ->GetEntries());
h1_ClustSat_MyCorrModulGeom_DeltaQ->Write();
profile_ClustSat_MyCorrModulGeom_DeltaQ->Write();

h1_ClustSat_MyCorrModulGeom_DeltaQ40->Scale(1./h1_ClustSat_MyCorrModulGeom_DeltaQ40->GetEntries());
h1_ClustSat_MyCorrModulGeom_DeltaQ40->Write();
profile_ClustSat_MyCorrModulGeom_DeltaQ40->Write();
h2_ClustSat_MyCorrModulGeom_DeltaQ40->Write();
float p0=-1,p1=-1,Chi2=-1;
ProfileMPV(h2_ClustSat_MyCorrModulGeom_DeltaQ40,p0,p1,Chi2)->Write();
cout<<p0<<" "<<p1<<" "<<Chi2<<endl;

h1_ClustSat_MyCorrModulGeom_DeltaQ30->Scale(1./h1_ClustSat_MyCorrModulGeom_DeltaQ30->GetEntries());
h1_ClustSat_MyCorrModulGeom_DeltaQ30->Write();
profile_ClustSat_MyCorrModulGeom_DeltaQ30->Write();

h1_ClustSat_MyCorrModulGeom_DeltaQ50->Scale(1./h1_ClustSat_MyCorrModulGeom_DeltaQ50->GetEntries());
h1_ClustSat_MyCorrModulGeom_DeltaQ50->Write();
profile_ClustSat_MyCorrModulGeom_DeltaQ50->Write();

h1_ClustSat_MyCorrModulGeom_DeltaQ60->Scale(1./h1_ClustSat_MyCorrModulGeom_DeltaQ60->GetEntries());
h1_ClustSat_MyCorrModulGeom_DeltaQ60->Write();
profile_ClustSat_MyCorrModulGeom_DeltaQ60->Write();

h1_ClustSat_MyCorrModulGeom_DeltaQ20->Scale(1./h1_ClustSat_MyCorrModulGeom_DeltaQ20->GetEntries());
h1_ClustSat_MyCorrModulGeom_DeltaQ20->Write();
profile_ClustSat_MyCorrModulGeom_DeltaQ20->Write();

h1_ClustSat_MyCorrModulGeom_DeltaQ10->Scale(1./h1_ClustSat_MyCorrModulGeom_DeltaQ10->GetEntries());
h1_ClustSat_MyCorrModulGeom_DeltaQ10->Write();
profile_ClustSat_MyCorrModulGeom_DeltaQ10->Write();

h1_ClustSat_MyCorrModulGeom_DeltaQ15->Scale(1./h1_ClustSat_MyCorrModulGeom_DeltaQ15->GetEntries());
h1_ClustSat_MyCorrModulGeom_DeltaQ15->Write();
profile_ClustSat_MyCorrModulGeom_DeltaQ15->Write();

h1_ClustSat_MyCorrModulGeom_DeltaQ5->Scale(1./h1_ClustSat_MyCorrModulGeom_DeltaQ5->GetEntries());
h1_ClustSat_MyCorrModulGeom_DeltaQ5->Write();
profile_ClustSat_MyCorrModulGeom_DeltaQ5->Write();

h1_ClustSat_MyCorrModulGeom_ThresholdSat25->Scale(1./h1_ClustSat_MyCorrModulGeom_ThresholdSat25->GetEntries());
h1_ClustSat_MyCorrModulGeom_ThresholdSat25->Write();
profile_ClustSat_MyCorrModulGeom_ThresholdSat25->Write();

h1_ClustSat_MyCorrModulGeom_RatioSat60->Scale(1./h1_ClustSat_MyCorrModulGeom_RatioSat60->GetEntries());
h1_ClustSat_MyCorrModulGeom_RatioSat60->Write();
profile_ClustSat_MyCorrModulGeom_RatioSat60->Write();

h1_ClustSat_MyCorrModulGeom_RatioSat90->Scale(1./h1_ClustSat_MyCorrModulGeom_RatioSat90->GetEntries());
h1_ClustSat_MyCorrModulGeom_RatioSat90->Write();
profile_ClustSat_MyCorrModulGeom_RatioSat90->Write();

h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat60->Scale(1./h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat60->GetEntries());
h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat60->Write();
profile_ClustSat_MyCorrModulGeom_DeltaQRatioSat60->Write();

h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat90->Scale(1./h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat90->GetEntries());
h1_ClustSat_MyCorrModulGeom_DeltaQRatioSat90->Write();
profile_ClustSat_MyCorrModulGeom_DeltaQRatioSat90->Write();

h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat60->Scale(1./h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat60->GetEntries());
h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat60->Write();
profile_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat60->Write();

h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat90->Scale(1./h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat90->GetEntries());
h1_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat90->Write();
profile_ClustSat_MyCorrModulGeom_ThresholdSat25RatioSat90->Write();

h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQ->Scale(1./h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQ->GetEntries());
h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQ->Write();
profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQ->Write();

h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60->Scale(1./h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60->GetEntries());
h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60->Write();
profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat60->Write();

h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat90->Scale(1./h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat90->GetEntries());
h1_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat90->Write();
profile_ClustSat_MyCorrModulGeom_ThresholdSat25DeltaQRatioSat90->Write();

h1_ClustSat_MyCorrModulGeom_ChargeMinEstimated->Scale(1./h1_ClustSat_MyCorrModulGeom_ChargeMinEstimated->GetEntries());
h1_ClustSat_MyCorrModulGeom_ChargeMinEstimated->Write();
profile_ClustSat_MyCorrModulGeom_ChargeMinEstimated->Write();


/*for(int i=0;i<h2_ClustSat_ThresholdSatCrit->GetNbinsY()+1;i++)
{
    int count=0;
    for(int j=0;j<h2_ClustSat_ThresholdSatCrit->GetNbinsX()+1;j++) count+=h2_ClustSat_ThresholdSatCrit->GetBinContent(j,i);
    for(int j=0;j<h2_ClustSat_ThresholdSatCrit->GetNbinsX()+1;j++) h2_ClustSat_ThresholdSatCrit->SetBinContent(j,i,(float)h2_ClustSat_ThresholdSatCrit->GetBinContent(j,i)/(float)count);
}*/
h2_ClustSat_ThresholdSatCrit->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
h2_ClustSat_ThresholdSatCrit->GetYaxis()->SetTitle("Q^{#pm1} > y");
h2_ClustSat_ThresholdSatCrit->GetZaxis()->SetTitle("[u.a.]");

/*for(int i=0;i<h2_ClustSat_DeltaQ->GetNbinsY()+1;i++)
{
    int count=0;
    for(int j=0;j<h2_ClustSat_DeltaQ->GetNbinsX()+1;j++) count+=h2_ClustSat_DeltaQ->GetBinContent(j,i);
    for(int j=0;j<h2_ClustSat_DeltaQ->GetNbinsX()+1;j++) h2_ClustSat_DeltaQ->SetBinContent(j,i,(float)h2_ClustSat_DeltaQ->GetBinContent(j,i)/(float)count);
}*/
h2_ClustSat_DeltaQ->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
h2_ClustSat_DeltaQ->GetYaxis()->SetTitle("|Q^{+1}-Q^{-1}| < y");
h2_ClustSat_DeltaQ->GetZaxis()->SetTitle("[u.a.]");

h2_ClustSat_DeltaQ->Write();
h2_ClustSat_ThresholdSatCrit->Write();


h1_ClustSat_DeltaQ->Write();
h1_ClustSat_DeltaQ35->Write();
h1_ClustSat_DeltaQ30->Write();
h1_ClustSat_DeltaQ25->Write();
h1_ClustSat_DeltaQ20->Write();
h1_ClustSat_DeltaQ15->Write();




h1_DeltaQ_True->Write();
h1_DeltaQ_False->Write();



h1_FakeClusters_testThresholdDiff->Write();
h1_FakeClusters_testThresholdSat->Write();
h1_FakeClusters->Write();

h1_ClustSat_testThresholdSat_HscpCorr->Write();
h1_ClustSat_testThresholdDiff_HscpCorr->Write();
h1_ClustSat_testBothThreshold_HscpCorr->Write();


h_ClustSat_MyFullCorr->Write();
h_ClustSat_HscpCorrNoInv->Write();

h_ClustSatAndHscpTest_WithoutCorr->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
h_ClustSatAndHscpTest_MyFullCorr->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
h_ClustSatAndHscpTest_HscpCorrNoInv->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
h_ClustSatAndHscpTest_MyFullCorrModulGeom->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
h_ClustSatAndHscpTest_MyFullCorr->Write();
h_ClustSatAndHscpTest_HscpCorrNoInv->Write();
h_ClustSatAndHscpTest_WithoutCorr->Write();
h_ClustSatAndHscpTest_MyFullCorrModulGeom->Write();

h_ClustSat_Esim->Scale(1./h_ClustSat_Esim->Integral());
h_ClustSat_Erec->Scale(1./h_ClustSat_Erec->Integral());
h_ClustSat_ClusterSize->Scale(1./h_ClustSat_ClusterSize->Integral());
h_ClustSat_Sat254->Scale(1./h_ClustSat_Sat254->Integral());
h_ClustSat_Sat255->Scale(1./h_ClustSat_Sat255->Integral());
h_ClustSat_NSat->Scale(1./h_ClustSat_NSat->Integral());
h_ClustSat_layerlabel->Scale(1./h_ClustSat_layerlabel->Integral());
h_ClustSat_Eta->Scale(1./h_ClustSat_Eta->Integral());

h_ClustSat_SurCorr_Esim->Scale(1./h_ClustSat_SurCorr_Esim->Integral());
h_ClustSat_SurCorr_Erec->Scale(1./h_ClustSat_SurCorr_Erec->Integral());
h_ClustSat_SurCorr_ClusterSize->Scale(1./h_ClustSat_SurCorr_ClusterSize->Integral());
h_ClustSat_SurCorr_Sat254->Scale(1./h_ClustSat_SurCorr_Sat254->Integral());
h_ClustSat_SurCorr_Sat255->Scale(1./h_ClustSat_SurCorr_Sat255->Integral());
h_ClustSat_SurCorr_NSat->Scale(1./h_ClustSat_SurCorr_NSat->Integral());
h_ClustSat_SurCorr_layerlabel->Scale(1./h_ClustSat_SurCorr_layerlabel->Integral());
h_ClustSat_SurCorr_Eta->Scale(1./h_ClustSat_SurCorr_Eta->Integral());

h_ClustSat_Esim->Write();
h_ClustSat_Erec->Write();
h_ClustSat_ClusterSize->Write();
h_ClustSat_Sat254->Write();
h_ClustSat_Sat255->Write();
h_ClustSat_NSat->Write();
h_ClustSat_layerlabel->Write();
h_ClustSat_Eta->Write();

h_ClustSat_SurCorr_Esim->Write();
h_ClustSat_SurCorr_Erec->Write();
h_ClustSat_SurCorr_ClusterSize->Write();
h_ClustSat_SurCorr_Sat254->Write();
h_ClustSat_SurCorr_Sat255->Write();
h_ClustSat_SurCorr_NSat->Write();
h_ClustSat_SurCorr_layerlabel->Write();
h_ClustSat_SurCorr_Eta->Write();

h_ClustSat_SurCorr_Esim->Divide(h_ClustSat_Esim);
h_ClustSat_SurCorr_Esim->Write();

h_ClustSat_SurCorr_Erec->Divide(h_ClustSat_Erec);
h_ClustSat_SurCorr_Erec->Write();

h_ClustSat_SurCorr_layerlabel->Divide(h_ClustSat_layerlabel);
h_ClustSat_SurCorr_layerlabel->Write();

TCanvas* cEsimSurCorr = new TCanvas();
cEsimSurCorr->Divide(1,2);
cEsimSurCorr->cd(1);


h2_ClustSat_EsimEMyCorrRatioSat->GetXaxis()->SetTitle("E_{sim} [GeV]");
h2_ClustSat_EsimEMyCorrRatioSat->GetYaxis()->SetTitle("E_{corr} [GeV]");
h2_ClustSat_EsimECorrHscp->GetXaxis()->SetTitle("E_{sim} [GeV]");
h2_ClustSat_EsimECorrHscp->GetYaxis()->SetTitle("E_{corr} [GeV]");
h2_ClustSat_IhSimIhMyCorrRatioSat->GetXaxis()->SetTitle("I_{h}^{sim}");
h2_ClustSat_IhSimIhMyCorrRatioSat->GetYaxis()->SetTitle("I_{h}^{corr}");
h2_ClustSat_IhSimIhHscpCorr->GetXaxis()->SetTitle("I_{h}^{sim}");
h2_ClustSat_IhSimIhHscpCorr->GetYaxis()->SetTitle("I_{h}^{corr}");

h2_ClustSat_EsimEMyCorrRatioSat->Write();
h2_ClustSat_EsimECorrHscp->Write();
prof_ClustSat_EsimEMyCorrRatioSat->Write();
prof_ClustSat_EsimECorrHscp->Write();
h2_ClustSat_IhSimIhMyCorrRatioSat->Write();
h2_ClustSat_IhSimIhHscpCorr->Write();
prof_ClustSat_IhSimIhMyCorrRatioSat->Write();
prof_ClustSat_IhSimIhHscpCorr->Write();

h2_ClustSatSurCorr_EsimEMyCorrRatioSat->GetXaxis()->SetTitle("E_{sim} [GeV]");
h2_ClustSatSurCorr_EsimEMyCorrRatioSat->GetYaxis()->SetTitle("E_{corr} [GeV]");
h2_ClustSatSurCorr_EsimECorrHscp->GetXaxis()->SetTitle("E_{sim} [GeV]");
h2_ClustSatSurCorr_EsimECorrHscp->GetYaxis()->SetTitle("E_{corr} [GeV]");
h2_ClustSatSurCorr_IhSimIhMyCorrRatioSat->GetXaxis()->SetTitle("I_{h}^{sim}");
h2_ClustSatSurCorr_IhSimIhMyCorrRatioSat->GetYaxis()->SetTitle("I_{h}^{corr}");
h2_ClustSatSurCorr_IhSimIhHscpCorr->GetXaxis()->SetTitle("I_{h}^{sim}");
h2_ClustSatSurCorr_IhSimIhHscpCorr->GetYaxis()->SetTitle("I_{h}^{corr}");

h2_ClustSatSurCorr_EsimEMyCorrRatioSat->Write();
h2_ClustSatSurCorr_EsimECorrHscp->Write();
prof_ClustSatSurCorr_EsimEMyCorrRatioSat->Write();
prof_ClustSatSurCorr_EsimECorrHscp->Write();
h2_ClustSatSurCorr_IhSimIhMyCorrRatioSat->Write();
h2_ClustSatSurCorr_IhSimIhHscpCorr->Write();
prof_ClustSatSurCorr_IhSimIhMyCorrRatioSat->Write();
prof_ClustSatSurCorr_IhSimIhHscpCorr->Write();

h2_ClustSatNoSurCorr_EsimEMyCorrRatioSat->GetXaxis()->SetTitle("E_{sim} [GeV]");
h2_ClustSatNoSurCorr_EsimEMyCorrRatioSat->GetYaxis()->SetTitle("E_{corr} [GeV]");
h2_ClustSatNoSurCorr_EsimECorrHscp->GetXaxis()->SetTitle("E_{sim} [GeV]");
h2_ClustSatNoSurCorr_EsimECorrHscp->GetYaxis()->SetTitle("E_{corr} [GeV]");
h2_ClustSatNoSurCorr_IhSimIhMyCorrRatioSat->GetXaxis()->SetTitle("I_{h}^{sim}");
h2_ClustSatNoSurCorr_IhSimIhMyCorrRatioSat->GetYaxis()->SetTitle("I_{h}^{corr}");
h2_ClustSatNoSurCorr_IhSimIhHscpCorr->GetXaxis()->SetTitle("I_{h}^{sim}");
h2_ClustSatNoSurCorr_IhSimIhHscpCorr->GetYaxis()->SetTitle("I_{h}^{corr}");

h2_ClustSatNoSurCorr_EsimEMyCorrRatioSat->Write();
h2_ClustSatNoSurCorr_EsimECorrHscp->Write();
prof_ClustSatNoSurCorr_EsimEMyCorrRatioSat->Write();
prof_ClustSatNoSurCorr_EsimECorrHscp->Write();
h2_ClustSatNoSurCorr_IhSimIhMyCorrRatioSat->Write();
h2_ClustSatNoSurCorr_IhSimIhHscpCorr->Write();
prof_ClustSatNoSurCorr_IhSimIhMyCorrRatioSat->Write();
prof_ClustSatNoSurCorr_IhSimIhHscpCorr->Write();

h2_ClustSatMoreThanOneSat_EsimEMyCorrRatioSat->GetXaxis()->SetTitle("E_{sim} [GeV]");
h2_ClustSatMoreThanOneSat_EsimEMyCorrRatioSat->GetYaxis()->SetTitle("E_{corr} [GeV]");
h2_ClustSatMoreThanOneSat_EsimECorrHscp->GetXaxis()->SetTitle("E_{sim} [GeV]");
h2_ClustSatMoreThanOneSat_EsimECorrHscp->GetYaxis()->SetTitle("E_{corr} [GeV]");
h2_ClustSatMoreThanOneSat_IhSimIhMyCorrRatioSat->GetXaxis()->SetTitle("I_{h}^{sim}");
h2_ClustSatMoreThanOneSat_IhSimIhMyCorrRatioSat->GetYaxis()->SetTitle("I_{h}^{corr}");
h2_ClustSatMoreThanOneSat_IhSimIhHscpCorr->GetXaxis()->SetTitle("I_{h}^{sim}");
h2_ClustSatMoreThanOneSat_IhSimIhHscpCorr->GetYaxis()->SetTitle("I_{h}^{corr}");

h2_ClustSatMoreThanOneSat_EsimEMyCorrRatioSat->Write();
h2_ClustSatMoreThanOneSat_EsimECorrHscp->Write();
prof_ClustSatMoreThanOneSat_EsimEMyCorrRatioSat->Write();
prof_ClustSatMoreThanOneSat_EsimECorrHscp->Write();
h2_ClustSatMoreThanOneSat_IhSimIhMyCorrRatioSat->Write();
h2_ClustSatMoreThanOneSat_IhSimIhHscpCorr->Write();
prof_ClustSatMoreThanOneSat_IhSimIhMyCorrRatioSat->Write();
prof_ClustSatMoreThanOneSat_IhSimIhHscpCorr->Write();

h2_ClustSatOneSat_EsimEMyCorrRatioSat->GetXaxis()->SetTitle("E_{sim} [GeV]");
h2_ClustSatOneSat_EsimEMyCorrRatioSat->GetYaxis()->SetTitle("E_{corr} [GeV]");
h2_ClustSatOneSat_EsimECorrHscp->GetXaxis()->SetTitle("E_{sim} [GeV]");
h2_ClustSatOneSat_EsimECorrHscp->GetYaxis()->SetTitle("E_{corr} [GeV]");
h2_ClustSatOneSat_IhSimIhMyCorrRatioSat->GetXaxis()->SetTitle("I_{h}^{sim}");
h2_ClustSatOneSat_IhSimIhMyCorrRatioSat->GetYaxis()->SetTitle("I_{h}^{corr}");
h2_ClustSatOneSat_IhSimIhHscpCorr->GetXaxis()->SetTitle("I_{h}^{sim}");
h2_ClustSatOneSat_IhSimIhHscpCorr->GetYaxis()->SetTitle("I_{h}^{corr}");

h2_ClustSatOneSat_EsimEMyCorrRatioSat->Write();
h2_ClustSatOneSat_EsimECorrHscp->Write();
prof_ClustSatOneSat_EsimEMyCorrRatioSat->Write();
prof_ClustSatOneSat_EsimECorrHscp->Write();
h2_ClustSatOneSat_IhSimIhMyCorrRatioSat->Write();
h2_ClustSatOneSat_IhSimIhHscpCorr->Write();
prof_ClustSatOneSat_IhSimIhMyCorrRatioSat->Write();
prof_ClustSatOneSat_IhSimIhHscpCorr->Write();

h2_ClustSatPassingHscpTest_EsimEMyCorrRatioSat->GetXaxis()->SetTitle("E_{sim} [GeV]");
h2_ClustSatPassingHscpTest_EsimEMyCorrRatioSat->GetYaxis()->SetTitle("E_{corr} [GeV]");
h2_ClustSatPassingHscpTest_EsimECorrHscp->GetXaxis()->SetTitle("E_{sim} [GeV]");
h2_ClustSatPassingHscpTest_EsimECorrHscp->GetYaxis()->SetTitle("E_{corr} [GeV]");
h2_ClustSatPassingHscpTest_IhSimIhMyCorrRatioSat->GetXaxis()->SetTitle("I_{h}^{sim}");
h2_ClustSatPassingHscpTest_IhSimIhMyCorrRatioSat->GetYaxis()->SetTitle("I_{h}^{corr}");
h2_ClustSatPassingHscpTest_IhSimIhHscpCorr->GetXaxis()->SetTitle("I_{h}^{sim}");
h2_ClustSatPassingHscpTest_IhSimIhHscpCorr->GetYaxis()->SetTitle("I_{h}^{corr}");

h2_ClustSatPassingHscpTest_EsimEMyCorrRatioSat->Write();
h2_ClustSatPassingHscpTest_EsimECorrHscp->Write();
prof_ClustSatPassingHscpTest_EsimEMyCorrRatioSat->Write();
prof_ClustSatPassingHscpTest_EsimECorrHscp->Write();
h2_ClustSatPassingHscpTest_IhSimIhMyCorrRatioSat->Write();
h2_ClustSatPassingHscpTest_IhSimIhHscpCorr->Write();
prof_ClustSatPassingHscpTest_IhSimIhMyCorrRatioSat->Write();
prof_ClustSatPassingHscpTest_IhSimIhHscpCorr->Write();





h_ClustSat_MyCorrRatioSat0->Write();
h_ClustSat_MyCorrRatioSat10->Write();
h_ClustSat_MyCorrRatioSat20->Write();
h_ClustSat_MyCorrRatioSat30->Write();
h_ClustSat_MyCorrRatioSat40->Write();
h_ClustSat_MyCorrRatioSat50->Write();
h_ClustSat_MyCorrRatioSat60->Write();
h_ClustSat_MyCorrRatioSat70->Write();
h_ClustSat_MyCorrRatioSat80->Write();
h_ClustSat_MyCorrRatioSat90->Write();
h_ClustSat_MyCorrRatioSat100->Write();



h2_IhHarmonic_sim->Write();
h2_IhHarmonic_MyCorr->Write();
h2_IhHarmonic_HscpCorrWOinv->Write();
h2_IhHarmonic_HscpCorr->Write();

h_ClustSat_HscpCorrWinv_firstThreshold20->Write();
h_ClustSat_HscpCorrWinv_firstThreshold50->Write();
h_ClustSat_HscpCorrWinv_firstThreshold100->Write();
h_ClustSat_HscpCorrWinv_firstThreshold150->Write();
h_ClustSat_HscpCorrWinv_firstThreshold200->Write();

h_ClustSat_HscpCorrWinv_secondThreshold25->Write();
h_ClustSat_HscpCorrWinv_secondThreshold50->Write();
h_ClustSat_HscpCorrWinv_secondThreshold100->Write();
h_ClustSat_HscpCorrWinv_secondThreshold150->Write();
h_ClustSat_HscpCorrWinv_secondThreshold200->Write();

h_ClustSat_HscpCorrWinv_firstXTalk10->Write();
h_ClustSat_HscpCorrWinv_firstXTalk11->Write();
h_ClustSat_HscpCorrWinv_firstXTalk12->Write();
h_ClustSat_HscpCorrWinv_firstXTalk9->Write();
h_ClustSat_HscpCorrWinv_firstXTalk8->Write();

h_ClustSat_HscpCorrWinv_secondXTalk4->Write();
h_ClustSat_HscpCorrWinv_secondXTalk5->Write();
h_ClustSat_HscpCorrWinv_secondXTalk6->Write();
h_ClustSat_HscpCorrWinv_secondXTalk3->Write();
h_ClustSat_HscpCorrWinv_secondXTalk2->Write();


h_ClustSat_LowP_ClusterSize->Write();
h_ClustSat_LowP_Sat254->Write();
h_ClustSat_LowP_Sat255->Write();
h_ClustSat_LowP_NSat->Write();
h_ClustSat_LowP_layerlabel->Write();
h_ClustSat_LowP_Eta->Write();
h_ClustSat_LowP_Esim->Write();
h_ClustSat_LowP_Emycorr->Write();
h_ClustSat_LowP_EcorrHscp->Write();

h_ClustSat_MyCorrWOinv_TIB->Write();
h_ClustSat_MyCorrWOinv_TID->Write();
h_ClustSat_MyCorrWOinv_TOB->Write();
h_ClustSat_MyCorrWOinv_TEC->Write();

h_ClustSat_HscpCorrWOinv_TIB->Write();
h_ClustSat_HscpCorrWOinv_TID->Write();
h_ClustSat_HscpCorrWOinv_TOB->Write();
h_ClustSat_HscpCorrWOinv_TEC->Write();

h_ClustSat_MyCorrWOinv_Layer1->Write();
h_ClustSat_MyCorrWOinv_Layer2->Write();
h_ClustSat_MyCorrWOinv_Layer3->Write();
h_ClustSat_MyCorrWOinv_Layer4->Write();
h_ClustSat_MyCorrWOinv_Layer5->Write();
h_ClustSat_MyCorrWOinv_Layer6->Write();
h_ClustSat_MyCorrWOinv_Layer7->Write();
h_ClustSat_MyCorrWOinv_Layer8->Write();
h_ClustSat_MyCorrWOinv_Layer9->Write();
h_ClustSat_MyCorrWOinv_Layer10->Write();

h_ClustSat_HscpCorrWOinv_Layer1->Write();
h_ClustSat_HscpCorrWOinv_Layer2->Write();
h_ClustSat_HscpCorrWOinv_Layer3->Write();
h_ClustSat_HscpCorrWOinv_Layer4->Write();
h_ClustSat_HscpCorrWOinv_Layer5->Write();
h_ClustSat_HscpCorrWOinv_Layer6->Write();
h_ClustSat_HscpCorrWOinv_Layer7->Write();
h_ClustSat_HscpCorrWOinv_Layer8->Write();
h_ClustSat_HscpCorrWOinv_Layer9->Write();
h_ClustSat_HscpCorrWOinv_Layer10->Write();

h_ClustSat_MyCorrWOinv_size1->Write();
h_ClustSat_MyCorrWOinv_size2->Write();
h_ClustSat_MyCorrWOinv_size3->Write();
h_ClustSat_MyCorrWOinv_size4->Write();
h_ClustSat_MyCorrWOinv_size5->Write();
h_ClustSat_MyCorrWOinv_size6->Write();
h_ClustSat_MyCorrWOinv_size7->Write();
h_ClustSat_MyCorrWOinv_size8->Write();

h_ClustSat_HscpCorrWOinv_size1->Write();
h_ClustSat_HscpCorrWOinv_size2->Write();
h_ClustSat_HscpCorrWOinv_size3->Write();
h_ClustSat_HscpCorrWOinv_size4->Write();
h_ClustSat_HscpCorrWOinv_size5->Write();
h_ClustSat_HscpCorrWOinv_size6->Write();
h_ClustSat_HscpCorrWOinv_size7->Write();
h_ClustSat_HscpCorrWOinv_size8->Write();

h_ClustSat_MyCorrWOinv_sat2->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");

h_ClustSat_MyCorrWOinv_sat1->Write();
h_ClustSat_MyCorrWOinv_sat2->Write();
h_ClustSat_MyCorrWOinv_sat3->Write();

h_ClustSat_HscpCorrWOinv_sat1->Write();
h_ClustSat_HscpCorrWOinv_sat2->Write();
h_ClustSat_HscpCorrWOinv_sat3->Write();


h_ClustSat_MyCorrWOinv_0p500->Write();
h_ClustSat_MyCorrWOinv_500p1000->Write();
h_ClustSat_MyCorrWOinv_1000p1500->Write();
h_ClustSat_MyCorrWOinv_1500p2000->Write();
h_ClustSat_MyCorrWOinv_2000p4000->Write();
h_ClustSat_MyCorrWOinv_0eta05->Write();
h_ClustSat_MyCorrWOinv_05eta10->Write();
h_ClustSat_MyCorrWOinv_10eta15->Write();
h_ClustSat_MyCorrWOinv_15eta25->Write();
h_ClustSat_MyCorrWOinv_20eta25->Write();

h_ClustSat_HscpCorrWOinv_0p500->Write();
h_ClustSat_HscpCorrWOinv_500p1000->Write();
h_ClustSat_HscpCorrWOinv_1000p1500->Write();
h_ClustSat_HscpCorrWOinv_1500p2000->Write();
h_ClustSat_HscpCorrWOinv_2000p4000->Write();
h_ClustSat_HscpCorrWOinv_0eta05->Write();
h_ClustSat_HscpCorrWOinv_05eta10->Write();
h_ClustSat_HscpCorrWOinv_10eta15->Write();
h_ClustSat_HscpCorrWOinv_15eta25->Write();
h_ClustSat_HscpCorrWOinv_20eta25->Write();


TCanvas* cMyCorrDiffP = new TCanvas();
TLegend* legendMyCorrP = new TLegend(0.7,0.7,0.9,0.9);
legendMyCorrP->AddEntry(h_ClustSat_MyCorrWOinv_0p500,"0<p<500","l");
legendMyCorrP->AddEntry(h_ClustSat_MyCorrWOinv_500p1000,"500<p<1000","l");
legendMyCorrP->AddEntry(h_ClustSat_MyCorrWOinv_1000p1500,"1000<p<1500","l");
legendMyCorrP->AddEntry(h_ClustSat_MyCorrWOinv_1500p2000,"1500<p<2000","l");
legendMyCorrP->AddEntry(h_ClustSat_MyCorrWOinv_2000p4000,"2000<p<4000","l");
h_ClustSat_MyCorrWOinv_0p500->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_500p1000->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_1000p1500->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_1500p2000->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_2000p4000->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_0p500->SetLineColor(30);
h_ClustSat_MyCorrWOinv_500p1000->SetLineColor(38);
h_ClustSat_MyCorrWOinv_1000p1500->SetLineColor(42);
h_ClustSat_MyCorrWOinv_1500p2000->SetLineColor(46);
h_ClustSat_MyCorrWOinv_2000p4000->SetLineColor(1);
h_ClustSat_MyCorrWOinv_0p500->Scale(1./h_ClustSat_MyCorrWOinv_0p500->GetEntries());
h_ClustSat_MyCorrWOinv_500p1000->Scale(1./h_ClustSat_MyCorrWOinv_500p1000->GetEntries());
h_ClustSat_MyCorrWOinv_1000p1500->Scale(1./h_ClustSat_MyCorrWOinv_1000p1500->GetEntries());
h_ClustSat_MyCorrWOinv_1500p2000->Scale(1./h_ClustSat_MyCorrWOinv_1500p2000->GetEntries());
h_ClustSat_MyCorrWOinv_2000p4000->Scale(1./h_ClustSat_MyCorrWOinv_2000p4000->GetEntries());
h_ClustSat_MyCorrWOinv_0p500->Draw();
h_ClustSat_MyCorrWOinv_500p1000->Draw("same");
h_ClustSat_MyCorrWOinv_1000p1500->Draw("same");
h_ClustSat_MyCorrWOinv_1500p2000->Draw("same");
h_ClustSat_MyCorrWOinv_2000p4000->Draw("same");
legendMyCorrP->Draw("same");
cMyCorrDiffP->Write();

TCanvas* cMyCorrDiffeta = new TCanvas();
TLegend* legendMyCorreta = new TLegend(0.7,0.7,0.9,0.9);
legendMyCorreta->AddEntry(h_ClustSat_MyCorrWOinv_0eta05,"|#eta|<0.5","l");
legendMyCorreta->AddEntry(h_ClustSat_MyCorrWOinv_05eta10,"0.5<|#eta|<1.0","l");
legendMyCorreta->AddEntry(h_ClustSat_MyCorrWOinv_10eta15,"1.0<|#eta|<1.5","l");
legendMyCorreta->AddEntry(h_ClustSat_MyCorrWOinv_15eta25,"1.5<|#eta|<2.5","l");
//legendMyCorreta->AddEntry(h_ClustSat_MyCorrWOinv_20eta25,"2.0<|#eta|<2.5","l");
h_ClustSat_MyCorrWOinv_0eta05->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_05eta10->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_10eta15->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_15eta25->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_20eta25->SetLineWidth(2);
h_ClustSat_MyCorrWOinv_0eta05->SetLineColor(30);
h_ClustSat_MyCorrWOinv_05eta10->SetLineColor(38);
h_ClustSat_MyCorrWOinv_10eta15->SetLineColor(42);
h_ClustSat_MyCorrWOinv_15eta25->SetLineColor(46);
h_ClustSat_MyCorrWOinv_20eta25->SetLineColor(1);
h_ClustSat_MyCorrWOinv_0eta05->Scale(1./h_ClustSat_MyCorrWOinv_0eta05->GetEntries());
h_ClustSat_MyCorrWOinv_05eta10->Scale(1./h_ClustSat_MyCorrWOinv_05eta10->GetEntries());
h_ClustSat_MyCorrWOinv_10eta15->Scale(1./h_ClustSat_MyCorrWOinv_10eta15->GetEntries());
h_ClustSat_MyCorrWOinv_15eta25->Scale(1./h_ClustSat_MyCorrWOinv_15eta25->GetEntries());
h_ClustSat_MyCorrWOinv_20eta25->Scale(1./h_ClustSat_MyCorrWOinv_20eta25->GetEntries());
h_ClustSat_MyCorrWOinv_0eta05->Draw();
h_ClustSat_MyCorrWOinv_05eta10->Draw("same");
h_ClustSat_MyCorrWOinv_10eta15->Draw("same");
h_ClustSat_MyCorrWOinv_15eta25->Draw("same");
//h_ClustSat_MyCorrWOinv_20eta25->Draw("same");
legendMyCorreta->Draw("same");
cMyCorrDiffeta->Write();

TCanvas* cHscpCorrDiffP = new TCanvas();
TLegend* legendHscpCorrP = new TLegend(0.7,0.7,0.9,0.9);
legendHscpCorrP->AddEntry(h_ClustSat_HscpCorrWOinv_0p500,"0<p<500","l");
legendHscpCorrP->AddEntry(h_ClustSat_HscpCorrWOinv_500p1000,"500<p<1000","l");
legendHscpCorrP->AddEntry(h_ClustSat_HscpCorrWOinv_1000p1500,"1000<p<1500","l");
legendHscpCorrP->AddEntry(h_ClustSat_HscpCorrWOinv_1500p2000,"1500<p<2000","l");
legendHscpCorrP->AddEntry(h_ClustSat_HscpCorrWOinv_2000p4000,"2000<p<4000","l");
h_ClustSat_HscpCorrWOinv_0p500->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_500p1000->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_1000p1500->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_1500p2000->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_2000p4000->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_0p500->SetLineColor(30);
h_ClustSat_HscpCorrWOinv_500p1000->SetLineColor(38);
h_ClustSat_HscpCorrWOinv_1000p1500->SetLineColor(42);
h_ClustSat_HscpCorrWOinv_1500p2000->SetLineColor(46);
h_ClustSat_HscpCorrWOinv_2000p4000->SetLineColor(1);
h_ClustSat_HscpCorrWOinv_0p500->Scale(1./h_ClustSat_HscpCorrWOinv_0p500->GetEntries());
h_ClustSat_HscpCorrWOinv_500p1000->Scale(1./h_ClustSat_HscpCorrWOinv_500p1000->GetEntries());
h_ClustSat_HscpCorrWOinv_1000p1500->Scale(1./h_ClustSat_HscpCorrWOinv_1000p1500->GetEntries());
h_ClustSat_HscpCorrWOinv_1500p2000->Scale(1./h_ClustSat_HscpCorrWOinv_1500p2000->GetEntries());
h_ClustSat_HscpCorrWOinv_2000p4000->Scale(1./h_ClustSat_HscpCorrWOinv_2000p4000->GetEntries());
h_ClustSat_HscpCorrWOinv_0p500->Draw();
h_ClustSat_HscpCorrWOinv_500p1000->Draw("same");
h_ClustSat_HscpCorrWOinv_1000p1500->Draw("same");
h_ClustSat_HscpCorrWOinv_1500p2000->Draw("same");
h_ClustSat_HscpCorrWOinv_2000p4000->Draw("same");
legendHscpCorrP->Draw("same");
cHscpCorrDiffP->Write();

TCanvas* cHscpCorrDiffeta = new TCanvas();
TLegend* legendHscpCorreta = new TLegend(0.7,0.7,0.9,0.9);
legendHscpCorreta->AddEntry(h_ClustSat_HscpCorrWOinv_0eta05,"|#eta|<0.5","l");
legendHscpCorreta->AddEntry(h_ClustSat_HscpCorrWOinv_05eta10,"0.5<|#eta|<1.0","l");
legendHscpCorreta->AddEntry(h_ClustSat_HscpCorrWOinv_10eta15,"1.0<|#eta|<1.5","l");
legendHscpCorreta->AddEntry(h_ClustSat_HscpCorrWOinv_15eta25,"1.5<|#eta|<2.5","l");
//legendHscpCorreta->AddEntry(h_ClustSat_HscpCorrWOinv_20eta25,"2.0<|#eta|<2.5","l");
h_ClustSat_HscpCorrWOinv_0eta05->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_05eta10->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_10eta15->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_15eta25->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_20eta25->SetLineWidth(2);
h_ClustSat_HscpCorrWOinv_0eta05->SetLineColor(30);
h_ClustSat_HscpCorrWOinv_05eta10->SetLineColor(38);
h_ClustSat_HscpCorrWOinv_10eta15->SetLineColor(42);
h_ClustSat_HscpCorrWOinv_15eta25->SetLineColor(46);
h_ClustSat_HscpCorrWOinv_20eta25->SetLineColor(1);
h_ClustSat_HscpCorrWOinv_0eta05->Scale(1./h_ClustSat_HscpCorrWOinv_0eta05->GetEntries());
h_ClustSat_HscpCorrWOinv_05eta10->Scale(1./h_ClustSat_HscpCorrWOinv_05eta10->GetEntries());
h_ClustSat_HscpCorrWOinv_10eta15->Scale(1./h_ClustSat_HscpCorrWOinv_10eta15->GetEntries());
h_ClustSat_HscpCorrWOinv_15eta25->Scale(1./h_ClustSat_HscpCorrWOinv_15eta25->GetEntries());
h_ClustSat_HscpCorrWOinv_20eta25->Scale(1./h_ClustSat_HscpCorrWOinv_20eta25->GetEntries());
h_ClustSat_HscpCorrWOinv_0eta05->Draw();
h_ClustSat_HscpCorrWOinv_05eta10->Draw("same");
h_ClustSat_HscpCorrWOinv_10eta15->Draw("same");
h_ClustSat_HscpCorrWOinv_15eta25->Draw("same");
//h_ClustSat_HscpCorrWOinv_20eta25->Draw("same");
legendHscpCorreta->Draw("same");
cHscpCorrDiffeta->Write();

/*h_ClustSat_DistribP->Scale(1/h_ClustSat_DistribP->Integral());
h_ClustSat_DistribEta->Scale(1/h_ClustSat_DistribEta->Integral());
h_ClustSat_DistribSize->Scale((float)(1./(float)h_ClustSat_DistribSize->GetEntries()));
h_ClustSat_DistribNSat->Scale((float)(1./(float)h_ClustSat_DistribNSat->GetEntries()));

h_DistribP->Scale(1/h_DistribP->Integral());
h_DistribEta->Scale(1/h_DistribEta->Integral());
h_DistribSize->Scale((float)(1./(float)h_DistribSize->GetEntries()));
h_DistribNSat->Scale((float)(1./(float)h_DistribNSat->GetEntries()));*/

h_ClustSat_DistribP->Write();
h_ClustSat_DistribEta->Write();
h_ClustSat_DistribSize->Write();
h_ClustSat_DistribNSat->Write();

h_DistribP->Write();
h_DistribEta->Write();
h_DistribSize->Write();
h_DistribNSat->Write();

h2_ClustSat_MyCorrWOinv_p->Write();
h2_ClustSat_MyCorrWOinv_eta->Write();

hWithoutCorr->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
hWithoutCorr->Write();
hCorrHSCP_WithoutInversion->Write();
hCorrHSCP->Write();
hMyCorr_WithoutInversion->Write();
hMyCorr->Write();
hMyCorr_WithoutInversion_RatioSat->Write();
hMyCorr_RatioSat->Write();

h_ClustSat_WithoutCorr->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
h_ClustSat_CorrHSCP_WithoutInversion->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
h_ClustSat_MyCorr_WithoutInversion->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
h_ClustSat_MyCorr_WithoutInversion_RatioSat->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
h_ClustSat_WithoutCorr->Write();
h_ClustSat_CorrHSCP_WithoutInversion->Write();
h_ClustSat_CorrHSCP->Write();
h_ClustSat_MyCorr_WithoutInversion->Write();
h_ClustSat_MyCorr->Write();
h_ClustSat_MyCorr_WithoutInversion_RatioSat->Write();
h_ClustSat_MyCorr_RatioSat->Write();

h_ClustNoSat_WithoutCorr->Write();
h_ClustNoSat_CorrHSCP_WithoutInversion->Write();
h_ClustNoSat_CorrHSCP->Write();
h_ClustNoSat_MyCorr_WithoutInversion->Write();
h_ClustNoSat_MyCorr->Write();
h_ClustNoSat_MyCorr_WithoutInversion_RatioSat->Write();
h_ClustNoSat_MyCorr_RatioSat->Write();

hMyCorrModulGoem->Write();

TCanvas* cStudyCorr = new TCanvas();
TLegend* legendStudyCorr = new TLegend(0.7,0.7,0.9,0.9);
hWithoutCorr->SetLineColor(1);
hWithoutCorr->SetLineWidth(2);
legendStudyCorr->AddEntry(hWithoutCorr,"W/O correction & W/O inversion","l");
hCorrHSCP_WithoutInversion->SetLineColor(40);
hCorrHSCP_WithoutInversion->SetLineWidth(2);
legendStudyCorr->AddEntry(hCorrHSCP_WithoutInversion,"Hscp correction W/O inversion","l");
hMyCorr_WithoutInversion_RatioSat->SetLineColor(46);
hMyCorr_WithoutInversion_RatioSat->SetLineWidth(2);
legendStudyCorr->AddEntry(hMyCorr_WithoutInversion_RatioSat,"Dylan's correction W/O inversion","l");
hWithoutCorr->GetYaxis()->SetRangeUser(0,15000);
hWithoutCorr->Draw();
hCorrHSCP_WithoutInversion->Draw("same");
hMyCorr_WithoutInversion_RatioSat->Draw("same");
legendStudyCorr->Draw("same");
cStudyCorr->Write();



TCanvas* cStudyCorrOnlySat = new TCanvas("cStudyCorrOnlySat","");

TLegend* legendStudyCorrOnlySat = new TLegend(0.7,0.7,0.9,0.9);
legendStudyCorrOnlySat->AddEntry(h_ClustSat_WithoutCorr,"No correction and no inversion","l");
legendStudyCorrOnlySat->AddEntry(h_ClustSatAndHscpTest_WithoutCorr,"No correction and no inversion - passing the correction criteria","l");
legendStudyCorrOnlySat->AddEntry(h_ClustSat_MyFullCorr,"Dylan's correction","l");
legendStudyCorrOnlySat->AddEntry(h_ClustSatAndHscpTest_MyFullCorr,"Dylan's correction - passing the correction criteria","l");
legendStudyCorrOnlySat->AddEntry(h_ClustSat_HscpCorrNoInv,"Hscp correction and no inversion","l");
legendStudyCorrOnlySat->AddEntry(h_ClustSatAndHscpTest_HscpCorrNoInv,"Hscp correction and no inversion - passing the correction criteria","l");

h_ClustSat_WithoutCorr->SetLineColor(1);
h_ClustSatAndHscpTest_WithoutCorr->SetLineStyle(2);
h_ClustSatAndHscpTest_WithoutCorr->SetLineColor(1);

h_ClustSat_MyFullCorr->SetLineColor(46);
h_ClustSatAndHscpTest_MyFullCorr->SetLineStyle(2);
h_ClustSatAndHscpTest_MyFullCorr->SetLineColor(46);

h_ClustSat_HscpCorrNoInv->SetLineColor(40);
h_ClustSatAndHscpTest_HscpCorrNoInv->SetLineStyle(2);
h_ClustSatAndHscpTest_HscpCorrNoInv->SetLineColor(40);


h_ClustSat_MyFullCorr->Draw();
h_ClustSatAndHscpTest_MyFullCorr->Draw("same");
h_ClustSat_HscpCorrNoInv->Draw("same");
h_ClustSatAndHscpTest_HscpCorrNoInv->Draw("same");
h_ClustSat_WithoutCorr->Draw("same");
h_ClustSatAndHscpTest_WithoutCorr->Draw("same");
legendStudyCorrOnlySat->Draw("same");


cStudyCorrOnlySat->Write();








TCanvas* cmass = new TCanvas();
TLegend* legendmass = new TLegend(0.7,0.7,0.9,0.9);
TLine* line2400 = new TLine(2400,0,2400,100);
legendmass->AddEntry(massEsim,"M(E_{sim})","l");
legendmass->AddEntry(massEfullLayer,"M(E_{full}^{Layer})","l");
massEsim->GetXaxis()->SetTitle("m [GeV/c^{2}]");
massEsim->GetYaxis()->SetTitle("[u.a.]");
massEsim->Draw();
massEfullLayer->Draw("SAME");
legendmass->Draw("same");
line2400->Draw("same");
cmass->Write();

TCanvas* cmassModul = new TCanvas();
TLegend* legendmassModul = new TLegend(0.7,0.7,0.9,0.9);
legendmassModul->AddEntry(massEsim,"M(E_{sim})","l");
legendmassModul->AddEntry(massEfullLayer,"M(E_{full}^{ModulGeom})","l");
massEsim->GetXaxis()->SetTitle("m [GeV/c^{2}]");
massEsim->GetYaxis()->SetTitle("[u.a.]");
massEsim->Draw();
massEfullModulGeom->Draw("SAME");
legendmassModul->Draw("same");
line2400->Draw("same");
cmassModul->Write();


//cout<<massEfullLayer->GetMean()<<endl;
//cout<<massEfullLayer->GetStdDev()<<endl;
//cout<<massEfullModulGeom->GetMean()<<endl;
//cout<<massEfullModulGeom->GetStdDev()<<endl;





massEsim->Write();
massErec->Write();
hmassEcorr->Write();
massEfullLayer->Write();
massEfullModulGeom->Write();


pVsdEsimdx->Write();
pVsdErecdx->Write();
pVsdEcorrdx->Write();
pVsdEfulldx_Layer->Write();
pVsdEfulldx_ModulGeom->Write();




distribmassprotonEloss->Write();
distribmassprotonEcorr->Write();
distribmassprotonEcorrNew->Write();


TCanvas* canvasKCecorrnew = new TCanvas();
h2PvsDecorrDxProtonKCnewSatLayer->Draw();
historececorrnew->Draw("same");
canvasKCecorrnew->Write();

TCanvas* canvasKCecorr = new TCanvas();
h2PvsDecorrDxProtonKC->Draw();
historececorr->Draw("same");
canvasKCecorr->Write();

TCanvas* canvasKC = new TCanvas();
h2PvsDelossDxProtonKC->Draw();
historeceloss->Draw("same");
canvasKC->Write();


h2PoverMDecorrDxLayer->Write();
h2PoverMDecorrDxModulGeom->Write();

h2PvsDecorrDxLayer->Write();
h2PvsDecorrDxModulGeom->Write();

hmodulgeom->Write();
hlayerlabel->Write();
SetHistoLabel(hmodulgeomvslayer);
hmodulgeomvslayer->Write();

TCanvas* cLayer = new TCanvas("layer","");
hEcorrLayer->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
hEcorrLayer->Draw();
hErecLayer->SetLineColor(2);
hErecLayer->Draw("same");
cLayer->Write();

TCanvas* cModulGeom = new TCanvas("ModulGeom","");
hEcorrModulGeom->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
hEcorrModulGeom->Draw();
hErecModulGeom->SetLineColor(2);
hErecModulGeom->Draw("same");
cModulGeom->Write();

TCanvas* cModulGeomNoSat = new TCanvas("ModulGeomNoSat","");
hEcorrModulGeomNoSat->Draw();
hErecModulGeomNoSat->SetLineColor(2);
hErecModulGeomNoSat->Draw("same");
cModulGeomNoSat->Write();

TCanvas* cLayerModulGeom = new TCanvas("LayerModulGeom","");
TLegend* legEcorrLayerModulGeom = new TLegend(0.7,0.7,0.9,0.9);
legEcorrLayerModulGeom->AddEntry(hEcorrLayer,"Layer","l");
legEcorrLayerModulGeom->AddEntry(hEcorrModulGeom,"ModulGeom","l");
hEcorrLayer->GetXaxis()->SetTitle("#frac{E_{corr}-E_{sim}}{E_{sim}}");
hEcorrLayer->Draw();
hEcorrModulGeom->SetLineColor(2);
hEcorrModulGeom->Draw("same");
legEcorrLayerModulGeom->Draw("same");
cLayerModulGeom->Write();

h2ElossErec->GetXaxis()->SetTitle("E_{sim} [GeV]");
h2ElossErec->GetYaxis()->SetTitle("E_{rec} [GeV]");

h2ElossErec->Write();
h2ElossErecLayer->Write();
h2PvsDecorrDxProtonKC->Write();
h2ElossErecModulGeom->Write();
h2ElossEcorrModulGeom->Write();
h2ElossErecModulGeomNoSat->Write();
h2ElossEcorrModulGeomNoSat->Write();


h2EcorrLayerVsEcorrModulGeom->Write();
h1EcorrLayerVsEcorrModulGeom->Write();


h2ElossEcorrSatLayer->Write();







outputfile->Write();
outputfile->Close();

delete outputfile;




return 0;

}