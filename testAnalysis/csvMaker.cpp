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






//CONSTANTES K ET C

const float CsusyCorr = 1.14;	//pour gluino
const float KsusyCorr = 6.7;



const float CprotonEloss = 3.27;
const float KprotonEloss = 2.49;

const float CprotonEcorrNew = 6.73;
const float KprotonEcorrNew = 1.18;

const float CprotonEcorr = 6.15;
const float KprotonEcorr = 2.91;

const float MGluino = 2400.;







//FUNCTIONS -- LABEL HISTO

void StudyGraph (ofstream & ofile, TH1F * histo)
{
  float DistribMean = histo->GetMean ();
  float DistribStdDev = histo->GetStdDev ();
  float DistribMeanErr = histo->GetMeanError ();
  float DistribStdDevErr = histo->GetStdDevError ();
  int FitStatus = histo->Fit ("gaus", "QRN", "", -0.1, 0.1);
  //TFitResultPtr fit = histo->Fit("gaus","QRS","",DistribMean-DistribStdDev,DistribMean+DistribStdDev);
  TFitResultPtr fit = histo->Fit ("gaus", "QRSN", "", -0.1, 0.1);
  float FitMean = 0;
  float FitSigma = 0;
  float FitMeanErr = 0;
  float FitSigmaErr = 0;
  if (FitStatus > -1) {
    FitMean = fit->Parameter (1);
    FitSigma = fit->Parameter (2);
    FitMeanErr = fit->Error (1);
    FitSigmaErr = fit->Error (2);
  }
  int CountPass = 0;
  int CountPassGlissantSigma = 0;
  int CountPassGlissant = 0;
  int CountPassRight = 0;
  int CountPassLeft = 0;
  for (int i = 0; i < histo->GetNbinsX (); i++) {
    if (histo->GetBinCenter (i) >= FitMean - FitSigma && histo->GetBinCenter (i) <= FitMean + FitSigma) {
      CountPassGlissantSigma += histo->GetBinContent (i);
    }
    if (histo->GetBinCenter (i) >= FitMean - 0.34 && histo->GetBinCenter (i) <= FitMean + 0.34) {
      CountPassGlissant += histo->GetBinContent (i);
    }
    if (histo->GetBinCenter (i) >= -0.25 && histo->GetBinCenter (i) <= 0.25) {
      CountPass += histo->GetBinContent (i);
    }
    if (histo->GetBinCenter (i) >= 0.4)
      CountPassRight += histo->GetBinContent (i);
    if (histo->GetBinCenter (i) <= -0.4)
      CountPassLeft += histo->GetBinContent (i);
  }
  float FracInterval = (float) CountPass / (float) histo->GetEntries ();
  float FracIntervalGlissantSigma = (float) CountPassGlissantSigma / (float) histo->GetEntries ();
  float FracIntervalGlissant = (float) CountPassGlissant / (float) histo->GetEntries ();
  float FracRight = (float) CountPassRight / (float) histo->GetEntries ();
  float FracLeft = (float) CountPassLeft / (float) histo->GetEntries ();
  //ofile<<"Name"<<"\t\t"<<"DistribMean"<<"\t"<<"DistribMeanErr"<<"\t"<<"DistribStdDev"<<"\t"<<"DistribStdDevErr"<<"\t"<<"FitMean"<<"\t"<<"FitMeanErr"<<"\t"<<"FitSigma"<<"\t"<<"FitSigmaErr"<<"\t"<<"FracIntervalGlissantSigma"<<"\t"<<"FracInterval"<<endl;
  ofile << endl;
  ofile << histo->GetName () << "\t\t" << histo->
    GetMaximum () << "\t" << DistribMean << "\t" << DistribMeanErr << "\t" << DistribStdDev << "\t" << DistribStdDevErr << "\t FitMean : " << FitMean << "\t" << FitMeanErr << "\t FitSigma : " <<
    FitSigma << "\t" << FitSigmaErr << "\t" << FracIntervalGlissantSigma << "\t" << FracIntervalGlissant << "\t" << FracInterval << "\t FracRight : " << FracRight << "\t FracLeft : " << FracLeft <<
    endl;

}

void SetHistoLabel (TH2I * histo)
{
  histo->GetXaxis ()->SetBinLabel (1, "TIB 1");
  histo->GetXaxis ()->SetBinLabel (2, "TIB 2");
  histo->GetXaxis ()->SetBinLabel (3, "TIB 3");
  histo->GetXaxis ()->SetBinLabel (4, "TIB 4");
  histo->GetXaxis ()->SetBinLabel (5, "TOB 1");
  histo->GetXaxis ()->SetBinLabel (6, "TOB 2");
  histo->GetXaxis ()->SetBinLabel (7, "TOB 3");
  histo->GetXaxis ()->SetBinLabel (8, "TOB 4");
  histo->GetXaxis ()->SetBinLabel (9, "TOB 5");
  histo->GetXaxis ()->SetBinLabel (10, "TOB 6");
  histo->GetXaxis ()->SetBinLabel (11, "TID 1");
  histo->GetXaxis ()->SetBinLabel (12, "TID 2");
  histo->GetXaxis ()->SetBinLabel (13, "TEC 1");
  histo->GetXaxis ()->SetBinLabel (14, "TEC 2");
  histo->GetXaxis ()->SetBinLabel (15, "TEC 3");
  histo->GetXaxis ()->SetBinLabel (16, "TEC 4");
  histo->GetXaxis ()->SetBinLabel (17, "TEC 5");
  histo->GetXaxis ()->SetBinLabel (18, "TEC 6");
  histo->GetXaxis ()->SetBinLabel (19, "TEC 7");
  histo->GetXaxis ()->SetBinLabel (20, "TEC 8");
  histo->GetXaxis ()->SetBinLabel (21, "TEC 9");
  histo->GetXaxis ()->LabelsOption ("v");

  histo->GetYaxis ()->SetBinLabel (1, "IB1");
  histo->GetYaxis ()->SetBinLabel (2, "IB2");
  histo->GetYaxis ()->SetBinLabel (3, "OB1");
  histo->GetYaxis ()->SetBinLabel (4, "OB2");
  histo->GetYaxis ()->SetBinLabel (5, "W1A");
  histo->GetYaxis ()->SetBinLabel (6, "W2A");
  histo->GetYaxis ()->SetBinLabel (7, "W3A");
  histo->GetYaxis ()->SetBinLabel (8, "W1B");
  histo->GetYaxis ()->SetBinLabel (9, "W2B");
  histo->GetYaxis ()->SetBinLabel (10, "W3B");
  histo->GetYaxis ()->SetBinLabel (11, "W4");
  histo->GetYaxis ()->SetBinLabel (12, "W5");
  histo->GetYaxis ()->SetBinLabel (13, "W6");
  histo->GetYaxis ()->SetBinLabel (14, "W7");

}

void SetHistoLabel (TCanvas * canvas, TH1D * histo)
{

  int Wsize = histo->GetMaximum () * 1.05;
  histo->GetYaxis ()->SetRangeUser (0, histo->GetMaximum () * 1.05);

  TLine *line1 = new TLine (4, 0, 4, 0.25);
  TLine *line2 = new TLine (10, 0, 10, 0.25);
  TLine *line3 = new TLine (12, 0, 12, 0.25);

  if (Wsize != 0) {
    line1->SetY2 (Wsize);
    line2->SetY2 (Wsize);
    line3->SetY2 (Wsize);
  }


  line1->SetLineStyle (2);
  line2->SetLineStyle (2);
  line3->SetLineStyle (2);

  histo->GetXaxis ()->SetBinLabel (1, "TIB 1");
  histo->GetXaxis ()->SetBinLabel (2, "TIB 2");
  histo->GetXaxis ()->SetBinLabel (3, "TIB 3");
  histo->GetXaxis ()->SetBinLabel (4, "TIB 4");
  histo->GetXaxis ()->SetBinLabel (5, "TOB 1");
  histo->GetXaxis ()->SetBinLabel (6, "TOB 2");
  histo->GetXaxis ()->SetBinLabel (7, "TOB 3");
  histo->GetXaxis ()->SetBinLabel (8, "TOB 4");
  histo->GetXaxis ()->SetBinLabel (9, "TOB 5");
  histo->GetXaxis ()->SetBinLabel (10, "TOB 6");
  histo->GetXaxis ()->SetBinLabel (11, "TID 1");
  histo->GetXaxis ()->SetBinLabel (12, "TID 2");
  histo->GetXaxis ()->SetBinLabel (13, "TEC 1");
  histo->GetXaxis ()->SetBinLabel (14, "TEC 2");
  histo->GetXaxis ()->SetBinLabel (15, "TEC 3");
  histo->GetXaxis ()->SetBinLabel (16, "TEC 4");
  histo->GetXaxis ()->SetBinLabel (17, "TEC 5");
  histo->GetXaxis ()->SetBinLabel (18, "TEC 6");
  histo->GetXaxis ()->SetBinLabel (19, "TEC 7");
  histo->GetXaxis ()->SetBinLabel (20, "TEC 8");
  histo->GetXaxis ()->SetBinLabel (21, "TEC 9");
  histo->GetXaxis ()->LabelsOption ("v");
  canvas->cd ();
  histo->Draw ();
  line1->Draw ();
  line2->Draw ();
  line3->Draw ();

}


string LabelModulGeom (int modulgeom)
{
  if (modulgeom == 1)
    return "IB1";
  if (modulgeom == 2)
    return "IB2";
  if (modulgeom == 3)
    return "OB1";
  if (modulgeom == 4)
    return "OB2";
  if (modulgeom == 5)
    return "W1A";
  if (modulgeom == 6)
    return "W2A";
  if (modulgeom == 7)
    return "W3A";
  if (modulgeom == 8)
    return "W1B";
  if (modulgeom == 9)
    return "W2B";
  if (modulgeom == 10)
    return "W3B";
  if (modulgeom == 11)
    return "W4";
  if (modulgeom == 12)
    return "W5";
  if (modulgeom == 13)
    return "W6";
  if (modulgeom == 14)
    return "W7";
}

string LabelLayer (int i)
{
  if (i == 1)
    return "TIB1";
  if (i == 2)
    return "TIB2";
  if (i == 3)
    return "TIB3";
  if (i == 4)
    return "TIB4";
  if (i == 5)
    return "TOB1";
  if (i == 6)
    return "TOB2";
  if (i == 7)
    return "TOB3";
  if (i == 8)
    return "TOB4";
  if (i == 9)
    return "TOB5";
  if (i == 10)
    return "TOB6";
  if (i == 11)
    return "TID1";
  if (i == 12)
    return "TID2";
  if (i == 13)
    return "TEC1";
  if (i == 14)
    return "TEC2";
  if (i == 15)
    return "TEC3";
  if (i == 16)
    return "TEC4";
  if (i == 17)
    return "TEC5";
  if (i == 18)
    return "TEC6";
  if (i == 19)
    return "TEC7";
  if (i == 20)
    return "TEC8";
  if (i == 21)
    return "TEC9";
}

vector < float >FactorKC (TH2F * histo, float mass)	//histo pre-filtre sur le particleID --> il faut donner la masse correspondante 
{

  TProfile profile;
  TH2F *histoclone = (TH2F *) histo->Clone ();
  histoclone->Reset ();
  TFitResultPtr fit_res;
  TH1D *projY;
  TH1F *historec = new TH1F ("", "", 1000, 0, 5);
  int divider = 50;
  for (int i = 1; i <= histo->GetNbinsX (); i++) {
    for (int j = 2; j <= histo->GetNbinsY (); j++) {
      if (histo->GetBinContent (i, j) > 0)
	histoclone->SetBinContent (i, j, histo->GetBinContent (i, j));
    }
    if (i % divider == 0) {
      projY = histoclone->ProjectionY ();
      projY->Rebin (5);
      if (projY->GetEntries () > 0) {
	fit_res = projY->Fit ("gaus", "0QS");
	TH1D *projX = histoclone->ProjectionX ();
	//projX->Rebin(5);
	//TFitResultPtr fit_resX = projX->Fit("gaus","QS");
	//historec.SetBinContent(projX->GetMean(),fit_res->Parameter(1));
	//historec.SetBinError(projX->GetMean(),fit_res->Error(1));
	historec->SetBinContent (i - (float) divider / 2, fit_res->Parameter (1));
	historec->SetBinError (i - (float) divider / 2, fit_res->Error (1));
	historec->SetMarkerStyle (3);
	historec->SetMarkerColor (4);

      }
      histoclone->Reset ();
    }
  }
  TF1 *InvSquare = new TF1 ("InvSquare", "[0]+[1]*(1/pow(x,2))");
  TFitResultPtr fit = historec->Fit ("InvSquare", "QRS", "", 0, 3000);
  vector < float >vectres;
  vectres.push_back (fit->Parameter (1) / (mass * mass));	//K
  vectres.push_back (fit->Parameter (0));	//C
  vectres.push_back (fit->Error (1) / (mass * mass));
  vectres.push_back (fit->Error (0));
  return vectres;
}


vector < int >CrossTalkInv (Correction corr, const std::vector < int >&Q, const float x1, const float x2, float threshold, float thresholdSat, int label, int correctType, float RatioSat254, float RatioSat255, float ThresholdRatioSat)	//Methode de correction + XTalkInv
{

  bool DylanCorr = true;
  bool inversion = true;
  if (correctType == 1) {
    DylanCorr = false;
    inversion = false;
  }
  if (correctType == 2) {
    DylanCorr = false;
  }
  if (correctType == 3) {
    inversion = false;
  }

  float charge = 0.;

  const unsigned N = Q.size ();
  std::vector < int >QII;
  std::vector < float >QI (N, 0);
  Double_t a = 1 - 2 * x1 - 2 * x2;
  //  bool debugbool=false;
  TMatrix A (N, N);

  float RatioSatTotal = RatioSat254 + RatioSat255;

  vector < int >::const_iterator mQ = max_element (Q.begin (), Q.end ());

  //----------------------

  if (!DylanCorr) {
    //---  que pour 1 max bien net 
    if (Q.size () < 2 || Q.size () > 8) {
      for (unsigned int i = 0; i < Q.size (); i++) {
	QII.push_back (Q[i]);
      }
      return QII;
    }

    if (*mQ > 253) {

      if (*mQ == 255 && *(mQ - 1) > 253 && *(mQ + 1) > 253)
	return Q;

      if (*(mQ - 1) > thresholdSat && *(mQ + 1) > thresholdSat && *(mQ - 1) < 254 && *(mQ + 1) < 254 && abs (*(mQ - 1) - *(mQ + 1)) < 40)
	//if(*(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40)
	//if(RatioSatTotal>=ThresholdRatioSat)
      {
	QII.push_back ((10 * (*(mQ - 1)) + 10 * (*(mQ + 1))) / 2);
	return QII;
      }
    }
  }

  //----------------------

  if (DylanCorr && RatioSatTotal >= ThresholdRatioSat && *mQ > 253 && (*(mQ - 1) > thresholdSat && *(mQ + 1) > thresholdSat && *(mQ - 1) < 254 && *(mQ + 1) < 254 && abs (*(mQ - 1) - *(mQ + 1)) < 40)) {
    int NSat254 = 0;
    int NSat255 = 0;
    float ClusterCharge = 0;

    for (int i = 0; i < Q.size (); i++) {
      if (Q[i] == 254)
	NSat254++;
      if (Q[i] == 255)
	NSat255++;
      ClusterCharge += (Q[i]);
    }

    //float ClusterChargeCorr=corr.ChargeCorr(ClusterCharge*(3.61*pow(10,-9)*247),label,Q.size(),NSat254,NSat255)/(3.61*pow(10,-9)*247);
    float ClusterChargeCorr = corr.ChargeCorr (ClusterCharge, label, Q.size (), NSat254, NSat255) * 1.038;
    float DiffClusterCharge = ClusterChargeCorr - ClusterCharge;
    for (unsigned int i = 0; i < Q.size (); i++) {
      QII.push_back (Q[i]);
    }
    //Incomplet pour le moment --> a voir au moment de ClusterCleaning
    //Le surplus de charge ne se fait pas que sur une piste 

    vector < int >::iterator maxQ = max_element (QII.begin (), QII.end ());
    if (DiffClusterCharge >= 0)
      QII.at (std::distance (QII.begin (), maxQ)) += DiffClusterCharge;
    return QII;
  }

  //----------------------

  if (inversion) {
    for (unsigned int i = 0; i < N; i++) {
      A (i, i) = a;
      if (i < N - 1) {
	A (i + 1, i) = x1;
	A (i, i + 1) = x1;
      }
      else
	continue;
      if (i < N - 2) {
	A (i + 2, i) = x2;
	A (i, i + 2) = x2;
      }
    }

    if (N == 1)
      A (0, 0) = 1 / a;
    else
      A.InvertFast ();

    for (unsigned int i = 0; i < N; i++) {
      for (unsigned int j = 0; j < N; j++) {
	QI[i] += A (i, j) * (float) Q[j];
      }
    }

    for (unsigned int i = 0; i < QI.size (); i++) {
      if (QI[i] < threshold)
	QI[i] = 0;
      QII.push_back ((int) QI[i]);
    }
    return QII;

  }

  if (!inversion)
    return Q;

}

vector < int >CrossTalkInvStudy (Correction corr, const std::vector < int >&Q, const float x1, const float x2, float threshold, float thresholdSat, float thresholdDiff, int label, int correctType, float RatioSat254, float RatioSat255, float ThresholdRatioSat, bool & testThresholdSat, bool & testThresholdDiff, bool & testBothThreshold)	//Methode de correction + XTalkInv
{
  testThresholdSat = false;
  testThresholdDiff = false;
  testBothThreshold = false;

  bool DylanCorr = true;
  bool inversion = true;
  if (correctType == 1) {
    DylanCorr = false;
    inversion = false;
  }
  if (correctType == 2) {
    DylanCorr = false;
  }
  if (correctType == 3) {
    inversion = false;
  }

  float charge = 0.;

  const unsigned N = Q.size ();
  std::vector < int >QII;
  std::vector < float >QI (N, 0);
  Double_t a = 1 - 2 * x1 - 2 * x2;
  //  bool debugbool=false;
  TMatrix A (N, N);

  float RatioSatTotal = RatioSat254 + RatioSat255;

  vector < int >::const_iterator mQ = max_element (Q.begin (), Q.end ());

  //----------------------

  if (!DylanCorr) {
    //---  que pour 1 max bien net 
    if (Q.size () < 2 || Q.size () > 8) {
      for (unsigned int i = 0; i < Q.size (); i++) {
	QII.push_back (Q[i]);
      }
      return QII;
    }

    if (*mQ > 253) {

      if (*mQ == 255 && *(mQ - 1) > 253 && *(mQ + 1) > 253)
	return Q;

      if (*(mQ - 1) > thresholdSat && *(mQ + 1) > thresholdSat && *(mQ - 1) < 254 && *(mQ + 1) < 254) {
	testThresholdSat = true;
      }

      if (abs (*(mQ - 1) - *(mQ + 1)) < thresholdDiff && *(mQ - 1) < 254 && *(mQ + 1) < 254) {
	testThresholdDiff = true;
      }

      if (*(mQ - 1) > thresholdSat && *(mQ + 1) > thresholdSat && *(mQ - 1) < 254 && *(mQ + 1) < 254 && abs (*(mQ - 1) - *(mQ + 1)) < thresholdDiff) {
	testBothThreshold = true;
      }
      if (testThresholdSat) {
	QII.push_back ((10 * (*(mQ - 1)) + 10 * (*(mQ + 1))) / 2);
	return QII;
      }
      if (testThresholdDiff) {
	QII.push_back ((10 * (*(mQ - 1)) + 10 * (*(mQ + 1))) / 2);
	return QII;
      }
      if (testBothThreshold) {
	QII.push_back ((10 * (*(mQ - 1)) + 10 * (*(mQ + 1))) / 2);
	return QII;
      }

    }
  }

  //----------------------

  if (DylanCorr && RatioSatTotal >= ThresholdRatioSat && *mQ > 253
      && (*(mQ - 1) > thresholdSat && *(mQ + 1) > thresholdSat && *(mQ - 1) < 254 && *(mQ + 1) < 254 && abs (*(mQ - 1) - *(mQ + 1)) < thresholdDiff)) {
    int NSat254 = 0;
    int NSat255 = 0;
    float ClusterCharge = 0;

    for (int i = 0; i < Q.size (); i++) {
      if (Q[i] == 254)
	NSat254++;
      if (Q[i] == 255)
	NSat255++;
      ClusterCharge += (Q[i]);
    }

    //float ClusterChargeCorr=corr.ChargeCorr(ClusterCharge*(3.61*pow(10,-9)*247),label,Q.size(),NSat254,NSat255)/(3.61*pow(10,-9)*247);
    float ClusterChargeCorr = corr.ChargeCorr (ClusterCharge, label, Q.size (), NSat254, NSat255);
    float DiffClusterCharge = ClusterChargeCorr - ClusterCharge;
    for (unsigned int i = 0; i < Q.size (); i++) {
      QII.push_back (Q[i]);
    }
    //Incomplet pour le moment --> a voir au moment de ClusterCleaning
    //Le surplus de charge ne se fait pas que sur une piste 

    vector < int >::iterator maxQ = max_element (QII.begin (), QII.end ());
    if (DiffClusterCharge >= 0)
      QII.at (std::distance (QII.begin (), maxQ)) += DiffClusterCharge;
    return QII;
  }


  if (!inversion)
    return Q;

  //----------------------

  if (inversion) {
    for (unsigned int i = 0; i < N; i++) {
      A (i, i) = a;
      if (i < N - 1) {
	A (i + 1, i) = x1;
	A (i, i + 1) = x1;
      }
      else
	continue;
      if (i < N - 2) {
	A (i + 2, i) = x2;
	A (i, i + 2) = x2;
      }
    }

    if (N == 1)
      A (0, 0) = 1 / a;
    else
      A.InvertFast ();

    for (unsigned int i = 0; i < N; i++) {
      for (unsigned int j = 0; j < N; j++) {
	QI[i] += A (i, j) * (float) Q[j];
      }
    }

    for (unsigned int i = 0; i < QI.size (); i++) {
      if (QI[i] < threshold)
	QI[i] = 0;
      QII.push_back ((int) QI[i]);
    }
    return QII;

  }


}







//MAIN



int main (int argc, char **argv)
{







  TChain chain;
  chain.SetName ("stage/ttree");
  for (int i = 1; i < argc - 1; i++)
    chain.Add (argv[i]);
  string s1 = argv[1];
  string s2 = s1.substr (0, s1.find ('.')) + "_results.root";

  Builder *b1 = new Builder (chain);
  b1->SetBranchAdd();


  int nentries = b1->GetEntries ();
  cout<<"nentries: "<<nentries<<endl;
  b1->GetEntry(1);
  cout<<b1->GetNtracks()<<endl;
  //return 0;

//OUVERTURE FICHIERS DE CONFIGS POUR LA METHODE

  Correction EcorrLayer;
  TFile *fileLayer = TFile::Open ("testMethodLayer.root");
  //TFile* fileLayer = TFile::Open("testMethodModulGeom.root");
  TTree *treeLayer = (TTree *) fileLayer->Get ("tree");
  EcorrLayer.SetTree (*treeLayer);
  EcorrLayer.ReadLayer ();
  //EcorrLayer.ReadModulGeom();

  Correction EcorrModulGeom;
  TFile *fileModulGeom = TFile::Open ("testMethodModulGeom.root");
  //TFile* fileModulGeom = TFile::Open("testMethodModulGeom.root");
  TTree *treeModulGeom = (TTree *) fileModulGeom->Get ("tree");
  EcorrModulGeom.SetTree (*treeModulGeom);
  EcorrModulGeom.ReadModulGeom ();



// OUVERTURE FICHIER CONFIG FACTEURS K ET C

  ifstream inputfileFactorKC ("factKC.txt");
  float K_Esim, C_Esim;
  float K_Erec, C_Erec;
  float K_Ecorr, C_Ecorr;
  float K_EfullLayer, C_EfullLayer;
  float K_EfullModulGeom, C_EfullModulGeom;

  inputfileFactorKC >> K_Esim;
  inputfileFactorKC >> C_Esim;
  inputfileFactorKC >> K_Erec;
  inputfileFactorKC >> C_Erec;
  inputfileFactorKC >> K_Ecorr;
  inputfileFactorKC >> C_Ecorr;
  inputfileFactorKC >> K_EfullLayer;
  inputfileFactorKC >> C_EfullLayer;
  inputfileFactorKC >> K_EfullModulGeom;
  inputfileFactorKC >> C_EfullModulGeom;




  bool method = false;
  bool corr = false;
  bool FactKC = false;
  bool XTalkInvStudy = false;

  bool gluino = false;

  float idlow = 2212, idup = 2212;

  float ratsat = 0.;

//CHOIX OPTION



  if (atof (argv[argc - 1]) == 1)
    method = true;
  if (atof (argv[argc - 1]) == 2)
    corr = true;
  if (atof (argv[argc - 1]) == 3) {
    corr = true;
    FactKC = true;
  }
  if (atof (argv[argc - 1]) == 4)
    XTalkInvStudy = true;


  if (atof (argv[argc - 2]) == 1) {
    gluino = true;
    idlow = 1000001;
    idup = 2000015;
  }

  ratsat = atof (argv[argc - 3]);


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

  vector < vector < vector < vector < TH2F * >>>>VectLayerVectNStripVectNStripSat254VectNStripSat255Histo;
  vector < vector < vector < vector < TH2F * >>>>VectModulGeomVectNStripVectNStripSat254VectNStripSat255Histo;
  vector < vector < vector < vector < TH2F * >>>>VectModulGeomVectNStripVectNStripSat254VectNStripSat255HistoNoSat;
  if (method) {
//METHODE 1 -- CATEGORIES LAYER
    for (int layer = 1; layer < 22; layer++) {
      vector < vector < vector < TH2F * >>>VectNStripVectNStripSat254VectNStripSat255Histo;
      for (int nstrip = 3; nstrip < 7; nstrip++) {
	vector < vector < TH2F * >>VectNStripSat254VectNStripSat255Histo;
	for (int nstripsat254 = 0; nstripsat254 < nstrip + 1; nstripsat254++) {
	  vector < TH2F * >VectNStripSat255Histo;
	  for (int nstripsat255 = 0; nstripsat255 < nstrip - nstripsat254 + 1; nstripsat255++) {
	    string title = LabelLayer (layer) + " NStrip=" + to_string (nstrip) + " NStripSat254=" + to_string (nstripsat254) + " NStripSat255=" + to_string (nstripsat255);
	    //VectNStripSat255Histo.push_back(new TH2F(title.c_str(),title.c_str(),300,0,6000*pow(10,-6),300,0,6000*pow(10,-6)));
	    VectNStripSat255Histo.push_back (new TH2F (title.c_str (), title.c_str (), 200, 0, 5000, 200, 0, 5000));
	  }
	  VectNStripSat254VectNStripSat255Histo.push_back (VectNStripSat255Histo);
	}
	VectNStripVectNStripSat254VectNStripSat255Histo.push_back (VectNStripSat254VectNStripSat255Histo);
      }
      VectLayerVectNStripVectNStripSat254VectNStripSat255Histo.push_back (VectNStripVectNStripSat254VectNStripSat255Histo);
    }

//METHODE 2 -- CATEGORIES XTALK

    for (int ModulGeom = 1; ModulGeom < 15; ModulGeom++) {
      vector < vector < vector < TH2F * >>>VectNStripVectNStripSat254VectNStripSat255Histo;
      for (int nstrip = 3; nstrip < 7; nstrip++) {
	vector < vector < TH2F * >>VectNStripSat254VectNStripSat255Histo;
	for (int nstripsat254 = 0; nstripsat254 < nstrip + 1; nstripsat254++) {
	  vector < TH2F * >VectNStripSat255Histo;
	  for (int nstripsat255 = 0; nstripsat255 < nstrip - nstripsat254 + 1; nstripsat255++) {
	    string title = LabelModulGeom (ModulGeom) + " NStrip=" + to_string (nstrip) + " NStripSat254=" + to_string (nstripsat254) + " NStripSat255=" + to_string (nstripsat255);
	    //VectNStripSat255Histo.push_back(new TH2F(title.c_str(),title.c_str(),300,0,6000*pow(10,-6),300,0,6000*pow(10,-6)));
	    VectNStripSat255Histo.push_back (new TH2F (title.c_str (), title.c_str (), 200, 0, 5000, 200, 0, 5000));
	  }
	  VectNStripSat254VectNStripSat255Histo.push_back (VectNStripSat255Histo);
	}
	VectNStripVectNStripSat254VectNStripSat255Histo.push_back (VectNStripSat254VectNStripSat255Histo);
      }
      VectModulGeomVectNStripVectNStripSat254VectNStripSat255Histo.push_back (VectNStripVectNStripSat254VectNStripSat255Histo);
    }



  }


//DECLARATION HISTOGRAMMES POUR ETUDES

ofstream ofile("data.csv");


//Ratio nombre de pistes satures 

  int ncluster = 0;
  int sat1 = 0;
  int sat2 = 0;
  int sat3 = 0;
  int sat1HscpTest = 0;
  int nSat = 0;
  int nHscpTest = 0;
  int nRatSatTest = 0;
  int nHscpAndRatSatTest = 0;
  int sat2RatSatTest = 0;
  int surcorr = 0;
  int surcorrHSCP = 0;

  int nParticle = 0;
  int nPatriclePassingMyCrit = 0;

  int counterpass = 0;

  int countThresholdSat = 0;
  int countThresholdDiff = 0;
  int countBothThreshold = 0;
  int countSatCluster = 0;


  b1->SetThresholdPartId(0.7);
  b1->SetThresholdPt(0);
  b1->SetThresholdP(0);
  b1->SetThresholdEta(5);

//BOUCLE SUR LES EVENTS

  int entries = nentries;

  for (int i = 0; i < entries; i++) {
    if (i % 1000 == 0)
      cout << "Event " << i << endl;
    b1->GetEntry (i);

    //BOUCLE SUR LES TRACES

    cout<<b1->GetNtracks()<<endl;
    for (int track = 0; track < b1->GetNtracks (); track++) {
      vector < int >vect_partID;
      vector < float >vect_Eloss;
      vector < float >vect_charge;
      vector < float >vect_dedx;
      vector < float >vect_dqdx;
      vector < float >vect_dqcorrdx;
      vector < float >vect_dqcorrdx2;
      vector < float >vect_plength;

      vector < float >vect_decorrdxnewSatLayer;
      vector < float >vect_decorrdx;

      vector < float >dEsimdx;
      vector < float >dErecdx;
      vector < float >dEcorrdx;
      vector < float >dEfulldx_Layer;
      vector < float >dEfulldx_ModulGeom;



      vector < float >VectDecorrDxLayer;
      vector < float >VectDecorrDxModulGeom;
      float pt = b1->GetVectTrack ()[track].GetPt ();
      float p = b1->GetVectTrack ()[track].GetP ();
      float eta = b1->GetVectTrack ()[track].GetEta ();
      float phi = b1->GetVectTrack ()[track].GetPhi ();
      float ias_ampl = b1->GetVectTrack ()[track].GetIasAmpl ();
      int NCluster = b1->GetVectTrack ()[track].GetNCluster ();
      int NClustSat254 = b1->GetVectTrack ()[track].GetNSatCluster (254);
      int NClustSat255 = b1->GetVectTrack ()[track].GetNSatCluster (255);
      float RatioNClusterSat254 = (double) NClustSat254 / (double) NCluster;
      float RatioNClusterSat255 = (double) NClustSat255 / (double) NCluster;
      int id = b1->GetVectTrack ()[track].GetPartId ();
      float PoverM = GetPoverM (p, id);



      nParticle++;
      if (RatioNClusterSat254 + RatioNClusterSat255 >= ratsat)
	nPatriclePassingMyCrit++;


      bool testSurCorr = false;
      bool testMoreThanOne = false;
      bool testHscptrack = false;


      std::vector < float >dEdxsim;
      std::vector < float >dEdxMyCorr;
      std::vector < float >dEdxHscpCorrWOinv;
      std::vector < float >dEdxHscpCorr;


      cout<<idlow<<" "<<id<<" "<<idup<<endl;
      if (idlow <= abs (id) && abs (id) <= idup)
      {
       cout<<"ca passe"<<endl;
	//BOUCLE SUR LES CLUSTERS

	for (int cluster = 0; cluster < b1->GetVectTrack ()[track].GetNCluster (); cluster++) {
	  int fakecluster = 0;

	  float charge = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetSclusCharge ();
	  float charge_corr = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetSclusChargeCorr ();
	  float EcorrBef = charge_corr * (3.61 * pow (10, -9) * 247);
	  float Erec = charge * (3.61 * pow (10, -9) * 247);
	  float chargeNoSat = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetChargeWithoutSaturation ();
	  float ErecNoSat = chargeNoSat * (3.61 * pow (10, -9) * 247);
	  float Eloss = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetEloss ();
	  float pathlength = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetPathLength ();
	  int layer = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetLayer ();
	  int layerLabel = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetLayerLabel ();
	  bool sat254 = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetSat254 ();
	  bool sat255 = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetSat255 ();
	  bool shape = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetShape ();
	  int nsimhits = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetNSimHits ();
	  int nstrips = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetNStrip ();
	  int nsat254 = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetNSatStrip (254);
	  int nsat255 = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetNSatStrip (255);
	  bool edge = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].Edge ();
	  bool cut = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].Cut ();
	  int detid = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetDetId ();
	  int subdetid = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetSubDetId ();
	  int modulgeom = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetModulGeom ();
	  int maxstrip = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetMaxStrip ();
	  int cluster_id = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetPartId ();
	  float EcorrSatLayer = Erec;
	  float EcorrSatModul = Erec;
	  float EcorrNoSat = ErecNoSat;
	  int charge_sim = Eloss / (3.61 * pow (10, -9) * 247);

	  bool fakeclust = true;

	  for (int simhit = 0; simhit < b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetNSimHits (); simhit++) {
	    int simhit_id = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectSimHits ()[simhit].GetPartId ();
	    if (simhit_id == id)
	      fakeclust = false;
	  }

	  float Ratio_fakeclust = (float) fakecluster / (float) NCluster;

	  bool testThresholdSat = false;
	  bool testThresholdDiff = false;
	  bool testBothThreshold = false;

/*
	    h2ElossErec->Fill(Eloss,Erec);

			hmodulgeom->Fill(modulgeom);
			hlayerlabel->Fill(layerLabel);
			hmodulgeomvslayer->Fill(layerLabel-1,modulgeom-1);
*/
	  vect_dedx.push_back (Eloss / pathlength);

	  //BOUCLE SUR LES STRIPS D'UN CLUSTER

	  int ChargeCluster = 0;
	  vector < int >VectChargeStrip;

	  int CentralChargeAndNeighbours = 0;



	  float Qminus = 0.;
	  float Qplus = 0.;
	  bool testQminus = true;

	  float StripNoSat = 0;
	  float StripSat = 0;


	  for (int strip = 0; strip < b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetNStrip (); strip++) {
	    float ChargeStrip = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectStrips ()[strip].GetAmpl ();
	    VectChargeStrip.push_back (b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectStrips ()[strip].GetAmpl ());
	    ChargeCluster += b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectStrips ()[strip].GetAmpl ();

	    if (ChargeStrip > 253 && testQminus) {
	      testQminus = false;
	      Qminus = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectStrips ()[strip - 1].GetAmpl ();
	    }
	    if (ChargeStrip > 253) {
	      Qplus = b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectStrips ()[strip + 1].GetAmpl ();
	    }
	    if (ChargeStrip < 254)
	      StripNoSat += ChargeStrip;
	    if (ChargeStrip > 253)
	      StripSat += ChargeStrip;

	    if (ChargeStrip > 253) {
	      CentralChargeAndNeighbours += ChargeStrip;
	      CentralChargeAndNeighbours += b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectStrips ()[strip - 1].GetAmpl ();
	      if (b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectStrips ()[strip + 1].GetAmpl () > 253) {
		CentralChargeAndNeighbours += b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectStrips ()[strip + 1].GetAmpl ();
		if (b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectStrips ()[strip + 2].GetAmpl () > 253) {
		  CentralChargeAndNeighbours += b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectStrips ()[strip + 2].GetAmpl ();
		  CentralChargeAndNeighbours += b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectStrips ()[strip + 3].GetAmpl ();
		}
		else
		  CentralChargeAndNeighbours += b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectStrips ()[strip + 2].GetAmpl ();
	      }
	      else
		CentralChargeAndNeighbours += b1->GetVectTrack ()[track].GetVectClusters ()[cluster].GetVectStrips ()[strip + 1].GetAmpl ();
	    }

	  }

	  bool testCorrHscp = false;
	  vector < int >::const_iterator mQ = max_element (VectChargeStrip.begin (), VectChargeStrip.end ());
	  if (*mQ > 253 && *(mQ - 1) > 25 && *(mQ + 1) > 25 && *(mQ - 1) < 254 && *(mQ + 1) < 254 && abs (*(mQ - 1) - *(mQ + 1)) < 40)
	    testCorrHscp = true;

	  float DeltaQ = abs (*(mQ - 1) - *(mQ + 1));
	  float NewTestEnergy = CentralChargeAndNeighbours * (3.61 * pow (10, -9) * 247);
	  float DeltaQCalc = abs (Qplus - Qminus);


	  float ChargeMinEstimated = StripNoSat + 1.2 * StripSat;

	  //EcorrSatModul = EcorrModulGeom.ChargeCorr(Erec,modulgeom,nstrips,nsat254,nsat255);

	  //--------------------------------



	  if (corr) {
	    if (layerLabel <= 22)
	      dEsimdx.push_back (Eloss / pathlength);
	    dErecdx.push_back (Erec / pathlength);
	    dEcorrdx.push_back (EcorrBef / pathlength);
	    if (nstrips >= 3) {
	      if ((RatioNClusterSat254 + RatioNClusterSat255 >= ratsat)) {
		EcorrSatLayer = EcorrLayer.ChargeCorr (Erec, layerLabel, nstrips, nsat254, nsat255);
		dEfulldx_Layer.push_back (EcorrSatLayer / pathlength);
		EcorrSatModul = EcorrModulGeom.ChargeCorr (Erec, modulgeom, nstrips, nsat254, nsat255);
		dEfulldx_ModulGeom.push_back (EcorrSatModul / pathlength);



	      }
	      else {
		if (!(sat254 || sat255)) {
		  dEfulldx_Layer.push_back (Erec / pathlength);
		  dEfulldx_ModulGeom.push_back (Erec / pathlength);
		}
	      }
	    }
	  }









	}			//fin boucle cluster

	Estimator estimLayer (VectDecorrDxLayer);
	Estimator estimModulGeom (VectDecorrDxModulGeom);
	Estimator estimDelossDx (vect_dedx);


        cout<<p<<","<<eta<<","<<NCluster<<","<<NClustSat254<<endl;
        ofile<<p<<","<<eta<<","<<NCluster<<","<<NClustSat254<<endl;


	float K = KsusyCorr;
	float C = CsusyCorr;
	//float dedx = estimLayer.GetHarmonic2()*pow(10,3);
	float dedx = estimDelossDx.GetHarmonic2 () * pow (10, 3);

	float mass_dedx = sqrt ((dedx - C) * pow (p, 2) / K);





	Estimator estimDecorrDx (vect_decorrdx);
	Estimator estimDecorrDxnewSatLayer (vect_decorrdxnewSatLayer);
	float decorrdx = estimDecorrDx.GetHarmonic2 () * pow (10, 3);
	float decorrdxnewsatlayer = estimDecorrDxnewSatLayer.GetHarmonic2 () * pow (10, 3);

	//if(decorrdxnewsatlayer>0) h2PvsDecorrDxProtonKCnewSatLayer->Fill(p,decorrdxnewsatlayer);





	float massEloss = sqrt ((dedx - CprotonEloss) * pow (p, 2) / KprotonEloss);
	float massEcorr = sqrt ((decorrdx - CprotonEcorr) * pow (p, 2) / KprotonEcorr);
	float massEcorrNew = sqrt ((decorrdxnewsatlayer - CprotonEcorrNew) * pow (p, 2) / KprotonEcorrNew);
	//float massEloss = sqrt((dedx-CprotonEloss)*pow(p,2)/KprotonEloss);
	//float massEcorr = sqrt((decorrdx-CprotonEcorr)*pow(p,2)/KprotonEcorr);
	//float massEcorrNew = sqrt((decorrdxnewsatlayer-CprotonEcorrNew)*pow(p,2)/KprotonEcorrNew);


	float CsusySim = 4.53;	//pour gluino
	float KsusySim = 1.76;




	Estimator estimDeDx (dEsimdx);
	//pVsdEsimdx->Fill(p,estimDeDx.GetHarmonic2()*pow(10,3));
	//massEsim->Fill(sqrt((estimDeDx.GetHarmonic2()*pow(10,3)-C_Esim)*pow(p,2)/K_Esim));

	estimDeDx.SetVect (dErecdx);

	estimDeDx.SetVect (dEcorrdx);

	estimDeDx.SetVect (dEfulldx_Layer);

	estimDeDx.SetVect (dEfulldx_ModulGeom);




	Estimator estimdEdxsim (dEdxsim);
	Estimator estimdEdxMyCorr (dEdxMyCorr);
	Estimator estimdEdxHscpCorrWOinv (dEdxHscpCorrWOinv);
	Estimator estimdEdxHscpCorr (dEdxHscpCorr);








      }

    }
  }				//END LOOP ENTRIES

//Ratio number of saturated strip

ofile.close();

//ECRITURE RESULTATS ETUDES

}
