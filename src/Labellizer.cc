#include "../interface/Labellizer.h"

using namespace std;


const float m_electron = 511*pow(10,-6); //masse de l'electron en GeV
const float m_pion = 139.57*pow(10,-3); //masse du pion en GeV
const float m_kaon = 439.68*pow(10,-3); //masse du kaon en GeV
const float m_proton = 938.27*pow(10,-3); //masse du proton en GeV
const float m_Rhadrons = 2400; //masse des R-hadrons (proche masse gluino) en GeV 

const char* StringToChar(string str)
{
    return str.c_str();
}

void SetHistoLabel(TCanvas* canvas,TH1F* histo)
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

void SetHistoLabelPartID(TCanvas* canvas,TH2F* histo)
{
    histo->GetXaxis()->SetBinLabel(1,"e^{+/-}");
    histo->GetXaxis()->SetBinLabel(2,"#pi^{+/-}");
    histo->GetXaxis()->SetBinLabel(3,"K^{+/-}");
    histo->GetXaxis()->SetBinLabel(4,"p/#bar{p}");
    histo->GetXaxis()->SetBinLabel(5,"SM particle");
    histo->GetXaxis()->SetBinLabel(6,"R^{+/-}_{#tilde{g}u#bar{d}}");
    histo->GetXaxis()->SetBinLabel(7,"SUSY particle");
    canvas->cd();
    histo->Draw();
}

string LabelParticle(int i)
{
    if(i==11 || i==-11) return "e^{+/-}";
    if(i==321 || i==-321) return "K^{+/-}";
    if(i==211 || i==-211) return "#pi^{+/-}";
    if(i==2212 || i==-2212) return "p/#bar{p}";
    if(i==1009213 || i==-1009213) return "R^{+}_{#tilde{g}u#bar{d}}";
    //if(((int)i/1000==1009 || (int)i/1000==-1009) && !(i==1009213 || i==-1009213)) return "R^{+}";
    if((1000001<=i && i<=2000015) || (-2000015<=i && i<=-1000001)) return "SUSY particle";
    else return "SM particle";
}

int GetPartID(const vector<Cluster> &VectClust,float &threshold)
{
    int PartID=0;
    vector<int> VectPartID_Cluster;
    for(int cluster=0;cluster<VectClust.size();cluster++)
    {
        for(int simhit=0;simhit<VectClust[cluster].GetVectSimHits().size();simhit++)
        {
            PartID=VectClust[cluster].GetVectSimHits()[simhit].GetPartId();
            VectPartID_Cluster.push_back(PartID);
        }
    }
    sort(VectPartID_Cluster.begin(),VectPartID_Cluster.end());
    vector<int> VectPartID;
    vector<int> counter;
    VectPartID.push_back(VectPartID_Cluster[0]);
    counter.push_back(1);
    int j=0;
    for(int i=1;i<VectPartID_Cluster.size();i++)
    {
        if(VectPartID_Cluster[i]==VectPartID[j])
        {
            counter[j]++;
        }
        else
        {
            VectPartID.push_back(VectPartID_Cluster[i]);
            counter.push_back(1);
            j++;
        }
    }
    vector<float> ratio;
    for(int i=0;i<counter.size();i++)
    {
        ratio.push_back((float)counter[i]/(float)VectPartID_Cluster.size());
    }
    float max=0;
    int indice=0;
    for(int i=0;i<ratio.size();i++)
    {
        if(ratio[i]>max)
        {
            max=ratio[i];
            indice=i;
        }
    }
    threshold=max;
    return VectPartID[indice];
}

int ReBinPartID(int i)
{
    if(i==11 || i==-11) return 0; //electron
    if(i==211 || i==-211) return 1; //pion
    if(i==321 || i==-321) return 2; //kaon
    if(i==2212 || i==-2212) return 3; //proton
    if(i==1009213 || i==-1009213) return 5; //R-hadron 
    //if(((int)i/1000==1009 || (int)i/1000==-1009) && !(i==1009213 || i==-1009213)) return 6;
    if((1000001<=i && i<=2000015) || (-2000015<=i && i<=-1000001)) return 6; //SUSY
    else return 4; //SM
}

float GetPoverM(float p,int i)
{
    if(i==11 || i==-11) return p/m_electron;
    if(i==321 || i==-321) return p/m_kaon;
    if(i==211 || i==-211) return p/m_pion;
    if(i==2212 || i==-2212) return p/m_proton;
    if(i==1009213 || i==-1009213) return p/m_Rhadrons;
    if((1000001<=i && i<=2000015) || (-2000015<=i && i<=-1000001)) return p/m_Rhadrons;
}

string LoopPartID(int i)
{
    if(i==0) return "e^{+/-}";
    if(i==1) return "#pi^{+/-}";
    if(i==2) return "K^{+/-}";
    if(i==3) return "p/#bar{p}";
    if(i==4) return "SM";
    if(i==5) return "R^{+/-}_{#tilde{g}u#bar{d}}";
    if(i==6) return "SUSY";
}

string Label(int i)
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

int GetThickness(int layerLabel)
{
    if(layerLabel==1 || layerLabel==2) return 80; //en um
    if(layerLabel==3 || layerLabel==4) return 120;
    if(layerLabel==5 || layerLabel==6 || layerLabel==7 || layerLabel==8) return 183;
    if(layerLabel==9 || layerLabel==10) return 122;
    else return 100;        //attention, a modifier lorsque je m'interesserai aux valeurs pour le TEC 
}

int GetPitch(int layerLabel)
{
    if(layerLabel==1 || layerLabel==2 || layerLabel==3 || layerLabel==4) return 320; //en um
    if(layerLabel==5 || layerLabel==6 || layerLabel==7 || layerLabel==8 || layerLabel==9 || layerLabel==10) return 500;
    else return 100;        //attention, a modifier lorsque je m'interesserai aux valeurs pour le TEC 
}