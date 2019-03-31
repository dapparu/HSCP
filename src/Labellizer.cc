#include "../interface/Labellizer.h"

using namespace std;


const float m_proton = 938.27*pow(10,-3); //masse du proton en MeV
const float m_pion = 139.57*pow(10,-3); //masse du pion en MeV
const float m_Rhadrons = 2400; //masse des R-hadrons (proche masse gluino) en MeV 

char StirngToChar(string str)
{
    char c[str.length()+1];
    strcpy(c,str);
    return c;
}

void SetHistoLabel(TCanvas* canvas,TH1F* histo)
{

    int Wsize = histo->GetMaximum()*1.05;
    histo->GetYaxis()->SetRangeUser(0,histo->GetMaximum()*1.05);

    TLine* line1 = new TLine(5,0,5,0.25);
    TLine* line2 = new TLine(11,0,11,0.25);
    TLine* line3 = new TLine(13,0,13,0.25);

    if(Wsize!=0)
    {
        line1->SetY2(Wsize);
        line2->SetY2(Wsize);
        line3->SetY2(Wsize);
    }


    line1->SetLineStyle(2);
    line2->SetLineStyle(2);
    line3->SetLineStyle(2);
    
    histo->GetXaxis()->SetBinLabel(1+1,"TIB 1");
    histo->GetXaxis()->SetBinLabel(1+2,"TIB 2");
    histo->GetXaxis()->SetBinLabel(1+3,"TIB 3");
    histo->GetXaxis()->SetBinLabel(1+4,"TIB 4");
    histo->GetXaxis()->SetBinLabel(1+5,"TOB 1");
    histo->GetXaxis()->SetBinLabel(1+6,"TOB 2");
    histo->GetXaxis()->SetBinLabel(1+7,"TOB 3");
    histo->GetXaxis()->SetBinLabel(1+8,"TOB 4");
    histo->GetXaxis()->SetBinLabel(1+9,"TOB 5");
    histo->GetXaxis()->SetBinLabel(1+10,"TOB 6");
    histo->GetXaxis()->SetBinLabel(1+11,"TID 1");
    histo->GetXaxis()->SetBinLabel(1+12,"TID 2");
    histo->GetXaxis()->SetBinLabel(1+13,"TEC 1");
    histo->GetXaxis()->SetBinLabel(1+14,"TEC 2");
    histo->GetXaxis()->SetBinLabel(1+15,"TEC 3");
    histo->GetXaxis()->SetBinLabel(1+16,"TEC 4");
    histo->GetXaxis()->SetBinLabel(1+17,"TEC 5");
    histo->GetXaxis()->SetBinLabel(1+18,"TEC 6");
    histo->GetXaxis()->SetBinLabel(1+19,"TEC 7");
    histo->GetXaxis()->SetBinLabel(1+20,"TEC 8");
    histo->GetXaxis()->SetBinLabel(1+21,"TEC 9");
    histo->GetXaxis()->LabelsOption("v");
    canvas->cd();
    histo->Draw();
    line1->Draw();
    line2->Draw();
    line3->Draw();
    
}

void SetHistoLabelPartID(TCanvas* canvas,TH2F* histo)
{
    histo->GetXaxis()->SetBinLabel(1+1,"#pi^{+/-}");
    histo->GetXaxis()->SetBinLabel(1+2,"p/#bar{p}");
    histo->GetXaxis()->SetBinLabel(1+3,"R^{+/-}_{#tilde{g}u#bar{d}}");
    histo->GetXaxis()->SetBinLabel(1+4,"R^{+/-}");
    canvas->cd();
    histo->Draw();
}

string Label(int i)
{
    if(i==1)        return "TIB1";
    else if(i==2)   return "TIB2";
    else if(i==3)   return "TIB3";
    else if(i==4)   return "TIB4";
    else if(i==5)   return "TOB1";
    else if(i==6)   return "TOB2";
    else if(i==7)   return "TOB3";
    else if(i==8)   return "TOB4";
    else if(i==9)   return "TOB5";
    else if(i==10)  return "TOB6";
    else if(i==11)  return "TID1";
    else if(i==12)  return "TID2";
    else if(i==13)  return "TEC1";
    else if(i==14)  return "TEC2";
    else if(i==15)  return "TEC3";
    else if(i==16)  return "TEC4";
    else if(i==17)  return "TEC5";
    else if(i==18)  return "TEC6";
    else if(i==19)  return "TEC7";
    else if(i==20)  return "TEC8";
    else if(i==21)  return "TEC9";
}

string LabelParticle(int i)
{
    if(i==211 || i==-211) return "#pi";
    if(i==2212 || i==-2212) return "p";
    if(i==1009213 || i==-1009213) return "R^{+}_{#tilde{g}u#bar{d}}";
    else if((int)i/1000==1009 || (int)i/1000==-1009) return "R^{+}";
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
    if(i==211 || i==-211) return 1;
    else if(i==2212 || i==-2212) return 2;
    else if(i==1009213 || i==-1009213) return 3;
    else if((int)i/1000==1009 || (int)i/1000==-1009) return 4;
}

float GetPoverM(float p,int i)
{
    if(i==211 || i==-211) return p/m_pion;
    if(i==2212 || i==2212) return p/m_proton;
    if((int)i/1000==1009 || (int)i/1000==-1009) return p/m_Rhadrons;
}
