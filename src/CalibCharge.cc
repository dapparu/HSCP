#include "../interface/CalibCharge.h"

using namespace std;


CalibCharge::CalibCharge()
{
    float p0_TOB1NStrip3NStripSat1=0,       p1_TOB1NStrip3NStripSat1=1;
    float p0_TOB1NStrip4NStripSat1=0,       p1_TOB1NStrip4NStripSat1=1;
    float p0_TOB1NStrip5NStripSat1=0,       p1_TOB1NStrip5NStripSat1=1;
    float p0_TOB1NStrip6orMoreNStripSat1=0, p1_TOB1NStrip6orMoreNStripSat1=1;
}

CalibCharge::~CalibCharge()
{

}

bool CalibCharge::SelectedArea(float x1,float y1,float x2,float y2,float x3,float y3,float x4,float y4,float x,float y)
{
    bool res=false;
    float coef1 = (y3-y1)/(x3-x1);
    float coef2 = (y4-y2)/(x4-x2);
    float coef3 = (y2-y1)/(x2-x1);
    float coef4 = (y4-y3)/(x4-x3);
    float orign1 = y1-coef1*x1;
    float orign2 = y4-coef2*x4;
    //float orign3 = y2-coef3*x2;
    //float orign4 = y3-coef4*x3;
    //if(x>=(y-orign3)/coef3 && x<=(y-orign4)/coef4 && y>=coef1*x+orign1 && y<=coef2*x+orign2) res=true;
    if(x>=x1 && x>=x2 && x<=x3 && x<=x4 && y>=coef1*x+orign1 && y<=coef2*x+orign2) res=true;
    return res;    
}

bool CalibCharge::Area(int layerLabel,int NStrip,float eloss,float charge)
{
    if(layerLabel==5 && NStrip==3) return SelectedArea(221*pow(10,-6),252*pow(10,-6),221*pow(10,-6),294*pow(10,-6),618*pow(10,-6),339*pow(10,-6),618*pow(10,-6),386*pow(10,-6),eloss,charge);
    if(layerLabel==5 && NStrip==4) return SelectedArea(232*pow(10,-6),268*pow(10,-6),232*pow(10,-6),303*pow(10,-6),893*pow(10,-6),411*pow(10,-6),893*pow(10,-6),458*pow(10,-6),eloss,charge);
    if(layerLabel==5 && NStrip==5) return SelectedArea(402*pow(10,-6),314*pow(10,-6),402*pow(10,-6),354*pow(10,-6),1297*pow(10,-6),542*pow(10,-6),1297*pow(10,-6),597*pow(10,-6),eloss,charge);
}

void CalibCharge::SetParameters()
{
    string s_fit = "./data/fit_res.txt";

    ifstream ifile_fit(s_fit.c_str(),ios::in);

    if(ifile_fit)
    {
        while(!ifile_fit.eof())
        {
            ifile_fit>>TOB1NStrip3NStripSat1>>p0_TOB1NStrip3NStripSat1>>p1_TOB1NStrip3NStripSat1;
            ifile_fit>>TOB1NStrip4NStripSat1>>p0_TOB1NStrip4NStripSat1>>p1_TOB1NStrip4NStripSat1;
            ifile_fit>>TOB1NStrip5NStripSat1>>p0_TOB1NStrip5NStripSat1>>p1_TOB1NStrip5NStripSat1;
            ifile_fit>>TOB1NStrip6orMoreNStripSat1>>p0_TOB1NStrip6orMoreNStripSat1>>p1_TOB1NStrip6orMoreNStripSat1;
        }
        ifile_fit.close();
    }
    else cerr<<"Fichier de fit non ouvert"<<endl;
}

float CalibCharge::Charge(int layerLabel,int NStrip,float NStripSat,float charge)
{
    if(layerLabel==5 && NStrip==3 && NStripSat==1) return -(p0_TOB1NStrip3NStripSat1-charge)/p1_TOB1NStrip3NStripSat1;
    if(layerLabel==5 && NStrip==4 && NStripSat==1) return -(p0_TOB1NStrip4NStripSat1-charge)/p1_TOB1NStrip4NStripSat1;
    if(layerLabel==5 && NStrip==5 && NStripSat==1) return -(p0_TOB1NStrip5NStripSat1-charge)/p1_TOB1NStrip5NStripSat1;
    if(layerLabel==5 && NStrip>=6 && NStripSat==1) return -(p0_TOB1NStrip6orMoreNStripSat1-charge)/p1_TOB1NStrip6orMoreNStripSat1;
    else return charge;
}