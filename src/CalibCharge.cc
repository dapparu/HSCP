#include "../interface/CalibCharge.h"

using namespace std;


CalibCharge::CalibCharge()
{
    p0_TIB1_254_NStrip3_NStripSat1=0, p1_TIB1_254_NStrip3_NStripSat1=1;
    p0_TIB1_254_NStrip4_NStripSat1=0, p1_TIB1_254_NStrip4_NStripSat1=1;
    p0_TIB1_254_NStrip4_NStripSat2=0, p1_TIB1_254_NStrip4_NStripSat2=1;
    p0_TIB1_254_NStrip5_NStripSat1=0, p1_TIB1_254_NStrip5_NStripSat1=1;
    p0_TIB1_254_NStrip5_NStripSat2=0, p1_TIB1_254_NStrip5_NStripSat2=1;
    p0_TIB1_254_NStrip6_NStripSat1=0, p1_TIB1_254_NStrip6_NStripSat1=1;

    p0_TIB1_255_NStrip3_NStripSat1=0, p1_TIB1_255_NStrip3_NStripSat1=1;
    p0_TIB1_255_NStrip4_NStripSat1=0, p1_TIB1_255_NStrip4_NStripSat1=1;
    p0_TIB1_255_NStrip4_NStripSat2=0, p1_TIB1_255_NStrip4_NStripSat2=1;
    p0_TIB1_255_NStrip5_NStripSat1=0, p1_TIB1_255_NStrip5_NStripSat1=1;
    p0_TIB1_255_NStrip5_NStripSat2=0, p1_TIB1_255_NStrip5_NStripSat2=1;
    p0_TIB1_255_NStrip6_NStripSat1=0, p1_TIB1_255_NStrip6_NStripSat1=1;

    p0_TIB2_254_NStrip3_NStripSat1=0, p1_TIB2_254_NStrip3_NStripSat1=1;
    p0_TIB2_254_NStrip4_NStripSat1=0, p1_TIB2_254_NStrip4_NStripSat1=1;
    p0_TIB2_254_NStrip4_NStripSat2=0, p1_TIB2_254_NStrip4_NStripSat2=1;
    p0_TIB2_254_NStrip5_NStripSat1=0, p1_TIB2_254_NStrip5_NStripSat1=1;
    p0_TIB2_254_NStrip5_NStripSat2=0, p1_TIB2_254_NStrip5_NStripSat2=1;
    p0_TIB2_254_NStrip6_NStripSat1=0, p1_TIB2_254_NStrip6_NStripSat1=1;

    p0_TIB2_255_NStrip3_NStripSat1=0, p1_TIB2_255_NStrip3_NStripSat1=1;
    p0_TIB2_255_NStrip4_NStripSat1=0, p1_TIB2_255_NStrip4_NStripSat1=1;
    p0_TIB2_255_NStrip4_NStripSat2=0, p1_TIB2_255_NStrip4_NStripSat2=1;
    p0_TIB2_255_NStrip5_NStripSat1=0, p1_TIB2_255_NStrip5_NStripSat1=1;
    p0_TIB2_255_NStrip5_NStripSat2=0, p1_TIB2_255_NStrip5_NStripSat2=1;
    p0_TIB2_255_NStrip6_NStripSat1=0, p1_TIB2_255_NStrip6_NStripSat1=1;

    p0_TIB3_254_NStrip3_NStripSat1=0, p1_TIB3_254_NStrip3_NStripSat1=1;
    p0_TIB3_254_NStrip4_NStripSat1=0, p1_TIB3_254_NStrip4_NStripSat1=1;
    p0_TIB3_254_NStrip4_NStripSat2=0, p1_TIB3_254_NStrip4_NStripSat2=1;
    p0_TIB3_254_NStrip5_NStripSat1=0, p1_TIB3_254_NStrip5_NStripSat1=1;
    p0_TIB3_254_NStrip5_NStripSat2=0, p1_TIB3_254_NStrip5_NStripSat2=1;
    p0_TIB3_254_NStrip6_NStripSat1=0, p1_TIB3_254_NStrip6_NStripSat1=1;

    p0_TIB3_255_NStrip3_NStripSat1=0, p1_TIB3_255_NStrip3_NStripSat1=1;
    p0_TIB3_255_NStrip4_NStripSat1=0, p1_TIB3_255_NStrip4_NStripSat1=1;
    p0_TIB3_255_NStrip4_NStripSat2=0, p1_TIB3_255_NStrip4_NStripSat2=1;
    p0_TIB3_255_NStrip5_NStripSat1=0, p1_TIB3_255_NStrip5_NStripSat1=1;
    p0_TIB3_255_NStrip5_NStripSat2=0, p1_TIB3_255_NStrip5_NStripSat2=1;
    p0_TIB3_255_NStrip6_NStripSat1=0, p1_TIB3_255_NStrip6_NStripSat1=1;

    p0_TIB4_254_NStrip3_NStripSat1=0, p1_TIB4_254_NStrip3_NStripSat1=1;
    p0_TIB4_254_NStrip4_NStripSat1=0, p1_TIB4_254_NStrip4_NStripSat1=1;
    p0_TIB4_254_NStrip4_NStripSat2=0, p1_TIB4_254_NStrip4_NStripSat2=1;
    p0_TIB4_254_NStrip5_NStripSat1=0, p1_TIB4_254_NStrip5_NStripSat1=1;
    p0_TIB4_254_NStrip5_NStripSat2=0, p1_TIB4_254_NStrip5_NStripSat2=1;
    p0_TIB4_254_NStrip6_NStripSat1=0, p1_TIB4_254_NStrip6_NStripSat1=1;

    p0_TIB4_255_NStrip3_NStripSat1=0, p1_TIB4_255_NStrip3_NStripSat1=1;
    p0_TIB4_255_NStrip4_NStripSat1=0, p1_TIB4_255_NStrip4_NStripSat1=1;
    p0_TIB4_255_NStrip4_NStripSat2=0, p1_TIB4_255_NStrip4_NStripSat2=1;
    p0_TIB4_255_NStrip5_NStripSat1=0, p1_TIB4_255_NStrip5_NStripSat1=1;
    p0_TIB4_255_NStrip5_NStripSat2=0, p1_TIB4_255_NStrip5_NStripSat2=1;
    p0_TIB4_255_NStrip6_NStripSat1=0, p1_TIB4_255_NStrip6_NStripSat1=1;

    p0_TOB1_254_NStrip3_NStripSat1=0, p1_TOB1_254_NStrip3_NStripSat1=1;
    p0_TOB1_254_NStrip4_NStripSat1=0, p1_TOB1_254_NStrip4_NStripSat1=1;
    p0_TOB1_254_NStrip4_NStripSat2=0, p1_TOB1_254_NStrip4_NStripSat2=1;
    p0_TOB1_254_NStrip5_NStripSat1=0, p1_TOB1_254_NStrip5_NStripSat1=1;
    p0_TOB1_254_NStrip5_NStripSat2=0, p1_TOB1_254_NStrip5_NStripSat2=1;
    p0_TOB1_254_NStrip6_NStripSat1=0, p1_TOB1_254_NStrip6_NStripSat1=1;

    p0_TOB1_255_NStrip3_NStripSat1=0, p1_TOB1_255_NStrip3_NStripSat1=1;
    p0_TOB1_255_NStrip4_NStripSat1=0, p1_TOB1_255_NStrip4_NStripSat1=1;
    p0_TOB1_255_NStrip4_NStripSat2=0, p1_TOB1_255_NStrip4_NStripSat2=1;
    p0_TOB1_255_NStrip5_NStripSat1=0, p1_TOB1_255_NStrip5_NStripSat1=1;
    p0_TOB1_255_NStrip5_NStripSat2=0, p1_TOB1_255_NStrip5_NStripSat2=1;
    p0_TOB1_255_NStrip6_NStripSat1=0, p1_TOB1_255_NStrip6_NStripSat1=1;

    p0_TOB2_254_NStrip3_NStripSat1=0, p1_TOB2_254_NStrip3_NStripSat1=1;
    p0_TOB2_254_NStrip4_NStripSat1=0, p1_TOB2_254_NStrip4_NStripSat1=1;
    p0_TOB2_254_NStrip4_NStripSat2=0, p1_TOB2_254_NStrip4_NStripSat2=1;
    p0_TOB2_254_NStrip5_NStripSat1=0, p1_TOB2_254_NStrip5_NStripSat1=1;
    p0_TOB2_254_NStrip5_NStripSat2=0, p1_TOB2_254_NStrip5_NStripSat2=1;
    p0_TOB2_254_NStrip6_NStripSat1=0, p1_TOB2_254_NStrip6_NStripSat1=1;

    p0_TOB2_255_NStrip3_NStripSat1=0, p1_TOB2_255_NStrip3_NStripSat1=1;
    p0_TOB2_255_NStrip4_NStripSat1=0, p1_TOB2_255_NStrip4_NStripSat1=1;
    p0_TOB2_255_NStrip4_NStripSat2=0, p1_TOB2_255_NStrip4_NStripSat2=1;
    p0_TOB2_255_NStrip5_NStripSat1=0, p1_TOB2_255_NStrip5_NStripSat1=1;
    p0_TOB2_255_NStrip5_NStripSat2=0, p1_TOB2_255_NStrip5_NStripSat2=1;
    p0_TOB2_255_NStrip6_NStripSat1=0, p1_TOB2_255_NStrip6_NStripSat1=1;

    p0_TOB3_254_NStrip3_NStripSat1=0, p1_TOB3_254_NStrip3_NStripSat1=1;
    p0_TOB3_254_NStrip4_NStripSat1=0, p1_TOB3_254_NStrip4_NStripSat1=1;
    p0_TOB3_254_NStrip4_NStripSat2=0, p1_TOB3_254_NStrip4_NStripSat2=1;
    p0_TOB3_254_NStrip5_NStripSat1=0, p1_TOB3_254_NStrip5_NStripSat1=1;
    p0_TOB3_254_NStrip5_NStripSat2=0, p1_TOB3_254_NStrip5_NStripSat2=1;
    p0_TOB3_254_NStrip6_NStripSat1=0, p1_TOB3_254_NStrip6_NStripSat1=1;

    p0_TOB3_255_NStrip3_NStripSat1=0, p1_TOB3_255_NStrip3_NStripSat1=1;
    p0_TOB3_255_NStrip4_NStripSat1=0, p1_TOB3_255_NStrip4_NStripSat1=1;
    p0_TOB3_255_NStrip4_NStripSat2=0, p1_TOB3_255_NStrip4_NStripSat2=1;
    p0_TOB3_255_NStrip5_NStripSat1=0, p1_TOB3_255_NStrip5_NStripSat1=1;
    p0_TOB3_255_NStrip5_NStripSat2=0, p1_TOB3_255_NStrip5_NStripSat2=1;
    p0_TOB3_255_NStrip6_NStripSat1=0, p1_TOB3_255_NStrip6_NStripSat1=1;

    p0_TOB4_254_NStrip3_NStripSat1=0, p1_TOB4_254_NStrip3_NStripSat1=1;
    p0_TOB4_254_NStrip4_NStripSat1=0, p1_TOB4_254_NStrip4_NStripSat1=1;
    p0_TOB4_254_NStrip4_NStripSat2=0, p1_TOB4_254_NStrip4_NStripSat2=1;
    p0_TOB4_254_NStrip5_NStripSat1=0, p1_TOB4_254_NStrip5_NStripSat1=1;
    p0_TOB4_254_NStrip5_NStripSat2=0, p1_TOB4_254_NStrip5_NStripSat2=1;
    p0_TOB4_254_NStrip6_NStripSat1=0, p1_TOB4_254_NStrip6_NStripSat1=1;

    p0_TOB4_255_NStrip3_NStripSat1=0, p1_TOB4_255_NStrip3_NStripSat1=1;
    p0_TOB4_255_NStrip4_NStripSat1=0, p1_TOB4_255_NStrip4_NStripSat1=1;
    p0_TOB4_255_NStrip4_NStripSat2=0, p1_TOB4_255_NStrip4_NStripSat2=1;
    p0_TOB4_255_NStrip5_NStripSat1=0, p1_TOB4_255_NStrip5_NStripSat1=1;
    p0_TOB4_255_NStrip5_NStripSat2=0, p1_TOB4_255_NStrip5_NStripSat2=1;
    p0_TOB4_255_NStrip6_NStripSat1=0, p1_TOB4_255_NStrip6_NStripSat1=1;

    p0_TOB5_254_NStrip3_NStripSat1=0, p1_TOB5_254_NStrip3_NStripSat1=1;
    p0_TOB5_254_NStrip4_NStripSat1=0, p1_TOB5_254_NStrip4_NStripSat1=1;
    p0_TOB5_254_NStrip4_NStripSat2=0, p1_TOB5_254_NStrip4_NStripSat2=1;
    p0_TOB5_254_NStrip5_NStripSat1=0, p1_TOB5_254_NStrip5_NStripSat1=1;
    p0_TOB5_254_NStrip5_NStripSat2=0, p1_TOB5_254_NStrip5_NStripSat2=1;
    p0_TOB5_254_NStrip6_NStripSat1=0, p1_TOB5_254_NStrip6_NStripSat1=1;

    p0_TOB5_255_NStrip3_NStripSat1=0, p1_TOB5_255_NStrip3_NStripSat1=1;
    p0_TOB5_255_NStrip4_NStripSat1=0, p1_TOB5_255_NStrip4_NStripSat1=1;
    p0_TOB5_255_NStrip4_NStripSat2=0, p1_TOB5_255_NStrip4_NStripSat2=1;
    p0_TOB5_255_NStrip5_NStripSat1=0, p1_TOB5_255_NStrip5_NStripSat1=1;
    p0_TOB5_255_NStrip5_NStripSat2=0, p1_TOB5_255_NStrip5_NStripSat2=1;
    p0_TOB5_255_NStrip6_NStripSat1=0, p1_TOB5_255_NStrip6_NStripSat1=1;

    p0_TOB6_254_NStrip3_NStripSat1=0, p1_TOB6_254_NStrip3_NStripSat1=1;
    p0_TOB6_254_NStrip4_NStripSat1=0, p1_TOB6_254_NStrip4_NStripSat1=1;
    p0_TOB6_254_NStrip4_NStripSat2=0, p1_TOB6_254_NStrip4_NStripSat2=1;
    p0_TOB6_254_NStrip5_NStripSat1=0, p1_TOB6_254_NStrip5_NStripSat1=1;
    p0_TOB6_254_NStrip5_NStripSat2=0, p1_TOB6_254_NStrip5_NStripSat2=1;
    p0_TOB6_254_NStrip6_NStripSat1=0, p1_TOB6_254_NStrip6_NStripSat1=1;

    p0_TOB6_255_NStrip3_NStripSat1=0, p1_TOB6_255_NStrip3_NStripSat1=1;
    p0_TOB6_255_NStrip4_NStripSat1=0, p1_TOB6_255_NStrip4_NStripSat1=1;
    p0_TOB6_255_NStrip4_NStripSat2=0, p1_TOB6_255_NStrip4_NStripSat2=1;
    p0_TOB6_255_NStrip5_NStripSat1=0, p1_TOB6_255_NStrip5_NStripSat1=1;
    p0_TOB6_255_NStrip5_NStripSat2=0, p1_TOB6_255_NStrip5_NStripSat2=1;
    p0_TOB6_255_NStrip6_NStripSat1=0, p1_TOB6_255_NStrip6_NStripSat1=1;
}

CalibCharge::~CalibCharge()
{

}

bool CalibCharge::SelectedArea(float x1,float y1,float x2,float y2,float x3,float y3,float x4,float y4,float x,float y)
{
    x1*=pow(10,-6),x2*=pow(10,-6),x3*=pow(10,-6),x4*=pow(10,-6);
    y1*=pow(10,-6),y2*=pow(10,-6),y3*=pow(10,-6),y4*=pow(10,-6);
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

bool CalibCharge::Area(int layerLabel,int NStrip,int NStripSat,bool sat254,float eloss,float charge)
{
    if(layerLabel==5 && NStrip==3 && NStripSat==1 && sat254==true) return SelectedArea(221,252,221,294,618,339,618,386,eloss,charge);
    if(layerLabel==5 && NStrip==4 && NStripSat==1 && sat254==true) return SelectedArea(232,268,232,303,893,411,893,458,eloss,charge);
    if(layerLabel==5 && NStrip==4 && NStripSat==2 && sat254==true) return SelectedArea(465,497,465,539,1082,572,1082,607,eloss,charge);
    if(layerLabel==5 && NStrip==5 && NStripSat==1 && sat254==true) return SelectedArea(402,314,402,354,1297,542,1297,597,eloss,charge);
    if(layerLabel==5 && NStrip==5 && NStripSat==2 && sat254==true) return SelectedArea(556,509,556,575,1444,631,1444,685,eloss,charge);
    if(layerLabel==5 && NStrip==6 && NStripSat==1 && sat254==true) return SelectedArea(342,325,342,372,1309,535,1309,590,eloss,charge);

    if(layerLabel==5 && NStrip==3 && NStripSat==1 && sat254==false) return SelectedArea(231,252,231,306,605,330,605,384,eloss,charge);
    if(layerLabel==5 && NStrip==4 && NStripSat==1 && sat254==false) return SelectedArea(247,265,247,338,803,390,803,456,eloss,charge);
    if(layerLabel==5 && NStrip==4 && NStripSat==2 && sat254==false) return true;
    if(layerLabel==5 && NStrip==5 && NStripSat==1 && sat254==false) return SelectedArea(263,276,263,354,1330,522,1330,617,eloss,charge);
    if(layerLabel==5 && NStrip==5 && NStripSat==2 && sat254==false) return true;
    if(layerLabel==5 && NStrip==6 && NStripSat==1 && sat254==false) return true;

    else return false;

    /*if(layerLabel==6 && NStrip==3 && NStripSat==1 && sat254==true) return SelectedArea(,eloss,charge);
    if(layerLabel==6 && NStrip==4 && NStripSat==1 && sat254==true) return SelectedArea(,eloss,charge);
    if(layerLabel==6 && NStrip==4 && NStripSat==2 && sat254==true) return SelectedArea(,eloss,charge);
    if(layerLabel==6 && NStrip==5 && NStripSat==1 && sat254==true) return SelectedArea(,eloss,charge);
    if(layerLabel==6 && NStrip==5 && NStripSat==2 && sat254==true) return SelectedArea(,eloss,charge);
    if(layerLabel==6 && NStrip==6 && NStripSat==1 && sat254==true) return SelectedArea(,eloss,charge);

    if(layerLabel==6 && NStrip==3 && NStripSat==1 && sat254==false) return SelectedArea(,eloss,charge);
    if(layerLabel==6 && NStrip==4 && NStripSat==1 && sat254==false) return SelectedArea(,eloss,charge);
    if(layerLabel==6 && NStrip==4 && NStripSat==2 && sat254==false) return SelectedArea(,eloss,charge);
    if(layerLabel==6 && NStrip==5 && NStripSat==1 && sat254==false) return SelectedArea(,eloss,charge);
    if(layerLabel==6 && NStrip==5 && NStripSat==2 && sat254==false) return SelectedArea(,eloss,charge);
    if(layerLabel==6 && NStrip==6 && NStripSat==1 && sat254==false) return SelectedArea(,eloss,charge);*/
}

void CalibCharge::SetParameters()
{
    string s_fit = "./data/fit_res.txt";

    ifstream ifile_fit(s_fit.c_str(),ios::in);

    if(ifile_fit)
    {
        while(!ifile_fit.eof())
        {
            ifile_fit >> TIB1_254_NStrip3_NStripSat1 >> p0_TIB1_254_NStrip3_NStripSat1 >> p1_TIB1_254_NStrip3_NStripSat1;
            ifile_fit >> TIB1_254_NStrip4_NStripSat1 >> p0_TIB1_254_NStrip4_NStripSat1 >> p1_TIB1_254_NStrip4_NStripSat1;
            ifile_fit >> TIB1_254_NStrip4_NStripSat2 >> p0_TIB1_254_NStrip4_NStripSat2 >> p1_TIB1_254_NStrip4_NStripSat2;
            ifile_fit >> TIB1_254_NStrip5_NStripSat1 >> p0_TIB1_254_NStrip5_NStripSat1 >> p1_TIB1_254_NStrip5_NStripSat1;
            ifile_fit >> TIB1_254_NStrip5_NStripSat2 >> p0_TIB1_254_NStrip5_NStripSat2 >> p1_TIB1_254_NStrip5_NStripSat2;
            ifile_fit >> TIB1_254_NStrip6_NStripSat1 >> p0_TIB1_254_NStrip6_NStripSat1 >> p1_TIB1_254_NStrip6_NStripSat1;

            ifile_fit >> TIB1_255_NStrip3_NStripSat1 >> p0_TIB1_255_NStrip3_NStripSat1 >> p1_TIB1_255_NStrip3_NStripSat1;
            ifile_fit >> TIB1_255_NStrip4_NStripSat1 >> p0_TIB1_255_NStrip4_NStripSat1 >> p1_TIB1_255_NStrip4_NStripSat1;
            ifile_fit >> TIB1_255_NStrip4_NStripSat2 >> p0_TIB1_255_NStrip4_NStripSat2 >> p1_TIB1_255_NStrip4_NStripSat2;
            ifile_fit >> TIB1_255_NStrip5_NStripSat1 >> p0_TIB1_255_NStrip5_NStripSat1 >> p1_TIB1_255_NStrip5_NStripSat1;
            ifile_fit >> TIB1_255_NStrip5_NStripSat2 >> p0_TIB1_255_NStrip5_NStripSat2 >> p1_TIB1_255_NStrip5_NStripSat2;
            ifile_fit >> TIB1_255_NStrip6_NStripSat1 >> p0_TIB1_255_NStrip6_NStripSat1 >> p1_TIB1_255_NStrip6_NStripSat1;

            ifile_fit >> TIB2_254_NStrip3_NStripSat1 >> p0_TIB2_254_NStrip3_NStripSat1 >> p1_TIB2_254_NStrip3_NStripSat1;
            ifile_fit >> TIB2_254_NStrip4_NStripSat1 >> p0_TIB2_254_NStrip4_NStripSat1 >> p1_TIB2_254_NStrip4_NStripSat1;
            ifile_fit >> TIB2_254_NStrip4_NStripSat2 >> p0_TIB2_254_NStrip4_NStripSat2 >> p1_TIB2_254_NStrip4_NStripSat2;
            ifile_fit >> TIB2_254_NStrip5_NStripSat1 >> p0_TIB2_254_NStrip5_NStripSat1 >> p1_TIB2_254_NStrip5_NStripSat1;
            ifile_fit >> TIB2_254_NStrip5_NStripSat2 >> p0_TIB2_254_NStrip5_NStripSat2 >> p1_TIB2_254_NStrip5_NStripSat2;
            ifile_fit >> TIB2_254_NStrip6_NStripSat1 >> p0_TIB2_254_NStrip6_NStripSat1 >> p1_TIB2_254_NStrip6_NStripSat1;

            ifile_fit >> TIB2_255_NStrip3_NStripSat1 >> p0_TIB2_255_NStrip3_NStripSat1 >> p1_TIB2_255_NStrip3_NStripSat1;
            ifile_fit >> TIB2_255_NStrip4_NStripSat1 >> p0_TIB2_255_NStrip4_NStripSat1 >> p1_TIB2_255_NStrip4_NStripSat1;
            ifile_fit >> TIB2_255_NStrip4_NStripSat2 >> p0_TIB2_255_NStrip4_NStripSat2 >> p1_TIB2_255_NStrip4_NStripSat2;
            ifile_fit >> TIB2_255_NStrip5_NStripSat1 >> p0_TIB2_255_NStrip5_NStripSat1 >> p1_TIB2_255_NStrip5_NStripSat1;
            ifile_fit >> TIB2_255_NStrip5_NStripSat2 >> p0_TIB2_255_NStrip5_NStripSat2 >> p1_TIB2_255_NStrip5_NStripSat2;
            ifile_fit >> TIB2_255_NStrip6_NStripSat1 >> p0_TIB2_255_NStrip6_NStripSat1 >> p1_TIB2_255_NStrip6_NStripSat1;

            ifile_fit >> TIB3_254_NStrip3_NStripSat1 >> p0_TIB3_254_NStrip3_NStripSat1 >> p1_TIB3_254_NStrip3_NStripSat1;
            ifile_fit >> TIB3_254_NStrip4_NStripSat1 >> p0_TIB3_254_NStrip4_NStripSat1 >> p1_TIB3_254_NStrip4_NStripSat1;
            ifile_fit >> TIB3_254_NStrip4_NStripSat2 >> p0_TIB3_254_NStrip4_NStripSat2 >> p1_TIB3_254_NStrip4_NStripSat2;
            ifile_fit >> TIB3_254_NStrip5_NStripSat1 >> p0_TIB3_254_NStrip5_NStripSat1 >> p1_TIB3_254_NStrip5_NStripSat1;
            ifile_fit >> TIB3_254_NStrip5_NStripSat2 >> p0_TIB3_254_NStrip5_NStripSat2 >> p1_TIB3_254_NStrip5_NStripSat2;
            ifile_fit >> TIB3_254_NStrip6_NStripSat1 >> p0_TIB3_254_NStrip6_NStripSat1 >> p1_TIB3_254_NStrip6_NStripSat1;

            ifile_fit >> TIB3_255_NStrip3_NStripSat1 >> p0_TIB3_255_NStrip3_NStripSat1 >> p1_TIB3_255_NStrip3_NStripSat1;
            ifile_fit >> TIB3_255_NStrip4_NStripSat1 >> p0_TIB3_255_NStrip4_NStripSat1 >> p1_TIB3_255_NStrip4_NStripSat1;
            ifile_fit >> TIB3_255_NStrip4_NStripSat2 >> p0_TIB3_255_NStrip4_NStripSat2 >> p1_TIB3_255_NStrip4_NStripSat2;
            ifile_fit >> TIB3_255_NStrip5_NStripSat1 >> p0_TIB3_255_NStrip5_NStripSat1 >> p1_TIB3_255_NStrip5_NStripSat1;
            ifile_fit >> TIB3_255_NStrip5_NStripSat2 >> p0_TIB3_255_NStrip5_NStripSat2 >> p1_TIB3_255_NStrip5_NStripSat2;
            ifile_fit >> TIB3_255_NStrip6_NStripSat1 >> p0_TIB3_255_NStrip6_NStripSat1 >> p1_TIB3_255_NStrip6_NStripSat1;

            ifile_fit >> TIB4_254_NStrip3_NStripSat1 >> p0_TIB4_254_NStrip3_NStripSat1 >> p1_TIB4_254_NStrip3_NStripSat1;
            ifile_fit >> TIB4_254_NStrip4_NStripSat1 >> p0_TIB4_254_NStrip4_NStripSat1 >> p1_TIB4_254_NStrip4_NStripSat1;
            ifile_fit >> TIB4_254_NStrip4_NStripSat2 >> p0_TIB4_254_NStrip4_NStripSat2 >> p1_TIB4_254_NStrip4_NStripSat2;
            ifile_fit >> TIB4_254_NStrip5_NStripSat1 >> p0_TIB4_254_NStrip5_NStripSat1 >> p1_TIB4_254_NStrip5_NStripSat1;
            ifile_fit >> TIB4_254_NStrip5_NStripSat2 >> p0_TIB4_254_NStrip5_NStripSat2 >> p1_TIB4_254_NStrip5_NStripSat2;
            ifile_fit >> TIB4_254_NStrip6_NStripSat1 >> p0_TIB4_254_NStrip6_NStripSat1 >> p1_TIB4_254_NStrip6_NStripSat1;

            ifile_fit >> TIB4_255_NStrip3_NStripSat1 >> p0_TIB4_255_NStrip3_NStripSat1 >> p1_TIB4_255_NStrip3_NStripSat1;
            ifile_fit >> TIB4_255_NStrip4_NStripSat1 >> p0_TIB4_255_NStrip4_NStripSat1 >> p1_TIB4_255_NStrip4_NStripSat1;
            ifile_fit >> TIB4_255_NStrip4_NStripSat2 >> p0_TIB4_255_NStrip4_NStripSat2 >> p1_TIB4_255_NStrip4_NStripSat2;
            ifile_fit >> TIB4_255_NStrip5_NStripSat1 >> p0_TIB4_255_NStrip5_NStripSat1 >> p1_TIB4_255_NStrip5_NStripSat1;
            ifile_fit >> TIB4_255_NStrip5_NStripSat2 >> p0_TIB4_255_NStrip5_NStripSat2 >> p1_TIB4_255_NStrip5_NStripSat2;
            ifile_fit >> TIB4_255_NStrip6_NStripSat1 >> p0_TIB4_255_NStrip6_NStripSat1 >> p1_TIB4_255_NStrip6_NStripSat1;

            ifile_fit >> TOB1_254_NStrip3_NStripSat1 >> p0_TOB1_254_NStrip3_NStripSat1 >> p1_TOB1_254_NStrip3_NStripSat1;
            ifile_fit >> TOB1_254_NStrip4_NStripSat1 >> p0_TOB1_254_NStrip4_NStripSat1 >> p1_TOB1_254_NStrip4_NStripSat1;
            ifile_fit >> TOB1_254_NStrip4_NStripSat2 >> p0_TOB1_254_NStrip4_NStripSat2 >> p1_TOB1_254_NStrip4_NStripSat2;
            ifile_fit >> TOB1_254_NStrip5_NStripSat1 >> p0_TOB1_254_NStrip5_NStripSat1 >> p1_TOB1_254_NStrip5_NStripSat1;
            ifile_fit >> TOB1_254_NStrip5_NStripSat2 >> p0_TOB1_254_NStrip5_NStripSat2 >> p1_TOB1_254_NStrip5_NStripSat2;
            ifile_fit >> TOB1_254_NStrip6_NStripSat1 >> p0_TOB1_254_NStrip6_NStripSat1 >> p1_TOB1_254_NStrip6_NStripSat1;

            ifile_fit >> TOB1_255_NStrip3_NStripSat1 >> p0_TOB1_255_NStrip3_NStripSat1 >> p1_TOB1_255_NStrip3_NStripSat1;
            ifile_fit >> TOB1_255_NStrip4_NStripSat1 >> p0_TOB1_255_NStrip4_NStripSat1 >> p1_TOB1_255_NStrip4_NStripSat1;
            ifile_fit >> TOB1_255_NStrip4_NStripSat2 >> p0_TOB1_255_NStrip4_NStripSat2 >> p1_TOB1_255_NStrip4_NStripSat2;
            ifile_fit >> TOB1_255_NStrip5_NStripSat1 >> p0_TOB1_255_NStrip5_NStripSat1 >> p1_TOB1_255_NStrip5_NStripSat1;
            ifile_fit >> TOB1_255_NStrip5_NStripSat2 >> p0_TOB1_255_NStrip5_NStripSat2 >> p1_TOB1_255_NStrip5_NStripSat2;
            ifile_fit >> TOB1_255_NStrip6_NStripSat1 >> p0_TOB1_255_NStrip6_NStripSat1 >> p1_TOB1_255_NStrip6_NStripSat1;

            ifile_fit >> TOB2_254_NStrip3_NStripSat1 >> p0_TOB2_254_NStrip3_NStripSat1 >> p1_TOB2_254_NStrip3_NStripSat1;
            ifile_fit >> TOB2_254_NStrip4_NStripSat1 >> p0_TOB2_254_NStrip4_NStripSat1 >> p1_TOB2_254_NStrip4_NStripSat1;
            ifile_fit >> TOB2_254_NStrip4_NStripSat2 >> p0_TOB2_254_NStrip4_NStripSat2 >> p1_TOB2_254_NStrip4_NStripSat2;
            ifile_fit >> TOB2_254_NStrip5_NStripSat1 >> p0_TOB2_254_NStrip5_NStripSat1 >> p1_TOB2_254_NStrip5_NStripSat1;
            ifile_fit >> TOB2_254_NStrip5_NStripSat2 >> p0_TOB2_254_NStrip5_NStripSat2 >> p1_TOB2_254_NStrip5_NStripSat2;
            ifile_fit >> TOB2_254_NStrip6_NStripSat1 >> p0_TOB2_254_NStrip6_NStripSat1 >> p1_TOB2_254_NStrip6_NStripSat1;

            ifile_fit >> TOB2_255_NStrip3_NStripSat1 >> p0_TOB2_255_NStrip3_NStripSat1 >> p1_TOB2_255_NStrip3_NStripSat1;
            ifile_fit >> TOB2_255_NStrip4_NStripSat1 >> p0_TOB2_255_NStrip4_NStripSat1 >> p1_TOB2_255_NStrip4_NStripSat1;
            ifile_fit >> TOB2_255_NStrip4_NStripSat2 >> p0_TOB2_255_NStrip4_NStripSat2 >> p1_TOB2_255_NStrip4_NStripSat2;
            ifile_fit >> TOB2_255_NStrip5_NStripSat1 >> p0_TOB2_255_NStrip5_NStripSat1 >> p1_TOB2_255_NStrip5_NStripSat1;
            ifile_fit >> TOB2_255_NStrip5_NStripSat2 >> p0_TOB2_255_NStrip5_NStripSat2 >> p1_TOB2_255_NStrip5_NStripSat2;
            ifile_fit >> TOB2_255_NStrip6_NStripSat1 >> p0_TOB2_255_NStrip6_NStripSat1 >> p1_TOB2_255_NStrip6_NStripSat1;

            ifile_fit >> TOB3_254_NStrip3_NStripSat1 >> p0_TOB3_254_NStrip3_NStripSat1 >> p1_TOB3_254_NStrip3_NStripSat1;
            ifile_fit >> TOB3_254_NStrip4_NStripSat1 >> p0_TOB3_254_NStrip4_NStripSat1 >> p1_TOB3_254_NStrip4_NStripSat1;
            ifile_fit >> TOB3_254_NStrip4_NStripSat2 >> p0_TOB3_254_NStrip4_NStripSat2 >> p1_TOB3_254_NStrip4_NStripSat2;
            ifile_fit >> TOB3_254_NStrip5_NStripSat1 >> p0_TOB3_254_NStrip5_NStripSat1 >> p1_TOB3_254_NStrip5_NStripSat1;
            ifile_fit >> TOB3_254_NStrip5_NStripSat2 >> p0_TOB3_254_NStrip5_NStripSat2 >> p1_TOB3_254_NStrip5_NStripSat2;
            ifile_fit >> TOB3_254_NStrip6_NStripSat1 >> p0_TOB3_254_NStrip6_NStripSat1 >> p1_TOB3_254_NStrip6_NStripSat1;

            ifile_fit >> TOB3_255_NStrip3_NStripSat1 >> p0_TOB3_255_NStrip3_NStripSat1 >> p1_TOB3_255_NStrip3_NStripSat1;
            ifile_fit >> TOB3_255_NStrip4_NStripSat1 >> p0_TOB3_255_NStrip4_NStripSat1 >> p1_TOB3_255_NStrip4_NStripSat1;
            ifile_fit >> TOB3_255_NStrip4_NStripSat2 >> p0_TOB3_255_NStrip4_NStripSat2 >> p1_TOB3_255_NStrip4_NStripSat2;
            ifile_fit >> TOB3_255_NStrip5_NStripSat1 >> p0_TOB3_255_NStrip5_NStripSat1 >> p1_TOB3_255_NStrip5_NStripSat1;
            ifile_fit >> TOB3_255_NStrip5_NStripSat2 >> p0_TOB3_255_NStrip5_NStripSat2 >> p1_TOB3_255_NStrip5_NStripSat2;
            ifile_fit >> TOB3_255_NStrip6_NStripSat1 >> p0_TOB3_255_NStrip6_NStripSat1 >> p1_TOB3_255_NStrip6_NStripSat1;

            ifile_fit >> TOB4_254_NStrip3_NStripSat1 >> p0_TOB4_254_NStrip3_NStripSat1 >> p1_TOB4_254_NStrip3_NStripSat1;
            ifile_fit >> TOB4_254_NStrip4_NStripSat1 >> p0_TOB4_254_NStrip4_NStripSat1 >> p1_TOB4_254_NStrip4_NStripSat1;
            ifile_fit >> TOB4_254_NStrip4_NStripSat2 >> p0_TOB4_254_NStrip4_NStripSat2 >> p1_TOB4_254_NStrip4_NStripSat2;
            ifile_fit >> TOB4_254_NStrip5_NStripSat1 >> p0_TOB4_254_NStrip5_NStripSat1 >> p1_TOB4_254_NStrip5_NStripSat1;
            ifile_fit >> TOB4_254_NStrip5_NStripSat2 >> p0_TOB4_254_NStrip5_NStripSat2 >> p1_TOB4_254_NStrip5_NStripSat2;
            ifile_fit >> TOB4_254_NStrip6_NStripSat1 >> p0_TOB4_254_NStrip6_NStripSat1 >> p1_TOB4_254_NStrip6_NStripSat1;

            ifile_fit >> TOB4_255_NStrip3_NStripSat1 >> p0_TOB4_255_NStrip3_NStripSat1 >> p1_TOB4_255_NStrip3_NStripSat1;
            ifile_fit >> TOB4_255_NStrip4_NStripSat1 >> p0_TOB4_255_NStrip4_NStripSat1 >> p1_TOB4_255_NStrip4_NStripSat1;
            ifile_fit >> TOB4_255_NStrip4_NStripSat2 >> p0_TOB4_255_NStrip4_NStripSat2 >> p1_TOB4_255_NStrip4_NStripSat2;
            ifile_fit >> TOB4_255_NStrip5_NStripSat1 >> p0_TOB4_255_NStrip5_NStripSat1 >> p1_TOB4_255_NStrip5_NStripSat1;
            ifile_fit >> TOB4_255_NStrip5_NStripSat2 >> p0_TOB4_255_NStrip5_NStripSat2 >> p1_TOB4_255_NStrip5_NStripSat2;
            ifile_fit >> TOB4_255_NStrip6_NStripSat1 >> p0_TOB4_255_NStrip6_NStripSat1 >> p1_TOB4_255_NStrip6_NStripSat1;

            ifile_fit >> TOB5_254_NStrip3_NStripSat1 >> p0_TOB5_254_NStrip3_NStripSat1 >> p1_TOB5_254_NStrip3_NStripSat1;
            ifile_fit >> TOB5_254_NStrip4_NStripSat1 >> p0_TOB5_254_NStrip4_NStripSat1 >> p1_TOB5_254_NStrip4_NStripSat1;
            ifile_fit >> TOB5_254_NStrip4_NStripSat2 >> p0_TOB5_254_NStrip4_NStripSat2 >> p1_TOB5_254_NStrip4_NStripSat2;
            ifile_fit >> TOB5_254_NStrip5_NStripSat1 >> p0_TOB5_254_NStrip5_NStripSat1 >> p1_TOB5_254_NStrip5_NStripSat1;
            ifile_fit >> TOB5_254_NStrip5_NStripSat2 >> p0_TOB5_254_NStrip5_NStripSat2 >> p1_TOB5_254_NStrip5_NStripSat2;
            ifile_fit >> TOB5_254_NStrip6_NStripSat1 >> p0_TOB5_254_NStrip6_NStripSat1 >> p1_TOB5_254_NStrip6_NStripSat1;

            ifile_fit >> TOB5_255_NStrip3_NStripSat1 >> p0_TOB5_255_NStrip3_NStripSat1 >> p1_TOB5_255_NStrip3_NStripSat1;
            ifile_fit >> TOB5_255_NStrip4_NStripSat1 >> p0_TOB5_255_NStrip4_NStripSat1 >> p1_TOB5_255_NStrip4_NStripSat1;
            ifile_fit >> TOB5_255_NStrip4_NStripSat2 >> p0_TOB5_255_NStrip4_NStripSat2 >> p1_TOB5_255_NStrip4_NStripSat2;
            ifile_fit >> TOB5_255_NStrip5_NStripSat1 >> p0_TOB5_255_NStrip5_NStripSat1 >> p1_TOB5_255_NStrip5_NStripSat1;
            ifile_fit >> TOB5_255_NStrip5_NStripSat2 >> p0_TOB5_255_NStrip5_NStripSat2 >> p1_TOB5_255_NStrip5_NStripSat2;
            ifile_fit >> TOB5_255_NStrip6_NStripSat1 >> p0_TOB5_255_NStrip6_NStripSat1 >> p1_TOB5_255_NStrip6_NStripSat1;

            ifile_fit >> TOB6_254_NStrip3_NStripSat1 >> p0_TOB6_254_NStrip3_NStripSat1 >> p1_TOB6_254_NStrip3_NStripSat1;
            ifile_fit >> TOB6_254_NStrip4_NStripSat1 >> p0_TOB6_254_NStrip4_NStripSat1 >> p1_TOB6_254_NStrip4_NStripSat1;
            ifile_fit >> TOB6_254_NStrip4_NStripSat2 >> p0_TOB6_254_NStrip4_NStripSat2 >> p1_TOB6_254_NStrip4_NStripSat2;
            ifile_fit >> TOB6_254_NStrip5_NStripSat1 >> p0_TOB6_254_NStrip5_NStripSat1 >> p1_TOB6_254_NStrip5_NStripSat1;
            ifile_fit >> TOB6_254_NStrip5_NStripSat2 >> p0_TOB6_254_NStrip5_NStripSat2 >> p1_TOB6_254_NStrip5_NStripSat2;
            ifile_fit >> TOB6_254_NStrip6_NStripSat1 >> p0_TOB6_254_NStrip6_NStripSat1 >> p1_TOB6_254_NStrip6_NStripSat1;

            ifile_fit >> TOB6_255_NStrip3_NStripSat1 >> p0_TOB6_255_NStrip3_NStripSat1 >> p1_TOB6_255_NStrip3_NStripSat1;
            ifile_fit >> TOB6_255_NStrip4_NStripSat1 >> p0_TOB6_255_NStrip4_NStripSat1 >> p1_TOB6_255_NStrip4_NStripSat1;
            ifile_fit >> TOB6_255_NStrip4_NStripSat2 >> p0_TOB6_255_NStrip4_NStripSat2 >> p1_TOB6_255_NStrip4_NStripSat2;
            ifile_fit >> TOB6_255_NStrip5_NStripSat1 >> p0_TOB6_255_NStrip5_NStripSat1 >> p1_TOB6_255_NStrip5_NStripSat1;
            ifile_fit >> TOB6_255_NStrip5_NStripSat2 >> p0_TOB6_255_NStrip5_NStripSat2 >> p1_TOB6_255_NStrip5_NStripSat2;
            ifile_fit >> TOB6_255_NStrip6_NStripSat1 >> p0_TOB6_255_NStrip6_NStripSat1 >> p1_TOB6_255_NStrip6_NStripSat1;
        }
        ifile_fit.close();
    }
    else cerr<<"Fichier de fit non ouvert"<<endl;
}

float CalibCharge::Charge(int layerLabel,int NStrip,float NStripSat,bool sat254,float charge)
{

    if(layerLabel==1 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TIB1_254_NStrip3_NStripSat1-charge)/p1_TIB1_254_NStrip3_NStripSat1;
    if(layerLabel==1 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TIB1_254_NStrip4_NStripSat1-charge)/p1_TIB1_254_NStrip4_NStripSat1;
    if(layerLabel==1 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TIB1_254_NStrip4_NStripSat2-charge)/p1_TIB1_254_NStrip4_NStripSat2;
    if(layerLabel==1 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TIB1_254_NStrip5_NStripSat1-charge)/p1_TIB1_254_NStrip5_NStripSat1;
    if(layerLabel==1 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TIB1_254_NStrip5_NStripSat2-charge)/p1_TIB1_254_NStrip5_NStripSat2;
    if(layerLabel==1 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TIB1_254_NStrip6_NStripSat1-charge)/p1_TIB1_254_NStrip6_NStripSat1;
    if(layerLabel==1 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TIB1_255_NStrip3_NStripSat1-charge)/p1_TIB1_255_NStrip3_NStripSat1;
    if(layerLabel==1 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TIB1_255_NStrip4_NStripSat1-charge)/p1_TIB1_255_NStrip4_NStripSat1;
    if(layerLabel==1 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TIB1_255_NStrip4_NStripSat2-charge)/p1_TIB1_255_NStrip4_NStripSat2;
    if(layerLabel==1 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TIB1_255_NStrip5_NStripSat1-charge)/p1_TIB1_255_NStrip5_NStripSat1;
    if(layerLabel==1 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TIB1_255_NStrip5_NStripSat2-charge)/p1_TIB1_255_NStrip5_NStripSat2;
    if(layerLabel==1 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TIB1_255_NStrip6_NStripSat1-charge)/p1_TIB1_255_NStrip6_NStripSat1;

    if(layerLabel==2 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TIB2_254_NStrip3_NStripSat1-charge)/p1_TIB2_254_NStrip3_NStripSat1;
    if(layerLabel==2 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TIB2_254_NStrip4_NStripSat1-charge)/p1_TIB2_254_NStrip4_NStripSat1;
    if(layerLabel==2 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TIB2_254_NStrip4_NStripSat2-charge)/p1_TIB2_254_NStrip4_NStripSat2;
    if(layerLabel==2 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TIB2_254_NStrip5_NStripSat1-charge)/p1_TIB2_254_NStrip5_NStripSat1;
    if(layerLabel==2 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TIB2_254_NStrip5_NStripSat2-charge)/p1_TIB2_254_NStrip5_NStripSat2;
    if(layerLabel==2 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TIB2_254_NStrip6_NStripSat1-charge)/p1_TIB2_254_NStrip6_NStripSat1;
    if(layerLabel==2 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TIB2_255_NStrip3_NStripSat1-charge)/p1_TIB2_255_NStrip3_NStripSat1;
    if(layerLabel==2 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TIB2_255_NStrip4_NStripSat1-charge)/p1_TIB2_255_NStrip4_NStripSat1;
    if(layerLabel==2 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TIB2_255_NStrip4_NStripSat2-charge)/p1_TIB2_255_NStrip4_NStripSat2;
    if(layerLabel==2 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TIB2_255_NStrip5_NStripSat1-charge)/p1_TIB2_255_NStrip5_NStripSat1;
    if(layerLabel==2 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TIB2_255_NStrip5_NStripSat2-charge)/p1_TIB2_255_NStrip5_NStripSat2;
    if(layerLabel==2 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TIB2_255_NStrip6_NStripSat1-charge)/p1_TIB2_255_NStrip6_NStripSat1;

    if(layerLabel==3 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TIB3_254_NStrip3_NStripSat1-charge)/p1_TIB3_254_NStrip3_NStripSat1;
    if(layerLabel==3 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TIB3_254_NStrip4_NStripSat1-charge)/p1_TIB3_254_NStrip4_NStripSat1;
    if(layerLabel==3 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TIB3_254_NStrip4_NStripSat2-charge)/p1_TIB3_254_NStrip4_NStripSat2;
    if(layerLabel==3 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TIB3_254_NStrip5_NStripSat1-charge)/p1_TIB3_254_NStrip5_NStripSat1;
    if(layerLabel==3 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TIB3_254_NStrip5_NStripSat2-charge)/p1_TIB3_254_NStrip5_NStripSat2;
    if(layerLabel==3 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TIB3_254_NStrip6_NStripSat1-charge)/p1_TIB3_254_NStrip6_NStripSat1;
    if(layerLabel==3 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TIB3_255_NStrip3_NStripSat1-charge)/p1_TIB3_255_NStrip3_NStripSat1;
    if(layerLabel==3 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TIB3_255_NStrip4_NStripSat1-charge)/p1_TIB3_255_NStrip4_NStripSat1;
    if(layerLabel==3 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TIB3_255_NStrip4_NStripSat2-charge)/p1_TIB3_255_NStrip4_NStripSat2;
    if(layerLabel==3 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TIB3_255_NStrip5_NStripSat1-charge)/p1_TIB3_255_NStrip5_NStripSat1;
    if(layerLabel==3 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TIB3_255_NStrip5_NStripSat2-charge)/p1_TIB3_255_NStrip5_NStripSat2;
    if(layerLabel==3 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TIB3_255_NStrip6_NStripSat1-charge)/p1_TIB3_255_NStrip6_NStripSat1;

    if(layerLabel==4 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TIB4_254_NStrip3_NStripSat1-charge)/p1_TIB4_254_NStrip3_NStripSat1;
    if(layerLabel==4 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TIB4_254_NStrip4_NStripSat1-charge)/p1_TIB4_254_NStrip4_NStripSat1;
    if(layerLabel==4 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TIB4_254_NStrip4_NStripSat2-charge)/p1_TIB4_254_NStrip4_NStripSat2;
    if(layerLabel==4 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TIB4_254_NStrip5_NStripSat1-charge)/p1_TIB4_254_NStrip5_NStripSat1;
    if(layerLabel==4 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TIB4_254_NStrip5_NStripSat2-charge)/p1_TIB4_254_NStrip5_NStripSat2;
    if(layerLabel==4 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TIB4_254_NStrip6_NStripSat1-charge)/p1_TIB4_254_NStrip6_NStripSat1;
    if(layerLabel==4 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TIB4_255_NStrip3_NStripSat1-charge)/p1_TIB4_255_NStrip3_NStripSat1;
    if(layerLabel==4 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TIB4_255_NStrip4_NStripSat1-charge)/p1_TIB4_255_NStrip4_NStripSat1;
    if(layerLabel==4 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TIB4_255_NStrip4_NStripSat2-charge)/p1_TIB4_255_NStrip4_NStripSat2;
    if(layerLabel==4 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TIB4_255_NStrip5_NStripSat1-charge)/p1_TIB4_255_NStrip5_NStripSat1;
    if(layerLabel==4 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TIB4_255_NStrip5_NStripSat2-charge)/p1_TIB4_255_NStrip5_NStripSat2;
    if(layerLabel==4 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TIB4_255_NStrip6_NStripSat1-charge)/p1_TIB4_255_NStrip6_NStripSat1;

    if(layerLabel==5 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TOB1_254_NStrip3_NStripSat1-charge)/p1_TOB1_254_NStrip3_NStripSat1;
    if(layerLabel==5 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TOB1_254_NStrip4_NStripSat1-charge)/p1_TOB1_254_NStrip4_NStripSat1;
    if(layerLabel==5 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TOB1_254_NStrip4_NStripSat2-charge)/p1_TOB1_254_NStrip4_NStripSat2;
    if(layerLabel==5 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TOB1_254_NStrip5_NStripSat1-charge)/p1_TOB1_254_NStrip5_NStripSat1;
    if(layerLabel==5 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TOB1_254_NStrip5_NStripSat2-charge)/p1_TOB1_254_NStrip5_NStripSat2;
    if(layerLabel==5 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TOB1_254_NStrip6_NStripSat1-charge)/p1_TOB1_254_NStrip6_NStripSat1;
    if(layerLabel==5 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TOB1_255_NStrip3_NStripSat1-charge)/p1_TOB1_255_NStrip3_NStripSat1;
    if(layerLabel==5 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TOB1_255_NStrip4_NStripSat1-charge)/p1_TOB1_255_NStrip4_NStripSat1;
    if(layerLabel==5 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TOB1_255_NStrip4_NStripSat2-charge)/p1_TOB1_255_NStrip4_NStripSat2;
    if(layerLabel==5 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TOB1_255_NStrip5_NStripSat1-charge)/p1_TOB1_255_NStrip5_NStripSat1;
    if(layerLabel==5 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TOB1_255_NStrip5_NStripSat2-charge)/p1_TOB1_255_NStrip5_NStripSat2;
    if(layerLabel==5 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TOB1_255_NStrip6_NStripSat1-charge)/p1_TOB1_255_NStrip6_NStripSat1;

    if(layerLabel==6 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TOB2_254_NStrip3_NStripSat1-charge)/p1_TOB2_254_NStrip3_NStripSat1;
    if(layerLabel==6 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TOB2_254_NStrip4_NStripSat1-charge)/p1_TOB2_254_NStrip4_NStripSat1;
    if(layerLabel==6 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TOB2_254_NStrip4_NStripSat2-charge)/p1_TOB2_254_NStrip4_NStripSat2;
    if(layerLabel==6 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TOB2_254_NStrip5_NStripSat1-charge)/p1_TOB2_254_NStrip5_NStripSat1;
    if(layerLabel==6 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TOB2_254_NStrip5_NStripSat2-charge)/p1_TOB2_254_NStrip5_NStripSat2;
    if(layerLabel==6 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TOB2_254_NStrip6_NStripSat1-charge)/p1_TOB2_254_NStrip6_NStripSat1;
    if(layerLabel==6 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TOB2_255_NStrip3_NStripSat1-charge)/p1_TOB2_255_NStrip3_NStripSat1;
    if(layerLabel==6 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TOB2_255_NStrip4_NStripSat1-charge)/p1_TOB2_255_NStrip4_NStripSat1;
    if(layerLabel==6 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TOB2_255_NStrip4_NStripSat2-charge)/p1_TOB2_255_NStrip4_NStripSat2;
    if(layerLabel==6 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TOB2_255_NStrip5_NStripSat1-charge)/p1_TOB2_255_NStrip5_NStripSat1;
    if(layerLabel==6 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TOB2_255_NStrip5_NStripSat2-charge)/p1_TOB2_255_NStrip5_NStripSat2;
    if(layerLabel==6 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TOB2_255_NStrip6_NStripSat1-charge)/p1_TOB2_255_NStrip6_NStripSat1;

    if(layerLabel==7 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TOB3_254_NStrip3_NStripSat1-charge)/p1_TOB3_254_NStrip3_NStripSat1;
    if(layerLabel==7 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TOB3_254_NStrip4_NStripSat1-charge)/p1_TOB3_254_NStrip4_NStripSat1;
    if(layerLabel==7 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TOB3_254_NStrip4_NStripSat2-charge)/p1_TOB3_254_NStrip4_NStripSat2;
    if(layerLabel==7 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TOB3_254_NStrip5_NStripSat1-charge)/p1_TOB3_254_NStrip5_NStripSat1;
    if(layerLabel==7 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TOB3_254_NStrip5_NStripSat2-charge)/p1_TOB3_254_NStrip5_NStripSat2;
    if(layerLabel==7 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TOB3_254_NStrip6_NStripSat1-charge)/p1_TOB3_254_NStrip6_NStripSat1;
    if(layerLabel==7 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TOB3_255_NStrip3_NStripSat1-charge)/p1_TOB3_255_NStrip3_NStripSat1;
    if(layerLabel==7 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TOB3_255_NStrip4_NStripSat1-charge)/p1_TOB3_255_NStrip4_NStripSat1;
    if(layerLabel==7 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TOB3_255_NStrip4_NStripSat2-charge)/p1_TOB3_255_NStrip4_NStripSat2;
    if(layerLabel==7 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TOB3_255_NStrip5_NStripSat1-charge)/p1_TOB3_255_NStrip5_NStripSat1;
    if(layerLabel==7 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TOB3_255_NStrip5_NStripSat2-charge)/p1_TOB3_255_NStrip5_NStripSat2;
    if(layerLabel==7 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TOB3_255_NStrip6_NStripSat1-charge)/p1_TOB3_255_NStrip6_NStripSat1;

    if(layerLabel==8 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TOB4_254_NStrip3_NStripSat1-charge)/p1_TOB4_254_NStrip3_NStripSat1;
    if(layerLabel==8 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TOB4_254_NStrip4_NStripSat1-charge)/p1_TOB4_254_NStrip4_NStripSat1;
    if(layerLabel==8 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TOB4_254_NStrip4_NStripSat2-charge)/p1_TOB4_254_NStrip4_NStripSat2;
    if(layerLabel==8 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TOB4_254_NStrip5_NStripSat1-charge)/p1_TOB4_254_NStrip5_NStripSat1;
    if(layerLabel==8 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TOB4_254_NStrip5_NStripSat2-charge)/p1_TOB4_254_NStrip5_NStripSat2;
    if(layerLabel==8 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TOB4_254_NStrip6_NStripSat1-charge)/p1_TOB4_254_NStrip6_NStripSat1;
    if(layerLabel==8 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TOB4_255_NStrip3_NStripSat1-charge)/p1_TOB4_255_NStrip3_NStripSat1;
    if(layerLabel==8 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TOB4_255_NStrip4_NStripSat1-charge)/p1_TOB4_255_NStrip4_NStripSat1;
    if(layerLabel==8 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TOB4_255_NStrip4_NStripSat2-charge)/p1_TOB4_255_NStrip4_NStripSat2;
    if(layerLabel==8 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TOB4_255_NStrip5_NStripSat1-charge)/p1_TOB4_255_NStrip5_NStripSat1;
    if(layerLabel==8 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TOB4_255_NStrip5_NStripSat2-charge)/p1_TOB4_255_NStrip5_NStripSat2;
    if(layerLabel==8 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TOB4_255_NStrip6_NStripSat1-charge)/p1_TOB4_255_NStrip6_NStripSat1;

    if(layerLabel==9 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TOB5_254_NStrip3_NStripSat1-charge)/p1_TOB5_254_NStrip3_NStripSat1;
    if(layerLabel==9 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TOB5_254_NStrip4_NStripSat1-charge)/p1_TOB5_254_NStrip4_NStripSat1;
    if(layerLabel==9 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TOB5_254_NStrip4_NStripSat2-charge)/p1_TOB5_254_NStrip4_NStripSat2;
    if(layerLabel==9 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TOB5_254_NStrip5_NStripSat1-charge)/p1_TOB5_254_NStrip5_NStripSat1;
    if(layerLabel==9 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TOB5_254_NStrip5_NStripSat2-charge)/p1_TOB5_254_NStrip5_NStripSat2;
    if(layerLabel==9 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TOB5_254_NStrip6_NStripSat1-charge)/p1_TOB5_254_NStrip6_NStripSat1;
    if(layerLabel==9 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TOB5_255_NStrip3_NStripSat1-charge)/p1_TOB5_255_NStrip3_NStripSat1;
    if(layerLabel==9 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TOB5_255_NStrip4_NStripSat1-charge)/p1_TOB5_255_NStrip4_NStripSat1;
    if(layerLabel==9 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TOB5_255_NStrip4_NStripSat2-charge)/p1_TOB5_255_NStrip4_NStripSat2;
    if(layerLabel==9 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TOB5_255_NStrip5_NStripSat1-charge)/p1_TOB5_255_NStrip5_NStripSat1;
    if(layerLabel==9 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TOB5_255_NStrip5_NStripSat2-charge)/p1_TOB5_255_NStrip5_NStripSat2;
    if(layerLabel==9 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TOB5_255_NStrip6_NStripSat1-charge)/p1_TOB5_255_NStrip6_NStripSat1;

    if(layerLabel==10 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TOB6_254_NStrip3_NStripSat1-charge)/p1_TOB6_254_NStrip3_NStripSat1;
    if(layerLabel==10 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TOB6_254_NStrip4_NStripSat1-charge)/p1_TOB6_254_NStrip4_NStripSat1;
    if(layerLabel==10 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TOB6_254_NStrip4_NStripSat2-charge)/p1_TOB6_254_NStrip4_NStripSat2;
    if(layerLabel==10 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TOB6_254_NStrip5_NStripSat1-charge)/p1_TOB6_254_NStrip5_NStripSat1;
    if(layerLabel==10 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TOB6_254_NStrip5_NStripSat2-charge)/p1_TOB6_254_NStrip5_NStripSat2;
    if(layerLabel==10 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TOB6_254_NStrip6_NStripSat1-charge)/p1_TOB6_254_NStrip6_NStripSat1;
    if(layerLabel==10 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TOB6_255_NStrip3_NStripSat1-charge)/p1_TOB6_255_NStrip3_NStripSat1;
    if(layerLabel==10 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TOB6_255_NStrip4_NStripSat1-charge)/p1_TOB6_255_NStrip4_NStripSat1;
    if(layerLabel==10 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TOB6_255_NStrip4_NStripSat2-charge)/p1_TOB6_255_NStrip4_NStripSat2;
    if(layerLabel==10 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TOB6_255_NStrip5_NStripSat1-charge)/p1_TOB6_255_NStrip5_NStripSat1;
    if(layerLabel==10 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TOB6_255_NStrip5_NStripSat2-charge)/p1_TOB6_255_NStrip5_NStripSat2;
    if(layerLabel==10 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TOB6_255_NStrip6_NStripSat1-charge)/p1_TOB6_255_NStrip6_NStripSat1;

    else return charge;
}