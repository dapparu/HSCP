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
    p0_TIB1_254_NStrip6_NStripSat2=0, p1_TIB1_254_NStrip6_NStripSat2=1;

    p0_TIB1_255_NStrip3_NStripSat1=0, p1_TIB1_255_NStrip3_NStripSat1=1;
    p0_TIB1_255_NStrip4_NStripSat1=0, p1_TIB1_255_NStrip4_NStripSat1=1;
    p0_TIB1_255_NStrip4_NStripSat2=0, p1_TIB1_255_NStrip4_NStripSat2=1;
    p0_TIB1_255_NStrip5_NStripSat1=0, p1_TIB1_255_NStrip5_NStripSat1=1;
    p0_TIB1_255_NStrip5_NStripSat2=0, p1_TIB1_255_NStrip5_NStripSat2=1;
    p0_TIB1_255_NStrip6_NStripSat1=0, p1_TIB1_255_NStrip6_NStripSat1=1;
    p0_TIB1_255_NStrip6_NStripSat2=0, p1_TIB1_255_NStrip6_NStripSat2=1;

    p0_TIB2_254_NStrip3_NStripSat1=0, p1_TIB2_254_NStrip3_NStripSat1=1;
    p0_TIB2_254_NStrip4_NStripSat1=0, p1_TIB2_254_NStrip4_NStripSat1=1;
    p0_TIB2_254_NStrip4_NStripSat2=0, p1_TIB2_254_NStrip4_NStripSat2=1;
    p0_TIB2_254_NStrip5_NStripSat1=0, p1_TIB2_254_NStrip5_NStripSat1=1;
    p0_TIB2_254_NStrip5_NStripSat2=0, p1_TIB2_254_NStrip5_NStripSat2=1;
    p0_TIB2_254_NStrip6_NStripSat1=0, p1_TIB2_254_NStrip6_NStripSat1=1;
    p0_TIB2_254_NStrip6_NStripSat2=0, p1_TIB2_254_NStrip6_NStripSat2=1;

    p0_TIB2_255_NStrip3_NStripSat1=0, p1_TIB2_255_NStrip3_NStripSat1=1;
    p0_TIB2_255_NStrip4_NStripSat1=0, p1_TIB2_255_NStrip4_NStripSat1=1;
    p0_TIB2_255_NStrip4_NStripSat2=0, p1_TIB2_255_NStrip4_NStripSat2=1;
    p0_TIB2_255_NStrip5_NStripSat1=0, p1_TIB2_255_NStrip5_NStripSat1=1;
    p0_TIB2_255_NStrip5_NStripSat2=0, p1_TIB2_255_NStrip5_NStripSat2=1;
    p0_TIB2_255_NStrip6_NStripSat1=0, p1_TIB2_255_NStrip6_NStripSat1=1;
    p0_TIB2_255_NStrip6_NStripSat2=0, p1_TIB2_255_NStrip6_NStripSat2=1;

    p0_TIB3_254_NStrip3_NStripSat1=0, p1_TIB3_254_NStrip3_NStripSat1=1;
    p0_TIB3_254_NStrip4_NStripSat1=0, p1_TIB3_254_NStrip4_NStripSat1=1;
    p0_TIB3_254_NStrip4_NStripSat2=0, p1_TIB3_254_NStrip4_NStripSat2=1;
    p0_TIB3_254_NStrip5_NStripSat1=0, p1_TIB3_254_NStrip5_NStripSat1=1;
    p0_TIB3_254_NStrip5_NStripSat2=0, p1_TIB3_254_NStrip5_NStripSat2=1;
    p0_TIB3_254_NStrip6_NStripSat1=0, p1_TIB3_254_NStrip6_NStripSat1=1;
    p0_TIB3_254_NStrip6_NStripSat2=0, p1_TIB3_254_NStrip6_NStripSat2=1;

    p0_TIB3_255_NStrip3_NStripSat1=0, p1_TIB3_255_NStrip3_NStripSat1=1;
    p0_TIB3_255_NStrip4_NStripSat1=0, p1_TIB3_255_NStrip4_NStripSat1=1;
    p0_TIB3_255_NStrip4_NStripSat2=0, p1_TIB3_255_NStrip4_NStripSat2=1;
    p0_TIB3_255_NStrip5_NStripSat1=0, p1_TIB3_255_NStrip5_NStripSat1=1;
    p0_TIB3_255_NStrip5_NStripSat2=0, p1_TIB3_255_NStrip5_NStripSat2=1;
    p0_TIB3_255_NStrip6_NStripSat1=0, p1_TIB3_255_NStrip6_NStripSat1=1;
    p0_TIB3_255_NStrip6_NStripSat2=0, p1_TIB3_255_NStrip6_NStripSat2=1;

    p0_TIB4_254_NStrip3_NStripSat1=0, p1_TIB4_254_NStrip3_NStripSat1=1;
    p0_TIB4_254_NStrip4_NStripSat1=0, p1_TIB4_254_NStrip4_NStripSat1=1;
    p0_TIB4_254_NStrip4_NStripSat2=0, p1_TIB4_254_NStrip4_NStripSat2=1;
    p0_TIB4_254_NStrip5_NStripSat1=0, p1_TIB4_254_NStrip5_NStripSat1=1;
    p0_TIB4_254_NStrip5_NStripSat2=0, p1_TIB4_254_NStrip5_NStripSat2=1;
    p0_TIB4_254_NStrip6_NStripSat1=0, p1_TIB4_254_NStrip6_NStripSat1=1;
    p0_TIB4_254_NStrip6_NStripSat2=0, p1_TIB4_254_NStrip6_NStripSat2=1;

    p0_TIB4_255_NStrip3_NStripSat1=0, p1_TIB4_255_NStrip3_NStripSat1=1;
    p0_TIB4_255_NStrip4_NStripSat1=0, p1_TIB4_255_NStrip4_NStripSat1=1;
    p0_TIB4_255_NStrip4_NStripSat2=0, p1_TIB4_255_NStrip4_NStripSat2=1;
    p0_TIB4_255_NStrip5_NStripSat1=0, p1_TIB4_255_NStrip5_NStripSat1=1;
    p0_TIB4_255_NStrip5_NStripSat2=0, p1_TIB4_255_NStrip5_NStripSat2=1;
    p0_TIB4_255_NStrip6_NStripSat1=0, p1_TIB4_255_NStrip6_NStripSat1=1;
    p0_TIB4_255_NStrip6_NStripSat2=0, p1_TIB4_255_NStrip6_NStripSat2=1;

    p0_TOB1_254_NStrip3_NStripSat1=0, p1_TOB1_254_NStrip3_NStripSat1=1;
    p0_TOB1_254_NStrip4_NStripSat1=0, p1_TOB1_254_NStrip4_NStripSat1=1;
    p0_TOB1_254_NStrip4_NStripSat2=0, p1_TOB1_254_NStrip4_NStripSat2=1;
    p0_TOB1_254_NStrip5_NStripSat1=0, p1_TOB1_254_NStrip5_NStripSat1=1;
    p0_TOB1_254_NStrip5_NStripSat2=0, p1_TOB1_254_NStrip5_NStripSat2=1;
    p0_TOB1_254_NStrip6_NStripSat1=0, p1_TOB1_254_NStrip6_NStripSat1=1;
    p0_TOB1_254_NStrip6_NStripSat2=0, p1_TOB1_254_NStrip6_NStripSat2=1;

    p0_TOB1_255_NStrip3_NStripSat1=0, p1_TOB1_255_NStrip3_NStripSat1=1;
    p0_TOB1_255_NStrip4_NStripSat1=0, p1_TOB1_255_NStrip4_NStripSat1=1;
    p0_TOB1_255_NStrip4_NStripSat2=0, p1_TOB1_255_NStrip4_NStripSat2=1;
    p0_TOB1_255_NStrip5_NStripSat1=0, p1_TOB1_255_NStrip5_NStripSat1=1;
    p0_TOB1_255_NStrip5_NStripSat2=0, p1_TOB1_255_NStrip5_NStripSat2=1;
    p0_TOB1_255_NStrip6_NStripSat1=0, p1_TOB1_255_NStrip6_NStripSat1=1;
    p0_TOB1_255_NStrip6_NStripSat2=0, p1_TOB1_255_NStrip6_NStripSat2=1;

    p0_TOB2_254_NStrip3_NStripSat1=0, p1_TOB2_254_NStrip3_NStripSat1=1;
    p0_TOB2_254_NStrip4_NStripSat1=0, p1_TOB2_254_NStrip4_NStripSat1=1;
    p0_TOB2_254_NStrip4_NStripSat2=0, p1_TOB2_254_NStrip4_NStripSat2=1;
    p0_TOB2_254_NStrip5_NStripSat1=0, p1_TOB2_254_NStrip5_NStripSat1=1;
    p0_TOB2_254_NStrip5_NStripSat2=0, p1_TOB2_254_NStrip5_NStripSat2=1;
    p0_TOB2_254_NStrip6_NStripSat1=0, p1_TOB2_254_NStrip6_NStripSat1=1;
    p0_TOB2_254_NStrip6_NStripSat2=0, p1_TOB2_254_NStrip6_NStripSat2=1;

    p0_TOB2_255_NStrip3_NStripSat1=0, p1_TOB2_255_NStrip3_NStripSat1=1;
    p0_TOB2_255_NStrip4_NStripSat1=0, p1_TOB2_255_NStrip4_NStripSat1=1;
    p0_TOB2_255_NStrip4_NStripSat2=0, p1_TOB2_255_NStrip4_NStripSat2=1;
    p0_TOB2_255_NStrip5_NStripSat1=0, p1_TOB2_255_NStrip5_NStripSat1=1;
    p0_TOB2_255_NStrip5_NStripSat2=0, p1_TOB2_255_NStrip5_NStripSat2=1;
    p0_TOB2_255_NStrip6_NStripSat1=0, p1_TOB2_255_NStrip6_NStripSat1=1;
    p0_TOB2_255_NStrip6_NStripSat2=0, p1_TOB2_255_NStrip6_NStripSat2=1;

    p0_TOB3_254_NStrip3_NStripSat1=0, p1_TOB3_254_NStrip3_NStripSat1=1;
    p0_TOB3_254_NStrip4_NStripSat1=0, p1_TOB3_254_NStrip4_NStripSat1=1;
    p0_TOB3_254_NStrip4_NStripSat2=0, p1_TOB3_254_NStrip4_NStripSat2=1;
    p0_TOB3_254_NStrip5_NStripSat1=0, p1_TOB3_254_NStrip5_NStripSat1=1;
    p0_TOB3_254_NStrip5_NStripSat2=0, p1_TOB3_254_NStrip5_NStripSat2=1;
    p0_TOB3_254_NStrip6_NStripSat1=0, p1_TOB3_254_NStrip6_NStripSat1=1;
    p0_TOB3_254_NStrip6_NStripSat2=0, p1_TOB3_254_NStrip6_NStripSat2=1;

    p0_TOB3_255_NStrip3_NStripSat1=0, p1_TOB3_255_NStrip3_NStripSat1=1;
    p0_TOB3_255_NStrip4_NStripSat1=0, p1_TOB3_255_NStrip4_NStripSat1=1;
    p0_TOB3_255_NStrip4_NStripSat2=0, p1_TOB3_255_NStrip4_NStripSat2=1;
    p0_TOB3_255_NStrip5_NStripSat1=0, p1_TOB3_255_NStrip5_NStripSat1=1;
    p0_TOB3_255_NStrip5_NStripSat2=0, p1_TOB3_255_NStrip5_NStripSat2=1;
    p0_TOB3_255_NStrip6_NStripSat1=0, p1_TOB3_255_NStrip6_NStripSat1=1;
    p0_TOB3_255_NStrip6_NStripSat2=0, p1_TOB3_255_NStrip6_NStripSat2=1;

    p0_TOB4_254_NStrip3_NStripSat1=0, p1_TOB4_254_NStrip3_NStripSat1=1;
    p0_TOB4_254_NStrip4_NStripSat1=0, p1_TOB4_254_NStrip4_NStripSat1=1;
    p0_TOB4_254_NStrip4_NStripSat2=0, p1_TOB4_254_NStrip4_NStripSat2=1;
    p0_TOB4_254_NStrip5_NStripSat1=0, p1_TOB4_254_NStrip5_NStripSat1=1;
    p0_TOB4_254_NStrip5_NStripSat2=0, p1_TOB4_254_NStrip5_NStripSat2=1;
    p0_TOB4_254_NStrip6_NStripSat1=0, p1_TOB4_254_NStrip6_NStripSat1=1;
    p0_TOB4_254_NStrip6_NStripSat2=0, p1_TOB4_254_NStrip6_NStripSat2=1;

    p0_TOB4_255_NStrip3_NStripSat1=0, p1_TOB4_255_NStrip3_NStripSat1=1;
    p0_TOB4_255_NStrip4_NStripSat1=0, p1_TOB4_255_NStrip4_NStripSat1=1;
    p0_TOB4_255_NStrip4_NStripSat2=0, p1_TOB4_255_NStrip4_NStripSat2=1;
    p0_TOB4_255_NStrip5_NStripSat1=0, p1_TOB4_255_NStrip5_NStripSat1=1;
    p0_TOB4_255_NStrip5_NStripSat2=0, p1_TOB4_255_NStrip5_NStripSat2=1;
    p0_TOB4_255_NStrip6_NStripSat1=0, p1_TOB4_255_NStrip6_NStripSat1=1;
    p0_TOB4_255_NStrip6_NStripSat2=0, p1_TOB4_255_NStrip6_NStripSat2=1;

    p0_TOB5_254_NStrip3_NStripSat1=0, p1_TOB5_254_NStrip3_NStripSat1=1;
    p0_TOB5_254_NStrip4_NStripSat1=0, p1_TOB5_254_NStrip4_NStripSat1=1;
    p0_TOB5_254_NStrip4_NStripSat2=0, p1_TOB5_254_NStrip4_NStripSat2=1;
    p0_TOB5_254_NStrip5_NStripSat1=0, p1_TOB5_254_NStrip5_NStripSat1=1;
    p0_TOB5_254_NStrip5_NStripSat2=0, p1_TOB5_254_NStrip5_NStripSat2=1;
    p0_TOB5_254_NStrip6_NStripSat1=0, p1_TOB5_254_NStrip6_NStripSat1=1;
    p0_TOB5_254_NStrip6_NStripSat2=0, p1_TOB5_254_NStrip6_NStripSat2=1;

    p0_TOB5_255_NStrip3_NStripSat1=0, p1_TOB5_255_NStrip3_NStripSat1=1;
    p0_TOB5_255_NStrip4_NStripSat1=0, p1_TOB5_255_NStrip4_NStripSat1=1;
    p0_TOB5_255_NStrip4_NStripSat2=0, p1_TOB5_255_NStrip4_NStripSat2=1;
    p0_TOB5_255_NStrip5_NStripSat1=0, p1_TOB5_255_NStrip5_NStripSat1=1;
    p0_TOB5_255_NStrip5_NStripSat2=0, p1_TOB5_255_NStrip5_NStripSat2=1;
    p0_TOB5_255_NStrip6_NStripSat1=0, p1_TOB5_255_NStrip6_NStripSat1=1;
    p0_TOB5_255_NStrip6_NStripSat2=0, p1_TOB5_255_NStrip6_NStripSat2=1;

    p0_TOB6_254_NStrip3_NStripSat1=0, p1_TOB6_254_NStrip3_NStripSat1=1;
    p0_TOB6_254_NStrip4_NStripSat1=0, p1_TOB6_254_NStrip4_NStripSat1=1;
    p0_TOB6_254_NStrip4_NStripSat2=0, p1_TOB6_254_NStrip4_NStripSat2=1;
    p0_TOB6_254_NStrip5_NStripSat1=0, p1_TOB6_254_NStrip5_NStripSat1=1;
    p0_TOB6_254_NStrip5_NStripSat2=0, p1_TOB6_254_NStrip5_NStripSat2=1;
    p0_TOB6_254_NStrip6_NStripSat1=0, p1_TOB6_254_NStrip6_NStripSat1=1;
    p0_TOB6_254_NStrip6_NStripSat2=0, p1_TOB6_254_NStrip6_NStripSat2=1;

    p0_TOB6_255_NStrip3_NStripSat1=0, p1_TOB6_255_NStrip3_NStripSat1=1;
    p0_TOB6_255_NStrip4_NStripSat1=0, p1_TOB6_255_NStrip4_NStripSat1=1;
    p0_TOB6_255_NStrip4_NStripSat2=0, p1_TOB6_255_NStrip4_NStripSat2=1;
    p0_TOB6_255_NStrip5_NStripSat1=0, p1_TOB6_255_NStrip5_NStripSat1=1;
    p0_TOB6_255_NStrip5_NStripSat2=0, p1_TOB6_255_NStrip5_NStripSat2=1;
    p0_TOB6_255_NStrip6_NStripSat1=0, p1_TOB6_255_NStrip6_NStripSat1=1;
    p0_TOB6_255_NStrip6_NStripSat2=0, p1_TOB6_255_NStrip6_NStripSat2=1;

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
    float orign1 = y1-coef1*x1;
    float orign2 = y4-coef2*x4;
    /*float coef3 = (y2-y1)/(x2-x1);
    float coef4 = (y4-y3)/(x4-x3);
    float orign3 = y2-coef3*x2;
    float orign4 = y3-coef4*x3;*/
    //if(x>=(y-orign3)/coef3 && x<=(y-orign4)/coef4 && y>=coef1*x+orign1 && y<=coef2*x+orign2) res=true;
    if(x>=x1 && x>=x2 && x<=x3 && x<=x4 && y>=coef1*x+orign1 && y<=coef2*x+orign2) res=true;
    return res;    
}

bool CalibCharge::Area(int layerLabel,int NStrip,int NStripSat,bool sat254,float eloss,float charge)
{
    if(layerLabel==1 && NStrip==3 && NStripSat==1 && sat254==true) return SelectedArea(245,261,245,293,595,332,595,365,eloss,charge);
    if(layerLabel==1 && NStrip==4 && NStripSat==1 && sat254==true) return SelectedArea(254,270,254,317,825,392,825,422,eloss,charge);
    if(layerLabel==1 && NStrip==4 && NStripSat==2 && sat254==true) return SelectedArea(501,497,501,530,1119,566,1119,604,eloss,charge);
    if(layerLabel==1 && NStrip==5 && NStripSat==1 && sat254==true) return SelectedArea(323,296,323,338,1131,452,1131,509,eloss,charge);
    if(layerLabel==1 && NStrip==5 && NStripSat==2 && sat254==true) return SelectedArea(545,503,545,560,1425,604,1425,649,eloss,charge);
    if(layerLabel==1 && NStrip==6 && NStripSat==1 && sat254==true) return SelectedArea(341,294,341,369,1069,444,1069,534,eloss,charge);
    if(layerLabel==1 && NStrip==6 && NStripSat==2 && sat254==true) return SelectedArea(678,525,678,579,1472,619,1472,680,eloss,charge);

    if(layerLabel==1 && NStrip==3 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==1 && NStrip==4 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==1 && NStrip==4 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==1 && NStrip==5 && NStripSat==1 && sat254==false) return SelectedArea(1036,438,1036,510,2176,665,2176,743,eloss,charge);
    if(layerLabel==1 && NStrip==5 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==1 && NStrip==6 && NStripSat==1 && sat254==false) return SelectedArea(1021,430,1021,564,2203,627,2203,809,eloss,charge);
    if(layerLabel==1 && NStrip==6 && NStripSat==2 && sat254==false) return SelectedArea(1036,438,1036,510,2176,665,2176,743,eloss,charge);

//---------------------------------------------------------------------------------------------

    if(layerLabel==2 && NStrip==3 && NStripSat==1 && sat254==true) return SelectedArea(242,261,242,296,618,332,618,368,eloss,charge);
    if(layerLabel==2 && NStrip==4 && NStripSat==1 && sat254==true) return SelectedArea(253,277,253,306,784,377,784,409,eloss,charge);
    if(layerLabel==2 && NStrip==4 && NStripSat==2 && sat254==true) return SelectedArea(477,500,477,533,1059,557,1059,590,eloss,charge);
    if(layerLabel==2 && NStrip==5 && NStripSat==1 && sat254==true) return SelectedArea(275,277,275,303,877,394,877,440,eloss,charge);
    if(layerLabel==2 && NStrip==5 && NStripSat==2 && sat254==true) return SelectedArea(530,506,530,548,1380,607,1380,652,eloss,charge);
    if(layerLabel==2 && NStrip==6 && NStripSat==1 && sat254==true) return SelectedArea(362,306,362,352,1148,470,1148,520,eloss,charge);
    if(layerLabel==2 && NStrip==6 && NStripSat==2 && sat254==true) return SelectedArea(688,521,688,585,1456,610,1456,690,eloss,charge);

    if(layerLabel==2 && NStrip==3 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==2 && NStrip==4 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==2 && NStrip==4 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==2 && NStrip==5 && NStripSat==1 && sat254==false) return SelectedArea(1000,426,1000,492,2300,665,2300,749,eloss,charge);
    if(layerLabel==2 && NStrip==5 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==2 && NStrip==6 && NStripSat==1 && sat254==false) return SelectedArea(916,424,916,570,2533,700,2533,900,eloss,charge);
    if(layerLabel==2 && NStrip==6 && NStripSat==1 && sat254==false) return SelectedArea(2139,645,2139,900,4100,670,4100,1000,eloss,charge);

//---------------------------------------------------------------------------------------------

    if(layerLabel==3 && NStrip==3 && NStripSat==1 && sat254==true) return SelectedArea(231,243,231,305,878,340,878,400,eloss,charge);
    if(layerLabel==3 && NStrip==4 && NStripSat==1 && sat254==true) return SelectedArea(258,258,258,293,820,353,820,380,eloss,charge);
    if(layerLabel==3 && NStrip==4 && NStripSat==2 && sat254==true) return SelectedArea(444,433,444,542,1222,530,1222,616,eloss,charge);
    if(layerLabel==3 && NStrip==5 && NStripSat==1 && sat254==true) return SelectedArea(251,258,251,309,1122,404,1122,446,eloss,charge);
    if(layerLabel==3 && NStrip==5 && NStripSat==2 && sat254==true) return SelectedArea(491,482,491,550,1400,550,1400,650,eloss,charge);
    if(layerLabel==3 && NStrip==6 && NStripSat==1 && sat254==true) return false;
    if(layerLabel==3 && NStrip==6 && NStripSat==2 && sat254==true) return false;

    if(layerLabel==3 && NStrip==3 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==3 && NStrip==4 && NStripSat==1 && sat254==false) return SelectedArea(860,290,860,480,1600,420,1600,620,eloss,charge);
    if(layerLabel==3 && NStrip==4 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==3 && NStrip==5 && NStripSat==1 && sat254==false) return SelectedArea(912,320,912,550,3100,540,3100,830,eloss,charge);
    if(layerLabel==3 && NStrip==5 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==3 && NStrip==6 && NStripSat==1 && sat254==false) return SelectedArea(912,320,912,550,3100,540,3100,830,eloss,charge);
    if(layerLabel==3 && NStrip==6 && NStripSat==2 && sat254==false) return false;

//---------------------------------------------------------------------------------------------

    if(layerLabel==4 && NStrip==3 && NStripSat==1 && sat254==true) return true;
    if(layerLabel==4 && NStrip==4 && NStripSat==1 && sat254==true) return SelectedArea(238,258,238,290,997,374,997,407,eloss,charge);
    if(layerLabel==4 && NStrip==4 && NStripSat==2 && sat254==true) return true;
    if(layerLabel==4 && NStrip==5 && NStripSat==1 && sat254==true) return SelectedArea(396,279,396,3323,1170,401,1170,440,eloss,charge);
    if(layerLabel==4 && NStrip==5 && NStripSat==2 && sat254==true) return true;
    if(layerLabel==4 && NStrip==6 && NStripSat==1 && sat254==true) return false;

    if(layerLabel==4 && NStrip==3 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==4 && NStrip==4 && NStripSat==1 && sat254==false) return true;
    if(layerLabel==4 && NStrip==4 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==4 && NStrip==5 && NStripSat==1 && sat254==false) return SelectedArea(990,354,990,432,2459,504,2459,593,eloss,charge);
    if(layerLabel==4 && NStrip==5 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==4 && NStrip==6 && NStripSat==1 && sat254==false) return false;

//---------------------------------------------------------------------------------------------

    if(layerLabel==5 && NStrip==3 && NStripSat==1 && sat254==true) return SelectedArea(221,252,221,294,618,339,618,386,eloss,charge);
    if(layerLabel==5 && NStrip==4 && NStripSat==1 && sat254==true) return SelectedArea(232,268,232,303,893,411,893,458,eloss,charge);
    if(layerLabel==5 && NStrip==4 && NStripSat==2 && sat254==true) return true;
    if(layerLabel==5 && NStrip==5 && NStripSat==1 && sat254==true) return SelectedArea(402,314,402,354,1297,542,1297,597,eloss,charge);
    if(layerLabel==5 && NStrip==5 && NStripSat==2 && sat254==true) return SelectedArea(556,509,556,575,1444,631,1444,685,eloss,charge);
    if(layerLabel==5 && NStrip==6 && NStripSat==1 && sat254==true) return SelectedArea(342,325,342,372,1309,535,1309,590,eloss,charge);

    if(layerLabel==5 && NStrip==3 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==5 && NStrip==4 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==5 && NStrip==4 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==5 && NStrip==5 && NStripSat==1 && sat254==false) return true;
    if(layerLabel==5 && NStrip==5 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==5 && NStrip==6 && NStripSat==1 && sat254==false) return true;

//---------------------------------------------------------------------------------------------

    if(layerLabel==6 && NStrip==3 && NStripSat==1 && sat254==true) return true;
    if(layerLabel==6 && NStrip==4 && NStripSat==1 && sat254==true) return SelectedArea(250,271,250,330,693,375,693,439,eloss,charge);
    if(layerLabel==6 && NStrip==4 && NStripSat==2 && sat254==true) return true;
    if(layerLabel==6 && NStrip==5 && NStripSat==1 && sat254==true) return SelectedArea(321,298,321,352,1195,502,1195,557,eloss,charge);
    if(layerLabel==6 && NStrip==5 && NStripSat==2 && sat254==true) return true;
    if(layerLabel==6 && NStrip==6 && NStripSat==1 && sat254==true) return SelectedArea(367,320,367,370,1101,502,1101,552,eloss,charge);

    if(layerLabel==6 && NStrip==3 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==6 && NStrip==4 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==6 && NStrip==4 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==6 && NStrip==5 && NStripSat==1 && sat254==false) return true;
    if(layerLabel==6 && NStrip==5 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==6 && NStrip==6 && NStripSat==1 && sat254==false) return true;

//---------------------------------------------------------------------------------------------

    if(layerLabel==7 && NStrip==3 && NStripSat==1 && sat254==true) return true;
    if(layerLabel==7 && NStrip==4 && NStripSat==1 && sat254==true) return SelectedArea(235,278,235,315,689,382,689,412,eloss,charge);
    if(layerLabel==7 && NStrip==4 && NStripSat==2 && sat254==true) return true;
    if(layerLabel==7 && NStrip==5 && NStripSat==1 && sat254==true) return SelectedArea(326,304,326,345,1130,490,1130,531,eloss,charge);
    if(layerLabel==7 && NStrip==5 && NStripSat==2 && sat254==true) return true;
    if(layerLabel==7 && NStrip==6 && NStripSat==1 && sat254==true) return false;

    if(layerLabel==7 && NStrip==3 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==7 && NStrip==4 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==7 && NStrip==4 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==7 && NStrip==5 && NStripSat==1 && sat254==false) return true;
    if(layerLabel==7 && NStrip==5 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==7 && NStrip==6 && NStripSat==1 && sat254==false) return false;

//---------------------------------------------------------------------------------------------

    if(layerLabel==8 && NStrip==3 && NStripSat==1 && sat254==true) return true;
    if(layerLabel==8 && NStrip==4 && NStripSat==1 && sat254==true) return SelectedArea(261,270,261,324,748,388,748,419,eloss,charge);
    if(layerLabel==8 && NStrip==4 && NStripSat==2 && sat254==true) return true;
    if(layerLabel==8 && NStrip==5 && NStripSat==1 && sat254==true) return SelectedArea(287,290,287,334,1130,480,1130,531,eloss,charge);
    if(layerLabel==8 && NStrip==5 && NStripSat==2 && sat254==true) return true;
    if(layerLabel==8 && NStrip==6 && NStripSat==1 && sat254==true) return false;

    if(layerLabel==8 && NStrip==3 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==8 && NStrip==4 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==8 && NStrip==4 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==8 && NStrip==5 && NStripSat==1 && sat254==false) return true;
    if(layerLabel==8 && NStrip==5 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==8 && NStrip==6 && NStripSat==1 && sat254==false) return true;

//---------------------------------------------------------------------------------------------

    if(layerLabel==9 && NStrip==3 && NStripSat==1 && sat254==true) return true;
    if(layerLabel==9 && NStrip==4 && NStripSat==1 && sat254==true) return SelectedArea(246,267,246,317,570,359,570,401,eloss,charge);
    if(layerLabel==9 && NStrip==4 && NStripSat==2 && sat254==true) return false;
    if(layerLabel==9 && NStrip==5 && NStripSat==1 && sat254==true) return SelectedArea(288,299,288,350,1082,497,1082,539,eloss,charge);
    if(layerLabel==9 && NStrip==5 && NStripSat==2 && sat254==true) return true;
    if(layerLabel==9 && NStrip==6 && NStripSat==1 && sat254==true) return false;

    if(layerLabel==9 && NStrip==3 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==9 && NStrip==4 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==9 && NStrip==4 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==9 && NStrip==5 && NStripSat==1 && sat254==false) return true;
    if(layerLabel==9 && NStrip==5 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==9 && NStrip==6 && NStripSat==1 && sat254==false) return false;  

//---------------------------------------------------------------------------------------------

    if(layerLabel==10 && NStrip==3 && NStripSat==1 && sat254==true) return true;
    if(layerLabel==10 && NStrip==4 && NStripSat==1 && sat254==true) return SelectedArea(250,287,250,320,561,353,561,398,eloss,charge);
    if(layerLabel==10 && NStrip==4 && NStripSat==2 && sat254==true) return false;
    if(layerLabel==10 && NStrip==5 && NStripSat==1 && sat254==true) return SelectedArea(290,302,290,347,967,470,967,509,eloss,charge);
    if(layerLabel==10 && NStrip==5 && NStripSat==2 && sat254==true) return true;
    if(layerLabel==10 && NStrip==6 && NStripSat==1 && sat254==true) return false;

    if(layerLabel==10 && NStrip==3 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==10 && NStrip==4 && NStripSat==1 && sat254==false) return false;
    if(layerLabel==10 && NStrip==4 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==10 && NStrip==5 && NStripSat==1 && sat254==false) return true;
    if(layerLabel==10 && NStrip==5 && NStripSat==2 && sat254==false) return false;
    if(layerLabel==10 && NStrip==6 && NStripSat==1 && sat254==false) return false;        

    else return false;
}

void CalibCharge::SetParameters()
{
    string s_fit = "./data/fit_res.txt";

    ifstream ifile_fit(s_fit.c_str(),ios::in);

    if(ifile_fit)
    {
        ifile_fit >> TIB1_254_NStrip3_NStripSat1 >> p0_TIB1_254_NStrip3_NStripSat1 >> p1_TIB1_254_NStrip3_NStripSat1;
        ifile_fit >> TIB1_254_NStrip4_NStripSat1 >> p0_TIB1_254_NStrip4_NStripSat1 >> p1_TIB1_254_NStrip4_NStripSat1;
        ifile_fit >> TIB1_254_NStrip4_NStripSat2 >> p0_TIB1_254_NStrip4_NStripSat2 >> p1_TIB1_254_NStrip4_NStripSat2;
        ifile_fit >> TIB1_254_NStrip5_NStripSat1 >> p0_TIB1_254_NStrip5_NStripSat1 >> p1_TIB1_254_NStrip5_NStripSat1;
        ifile_fit >> TIB1_254_NStrip5_NStripSat2 >> p0_TIB1_254_NStrip5_NStripSat2 >> p1_TIB1_254_NStrip5_NStripSat2;
        ifile_fit >> TIB1_254_NStrip6_NStripSat1 >> p0_TIB1_254_NStrip6_NStripSat1 >> p1_TIB1_254_NStrip6_NStripSat1;
        ifile_fit >> TIB1_254_NStrip6_NStripSat2 >> p0_TIB1_254_NStrip6_NStripSat2 >> p1_TIB1_254_NStrip6_NStripSat2;

        ifile_fit >> TIB1_255_NStrip3_NStripSat1 >> p0_TIB1_255_NStrip3_NStripSat1 >> p1_TIB1_255_NStrip3_NStripSat1;
        ifile_fit >> TIB1_255_NStrip4_NStripSat1 >> p0_TIB1_255_NStrip4_NStripSat1 >> p1_TIB1_255_NStrip4_NStripSat1;
        ifile_fit >> TIB1_255_NStrip4_NStripSat2 >> p0_TIB1_255_NStrip4_NStripSat2 >> p1_TIB1_255_NStrip4_NStripSat2;
        ifile_fit >> TIB1_255_NStrip5_NStripSat1 >> p0_TIB1_255_NStrip5_NStripSat1 >> p1_TIB1_255_NStrip5_NStripSat1;
        ifile_fit >> TIB1_255_NStrip5_NStripSat2 >> p0_TIB1_255_NStrip5_NStripSat2 >> p1_TIB1_255_NStrip5_NStripSat2;
        ifile_fit >> TIB1_255_NStrip6_NStripSat1 >> p0_TIB1_255_NStrip6_NStripSat1 >> p1_TIB1_255_NStrip6_NStripSat1;
        ifile_fit >> TIB1_255_NStrip6_NStripSat2 >> p0_TIB1_255_NStrip6_NStripSat2 >> p1_TIB1_255_NStrip6_NStripSat2;

        ifile_fit >> TIB2_254_NStrip3_NStripSat1 >> p0_TIB2_254_NStrip3_NStripSat1 >> p1_TIB2_254_NStrip3_NStripSat1;
        ifile_fit >> TIB2_254_NStrip4_NStripSat1 >> p0_TIB2_254_NStrip4_NStripSat1 >> p1_TIB2_254_NStrip4_NStripSat1;
        ifile_fit >> TIB2_254_NStrip4_NStripSat2 >> p0_TIB2_254_NStrip4_NStripSat2 >> p1_TIB2_254_NStrip4_NStripSat2;
        ifile_fit >> TIB2_254_NStrip5_NStripSat1 >> p0_TIB2_254_NStrip5_NStripSat1 >> p1_TIB2_254_NStrip5_NStripSat1;
        ifile_fit >> TIB2_254_NStrip5_NStripSat2 >> p0_TIB2_254_NStrip5_NStripSat2 >> p1_TIB2_254_NStrip5_NStripSat2;
        ifile_fit >> TIB2_254_NStrip6_NStripSat1 >> p0_TIB2_254_NStrip6_NStripSat1 >> p1_TIB2_254_NStrip6_NStripSat1;
        ifile_fit >> TIB2_254_NStrip6_NStripSat2 >> p0_TIB2_254_NStrip6_NStripSat2 >> p1_TIB2_254_NStrip6_NStripSat2;

        ifile_fit >> TIB2_255_NStrip3_NStripSat1 >> p0_TIB2_255_NStrip3_NStripSat1 >> p1_TIB2_255_NStrip3_NStripSat1;
        ifile_fit >> TIB2_255_NStrip4_NStripSat1 >> p0_TIB2_255_NStrip4_NStripSat1 >> p1_TIB2_255_NStrip4_NStripSat1;
        ifile_fit >> TIB2_255_NStrip4_NStripSat2 >> p0_TIB2_255_NStrip4_NStripSat2 >> p1_TIB2_255_NStrip4_NStripSat2;
        ifile_fit >> TIB2_255_NStrip5_NStripSat1 >> p0_TIB2_255_NStrip5_NStripSat1 >> p1_TIB2_255_NStrip5_NStripSat1;
        ifile_fit >> TIB2_255_NStrip5_NStripSat2 >> p0_TIB2_255_NStrip5_NStripSat2 >> p1_TIB2_255_NStrip5_NStripSat2;
        ifile_fit >> TIB2_255_NStrip6_NStripSat1 >> p0_TIB2_255_NStrip6_NStripSat1 >> p1_TIB2_255_NStrip6_NStripSat1;
        ifile_fit >> TIB2_255_NStrip6_NStripSat2 >> p0_TIB2_255_NStrip6_NStripSat2 >> p1_TIB2_255_NStrip6_NStripSat2;

        ifile_fit >> TIB3_254_NStrip3_NStripSat1 >> p0_TIB3_254_NStrip3_NStripSat1 >> p1_TIB3_254_NStrip3_NStripSat1;
        ifile_fit >> TIB3_254_NStrip4_NStripSat1 >> p0_TIB3_254_NStrip4_NStripSat1 >> p1_TIB3_254_NStrip4_NStripSat1;
        ifile_fit >> TIB3_254_NStrip4_NStripSat2 >> p0_TIB3_254_NStrip4_NStripSat2 >> p1_TIB3_254_NStrip4_NStripSat2;
        ifile_fit >> TIB3_254_NStrip5_NStripSat1 >> p0_TIB3_254_NStrip5_NStripSat1 >> p1_TIB3_254_NStrip5_NStripSat1;
        ifile_fit >> TIB3_254_NStrip5_NStripSat2 >> p0_TIB3_254_NStrip5_NStripSat2 >> p1_TIB3_254_NStrip5_NStripSat2;
        ifile_fit >> TIB3_254_NStrip6_NStripSat1 >> p0_TIB3_254_NStrip6_NStripSat1 >> p1_TIB3_254_NStrip6_NStripSat1;
        ifile_fit >> TIB3_254_NStrip6_NStripSat2 >> p0_TIB3_254_NStrip6_NStripSat2 >> p1_TIB3_254_NStrip6_NStripSat2;

        ifile_fit >> TIB3_255_NStrip3_NStripSat1 >> p0_TIB3_255_NStrip3_NStripSat1 >> p1_TIB3_255_NStrip3_NStripSat1;
        ifile_fit >> TIB3_255_NStrip4_NStripSat1 >> p0_TIB3_255_NStrip4_NStripSat1 >> p1_TIB3_255_NStrip4_NStripSat1;
        ifile_fit >> TIB3_255_NStrip4_NStripSat2 >> p0_TIB3_255_NStrip4_NStripSat2 >> p1_TIB3_255_NStrip4_NStripSat2;
        ifile_fit >> TIB3_255_NStrip5_NStripSat1 >> p0_TIB3_255_NStrip5_NStripSat1 >> p1_TIB3_255_NStrip5_NStripSat1;
        ifile_fit >> TIB3_255_NStrip5_NStripSat2 >> p0_TIB3_255_NStrip5_NStripSat2 >> p1_TIB3_255_NStrip5_NStripSat2;
        ifile_fit >> TIB3_255_NStrip6_NStripSat1 >> p0_TIB3_255_NStrip6_NStripSat1 >> p1_TIB3_255_NStrip6_NStripSat1;
        ifile_fit >> TIB3_255_NStrip6_NStripSat2 >> p0_TIB3_255_NStrip6_NStripSat2 >> p1_TIB3_255_NStrip6_NStripSat2;

        ifile_fit >> TIB4_254_NStrip3_NStripSat1 >> p0_TIB4_254_NStrip3_NStripSat1 >> p1_TIB4_254_NStrip3_NStripSat1;
        ifile_fit >> TIB4_254_NStrip4_NStripSat1 >> p0_TIB4_254_NStrip4_NStripSat1 >> p1_TIB4_254_NStrip4_NStripSat1;
        ifile_fit >> TIB4_254_NStrip4_NStripSat2 >> p0_TIB4_254_NStrip4_NStripSat2 >> p1_TIB4_254_NStrip4_NStripSat2;
        ifile_fit >> TIB4_254_NStrip5_NStripSat1 >> p0_TIB4_254_NStrip5_NStripSat1 >> p1_TIB4_254_NStrip5_NStripSat1;
        ifile_fit >> TIB4_254_NStrip5_NStripSat2 >> p0_TIB4_254_NStrip5_NStripSat2 >> p1_TIB4_254_NStrip5_NStripSat2;
        ifile_fit >> TIB4_254_NStrip6_NStripSat1 >> p0_TIB4_254_NStrip6_NStripSat1 >> p1_TIB4_254_NStrip6_NStripSat1;
        ifile_fit >> TIB4_254_NStrip6_NStripSat2 >> p0_TIB4_254_NStrip6_NStripSat2 >> p1_TIB4_254_NStrip6_NStripSat2;

        ifile_fit >> TIB4_255_NStrip3_NStripSat1 >> p0_TIB4_255_NStrip3_NStripSat1 >> p1_TIB4_255_NStrip3_NStripSat1;
        ifile_fit >> TIB4_255_NStrip4_NStripSat1 >> p0_TIB4_255_NStrip4_NStripSat1 >> p1_TIB4_255_NStrip4_NStripSat1;
        ifile_fit >> TIB4_255_NStrip4_NStripSat2 >> p0_TIB4_255_NStrip4_NStripSat2 >> p1_TIB4_255_NStrip4_NStripSat2;
        ifile_fit >> TIB4_255_NStrip5_NStripSat1 >> p0_TIB4_255_NStrip5_NStripSat1 >> p1_TIB4_255_NStrip5_NStripSat1;
        ifile_fit >> TIB4_255_NStrip5_NStripSat2 >> p0_TIB4_255_NStrip5_NStripSat2 >> p1_TIB4_255_NStrip5_NStripSat2;
        ifile_fit >> TIB4_255_NStrip6_NStripSat1 >> p0_TIB4_255_NStrip6_NStripSat1 >> p1_TIB4_255_NStrip6_NStripSat1;
        ifile_fit >> TIB4_255_NStrip6_NStripSat2 >> p0_TIB4_255_NStrip6_NStripSat2 >> p1_TIB4_255_NStrip6_NStripSat2;

        ifile_fit >> TOB1_254_NStrip3_NStripSat1 >> p0_TOB1_254_NStrip3_NStripSat1 >> p1_TOB1_254_NStrip3_NStripSat1;
        ifile_fit >> TOB1_254_NStrip4_NStripSat1 >> p0_TOB1_254_NStrip4_NStripSat1 >> p1_TOB1_254_NStrip4_NStripSat1;
        ifile_fit >> TOB1_254_NStrip4_NStripSat2 >> p0_TOB1_254_NStrip4_NStripSat2 >> p1_TOB1_254_NStrip4_NStripSat2;
        ifile_fit >> TOB1_254_NStrip5_NStripSat1 >> p0_TOB1_254_NStrip5_NStripSat1 >> p1_TOB1_254_NStrip5_NStripSat1;
        ifile_fit >> TOB1_254_NStrip5_NStripSat2 >> p0_TOB1_254_NStrip5_NStripSat2 >> p1_TOB1_254_NStrip5_NStripSat2;
        ifile_fit >> TOB1_254_NStrip6_NStripSat1 >> p0_TOB1_254_NStrip6_NStripSat1 >> p1_TOB1_254_NStrip6_NStripSat1;
        ifile_fit >> TOB1_254_NStrip6_NStripSat2 >> p0_TOB1_254_NStrip6_NStripSat2 >> p1_TOB1_254_NStrip6_NStripSat2;

        ifile_fit >> TOB1_255_NStrip3_NStripSat1 >> p0_TOB1_255_NStrip3_NStripSat1 >> p1_TOB1_255_NStrip3_NStripSat1;
        ifile_fit >> TOB1_255_NStrip4_NStripSat1 >> p0_TOB1_255_NStrip4_NStripSat1 >> p1_TOB1_255_NStrip4_NStripSat1;
        ifile_fit >> TOB1_255_NStrip4_NStripSat2 >> p0_TOB1_255_NStrip4_NStripSat2 >> p1_TOB1_255_NStrip4_NStripSat2;
        ifile_fit >> TOB1_255_NStrip5_NStripSat1 >> p0_TOB1_255_NStrip5_NStripSat1 >> p1_TOB1_255_NStrip5_NStripSat1;
        ifile_fit >> TOB1_255_NStrip5_NStripSat2 >> p0_TOB1_255_NStrip5_NStripSat2 >> p1_TOB1_255_NStrip5_NStripSat2;
        ifile_fit >> TOB1_255_NStrip6_NStripSat1 >> p0_TOB1_255_NStrip6_NStripSat1 >> p1_TOB1_255_NStrip6_NStripSat1;
        ifile_fit >> TOB1_255_NStrip6_NStripSat2 >> p0_TOB1_255_NStrip6_NStripSat2 >> p1_TOB1_255_NStrip6_NStripSat2;

        ifile_fit >> TOB2_254_NStrip3_NStripSat1 >> p0_TOB2_254_NStrip3_NStripSat1 >> p1_TOB2_254_NStrip3_NStripSat1;
        ifile_fit >> TOB2_254_NStrip4_NStripSat1 >> p0_TOB2_254_NStrip4_NStripSat1 >> p1_TOB2_254_NStrip4_NStripSat1;
        ifile_fit >> TOB2_254_NStrip4_NStripSat2 >> p0_TOB2_254_NStrip4_NStripSat2 >> p1_TOB2_254_NStrip4_NStripSat2;
        ifile_fit >> TOB2_254_NStrip5_NStripSat1 >> p0_TOB2_254_NStrip5_NStripSat1 >> p1_TOB2_254_NStrip5_NStripSat1;
        ifile_fit >> TOB2_254_NStrip5_NStripSat2 >> p0_TOB2_254_NStrip5_NStripSat2 >> p1_TOB2_254_NStrip5_NStripSat2;
        ifile_fit >> TOB2_254_NStrip6_NStripSat1 >> p0_TOB2_254_NStrip6_NStripSat1 >> p1_TOB2_254_NStrip6_NStripSat1;
        ifile_fit >> TOB2_254_NStrip6_NStripSat2 >> p0_TOB2_254_NStrip6_NStripSat2 >> p1_TOB2_254_NStrip6_NStripSat2;

        ifile_fit >> TOB2_255_NStrip3_NStripSat1 >> p0_TOB2_255_NStrip3_NStripSat1 >> p1_TOB2_255_NStrip3_NStripSat1;
        ifile_fit >> TOB2_255_NStrip4_NStripSat1 >> p0_TOB2_255_NStrip4_NStripSat1 >> p1_TOB2_255_NStrip4_NStripSat1;
        ifile_fit >> TOB2_255_NStrip4_NStripSat2 >> p0_TOB2_255_NStrip4_NStripSat2 >> p1_TOB2_255_NStrip4_NStripSat2;
        ifile_fit >> TOB2_255_NStrip5_NStripSat1 >> p0_TOB2_255_NStrip5_NStripSat1 >> p1_TOB2_255_NStrip5_NStripSat1;
        ifile_fit >> TOB2_255_NStrip5_NStripSat2 >> p0_TOB2_255_NStrip5_NStripSat2 >> p1_TOB2_255_NStrip5_NStripSat2;
        ifile_fit >> TOB2_255_NStrip6_NStripSat1 >> p0_TOB2_255_NStrip6_NStripSat1 >> p1_TOB2_255_NStrip6_NStripSat1;
        ifile_fit >> TOB2_255_NStrip6_NStripSat2 >> p0_TOB2_255_NStrip6_NStripSat2 >> p1_TOB2_255_NStrip6_NStripSat2;

        ifile_fit >> TOB3_254_NStrip3_NStripSat1 >> p0_TOB3_254_NStrip3_NStripSat1 >> p1_TOB3_254_NStrip3_NStripSat1;
        ifile_fit >> TOB3_254_NStrip4_NStripSat1 >> p0_TOB3_254_NStrip4_NStripSat1 >> p1_TOB3_254_NStrip4_NStripSat1;
        ifile_fit >> TOB3_254_NStrip4_NStripSat2 >> p0_TOB3_254_NStrip4_NStripSat2 >> p1_TOB3_254_NStrip4_NStripSat2;
        ifile_fit >> TOB3_254_NStrip5_NStripSat1 >> p0_TOB3_254_NStrip5_NStripSat1 >> p1_TOB3_254_NStrip5_NStripSat1;
        ifile_fit >> TOB3_254_NStrip5_NStripSat2 >> p0_TOB3_254_NStrip5_NStripSat2 >> p1_TOB3_254_NStrip5_NStripSat2;
        ifile_fit >> TOB3_254_NStrip6_NStripSat1 >> p0_TOB3_254_NStrip6_NStripSat1 >> p1_TOB3_254_NStrip6_NStripSat1;
        ifile_fit >> TOB3_254_NStrip6_NStripSat2 >> p0_TOB3_254_NStrip6_NStripSat2 >> p1_TOB3_254_NStrip6_NStripSat2;

        ifile_fit >> TOB3_255_NStrip3_NStripSat1 >> p0_TOB3_255_NStrip3_NStripSat1 >> p1_TOB3_255_NStrip3_NStripSat1;
        ifile_fit >> TOB3_255_NStrip4_NStripSat1 >> p0_TOB3_255_NStrip4_NStripSat1 >> p1_TOB3_255_NStrip4_NStripSat1;
        ifile_fit >> TOB3_255_NStrip4_NStripSat2 >> p0_TOB3_255_NStrip4_NStripSat2 >> p1_TOB3_255_NStrip4_NStripSat2;
        ifile_fit >> TOB3_255_NStrip5_NStripSat1 >> p0_TOB3_255_NStrip5_NStripSat1 >> p1_TOB3_255_NStrip5_NStripSat1;
        ifile_fit >> TOB3_255_NStrip5_NStripSat2 >> p0_TOB3_255_NStrip5_NStripSat2 >> p1_TOB3_255_NStrip5_NStripSat2;
        ifile_fit >> TOB3_255_NStrip6_NStripSat1 >> p0_TOB3_255_NStrip6_NStripSat1 >> p1_TOB3_255_NStrip6_NStripSat1;
        ifile_fit >> TOB3_255_NStrip6_NStripSat2 >> p0_TOB3_255_NStrip6_NStripSat2 >> p1_TOB3_255_NStrip6_NStripSat2;

        ifile_fit >> TOB4_254_NStrip3_NStripSat1 >> p0_TOB4_254_NStrip3_NStripSat1 >> p1_TOB4_254_NStrip3_NStripSat1;
        ifile_fit >> TOB4_254_NStrip4_NStripSat1 >> p0_TOB4_254_NStrip4_NStripSat1 >> p1_TOB4_254_NStrip4_NStripSat1;
        ifile_fit >> TOB4_254_NStrip4_NStripSat2 >> p0_TOB4_254_NStrip4_NStripSat2 >> p1_TOB4_254_NStrip4_NStripSat2;
        ifile_fit >> TOB4_254_NStrip5_NStripSat1 >> p0_TOB4_254_NStrip5_NStripSat1 >> p1_TOB4_254_NStrip5_NStripSat1;
        ifile_fit >> TOB4_254_NStrip5_NStripSat2 >> p0_TOB4_254_NStrip5_NStripSat2 >> p1_TOB4_254_NStrip5_NStripSat2;
        ifile_fit >> TOB4_254_NStrip6_NStripSat1 >> p0_TOB4_254_NStrip6_NStripSat1 >> p1_TOB4_254_NStrip6_NStripSat1;
        ifile_fit >> TOB4_254_NStrip6_NStripSat2 >> p0_TOB4_254_NStrip6_NStripSat2 >> p1_TOB4_254_NStrip6_NStripSat2;

        ifile_fit >> TOB4_255_NStrip3_NStripSat1 >> p0_TOB4_255_NStrip3_NStripSat1 >> p1_TOB4_255_NStrip3_NStripSat1;
        ifile_fit >> TOB4_255_NStrip4_NStripSat1 >> p0_TOB4_255_NStrip4_NStripSat1 >> p1_TOB4_255_NStrip4_NStripSat1;
        ifile_fit >> TOB4_255_NStrip4_NStripSat2 >> p0_TOB4_255_NStrip4_NStripSat2 >> p1_TOB4_255_NStrip4_NStripSat2;
        ifile_fit >> TOB4_255_NStrip5_NStripSat1 >> p0_TOB4_255_NStrip5_NStripSat1 >> p1_TOB4_255_NStrip5_NStripSat1;
        ifile_fit >> TOB4_255_NStrip5_NStripSat2 >> p0_TOB4_255_NStrip5_NStripSat2 >> p1_TOB4_255_NStrip5_NStripSat2;
        ifile_fit >> TOB4_255_NStrip6_NStripSat1 >> p0_TOB4_255_NStrip6_NStripSat1 >> p1_TOB4_255_NStrip6_NStripSat1;
        ifile_fit >> TOB4_255_NStrip6_NStripSat2 >> p0_TOB4_255_NStrip6_NStripSat2 >> p1_TOB4_255_NStrip6_NStripSat2;

        ifile_fit >> TOB5_254_NStrip3_NStripSat1 >> p0_TOB5_254_NStrip3_NStripSat1 >> p1_TOB5_254_NStrip3_NStripSat1;
        ifile_fit >> TOB5_254_NStrip4_NStripSat1 >> p0_TOB5_254_NStrip4_NStripSat1 >> p1_TOB5_254_NStrip4_NStripSat1;
        ifile_fit >> TOB5_254_NStrip4_NStripSat2 >> p0_TOB5_254_NStrip4_NStripSat2 >> p1_TOB5_254_NStrip4_NStripSat2;
        ifile_fit >> TOB5_254_NStrip5_NStripSat1 >> p0_TOB5_254_NStrip5_NStripSat1 >> p1_TOB5_254_NStrip5_NStripSat1;
        ifile_fit >> TOB5_254_NStrip5_NStripSat2 >> p0_TOB5_254_NStrip5_NStripSat2 >> p1_TOB5_254_NStrip5_NStripSat2;
        ifile_fit >> TOB5_254_NStrip6_NStripSat1 >> p0_TOB5_254_NStrip6_NStripSat1 >> p1_TOB5_254_NStrip6_NStripSat1;
        ifile_fit >> TOB5_254_NStrip6_NStripSat2 >> p0_TOB5_254_NStrip6_NStripSat2 >> p1_TOB5_254_NStrip6_NStripSat2;

        ifile_fit >> TOB5_255_NStrip3_NStripSat1 >> p0_TOB5_255_NStrip3_NStripSat1 >> p1_TOB5_255_NStrip3_NStripSat1;
        ifile_fit >> TOB5_255_NStrip4_NStripSat1 >> p0_TOB5_255_NStrip4_NStripSat1 >> p1_TOB5_255_NStrip4_NStripSat1;
        ifile_fit >> TOB5_255_NStrip4_NStripSat2 >> p0_TOB5_255_NStrip4_NStripSat2 >> p1_TOB5_255_NStrip4_NStripSat2;
        ifile_fit >> TOB5_255_NStrip5_NStripSat1 >> p0_TOB5_255_NStrip5_NStripSat1 >> p1_TOB5_255_NStrip5_NStripSat1;
        ifile_fit >> TOB5_255_NStrip5_NStripSat2 >> p0_TOB5_255_NStrip5_NStripSat2 >> p1_TOB5_255_NStrip5_NStripSat2;
        ifile_fit >> TOB5_255_NStrip6_NStripSat1 >> p0_TOB5_255_NStrip6_NStripSat1 >> p1_TOB5_255_NStrip6_NStripSat1;
        ifile_fit >> TOB5_255_NStrip6_NStripSat2 >> p0_TOB5_255_NStrip6_NStripSat2 >> p1_TOB5_255_NStrip6_NStripSat2;

        ifile_fit >> TOB6_254_NStrip3_NStripSat1 >> p0_TOB6_254_NStrip3_NStripSat1 >> p1_TOB6_254_NStrip3_NStripSat1;
        ifile_fit >> TOB6_254_NStrip4_NStripSat1 >> p0_TOB6_254_NStrip4_NStripSat1 >> p1_TOB6_254_NStrip4_NStripSat1;
        ifile_fit >> TOB6_254_NStrip4_NStripSat2 >> p0_TOB6_254_NStrip4_NStripSat2 >> p1_TOB6_254_NStrip4_NStripSat2;
        ifile_fit >> TOB6_254_NStrip5_NStripSat1 >> p0_TOB6_254_NStrip5_NStripSat1 >> p1_TOB6_254_NStrip5_NStripSat1;
        ifile_fit >> TOB6_254_NStrip5_NStripSat2 >> p0_TOB6_254_NStrip5_NStripSat2 >> p1_TOB6_254_NStrip5_NStripSat2;
        ifile_fit >> TOB6_254_NStrip6_NStripSat1 >> p0_TOB6_254_NStrip6_NStripSat1 >> p1_TOB6_254_NStrip6_NStripSat1;
        ifile_fit >> TOB6_254_NStrip6_NStripSat2 >> p0_TOB6_254_NStrip6_NStripSat2 >> p1_TOB6_254_NStrip6_NStripSat2;

        ifile_fit >> TOB6_255_NStrip3_NStripSat1 >> p0_TOB6_255_NStrip3_NStripSat1 >> p1_TOB6_255_NStrip3_NStripSat1;
        ifile_fit >> TOB6_255_NStrip4_NStripSat1 >> p0_TOB6_255_NStrip4_NStripSat1 >> p1_TOB6_255_NStrip4_NStripSat1;
        ifile_fit >> TOB6_255_NStrip4_NStripSat2 >> p0_TOB6_255_NStrip4_NStripSat2 >> p1_TOB6_255_NStrip4_NStripSat2;
        ifile_fit >> TOB6_255_NStrip5_NStripSat1 >> p0_TOB6_255_NStrip5_NStripSat1 >> p1_TOB6_255_NStrip5_NStripSat1;
        ifile_fit >> TOB6_255_NStrip5_NStripSat2 >> p0_TOB6_255_NStrip5_NStripSat2 >> p1_TOB6_255_NStrip5_NStripSat2;
        ifile_fit >> TOB6_255_NStrip6_NStripSat1 >> p0_TOB6_255_NStrip6_NStripSat1 >> p1_TOB6_255_NStrip6_NStripSat1;
        ifile_fit >> TOB6_255_NStrip6_NStripSat2 >> p0_TOB6_255_NStrip6_NStripSat2 >> p1_TOB6_255_NStrip6_NStripSat2;
    
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
    if(layerLabel==1 && NStrip>=6 && NStripSat==2 && sat254==true) return -(p0_TIB1_254_NStrip6_NStripSat2-charge)/p1_TIB1_254_NStrip6_NStripSat2;
    if(layerLabel==1 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TIB1_255_NStrip3_NStripSat1-charge)/p1_TIB1_255_NStrip3_NStripSat1;
    if(layerLabel==1 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TIB1_255_NStrip4_NStripSat1-charge)/p1_TIB1_255_NStrip4_NStripSat1;
    if(layerLabel==1 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TIB1_255_NStrip4_NStripSat2-charge)/p1_TIB1_255_NStrip4_NStripSat2;
    if(layerLabel==1 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TIB1_255_NStrip5_NStripSat1-charge)/p1_TIB1_255_NStrip5_NStripSat1;
    if(layerLabel==1 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TIB1_255_NStrip5_NStripSat2-charge)/p1_TIB1_255_NStrip5_NStripSat2;
    if(layerLabel==1 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TIB1_255_NStrip6_NStripSat1-charge)/p1_TIB1_255_NStrip6_NStripSat1;
    if(layerLabel==1 && NStrip>=6 && NStripSat==2 && sat254==false) return -(p0_TIB1_255_NStrip6_NStripSat2-charge)/p1_TIB1_255_NStrip6_NStripSat2;

    if(layerLabel==2 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TIB2_254_NStrip3_NStripSat1-charge)/p1_TIB2_254_NStrip3_NStripSat1;
    if(layerLabel==2 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TIB2_254_NStrip4_NStripSat1-charge)/p1_TIB2_254_NStrip4_NStripSat1;
    if(layerLabel==2 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TIB2_254_NStrip4_NStripSat2-charge)/p1_TIB2_254_NStrip4_NStripSat2;
    if(layerLabel==2 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TIB2_254_NStrip5_NStripSat1-charge)/p1_TIB2_254_NStrip5_NStripSat1;
    if(layerLabel==2 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TIB2_254_NStrip5_NStripSat2-charge)/p1_TIB2_254_NStrip5_NStripSat2;
    if(layerLabel==2 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TIB2_254_NStrip6_NStripSat1-charge)/p1_TIB2_254_NStrip6_NStripSat1;
    if(layerLabel==2 && NStrip>=6 && NStripSat==2 && sat254==true) return -(p0_TIB2_254_NStrip6_NStripSat2-charge)/p1_TIB2_254_NStrip6_NStripSat2;
    if(layerLabel==2 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TIB2_255_NStrip3_NStripSat1-charge)/p1_TIB2_255_NStrip3_NStripSat1;
    if(layerLabel==2 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TIB2_255_NStrip4_NStripSat1-charge)/p1_TIB2_255_NStrip4_NStripSat1;
    if(layerLabel==2 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TIB2_255_NStrip4_NStripSat2-charge)/p1_TIB2_255_NStrip4_NStripSat2;
    if(layerLabel==2 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TIB2_255_NStrip5_NStripSat1-charge)/p1_TIB2_255_NStrip5_NStripSat1;
    if(layerLabel==2 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TIB2_255_NStrip5_NStripSat2-charge)/p1_TIB2_255_NStrip5_NStripSat2;
    if(layerLabel==2 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TIB2_255_NStrip6_NStripSat1-charge)/p1_TIB2_255_NStrip6_NStripSat1;
    if(layerLabel==2 && NStrip>=6 && NStripSat==2 && sat254==false) return -(p0_TIB2_255_NStrip6_NStripSat2-charge)/p1_TIB2_255_NStrip6_NStripSat2;

    if(layerLabel==3 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TIB3_254_NStrip3_NStripSat1-charge)/p1_TIB3_254_NStrip3_NStripSat1;
    if(layerLabel==3 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TIB3_254_NStrip4_NStripSat1-charge)/p1_TIB3_254_NStrip4_NStripSat1;
    if(layerLabel==3 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TIB3_254_NStrip4_NStripSat2-charge)/p1_TIB3_254_NStrip4_NStripSat2;
    if(layerLabel==3 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TIB3_254_NStrip5_NStripSat1-charge)/p1_TIB3_254_NStrip5_NStripSat1;
    if(layerLabel==3 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TIB3_254_NStrip5_NStripSat2-charge)/p1_TIB3_254_NStrip5_NStripSat2;
    if(layerLabel==3 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TIB3_254_NStrip6_NStripSat1-charge)/p1_TIB3_254_NStrip6_NStripSat1;
    if(layerLabel==3 && NStrip>=6 && NStripSat==2 && sat254==true) return -(p0_TIB3_254_NStrip6_NStripSat2-charge)/p1_TIB3_254_NStrip6_NStripSat2;
    if(layerLabel==3 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TIB3_255_NStrip3_NStripSat1-charge)/p1_TIB3_255_NStrip3_NStripSat1;
    if(layerLabel==3 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TIB3_255_NStrip4_NStripSat1-charge)/p1_TIB3_255_NStrip4_NStripSat1;
    if(layerLabel==3 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TIB3_255_NStrip4_NStripSat2-charge)/p1_TIB3_255_NStrip4_NStripSat2;
    if(layerLabel==3 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TIB3_255_NStrip5_NStripSat1-charge)/p1_TIB3_255_NStrip5_NStripSat1;
    if(layerLabel==3 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TIB3_255_NStrip5_NStripSat2-charge)/p1_TIB3_255_NStrip5_NStripSat2;
    if(layerLabel==3 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TIB3_255_NStrip6_NStripSat1-charge)/p1_TIB3_255_NStrip6_NStripSat1;
    if(layerLabel==3 && NStrip>=6 && NStripSat==2 && sat254==false) return -(p0_TIB3_255_NStrip6_NStripSat2-charge)/p1_TIB3_255_NStrip6_NStripSat2;

    if(layerLabel==4 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TIB4_254_NStrip3_NStripSat1-charge)/p1_TIB4_254_NStrip3_NStripSat1;
    if(layerLabel==4 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TIB4_254_NStrip4_NStripSat1-charge)/p1_TIB4_254_NStrip4_NStripSat1;
    if(layerLabel==4 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TIB4_254_NStrip4_NStripSat2-charge)/p1_TIB4_254_NStrip4_NStripSat2;
    if(layerLabel==4 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TIB4_254_NStrip5_NStripSat1-charge)/p1_TIB4_254_NStrip5_NStripSat1;
    if(layerLabel==4 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TIB4_254_NStrip5_NStripSat2-charge)/p1_TIB4_254_NStrip5_NStripSat2;
    if(layerLabel==4 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TIB4_254_NStrip6_NStripSat1-charge)/p1_TIB4_254_NStrip6_NStripSat1;
    if(layerLabel==4 && NStrip>=6 && NStripSat==2 && sat254==true) return -(p0_TIB4_254_NStrip6_NStripSat2-charge)/p1_TIB4_254_NStrip6_NStripSat2;
    if(layerLabel==4 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TIB4_255_NStrip3_NStripSat1-charge)/p1_TIB4_255_NStrip3_NStripSat1;
    if(layerLabel==4 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TIB4_255_NStrip4_NStripSat1-charge)/p1_TIB4_255_NStrip4_NStripSat1;
    if(layerLabel==4 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TIB4_255_NStrip4_NStripSat2-charge)/p1_TIB4_255_NStrip4_NStripSat2;
    if(layerLabel==4 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TIB4_255_NStrip5_NStripSat1-charge)/p1_TIB4_255_NStrip5_NStripSat1;
    if(layerLabel==4 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TIB4_255_NStrip5_NStripSat2-charge)/p1_TIB4_255_NStrip5_NStripSat2;
    if(layerLabel==4 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TIB4_255_NStrip6_NStripSat1-charge)/p1_TIB4_255_NStrip6_NStripSat1;
    if(layerLabel==4 && NStrip>=6 && NStripSat==2 && sat254==false) return -(p0_TIB4_255_NStrip6_NStripSat2-charge)/p1_TIB4_255_NStrip6_NStripSat2;

    if(layerLabel==5 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TOB1_254_NStrip3_NStripSat1-charge)/p1_TOB1_254_NStrip3_NStripSat1;
    if(layerLabel==5 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TOB1_254_NStrip4_NStripSat1-charge)/p1_TOB1_254_NStrip4_NStripSat1;
    if(layerLabel==5 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TOB1_254_NStrip4_NStripSat2-charge)/p1_TOB1_254_NStrip4_NStripSat2;
    if(layerLabel==5 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TOB1_254_NStrip5_NStripSat1-charge)/p1_TOB1_254_NStrip5_NStripSat1;
    if(layerLabel==5 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TOB1_254_NStrip5_NStripSat2-charge)/p1_TOB1_254_NStrip5_NStripSat2;
    if(layerLabel==5 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TOB1_254_NStrip6_NStripSat1-charge)/p1_TOB1_254_NStrip6_NStripSat1;
    if(layerLabel==5 && NStrip>=6 && NStripSat==2 && sat254==true) return -(p0_TOB1_254_NStrip6_NStripSat2-charge)/p1_TOB1_254_NStrip6_NStripSat2;
    if(layerLabel==5 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TOB1_255_NStrip3_NStripSat1-charge)/p1_TOB1_255_NStrip3_NStripSat1;
    if(layerLabel==5 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TOB1_255_NStrip4_NStripSat1-charge)/p1_TOB1_255_NStrip4_NStripSat1;
    if(layerLabel==5 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TOB1_255_NStrip4_NStripSat2-charge)/p1_TOB1_255_NStrip4_NStripSat2;
    if(layerLabel==5 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TOB1_255_NStrip5_NStripSat1-charge)/p1_TOB1_255_NStrip5_NStripSat1;
    if(layerLabel==5 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TOB1_255_NStrip5_NStripSat2-charge)/p1_TOB1_255_NStrip5_NStripSat2;
    if(layerLabel==5 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TOB1_255_NStrip6_NStripSat1-charge)/p1_TOB1_255_NStrip6_NStripSat1;
    if(layerLabel==5 && NStrip>=6 && NStripSat==2 && sat254==false) return -(p0_TOB1_255_NStrip6_NStripSat2-charge)/p1_TOB1_255_NStrip6_NStripSat2;

    if(layerLabel==6 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TOB2_254_NStrip3_NStripSat1-charge)/p1_TOB2_254_NStrip3_NStripSat1;
    if(layerLabel==6 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TOB2_254_NStrip4_NStripSat1-charge)/p1_TOB2_254_NStrip4_NStripSat1;
    if(layerLabel==6 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TOB2_254_NStrip4_NStripSat2-charge)/p1_TOB2_254_NStrip4_NStripSat2;
    if(layerLabel==6 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TOB2_254_NStrip5_NStripSat1-charge)/p1_TOB2_254_NStrip5_NStripSat1;
    if(layerLabel==6 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TOB2_254_NStrip5_NStripSat2-charge)/p1_TOB2_254_NStrip5_NStripSat2;
    if(layerLabel==6 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TOB2_254_NStrip6_NStripSat1-charge)/p1_TOB2_254_NStrip6_NStripSat1;
    if(layerLabel==6 && NStrip>=6 && NStripSat==2 && sat254==true) return -(p0_TOB2_254_NStrip6_NStripSat2-charge)/p1_TOB2_254_NStrip6_NStripSat2;
    if(layerLabel==6 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TOB2_255_NStrip3_NStripSat1-charge)/p1_TOB2_255_NStrip3_NStripSat1;
    if(layerLabel==6 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TOB2_255_NStrip4_NStripSat1-charge)/p1_TOB2_255_NStrip4_NStripSat1;
    if(layerLabel==6 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TOB2_255_NStrip4_NStripSat2-charge)/p1_TOB2_255_NStrip4_NStripSat2;
    if(layerLabel==6 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TOB2_255_NStrip5_NStripSat1-charge)/p1_TOB2_255_NStrip5_NStripSat1;
    if(layerLabel==6 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TOB2_255_NStrip5_NStripSat2-charge)/p1_TOB2_255_NStrip5_NStripSat2;
    if(layerLabel==6 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TOB2_255_NStrip6_NStripSat1-charge)/p1_TOB2_255_NStrip6_NStripSat1;
    if(layerLabel==6 && NStrip>=6 && NStripSat==2 && sat254==false) return -(p0_TOB2_255_NStrip6_NStripSat2-charge)/p1_TOB2_255_NStrip6_NStripSat2;

    if(layerLabel==7 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TOB3_254_NStrip3_NStripSat1-charge)/p1_TOB3_254_NStrip3_NStripSat1;
    if(layerLabel==7 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TOB3_254_NStrip4_NStripSat1-charge)/p1_TOB3_254_NStrip4_NStripSat1;
    if(layerLabel==7 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TOB3_254_NStrip4_NStripSat2-charge)/p1_TOB3_254_NStrip4_NStripSat2;
    if(layerLabel==7 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TOB3_254_NStrip5_NStripSat1-charge)/p1_TOB3_254_NStrip5_NStripSat1;
    if(layerLabel==7 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TOB3_254_NStrip5_NStripSat2-charge)/p1_TOB3_254_NStrip5_NStripSat2;
    if(layerLabel==7 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TOB3_254_NStrip6_NStripSat1-charge)/p1_TOB3_254_NStrip6_NStripSat1;
    if(layerLabel==7 && NStrip>=6 && NStripSat==2 && sat254==true) return -(p0_TOB3_254_NStrip6_NStripSat2-charge)/p1_TOB3_254_NStrip6_NStripSat2;
    if(layerLabel==7 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TOB3_255_NStrip3_NStripSat1-charge)/p1_TOB3_255_NStrip3_NStripSat1;
    if(layerLabel==7 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TOB3_255_NStrip4_NStripSat1-charge)/p1_TOB3_255_NStrip4_NStripSat1;
    if(layerLabel==7 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TOB3_255_NStrip4_NStripSat2-charge)/p1_TOB3_255_NStrip4_NStripSat2;
    if(layerLabel==7 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TOB3_255_NStrip5_NStripSat1-charge)/p1_TOB3_255_NStrip5_NStripSat1;
    if(layerLabel==7 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TOB3_255_NStrip5_NStripSat2-charge)/p1_TOB3_255_NStrip5_NStripSat2;
    if(layerLabel==7 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TOB3_255_NStrip6_NStripSat1-charge)/p1_TOB3_255_NStrip6_NStripSat1;
    if(layerLabel==7 && NStrip>=6 && NStripSat==2 && sat254==false) return -(p0_TOB3_255_NStrip6_NStripSat2-charge)/p1_TOB3_255_NStrip6_NStripSat2;

    if(layerLabel==8 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TOB4_254_NStrip3_NStripSat1-charge)/p1_TOB4_254_NStrip3_NStripSat1;
    if(layerLabel==8 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TOB4_254_NStrip4_NStripSat1-charge)/p1_TOB4_254_NStrip4_NStripSat1;
    if(layerLabel==8 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TOB4_254_NStrip4_NStripSat2-charge)/p1_TOB4_254_NStrip4_NStripSat2;
    if(layerLabel==8 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TOB4_254_NStrip5_NStripSat1-charge)/p1_TOB4_254_NStrip5_NStripSat1;
    if(layerLabel==8 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TOB4_254_NStrip5_NStripSat2-charge)/p1_TOB4_254_NStrip5_NStripSat2;
    if(layerLabel==8 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TOB4_254_NStrip6_NStripSat1-charge)/p1_TOB4_254_NStrip6_NStripSat1;
    if(layerLabel==8 && NStrip>=6 && NStripSat==2 && sat254==true) return -(p0_TOB4_254_NStrip6_NStripSat2-charge)/p1_TOB4_254_NStrip6_NStripSat2;
    if(layerLabel==8 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TOB4_255_NStrip3_NStripSat1-charge)/p1_TOB4_255_NStrip3_NStripSat1;
    if(layerLabel==8 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TOB4_255_NStrip4_NStripSat1-charge)/p1_TOB4_255_NStrip4_NStripSat1;
    if(layerLabel==8 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TOB4_255_NStrip4_NStripSat2-charge)/p1_TOB4_255_NStrip4_NStripSat2;
    if(layerLabel==8 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TOB4_255_NStrip5_NStripSat1-charge)/p1_TOB4_255_NStrip5_NStripSat1;
    if(layerLabel==8 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TOB4_255_NStrip5_NStripSat2-charge)/p1_TOB4_255_NStrip5_NStripSat2;
    if(layerLabel==8 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TOB4_255_NStrip6_NStripSat1-charge)/p1_TOB4_255_NStrip6_NStripSat1;
    if(layerLabel==8 && NStrip>=6 && NStripSat==2 && sat254==false) return -(p0_TOB4_255_NStrip6_NStripSat2-charge)/p1_TOB4_255_NStrip6_NStripSat2;

    if(layerLabel==9 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TOB5_254_NStrip3_NStripSat1-charge)/p1_TOB5_254_NStrip3_NStripSat1;
    if(layerLabel==9 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TOB5_254_NStrip4_NStripSat1-charge)/p1_TOB5_254_NStrip4_NStripSat1;
    if(layerLabel==9 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TOB5_254_NStrip4_NStripSat2-charge)/p1_TOB5_254_NStrip4_NStripSat2;
    if(layerLabel==9 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TOB5_254_NStrip5_NStripSat1-charge)/p1_TOB5_254_NStrip5_NStripSat1;
    if(layerLabel==9 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TOB5_254_NStrip5_NStripSat2-charge)/p1_TOB5_254_NStrip5_NStripSat2;
    if(layerLabel==9 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TOB5_254_NStrip6_NStripSat1-charge)/p1_TOB5_254_NStrip6_NStripSat1;
    if(layerLabel==9 && NStrip>=6 && NStripSat==2 && sat254==true) return -(p0_TOB5_254_NStrip6_NStripSat2-charge)/p1_TOB5_254_NStrip6_NStripSat2;
    if(layerLabel==9 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TOB5_255_NStrip3_NStripSat1-charge)/p1_TOB5_255_NStrip3_NStripSat1;
    if(layerLabel==9 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TOB5_255_NStrip4_NStripSat1-charge)/p1_TOB5_255_NStrip4_NStripSat1;
    if(layerLabel==9 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TOB5_255_NStrip4_NStripSat2-charge)/p1_TOB5_255_NStrip4_NStripSat2;
    if(layerLabel==9 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TOB5_255_NStrip5_NStripSat1-charge)/p1_TOB5_255_NStrip5_NStripSat1;
    if(layerLabel==9 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TOB5_255_NStrip5_NStripSat2-charge)/p1_TOB5_255_NStrip5_NStripSat2;
    if(layerLabel==9 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TOB5_255_NStrip6_NStripSat1-charge)/p1_TOB5_255_NStrip6_NStripSat1;
    if(layerLabel==9 && NStrip>=6 && NStripSat==2 && sat254==false) return -(p0_TOB5_255_NStrip6_NStripSat2-charge)/p1_TOB5_255_NStrip6_NStripSat2;

    if(layerLabel==10 && NStrip==3 && NStripSat==1 && sat254==true) return -(p0_TOB6_254_NStrip3_NStripSat1-charge)/p1_TOB6_254_NStrip3_NStripSat1;
    if(layerLabel==10 && NStrip==4 && NStripSat==1 && sat254==true) return -(p0_TOB6_254_NStrip4_NStripSat1-charge)/p1_TOB6_254_NStrip4_NStripSat1;
    if(layerLabel==10 && NStrip==4 && NStripSat==2 && sat254==true) return -(p0_TOB6_254_NStrip4_NStripSat2-charge)/p1_TOB6_254_NStrip4_NStripSat2;
    if(layerLabel==10 && NStrip==5 && NStripSat==1 && sat254==true) return -(p0_TOB6_254_NStrip5_NStripSat1-charge)/p1_TOB6_254_NStrip5_NStripSat1;
    if(layerLabel==10 && NStrip==5 && NStripSat==2 && sat254==true) return -(p0_TOB6_254_NStrip5_NStripSat2-charge)/p1_TOB6_254_NStrip5_NStripSat2;
    if(layerLabel==10 && NStrip>=6 && NStripSat==1 && sat254==true) return -(p0_TOB6_254_NStrip6_NStripSat1-charge)/p1_TOB6_254_NStrip6_NStripSat1;
    if(layerLabel==10 && NStrip>=6 && NStripSat==2 && sat254==true) return -(p0_TOB6_254_NStrip6_NStripSat2-charge)/p1_TOB6_254_NStrip6_NStripSat2;
    if(layerLabel==10 && NStrip==3 && NStripSat==1 && sat254==false) return -(p0_TOB6_255_NStrip3_NStripSat1-charge)/p1_TOB6_255_NStrip3_NStripSat1;
    if(layerLabel==10 && NStrip==4 && NStripSat==1 && sat254==false) return -(p0_TOB6_255_NStrip4_NStripSat1-charge)/p1_TOB6_255_NStrip4_NStripSat1;
    if(layerLabel==10 && NStrip==4 && NStripSat==2 && sat254==false) return -(p0_TOB6_255_NStrip4_NStripSat2-charge)/p1_TOB6_255_NStrip4_NStripSat2;
    if(layerLabel==10 && NStrip==5 && NStripSat==1 && sat254==false) return -(p0_TOB6_255_NStrip5_NStripSat1-charge)/p1_TOB6_255_NStrip5_NStripSat1;
    if(layerLabel==10 && NStrip==5 && NStripSat==2 && sat254==false) return -(p0_TOB6_255_NStrip5_NStripSat2-charge)/p1_TOB6_255_NStrip5_NStripSat2;
    if(layerLabel==10 && NStrip>=6 && NStripSat==1 && sat254==false) return -(p0_TOB6_255_NStrip6_NStripSat1-charge)/p1_TOB6_255_NStrip6_NStripSat1;
    if(layerLabel==10 && NStrip>=6 && NStripSat==2 && sat254==false) return -(p0_TOB6_255_NStrip6_NStripSat2-charge)/p1_TOB6_255_NStrip6_NStripSat2;

    else return charge;
}