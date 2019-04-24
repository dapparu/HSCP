#include "../interface/Cluster.h"

using namespace std;

Cluster::Cluster()
{
	dedx_charge_		= .0;
	sclus_charge_		= .0;
	pathlength_			= .0;
	eloss_ 				= .0;
	nstrips_			= 0;
	nsimhits_			= 0;
	detid_				= 0;
	subdetid_			= 0;
	firstsclus_			= 0;
	sat254_				= false;
	sat255_				= false;
	shape_				= false;
	partid_				= 0;
}

Cluster::Cluster(float dedx_charge,float sclus_charge,float pathlength,float eloss,int nstrips,int nsimhits,int detid,int subdetid,bool sat254,bool sat255,bool shape,int firstsclus,int partId,const vector<ClusterStrip> &VectStrips,const vector<SimHit> &VectSimHits)
{
	dedx_charge_		= dedx_charge;
	sclus_charge_		= sclus_charge;
	pathlength_			= pathlength;
	eloss_				= eloss;
	nstrips_			= nstrips;
	nsimhits_			= nsimhits;
	detid_				= detid;
	subdetid_			= subdetid;
	sat254_				= sat254;
	sat255_				= sat255;
	shape_				= shape;
	firstsclus_			= firstsclus;
	VectStrips_			= VectStrips;
	VectSimHits_ 		= VectSimHits;
	partid_				= partId;
}

Cluster::~Cluster()
{
}

float Cluster::GetDedxCharge() const
{
	return dedx_charge_;
}

float Cluster::GetSclusCharge() const
{
	return sclus_charge_;
}

float Cluster::GetPathLength() const
{
	return pathlength_;
}

float Cluster::GetEloss() const
{
	return eloss_;
}

int Cluster::GetFirstsclus() const
{
	return firstsclus_;
}

int Cluster::GetLayer() const
{
	if(subdetid_==3) return 3*10+((detid_>>14)&0x7)-1;			//TIB
	else if(subdetid_==4) return 4*10+((detid_>>14)&0x3);		//TID
	else if(subdetid_==5) return 5*10+((detid_>>14)&0x7)-1;		//TOB
	else if(subdetid_==6) return 6*10+((detid_>>14)&0xF)-1;		//TEC
}

int Cluster::GetLayerLabel() const
{
	if(subdetid_==3)
	{
		if(((detid_>>14)&0x7)==1) return 1;
		else if(((detid_>>14)&0x7)==2) return 2;
		else if(((detid_>>14)&0x7)==3) return 3;
		else if(((detid_>>14)&0x7)==4) return 4; 
	}
	else if(subdetid_==5)
	{
		if(((detid_>>14)&0x7)==1) return 5;
		else if(((detid_>>14)&0x7)==2) return 6;
		else if(((detid_>>14)&0x7)==3) return 7;
		else if(((detid_>>14)&0x7)==4) return 8;
		else if(((detid_>>14)&0x7)==5) return 9;
		else if(((detid_>>14)&0x7)==6) return 10;
	}
	else if(subdetid_==4)
	{
		if(((detid_>>14)&0x3)==0) return 11;
		else if(((detid_>>14)&0x3)==1) return 12;
	}
	else if(subdetid_==6)
	{
		if(((detid_>>14)&0xF)==1) return 13;
		else if(((detid_>>14)&0xF)==2) return 14;
		else if(((detid_>>14)&0xF)==3) return 15;
		else if(((detid_>>14)&0xF)==4) return 16;
		else if(((detid_>>14)&0xF)==5) return 17;
		else if(((detid_>>14)&0xF)==6) return 18;
		else if(((detid_>>14)&0xF)==7) return 19;
		else if(((detid_>>14)&0xF)==8) return 20;
		else if(((detid_>>14)&0xF)==9) return 21;
	}
}

int Cluster::GetNStripLayer(int layerlabel) const
{
	if(layerlabel==1 || layerlabel==2 || layerlabel==9 || layerlabel==10) return 768;
	if(layerlabel==3 || layerlabel==4 || layerlabel==7 || layerlabel==8) return 512; //bien traite pour le barrel --> refaire pour les caps
	else return 0;
}

bool Cluster::GetSat254() const
{
	return sat254_;
}

bool Cluster::GetSat255() const
{
	return sat255_;
}

bool Cluster::GetShape() const
{
	return shape_;
}

int Cluster::GetNSatStrip(int sat) const
{
	int nsat=0;
	for(int i=0;i<VectStrips_.size();i++)
	{
		if(VectStrips_[i].GetAmpl()==sat) nsat++;
	}
	return nsat;
}

int Cluster::GetNSatStripBoth() const
{
	int nsat=0;
	for(int i=0;i<VectStrips_.size();i++)
	{
		if(VectStrips_[i].GetAmpl()==254 || VectStrips_[i].GetAmpl()==255) nsat++;
	}
	return nsat;
}

int Cluster::GetNStrip() const
{
	return nstrips_;
}	

int Cluster::GetNSimHits() const
{
	return nsimhits_;
}

const vector<ClusterStrip>& Cluster::GetVectStrips() const
{
	return VectStrips_;
}

const vector<SimHit>& Cluster::GetVectSimHits() const
{
	return VectSimHits_;
}

TProfile& Cluster::GetDistribStrip() const
{
	TProfile* profDistribStrip = new TProfile("DistripStrip","DistribStrip",16,0,15);
	profDistribStrip->SetBins(VectStrips_.size()+2,0,VectStrips_.size()+2);
	for(int i=0;i<VectStrips_.size();i++)
	{
		profDistribStrip->Fill(i+1,VectStrips_[i].GetAmpl());
	}
	if(sat254_ && !sat255_) profDistribStrip->SetLineColor(2);
	if(sat255_) profDistribStrip->SetLineColor(3);
	profDistribStrip->GetXaxis()->SetTitle("Strip");
	profDistribStrip->GetYaxis()->SetTitle("Charge");
	return *profDistribStrip;
}
	
bool Cluster::Edge() const
{
	bool Edge_res=false;
	int nstriplayer=this->GetNStripLayer(this->GetLayerLabel());
	if(nstriplayer==0 && (firstsclus_==0 || firstsclus_+nstrips_==512 || firstsclus_+nstrips_==768)) Edge_res=true;
	else if(nstriplayer!=0 && (firstsclus_==0 || firstsclus_+nstrips_==nstriplayer)) Edge_res=true;
	return Edge_res;
}

bool Cluster::Cut() const
{
	bool Cut_res = false;
	if(VectStrips_.size()>=3)
	{
		for(int i=1;i<VectStrips_.size()-1;i++)
		{
			if(VectStrips_[i].GetAmpl()==0) Cut_res=true;
		}
	}
	return Cut_res;
}

int Cluster::GetPartId() const
{
	return partid_;
}

void Cluster::SetPartId(int partid)
{
	partid_ = partid;
}

int Cluster::GetDetId() const
{
	return detid_;
}

int Cluster::GetSubDetId() const
{
	return subdetid_;
}