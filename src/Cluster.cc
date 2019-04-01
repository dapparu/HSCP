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
	sat254_				= false;
	sat255_				= false;
}

Cluster::Cluster(float dedx_charge,float sclus_charge,float pathlength,float eloss,int nstrips,int nsimhits,int detid,int subdetid,bool sat254,bool sat255,const vector<ClusterStrip> &VectStrips,const vector<SimHit> &VectSimHits)
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
	VectStrips_			= VectStrips;
	VectSimHits_ 		= VectSimHits;
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

bool Cluster::GetSat254() const
{
	return sat254_;
}

bool Cluster::GetSat255() const
{
	return sat255_;
}

int Cluster::GetNSatStrip(int sat) const
{
	if(sat==254)
	{
		int nsat254=0;
		for(int i=0;i<VectStrips_.size();i++)
		{
			if(VectStrips_[i].GetAmpl()==254) nsat254++;
		}
		return nsat254;
	}
	if(sat==255)
	{
		int nsat255=0;
		for(int i=0;i<VectStrips_.size();i++)
		{
			if(VectStrips_[i].GetAmpl()==255) nsat255++;
		}
		return nsat255;
	}
	else return 0;
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
	
