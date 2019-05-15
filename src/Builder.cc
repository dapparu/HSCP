#include <iostream>
#include <vector>

#include "../interface/Builder.h"


using namespace std;

Builder::Builder()
{
    cerr<<"Input a TTree"<<endl;
}

Builder::Builder(TChain &chain)
{
    chain_  =   &chain;

    ntracks=0;

	ThresholdPartId_ = 0.;
	ThresholdPt_	 = 0.;

	for(int i=0;i<ngenpart_max;i++){
		gen_pdg[i]=0;
		gen_pt[i]=0;
		gen_eta[i]=0;
		gen_phi[i]=0;
		gen_isHardProcess[i]=0;
		gen_status[i]=0;
		gen_moth_pdg[i]=0;
		gen_ndaughter[i]=0;
		gen_daughter_pdg[i]=0;
	}
	for(int i=0;i<ntracks_max;i++){
		track_pt[i]=0;
		track_p[i]=0;
		track_nhits[i]=0;
		track_index_hit[i]=0;
		track_eta[i]=0;
	}
	for(int i=0;i<ndedxhits_max;i++){
		dedx_detid[i]=0;
		dedx_subdetid[i]=0;
		dedx_charge[i]=0;
		dedx_pathlength[i]=0;
		dedx_isstrip[i]=0;
		sclus_sat254[i]=0;
		sclus_sat255[i]=0;
		sclus_index_strip[i]=0;
		sclus_nstrip[i]=0;
		sclus_index_simhit[i]=0;
		sclus_nsimhit[i]=0;
		sclus_eloss[i]=0;
		sclus_firstsclus[i]=0;
	}
	for(int i=0;i<nstrips_max;i++){
		strip_ampl[i]=0;
	}
	for(int i=0;i<nsimhits_max;i++){
		simhit_pid[i]=0;
		simhit_p[i]=0;
		simhit_eloss[i]=0;
		simhit_tof[i]=0;
		simhit_xentry[i]=0;
		simhit_xexit[i]=0;
	}
}

Builder::~Builder()
{
	delete chain_;
	VectTrack_.clear();
}

void Builder::SetBranchAdd()
{
    chain_->SetBranchAddress("runNumber", &runNumber);
	chain_->SetBranchAddress("event", &event);
	chain_->SetBranchAddress("ngenpart", &ngenpart);
	chain_->SetBranchAddress("gen_pdg", &gen_pdg);
	chain_->SetBranchAddress("gen_pt", &gen_pt);
	chain_->SetBranchAddress("gen_eta", &gen_eta);
	chain_->SetBranchAddress("gen_phi", &gen_phi);
	chain_->SetBranchAddress("gen_isHardProcess", &gen_isHardProcess);
	chain_->SetBranchAddress("gen_status", &gen_status);
	chain_->SetBranchAddress("gen_moth_pdg", &gen_moth_pdg);
	chain_->SetBranchAddress("gen_ndaughter", &gen_ndaughter);
	chain_->SetBranchAddress("gen_daughter_pdg", &gen_daughter_pdg);
	chain_->SetBranchAddress("ntracks", &ntracks);
	chain_->SetBranchAddress("track_pt", &track_pt);
	chain_->SetBranchAddress("track_p", &track_p);
	chain_->SetBranchAddress("track_eta", &track_eta);
	chain_->SetBranchAddress("track_phi", &track_phi);
	chain_->SetBranchAddress("track_chi2", &track_chi2);
	chain_->SetBranchAddress("track_nvalidhits", &track_nvalidhits);
	chain_->SetBranchAddress("track_index_hit", &track_index_hit);
	chain_->SetBranchAddress("track_nhits", &track_nhits);
	chain_->SetBranchAddress("ndedxhits", &ndedxhits);
	chain_->SetBranchAddress("dedx_detid", &dedx_detid);
	chain_->SetBranchAddress("dedx_subdetid", &dedx_subdetid);
	chain_->SetBranchAddress("dedx_charge", &dedx_charge);
	chain_->SetBranchAddress("dedx_pathlength", &dedx_pathlength);
	chain_->SetBranchAddress("dedx_posx", &dedx_posx);
	chain_->SetBranchAddress("dedx_posy", &dedx_posy);
	chain_->SetBranchAddress("dedx_posz", &dedx_posz);
	chain_->SetBranchAddress("dedx_isstrip", &dedx_isstrip);
	chain_->SetBranchAddress("dedx_ispixel", &dedx_ispixel);
	chain_->SetBranchAddress("sclus_firstsclus", &sclus_firstsclus);
	chain_->SetBranchAddress("sclus_barycenter", &sclus_barycenter);
	chain_->SetBranchAddress("sclus_charge", &sclus_charge);
	chain_->SetBranchAddress("sclus_errorclus", &sclus_errorclus);
	chain_->SetBranchAddress("sclus_ismerged", &sclus_ismerged);
	chain_->SetBranchAddress("sclus_index_strip", &sclus_index_strip);
	chain_->SetBranchAddress("sclus_nstrip", &sclus_nstrip);
	chain_->SetBranchAddress("sclus_sat254", &sclus_sat254);
	chain_->SetBranchAddress("sclus_sat255", &sclus_sat255);
	chain_->SetBranchAddress("sclus_shape", &sclus_shape);
	chain_->SetBranchAddress("sclus_index_simhit", &sclus_index_simhit);
	chain_->SetBranchAddress("sclus_nsimhit", &sclus_nsimhit);
	chain_->SetBranchAddress("sclus_eloss", &sclus_eloss);
	chain_->SetBranchAddress("nstrips", &nstrips);
	chain_->SetBranchAddress("strip_ampl", &strip_ampl);
	chain_->SetBranchAddress("nsimhits", &nsimhits);
	chain_->SetBranchAddress("simhit_pid", &simhit_pid);
	chain_->SetBranchAddress("simhit_process", &simhit_process);
	chain_->SetBranchAddress("simhit_p", &simhit_p);
	chain_->SetBranchAddress("simhit_eloss", &simhit_eloss);
	chain_->SetBranchAddress("simhit_tof", &simhit_tof);
	chain_->SetBranchAddress("simhit_segment", &simhit_segment);
	chain_->SetBranchAddress("simhit_xentry", &simhit_xentry);
	chain_->SetBranchAddress("simhit_yentry", &simhit_yentry);
	chain_->SetBranchAddress("simhit_zentry", &simhit_zentry);
	chain_->SetBranchAddress("simhit_xexit", &simhit_xexit);
	chain_->SetBranchAddress("simhit_yexit", &simhit_yexit);
	chain_->SetBranchAddress("simhit_zexit", &simhit_zexit);
}

void Builder::GetEntry(int i)
{
	VectTrack_.clear();
    chain_->GetEntry(i);
	vector<int> 	VectPartID_Cluster;
	int SizeSimHit=0;
    if(ntracks>0)
    {
        for(int itrack=0;itrack<ntracks;itrack++)
        {
			//cout<<"track "<<itrack<<endl;
            if(track_pt[itrack]>ThresholdPt_ && track_p[itrack]>0 && track_chi2[itrack]<=5 && track_nvalidhits[itrack]>=8 && track_nhits[itrack]<=40 && -0.8<=track_eta[itrack] && track_eta[itrack]<=0.8) //on applique des selections sur nos traces, des criteres de qualite. On veut minimum pt>threshold GeV 
            {	
				float ThresholdPartId=1;
                vector<Cluster> VectClust;
                for(int iclust=track_index_hit[itrack];iclust<track_index_hit[itrack]+track_nhits[itrack];iclust++)
                {
					//cout<<"clust "<<iclust<<endl;
                    if(dedx_isstrip[iclust]==true && dedx_charge[iclust]>0 && dedx_pathlength[iclust]>0)
                    {
                        vector<ClusterStrip> VectStrips;
                        vector<SimHit> VectSimHits;
                        for(int istrip=sclus_index_strip[iclust];istrip<sclus_index_strip[iclust]+sclus_nstrip[iclust];istrip++)
                        {
							//cout<<"strip "<<istrip<<endl;
                            ClusterStrip strip1(strip_ampl[istrip]);
                            VectStrips.push_back(strip1);
                        }
                        for(int isimhit=sclus_index_simhit[iclust];isimhit<sclus_index_simhit[iclust]+sclus_nsimhit[iclust];isimhit++)
                        {
							//cout<<"simhit "<<isimhit<<endl;
                            SimHit simhit1(simhit_pid[isimhit],simhit_p[isimhit],simhit_eloss[isimhit],simhit_tof[isimhit],simhit_xentry[isimhit],simhit_xexit[isimhit]);
                            VectSimHits.push_back(simhit1);
							VectPartID_Cluster.push_back(simhit_pid[isimhit]);
                        }
                        Cluster clust1(dedx_charge[iclust],sclus_charge[iclust]*(3.61*pow(10,-9)*247),dedx_pathlength[iclust],sclus_eloss[iclust],sclus_nstrip[iclust],sclus_nsimhit[iclust],dedx_detid[iclust],dedx_subdetid[iclust],sclus_sat254[iclust],sclus_sat255[iclust],sclus_shape[iclust],sclus_firstsclus[iclust],VectPartID_Cluster[0],VectStrips,VectSimHits);
						if(clust1.Cut()==false && clust1.Edge()==false) VectClust.push_back(clust1); //filtre sur les clusters
						SizeSimHit=VectSimHits.size();
                    }
                }
				if(VectClust.size()>=3)
				{
					if(SizeSimHit!=VectClust.size())
					{
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
    					ThresholdPartId=max;
						for(int iclust=0;iclust<VectClust.size();iclust++) VectClust[iclust].SetPartId(VectPartID[indice]);
						VectPartID.clear();
					}
					Track track1(track_pt[itrack],track_p[itrack],track_eta[itrack],track_phi[itrack],track_nhits[itrack],ndedxhits,VectClust);
					track1.SetPartId(VectClust[0].GetPartId());
					if(ThresholdPartId>=ThresholdPartId_) VectTrack_.push_back(track1); //on ne garde que les traces qui verifient la condition d'au minimum threshold% de partId identique
					VectPartID_Cluster.clear();
				}	
            }
        }
    }
}

int Builder::GetEntries() const
{
    return chain_->GetEntries();
}

int Builder::GetNtracks() const
{
	return VectTrack_.size();
}

const vector<Track>& Builder::GetVectTrack() const
{
	return VectTrack_;
}

void Builder::SetThresholdPartId(float threshold)
{
	ThresholdPartId_ = threshold;
}

float Builder::GetThresholdPartId() const
{
	return ThresholdPartId_;
}

void Builder::SetThresholdPt(float threshold)
{
	ThresholdPt_ = threshold;
}

float Builder::GetThresholdPt() const
{
	return ThresholdPt_;
}

/*void Builder::SetCalibration(float factor,int entries)
{
	chain_->GetEntry(entries);
    if(ntracks>0)
    {
        for(int itrack=0;itrack<ntracks;itrack++)
        {
			//cout<<"track "<<itrack<<endl;
            if(track_pt[itrack]>0 && track_p[itrack]>0)
            {
                vector<Cluster> VectClust;
                for(int iclust=track_index_hit[itrack];iclust<track_index_hit[itrack]+track_nhits[itrack];iclust++)
                {
					//cout<<"clust "<<iclust<<endl;
                    if(dedx_isstrip[iclust]==true && dedx_charge[iclust]>0 && dedx_pathlength[iclust]>0)
                    {
                        vector<ClusterStrip> VectStrips;
                        vector<SimHit> VectSimHits;
                        for(int istrip=sclus_index_strip[iclust];istrip<sclus_index_strip[iclust]+sclus_nstrip[iclust];istrip++)
                        {
                            ClusterStrip strip1(strip_ampl[istrip]);
                            VectStrips.push_back(strip1);
                        }
                        for(int isimhit=sclus_index_simhit[iclust];isimhit<sclus_index_simhit[iclust]+sclus_nsimhit[iclust];isimhit++)
                        {
                            SimHit simhit1(simhit_pid[isimhit],simhit_p[isimhit],simhit_eloss[isimhit],simhit_tof[isimhit]);
                            VectSimHits.push_back(simhit1);
                        }
                        Cluster clust1(dedx_charge[iclust],sclus_charge[iclust]*factor*(3.61*pow(10,-9)*247),dedx_pathlength[iclust],sclus_eloss[iclust],sclus_nstrip[iclust],sclus_nsimhit[iclust],dedx_detid[iclust],dedx_subdetid[iclust],sclus_sat254[iclust],sclus_sat255[iclust],sclus_firstsclus[iclust],VectStrips,VectSimHits);
                        VectClust.push_back(clust1);
                    }
                }
                Track track1(track_pt[itrack],track_p[itrack],track_nhits[itrack],ndedxhits,VectClust);
                VectTrack_.push_back(track1);
            }
        }
    }
}*/
