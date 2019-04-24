#include <vector>
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"

#include "Cluster.h"

using namespace std;

class Track{

	private:
		
		float				pt_;
		float				p_;
		float 				eta_;
		float 				phi_;
		float				nhits_;
		int					ndedxhits_;
		vector<Cluster>		VectClusters_;
		int 				PartId_;

	public:

		Track();
		Track(float pt,float p,float eta,float phi,float nhits,int ndedxhits,const vector<Cluster> &VectClusters);
		~Track();
		float GetPt() const;
		float GetP() const;
		float GetEta() const;
		float GetPhi() const;
		int GetNhits() const;
		int GetNCluster() const;
		int GetNSatCluster(int sat) const;
		int GetNSatClusterBoth() const;
		const vector<Cluster>& GetVectClusters() const;
		TCanvas& GetProfCluster() const;
		int GetPartId() const;
		void SetPartId(int partid);

};
