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
		float				nhits_;
		int					ndedxhits_;
		vector<Cluster>		VectClusters_;
		int 				PartId_;

	public:

		Track();
		Track(float pt,float p,float nhits,int ndedxhits,const vector<Cluster> &VectClusters);
		~Track();
		float GetPt() const;
		float GetP() const;
		int GetNhits() const;
		int GetNCluster() const;
		int GetNSatCluster(int sat) const;
		int GetNSatClusterBoth() const;
		const vector<Cluster>& GetVectClusters() const;
		TCanvas& GetProfCluster() const;
		int GetPartId() const;
		void SetPartId(int partid);

};
