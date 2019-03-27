#include "TChain.h"

#include "Track.h"

const int ngenpart_max = 226; //taille obtenue grace au MakeClass de root, applique sur le fichier merge.root (fusion des trois ntuples)
const int ntracks_max = 216;//205;
const int ndedxhits_max = 3417;//3157;
const int nstrips_max = 8138;//8135;
const int nsimhits_max = 2847;//2712;
const int nentries = 15188;

class Builder{

	private:

		TChain* chain_;
		vector<Track> VectTrack_;

		Int_t           runNumber;
   		Int_t           event;
   		Int_t           ngenpart;
   		Int_t           gen_pdg[ngenpart_max];   //[ngenpart]
   		Float_t         gen_pt[ngenpart_max];   //[ngenpart]
   		Float_t         gen_eta[ngenpart_max];   //[ngenpart]
   		Float_t         gen_phi[ngenpart_max];   //[ngenpart]
   		Bool_t          gen_isHardProcess[ngenpart_max];   //[ngenpart]
   		Int_t           gen_status[ngenpart_max];   //[ngenpart]
   		Int_t           gen_moth_pdg[ngenpart_max];   //[ngenpart]
   		Int_t           gen_ndaughter[ngenpart_max];   //[ngenpart]
   		Int_t           gen_daughter_pdg[ngenpart_max];   //[ngenpart]
   		Int_t           ntracks;
   		Float_t         track_pt[ntracks_max];   //[ntracks]
   		Float_t         track_p[ntracks_max];   //[ntracks]
   		Float_t         track_eta[ntracks_max];   //[ntracks]
   		Float_t         track_phi[ntracks_max];   //[ntracks]
   		Float_t         track_chi2[ntracks_max];   //[ntracks]
   		Int_t           track_nvalidhits[ntracks_max];   //[ntracks]
   		Int_t           track_index_hit[ntracks_max];   //[ntracks]
   		Int_t           track_nhits[ntracks_max];   //[ntracks]
   		Int_t           ndedxhits;
   		Int_t           dedx_detid[ndedxhits_max];   //[ndedxhits]
   		Int_t           dedx_subdetid[ndedxhits_max];   //[ndedxhits]
   		Float_t         dedx_charge[ndedxhits_max];   //[ndedxhits]
   		Float_t         dedx_pathlength[ndedxhits_max];   //[ndedxhits]
   		Float_t         dedx_posx[ndedxhits_max];   //[ndedxhits]
   		Float_t         dedx_posy[ndedxhits_max];   //[ndedxhits]
   		Float_t         dedx_posz[ndedxhits_max];   //[ndedxhits]
   		Bool_t          dedx_isstrip[ndedxhits_max];   //[ndedxhits]
   		Bool_t          dedx_ispixel[ndedxhits_max];   //[ndedxhits]
   		Int_t           sclus_firstsclus[ndedxhits_max];   //[ndedxhits]
   		Float_t         sclus_barycenter[ndedxhits_max];   //[ndedxhits]
   		Float_t         sclus_charge[ndedxhits_max];   //[ndedxhits]
   		Float_t         sclus_errorclus[ndedxhits_max];   //[ndedxhits]
   		Bool_t          sclus_ismerged[ndedxhits_max];   //[ndedxhits]
   		Int_t           sclus_index_strip[ndedxhits_max];   //[ndedxhits]
   		Int_t           sclus_nstrip[ndedxhits_max];   //[ndedxhits]
   		Bool_t          sclus_sat254[ndedxhits_max];   //[ndedxhits]
   		Bool_t          sclus_sat255[ndedxhits_max];   //[ndedxhits]
   		Bool_t          sclus_shape[ndedxhits_max];   //[ndedxhits]
   		Int_t           sclus_index_simhit[ndedxhits_max];   //[ndedxhits]
   		Int_t           sclus_nsimhit[ndedxhits_max];   //[ndedxhits]
   		Float_t         sclus_eloss[ndedxhits_max];   //[ndedxhits]
   		Int_t           nstrips;
   		Int_t           strip_ampl[nstrips_max];   //[nstrips]
   		Int_t           nsimhits;
   		Int_t           simhit_pid[nsimhits_max];   //[nsimhits]
   		Int_t           simhit_process[nsimhits_max];   //[nsimhits]
   		Float_t         simhit_p[nsimhits_max];   //[nsimhits]
   		Float_t         simhit_eloss[nsimhits_max];   //[nsimhits]
   		Float_t         simhit_tof[nsimhits_max];   //[nsimhits]
   		Float_t         simhit_segment[nsimhits_max];   //[nsimhits]
   		Float_t         simhit_xentry[nsimhits_max];   //[nsimhits]
   		Float_t         simhit_yentry[nsimhits_max];   //[nsimhits]
   		Float_t         simhit_zentry[nsimhits_max];   //[nsimhits]
   		Float_t         simhit_xexit[nsimhits_max];   //[nsimhits]
   		Float_t         simhit_yexit[nsimhits_max];   //[nsimhits]
   		Float_t         simhit_zexit[nsimhits_max];   //[nsimhits]

	public:

		Builder();///default constructor
		Builder(TChain &chain);///constructor using the input files 
		~Builder();///default destructor 
		void SetBranchAdd();///set the branch adresses
		void GetEntry(int i);///load entry i
		int GetEntries() const;///return the number of entries of the files 
		int GetNtracks() const;///return the number of tracks of the current event
		const vector<Track>& GetVectTrack() const;///return the tracks of the current event 
		void SetCalibration(float factor,int entries);///loop on the tree and recalibration 

};
