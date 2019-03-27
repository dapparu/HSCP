#include <vector>

using namespace std;

class SimHit{

	private:
		
		int 	partId_;
		float 	p_;
		float 	eloss_;
		float	tof_;
	
	public:

		SimHit();
		SimHit(int partId,float p,float eloss,float tof);
		~SimHit();
		int GetPartId() const;
		float GetP() const;
		float GetEloss() const;
		float GetToF() const;

};
