#include <vector>

using namespace std;

class SimHit{

	private:
		
		int 	partId_;
		float 	p_;
		float 	eloss_;
		float	tof_;
		float 	xentry_;
		float 	xexit_;
	
	public:

		SimHit();
		SimHit(int partId,float p,float eloss,float tof,float xentry,float xexit);
		~SimHit();
		int GetPartId() const;
		float GetP() const;
		float GetEloss() const;
		float GetToF() const;
		float GetXentry() const;
		float GetXexit() const;

};
