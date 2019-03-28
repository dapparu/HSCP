#include "../interface/SimHit.h"

using namespace std;

SimHit::SimHit()
{
	partId_	= 0;
	p_		= .0;
	eloss_	= .0;
	tof_	= .0;
}

SimHit::SimHit(int partId,float p,float eloss,float tof)
{
	partId_	= partId;
	p_	= p;
	eloss_	= eloss;
	tof_	= tof;
}

SimHit::~SimHit()
{
}

int SimHit::GetPartId() const
{
	return partId_;
}

float SimHit::GetP() const
{
	return p_;
}

float SimHit::GetEloss() const
{
	return eloss_;
}

float SimHit::GetToF() const
{
	return tof_;
}
