#include "../interface/ClusterStrip.h"

ClusterStrip::ClusterStrip()
{
	ampl_ = .0;
}

ClusterStrip::ClusterStrip(int ampl)
{
	ampl_ = ampl;
}

ClusterStrip::~ClusterStrip()
{
}

int ClusterStrip::GetAmpl() const
{
	return ampl_;
}
