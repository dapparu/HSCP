#include "../interface/Estimator.h"


#include <math.h>
#include <vector>
#include <iostream> 

using namespace std;

Estimator::Estimator()
{
	Vect_.push_back(1);
}

Estimator::Estimator(vector<float> Vect)
{
	Vect_ = Vect;
}

Estimator::Estimator(vector<int> Vect)
{
	for(int i=0;i<Vect.size();i++)
	{
		Vect_.push_back(Vect[i]);
	}
}

Estimator::~Estimator()
{
	Vect_.clear();
}

void Estimator::SetVect(vector<float> Vect)
{
	Vect_ = Vect;
}

vector<float> Estimator::GetVect()
{
	return Vect_;
}

void Estimator::DispVect()
{
	for(int i=0;i<Vect_.size();i++)
	{
		cout<<"Element "<<i<<" du vecteur : "<<Vect_[i]<<endl;
	}
}

float Estimator::GetMean()
{
	float mean = .0;
	float size = Vect_.size();
	for(int i=0;i<size;i++)
	{
		mean += Vect_[i];
	}
	return mean/size;
}

float Estimator::GetMean(vector<float> Vect)
{
	float mean = .0;
	float size = Vect.size();
	for(int i=0;i<size;i++)
	{
		mean+=Vect[i];
	}
	return mean/size;
}

float Estimator::GetMedian()
{
	sort(Vect_.begin(),Vect_.end());
	float size = Vect_.size();
	float med = Vect_[size/2.];
	return med;
}

float Estimator::GetTrunc40()
{
	sort(Vect_.begin(),Vect_.end());
	float size = Vect_.size();
	float size60 = size*0.6;
	vector<float> Vect60;
	for(int i=0; i<size60; i++)
	{
		Vect60.push_back(Vect_[i]);
	}
	//SetVect(Vect60);
	return this->GetMean(Vect60);
}

float Estimator::GetHarmonic2()
{
	float mean = .0;
	float size = Vect_.size();
	for(int i=0;i<size;i++)
	{
		mean+=pow(Vect_[i],-2);
	}
	return pow(mean/size,-0.5);
}

float Estimator::GetStdDev()
{
	float size = Vect_.size();
	float res=.0;
	for(int i=0; i<size; i++)
	{
		//cout<<"Vect_ = "<<Vect_[i]<<"    harmonic2 = "<<harmonic2()<<endl;
		res+=pow(Vect_[i]-this->GetMean(),2);
	}
	res=pow(res/(float)(size-1),0.5);
	return res;
}

