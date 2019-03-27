#include <vector>

using namespace std;

class Estimator
{

    private :
 
        vector<float> Vect_;

    public :
 
        Estimator();
        Estimator(vector<float>);
        Estimator(vector<int>);
        ~Estimator();
        void SetVect(vector<float>);
        vector<float> GetVect();
        void DispVect();
        float GetMean();
        float GetMean(vector<float>);
        float GetMedian();
        float GetTrunc40();
        float GetHarmonic2();
        float GetStdDev();
};
