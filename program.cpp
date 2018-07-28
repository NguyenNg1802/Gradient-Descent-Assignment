#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>
#include <sstream>
using namespace std;

#define NUM_SAMPLES 20000

struct Point2D
{
    double   x, t;
};

struct LinearModel
{
    double  Ermsd;
    double  a, b;
    double  start_a,start_b;
    int  errHistogram[10];
};

struct GD
{
    int  nIterations;
    int  nFolds;
    double  learningRate;
    Point2D startPoint;
};

struct DataIterations
{
	vector<Point2D> data;
	int size;
	int group;
};

int main()
{

	GD GD;
	LinearModel model;
	DataIterations dataIter;
	string str;
	
	Point2D P[NUM_SAMPLES];
    int i = 0;
    int line = 1;
    
	while (getline(cin,str))
	{
        stringstream ss(str);
        string aa;
        if ( line == 4 ) {
            ss >> aa >> GD.nIterations;
        }
        if ( line == 5) {
            ss >> aa >> GD.learningRate;
        }
        if ( line == 6) {
            ss >> aa >> model.start_a;
        }
        if ( line == 7) {
            ss >> aa >> model.start_b;
        }
        if ( line == 8) {
            ss >> aa >> GD.nFolds;
        }
        if ( line >= 12) {
            ss >> P[i].x >> P[i].t;
            dataIter.data.push_back( P[i] );
            i++;
        }
        line++;
	}

    

	int D = int(dataIter.data.size()) / GD.nFolds;
	dataIter.size = int(dataIter.data.size()) / GD.nFolds;
	dataIter.group = GD.nFolds;


	for (int i = 0; i < dataIter.group; i++)
	{
        double  La = 0, Lb = 0, eAVG = 0, sigma = 0, Vmax = 0, Vmin = 0, sumSigma = 0, sumErmsd = 0;
        int TST_start = dataIter.size*(i);
        int TST_end = (i+1)*dataIter.size - 1;
        if ((i == dataIter.group - 1) && (dataIter.data.size() % dataIter.group != 0))
        {
            TST_end = TST_end + (dataIter.data.size() % dataIter.group) ;
        }
		vector<Point2D> TRN, TST;
		for (int k = 0; k < dataIter.data.size(); k++)
		{
			if ((k >= TST_start) && (k <= TST_end)) TST.push_back(dataIter.data[k]);
			else TRN.push_back(dataIter.data[k]);
		}
		model.a = model.start_a;
		model.b = model.start_b;
		for (int m = 0; m < GD.nIterations; m++)
		{
            double  ga = 0.0, gb = 0.0;
            La = 0.0;
            Lb = 0.0;
			for (int j = 0; j < TRN.size(); j++)
			{
				La += (model.a*TRN[j].x + model.b - TRN[j].t)*TRN[j].x;
				Lb += model.a*TRN[j].x + model.b - TRN[j].t;
			}
			ga = La / sqrt(La*La + Lb*Lb);
			gb = Lb / sqrt(La*La + Lb*Lb);
			model.a = model.a - GD.learningRate*ga;
			model.b = model.b - GD.learningRate*gb;
		}
        double  sumei = 0;
        vector<float> ei(TST.size());
		for (int j = 0; j < TST.size(); j++)
		{
			sumErmsd += (model.a*TST[j].x + model.b - TST[j].t) * (model.a*TST[j].x + model.b - TST[j].t);
			ei[j] = model.a*TST[j].x + model.b - TST[j].t;
			sumei += ei[j];
		}
        eAVG = sumei / TST.size();
        for (int j = 0; j < TST.size(); j++)
        {
            sumSigma += pow(ei[j] - eAVG, 2);
        }
        
		sigma = sqrt(sumSigma / D);
		Vmax = eAVG + 3 * sigma;
		Vmin = eAVG - 3 * sigma;
		model.Ermsd = sqrt(sumErmsd / D);
        
        double  L[10], x = (Vmax - Vmin) / 10;
        L[0] = Vmin;
		for (int i = 1; i < 10; i++)
			L[i] = L[i - 1] + x;
		int s = 0;
        
        
		for (unsigned int i = 0; i < ei.size(); i++)
		{
			if ((ei[i] >= Vmin) && (ei[i] < L[1]))
			{
				model.errHistogram[0]++;
				s++;
			}
			if ((ei[i] >= L[1]) && (ei[i] < L[2]))
			{
				model.errHistogram[1]++;
				s++;
			}
			if ((ei[i] >= L[2]) && (ei[i] < L[3]))
			{
				model.errHistogram[2]++;
				s++;
			}
			if ((ei[i] >= L[3]) && (ei[i] < L[4]))
			{
				model.errHistogram[3]++;
				s++;
			}
			if ((ei[i] >= L[4]) && (ei[i] < L[5]))
			{
				model.errHistogram[4]++;
				s++;
			}
			if ((ei[i] >= L[5]) && (ei[i] < L[6]))
			{
				model.errHistogram[5]++;
				s++;
			}
			if ((ei[i] >= L[6]) && (ei[i] < L[7]))
			{
				model.errHistogram[6]++;
				s++;
			}
			if ((ei[i] >= L[7]) && (ei[i] < L[8]))
			{
				model.errHistogram[7]++;
				s++;
			}
			if ((ei[i] >= L[8]) && (ei[i] < L[9]))
			{
				model.errHistogram[8]++;
				s++;
			}
			if ((ei[i] >= L[9]) && (ei[i] <= Vmax))
			{
				model.errHistogram[9]++;
				s++;
			}
		}

		cout << fixed << setw(10) << setprecision(5) << right << model.a << fixed << setw(10) << setprecision(5) << right << model.b << fixed << setw(10) << setprecision(5) << right << model.Ermsd;
		for (int i = 0; i < 10; i++)
			cout << fixed << setw(10) << setprecision(5) << right << double(model.errHistogram[i]) / s;
		cout << endl;
        
        for (int i = 0; i < 10; i++)
        {
            model.errHistogram[i] = 0;
        }
	}

	return 0;
}
