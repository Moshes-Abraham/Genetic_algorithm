#include <string>
#include <iostream>
#include <iterator>
#include <cmath>
#include <random>
using namespace std;

#include <Eigen/Eigen>
#include <boost/timer/timer.hpp>
#define LEN 33		// 33 bytes
#define NUM 20		/* 20 chromosomes */
#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211704
template<typename T>
string toString(const T& t)
{
	ostringstream oss;
	oss << t;
	return oss.str();
}

class GA
{
	public:
		GA(int n, double (*obj)(vector<double>), int dim) : gen(n), objFunc(obj), dimension(dim)
	{
//		InitGroup();
	}
		~GA(){}

		

	public:
		void Solve();			// solve the problem
		void check();


	//	GA(string func, int n) :objFunc(func), gen(n)
	//	{
	//		
	//		
	//	}

		void InitGroup();
		string generateChromosome();

		int bin2dec(string);
		void adapt();
		Eigen::ArrayXd maxrecord();
		vector<double> bin_x(vector<string>);
		void chfather();		// choose father chromosom
		void opcrossover();		// one point crossover
		vector<string> onecross(string, string, int);
		string vari(string, int);
		void variation(); 		// variate operator

		void setRange(Eigen::ArrayX2d);

		// Advanced options
	//	void configChromosome(int, int, int);
		void setChromosomeLength(int len){chromosomeLEN = len;}
		void setChromosomeNumber(int num){chromosomeNUM = num;}
		//void setChromosomeBreakPoint(){breakPointPos = (chromosomeLEN / dimension)}
		void setCrossOverPossibility(double op){pc = op;}
		void setVariationPossibility(double vp){pm = vp;}
		void setBreakPoint();

		
	//	double objFunc(double, double);
	private:
		int gen; 			// repeat times
		int dimension;		// 2 by default
		//vector<double> input;
		double (*objFunc)(vector<double>); 		// container for the objective function.
		vector<string> v;		// group
		Eigen::ArrayXd record{chromosomeLEN};	// must init with "{}" not with "()", or it will recognized record as a function! 
		Eigen::ArrayXd maxrec;

		//double x1begin = -3.0, x1end = 3.0, x2begin = -3.0, x2end = 3.0;
		Eigen::ArrayX2d range{dimension,2};
		int chromosomeLEN = 33;		// by default
		int chromosomeNUM = 20;

		int breakPointStep; 	//= std::round(chromosomeLEN / dimension);
		vector<int> breakPointPos;
//		vector<int> breakPointPosPrev;
		//int breakPointNUM = dimension - 1;
		
		double pm = 0.01;
		double pc = 0.25;
};
