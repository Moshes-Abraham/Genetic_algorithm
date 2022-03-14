#include <string>
#include <iostream>
#include <iterator>
#include <cmath>
#include <random>
#include <Eigen/Eigen>
#include <boost/timer/timer.hpp>
#define LEN 33					// 33 bytes
#define NUM 20					/* 20 chromosomes */
#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211704
using namespace std;
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
		GA(int n, double (*obj)(vector<double>), int dim) : gen(n), objFunc(obj), dimension(dim){}
		~GA(){}
	public:
		void Solve();			// solve the problem
		void check();
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
		void setChromosomeLength(int len){chromosomeLEN = len;}
		void setChromosomeNumber(int num){chromosomeNUM = num;}
		void setCrossOverPossibility(double op){pc = op;}
		void setVariationPossibility(double vp){pm = vp;}
		void setBreakPoint();
		void setToFindMinimum(bool c){choice = c;}

	private:
		int gen; 				// repeat times
		int dimension;				// 2 by default
		double (*objFunc)(vector<double>); 	// container for the objective function.
		vector<string> v;			// group
		Eigen::ArrayXd record{0};		// must init with "{}" not with "()", or it will recognized record as a function! 
		Eigen::ArrayXd maxrec;
		Eigen::ArrayX2d range{dimension,2};
		int chromosomeLEN = 33;			// by default
		int chromosomeNUM = 20;
		int breakPointStep;		 	//= std::round(chromosomeLEN / dimension);
		vector<int> breakPointPos;
		double pm = 0.01;
		double pc = 0.25;
		bool choice = false;			// to find the maximum by default
};
