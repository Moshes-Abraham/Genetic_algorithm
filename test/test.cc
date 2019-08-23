#include "../src/GA.hh"
double MyobjFunc(vector<double>);
double MyobjFunc2(vector<double>);
int main()
{
	int dim = 3;			// dimension of the obj function
	GA g{2000,MyobjFunc2,dim};	// init with times of iteration, obj function and dimension of the obj function
	Eigen::ArrayX2d r{dim,2};
	r <<	-3.0, 12.1,		// range of x1
		4.1, 5.8,		// range of x2
		2.1, 6.7;		// range of x3

	g.setRange(r);
	g.setChromosomeLength(43);
	g.setChromosomeNumber(10);	 
	g.setCrossOverPossibility(0.61);
	g.setVariationPossibility(0.31);
	g.Solve();
}


double MyobjFunc(vector<double> x)
{
	return 21.5 + x[0] * sin (4 * PI * x[0]) + x[1] * sin(20 * PI * x[1]);
}

double MyobjFunc2(vector<double> x) 	
{
	return 21.5 + x[0] * sin (4 * PI * x[0]) + x[1] * sin(20 * PI * x[1]) - x[2] * sin(6 * PI * x[2]);
}
