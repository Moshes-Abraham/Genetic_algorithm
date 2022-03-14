#include "../src/GA.hh"
double MyobjFunc(vector<double>);
double MyobjFunc2(vector<double>);
int main()
{
	int dim = 7;			// dimension of the obj function
	GA g{20000,MyobjFunc2,dim};	// init with times of iteration, obj function and dimension of the obj function
	Eigen::ArrayX2d r{dim,2};
	r <<	-2.0, 2.0,		// range of x1
		4.1, 5.8,		// range of x2
		4.1, 5.8,		// range of x3
		4.1, 5.8,		// range of x4
		4.1, 5.8,		// range of x5
		4.1, 5.8,		// range of x6
		-2.0, 2.0;		// range of x7

	g.setRange(r);
	g.setChromosomeLength(43);
	g.setChromosomeNumber(10);	 
	g.setCrossOverPossibility(0.15);
	g.setVariationPossibility(0.21);
	g.setToFindMinimum(true);
	g.Solve();
}


double MyobjFunc(vector<double> x)
{
	return 21.5 + x[0] * sin (4 * PI * x[0]) + x[1] * sin(20 * PI * x[1]);
}

double MyobjFunc2(vector<double> x) 	
{
	//return 21.5 + x[0] * sin (4 * PI * x[0]) + x[1] * sin(20 * PI * x[1]) - x[2] * sin(6 * PI * x[2]);
	//return -(x[0] * exp(-x[0]*x[0]-x[1]*x[1]));
	//return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]-x[4]*x[4]+x[5]*x[5]-x[6]*x[6]);
	return 21.5 + x[0] * sin (4 * PI * x[0]) + x[1] * sin(20 * PI * x[1])+ x[2] * sin (4 * PI * x[2])+ x[3] * sin (4 * PI * x[3])+ x[4] * sin (4 * PI * x[4])+ x[5] * sin (4 * PI * x[5])+ x[6] * sin (4 * PI * x[6]);
}
