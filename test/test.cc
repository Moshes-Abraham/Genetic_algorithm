#include "GA.hh"
#include <cmath>
#include <Eigen/Eigen>
#include <iostream>
#include <random>
using namespace std;

double MyobjFunc(vector<double>);
double MyobjFunc2(vector<double>);
int main()
{
//	Eigen::ArrayXd temp;
//	string v;
//
//	temp = temp.Random(30);
//	temp = temp.abs();
//	temp = round(temp);
//
//	cout << temp << endl;
//
//	for(int i = 0; i != 30; ++i)
//	{
//		cout << toString(temp(i)) << "\t";
//	}
//	cout << endl;
	//	GA g;
	//	g.InitGroup();
//	string s("10010011");
//	int value = 0.0;	
//	for (int i = 7, j = 0; i >= 0; --i, ++j)	
//	{
//		if (s[j] == '1')
//			value += static_cast<int>(pow(2,i));
//	}
//	cout << value << endl;
//	cout << s.size() << endl;
	//		GA g;
	//		string st("10110011100101");
	//		g.bin2dec(st);
//	GA g;
//	g.InitGroup();
//	string temp1; 
//	string tmp;
//	vector<string> temp2, temp3;
//	Eigen::ArrayXXd record(NUM,1);
//
//	for (int i = 0; i != NUM; ++i)
//	{
//		for (int j = 0; j != NUM; ++j)
//		{
//			temp1 = g.v[j];
//			for (auto it = temp1.cbegin(); it != (temp1.cbegin() + 18); ++it)
//			{
//				tmp += *it;
//			}
//
//			temp2.push_back(tmp);		// stands for x1
//			tmp.clear();
//			for (auto it = (temp1.cbegin() + 18); it != temp1.cend(); ++it)
//			{
//				tmp += *it;
//			}
//
//			temp3.push_back(tmp);		// stands for x2
//			tmp.clear();
//			
//		}
//		record(i, 0) = objFunc(g.bin2dec(temp2[i]),g.bin2dec(temp3[i]));
//	}
//
//	for (int i = 0; i != NUM; ++i)
//	{
//		cout << temp2[i] << "\t" << temp3[i];
//		cout << endl;
//	}
//	
//	cout << "record = " << endl;
//	cout << record << endl;
//		GA g;
//		g.InitGroup();		// remember to Init or it'll cause segmentation fault
//		g.adapt(MyobjFunc);
//		cout << endl;
	//	for (auto &it : g.v)
	//		cout << it << endl;

	//	vector<string> a;
	//	a = g.onecross(g.v[0],g.v[1],4);
	//	cout << endl;
	//	for (auto &it : a)
	//		cout << it << endl;

//		cout << g.record << endl;
//		cout << endl;
//	      g.maxrecord();
//	      cout << endl << g.maxrec << endl << endl;
//	      g.variation();
//	      g.adapt(MyobjFunc);
//	      g.maxrecord();
//	      cout << endl << g.record << endl << endl << g.maxrec << endl;
//
//		g.opcrossover();
//		g.adapt(MyobjFunc);
//		cout << endl;
//		cout << g.record << endl;
//		cout << endl;
//	      g.maxrecord();
//	      cout << endl << g.maxrec << endl << endl;
		
			
//	std::random_device rd;		// generate real random numbers
//	Eigen::ArrayXd r{NUM};
//	Eigen::ArrayXi mk{0};
//	const float pc = 0.25;
//	int l = 1;
//	while (l == 1)
//	{
//	      for (int i = 0; i != NUM; ++i)	
//	      {
//	      	r(i) = (rd() % 100) / (double)100;
//	      }
//
//		for (int i = 0, j = 0; i != NUM; ++i)
//		{
//			if (r(i) < pc)
//			{
//				mk.conservativeResize(j + 1);
//				mk(j) = i;
//				++j;
//			}
//		}
//		l = mk.rows();
//	}
//
//		cout << l << endl << mk << endl;
	//Eigen::Array3d v(1,4212,23321), w(3,4,5);
	//std::cout << v.maxCoeff() << endl; 		// maxCoeff() returns the max value in an array
	//Eigen::Array3d::Index maxRow;
	//double max = v.maxCoeff(&maxRow);		// get the max value and its index (remember it starts from 0)
	//std :: cout << "Max: " << max << " at: " << maxRow << std :: endl;
	//std :: cout << v(2) << std :: endl;
		//	std::random_device myrd;
		//	for (int i = 0; i != 10; ++i)
		//		//std::cout << myrd() % 100 << std::endl;	// generate 10 real random numbers between 0 and 100
		//		std::cout << ((myrd() % 100) /(double)100)<< std::endl;	// generate 10 real random numbers between 0 and 1. Remember (double) ! 
	
	// real test
	//int gen = 0;
	int dim = 3;			// dimension of the obj function
	GA g{2000,MyobjFunc2,dim};	// init with times of iteration, obj function and dimension of the obj function
	Eigen::ArrayX2d r{dim,2};
	r <<	-3.0, 12.1,		// range of x1
		4.1, 5.8,		// range of x2
		2.1, 6.7;		// range of x3

	g.setRange(r);
	g.setChromosomeLength(43);
	g.setChromosomeNumber(10);	// Error! 
	g.setCrossOverPossibility(0.61);
	g.setVariationPossibility(0.31);
	g.Solve();
	//	GA g{100,MyobjFunc};
//		g.InitGroup();
//		for (auto &it : g.v)
//			cout << it << endl;
//		g.adapt();
//		Eigen::Array3d temp;
//		g.maxrec = g.maxrecord();
//		g.chfather();
//		g.opcrossover();
//		g.variation();
//		for (auto &it : g.v)
//			cout << it << endl;
	//	cout << endl << g.record << endl;

}


double MyobjFunc(vector<double> x)
{
	return 21.5 + x[0] * sin (4 * PI * x[0]) + x[1] * sin(20 * PI * x[1]);
//	//return sin(x1) + sin(x2);
//	return x1 * exp(-x1*x1 - x2*x2);
//	//return sqrt(x1*x1+x2*x2);
}

double MyobjFunc2(vector<double> x) 		// A genius idea!
{
	//return x[0] * exp(-x[0]*x[0] - x[1]*x[1] - x[1]*x[1]);
	return 21.5 + x[0] * sin (4 * PI * x[0]) + x[1] * sin(20 * PI * x[1]) - x[2] * sin(6 * PI * x[2]);
}
