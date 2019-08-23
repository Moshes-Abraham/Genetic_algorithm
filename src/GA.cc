#include "GA.hh"

string GA::generateChromosome()
{
	Eigen::ArrayXd temp{chromosomeLEN};
	string s;
	std::random_device rd;		// generate real random numbers
	for (int i = 0; i != chromosomeLEN; ++i)	
	{
		temp(i) = std::round((rd() % 100) / (double)100);
	}

//	temp = temp.Random(LEN);	// fake random numbers
//	temp = temp.abs();
//	temp = round(temp);

	for (short i = 0; i != chromosomeLEN; ++i)
	{
		s += toString(temp(i));
	}

	return s;
}

void GA::InitGroup()
{
	string a ;
	for (int i = 0; i != chromosomeNUM; ++i)
	{
		a = generateChromosome();
		v.push_back(a);
	}

	// test
//	for (auto id : v)
//	{
//		cout << id << endl;
//	}

}

int GA::bin2dec(string s)
{
	// TODO: check string
	int value = 0;	

	for (int i = s.size() - 1, j = 0; i >= 0; --i, ++j)	
	{
		if (s[j] == '1')
			value += static_cast<int>(pow(2,i));
	}
//	cout << value << endl;
	return value;
}

//Eigen::ArrayXXd GA::adapt(double (* objFunc)(double, double))
/*void GA::adapt(double (* objFunc)(double, double))	// my original method, runs well but not quite efficient
{
	string temp1; 
	string tmp;
	vector<string> temp2, temp3;

	for (int i = 0; i != NUM; ++i)
	{
		for (int j = 0; j != NUM; ++j)
		{
			temp1 = v[j];
			for (auto it = temp1.cbegin(); it != (temp1.cbegin() + 18); ++it)
			{
				tmp += *it;
			}

			temp2.push_back(tmp);		// stands for x1
			tmp.clear();
			for (auto it = (temp1.cbegin() + 18); it != temp1.cend(); ++it)
			{
				tmp += *it;
			}

			temp3.push_back(tmp);		// stands for x2
			tmp.clear();
			
		}
		// TODO: Use function pointer
		record(i, 0) = (* objFunc)(bin_x(temp2[i],'1'),bin_x(temp3[i],'2'));
	}
}*/
// This is a more efficient way
void GA::adapt()
{
	//record.conservativeResize(chromosomeNUM);
	string temp1; 
//	vector<string> temp2, temp3;
	vector<vector<string> > input;
	vector<string> tmp;
	int segment = 1;
	int pos = 0;

	for (int i = 0; i != chromosomeNUM; ++i)
	{
		for (int j = 0; j != chromosomeNUM; ++j)
		{
			temp1 = v[j];
			for (int k = 0; k != dimension; ++k)
			{
				if (segment != dimension)
				{
					tmp.push_back(temp1.substr(breakPointPos[pos],breakPointStep));
					pos += 1;
					//input.push_back(tmp);
					//tmp.clear();
				}
				else
				{
					tmp.push_back(temp1.substr(breakPointPos[pos]));
					//input.push_back(tmp);
					//tmp.clear();
				}
				segment += 1;
			}
			segment = 0;
			pos = 0;
			input.push_back(tmp);
			tmp.clear();
			//temp2.push_back(temp1.substr(0,breakPointPos - 1));	// stands for x1
			//temp3.push_back(temp1.substr(breakPointPos));	// stands for x2
		}


	//	record(i) = (* objFunc)(bin_x(temp2[i],'1'),bin_x(temp3[i],'2'));
	//	double d = (* objFunc)(bin_x(input[i]));
	//	cout << "d = " << d << endl;
	//	cout << "record: " << record(i) << endl;
	//	record(i) = d;
		record(i) = (* objFunc)(bin_x(input[i]));
	}
	
}

Eigen::ArrayXd GA::maxrecord()
{
	Eigen::ArrayXd max{dimension + 1};
	Eigen::ArrayXd::Index i;	// Index of m
	double m = record.maxCoeff(&i);
	max(0) = m;
	string temp1 = v[i];
	vector<string> tmp;
	vector<double> b_x;
	int segment = 1;
	int pos = 0;

	for (int k = 0; k != dimension; ++k)
	{
		if (segment != dimension)
		{
			tmp.push_back(temp1.substr(breakPointPos[pos],breakPointStep));
			pos += 1;
			//input.push_back(tmp);
			//tmp.clear();
		}
		else
		{
			tmp.push_back(temp1.substr(breakPointPos[pos]));
			//input.push_back(tmp);
			//tmp.clear();
		}
		segment += 1;
	}
	b_x = bin_x(tmp);
	
	for (int i = 1, j = 0; i <= dimension; ++i, ++j)
	{
		max(i) = b_x[j];
	}

//	string temp2 = temp1.substr(0,breakPointPos - 1);
//	string temp3 = temp1.substr(breakPointPos);
////	cout << endl << temp2 << "\t" << temp3 << endl << endl;	// test
//	max(1) = bin_x(temp2,'1');
//	max(2) = bin_x(temp3,'2');
	return max;
}

vector<double> GA::bin_x(vector<string> s)
{
	//double num = 0.0;
	vector<double> ret;
	//vector<string> s;
	for (int i = 0; i != dimension; ++i)
	{
		if (i == (dimension - 1))
			ret.push_back(range(i,0) + bin2dec(s[i]) * ((range(i,1) - range(i,0))/(pow(2,s[i].size()) - 1)));
		else
			ret.push_back(range(i,0) + bin2dec(s[i]) * ((range(i,1) - range(i,0))/(pow(2,breakPointStep) - 1)));
	}

	return ret;

//	double num = 0.0;
//
//	switch (opt)
//	{
//		case '1':
//			//return num = -3.0 + bin2dec(bin) * ((12.1 - (-3.0))/(pow(2,18) - 1));
//			return num = x1begin + bin2dec(bin) * ((x1end - (x1begin))/(pow(2,breakPointPos) - 1));
//			break;
//		case '2':
//			//return num = 4.1 + bin2dec(bin) * ((5.8 - 4.1)/(pow(2,15) - 1));
//			return num = x2begin + bin2dec(bin) * ((x2end - (x2begin))/(pow(2,chromosomeLEN - breakPointPos) - 1));
//			break;
//		default:
//			break;
//	}
//	return num;
	
}

void GA::chfather()
{ 
	double F = 0.0;
	int k;
	vector<string> temp;
	Eigen::ArrayXd pk{chromosomeNUM};
	Eigen::ArrayXd qk{chromosomeNUM};
	Eigen::ArrayXd r{chromosomeNUM};
	for (int i = 0; i != chromosomeNUM; ++i)
	{
		F += record(i);		// the sun of all values in record
	}
	pk = record / F;
	
	qk = qk.LinSpaced(chromosomeNUM,0,0); 
	//qk(0) = pk(0);
	for (int i = 0; i != chromosomeNUM; ++i)
	{
		for (int j = 0; j <= i; ++j)
		{
			qk(i) += pk(j); 
		}
	}
	// test
//	cout << endl << qk << endl;

	std::random_device rd;		// generate real random numbers
	for (int i = 0; i != chromosomeNUM; ++i)	
	{
		r(i) = (rd() % 100) / (double)100;
	}
	for (int i = 0; i != chromosomeNUM; ++i)
	{
		k = 1;
		while (r(i) > qk(k))
		{
			k += 1;
		}
		r(i) = k;
	}
	// test
//	cout << endl << r << endl;

	temp = v;
	for (int i = 0; i != chromosomeNUM; ++i)
	{
		v[i] = temp[r(i)];
	}
	// test
//	for (auto &it : v)
//		cout << it << endl;


}

void GA::opcrossover()
{
	vector<string> g;
	std::random_device rd;		// generate real random numbers
	Eigen::ArrayXd r{chromosomeNUM};
	Eigen::ArrayXi mk{0};
	//const float pc = 0.25;
	int l = 1;
	while (l == 1)
	{
		for (int i = 0; i != chromosomeNUM; ++i)	
		{
			r(i) = (rd() % 100) / (double)100;
		}
		for (int i = 0, j = 0; i != chromosomeNUM; ++i)	// ATTENTION! chromosomeNUM should not larger than chromosomeLEN, or mk may larger than chromosomeLEN and cause "std::out_of_ranger"
		{
			if (r(i) < pc)
			{
				mk.conservativeResize(j + 1);
				mk(j) = i;
				++j;
			}
		}
		l = mk.rows();
	}
	// test
//	cout << l << endl << mk << endl;

	if (l % 2 == 1)
	{
		mk.conservativeResize(l - 1);
		l -= 1;
	}
	// test
//	cout << endl << mk << endl;

	Eigen::ArrayXi r1{l/2};
	//int rt = 0;		// temp random number
	for (int i = 0; i != (l / 2); ++ i)
	{
		// TODO: change to 0 ~ 31
	//	rt = ((rd() % 100) / (double)100 ) * (NUM - 1);
	//	if (rt != 0) 
	//	{
	//		r1(r) = rt; 	// bug! it should be: r1(i) = rt;
	//	}
	//	else
	//		r1(i) = 1;
			r1(i) = ((rd() % 100) / (double)100 ) * (chromosomeNUM - 1);
	}

	// test
//	cout << endl << r1 << endl;

	for (int i = 0; i != (l / 2); ++i)
	{
		g = onecross(v[mk(2 * i)], v[mk(2 * i + 1)], r1(i));
		v[mk(2 * i)] = g[0];
		v[mk(2 * i + 1)] = g[1];
	}
}

vector<string> GA::onecross(string gene1, string gene2, int pos)
{
	//int len = gene1.size();
	vector<string> s;
	string g1 = gene1.substr(0,pos) + gene2.substr(pos);	// swap	
	string g2 = gene2.substr(0,pos) + gene1.substr(pos);
	s.push_back(g1);
	s.push_back(g2);
	return s;
}

string GA::vari(string gold, int pos)
{
	string gnew(gold);
	if (gnew[pos] == '1')
		gnew[pos] = '0';
	else
		gnew[pos] = '1';
	return gnew;
}

void GA::variation()
{
	std::random_device rd;
	Eigen::ArrayXd r{chromosomeLEN};
	Eigen::ArrayXi k{0};
	//const float pm = 0.01;
	for (int i = 0; i != chromosomeNUM; ++i)
	{
		for (int j = 0; j != chromosomeLEN; ++j)
		{
			r(j) = (rd() % 100) / (double)100;
		}

		// test
	//	cout << endl << r << endl;
		for (int j = 0, m = 0; j != chromosomeLEN; ++j)
		{
			if (r(j) < pm)
			{
				k.conservativeResize(m + 1);
				k(m) = j;
				++m;
			}
		}
//		cout << k << endl;
		
		for (int j = 0; j != k.rows(); ++j)
		{
			v[i] = vari(v[i], k(j));
		}
	}
}

void GA::Solve()
{
	check();
	InitGroup();
	boost::timer::cpu_timer t;
	Eigen::ArrayXd temp;
	int mark = 0;

	setBreakPoint();
	adapt();
	maxrec = maxrecord();

	t.start();
	for (int i = 0; i != gen; ++i)
	{
		chfather();		// chose father chromosome

		// genetic operator
		opcrossover();
		variation();

		adapt();
		temp = maxrecord();

		if (temp(0) > maxrec(0))	// keep the best value
		{
			maxrec = temp;
			mark = i;
		}

	}
	t.stop();

	// outputs
	cout << endl << "********************************************************" << endl << endl;
	cout << "  REPORT:          " << endl;
	for (int i = 1; i <= dimension; ++i)
	{
		cout << "               X" << i << ": " << maxrec(i) << endl;
	}
//	cout << "               X1: " << maxrec(1) << endl;
//	cout << "               X2: " << maxrec(2) << endl;
	cout << "         F(X1,X2): " << maxrec(0) << endl;
	cout << "             FROM: " << mark << "(TH) GENERATION " << endl;
	cout << "             TIME: " << t.format(6,"%w SECONDS") << endl;
	cout << endl << "********************************************************" << endl;

}

void GA::setRange(Eigen::ArrayX2d r)
{
	range = r;
//	x1begin	= x1b;
//	x1end	= x1e;
//	x2begin	= x2b;
//	x2end	= x2e;
}

void GA::check()
{
	if (chromosomeNUM >= chromosomeLEN /*|| breakPointPos >= chromosomeLEN*/)
	{
		cerr << "Error: Invalied values for chromosome! Please reset proper chromosome values." << endl; 
		exit(0);
	}
	if (chromosomeLEN <= 0 || chromosomeNUM <= 0 /*|| breakPointPos <= 0*/)
	{
		cerr << "Error: Invalied values for chromosome! Please reset proper chromosome values." << endl; 
		exit(0);
	}
	if (pc <= 0 || pm <= 0 || pc >= 1 || pm >= 1)
	{
		cerr << "Error: Invalied values for the possibility!";
		exit(0);
	}
	if (range.rows() == 0)
	{
		cerr << "Error: You didn't set the range!" << endl;
		exit(0);
	}
	if (dimension <= 0)
	{
		cerr << "Error: Invalied dimention!" << endl;
		exit(0);
	}
	if (chromosomeLEN <= dimension)
	{
		cerr << "Error: Invalied chromosome length. You should set it larger than dimension." << endl;
	}
	for (int i = 0; i != dimension; ++i)
	{
		if (range(i,0) >= range(i,1))
		{
			cerr << "Error: Invalied range! The lower range should be smaller than the higher range." << endl;
			exit(0);
		}
	}

}

void GA::setBreakPoint()
{
	breakPointStep = std::round(chromosomeLEN / dimension);
	int bpsfb = 0;		// break point step for the beginning
	for (int i = 0; i != dimension; ++i)
	{
		breakPointPos.push_back(bpsfb);
		bpsfb += breakPointStep;
	}

//	bpsfb = 0;
//	for (int i = 0; i != dimension - 1; ++i)
//	{
//		breakPointPosPrev.push_back(bpsfb);
//		bpsfb += breakPointStep;
//	}
}
