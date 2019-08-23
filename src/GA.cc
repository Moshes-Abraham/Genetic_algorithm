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
}

int GA::bin2dec(string s)
{
	int value = 0;	
	for (int i = s.size() - 1, j = 0; i >= 0; --i, ++j)	
	{
		if (s[j] == '1')
			value += static_cast<int>(pow(2,i));
	}
	return value;
}

void GA::adapt()
{
	string temp1; 
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
				}
				else
				{
					tmp.push_back(temp1.substr(breakPointPos[pos]));
				}
				segment += 1;
			}
			segment = 0;
			pos = 0;
			input.push_back(tmp);
			tmp.clear();
		}
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
		}
		else
		{
			tmp.push_back(temp1.substr(breakPointPos[pos]));
		}
		segment += 1;
	}
	b_x = bin_x(tmp);
	
	for (int i = 1, j = 0; i <= dimension; ++i, ++j)
	{
		max(i) = b_x[j];
	}
	return max;
}

vector<double> GA::bin_x(vector<string> s)
{
	vector<double> ret;
	for (int i = 0; i != dimension; ++i)
	{
		if (i == (dimension - 1))
			ret.push_back(range(i,0) + bin2dec(s[i]) * ((range(i,1) - range(i,0))/(pow(2,s[i].size()) - 1)));
		else
			ret.push_back(range(i,0) + bin2dec(s[i]) * ((range(i,1) - range(i,0))/(pow(2,breakPointStep) - 1)));
	}
	return ret;
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
	for (int i = 0; i != chromosomeNUM; ++i)
	{
		for (int j = 0; j <= i; ++j)
		{
			qk(i) += pk(j); 
		}
	}
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

	temp = v;
	for (int i = 0; i != chromosomeNUM; ++i)
	{
		v[i] = temp[r(i)];
	}
}

void GA::opcrossover()
{
	vector<string> g;
	std::random_device rd;		// generate real random numbers
	Eigen::ArrayXd r{chromosomeNUM};
	Eigen::ArrayXi mk{0};
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

	if (l % 2 == 1)
	{
		mk.conservativeResize(l - 1);
		l -= 1;
	}

	Eigen::ArrayXi r1{l/2};
	for (int i = 0; i != (l / 2); ++ i)
	{
			r1(i) = ((rd() % 100) / (double)100 ) * (chromosomeNUM - 1);
	}

	for (int i = 0; i != (l / 2); ++i)
	{
		g = onecross(v[mk(2 * i)], v[mk(2 * i + 1)], r1(i));
		v[mk(2 * i)] = g[0];
		v[mk(2 * i + 1)] = g[1];
	}
}

vector<string> GA::onecross(string gene1, string gene2, int pos)
{
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
	for (int i = 0; i != chromosomeNUM; ++i)
	{
		for (int j = 0; j != chromosomeLEN; ++j)
		{
			r(j) = (rd() % 100) / (double)100;
		}

		for (int j = 0, m = 0; j != chromosomeLEN; ++j)
		{
			if (r(j) < pm)
			{
				k.conservativeResize(m + 1);
				k(m) = j;
				++m;
			}
		}
		
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
	cout << "         F(X1,X2): " << maxrec(0) << endl;
	cout << "             FROM: " << mark << "(TH) GENERATION " << endl;
	cout << "             TIME: " << t.format(6,"%w SECONDS") << endl;
	cout << endl << "********************************************************" << endl;

}

void GA::setRange(Eigen::ArrayX2d r)
{
	range = r;
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
}
