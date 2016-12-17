#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include <limits>

#include "miehps.h"

// const double pi=3.1415926535897932384626433832795028841971693993751058209749;

using namespace std;

int main(){
	cout.precision(numeric_limits<double>::digits10 + 1);
	cout.precision(10);
	cout << endl;
	cout << "Calculating Qext and Qsca:" << endl;
	long double qext, qsca;
	double x;
	int n;
	complex<double> m;

	complex<double> s010,s020,s10,s20;

	Mie *mie1;

	x=10000;
	m=complex<double>(10.0,-10.0);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS((long double)0.0);
	s10=mie1->getS1();
	s20=mie1->getS2();
	cout << "Case m: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl << endl;
	
	delete mie1;
}

