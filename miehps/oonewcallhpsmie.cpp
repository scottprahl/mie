#include <iostream>
#include <cmath>
#include <complex>
#include <limits>

#include "miehps.h"

const long double pi=3.1415926535897932384626433832795028841971693993751058209749;

using namespace std;

int main(){
	cout.precision(numeric_limits<double>::digits10 + 1);
	cout.precision(7);
	
	cout << "------ Tests ------" << endl << endl;

	cout << "Calculating Qext:" << endl;
	long double qext,qsca;
	double x;
	int n;
	complex<double> m;
	complex<long double> s10(0,0);
	complex<long double> s20(0,0);
	complex<long double> s1pi(0,0);
	complex<long double> s2pi(0,0);
	
	Mie *mie1;

	x=0.099;
	m=complex<double>(0.75,0);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS(0);
	s10=mie1->getS1();
	s20=s10;
	mie1->calcS(pi);
	s1pi=mie1->getS1();
	s2pi=mie1->getS2();
	cout << "Case a: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl;
	cout << "Amplitudes: S1(pi)=" << s1pi << " S2(pi)=" << s2pi << endl << endl;
	delete mie1;
	
	x=0.101;
	m=complex<double>(0.75,0);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS(0);
	s10=mie1->getS1();
	s20=s10;
	mie1->calcS(pi);
	s1pi=mie1->getS1();
	s2pi=mie1->getS2();
	cout << "Case b: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl;
	cout << "Amplitudes: S1(pi)=" << s1pi << " S2(pi)=" << s2pi << endl << endl;
	delete mie1;

	x=10.0;
	m=complex<double>(0.75,0);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS(0);
	s10=mie1->getS1();
	s20=s10;
	mie1->calcS(pi);
	s1pi=mie1->getS1();
	s2pi=mie1->getS2();
	cout << "Case c: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl;
	cout << "Amplitudes: S1(pi)=" << s1pi << " S2(pi)=" << s2pi << endl << endl;
	delete mie1;

	x=1000.0;
	m=complex<double>(0.75,0);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS(0);
	s10=mie1->getS1();
	s20=s10;
	mie1->calcS(pi);
	s1pi=mie1->getS1();
	s2pi=mie1->getS2();
	cout << "Case d: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl;
	cout << "Amplitudes: S1(pi)=" << s1pi << " S2(pi)=" << s2pi << endl << endl;
	delete mie1;

	x=100.0;
	m=complex<double>(1.33,-1E-5);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS(0);
	s10=mie1->getS1();
	s20=s10;
	mie1->calcS(pi);
	s1pi=mie1->getS1();
	s2pi=mie1->getS2();
	cout << "Case e: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl;
	cout << "Amplitudes: S1(pi)=" << s1pi << " S2(pi)=" << s2pi << endl << endl;
	delete mie1;

	x=10000.0;
	m=complex<double>(1.33,-1E-5);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS(0);
	s10=mie1->getS1();
	s20=s10;
	mie1->calcS(pi);
	s1pi=mie1->getS1();
	s2pi=mie1->getS2();
	cout << "Case f: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl;
	cout << "Amplitudes: S1(pi)=" << s1pi << " S2(pi)=" << s2pi << endl << endl;
	delete mie1;

	x=0.055;
	m=complex<double>(1.5,-1);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS(0);
	s10=mie1->getS1();
	s20=s10;
	mie1->calcS(pi);
	s1pi=mie1->getS1();
	s2pi=mie1->getS2();
	cout << "Case g: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl;
	cout << "Amplitudes: S1(pi)=" << s1pi << " S2(pi)=" << s2pi << endl << endl;
	delete mie1;

	x=0.056;
	m=complex<double>(1.5,-1.0);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS(0);
	s10=mie1->getS1();
	s20=s10;
	mie1->calcS(pi);
	s1pi=mie1->getS1();
	s2pi=mie1->getS2();
	cout << "Case h: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl;
	cout << "Amplitudes: S1(pi)=" << s1pi << " S2(pi)=" << s2pi << endl << endl;
	delete mie1;

	x=100;
	m=complex<double>(1.5,-1.0);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS(0);
	s10=mie1->getS1();
	s20=s10;
	mie1->calcS(pi);
	s1pi=mie1->getS1();
	s2pi=mie1->getS2();
	cout << "Case i: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl;
	cout << "Amplitudes: S1(pi)=" << s1pi << " S2(pi)=" << s2pi << endl << endl;
	delete mie1;

	x=10000;
	m=complex<double>(1.5,-1.0);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS(0);
	s10=mie1->getS1();
	s20=s10;
	mie1->calcS(pi);
	s1pi=mie1->getS1();
	s2pi=mie1->getS2();
	cout << "Case j: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl;
	cout << "Amplitudes: S1(pi)=" << s1pi << " S2(pi)=" << s2pi << endl << endl;
	delete mie1;	
	
	x=1;
	m=complex<double>(10.0,-10.0);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS(0);
	s10=mie1->getS1();
	s20=s10;
	mie1->calcS(pi);
	s1pi=mie1->getS1();
	s2pi=mie1->getS2();
	cout << "Case k: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl;
	cout << "Amplitudes: S1(pi)=" << s1pi << " S2(pi)=" << s2pi << endl << endl;
	delete mie1;	

	x=100;
	m=complex<double>(10.0,-10.0);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS(0);
	s10=mie1->getS1();
	s20=s10;
	mie1->calcS(pi);
	s1pi=mie1->getS1();
	s2pi=mie1->getS2();
	cout << "Case l: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl;
	cout << "Amplitudes: S1(pi)=" << s1pi << " S2(pi)=" << s2pi << endl << endl;
	delete mie1;	

	x=10000;
	m=complex<double>(10.0,-10.0);
	mie1 = new Mie(x,m);
	mie1->calc();
	qext=mie1->getQext(); qsca=mie1->getQsca();
	mie1->calcS(0);
	s10=mie1->getS1();
	s20=s10;
	mie1->calcS(pi);
	s1pi=mie1->getS1();
	s2pi=mie1->getS2();
	cout << "Case m: x=" << x << " m=" << m << " n=" << mie1->getNmax() << " Qext=" << qext << " Qsca=" << qsca << endl;
	cout << "Amplitudes: S1(0)=" << s10 << " S2(0)=" << s20 << endl;
	cout << "Amplitudes: S1(pi)=" << s1pi << " S2(pi)=" << s2pi << endl << endl;
	delete mie1;
}

