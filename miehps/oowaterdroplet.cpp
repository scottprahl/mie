#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include <limits>
#include <string>

#include "miehps.h"

const long double pi=3.1415926535897932384626433832795028841971693993751058209749;

using namespace std;

int main(){
//	double x;
	int n;
	long double a, thetafrom, thetato, thetadelta;
	int ntheta;
	long double qext,qsca,qabs,i1,i2;
	complex<long double> s1(0,0);
	complex<long double> s2(0,0);
	complex<double> m1(0,0);
	complex<double> m2(0,0);
	complex<double> m3(0,0);
	complex<double> m4(0,0);
	
	string fwatername;
	
	fstream fwater;
	fwatername="water.data";
	fwater.open(fwatername.c_str(),ios_base:: out);
		
	Mie *mie1;
	Mie *mie2;
	Mie *mie3;
	Mie *mie4;
	
	cout.precision(numeric_limits<double>::digits10 + 1);
	cout.precision(7);

	cout << "------ Calculations of a water droplet with radius a=1mm ------" << endl << endl;

	// set the parameters of the droplet
	
	long  double lambda1=400E-9; // wavelength to calculate
	long  double lambda2=500E-9; // wavelength to calculate
	long  double lambda3=600E-9; // wavelength to calculate
	long  double lambda4=700E-9; // wavelength to calculate
	a=1E-3; // radius of the droplet in m
/*	m1=complex<double>(1.339,-1.68E-09); // complex refraction index, imaginary part has to be negative
	m2=complex<double>(1.339,-9.24E-10);
	m3=complex<double>(1.333,-9.63E-09);
	m4=complex<double>(1.329,-3.35E-08); */
	
	m1=complex<double>(1.339,-1.86E-09); // complex refraction index, imaginary part has to be negative
	m2=complex<double>(1.335,-1.00E-09);
	m3=complex<double>(1.332,-1.09E-08);
	m4=complex<double>(1.331,-3.35E-08);

	long  double x1=2*pi*a/lambda1; // size parameter
	long  double x2=2*pi*a/lambda2; // size parameter
	long  double x3=2*pi*a/lambda3; // size parameter
	long  double x4=2*pi*a/lambda4; // size parameter

	// range of theta to calculate	
	thetafrom=0;
// 	thetafrom=2.2; // wasser neben
// 	thetafrom=2.35; // wasser haupt
	thetato=pi;
//	thetato=2.325; // wasser neben
//	thetato=2.45; // wasser haupt
	thetadelta=0.0001;
	ntheta=(thetato-thetafrom)/thetadelta;

	//prepare the calculations
	mie1 = new Mie(x1,m1);
	mie2 = new Mie(x2,m2);
	mie3 = new Mie(x3,m3);
	mie4 = new Mie(x4,m4);

	mie1->calc();
	mie2->calc();
	mie3->calc();
	mie4->calc();

	fwater << "# theta | unpolarized scattering intensity" << endl;
	cout << "lambda1=" << lambda1 << " x=" << mie1->getX() << " m=" << mie1->getM() << " n=" << mie1->getNmax() << endl;
	cout << "lambda2=" << lambda2 << " x=" << mie2->getX() << " m=" << mie2->getM() << " n=" << mie2->getNmax() << endl;
	cout << "lambda3=" << lambda3 << " x=" << mie3->getX() << " m=" << mie3->getM() << " n=" << mie3->getNmax() << endl;
	cout << "lambda4=" << lambda4 << " x=" << mie4->getX() << " m=" << mie4->getM() << " n=" << mie4->getNmax() << endl;
	cout << "theta: from=" << thetafrom << " to=" << thetato << " delta=" << thetadelta << " Number of angles to calc=" << ntheta << endl << endl;

	qext=mie1->getQext();
	qsca=mie1->getQsca();
	qabs=mie1->getQabs();

	double N=1000.0;
	double dist=1000.0;
	
	double gamma1,gamma2,gamma3,gamma4;
	
	gamma1=mie1->getQext()*a*a*N*pi;
	gamma2=mie2->getQext()*a*a*N*pi;
	gamma3=mie3->getQext()*a*a*N*pi;
	gamma4=mie4->getQext()*a*a*N*pi;

	cout << "lambda=" << lambda1 << " Qext=" << mie1->getQext() << " Qsca=" << mie1->getQsca() 
		<< " Qabs=" << mie1->getQabs() << endl;
	cout << "extinction coefficient=" << mie1->getQext()*a*a*N*pi << endl;
	cout << "extinction in " << dist << "m distance=" << exp(-mie1->getQext()*a*a*N*pi*1000.0) << endl << endl;
	cout << "lambda=" << lambda2 << " Qext=" << mie2->getQext() << " Qsca=" << mie2->getQsca() 
		<< " Qabs=" << mie2->getQabs() << endl;
	cout << "extinction coefficient=" << mie2->getQext()*a*a*N*pi << endl;
	cout << "extinction in " << dist << "m distance=" << exp(-mie2->getQext()*a*a*N*pi*1000.0) << endl << endl;
	cout << "lambda=" << lambda3 << " Qext=" << mie3->getQext() << " Qsca=" << mie3->getQsca() 
		<< " Qabs=" << mie3->getQabs() << endl;		
	cout << "extinction coefficient=" << mie2->getQext()*a*a*N*pi << endl;
	cout << "extinction in " << dist << "m distance=" << exp(-mie3->getQext()*a*a*N*pi*1000.0) << endl << endl;
	cout << "lambda=" << lambda4 << " Qext=" << mie4->getQext() << " Qsca=" << mie4->getQsca() 
		<< " Qabs=" << mie4->getQabs() << endl;
	cout << "extinction coefficient=" << mie2->getQext()*a*a*N*pi << endl; 
	cout << "extinction in " << dist << "m distance=" << exp(-mie4->getQext()*a*a*N*pi*1000.0) << endl << endl;

	double I0=1300.0;

//	double i01,i02,i03,i04;
/*	if(false){	
	// treat theta=0 separately
	i01=abs(mie1->getSi0())*abs(mie1->getSi0()) + abs(mie1->getSi0())*abs(mie1->getSi0());
	i02=abs(mie2->getSi0())*abs(mie2->getSi0()) + abs(mie2->getSi0())*abs(mie2->getSi0());
	i03=abs(mie3->getSi0())*abs(mie3->getSi0()) + abs(mie3->getSi0())*abs(mie3->getSi0());
	i04=abs(mie4->getSi0())*abs(mie4->getSi0()) + abs(mie4->getSi0())*abs(mie4->getSi0());
	
	fwater << 0.0 
		<< " " << I0*i01*exp(-gamma1*dist)*2*pi/(2.0*dist*dist*lambda1)
		<< " " << I0*i02*exp(-gamma2*dist)*2*pi/(2.0*dist*dist*lambda2)
		<< " " << I0*i03*exp(-gamma3*dist)*2*pi/(2.0*dist*dist*lambda3)
		<< " " << I0*i04*exp(-gamma4*dist)*2*pi/(2.0*dist*dist*lambda4)
		<< endl;
	}*/
	double i1i21,i1i22,i1i23,i1i24;

	double norm=(2.0*dist*dist*2.0*pi*2.0*pi);
	
	for(int i=0;i<ntheta;i++){
		mie1->calcS(thetafrom + (double)i*thetadelta);
		mie2->calcS(thetafrom + (double)i*thetadelta);
		mie3->calcS(thetafrom + (double)i*thetadelta);
		mie4->calcS(thetafrom + (double)i*thetadelta);
	
		i1i21=(abs(mie1->getS1())*abs(mie1->getS1()) + abs(mie1->getS2())*abs(mie1->getS2()));
		i1i22=(abs(mie2->getS1())*abs(mie2->getS1()) + abs(mie2->getS2())*abs(mie2->getS2()));
		i1i23=(abs(mie3->getS1())*abs(mie3->getS1()) + abs(mie3->getS2())*abs(mie3->getS2()));
		i1i24=(abs(mie4->getS1())*abs(mie4->getS1()) + abs(mie4->getS2())*abs(mie4->getS2()));
		
	
		fwater << thetafrom + (double)i*thetadelta
			<< " " << I0*i1i21*exp(-gamma1*dist)*lambda1*lambda1/norm
			<< " " << I0*i1i22*exp(-gamma2*dist)*lambda2*lambda2/norm
			<< " " << I0*i1i23*exp(-gamma3*dist)*lambda3*lambda3/norm
			<< " " << I0*i1i24*exp(-gamma4*dist)*lambda4*lambda4/norm
			<< endl;
	}
	
	fwater.close();
}
