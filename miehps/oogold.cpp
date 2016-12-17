#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include <limits>
#include <string>
#include <sstream>

#include "miehps.h"

const long double pi=3.1415926535897932384626433832795028841971693993751058209749;

using namespace std;

int main(){
	long double a, thetafrom, thetato, thetadelta;
	int ntheta;
	long double qext,qsca,qabs,i1,i2;
	complex<long double> s1(0,0);
	complex<long double> s2(0,0);

	double m2=1.5;
	
	string ftmp;
	ostringstream oss;
	
	int s=1;
	
	fstream fint,fq;
        cout << "Enter filename to save intensity values (without .data): ";
	cin >> ftmp;

	cout << "Enter radius in m: ";
	cin >> a;

	cout.precision(numeric_limits<double>::digits10 + 1);
	cout.precision(7);

	Mie *mie;

	cout << "------ Calculation of a gold sphere with radius a=" << a << "m ------" << endl << endl;

	// set the parameters of the droplet
	
	const int nlambda=15;
	
	long  double lambda[nlambda]={387.5E-9,400.0E-9,413.3E-9,427.5E-9,442.8E-9,
				459.2E-9,476.9E-9,495.9E-9,516.6E-9,539.1E-9,
				563.6E-9,652.6E-9,688.8E-9,729.3E-9,774.9E-9};
		
	for(int i=0;i<nlambda;i++){
		lambda[i]=lambda[i]/m2;
	}
	
	complex<double> m[nlambda];
	
	m[0]=complex<double>(1.674,-1.936)/m2;
	m[1]=complex<double>(1.658,-1.956)/m2;
	m[2]=complex<double>(1.636,-1.958)/m2;
	m[3]=complex<double>(1.616,-1.940)/m2;
	m[4]=complex<double>(1.562,-1.904)/m2;
	m[5]=complex<double>(1.426,-1.846)/m2;
	m[6]=complex<double>(1.242,-1.796)/m2;
	m[7]=complex<double>(0.916,-1.840)/m2;
	m[8]=complex<double>(0.608,-2.120)/m2;
	m[9]=complex<double>(0.402,-2.540)/m2;
	m[10]=complex<double>(0.306,-2.88)/m2;
	m[11]=complex<double>(0.166,-3.15)/m2;
	m[12]=complex<double>(0.160,-3.80)/m2;
	m[13]=complex<double>(0.164,-4.35)/m2;
	m[14]=complex<double>(0.174,-4.86)/m2;
	
	long  double x[nlambda]; // size parameter

	for(int i=0;i<nlambda;i++){
		x[i]=2*pi*a/lambda[i]; // size parameter
	}
	
	long double i0=0;
	
	// range of theta to calculate	
	thetafrom=0;
	thetato=pi;
	thetadelta=0.005;
	ntheta=(thetato-thetafrom)/thetadelta;

	cout << "theta: from=" << thetafrom << " to=" << thetato << " delta=" << thetadelta << endl << endl;
	
	fq.open("q.data",ios_base::out);
		
	for(int i=0;i<nlambda;i++){
		oss.str("");
		oss << ftmp << i << ".data";
		string fname(oss.str());
		fint.open(fname.c_str(),ios_base:: out);

		mie = new Mie(x[i],m[i]);

		mie->calc();
	
		fint << "# a=" << a << " m=" << m[i] << " lambda=" << lambda[i]*m2 << endl;
		fint << "# theta | unpolarized scattering intensity" << endl;
		cout << "lambda=" << lambda[i]*m2  << " x=" << mie->getX() << " m=" << mie->getM() << " n=" << mie->getNmax() << endl;

		cout << "lambda=" << lambda[i]*m2  << " Qext=" << mie->getQext() << " Qsca=" << mie->getQsca() 
			<< " Qabs=" << mie->getQabs() << endl;

		mie->calcS(0.001);

		i0=abs(mie->getS1())*abs(mie->getS1()) + abs(mie->getS2())*abs(mie->getS2());
		cout << "Si(0)=" << i0 << endl;
		fint << 0.0
			<< " " << i0
			<< endl;

		i0=1.0;
		for(int j=1;j<ntheta;j++){
			mie->calcS(thetafrom + (double)j*thetadelta);
	
			fint << thetafrom + (double)j*thetadelta
				<< " " << (abs(mie->getS1())*abs(mie->getS1()) + abs(mie->getS2())*abs(mie->getS2()))/i0
				<< endl;
		}
		
		fq << lambda[i]*m2 << " " << mie->getQext() << " " << mie->getQsca() << " " << mie->getQext() - mie->getQsca() << endl;
		cout << endl;
	
		fint.close();
		delete mie;
	}
}
