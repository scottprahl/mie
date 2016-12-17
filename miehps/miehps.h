/*  Copyright:
    Copyright (C) 2011 Hans-Peter Schadler <hps@abyle.org>

    License:
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#ifndef _MIEHPS_H
#define _MIEHPS_H

#include <complex>

using namespace std;

class Mie
{
	public:
		Mie(double xt, complex<double> mt);
		~Mie();

		void calc();
		void calcS(double theta);

		double getQext();
		double getQsca();
		double getQabs();
		complex<double> getS1();
		complex<double> getS2();

		int getNmax();
		complex<double> getM();
		double getX();

	private:
		int nm,n;
		complex<double> m, mx;
		double x;
		int Nstar;

		complex<double> *an;
		complex<double> *bn;
		complex<double> *rn;
		complex<double> *rntest;
		complex<double> ccotmx;
		double *chi;
		double *psi;
		double *pi;
		double qext;
		double qsca;
		double qabs;
		complex<double> s1;
		complex<double> s2;
		
		complex<double> Cot(complex<double> mx);
		
		double calcPi(double mu);
		double Tau(double mu,int n);
		void calcChi(double x,int n);
		void calcPsi(double x,int n);
		complex<double> Zeta(double x,int n);
		
		void calcR();
		void calcAnBn();

		complex<double> calcSi(double theta);
		void calcQscaQext();
};

#endif /* _MIEHPS_H */
