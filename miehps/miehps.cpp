/*  Copyright:
    Copyright (C) 2011 Hans-Peter Schadler <hps@abyle.org>
    Based on the algorithm described in the paper
    Hong Du, ``Mie-Scattering Calculation'', Appl. Opt. 43, 1951-1956 (2004)
    
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

#include <cmath>
#include <complex>

#include "miehps.h"

using namespace std;

Mie::Mie(double xt, complex<double> mt){
	x=xt;
	m=mt;
	mx=m*x;
	Nstar=60000;
	n=x+4.0*pow(x,1.0/3.0)+2.0+10; // when to cut of the summations (term taken from Hong Du)
	an=new complex<double>[n];
	bn=new complex<double>[n];
	rn=new complex<double>[Nstar+1];
	chi=new double[n+1];
	psi=new double[n+1];
	pi=new double[n+1];
	calcChi(x,n);
	calcPsi(x,n);
	ccotmx=conj(Cot(conj(mx)));
}

Mie::~Mie(){

}

void Mie::calc(){
	calcR();
	calcAnBn();
	calcQscaQext();
}

double Mie::getQext(){
	return qext;
}

double Mie::getQsca(){
	return qsca;
}

double Mie::getQabs(){
	return qext-qsca;
}

void Mie::calcS(double theta){
	calcSi(theta);
}

complex<double> Mie::getS1(){
	return s1;
}

complex<double> Mie::getS2(){
	return s2;
}
		
int Mie::getNmax(){
	return n;
}

complex<double> Mie::getM(){
	return m;
}

double Mie::getX(){
	return x;
}

double Mie::calcPi(double mu){
	double s=0.0,t=0.0;
	
	pi[0]=0.0;
	pi[1]=1.0;
	
	for(int i=1;i<n;i++){
		s=mu*pi[i];
		t=s-pi[i-1];
		pi[i+1]=s+t+t/(double)i;
	}
}

double Mie::Tau(double mu,int n){
	if(n==0){
		return 0.0;
	}else{
		return n*(mu*pi[n]-pi[n-1])-pi[n-1];
	}
}

void Mie::calcChi(double x,int n){
	chi[0]=cos(x);
	chi[1]=chi[0]/x+sin(x);	

	for(int i=1;i<n;i++){
		chi[i+1]=(2.0*i+1.0)*chi[i]/x-chi[i-1];
	}
}

void Mie::calcPsi(double x,int n){
	psi[0]=sin(x);
	psi[1]=psi[0]/x-cos(x);

	for(int i=1;i<n;i++){
		psi[i+1]=(2.0*i+1.0)*psi[i]/x-psi[i-1];
	}
}

complex<double> Mie::Zeta(double x,int n){
	return complex<double>(psi[n],chi[n]);
}

complex<double> Mie::Cot(complex<double> mx){
	complex<double> cot(0,0),num,denom;
	double realnum,imagnum,realdenom,imagdenom,realcot,imagcot,div,ctan,cexp;
	
	ctan=tan(real(mx));
	cexp=exp(-2.0*imag(mx));
	
	num=complex<double>(0,1)+ctan-cexp*ctan+complex<double>(0,1.0)*cexp;
	denom=-(double)1.0+complex<double>(0,1.0)*ctan+complex<double>(0,1.0)*cexp*ctan+cexp;
	
	realnum=real(num);
	imagnum=imag(num);
	realdenom=real(denom);
	imagdenom=imag(denom);
	
	div=(realdenom*realdenom+imagdenom*imagdenom);
	
	realcot=((realnum*realdenom+imagnum*imagdenom)/div);
	imagcot=((imagnum*realdenom-realnum*imagdenom)/div);

	return complex<double>(realcot,imagcot);
}

void Mie::calcR(){ //downward recurrence
	rn[Nstar-1]=(double)(2.0*Nstar+1.0)/mx;
	
	for(int i=Nstar-1;i>=1;i--){
		rn[i-1]=(2.0*i+1.0)/mx-1.0/rn[i];
	}
}

void Mie::calcAnBn(){
	complex<double> rf;
	
	for(int nc=1;nc<=n;nc++){
		rf=(rn[nc-1]/m+(double)nc*(1.0-1.0/(m*m))/x);
		
		an[nc-1]=(rf*psi[nc]-psi[nc-1])/(rf*Zeta(x,nc)-Zeta(x,nc-1));
		bn[nc-1]=((complex<double>)rn[nc-1]*(complex<double>)m*psi[nc]-psi[nc-1])/((complex<double>)rn[nc-1]*(complex<double>)m*Zeta(x,nc)-Zeta(x,nc-1));
	}
}

complex<double> Mie::calcSi(double theta){
	double mu=cos(theta);
	double tau;
	double fn;

	s1=complex<double>(0,0);
	s2=complex<double>(0,0);
	
	calcPi(mu);
	
	for(int i=1;i<=n;i++){
		tau=Tau(mu,i);
		fn=(double)(2.0*i+1.0)/(double)(i*(i+1.0));
		s1=s1+fn*(an[i-1]*pi[i]+bn[i-1]*tau);
		s2=s2+fn*(an[i-1]*tau+bn[i-1]*pi[i]);
	}
}

void Mie::calcQscaQext(){
	qsca=0;
	qext=0;
		
	for(int i=1;i<=n;i++){
		qext=qext+(double)(2.0*i+1.0)*(an[i-1].real()+bn[i-1].real());
		qsca=qsca+(double)(2.0*i+1.0)*(pow(abs(an[i-1]),2)+pow(abs(bn[i-1]),2));
	}

	qext=qext*2.0/(x*x);
	qsca=qsca*2.0/(x*x);
}
