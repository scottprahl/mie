@** A Driver Program to test Mie Scattering.

@
@(test_mie.c@>=
@<Definition for |MieTest|@>@;

@*1 Mie Testing.

@ Here is the obligatory program to test the Mie Scattering Program.

@<Definition for |MieTest|@>=

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "mie_array.h"
#include "mie_complex.h"
#include "mie_lobatto.h"
#include "mie.h"
	
int main() 
{
	@<Test Logarithmic derivative@>@;
	@<First Mie Test@>@;
	@<Second Mie Test@>@;
	@<Third Mie Test@>@;
	@<Fourth Mie Test@>@;
	@<Fifth Mie Test@>@;
	@<Sixth Mie Test@>@;
	@<Test Small Mie@>@;
	@<Test Backscattering@>@;
	return 0;
}

@ First, I make sure that the logarithmic derivatives are 
calculated correctly.
@<Test Logarithmic derivative@>=
{
	struct c_complex y, z, *D;
	long nstop;
	
	printf("\n***********************************************\n");
	printf("Zeroth test for logarithmic derivative\n");
	printf("   The result should for D_9(1.0) = 9.95228198\n");
	z = c_set(1.0,0.0);
	y = Lentz_Dn(z,9L);
	printf("   The actual value               = %11.8f +i%12.8f\n\n",y.re, y.im);
	
	z= c_set(62*1.28,-62*1.37);
	nstop = 50;
	
	D = new_carray(nstop);
	printf("   For n = %ld \n", nstop);
	printf("   For j = %ld \n", 10L);
	printf("   For z = %10.5f +i %10.5f\n",z.re,z.im);
	
	printf("   Mathematica                 gives %10.6f +i %10.6f\n",0.004087,1.0002620);
	y = Lentz_Dn(z,10L);
	printf("   Dn[10] continued fraction   gives %10.6f +i %10.6f\n",y.re,y.im);
	
	Dn_up(z, nstop, D);
	printf("   Dn[10] upwards recurrence   gives %10.6f +i %10.6f\n",D[10].re,D[10].im);
	
	Dn_down(z, nstop, D);
	printf("   Dn[10] downwards recurrence gives %10.6f +i %10.6f\n",D[10].re,D[10].im);
	
	free_carray(D);
}

@ I did not get around to adding this test for a long time.  Finally, I did
a check and discovered that it failed for absorbing spheres!  I used Wiscombe's
test data and 
@<Test Small Mie@>=
{
	long nangles = 0;
	struct c_complex *s1 = NULL;
	struct c_complex *s2 = NULL;
	double *mu = NULL;
	struct c_complex m;
	double x= 20;
	double qext, qsca, qback, g;
	
	printf("\n***********************************************\n");
	printf("Small Mie Test\n");
	printf("           calc    Wiscombe       calc   Wiscombe\n");
	printf(" X         Qsca    Qsca            g         g\n");
	
	m=c_set(0.75, 0.0);
	x=0.099;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n",x,qext,0.000007,g,0.001448);
	x=0.101;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n",x,qext,0.000008,g,0.001507);

	m=c_set(1.5, -1.0);
	x=0.055;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n",x,qext,0.101491,g,0.000491);
	x=0.056;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n",x,qext,0.103347,g,0.000509);
	
	m=c_set(1e-10, -1e10);
	x=0.099;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n",x,qext,0.000321,g,-0.397357);
	x=0.101;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n",x,qext,0.000348,g,-0.397262);

	m=c_set(0.0, -1e10);
	x=0.099;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n",x,qext,0.000321,g,-0.397357);
	x=0.101;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.3f   % 8.6f % 8.6f   % 8.6f % 8.6f\n",x,qext,0.000348,g,-0.397262);
	
	printf("\n");
}

@ For the first Mie test, I chose the sample output values from the
book by Bohren and Huffman.  This is the only test that includes $Q_{\rm back}$.
The angle stuff is not used. 
@<First Mie Test@>=
{
    long nangles, i;
    struct c_complex *s1, *s2, m;
    double *w, *mu, x, rho, qext, qsca, qback, g;

	nangles = 10;
	m.re    = 1.55;
	m.im    = 0.0;
	x       = 5.213;
	rho     = 2 * x * (m.re - 1);
	
	printf("\n***********************************************\n");
	printf("First Mie Test -- cf. Bohren and Huffman pg 482\n");
	printf("    index of medium      %7.4f\n", 1.0);
	printf("    real index of sphere %7.4f\n", m.re);
	printf("    imag index of sphere %7.4f\n", m.im);
	printf("    quadrature angles    %ld\n\n", nangles);
	
	s1 = new_carray(nangles);
	s2 = new_carray(nangles);
	mu = new_darray(nangles);
	w  = new_darray(nangles);

	Lobatto (0.0, 3.1415926535, mu, w, nangles);
	for (i = 0; i < nangles; ++i)
	  mu[i] = cos(mu[i]);
	
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	
	printf("    x        Qsca       Qext      Qback    g\n");
	printf( "%7.3f %10.6f %10.6f %10.6f \n", 5.213, 3.10543, 3.10543, 2.92534);
	printf( "%7.3f %10.6f %10.6f %10.6f %10.6f \n", x, qsca, qext, qback, g);
	
	free_darray(mu);
	free_darray(w);
	free_carray(s1);
	free_carray(s2);
}

@ Another test that includes absorbing spheres is from the paper by Dave.
His table 2 tabulates the absorption efficiency.
@<Second Mie Test@>=
{
	double pi=3.14159265358979;
	long nangles = 0;
	struct c_complex *s1 = NULL;
	struct c_complex *s2 = NULL;
	double *mu = NULL;
	struct c_complex m;
	double x= 50.0 * pi;
	double qext, qsca, qback, g;
	
	printf("\n***********************************************\n");
	printf("Second Mie Test -- Dave Table 2\n");
	printf("          n                 Qa            Dave\n");
	
	m=c_set(1.342, 0.0);
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%10.5g%-+7.4fi    %10.5f    %10.5f\n", m.re,m.im, qext-qsca, 0.0);
	
	m=c_set(1.342, -0.0001);
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%10.5g%-+7.4fi    %10.5f    %10.5f\n", m.re,m.im, qext-qsca, 0.0535);
	
	m=c_set(1.342, -0.01);
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%10.5g%-+7.4fi    %10.5f    %10.5f\n", m.re,m.im, qext-qsca, 0.9649);
	
	m=c_set(1.342, -0.2);
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%10.5g%-+7.4fi    %10.5f    %10.5f\n", m.re,m.im, qext-qsca, 0.9542);
	
	m=c_set(1.342, -0.4);
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%10.5g%-+7.4fi    %10.5f    %10.5f\n", m.re,m.im, qext-qsca, 0.9221);
	
	m=c_set(1.342, -0.6);
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%10.5g%-+7.4fi    %10.5f    %10.5f\n", m.re,m.im, qext-qsca, 0.8808);
	
	m=c_set(1.342, -0.8);
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%10.5g%-+7.4fi    %10.5f    %10.5f\n", m.re,m.im, qext-qsca, 0.8369);
	
	m=c_set(1.342, -1.0);
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%10.5g%-+7.4fi    %10.5f    %10.5f\n", m.re,m.im, qext-qsca, 0.7910);
	
	printf("\n");
}


@ The third test uses tabulated values from van de Hulst to check the 
anisotropy calculation.
@<Third Mie Test@>=
{
	long nangles = 0;
	struct c_complex *s1 = NULL;
	struct c_complex *s2 = NULL;
	double *mu = NULL;
	struct c_complex m;
	double x= 20;
	double qext, qsca, qback, g;
	
	printf("\n***********************************************\n");
	printf("Third Mie Test -- van de Hulst page 161\n");
	printf(" x          Qs    vdH       Qs*g   vdH\n");
	
	m=c_set(0.0, 0.0);
	x=0.3;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.1f   %7.3f %7.3f   %7.3f %7.3f\n",x,qsca,0.028,qsca*g,-0.011);
	ez_Mie(x,0.0, &qsca, &g);
	printf("%4.1f   %7.3f %7.3f   %7.3f %7.3f\n",x,qsca,0.028,qsca*g,-0.011);
	
	x=1.0;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.1f   %7.3f %7.3f   %7.3f %7.3f\n",x,qsca,2.036,qsca*g,-0.385);
	
	x=1.5;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.1f   %7.3f %7.3f   %7.3f %7.3f\n",x,qsca,2.155,qsca*g,0.156);
	
	x=5.0;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.1f   %7.3f %7.3f   %7.3f %7.3f\n",x,qsca,2.116,qsca*g,0.965);
	
	printf("\n");
}

@ This test includes one value for totally reflecting spheres.
@<Fourth Mie Test@>=
{
	long nangles = 0;
	struct c_complex *s1 = NULL;
	struct c_complex *s2 = NULL;
	double *mu = NULL;
	struct c_complex m;
	double x= 20;
	double qext, qsca, qback, g;
	
	printf("\n***********************************************\n");
	printf("Fourth Mie Test -- van de Hulst page 277\n");
	printf(" x          Qs    vdH       Qs*g   vdH\n");
	
	m=c_set(3.41, -1.94);
	x=1.3;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.1f   %7.3f %7.3f   %7.2f %7.2f\n",x,qsca,1.669,qsca*g,0.30);
	
	m=c_set(7.20, -2.65);
	x=1.3;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.1f   %7.3f %7.3f   %7.2f %7.2f\n",x,qsca,1.860,qsca*g,0.31);
	
	m=c_set(0.0, 0.0);
	x=1.3;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%4.1f   %7.3f %7.3f   %7.2f %7.2f\n",x,qsca,2.266,qsca*g,-0.05);
	
	printf("\n");
}

@ Well I found a complete set of results posted by Wiscombe.  This
allows complete code coverage.  
@<Fifth Mie Test@>=
{
	long nangles = 0;
	struct c_complex *s1 = NULL;
	struct c_complex *s2 = NULL;
	double *mu = NULL;
	struct c_complex m;
	double x;
	double qext, qsca, qback, g;
	
	printf("\n***********************************************\n");
	printf("Fifth Mie Test -- Wiscombe\n");
	
	@<Wiscombe Non-absorbing spheres@>@;
	@<Wiscombe Absorbing water spheres@>@;
	@<Wiscombe Absorbing spheres@>@;
	@<Wiscombe Yet More Absorbing spheres@>@;
	@<Wiscombe perfectly conducting spheres@>@;
	printf("\n");
}

@ @<Wiscombe Non-absorbing spheres@>=
	printf("\nNon-Absorbing Spheres m=(0.75+0.0i)\n");
	printf("               Calc.     Wiscombe     Calc     Wiscombe\n");
	printf("   x            Qs          Qs          g          g\n");
	m=c_set(0.75, 0.0);
	x=0.099;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,0.000007,g,0.001448);
	ez_Mie(x,0.75, &qsca, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,0.000007,g,0.001448);
	
	m=c_set(0.75, 0.0);
	x=0.101;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,0.000008 ,g,0.001507);
	
	m=c_set(0.75, 0.0);
	x=10.0;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,2.232265,g,0.896473);
	
	m=c_set(0.75, 0.0);
	x=1000.0;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,1.997908,g,0.844944);
	ez_Mie(x,0.75, &qsca, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,1.997908,g,0.844944);

@ @<Wiscombe Absorbing water spheres@>=
	printf("\nAbsorbing Water Spheres m=(1.33-0.00001i)\n");
	printf("               Calc.     Wiscombe     Calc     Wiscombe\n");
	printf("   x            Qs          Qs          g          g\n");
	m=c_set(1.33, -0.00001);
	x=1.0;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,0.093923,g,0.184517);
	
	x=100.0;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,2.096594,g,0.868959);
	
	x=10000.0;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,1.723857,g,0.907840);
	
	@ @<Wiscombe Absorbing spheres@>=
	printf("\nAbsorbing Spheres m=(1.5-i)\n");
	printf("               Calc.     Wiscombe     Calc     Wiscombe\n");
	printf("   x            Qs          Qs          g          g\n");
	m=c_set(1.5, -1.00);
	
	x=0.055;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,0.000011,g,0.000491);
	
	x=0.056;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,0.000012,g,0.000509);
	
	x=1.0;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,0.6634538,g,0.192136);
	
	x=100.0;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,1.283697,g,0.850252);
	
	x=10000.0;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,1.236574,g,0.846310);

@ @<Wiscombe Yet More Absorbing spheres@>=
	printf("\n Yet More Absorbing Spheres m=(10-10i)\n");
	printf("               Calc.     Wiscombe     Calc     Wiscombe\n");
	printf("   x            Qs          Qs          g          g\n");
	m=c_set(10.0, -10.00);
	
	x=1.0;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,2.049405,g,-0.110664);
	
	x=100.0;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,1.836785,g,0.556215);
	
	x=10000.0;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,1.795393,g,0.548194);

@ @<Wiscombe perfectly conducting spheres@>=
	printf("\nPerfectly Conducting Spheres\n");
	printf("               Calc.     Wiscombe     Calc     Wiscombe\n");
	printf("   x            Qs          Qs          g          g\n");
	m=c_set(0.0, 0.0);
	x=0.099;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,0.000321,g,-0.397357);
	
	x=0.101;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,0.000348,g,-0.397262);
	
	x=100;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,2.008102,g,0.500926);
	ez_Mie(x,0.0, &qsca, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,2.008102,g,0.500926);
	
	x=10000;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%9.3f   %11.7f %11.7f   %11.7f %11.7f\n",x,qsca,2.000289,g,0.500070);


@ 	@<Sixth Mie Test@>=
{
	double pi=3.14159265358979;
	long nangles = 7;
	struct c_complex *s1 = NULL;
	struct c_complex *s2 = NULL;
	double *mu = NULL;
	struct c_complex m;
	double x;
	double qext, qsca, qback, g;
	char * form =  "%7.4f %8.5f%+-8.5fi    %8.5f%+-8.5fi  Calc\n";
	char * form2 = "%7.4f %8.5f%+-8.5fi    %8.5f%+-8.5fi   Wiscombe\n\n";
	long i;
	
	s1 = new_carray(nangles);
	s2 = new_carray(nangles);
	mu = new_darray(nangles);
	
	for (i=0; i<nangles; i++) mu[i] = cos(pi*i/6.0);

	printf("\n***********************************************\n");
	printf("Sixth Mie Test -- Wiscombe\n");
	printf("   angle       S1                    S2         \n");
	x=1.0;
	m=c_set(1.5,-1.0);
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	
	printf(form,mu[0],s1[0].re,s1[0].im,s2[0].re,s2[0].im);
	printf(form2,mu[0],5.84080E-01,1.90515E-01,5.84080E-01,1.90515E-01);
	printf(form,mu[1],s1[1].re,s1[1].im,s2[1].re,s2[1].im);
	printf(form2,mu[1],5.65702E-01,1.87200E-01,5.00161E-01,1.45611E-01);
	printf(form,mu[2],s1[2].re,s1[2].im,s2[2].re,s2[2].im);
	printf(form2,mu[2],5.17525E-01,1.78443E-01,2.87964E-01,4.10540E-02);
	printf(form,mu[3],s1[3].re,s1[3].im,s2[3].re,s2[3].im);
	printf(form2,mu[3],4.56340E-01,1.67167E-01,3.62285E-02,-6.18265E-02);
	printf(form,mu[4],s1[4].re,s1[4].im,s2[4].re,s2[4].im);
	printf(form2,mu[4],4.00212E-01,1.56643E-01,-1.74875E-01,-1.22959E-01);
	printf(form,mu[5],s1[5].re,s1[5].im,s2[5].re,s2[5].im);
	printf(form2,mu[5],3.62157E-01,1.49391E-01,-3.05682E-01,-1.43846E-01);
	printf(form,mu[6],s1[6].re,s1[6].im,s2[6].re,s2[6].im);
	printf(form2,mu[6],3.48844E-01,1.46829E-01,-3.48844E-01,-1.46829E-01);
	
	printf("\n");
	
	free_carray(s1);
	free_carray(s2);
	free_darray(mu);
}

@ 
@<Test Backscattering@>=
{
	long nangles = 0;
	struct c_complex *s1 = NULL;
	struct c_complex *s2 = NULL;
	double *mu = NULL;
	struct c_complex m;
	double x,qext, qsca, qback, g,ref;
	
	printf("\n***********************************************\n");
	printf("Backscattering Efficiency\n");
	printf("                                   Calc          reference\n");
	printf("    X         m.re     m.im        Qsca          qsca          ratio\n");
	
	m=c_set(1.55, 0.0);
	x   = 2*3.1415926535*0.525/0.6328;
	ref = 2.92534;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8g\n",x,m.re,m.im,qback,ref,qback/ref);

	m=c_set(0.0, -1000.0);
	x=0.099;
	ref = (4.77373E-07*4.77373E-07 +  1.45416E-03*1.45416E-03)/x/x*4;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.2f   % 8e % 8e %8.5f\n",x,m.re,m.im,qback,ref,qback/ref);
	x=0.101;
	ref = (5.37209E-07*5.37209E-07 +  1.54399E-03*1.54399E-03)/x/x*4;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.2f   % 8e % 8e %8.5f\n",x,m.re,m.im,qback,ref,qback/ref);
	x=100;
	ref = (4.35251E+01*4.35251E+01 +  2.45587E+01*2.45587E+01)/x/x*4;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.2f   % 8e % 8e %8.5f\n",x,m.re,m.im,qback,ref,qback/ref);

	m=c_set(0.75, 0.0);
	x=0.099;
	ref = (1.81756E-08*1.81756E-08 + 1.64810E-04*1.64810E-04)/x/x*4;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n",x,m.re,m.im,qback,ref,qback/ref);
	x=0.101;
	ref = (2.04875E-08*2.04875E-08 + 1.74965E-04*1.74965E-04)/x/x*4;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n",x,m.re,m.im,qback,ref,qback/ref);
	x=10.0;
	ref = (1.07857E+00*1.07857E+00 + 3.60881E-02*3.60881E-02)/x/x*4;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n",x,m.re,m.im,qback,ref,qback/ref);
	x=1000.0;
	ref = (1.70578E+01*1.70578E+01 +  4.84251E+02* 4.84251E+02)/x/x*4;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n",x,m.re,m.im,qback,ref,qback/ref);

	m=c_set(1.33, -0.00001);
	x=1.0;
	ref = (2.24362E-02*2.24362E-02 +  1.43711E-01*1.43711E-01)/x/x*4;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n",x,m.re,m.im,qback,ref,qback/ref);
	x=100.0;
	ref = (5.65921E+01*5.65921E+01 +  4.65097E+01*4.65097E+01)/x/x*4;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n",x,m.re,m.im,qback,ref,qback/ref);
	
	m=c_set(1.5, -1.0);
	x=0.055;
	ref = (7.66140E-05*7.66140E-05 +  8.33814E-05*8.33814E-05)/x/x*4;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n",x,m.re,m.im,qback,ref,qback/ref);
	x=0.056;
	ref = (8.08721E-05*8.08721E-05 +  8.80098E-05*8.80098E-05)/x/x*4;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n",x,m.re,m.im,qback,ref,qback/ref);
	x=1.0;
	ref = (3.48844E-01*3.48844E-01 +  1.46829E-01*1.46829E-01)/x/x*4;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n",x,m.re,m.im,qback,ref,qback/ref);
	x=100.0;
	ref = (2.02936E+01*2.02936E+01 +  4.38444E+00*4.38444E+00)/x/x*4;
	Mie (x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g);
	printf("%8.3f   % 8.4f % 8.4f   % 8e % 8e %8.5f\n",x,m.re,m.im,qback,ref,qback/ref);


	printf("\n");
}
