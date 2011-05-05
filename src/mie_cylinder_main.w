@** A driver program for the cylindrical Mie scattering code.

@ Here, then, is an overview of document structure
@(mie_cylinder_main.c@>=
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "mie_array.h"
#include "mie_complex.h"
#include "mie_cylinder.h"
#define  PI 3.14159265358979

@<Definition for |mie_cylinder_main|@>@;

@ @<Definition for |mie_cylinder_main|@>=
int main(int argc,char**argv)
{
	struct c_complex m, *t1, *t2, *t3, qexpar, qexper;
	double qscpar,qscper,x,zeta,lambda,diameter,qsca,t1norm,*theta;
	int i,nangles;
	long size;
	double n_medium, n_cylinder, k_cylinder;
	
	nangles=21;
	size = nangles;
	t1 = new_carray(size);
	t2 = new_carray(size);
	t3 = new_carray(size);
	theta = new_darray(size);
	
	for (i=0; i<nangles; i++)
		theta[i] = (double)i*PI/(nangles-1.0);

	lambda = 632.8;
	diameter = 1050;
	x = PI * diameter/lambda;
	n_cylinder = 1.55;
	n_medium   = 1.0;
	k_cylinder = 0.0;
	m=c_set(n_cylinder/n_medium,k_cylinder/n_medium);
	zeta = 90.0 * PI / 180.0;
	
	MieCylinder(x,m,zeta,theta,nangles,t1,t2,t3,&qexpar,&qexper,&qscpar,&qscper);

	qsca  = 0.5*(qscpar+qscper);
	
	printf("zeta       = %8.4f degrees\n", zeta * 180/PI);
	printf("n_medium   = %8.4f\n",n_medium);
	printf("n_cylinder = %8.4f + %8.4fi\n",n_cylinder,k_cylinder);
	printf("lambda     = %8.1f nm\n",lambda);
	printf("diameter   = %8.1f nm\n",diameter);
	printf("\n");
	printf("x          = %8.4f\n", x);
	printf("qexpar     = %8.4f + %8.4fi\n",qexpar.re,qexpar.im);
	printf("qscpar     = %8.4f\n",qscpar);
	printf("qexper     = %8.4f + %8.4fi\n",qexper.re,qexper.im);
	printf("qscper     = %8.4f\n",qscper);
	printf("\n");
	
	t1norm = 0.5 * c_norm(t1[0]) + 0.5*c_norm(t2[0]);
	for (i=0; i<nangles; i++) {
		double tpar,tper,t11,t12,pol,s1,s2,s3,t33,t34;
		tpar = c_norm(t1[i]);
		tper = c_norm(t2[i]);
		t11 = 0.5*(tpar+tper);
		t12 = 0.5*(tpar-tper);
		t33 = t1[i].re*t2[i].re+t1[i].im*t2[i].im;
		t34 = t1[i].im*t2[i].re-t1[i].re*t2[i].im;
		t33 /= t11;
		t34 /= t11;
		pol = t12/t11;
		
		s1=4*c_norm(t1[i])/qsca;
		s2=4*c_norm(t2[i])/qsca;
		s3=4*c_norm(t3[i])/qsca;
		printf("%8.3f \t %8.5f \t %8.5f \t %8.5f \t %8.5f\n",
		        theta[i]*180/PI,t11/t1norm,pol,t33,t34);
	}

	free_carray(t1);
	free_carray(t2);
	free_carray(t3);
	return(0);
}
