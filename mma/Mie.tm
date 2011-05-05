#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "libmie.h"
#include "mathlink.h"

static double * new_darray(long size)
{
  double *a;
  size += 2;
  a = (double *) malloc(sizeof(double)*size);
  return a+1;
}

static void free_darray(double *a)
{
  if (a!=0) free(a-1);
}

static void error(char *s)
{
	MLEvaluate(stdlink,s);
	MLNextPacket(stdlink);
	MLNewPacket(stdlink);
  	MLPutFunction(stdlink, "List", 5);
  	MLPutFunction(stdlink, "List", 4);
	MLPutSymbol(stdlink,"$Failed");
	MLPutSymbol(stdlink,"$Failed");
	MLPutSymbol(stdlink,"$Failed");
	MLPutSymbol(stdlink,"$Failed");
  	MLPutFunction(stdlink, "List", 1);
	MLPutSymbol(stdlink,"$Failed");
  	MLPutFunction(stdlink, "List", 1);
	MLPutSymbol(stdlink,"$Failed");
  	MLPutFunction(stdlink, "List", 1);
	MLPutSymbol(stdlink,"$Failed");
  	MLPutFunction(stdlink, "List", 1);
	MLPutSymbol(stdlink,"$Failed");
}

:Begin:
:Function:       MieMma
:Pattern:        MieMma[x_Real, mreal_Real, mimag_Real, mu_List]
:Arguments:      { x, mreal, mimag, mu}
:ArgumentTypes:  { Real, Real, Real, RealList }
:ReturnType:     Manual
:End:
:Evaluate:       Mie[x_, mreal_, mimag_, mu_] := MieMma[N[x], N[mreal], N[mimag], Flatten[{N[mu]}]]
:Evaluate:       Mie::badDiameter = "Size parameter x=pi*d/lambda must greater than zero but less than 10000."
:Evaluate:       Mie::badRealIndex = "Index of refraction must be positive and less than 10."
:Evaluate:       Mie::badImagIndex = "Imaginary part of index of refraction must be non-positive."
:Evaluate:       Mie::badAngle = "List of angles must be a list of cosines of the angles."

void MieMma(double x, double mreal, double mimag, double *mu, long mulen)
{
  double qsca, g, qext, qback;
  double *s1real, *s1imag, *s2real, *s2imag;
  double rr[4];
  long i;
  
  /* check the input variables */
  if (x<=0       || x>10000   ) {error("Message[Mie::badDiameter]"); return;}  
  if (mreal <= 0 || mreal > 10) {error("Message[Mie::badRealIndex]"); return;} 
  if (mimag > 0               ) {error("Message[Mie::badImagIndex]"); return;} 
  for (i=0; i<mulen; i++) {
  	if ((mu[i] < -1) || (mu[i] > 1)) {
  		error("Message[Mie::badAngle]"); 
  		return;
  	}
  }
  
  qsca = 0;
  g = 0;
  qext = 0;
  qback = 0;
  
  s1real=new_darray(mulen);
  s1imag=new_darray(mulen);
  s2real=new_darray(mulen);
  s2imag=new_darray(mulen);
  
  ez_Mie_Full(x, mreal, mimag, mulen, mu, s1real, s1imag, s2real, s2imag, &qext, &qsca, &qback, &g);

  rr[0] = qsca;
  rr[1] = g;
  rr[2] = qext;
  rr[3] = qback;
  
  MLPutFunction(stdlink, "List", 5);
  MLPutRealList(stdlink, rr, 4);
  MLPutRealList(stdlink, s1real, mulen);
  MLPutRealList(stdlink, s1imag, mulen);
  MLPutRealList(stdlink, s2real, mulen);
  MLPutRealList(stdlink, s2imag, mulen);
  
  free_darray(s1real);
  free_darray(s1imag);
  free_darray(s2real);
  free_darray(s2imag);
}

int main(int argc, char* argv[])
{
	return MLMain(argc, argv);
}


