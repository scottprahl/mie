@** Lobatto Quadrature.
These routines are useful in when both ends of the interval must
be calculated for some other reason.  Lobatto quadrature also works
well when both ends of the interval are equal to zero.

This global variable is needed because the degree of the
Legendre Polynomial must be known.  The routine |Lobatto| stores
the correct value in this.  I assume the Numerical Recipes version
of arrays in this program.

The |Lobatto| routine was tested for several values of |n| and
compared with the values from the paper by Michel.

@d NSLICES     1000
@d EPS 		1e-16

@(mie_lobatto.c@>=

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "mie_array.h"
#include "mie_legendre.h"
#include "mie_lobatto.h"
#include "mie_root.h"
	
	@h
    static long Lobatto_n_minus_1;

	@<Definition for |Lobatto_error|@>@;
	@<Definition for |Lobatto_fn1|@>@;
	@<Definition for |Lobatto_fn2|@>@;
	@<Definition for |Lobatto|@>@;
	
@ @(mie_lobatto.h@>=
	@<Prototype for |Lobatto|@>;

@ A simple error routine.
@<Prototype for |Lobatto_error|@>=
static void Lobatto_error(char *s)

@ @<Definition for |Lobatto_error|@>=
		@<Prototype for |Lobatto_error|@>
{
 printf("%s\n", s);
 exit(1);
 }

@*1 Lobatto functions.

These functions rely on the local variable |Lobatto_n_minus_1|.

@ This funtion is used to bracket all the roots.  It just
returns $P_n'(x)$.
@<Prototype for |Lobatto_fn1|@>=
static double Lobatto_fn1(double x)

@ @<Definition for |Lobatto_fn1|@>=
	@<Prototype for |Lobatto_fn1|@>
{
return LegendrePnd(Lobatto_n_minus_1, x);
}

@ This funtion is used to find each root.  Since Newton's
method is used both the function and the first derivative
are needed.  This routine returns both
$f=P_n'(x)$ and $df=P_n''(x)$.

@<Prototype for |Lobatto_fn2|@>=
static void Lobatto_fn2(double x, double *f, double *df)

@ @<Definition for |Lobatto_fn2|@>=
	@<Prototype for |Lobatto_fn2|@>
{
  *f= LegendrePnd(Lobatto_n_minus_1, x);
  *df= LegendrePndd(Lobatto_n_minus_1, x);
}


@*1 Lobatto Tables.

Here is a selection of commonly used number of quadrature points.

@ @<Values for |n==4|@>=
	x[2] = 0.4472135954999579;

	w[2] = 0.8333333333333333;
	break;

@ @<Values for |n==8|@>= 
	x[6] = 0.8717401485096066;
	x[5] = 0.5917001814331423;
	x[4] = 0.2092992179024789;
	
	w[6] = 0.2107042271435061;
	w[5] = 0.3411226924835043;
	w[4] = 0.4124587946587038;
	break;

@ @<Values for |n==16|@>=
	x[14] = 0.9695680462702180;
	x[13] = 0.8992005330934720;
	x[12] = 0.7920082918618151;
	x[11] = 0.6523887028824931;
	x[10] = 0.4860594218871376;
	x[ 9] = 0.2998304689007632;
	x[ 8] = 0.1013262735219495;
	
	w[14] = 0.0508503610059200;
	w[13] = 0.0893936973259308;
	w[12] = 0.1242553821325141;
	w[11] = 0.1540269808071643;
	w[10] = 0.1774919133917041;
	w[ 9] = 0.1936900238252036;
	w[ 8] = 0.2019583081782299;
	break;

@*1 Lobatto.
|Lobatto| calculates the |n| quadrature points $x_i$ and weights $w_i$
over the interval $(a,b)$.  The basis fo this is Lobatto quadrature
$$
\int_a^b f(x) dx = w_0 f(a) + \sum_{k=1}^{n-2} w_k f(x_k) + w_{n-1} f(b)
$$
where the quadrature points $x_k$ $(k=1,2,\cdots,n-2)$ are the zeros of
$$
P'_{n-1}(x_k) = 0
$$
and the weights are given by
$$
w_k = {2\over n(n-1)[P_{n-1}(x_k)]^2}
$$
Finally
$$
w_0 = w_{n-1} = {2\over n(n-1)}
$$

@<Prototype for |Lobatto|@>=
void Lobatto(double a, double b, double *x, double *w, long n)

@ @<Definition for |Lobatto|@>=
	@<Prototype for |Lobatto|@>
{
long nby2, n_odd, i;
double xm, xl, pnval;

	if (n<3) Lobatto_error("Number of Lobatto quadrature points less than 3");

	if (x==NULL) Lobatto_error("NULL value passed for x array to Lobatto");

	if (w==NULL) Lobatto_error("NULL value passed for w array to Lobatto");
	
	x[n-1] =  1.0;
	w[n-1] =  2.0 / n / (n-1);
    nby2 = n / 2 - 1;
	n_odd = n % 2;
	
	switch (n) {
		case  4: @<Values for |n==4|@>@;
		case  8: @<Values for |n==8|@>@;
		case 16: @<Values for |n==16|@>@;
		default: @<Values for arbitrary |n|@>@;
	}

	if (n_odd) @<Do middle value@>@;
	
	@<Do negative values@>@;
	@<Scale values@>@;
}

@ If the number of quadrature points in odd, then
the center value $x_i$ will always be zero.  Furthermore
$w_{n-i-1} = w_k$.

@<Do middle value@>=
	{
		i=nby2+1;
		x[i] = 0.0;
		pnval = LegendrePn(n-1,0.0);
		w[i] = 2/(n*(n-1)*pnval*pnval);
	}

@ The quadrature points are symmetric with $x_{n-i-1}=-x_i$
and $w_{n-i-1} = w_i$.  Just copy the top half of the array
into the bottom half.

@<Do negative values@>=
	for(i=0; i<=nby2; i++) {
		w[i] =  w[n-i-1];
		x[i] = -x[n-i-1];
	}

@ The code to scale values is easy.  Lobatto quadrature is
defined over the range $-1$ to 1.  To modify to the range $a$
to $b$, I just linearly scale
the width of each interval and weight as appropriate.
$$
x_i = {a+b\over2}+{a-b \over 2} x_i
$$
and
$$
w_i=   {b-a \over 2} w_i
$$

@<Scale values@>=
if ((a!=-1.0)|(b!=1.0)) {
	xm = (b + a) / 2.0;
	xl = (b - a) / 2.0;

	for (i = 0; i < n; i++) {
		x[i] = xm - xl * x[i];
		w[i] = xl * w[i];
	}
}

@ Here is the method for finding Lobatto quadrature points for
non-tabulated values.  There will be $n/2$ roots located between
0 and 1.  

The only strange part is that the local variable |Lobatto_n_minus_1|
is set here so that it will be set correctly when |Lobatto_fn1| and
|Lobatto_fn2| get called.

@<Values for arbitrary |n|@>=
{
	long nb, ndiv, size;
	double z, *xb1, *xb2;

	Lobatto_n_minus_1 = n-1;
	size = NSLICES;
	xb1 = new_darray(size);
	xb2 = new_darray(size);
	@<Bracket roots@>@;
	@<Find roots and weights@>@;
	free_darray(xb1);
	free_darray(xb2);
	break;
	}

@ Bracket $n/2$ roots, double |ndiv| if not enough roots are found.
I make every effort to find all the damn roots.

@<Bracket roots@>=
	
	ndiv = nby2;
	do{
		ndiv *= 2;
		if (ndiv>=NSLICES) ndiv = NSLICES-1;
		nb = nby2;
		bracketroot(Lobatto_fn1, 0.0, 1.0, ndiv, xb1, xb2, &nb);
	}
	while (nb < nby2 && ndiv < NSLICES-1);

	if (nb < nby2) 
		Lobatto_error("Cannot find enough roots for Lobatto quadrature");

@ Find the roots with an accuracy |EPS| and store them in the array |x|.
Put them in backwards so that |x[n-1]=-1| is in the correct spot.

@<Find roots and weights@>=
	for (i = 0; i < nby2; i++) {
		z = saferoot(Lobatto_fn2, xb1[i], xb2[i], EPS);
		x[n-nby2+i-1] = z;
		pnval = LegendrePn(n-1,z);
		w[n-nby2+i-1] = w[n-1] / (pnval*pnval);
	}

@*1 Testing Lobatto.

@(test_mie_lobatto.c@>=

#include <stddef.h>
#include <stdio.h>
#include "mie_array.h"
#include "mie_lobatto.h"
	
int main() 
{
	double *x, *w;
	double sum;
	long i, n;
	
	printf("testing n=10 --- sum should be 2\n");
	sum=0;
	n=10;
	x=new_darray(n);
	w=new_darray(n);	
	Lobatto(-1.0,1.0,x,w,n);
	printf("The x_i are\n");
	print_darray(x,n,0L,n-1);
	printf("The w_i are\n");
	print_darray(w,n,0L,n-1);
	for (i=0; i<n; i++) sum = sum+w[i];
	printf("sum   %20.15f\n",sum);
	free_darray(x);
	free_darray(w);
	printf("\n");
	
	printf("testing n=9 --- sum should be 2\n");
	sum=0;
	n=9;
	x=new_darray(n);
	w=new_darray(n);	
	Lobatto(-1.0,1.0,x,w,n);
	printf("The x_i are\n");
	print_darray(x,n,0L,n-1);
	printf("The w_i are\n");
	print_darray(w,n,0L,n-1);
	for (i=0; i<n; i++) sum = sum+w[i];
	printf("sum   %20.15f\n",sum);
	free_darray(x);
	free_darray(w);
	printf("\n");
	return 0;
}

	
