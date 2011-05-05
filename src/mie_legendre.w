@** Legendre Polynomials.

To calculate the values of the Legendre polynomial using
the recurrence relations given by  H. H. Michels, 
in ``Abscissas and weight coefficients for
Lobatto quadrature,'' {\it Math Comp\/}, {\bf 17}, 237-244 (1963).

These were checked for $n=10$ for $-1<=x<=1$ using the basic
differential equation
$$
(1-x^2) P_n''(x) - 2x P_n'(x)+ n(n+1) P_n(x) = 0
$$

$P_n(x)$ and $P_n'(x)$ were also checked against Abromowitz and Stegun.

@(mie_legendre.c@>=
	#include "mie_legendre.h"

	@<Definition for |LegendrePn|@>@;
	@<Definition for |LegendrePn_and_Pnm1|@>@;
	@<Definition for |LegendrePnd|@>@;
	@<Definition for |LegendrePndd|@>@;

@ @(mie_legendre.h@>=
	@<Prototype for |LegendrePn|@>;
	@<Prototype for |LegendrePn_and_Pnm1|@>;
	@<Prototype for |LegendrePnd|@>;
	@<Prototype for |LegendrePndd|@>;

@*1 Basic Legendre functions.

Returns the Legendre polynomial $P_n(x)$
using the recurrence formulas from Michel's paper,
$$
P_{n+1}(x) = \left({2n+1\over n+1}\right) x P_n(x) -
             \left({n   \over n+1}\right)   P_{n-1}(x)
$$
and
$$
P_0(x) = 1 \qquad\hbox{and}\qquad P_1(x) = x
$$
Also added special cases for $x=\pm1$.  These are $P_n(1)=1$ and
$P_n(-1)= (-1)^n$ and for what it is worth 
$$
P_n(0) = {1\cdot3\cdot5\cdots(2n-1)\over 2\cdot4\cdot6\cdots2n}
$$

If the argument is out of range ($\vert x\vert >1$) then instead
of complaining, I just return the value at $x=\pm1$ as appropriate.
If $n<0$ then I just return 1.0 and don't complain either.

Finally, I decided not to 
@<Prototype for |LegendrePn|@>=
double LegendrePn(long n, double x)

@ @<Definition for |LegendrePn|@>=
	@<Prototype for |LegendrePn|@>
{
	double pk, pkp1, pkm1;
	long k;

	if (n<=0) return 1.0;
	if (n==1) return x;
	
	if (x >= 1.0) return 1.0;	
	if (x <=-1.0) return (n % 2) ? -1.0 : 1.0;

	pk =  x;
	pkm1=1.0;
	for (k=1; k<n; k++) {
		pkp1 = ((2*k+1) * x * pk - k * pkm1)/(k+1) ;
		pkm1 = pk;
		pk = pkp1;
	}
	
	return pk;
}

@ |LegendrePn_and_Pnm1| returns $P_n(x)$ and $P_{n-1}(x)$

@<Prototype for |LegendrePn_and_Pnm1|@>=
void LegendrePn_and_Pnm1(long n, double x, double *Pnm1Val, double *PnVal)

@ @<Definition for |LegendrePn_and_Pnm1|@>=
	@<Prototype for |LegendrePn_and_Pnm1|@>
{
	long k;
	double Pk, Pkp1;
	double Pkm1=1.0;

	*Pnm1Val=1.0;
	*PnVal=1.0;
	if (x >= 1.0) x= 1.0;
	if (x <=-1.0) x=-1.0;

	Pk =  x;
	
	for (k=1; k<n; k++) {
		Pkp1 = ((2*k+1) * x * Pk - k * Pkm1) / (k+1);
		Pkm1 = Pk;
		Pk = Pkp1;
	}
	
	*Pnm1Val=Pkm1;
	*PnVal=Pk;
}

@*1 First derivative.

Returns the first derivative of the Legendre polynomial $P_n'(x)$
using the recurrence formulas from Michel's paper,
$$
P_{n+1}'(x) = \left({2n+1\over n}\right) x P_n'(x) -
              \left({n+1 \over n}\right) P_{n-1}'(x)
$$
and
$$
P_0'(x) = 0 \qquad\hbox{and}\qquad P_1'(x) = 1
$$

@<Prototype for |LegendrePnd|@>=
double LegendrePnd(long n, double x)

@ @<Definition for |LegendrePnd|@>=
	@<Prototype for |LegendrePnd|@>
{
	double p, pminus, pplus;
	long i;

	if (n <= 0) return 0;
	if (n == 1) return 1;

	if (x > 1.0) x= 1.0;
	if (x <-1.0) x=-1.0;

	pminus = 0;
	p =  1;

	for (i=1; i<n; i++) {
		pplus = ((2*i + 1) * x * p - (i + 1) * pminus)/i;
		pminus = p;
		p = pplus;
	}
	return p;
}


@*1 Second derivative.

Returns the second derivative of the Legendre polynomial $P_n''(x)$
using the recurrence formulas
$$
P_{n+1}''(x) = \left({2n+1\over n-1}\right) x P_n''(x) -
               \left({n+2 \over n-1}\right) P_{n-1}''(x)
$$
and
$$
P_1''(x) = 0 \qquad\hbox{and}\qquad P_2''(x) = 3
$$

@<Prototype for |LegendrePndd|@>=
double LegendrePndd(long n, double x)

@ @<Definition for |LegendrePndd|@>=
	@<Prototype for |LegendrePndd|@>
{
  double p, pminus, pplus;
  long m;

  if (n <= 1) return 0;
  if (n == 2) return 3;

  if (x > 1.0) x= 1.0;
  if (x <-1.0) x=-1.0;

  pminus = 0;
  p = 3;

  for (m=2; m < n; m++) {
      pplus = ((2*m+1) * x * p - (m + 2) * pminus)/ (m - 1);
      pminus = p;
      p = pplus;
    }
  return p;

}
