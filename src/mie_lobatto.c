

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "mie_array.h"
#include "mie_legendre.h"
#include "mie_lobatto.h"
#include "mie_root.h"

#define NSLICES 1000
#define EPS 1e-16 \



static long	    Lobatto_n_minus_1;



static void
Lobatto_error(char *s)
{


    printf("%s\n", s);
    exit(1);
}




static double
Lobatto_fn1(double x)
{


    return LegendrePnd(Lobatto_n_minus_1, x);
}




static void
Lobatto_fn2(double x, double *f, double *df)
{


    *f = LegendrePnd(Lobatto_n_minus_1, x);
    *df = LegendrePndd(Lobatto_n_minus_1, x);
}





void
Lobatto(double a, double b, double *x, double *w, long n)
{


    long		nby2, n_odd, i;
    double		xm, xl, pnval;

    if (n < 3)
	Lobatto_error("Number of Lobatto quadrature points less than 3");

    if (x == NULL)
	Lobatto_error("NULL value passed for x array to Lobatto");

    if (w == NULL)
	Lobatto_error("NULL value passed for w array to Lobatto");

    x[n - 1] = 1.0;
    w[n - 1] = 2.0 / n / (n - 1);
    nby2 = n / 2 - 1;
    n_odd = n % 2;

    switch (n) {
    case 4:
	x[2] = 0.4472135954999579;

	w[2] = 0.8333333333333333;
	break;


    case 8:
	x[6] = 0.8717401485096066;
	x[5] = 0.5917001814331423;
	x[4] = 0.2092992179024789;

	w[6] = 0.2107042271435061;
	w[5] = 0.3411226924835043;
	w[4] = 0.4124587946587038;
	break;


    case 16:
	x[14] = 0.9695680462702180;
	x[13] = 0.8992005330934720;
	x[12] = 0.7920082918618151;
	x[11] = 0.6523887028824931;
	x[10] = 0.4860594218871376;
	x[9] = 0.2998304689007632;
	x[8] = 0.1013262735219495;

	w[14] = 0.0508503610059200;
	w[13] = 0.0893936973259308;
	w[12] = 0.1242553821325141;
	w[11] = 0.1540269808071643;
	w[10] = 0.1774919133917041;
	w[9] = 0.1936900238252036;
	w[8] = 0.2019583081782299;
	break;


    default:
	{
	    long		nb, ndiv, size;
	    double		z, *xb1, *xb2;

	    Lobatto_n_minus_1 = n - 1;
	    size = NSLICES;
	    xb1 = new_darray(size);
	    xb2 = new_darray(size);


	    ndiv = nby2;
	    do {
		ndiv *= 2;
		if (ndiv >= NSLICES)
		    ndiv = NSLICES - 1;
		nb = nby2;
		bracketroot(Lobatto_fn1, 0.0, 1.0, ndiv, xb1, xb2, &nb);
	    }
	    while (nb < nby2 && ndiv < NSLICES - 1);

	    if (nb < nby2)
		Lobatto_error("Cannot find enough roots for Lobatto quadrature");



	    for (i = 0; i < nby2; i++) {
		z = saferoot(Lobatto_fn2, xb1[i], xb2[i], EPS);
		x[n - nby2 + i - 1] = z;
		pnval = LegendrePn(n - 1, z);
		w[n - nby2 + i - 1] = w[n - 1] / (pnval * pnval);
	    }


	    free_darray(xb1);
	    free_darray(xb2);
	    break;
	}


    }

    if (n_odd) {
	i = nby2 + 1;
	x[i] = 0.0;
	pnval = LegendrePn(n - 1, 0.0);
	w[i] = 2 / (n * (n - 1) * pnval * pnval);
    }




    for (i = 0; i <= nby2; i++) {
	w[i] = w[n - i - 1];
	x[i] = -x[n - i - 1];
    }



    if ((a != -1.0) | (b != 1.0)) {
	xm = (b + a) / 2.0;
	xl = (b - a) / 2.0;

	for (i = 0; i < n; i++) {
	    x[i] = xm - xl * x[i];
	    w[i] = xl * w[i];
	}
    }


}
