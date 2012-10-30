
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "mie_complex.h"



static void 
complex_error(char *s)
{
    printf("%s\n", s);
    exit(1);
}





struct c_complex 
c_set(double a, double b)
{
    struct c_complex    c;

    c.re = a;
    c.im = b;
    return c;
}




struct c_complex 
c_polarset(double r, double theta)
{
    return c_set(r * cos(theta), r * sin(theta));
}




double 
c_abs(struct c_complex z)
{
    double              x, y, temp;

    x = fabs(z.re);
    y = fabs(z.im);
    if (x == 0.0)
	return y;
    if (y == 0.0)
	return x;

    if (x > y) {
	temp = y / x;
	return (x * sqrt(1.0 + temp * temp));
    }
    temp = x / y;
    return (y * sqrt(1.0 + temp * temp));
}




double 
c_arg(struct c_complex z)
{
    return atan2(z.im, z.re);
}




double 
c_norm(struct c_complex z)
{
    return (z.re * z.re + z.im * z.im);
}




struct c_complex 
c_sqrt(struct c_complex z)
{
    double              a, b;

    if ((z.re == 0.0) && (z.im == 0.0))
	return c_set(0.0, 0.0);

    a = sqrt((fabs(z.re) + c_abs(z)) * 0.5);
    if (z.re >= 0)
	b = z.im / (a + a);
    else {
	b = z.im < 0 ? -a : a;
	a = z.im / (b + b);
    }

    return c_set(a, b);
}




struct c_complex 
c_sqr(struct c_complex z)
{
    return c_mul(z, z);
}




struct c_complex 
c_inv(struct c_complex w)
{
    double              r, d;

    if ((w.re == 0) && (w.im == 0))
	complex_error("Attempt to invert 0+0i");

    if (fabs(w.re) >= fabs(w.im)) {
	r = w.im / w.re;
	d = 1 / (w.re + r * w.im);
	return c_set(d, -r * d);
    }
    r = w.re / w.im;
    d = 1 / (w.im + r * w.re);
    return c_set(r * d, -d);
}




struct c_complex 
c_conj(struct c_complex z)
{
    return c_set(z.re, -z.im);
}




struct c_complex 
c_add(struct c_complex z, struct c_complex w)
{
    struct c_complex    c;

    c.im = z.im + w.im;
    c.re = z.re + w.re;
    return c;
}





struct c_complex 
c_sub(struct c_complex z, struct c_complex w)
{
    struct c_complex    c;

    c.im = z.im - w.im;
    c.re = z.re - w.re;
    return c;
}




struct c_complex 
c_mul(struct c_complex z, struct c_complex w)
{
    struct c_complex    c;

    c.re = z.re * w.re - z.im * w.im;
    c.im = z.im * w.re + z.re * w.im;
    return c;
}




struct c_complex 
c_div(struct c_complex z, struct c_complex w)
{
    struct c_complex    c;
    double              r, denom;

    if ((w.re == 0) && (w.im == 0))
	complex_error("Attempt to divide by 0+0i");

    if (fabs(w.re) >= fabs(w.im)) {
	r = w.im / w.re;
	denom = w.re + r * w.im;
	c.re = (z.re + r * z.im) / denom;
	c.im = (z.im - r * z.re) / denom;
    } else {
	r = w.re / w.im;
	denom = w.im + r * w.re;
	c.re = (z.re * r + z.im) / denom;
	c.im = (z.im * r - z.re) / denom;
    }
    return c;
}




double 
c_rdiv(struct c_complex z, struct c_complex w)
{
    double              r, c, denom;

    if ((w.re == 0) && (w.im == 0))
	complex_error("Attempt to find real part with divisor 0+0i");

    if (fabs(w.re) >= fabs(w.im)) {
	r = w.im / w.re;
	denom = w.re + r * w.im;
	c = (z.re + r * z.im) / denom;
    } else {
	r = w.re / w.im;
	denom = w.im + r * w.re;
	c = (z.re * r + z.im) / denom;
    }
    return c;
}




double 
c_rmul(struct c_complex z, struct c_complex w)
{
    return z.re * w.re - z.im * w.im;
}





struct c_complex 
c_sadd(double x, struct c_complex z)
{
    struct c_complex    c;

    c.re = x + z.re;
    c.im = z.im;
    return c;
}




struct c_complex 
c_sdiv(double x, struct c_complex w)
{
    struct c_complex    c;
    double              r, factor;

    if ((w.re == 0) && (w.im == 0))
	complex_error("Attempt to divide scalar by 0+0i");

    if (fabs(w.re) >= fabs(w.im)) {
	r = w.im / w.re;
	factor = x / (w.re + r * w.im);
	c.re = factor;
	c.im = -r * factor;
    } else {
	r = w.re / w.im;
	factor = x / (w.im + r * w.re);
	c.im = -factor;
	c.re = r * factor;
    }
    return c;
}




struct c_complex 
c_smul(double x, struct c_complex z)
{
    struct c_complex    c;

    c.re = z.re * x;
    c.im = z.im * x;
    return c;
}





struct c_complex 
c_sin(struct c_complex z)
{
    return c_set(sin(z.re) * cosh(z.im), cos(z.re) * sinh(z.im));
}




struct c_complex 
c_cos(struct c_complex z)
{
    return c_set(cos(z.re) * cosh(z.im), -(sin(z.re) * sinh(z.im)));
}




struct c_complex 
c_tan(struct c_complex z)
{
    double              t, x, y;

    if (z.im == 0)
	return c_set(tan(z.re), 0.0);
    if (z.im > DBL_MAX_10_EXP)
	return c_set(0.0, 1.0);
    if (z.im < -DBL_MAX_10_EXP)
	return c_set(0.0, -1.0);

    x = 2 * z.re;
    y = 2 * z.im;
    t = cos(x) + cosh(y);
    if (t == 0)
	complex_error("Complex tangent is infinite");

    return c_set(sin(x) / t, sinh(y) / t);
}




struct c_complex 
c_asin(struct c_complex z)
{
    struct c_complex    x;

    x = c_log(c_add(c_set(-z.im, z.re), c_sqrt(c_sub(c_set(1.0, 0.0), c_mul(z, z)))));
    return c_set(x.im, -x.re);
}




struct c_complex 
c_acos(struct c_complex z)
{
    struct c_complex    x;

    x = c_log(c_add(z, c_mul(c_set(0.0, 1.0), c_sqrt(c_sub(c_set(1.0, 0.0), c_sqr(z))))));
    return c_set(x.im, -x.re);
}




struct c_complex 
c_atan(struct c_complex z)
{
    struct c_complex    x;

    x = c_log(c_div(c_set(z.re, 1 + z.im), c_set(-z.re, 1 - z.im)));
    return c_set(-x.im / 2, x.re / 2);
}





struct c_complex 
c_sinh(struct c_complex z)
{
    return c_set(sinh(z.re) * cos(z.im), cosh(z.re) * sin(z.im));
}




struct c_complex 
c_cosh(struct c_complex z)
{
    return c_set(cosh(z.re) * cos(z.im), sinh(z.re) * sin(z.im));
}




struct c_complex 
c_tanh(struct c_complex z)
{
    double              x = 2 * z.re;
    double              y = 2 * z.im;
    double              t = 1.0 / (cosh(x) + cos(y));

    return c_set(t * sinh(x), t * sin(y));
}




struct c_complex 
c_atanh(struct c_complex z)
{
    return c_atan(c_set(-z.im, z.re));
}




struct c_complex 
c_asinh(struct c_complex z)
{
    return c_asin(c_set(-z.im, z.re));
}





struct c_complex 
c_exp(struct c_complex z)
{
    double              x = exp(z.re);

    return c_set(x * cos(z.im), x * sin(z.im));
}




struct c_complex 
c_log(struct c_complex z)
{
    return c_set(log(c_abs(z)), c_arg(z));
}




struct c_complex 
c_log10(struct c_complex z)
{
    return c_set(0.2171472409516259 * log(c_norm(z)), c_arg(z));
}





struct c_complex   *
new_carray(long size)
{
    struct c_complex   *a;

    if (size <= 0)
	complex_error("Non-positive complex array size chosen");

    a = (struct c_complex *) calloc(sizeof(struct c_complex), (unsigned long) size);

    if (a == NULL)
	complex_error("Can't allocate complex array");
    return a;
}




void 
free_carray(struct c_complex * a)
{
    if (a != NULL)
	free(a);
}




struct c_complex   *
copy_carray(struct c_complex * a, long size)
{
    struct c_complex   *b = NULL;

    if (a == NULL)
	complex_error("Can't duplicate a NULL complex array");

    b = new_carray(size);
    if (b != NULL)
	memcpy(b, a, size * sizeof(struct c_complex));
    return b;
}




void 
set_carray(struct c_complex * a, long size, struct c_complex z)
{
    long                j;

    if (a == NULL)
	complex_error("Can't operate on a NULL complex array");

    for (j = 0; j < size; j++)
	a[j] = z;
}
