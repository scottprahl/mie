@** Complex Number Routines.

Here are a bunch of routines to deal with complex numbers.  The functions
are pretty straightforward, but there are some subtle points in some of
the functions.  This could use some more error checking.

Changed names to not conflict with c++ routines

@ Here, then, is an overview of document structure
@(mie_complex.c@>=
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "mie_complex.h"

@<Definition for |complex_error|@>@;

@<Definition for |c_set|@>@;
@<Definition for |c_polarset|@>@;
@<Definition for |c_abs|@>@;
@<Definition for |c_arg|@>@;
@<Definition for |c_norm|@>@;
@<Definition for |c_sqrt|@>@;
@<Definition for |c_sqr|@>@;
@<Definition for |c_inv|@>@;
@<Definition for |c_conj|@>@;
@<Definition for |c_add|@>@;

@<Definition for |c_sub|@>@;
@<Definition for |c_mul|@>@;
@<Definition for |c_div|@>@;
@<Definition for |c_rdiv|@>@;
@<Definition for |c_rmul|@>@;

@<Definition for |c_sadd|@>@;
@<Definition for |c_sdiv|@>@;
@<Definition for |c_smul|@>@;

@<Definition for |c_sin|@>@;
@<Definition for |c_cos|@>@;
@<Definition for |c_tan|@>@;
@<Definition for |c_asin|@>@;
@<Definition for |c_acos|@>@;
@<Definition for |c_atan|@>@;

@<Definition for |c_sinh|@>@;
@<Definition for |c_cosh|@>@;
@<Definition for |c_tanh|@>@;
@<Definition for |c_atanh|@>@;
@<Definition for |c_asinh|@>@;

@<Definition for |c_exp|@>@;
@<Definition for |c_log|@>@;
@<Definition for |c_log10|@>@;

@<Definition for |new_carray|@>@;
@<Definition for |free_carray|@>@;
@<Definition for |copy_carray|@>@;
@<Definition for |set_carray|@>@;

@ Each function has its prototype exported to a header file along
with a couple of structure definitions.

@(mie_complex.h@>=
struct c_complex {double re; double im;};

@<Prototype for |c_set|@>;
@<Prototype for |c_polarset|@>;

@<Prototype for |c_abs|@>;
@<Prototype for |c_arg|@>;
@<Prototype for |c_sqr|@>;
@<Prototype for |c_conj|@>;
@<Prototype for |c_norm|@>;
@<Prototype for |c_sqrt|@>;
@<Prototype for |c_inv|@>;

@<Prototype for |c_add|@>;
@<Prototype for |c_sub|@>;
@<Prototype for |c_mul|@>;
@<Prototype for |c_div|@>;

@<Prototype for |c_rdiv|@>;
@<Prototype for |c_rmul|@>;
@<Prototype for |c_sadd|@>;
@<Prototype for |c_sdiv|@>;
@<Prototype for |c_smul|@>;

@<Prototype for |c_sin|@>;
@<Prototype for |c_cos|@>;
@<Prototype for |c_tan|@>;
@<Prototype for |c_asin|@>;
@<Prototype for |c_acos|@>;
@<Prototype for |c_atan|@>;

@<Prototype for |c_sinh|@>;
@<Prototype for |c_cosh|@>;
@<Prototype for |c_tanh|@>;
@<Prototype for |c_atanh|@>;
@<Prototype for |c_asinh|@>;

@<Prototype for |c_exp|@>;
@<Prototype for |c_log|@>;
@<Prototype for |c_log10|@>;

@<Prototype for |new_carray|@>;
@<Prototype for |free_carray|@>;
@<Prototype for |copy_carray|@>;
@<Prototype for |set_carray|@>;

@*1 Basic routines.

@ A simple error routine.
@<Prototype for |complex_error|@>=
static void complex_error(char *s)

@ @<Definition for |complex_error|@>=
		@<Prototype for |complex_error|@>
{
 printf("%s\n", s);
 exit(1);
 }

@ This is shorthand for setting a complex number.  It
just returns a complex equal to $a+b i$

@<Prototype for |c_set|@>=
struct c_complex c_set(double a, double b)  

@ @<Definition for |c_set|@>=
		@<Prototype for |c_set|@>
{
 struct c_complex c;
  c.re = a;
  c.im = b;
  return c;
}

@ A variation on |c_set| in which the complex
number is specifed using polar coordinates.

@<Prototype for |c_polarset|@>=
struct c_complex c_polarset (double r, double theta)

@ @<Definition for |c_polarset|@>=
		@<Prototype for |c_polarset|@>
{
  return c_set(r * cos(theta), r * sin(theta));
}

@ This routine returns the absolute value of a complex number
$\sqrt{zz^*}$.  To avoid unnecessary loss of accuracy as explained
in \S 5.4 of {\it Numerical Recipes in C}

@<Prototype for |c_abs|@>=
double c_abs(struct c_complex z)

@ @<Definition for |c_abs|@>=
		@<Prototype for |c_abs|@>
{
double x,y,temp;

x=fabs(z.re);
y=fabs(z.im);
if (x == 0.0)  return y;
if (y == 0.0)  return x;

if (x > y) {
	temp=y/x;
	return (x*sqrt(1.0+temp*temp));
} 
	
temp=x/y;
return (y*sqrt(1.0+temp*temp));
}

@ Returns the conjugate of the complex number $z$

@<Prototype for |c_conj|@>=
struct c_complex c_conj(struct c_complex z)

@ @<Definition for |c_conj|@>=
		@<Prototype for |c_conj|@>
{
return c_set(z.re,-z.im);
}

@ @<Prototype for |c_arg|@>=
double c_arg(struct c_complex z)

@ @<Definition for |c_arg|@>=
		@<Prototype for |c_arg|@>
{
  return atan2(z.im,z.re);
}

@ Returns the square of the modulus of the complex number $zz^*$

@<Prototype for |c_norm|@>=
double c_norm(struct c_complex z)

@ @<Definition for |c_norm|@>=
		@<Prototype for |c_norm|@>
{
  return (z.re*z.re + z.im*z.im);
}

@ @<Prototype for |c_sqrt|@>=
struct c_complex c_sqrt(struct c_complex z) 

@ @<Definition for |c_sqrt|@>=
		@<Prototype for |c_sqrt|@>
{
  double a, b;

  if ((z.re == 0.0) && (z.im == 0.0) ) 
  return c_set(0.0,0.0);

  a = sqrt ((fabs(z.re) + c_abs(z))*0.5);
  if (z.re >= 0 ) 
  	b = z.im / (a+a);
  else {
    b = z.im < 0 ? -a: a;
    a = z.im / (b+b);
  }

  return c_set(a,b);
}

@ Returns the product of a complex number with itself $z\cdot z$.
If you want $z\cdot z^*$ then use |c_norm|.

@<Prototype for |c_sqr|@>=
struct c_complex c_sqr(struct c_complex z) 

@ @<Definition for |c_sqr|@>=
		@<Prototype for |c_sqr|@>
{
  return c_mul(z,z);
}

@ Returns the reciprocal of |z|.

@<Prototype for |c_inv|@>=
struct c_complex c_inv(struct c_complex w)

@ @<Definition for |c_inv|@>=
		@<Prototype for |c_inv|@>
{
double r,d;
	
if ((w.re==0) && (w.im==0)) complex_error("Attempt to invert 0+0i");

if (fabs(w.re) >= fabs(w.im)) {
	r=w.im/w.re;
	d=1/(w.re+r*w.im);
	return c_set(d, -r*d);
}

r=w.re/w.im;
d=1/(w.im+r*w.re);
return c_set(r*d, -d);
}

@*1 Two complex numbers.

@ Returns the sum of the two complex numbers |a| and |b|

@<Prototype for |c_add|@>=
struct c_complex c_add(struct c_complex z, struct c_complex w)

@ @<Definition for |c_add|@>=
		@<Prototype for |c_add|@>
{
struct c_complex c;

c.im = z.im+w.im;
c.re = z.re+w.re;
return c;
}

@ Returns the difference of two complex numbers |z-w|

@<Prototype for |c_sub|@>=
struct c_complex c_sub(struct c_complex z, struct c_complex w)

@ @<Definition for |c_sub|@>=
		@<Prototype for |c_sub|@>
{
struct c_complex c;
c.im = z.im-w.im;
c.re = z.re-w.re;
return c;
}

@ Returns the product of two complex numbers $z\cdot w$

@<Prototype for |c_mul|@>=
struct c_complex c_mul(struct c_complex z, struct c_complex w)

@ @<Definition for |c_mul|@>=
		@<Prototype for |c_mul|@>
{
struct c_complex c;
c.re = z.re*w.re-z.im*w.im;
c.im = z.im*w.re+z.re*w.im;
return c;
}

@ Returns the quotient of two complex numbers |z/w|  
     see \S 5.4 of {\it Numerical Recipes in C}
 
@<Prototype for |c_div|@>=
struct c_complex c_div(struct c_complex z, struct c_complex w)

@ @<Definition for |c_div|@>=
		@<Prototype for |c_div|@>
{
struct c_complex c;
double r,denom;

if ((w.re==0) && (w.im==0)) complex_error("Attempt to divide by 0+0i"); 

if (fabs(w.re) >= fabs(w.im)) {
	r=w.im/w.re;
	denom=w.re+r*w.im;
	c.re=(z.re+r*z.im)/denom;
	c.im=(z.im-r*z.re)/denom;
} else {
	r=w.re/w.im;
	denom=w.im+r*w.re;
	c.re=(z.re*r+z.im)/denom;
	c.im=(z.im*r-z.re)/denom;
	}
return c;
}

@ Returns the real part of the quotient of two complex numbers $\hbox{Re}(z/w)$.
  Note how this is a special case of |c_div| above

@<Prototype for |c_rdiv|@>=
double c_rdiv(struct c_complex z, struct c_complex w)

@ @<Definition for |c_rdiv|@>=
		@<Prototype for |c_rdiv|@>
{
double r,c, denom;
	
if ((w.re==0) && (w.im==0)) complex_error("Attempt to find real part with divisor 0+0i"); 

if (fabs(w.re) >= fabs(w.im)) {
	r=w.im/w.re;
	denom=w.re+r*w.im;
	c=(z.re+r*z.im)/denom;
} else {
	r=w.re/w.im;
	denom=w.im+r*w.re;
	c=(z.re*r+z.im)/denom;
} 
return c;
}

@ Returns the real part of the product of two complex numbers $\hbox{Re}(z\cdot w)$

@<Prototype for |c_rmul|@>=
double c_rmul(struct c_complex z, struct c_complex w)

@ @<Definition for |c_rmul|@>=
		@<Prototype for |c_rmul|@>
{
  return z.re*w.re - z.im*w.im;
}

@*1 A scalar and a complex number.

@ Returns the product of a scalar with a complex number

@<Prototype for |c_smul|@>=
struct c_complex c_smul(double x, struct c_complex z)

@ @<Definition for |c_smul|@>=
		@<Prototype for |c_smul|@>
{
	struct c_complex c;
	c.re = z.re*x;
	c.im = z.im*x;
	return c;
}

@ Returns the sum of a scalar and a complex number

@<Prototype for |c_sadd|@>=
struct c_complex c_sadd(double x, struct c_complex z)

@ @<Definition for |c_sadd|@>=
		@<Prototype for |c_sadd|@>
{
	struct c_complex c;
	c.re = x+z.re;
	c.im = z.im;
	return c;
}

@ Returns the quotient of real number by a complex number |z|.
  Again a special case of |c_div|

@<Prototype for |c_sdiv|@>=
struct c_complex c_sdiv(double x, struct c_complex w)

@ @<Definition for |c_sdiv|@>=
		@<Prototype for |c_sdiv|@>
{
struct c_complex c;
double r,factor;
	
if ((w.re==0) && (w.im==0)) complex_error("Attempt to divide scalar by 0+0i"); 

if (fabs(w.re) >= fabs(w.im)) {
	r=w.im/w.re;
	factor=x/(w.re+r*w.im);
	c.re=factor;
	c.im=-r*factor;
} else {
	r=w.re/w.im;
	factor=x/(w.im+r*w.re);
	c.im=-factor;
	c.re=r*factor;
}
return c;
}

@*1 Trigonometric Functions.

@ The complex sine.

@<Prototype for |c_sin|@>=
struct c_complex c_sin(struct c_complex z)

@ @<Definition for |c_sin|@>=
		@<Prototype for |c_sin|@>
{
  return c_set (sin(z.re)*cosh(z.im), cos(z.re)*sinh(z.im));
}

@ The complex cosine.

@<Prototype for |c_cos|@>=
struct c_complex c_cos(struct c_complex z)

@ @<Definition for |c_cos|@>=
		@<Prototype for |c_cos|@>
{
  return c_set(cos(z.re) * cosh(z.im), -(sin(z.re) * sinh(z.im)));
}

@ The complex tangent.
$$
\tan (a+bi) = {\sin 2a + i \sinh 2b\over \cos 2a + \cosh 2b}
$$
or
$$
\tan (a+bi) = {2\sin 2a + i \exp(2b) - i \exp(-2b) \over 
               2\cos 2a + \exp(2b) + \exp(-2b)}
$$
it is easy to see that if $2b$ is large, then problems arise.

The number |DBL_MAX_10_EXP| is the value $c$ such that $10^c$ can
be represented by a double precision variable.  Now we are interested
in the maximum exponential, one would just multiply $c$ by $\ln 10 = 2.3$
to get such an exponential.  This could then be compared against the
value of $2b$ to figure out when an approximation should be used.  Slightly
more conservatively, one could just test to see when 
$$
  2b > 2 c
$$
and adjust accordingly.

@<Prototype for |c_tan|@>=
struct c_complex c_tan(struct c_complex z)

@ @<Definition for |c_tan|@>=
		@<Prototype for |c_tan|@>
{
  double t, x, y;

  if (z.im==0)  return c_set(tan(z.re), 0.0);
  if (z.im>DBL_MAX_10_EXP) return c_set(0.0, 1.0);
  if (z.im<-DBL_MAX_10_EXP) return c_set(0.0, -1.0);
  
  x = 2*z.re;
  y = 2*z.im;
  t = cos(x) +cosh(y);
  if (t==0) complex_error("Complex tangent is infinite");
  
  return c_set(sin(x)/t, sinh(y)/t );
}

@ The complex inverse sine.

@<Prototype for |c_asin|@>=
struct c_complex c_asin(struct c_complex z)    

@ @<Definition for |c_asin|@>=
		@<Prototype for |c_asin|@>
{
  struct c_complex x;
  x = c_log(c_add(c_set(-z.im,z.re),c_sqrt(c_sub(c_set(1.0,0.0),c_mul(z,z)))));
  return c_set ( x.im, -x.re );
}   

@ The complex inverse cosine

@<Prototype for |c_acos|@>=
struct c_complex c_acos(struct c_complex z)   

@ @<Definition for |c_acos|@>=
		@<Prototype for |c_acos|@>
{
  struct c_complex x;
  x = c_log (c_add(z, c_mul(c_set (0.0,1.0), c_sqrt(c_sub(c_set(1.0,0.0), c_sqr(z))))));
  return c_set ( x.im, -x.re );
}   

@ The complex inverse tangent

@<Prototype for |c_atan|@>=
struct c_complex c_atan(struct c_complex z)   

@ @<Definition for |c_atan|@>=
		@<Prototype for |c_atan|@>
{
  struct c_complex x;
  x = c_log(c_div(c_set(z.re,1+z.im),c_set(-z.re,1-z.im)));
  return c_set ( -x.im/2, x.re/2 );
}

@*1 Hyperbolic functions.

@ @<Prototype for |c_cosh|@>=
struct c_complex c_cosh(struct c_complex z )

@ @<Definition for |c_cosh|@>=
		@<Prototype for |c_cosh|@>
{
  return c_set (cosh(z.re)*cos(z.im), sinh(z.re)*sin(z.im));
}

@ @<Prototype for |c_sinh|@>=
struct c_complex c_sinh(struct c_complex z )

@ @<Definition for |c_sinh|@>=
		@<Prototype for |c_sinh|@>
{
  return c_set ( sinh(z.re) * cos(z.im), cosh(z.re) * sin(z.im) );
}

@ @<Prototype for |c_tanh|@>=
struct c_complex c_tanh(struct c_complex z )

@ @<Definition for |c_tanh|@>=
		@<Prototype for |c_tanh|@>
{
  double x = 2*z.re;
  double y = 2*z.im;
  double t = 1.0/(cosh(x) +cos(y));
  
  return c_set( t*sinh(x), t*sin(y) );
}    

@ @<Prototype for |c_atanh|@>=
struct c_complex c_atanh(struct c_complex z) 

@ @<Definition for |c_atanh|@>=
		@<Prototype for |c_atanh|@>
{
  return c_atan(c_set(-z.im, z.re));
}

@ @<Prototype for |c_asinh|@>=
struct c_complex c_asinh(struct c_complex z)

@ @<Definition for |c_asinh|@>=
		@<Prototype for |c_asinh|@>
{
  return c_asin(c_set(-z.im, z.re));
}

@*1 Exponentials and logarithms.

@ @<Prototype for |c_exp|@>=
struct c_complex c_exp(struct c_complex z )

@ @<Definition for |c_exp|@>=
		@<Prototype for |c_exp|@>
{
  double x = exp(z.re);
  return c_set( x*cos(z.im), x*sin(z.im) );
}    

@ @<Prototype for |c_log|@>=
struct c_complex c_log(struct c_complex z )

@ @<Definition for |c_log|@>=
		@<Prototype for |c_log|@>
{
  return c_set( log( c_abs(z) ), c_arg( z ) );
}    

@ @<Prototype for |c_log10|@>=
struct c_complex c_log10(struct c_complex z )

@ @<Definition for |c_log10|@>=
		@<Prototype for |c_log10|@>
{
  return c_set( 0.2171472409516259*log( c_norm(z) ), c_arg( z ) );
}    

@*1 Arrays of complex numbers.  

This assumes zero based arrays.

@ @<Prototype for |new_carray|@>=
struct c_complex * new_carray(long size)

@ @<Definition for |new_carray|@>=
		@<Prototype for |new_carray|@>
{
  struct c_complex *a;

  if (size<=0) complex_error("Non-positive complex array size chosen");
  
  a = (struct c_complex *) calloc (sizeof(struct c_complex), (unsigned long) size);

  if (a==NULL) complex_error("Can't allocate complex array");
  return a;
}

@ @<Prototype for |free_carray|@>=
void free_carray(struct c_complex *a)

@ @<Definition for |free_carray|@>=
		@<Prototype for |free_carray|@>
{
 if (a!=NULL) free(a);
}

@ This allocates a new complex array and 
copies the contents of |a| into it.

@<Prototype for |copy_carray|@>=
struct c_complex *copy_carray(struct c_complex *a, long size)

@ @<Definition for |copy_carray|@>=
		@<Prototype for |copy_carray|@>
{
  struct c_complex *b = NULL;
  
  if (a==NULL) complex_error("Can't duplicate a NULL complex array");
  
  b = new_carray(size);
  if (b!=NULL) memcpy(b, a, size * sizeof(struct c_complex));
   return b;
}

@ This puts |z| in all the entries in a complex array.

@<Prototype for |set_carray|@>=
void set_carray(struct c_complex * a, long size, struct c_complex z)

@ @<Definition for |set_carray|@>=
		@<Prototype for |set_carray|@>
{
  long j;

  if (a==NULL) complex_error("Can't operate on a NULL complex array");
  
  for (j=0; j< size; j++) a[j] = z;
}

