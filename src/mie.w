@** Mie Scattering Algorithms.

\def\ahat{\hat a}
\def\bhat{\hat b}
\def\Real{\mathop{\hbox{Re}}\nolimits}
\def\Imag{\mathop{\hbox{Im}}\nolimits}
\def\ds{\displaystyle}

This is a Mie scattering implementation.  Several resources were used
in creating this program.  First, the Fortran listing in Bohren and 
Huffman's book was used.  This listing was translated
into Pascal and refined using various suggestions by Wiscombe.  
This version was used for a couple of years and
later translated by me into C and then into CWeb with the documentation
you see here.  

Finally, consider using |ez_Mie| for problems that involve non-absorbing
spheres and you don't care about the scattering phase function.

A short to do list includes (1) use Wiscombe's trick to find the scattering functions,
(2) add code to deal with near zero entries in the Lentz routine, (3) 
allow calculation of extinction efficiencies with zero angles.

@ There are seven basic functions that are defined.  

@(mie.c@>=
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "mie_array.h"
#include "mie_complex.h"
#include "mie.h"

@<Definition for |mie_error|@>@;
@<Definition for |Lentz_Dn|@>@;
@<Definition for |Dn_down|@>@;
@<Definition for |Dn_up|@>@;
@<Definition for |small_Mie|@>@;
@<Definition for |small_conducting_Mie|@>@;
@<Definition for |Mie|@>@;
@<Definition for |ez_Mie|@>@;
@<Definition for |ez_Mie_Full|@>@;

@ Only the main function |Mie| is available for calling.
@(mie.h@>=
#define MIE_VERBOSE_ERROR_REPORTING 0
@<Prototype for |Lentz_Dn|@>;
@<Prototype for |Dn_down|@>;
@<Prototype for |Dn_up|@>;
@<Prototype for |small_Mie|@>;
@<Prototype for |small_conducting_Mie|@>;
@<Prototype for |Mie|@>;
@<Prototype for |ez_Mie|@>;
@<Prototype for |ez_Mie_Full|@>;

@ Some prototypes for the library interface.
@(libmie.h@>=
@<Prototype for |ez_Mie|@>;
@<Prototype for |ez_Mie_Full|@>;

@ A simple error routine that contains the only printf statement used
in the program.

@<Prototype for |mie_error|@>=
static void mie_error(char *s, int n)

@ @<Definition for |mie_error|@>=
		@<Prototype for |mie_error|@>
{
	if (MIE_VERBOSE_ERROR_REPORTING) {
		fprintf(stderr,"Mie Error %d -- %s\n", n, s);
		exit(n);
	}
}

@*1 The logarithmic derivative $D_n$.

@ This routine uses a continued fraction method to compute $D_n(z)$
proposed by Lentz.\footnote{*}{Lentz uses the notation $A_n$ 
instead of $D_n$, but I prefer the notation used by Bohren and Huffman.}
This method eliminates many weaknesses in previous algorithms using
forward recursion. 

I should add code to deal with $\alpha_{j,1}\approx 0$.

The logarithmic derivative $D_n$ is defined as
$$
D_n = -{n\over z} + {J_{n-1/2}(z)\over J_{n+1/2}(z)} 
$$
Equation (5) in Lentz's paper can be used to obtain
$$
{J_{n-1/2}(z)\over J_{n+1/2}(z)} =
{2n+1 \over z} + {1\over\ds -{2n+3\over z} +
                  {\strut 1 \over\ds {2n+5\over z} +
                      {\strut 1 \over\ds -{2n+7\over z} + \cdots}}}
$$

Now if
$$
\alpha_{i,j}=[a_i,a_{i-1},\ldots,a_j] = a_i + {1\over\ds a_{i-1} +
                  {\strut 1 \over\ds a_{i-2} + \cdots
                      {\strut 1 \over\ds a_j}}}
$$
we seek to create 
$$
\alpha = \alpha_{1,1}\,\alpha_{2,1}\cdots \alpha_{j,1} 
\qquad
\beta = \alpha_{2,2}\,\alpha_{3,2}\cdots \alpha_{j,2} 
$$
since Lentz showed that
$$
{J_{n-1/2}(z)\over J_{n+1/2}(z)} \approx {\alpha\over\beta}
$$

@ The whole goal is to iterate until the $\alpha$ and $\beta$
are identical to the number of digits desired. Once this is
achieved, then use equations this equation and the first equation for
the logarithmic derivative to calculate
$D_n(z)$.

@ @<Prototype for |Lentz_Dn|@>=
struct c_complex Lentz_Dn(struct c_complex z, long n)

@ @<Definition for |Lentz_Dn|@>=
		@<Prototype for |Lentz_Dn|@>
{
  struct c_complex alpha_j1, alpha_j2, zinv, aj;
  struct c_complex alpha, result,ratio,runratio;
  
@<Calculate first |alpha| and |beta|@>@;

do 
	@<Calculate next |ratio|@>@;

while (fabs(c_abs(ratio)-1.0) > 1e-12);

  result = c_add(c_sdiv((double)-n, z), runratio);
  return result;
}


@ Here I initialize for looping.  Of course it is kind of 
tricky, but what else would you expect.  The value of $a_j$ is
given by,
$$
a_j = (-1)^{j+1} {2n+2j-1\over z}
$$
The first terms for $\alpha$ and $beta$ are
$$
\alpha = a_1 \left(a_2 + {1\over a_1}\right)
\qquad
\beta = a_2
$$

@<Calculate first |alpha| and |beta|@>=

    zinv     = c_sdiv(2.0, z);
    alpha    = c_smul(n+0.5,zinv);
    aj       = c_smul(-n-1.5, zinv);
    alpha_j1 = c_add(aj, c_inv(alpha));
    alpha_j2 = aj;
    ratio    = c_div(alpha_j1,alpha_j2);
    runratio = c_mul(alpha,ratio);


@ To calculate the next $\alpha$ and $\beta$, I use
$$
a_{j+1} =  -a_j+(-1)^j{2\over z}
$$
to find the next $a_j$ and
$$
\alpha_{j+1} = a_j + {1\over\alpha_j},
\qquad\hbox{and}\qquad
\beta_{j+1} = a_j + {1\over\beta_j}
$$
and

@<Calculate next |ratio|@>=
{
    aj.re = zinv.re-aj.re;
    aj.im = zinv.im-aj.im;
    alpha_j1 = c_add(c_inv(alpha_j1), aj);
    alpha_j2 = c_add(c_inv(alpha_j2), aj);
    ratio = c_div(alpha_j1, alpha_j2);
    zinv.re *= -1;
    zinv.im *= -1;
    runratio = c_mul(ratio,runratio);
}

@*2 $D_n$ by upward recurrence.

Calculating the logarithmic derivative $D_n(\rho)$ using the upward
recurrence relation,
$$
D_n(z) = {1\over n/z - D_{n-1}(z)}-{n\over z}
$$

@ To calculate the initial value we must figure out
$D_0(z)$.  This is
$$
D_0(z) = {d\over dz} \ln\psi_0(z) = {d\over dz} \ln\sin(z)={\cos z\over\sin z}
$$
The only tricky part is finding the tangent of a complex number, but
this is all stuck in \.{complex.w}.

Finally, note that the returned array |*D| is set-up so that
$D_n(z)=$|D[n]|.  Therefore the first value for $D_1(z)$ will be
found not in |D[0]|, but rather in |D[1]|.
  

@ @<Prototype for |Dn_up|@>=
void Dn_up(struct c_complex z, long nstop, struct c_complex *D)

@ @<Definition for |Dn_up|@>=
		@<Prototype for |Dn_up|@>
{
  struct c_complex zinv, k_over_z;
  long k;
  
  D[0] = c_inv(c_tan(z));
  zinv = c_inv(z);

  for (k = 1; k < nstop; k++) {
    k_over_z = c_smul((double)k,zinv);
    D[k] = c_sub(c_inv(c_sub(k_over_z,D[k-1])),k_over_z);
  }
}

@*2 $D_n$ by downwards recurrence.

Start downwards recurrence using Lentz method, then find earlier
terms of the logarithmic derivative $D_n(z)$ using the recurrence relation,
$$
D_{n-1}(z) = {n\over z} - {1\over D_n(z) + n/z }
$$
This is a pretty straightforward procedure.
  
Finally, note that the returned array |*D| is set-up so that
$D_n(z)=$|D[n]|.  Therefore the first value for $D_1(z)$ will be
found not in |D[0]|, but rather in |D[1]|.
  
@ @<Prototype for |Dn_down|@>=
void Dn_down(struct c_complex z, long nstop, struct c_complex *D)

@ @<Definition for |Dn_down|@>=
		@<Prototype for |Dn_down|@>
{
  long k;
  struct c_complex zinv, k_over_z;

  D[nstop-1] = Lentz_Dn(z, nstop);
  zinv = c_inv(z);

  for (k = nstop-1; k >= 1; k--) {
    k_over_z = c_smul((double)k,zinv);
    D[k-1] = c_sub(k_over_z, c_inv(c_add(D[k], k_over_z)));
  }
}
@*1 Small Spheres.

This calculates everything accurately for small spheres.  This approximation
is necessary because in the small particle or Rayleigh limit $x\rightarrow0$ the
Mie formulas become ill-conditioned.  The method was taken from Wiscombe's paper
and has been tested for several complex indices of refraction.       
Wiscombe uses this when 
$$
x\vert m\vert\le0.1
$$ 
and says this routine should be accurate to six places.  

If |nangles==0| or |s1==NULL| or |s2==NULL| then this routine will do the
right thing---it will calculate the efficiencies and the anisotropy, but will
not calculate any of the scattering amplitudes.

Since it is not obvious |z0|$= i(m^2-1)$

@ @<Prototype for |small_Mie|@>=
void small_Mie(double x, struct c_complex m, double * mu, 
		      long nangles, struct c_complex *s1, 
		      struct c_complex *s2, double *qext, double *qsca, 
		      double *qback, double *g)

@ @<Definition for |small_Mie|@>=
		@<Prototype for |small_Mie|@>

{
  struct c_complex ahat1,ahat2,bhat1;
  struct c_complex z0,m2,m4;
  double x2,x3,x4;

  if ((s1==NULL) || (s2==NULL)) nangles=0;

  m2 = c_sqr(m);
  m4 = c_sqr(m2);
  x2 = x * x;
  x3 = x2 * x;
  x4 = x2 * x2;
  z0.re = - m2.im;
  z0.im =   m2.re-1;

  @<Calculate $\ahat_1$@>@;
  @<Calculate $\bhat_1$@>@;
  @<Calculate $\ahat_2$@>@;
  @<Calculate small Mie efficiencies and asymmetry@>@;
  @<Calculate small Mie scattering amplitudes@>@;
}

@ The formula for $\ahat_1$ is
$$
\ahat_1 = 2i{m^2-1\over3}{1-0.1x^2+{\ds4m^2+5\over\ds1400}x^4\over D}
$$
where 
$$
D=m^2+2+(1-0.7m^2)x^2-{8m^4-385m^2+350\over1400}x^4+2i{m^2-1\over3}x^3(1-0.1x^2)
$$
Note that I have disabled the case when the sphere has no index of refraction.
The perfectly conducting sphere equations are 
@<Calculate $\ahat_1$@>=
{ struct c_complex z1, z2, z3, z4, D;

  z1   = c_smul(2.0/3.0, z0);
  z2.re= 1.0-0.1*x2+(4.0*m2.re+5.0)*x4/1400.0;
  z2.im= 4.0*x4*m2.im/1400.0;
  z3   = c_mul(z1, z2);
  
  z4   = c_smul( x3*(1.0-0.1*x2), z1);
  D.re = 2.0+ m2.re + (1-0.7*m2.re)*x2 - (8.0*m4.re - 385.0*m2.re + 350.0)/1400.0*x4 + z4.re;
  D.im =      m2.im + ( -0.7*m2.im)*x2 - (8.0*m4.im - 385.0*m2.im        )/1400.0*x4 + z4.im;

  ahat1 = c_div(z3,D);

}

@ The formula for $\bhat_1$ is
$$
\bhat_1 = ix^2{m^2-1\over45}{1+{\ds2m^2-5\over\ds70}x^2\over1-{\ds2m^2-5\over\ds30}x^2}
$$
@<Calculate $\bhat_1$@>=
{
struct c_complex z2, z6, z7;
  z2   = c_smul(x2/45.0, z0);
  z6.re = 1.0 + (2.0 * m2.re - 5.0)*x2 / 70.0;
  z6.im = m2.im * x2 / 35.0;
  z7.re = 1.0 - (2.0 * m2.re - 5.0) *x2 / 30.0;
  z7.im = - m2.im * x2 / 15.0;
  bhat1 = c_mul( z2, c_div(z6,z7));
}

@ The formula for $\ahat_2$ is
$$
\ahat_2 = ix^2{m^2-1\over15}{1-{\ds1\over\ds14}x^2\over2m^2+3-{\ds2m^2-7\over\ds14}x^2}
$$
@<Calculate $\ahat_2$@>=
{ struct c_complex z3, z8;

  z3 = c_smul((1.0-x2/14.0)*x2/15.0,z0);
  z8.re = 2.0 * m2.re + 3.0 - (m2.re/7.0 - 0.5) * x2;
  z8.im = 2.0 * m2.im - m2.im / 7.0 * x2;
  ahat2 = c_div(z3, z8);

}

@ The scattering and extinction efficiencies are given by
$$
Q_{{\rm  ext}} = 6x \Real\left[\ahat_1+\bhat_1+{5\over3}\ahat_2\right]
$$
and
$$
Q_{{\rm  sca}} = 6x^4 T 
$$
with
$$
T         =\vert\ahat_1\vert^2+\vert\bhat_1\vert^2+{5\over3}\vert\ahat_2\vert^2
$$
and the anisotropy (average cosine of the phase function) is
$$
g          ={1\over T}\Real\left[\ahat_1(\ahat_2+\bhat_1)^*\right] 
$$

I also calculate the backscattering efficiency so that it
will be calculated correctly even when |nangles==0|.  
The backscattering efficiency $Q_{{\rm  back}}$ is defined as
$$
Q_{\rm  back} = {\sigma_{{\rm  back}}\over\pi a^2}
                  = {\vert S_1(-1)\vert^2 \over x^2}
$$
where $\sigma_{{\rm  back}}$ is the backscattering cross section.
The expression for $S_1(\mu)$ given in the chunk below yields 
$$
{S_1(-1)\over x}={3\over2}x^2\left[\ahat_1-\bhat_1-{5\over3}\ahat_2\right] 
$$
This only remains to be squared before the efficiency for backscattering
is obtained.

@ @<Calculate small Mie efficiencies and asymmetry@>=
{struct c_complex ss1;
 double T;

  T = c_norm(ahat1) + c_norm(bhat1) + (5.0/3.0)*c_norm(ahat2);
  *qsca = 6.0 * x4 * T;
  *qext = 6.0 * x * (ahat1.re + bhat1.re + (5.0/3.0) * ahat2.re);
  *g = (ahat1.re * (ahat2.re + bhat1.re) + ahat1.im * (ahat2.im + bhat1.im)) / T;
  ss1.re = 1.5*x2*(ahat1.re-bhat1.re-(5.0/3.0)*ahat2.re);
  ss1.im = 1.5*x2*(ahat1.im-bhat1.im-(5.0/3.0)*ahat2.im);
  *qback = 4*c_norm(ss1);  
}

@ Here is where the scattering functions get calculated according to
$$
S_1(\mu) = {3\over2}x^3\left[\ahat_1+\left(\bhat_1+{5\over3}\ahat_2\right)\mu\right] 
\qquad
S_2(\mu) = {3\over2}x^3\left[\bhat_1+\ahat_1\mu+ {5\over3}\ahat_2(2\mu^2-1)\right]
$$
Since this is the last thing to get calculated, I take the liberty of mucking around
with the variables $\ahat_1$, $\bhat_1$, $\ahat_2$, and $x^3$

@<Calculate small Mie scattering amplitudes@>=
{
double muj, angle;
long j;
  x3 *= 1.5;
  ahat1.re *= x3;
  ahat1.im *= x3;
  bhat1.re *= x3;
  bhat1.im *= x3;
  ahat2.re *= x3*(5.0/3.0);
  ahat2.im *= x3*(5.0/3.0);
  for (j = 0; j < nangles; j++) {
    muj = mu[j];
    angle = 2.0 *muj*muj- 1.0;
    s1[j].re = ahat1.re + (bhat1.re + ahat2.re) * muj ;
    s1[j].im = ahat1.im + (bhat1.im + ahat2.im) * muj;
    s2[j].re = bhat1.re + ahat1.re * muj + ahat2.re * angle;
    s2[j].im = bhat1.im + ahat1.im * muj + ahat2.im * angle;
  }
}

@*1 Small Perfectly Conducting Spheres.

@ @<Prototype for |small_conducting_Mie|@>=
void small_conducting_Mie(double x, struct c_complex m, double * mu, 
		      long nangles, struct c_complex *s1, 
		      struct c_complex *s2, double *qext, double *qsca, 
		      double *qback, double *g)

@ @<Definition for |small_conducting_Mie|@>=
		@<Prototype for |small_conducting_Mie|@>

{
  struct c_complex ahat1,ahat2,bhat1,bhat2;
  struct c_complex ss1;
  double x2,x3,x4,muj, angle;
  long j;

  if ((s1==NULL) || (s2==NULL)) nangles=0;

  m.re += 0.0;  /* suppress warning */
  x2 = x * x;
  x3 = x2 * x;
  x4 = x2 * x2;

  ahat1 = c_div(c_set(0.0,2.0/3.0*(1.0-0.2*x2)),c_set(1.0-0.5*x2,2.0/3.0*x3));
  bhat1 = c_div(c_set(0.0,(x2-10.0)/30.0), c_set(1+0.5*x2,-x3/3.0));
  ahat2 = c_set(0.0, x2/30.);
  bhat2 = c_set(0.0,-x2/45.);
  
  *qsca = 6.0 * x4 * (c_norm(ahat1) + c_norm(bhat1) + 
                      (5.0/3.0)*(c_norm(ahat2)+c_norm(bhat2)));
  *qext = *qsca;
  *g    = 6.0 * x4 * (ahat1.im * (        ahat2.im+bhat1.im) + 
                      bhat2.im * (5.0/9.0*ahat2.im+bhat1.im) +
                      ahat1.re * bhat1.re) / (*qsca);
                        
  ss1.re = 1.5*x2*(ahat1.re-bhat1.re);
  ss1.im = 1.5*x2*(ahat1.im-bhat1.im-(5.0/3.0)*(ahat2.im+bhat2.im));
  *qback = 4*c_norm(ss1);  

  x3 *= 1.5;
  ahat1.re *= x3;
  ahat1.im *= x3;
  bhat1.re *= x3;
  bhat1.im *= x3;
  ahat2.im *= x3*(5.0/3.0);
  bhat2.im *= x3*(5.0/3.0);
  for (j = 0; j < nangles; j++) {
    muj = mu[j];
    angle = 2.0 *muj*muj- 1.0;
    s1[j].re = ahat1.re + (bhat1.re           ) * muj;
    s1[j].im = ahat1.im + (bhat1.im + ahat2.im) * muj + bhat2.im * angle;;
    s2[j].re = bhat1.re + (ahat1.re           ) * muj;
    s2[j].im = bhat1.im + (ahat1.im + bhat2.im) * muj + ahat2.im * angle;
  }
}

@*1 Arbitrary Spheres.

Calculates the amplitude scattering matrix elements and efficiencies for 
extinction, total scattering and backscattering for a given size parameter 
and relative refractive index.  The basic algorithm follows Bohren and Huffman
originally written in Fortran.  The code was translated into
CWeb and documented by Scott Prahl. 

Many improvements suggested by Wiscombe have been incorporated.
In particular, either upward or downward iteration will be used to calculate the 
lograthmic derivative $D_n(z)$.

Routine preliminary checking suggests that everything is being calculated ok 
except $g$. 

Space must have been allocated for the scattering amplitude angles |s1| and
|s2| before this routine is called.

@ @<Prototype for |Mie|@>=
void Mie(double x, struct c_complex m, double *mu, long nangles, struct c_complex *s1, 
           struct c_complex *s2, double *qext, double *qsca, double *qback, double *g)

@ @<Definition for |Mie|@>=
		@<Prototype for |Mie|@>

{
  @<Declare variables for |Mie|@>@;
  
  @<Catch bogus input values@>@;
  @<Deal with small spheres@>@;

  @<Calculate |nstop|@>@;  

  @<Mie allocate and initialize angle arrays@>@;
  if (m.re>0)
  	@<Calculate the logarithmic derivatives@>@;

  @<Prepare to sum over all |nstop| terms@>@;

  for (n = 1; n <= nstop; n++) { 
    @<Establish $a_n$ and $b_n$@>@;
    @<Calculate phase function for each angle@>@;
    @<Increment cross sections @>@;
    @<Prepare for the next iteration@>@;
  }

 @<Calculate Efficiencies@>@;
 @<Free allocated memory@>@;
}

@   @<Declare variables for |Mie|@>=
  struct c_complex * D;
  struct c_complex  z1, an, bn, bnm1, anm1, qbcalc;
  double * pi0, * pi1, * tau;
  struct c_complex xi, xi0, xi1;
  double psi,psi0,psi1;
  double alpha, beta, factor;
  long n, k, nstop, sign;
  *qext=-1;
  *qsca=-1;
  *qback=-1;
  *g= -1;
  
@  @<Catch bogus input values@>=
  if (m.im>0.0) {
  	mie_error("This program requires m.im>=0",1);
  	return;
  }
  if (x<=0.0) {
  	mie_error("This program requires positive sphere sizes",2);
  	return;
  }
  if (nangles<0) {
  	mie_error("This program requires non-negative angle sizes",3);
  	return;
  }
  if (nangles<0) {
  	mie_error("This program requires non-negative angle sizes",4);
  	return;
  }
  if ((nangles > 0)&&(s1 == NULL)) {
  	mie_error("Space must be allocated for s1 if nangles!=0",5);
  	return;
  }
  if ((nangles > 0)&&(s2 == NULL)) {
  	mie_error("Space must be allocated for s2if nangles!=0",6);
  	return;
  }
  if (x>20000) { 
  	mie_error("Program not validated for spheres with x>20000",7);
  	return;
  }

@   @<Deal with small spheres@>=
  if ((m.re==0) && (x<0.1)) {
    small_conducting_Mie(x, m, mu, nangles, s1, s2, qext, qsca, qback, g);
    return;
  }

  if ((m.re>0.0) && (c_abs(m) * x < 0.1)) {
    small_Mie(x, m, mu, nangles, s1, s2, qext, qsca, qback, g);
    return;
  }

@ @<Mie allocate and initialize angle arrays@>=
  if (nangles > 0) {
	  set_carray(s1, nangles, c_set(0.0,0.0));
	  set_carray(s2, nangles, c_set(0.0,0.0));
 
	  pi0 = new_darray(nangles);                 
	  pi1 = new_darray(nangles);
	  tau = new_darray(nangles);
	  
	  set_darray (pi0, nangles, 0.0);
	  set_darray (tau, nangles, 0.0);
	  set_darray (pi1, nangles, 1.0);
 }
  
@ Calculate number of terms to be summed in series after Wiscombe

@<Calculate |nstop|@>=  
  nstop = floor(x + 4.05 * pow(x, 0.33333) + 2.0);

@  Allocate and initialize the space for the arrays.  One noteworthy aspect is that the 
complex array |D| is allocated from 0 to |nstop|.  This allows 
|D| to be a one-based array from 1 to |nstop| instead of a zero-based 
array from 0 to |nstop-1|.  Therefore |D[n]| will directly
correspond to $D_n$ in Bohren.  Furthermore, |an| and |bn| will
correspond to $a_n$ and $b_n$.  The angular arrays are still
zero-based.

Use formula 7 from Wiscombe's paper to figure out if upwards or
downwards recurrence should be used.  Namely if
$$
m_{\rm Im}x\le 13.78 m_{\rm Re}^2 - 10.8 m_{\rm Re} + 3.9
$$
the upward recurrence would be stable.

@<Calculate the logarithmic derivatives@>=
{
struct c_complex z;

  z = c_smul(x,m);
  
  D = new_carray(nstop+1);
  if (D==NULL) {
  	mie_error("Cannot allocate log array",8); 
  	return;
  }

  if (fabs(m.im * x) < ((13.78 * m.re - 10.8) * m.re + 3.9))                       
    Dn_up(z, nstop, D);
  else
    Dn_down(z, nstop, D);
  }
  
@ OK,  Here we go.  We need to start up the arrays.  First, recall
(page 128 Bohren and Huffman) that
$$
\psi_n(x) = x j_n(x)\qquad\hbox{and}\qquad \xi_n(x) = x j_n(x) + i x y_n(x)
$$
where $j_n$ and $y_n$ are spherical Bessel functions.  The first few terms
may be worked out as,
$$
\psi_0(x) = \sin x 
\qquad\hbox{and}\qquad
\psi_1(x) = {\sin x\over x} - \cos x
$$
and
$$
\xi_0(x) = \psi_0 + i \cos x
\qquad\hbox{and}\qquad
\xi_1(x) = \psi_1 + i \left[{\cos x\over x} + \sin x\right]
$$

@<Prepare to sum over all |nstop| terms@>=
  psi0 = sin(x);
  psi1 = psi0/x - cos(x);
  xi0 = c_set(psi0, cos(x));
  xi1 = c_set(psi1, cos(x)/x+sin(x));
  *qsca = 0.0;
  *g = 0.0;
  *qext = 0.0;
  sign = 1;
  qbcalc = c_set(0.0,0.0);
  anm1  = c_set(0.0,0.0);
  bnm1  = c_set(0.0,0.0);
 
@ The main equations for $a_n$ and $b_n$ 
in Bohren and Huffman Equation (4.88).
$$
a_n = {\Big[ D_n(mx)/m + n/x\Big] \psi_n(x)-\psi_{n-1}(x)\over
       \Big[ D_n(mx)/m + n/x\Big] \xi_n(x)- \xi_{n-1}(x)}
$$
and
$$
b_n = {\Big[m D_n(mx) + n/x\Big] \psi_n(x)-\psi_{n-1}(x)\over
       \Big[m D_n(mx) + n/x\Big] \xi_n(x)- \xi_{n-1}(x)}
$$

@<Establish $a_n$ and $b_n$@>=
    if (m.re==0.0) {
    	an = c_sdiv(n*psi1/x - psi0, c_sub(c_smul(n/x,xi1),xi0));
	bn = c_sdiv(psi1, xi1);
    } else if (m.im==0.0) {
	z1.re = D[n].re/m.re +n/x;
        an = c_sdiv(z1.re*psi1-psi0,c_sub(c_smul(z1.re,xi1),xi0));
	
	z1.re = D[n].re*m.re +n/x;
        bn = c_sdiv(z1.re*psi1-psi0,c_sub(c_smul(z1.re,xi1),xi0));
    } else {
	z1 = c_div(D[n], m);
        z1.re += n/x;
        an = c_div(c_set(z1.re*psi1-psi0, z1.im*psi1),c_sub(c_mul(z1,xi1),xi0));
	
        z1 = c_mul(D[n], m);
        z1.re += n/x;
        bn = c_div(c_set(z1.re*psi1-psi0, z1.im*psi1),c_sub(c_mul(z1,xi1),xi0));
    }
    
@ The scattering matrix is given by Equation 4.74 in Bohren and Huffman.
Namely,
$$
S_1 = \sum_n {2n+1\over n(n+1)} (a_n \pi_n+b_n\tau_n)
$$
and
$$
S_2 = \sum_n {2n+1\over n(n+1)} (a_n \tau_n+b_n\pi_n)
$$

Furthermore, equation 4.47 in Bohren and Huffman states
$$
\pi_n = {2n-1\over n-1}\mu \pi_{n-1} - {n\over n-1} \pi_{n-2}
$$
and
$$
\tau_n = n\mu\pi_n-(n+1)\pi_{n-1}
$$
@<Calculate phase function for each angle@>=
    for (k = 0; k < nangles; k++) {
      factor = (2.0*n+1.0)/(n+1.0)/n;
      tau[k] =  n * mu[k] * pi1[k] - (n + 1) * pi0[k];
      alpha = factor * pi1[k];
      beta  = factor * tau[k];
      s1[k].re += alpha * an.re + beta * bn.re;
      s1[k].im += alpha * an.im + beta * bn.im;
      s2[k].re += alpha * bn.re + beta * an.re;
      s2[k].im += alpha * bn.im + beta * an.im;
    } 
    
    for (k = 0; k < nangles; k++) {
      factor = pi1[k];
      pi1[k] = ((2.0 * n + 1.0) * mu[k] * pi1[k] - (n + 1.0) * pi0[k]) / n;
      pi0[k] = factor;
    }
    
@ From page 120 of Bohren and Huffman the anisotropy is given by
$$
Q_{\rm sca}\langle \cos\theta\rangle = {4\over x^2} \left[
\sum_{n=1}^{\infty} {n(n+2)\over n+1} \Real\lbrace a_na_{n+1}^*+b_nb_{n+1}^*\rbrace
+ \sum_{n=1}^{\infty} {2n+1\over n(n+1)} \Real\lbrace a_nb_n^*\rbrace\right]
$$
For computation purposes, this must be rewritten as
$$
Q_{\rm sca}\langle \cos\theta\rangle = {4\over x^2} \left[
\sum_{n=2}^{\infty} {(n^2-1)\over n} \Real\lbrace a_{n-1}a_n^*+b_{n-1}b_n^*\rbrace
+ \sum_{n=1}^{\infty} {2n+1\over n(n+1)} \Real\lbrace a_nb_n^*\rbrace\right]
$$
From page 122 we find an expression for the backscattering efficiency
$$
Q_{\rm back} = {\sigma_b\over\pi a^2} = {1\over x^2} \left\vert
\sum_{n=1}^{\infty} (2n+1)(-1)^n(a_n-b_n)\right\vert^2
$$
From page 103 we find an expression for the scattering cross section
$$
Q_{\rm sca} = {\sigma_s\over\pi a^2}
= {2\over x^2}\sum_{n=1}^{\infty} (2n+1)(\vert a_n\vert^2+\vert b_n\vert^2)
$$
The total extinction efficiency is also found on page 103
$$
Q_{\rm ext}= {\sigma_t\over\pi a^2}
= {2\over x^2}\sum_{n=1}^{\infty} (2n+1)\Real(a_n+b_n)
$$

 @<Increment cross sections @>=
   factor     = 2.0 * n +1.0;
   *g+=(n-1.0/n)*(anm1.re*an.re+anm1.im*an.im+bnm1.re*bn.re+bnm1.im*bn.im);
   *g         += factor/n/(n+1.0)* (an.re * bn.re + an.im * bn.im);
   *qsca      += factor * (c_norm(an) + c_norm(bn));
   *qext      += factor * (an.re + bn.re);
    sign      *= -1;
    qbcalc.re += sign * factor * (an.re - bn.re);
    qbcalc.im += sign * factor * (an.im - bn.im);
    
@ The
recurrence relations for $\psi$ and $\xi$ depend on the recursion relations
for the spherical Bessel functions (page 96 equation 4.11)
$$
z_{n-1}(x) + z_{n+1}(x) = {2n+1\over x} z_n(x)
$$
where $z_n$ might be either $j_n$ or $y_n$.   Thus
$$
\psi_{n+1}(x) = {2n+1\over x} \psi_n(x) - \psi_{n-1}(x)
\qquad\hbox{and}\qquad
\xi_{n+1}(x) = {2n+1\over x} \xi_n(x) - \xi_{n-1}(x)
$$
Furthermore, 

@<Prepare for the next iteration@>=
    factor = (2.0 *n + 1.0)/x;
    xi   = c_sub(c_smul(factor, xi1), xi0);
    xi0  = xi1;
    xi1  = xi;

    psi = factor * psi1 - psi0;
    psi0 = psi1;
    psi1 = xi1.re;

    anm1 = an;
    bnm1 = bn;

@ @<Calculate Efficiencies@>=
  *qsca  *= 2 / (x*x);
  *qext  *= 2 / (x*x);
  *g     *= 4 / (*qsca) / (x*x);
  *qback = c_norm(qbcalc) / (x*x);

@ @<Free allocated memory@>=
  if (m.re>0) free_carray(D);
  
  if (nangles > 0) {
    free_darray(pi0);
    free_darray(pi1);
    free_darray(tau);
  }

@*1 Easy Mie.

Given the size and real index of refraction, calculate the scattering efficiency and the 
anisotropy for a non-absorbing sphere.  If the sphere is totally reflecting, then let
the index of refraction be equal to zero.

To recover the scattering coefficient $\mu_s$ from the efficiency |qsca| just
multiply |qsca| by the geometric cross sectional area and the density of scatterers.

@*2 The function |ez_Mie|.
@<Prototype for |ez_Mie|@>=
void ez_Mie(double x, double n, double *qsca, double *g)

@ @<Definition for |ez_Mie|@>=
		@<Prototype for |ez_Mie|@>
{
long nangles=0;
double *mu=NULL;
struct c_complex *s1=NULL;
struct c_complex *s2=NULL;
struct c_complex m; 
double qext, qback;

m.re = n;
m.im = 0.0;

Mie (x, m, mu, nangles, s1, s2, &qext, qsca, &qback, g);
}


@*2 The function |ez_Mie_Full|.
This is a simple interface to that provides complete access to Mie calculations.
This function will return the scattering functions $S_1(\mu)$ and $S_2(\mu)$ for
each of the angles in the array |mu|.  Note that these are the cosines of the
angles and not the angles in radians.

This routine assumes that memory has been allocated for the arrays 
|mu|, |s1_real|, |s1_imag|, |s2_real| and, |s2_imag|.  The number of elements 
in these arrays is specified by |nangles|.

If you do not want to mess with angles then you probably want to
call |ez_Mie| instead of this function.

@<Prototype for |ez_Mie_Full|@>=
void ez_Mie_Full(double x, double m_real, double m_imag, long nangles,  double *mu,
         double *s1_real, double *s1_imag, double *s2_real, double *s2_imag, 
         double *qext, double *qsca, double *qback, double *g)

@ @<Definition for |ez_Mie_Full|@>=
		@<Prototype for |ez_Mie_Full|@>
{
struct c_complex *s1=NULL;
struct c_complex *s2=NULL;
struct c_complex m; 
int i;

	m.re = m_real;
	m.im = m_imag;
	
	s1 = new_carray(nangles);
	s2 = new_carray(nangles);
	
	Mie (x, m, mu, nangles, s1, s2, qext, qsca, qback, g);
	
	for (i=0; i<nangles; i++) {
		s1_imag[i] = s1[i].im;
		s1_real[i] = s1[i].re;
		s2_imag[i] = s2[i].im;
		s2_real[i] = s2[i].re;
	}
		
	free_carray(s1);
	free_carray(s2);
}
