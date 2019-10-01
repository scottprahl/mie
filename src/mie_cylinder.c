/*2:*/
#line 8 "mie_cylinder.w"

#include <math.h> 
#include <stdlib.h> 
#include <stdio.h> 
#include <stdarg.h> 
#include "mie_array.h"
#include "mie_complex.h"
#include "mie_cylinder.h"
#define  PI 3.14159265358979

/*5:*/
#line 38 "mie_cylinder.w"

/*4:*/
#line 35 "mie_cylinder.w"

static double bessj0(double x)

/*:4*/
#line 39 "mie_cylinder.w"

{
double ax,z;
double xx,y,ans,ans1,ans2;

if((ax= fabs(x))<8.0){
y= x*x;
ans1= 57568490574.0+y*(-13362590354.0+y*(651619640.7
+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
ans2= 57568490411.0+y*(1029532985.0+y*(9494680.718
+y*(59272.64853+y*(267.8532712+y*1.0))));
ans= ans1/ans2;
}else{
z= 8.0/ax;
y= z*z;
xx= ax-0.785398164;
ans1= 1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
+y*(-0.2073370639e-5+y*0.2093887211e-6)));
ans2= -0.1562499995e-1+y*(0.1430488765e-3
+y*(-0.6911147651e-5+y*(0.7621095161e-6
-y*0.934935152e-7)));
ans= sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
}
return ans;
}

/*:5*/
#line 18 "mie_cylinder.w"

/*9:*/
#line 102 "mie_cylinder.w"

/*8:*/
#line 99 "mie_cylinder.w"

static double bessj1(double x)

/*:8*/
#line 103 "mie_cylinder.w"


{
double ax,z;
double xx,y,ans,ans1,ans2;

if((ax= fabs(x))<8.0){
y= x*x;
ans1= x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
ans2= 144725228442.0+y*(2300535178.0+y*(18583304.74
+y*(99447.43394+y*(376.9991397+y*1.0))));
ans= ans1/ans2;
}else{
z= 8.0/ax;
y= z*z;
xx= ax-2.356194491;
ans1= 1.0+y*(0.183105e-2+y*(-0.3516396496e-4
+y*(0.2457520174e-5+y*(-0.240337019e-6))));
ans2= 0.04687499995+y*(-0.2002690873e-3
+y*(0.8449199096e-5+y*(-0.88228987e-6
+y*0.105787412e-6)));
ans= sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
if(x<0.0)ans= -ans;
}
return ans;
}

/*:9*/
#line 19 "mie_cylinder.w"

/*7:*/
#line 70 "mie_cylinder.w"

/*6:*/
#line 67 "mie_cylinder.w"

static double bessy0(double x)

/*:6*/
#line 71 "mie_cylinder.w"

{
double z;
double xx,y,ans,ans1,ans2;

if(x<8.0){
y= x*x;
ans1= -2957821389.0+y*(7062834065.0+y*(-512359803.6
+y*(10879881.29+y*(-86327.92757+y*228.4622733))));
ans2= 40076544269.0+y*(745249964.8+y*(7189466.438
+y*(47447.26470+y*(226.1030244+y*1.0))));
ans= (ans1/ans2)+0.636619772*bessj0(x)*log(x);
}else{
z= 8.0/x;
y= z*z;
xx= x-0.785398164;
ans1= 1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
+y*(-0.2073370639e-5+y*0.2093887211e-6)));
ans2= -0.1562499995e-1+y*(0.1430488765e-3
+y*(-0.6911147651e-5+y*(0.7621095161e-6
+y*(-0.934945152e-7))));
ans= sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
}
return ans;
}

/*:7*/
#line 20 "mie_cylinder.w"

/*11:*/
#line 136 "mie_cylinder.w"

/*10:*/
#line 133 "mie_cylinder.w"

static double bessy1(double x)

/*:10*/
#line 137 "mie_cylinder.w"

{
double z;
double xx,y,ans,ans1,ans2;

if(x<8.0){
y= x*x;
ans1= x*(-0.4900604943e13+y*(0.1275274390e13
+y*(-0.5153438139e11+y*(0.7349264551e9
+y*(-0.4237922726e7+y*0.8511937935e4)))));
ans2= 0.2499580570e14+y*(0.4244419664e12
+y*(0.3733650367e10+y*(0.2245904002e8
+y*(0.1020426050e6+y*(0.3549632885e3+y)))));
ans= (ans1/ans2)+0.636619772*(bessj1(x)*log(x)-1.0/x);
}else{
z= 8.0/x;
y= z*z;
xx= x-2.356194491;
ans1= 1.0+y*(0.183105e-2+y*(-0.3516396496e-4
+y*(0.2457520174e-5+y*(-0.240337019e-6))));
ans2= 0.04687499995+y*(-0.2002690873e-3
+y*(0.8449199096e-5+y*(-0.88228987e-6
+y*0.105787412e-6)));
ans= sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
}
return ans;
}

/*:11*/
#line 21 "mie_cylinder.w"


/*13:*/
#line 171 "mie_cylinder.w"

/*12:*/
#line 168 "mie_cylinder.w"

static void jn_real(double x,int n,double*BesselJn,double*BesselYn)

/*:12*/
#line 172 "mie_cylinder.w"

{
const double iacc= 40;
const double bigno= 1.0e+20;
const double bigni= 1.0e-20;
int j,m,jsum;
double tox,bj,bjp,bjm,sum;

BesselJn[0]= bessj0(x);
BesselJn[1]= bessj1(x);
BesselYn[0]= bessy0(x);
BesselYn[1]= bessy1(x);

tox= 2./x;
for(j= 1;j<n;j++)
BesselYn[j+1]= j*tox*BesselYn[j]-BesselYn[j-1];

if(x> n){
for(j= 1;j<n;j++)
BesselJn[j+1]= j*tox*BesselJn[j]-BesselJn[j-1];
return;
}

m= 2*floor((n+sqrt(iacc*n))/2);
jsum= 0;
sum= 0.0;
bjp= 0.0;
bj= 1.0;
for(j= m;j> 0;j--){
bjm= j*tox*bj-bjp;
bjp= bj;
bj= bjm;
if(fabs(bj)> bigno){
int i;
bj= bj*bigni;
bjp= bjp*bigni;
sum= sum*bigni;
for(i= j+1;i<=n;i++)
BesselJn[i]*= bigni;
}
if(jsum)
sum= sum+bj;
jsum= 1-jsum;

if(j<=n&&j>=2)
BesselJn[j]= bjp;
}

sum= 2.0*sum-bj;
for(j= 2;j<=n;j++)
BesselJn[j]/= sum;
}

/*:13*/
#line 23 "mie_cylinder.w"

/*15:*/
#line 230 "mie_cylinder.w"

/*14:*/
#line 227 "mie_cylinder.w"

static void jn_complex(struct c_complex z,int n,struct c_complex*Jn)

/*:14*/
#line 231 "mie_cylinder.w"

{
struct c_complex a,*JnTmp;
int nd,i;

nd= 2*floor((pow(101+c_abs(z),0.499)+n)/2);
JnTmp= new_carray((long)nd+1);
JnTmp[nd]= c_set(0.0,0.0);
JnTmp[nd-1]= c_set(1.0e-32,0.0);
a= c_set(0.0,0.0);

for(i= nd-1;i>=3;i-= 2){
JnTmp[i-1]= c_sub(c_smul(2.0*i,c_div(JnTmp[i],z)),JnTmp[i+1]);
JnTmp[i-2]= c_sub(c_smul(2.0*(i-1),c_div(JnTmp[i-1],z)),JnTmp[i]);
a.re+= JnTmp[i-1].re;
a.im+= JnTmp[i-1].im;
}

JnTmp[0]= c_sub(c_smul(2.0,c_div(JnTmp[1],z)),JnTmp[2]);
a.re= 2.0*a.re+JnTmp[0].re;
a.im= 2.0*a.im+JnTmp[0].im;
for(i= 0;i<n;i++)
Jn[i]= c_div(JnTmp[i],a);
free_carray(JnTmp);
}

/*:15*/
#line 24 "mie_cylinder.w"

/*17:*/
#line 276 "mie_cylinder.w"

/*16:*/
#line 269 "mie_cylinder.w"

void MieCylinderCoefficients(double x,struct c_complex m,double zeta,int n_terms,
struct c_complex an1[],
struct c_complex bn1[],
struct c_complex an2[],
struct c_complex bn2[])

/*:16*/
#line 277 "mie_cylinder.w"

{
struct c_complex*jn1,ci,eta,feta,m2xi;
double*BesselJn,*BesselYn,jnp,ynp,sin_zeta,xi,cos_zeta;
int i;
long size= n_terms+1;

jn1= new_carray(size);
BesselJn= new_darray(size);
BesselYn= new_darray(size);

cos_zeta= cos(zeta);
sin_zeta= sqrt((1.0+cos_zeta)*(1.0-cos_zeta));
ci= c_set(0.0,1.0);
eta= c_smul(x,c_sqrt(c_sub(c_sqr(m),c_set(cos_zeta*cos_zeta,0.0))));
xi= x*sin_zeta;
feta= c_mul(c_smul(cos_zeta,eta),c_sadd(-1.0,c_sdiv(xi*xi,c_sqr(eta))));

jn_real(xi,n_terms+1,BesselJn,BesselYn);
jn_complex(eta,n_terms+1,jn1);

m2xi= c_smul(xi,c_sqr(m));
for(i= 0;i<n_terms;i++){
struct c_complex hn,hnp,dn1,an,bn,cn,dn,vn,wn,den;
dn1= c_div(c_sub(c_div(c_smul((double)i,jn1[i]),eta),jn1[i+1]),jn1[i]);
jnp= i*BesselJn[i]/xi-BesselJn[i+1];
ynp= i*BesselYn[i]/xi-BesselYn[i+1];
hn= c_set(BesselJn[i],BesselYn[i]);
hnp= c_set(jnp,ynp);
dn= c_smul((double)i,c_mul(hn,feta));
cn= c_smul((double)i,c_smul(BesselJn[i],feta));
bn= c_smul(xi,c_sub(c_smul(BesselJn[i],c_mul(m2xi,dn1)),c_smul(jnp,eta)));
vn= c_smul(xi,c_sub(c_mul(hn,c_mul(m2xi,dn1)),c_mul(hnp,eta)));
wn= c_mul(c_smul(xi,ci),c_sub(c_mul(hnp,eta),c_smul(xi,c_mul(dn1,hn))));
an= c_mul(c_smul(-xi,ci),c_sub(c_smul(jnp,eta),c_smul(xi*BesselJn[i],dn1)));
den= c_add(c_mul(wn,vn),c_mul(ci,c_sqr(dn)));
an1[i]= c_div(c_sub(c_mul(cn,vn),c_mul(bn,dn)),den);
bn1[i]= c_div(c_add(c_mul(wn,bn),c_mul(c_mul(ci,dn),cn)),den);
an2[i]= c_div(c_sub(c_mul(c_mul(ci,dn),cn),c_mul(an,vn)),den);
bn2[i]= c_mul(c_add(c_mul(cn,wn),c_mul(an,dn)),c_div(ci,den));
bn2[i]= c_smul(-1.0,bn2[i]);
}

free_carray(jn1);
free_darray(BesselJn);
free_darray(BesselYn);
}

/*:17*/
#line 25 "mie_cylinder.w"

/*19:*/
#line 352 "mie_cylinder.w"

/*18:*/
#line 346 "mie_cylinder.w"

void MieCylinder(double x,struct c_complex m,double zeta,const double*theta,int nangles,
struct c_complex*t1,struct c_complex*t2,struct c_complex*t3,
struct c_complex*qexpar,struct c_complex*qexper,
double*qscpar,double*qscper)

/*:18*/
#line 353 "mie_cylinder.w"


{
struct c_complex*an1,*bn1,*an2,*bn2;
int i,k,n_stop_terms;
long size;

n_stop_terms= x+4.0*pow(x,1.0/3.0)+2.0;
size= n_stop_terms+1;
an1= new_carray(size);
bn1= new_carray(size);
an2= new_carray(size);
bn2= new_carray(size);
MieCylinderCoefficients(x,m,zeta,n_stop_terms,an1,bn1,an2,bn2);

*qexpar= c_smul(0.5,bn1[0]);
*qexper= c_smul(0.5,an2[0]);
*qscpar= 0.5*c_norm(bn1[0]);
*qscper= 0.5*c_norm(an2[0]);

for(i= 1;i<n_stop_terms;i++){
qexpar->re+= bn1[i].re;
qexpar->im+= bn1[i].im;
qexper->re+= an2[i].re;
qexper->im+= an2[i].im;
*qscpar+= c_norm(an1[i])+c_norm(bn1[i]);
*qscper+= c_norm(an2[i])+c_norm(bn2[i]);
}

*qscpar*= 4.0/x;
*qscper*= 4.0/x;
qexpar->re*= 4.0/x;
qexpar->im*= 4.0/x;
qexper->re*= 4.0/x;
qexper->im*= 4.0/x;

if(!theta)nangles= 0;

for(i= 0;i<nangles;i++){
double t= theta[i];
t1[i]= c_smul(0.5,bn1[0]);
t2[i]= c_smul(0.5,an2[0]);
t3[i]= c_set(0.0,0.0);
for(k= 1;k<=n_stop_terms;k++){
double ct= cos(k*t);
double st= sin(k*t);
t1[i].re+= ct*bn1[k].re;
t1[i].im+= ct*bn1[k].im;
t2[i].re+= ct*an2[k].re;
t2[i].im+= ct*an2[k].im;
t3[i].re+= st*an1[k].re;
t3[i].im+= st*an1[k].im;
}
}

free_carray(an1);
free_carray(an2);
free_carray(bn1);
free_carray(bn2);
}
/*:19*/
#line 26 "mie_cylinder.w"


/*:2*/
