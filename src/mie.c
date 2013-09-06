/*2:*/
#line 26 "./mie.w"

#include <math.h> 
#include <stddef.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include "mie_array.h"
#include "mie_complex.h"
#include "mie.h"

/*6:*/
#line 68 "./mie.w"

/*5:*/
#line 65 "./mie.w"

static void mie_error(char*s,int n)

/*:5*/
#line 69 "./mie.w"

{
if(MIE_VERBOSE_ERROR_REPORTING){
fprintf(stderr,"Mie Error %d -- %s\n",n,s);
exit(n);
}
}

/*:6*/
#line 35 "./mie.w"

/*11:*/
#line 125 "./mie.w"

/*10:*/
#line 122 "./mie.w"

struct c_complex Lentz_Dn(struct c_complex z,long n)

/*:10*/
#line 126 "./mie.w"

{
struct c_complex alpha_j1,alpha_j2,zinv,aj;
struct c_complex alpha,result,ratio,runratio;

/*12:*/
#line 156 "./mie.w"


zinv= c_sdiv(2.0,z);
alpha= c_smul(n+0.5,zinv);
aj= c_smul(-n-1.5,zinv);
alpha_j1= c_add(aj,c_inv(alpha));
alpha_j2= aj;
ratio= c_div(alpha_j1,alpha_j2);
runratio= c_mul(alpha,ratio);


/*:12*/
#line 131 "./mie.w"


do
/*13:*/
#line 179 "./mie.w"

{
aj.re= zinv.re-aj.re;
aj.im= zinv.im-aj.im;
alpha_j1= c_add(c_inv(alpha_j1),aj);
alpha_j2= c_add(c_inv(alpha_j2),aj);
ratio= c_div(alpha_j1,alpha_j2);
zinv.re*= -1;
zinv.im*= -1;
runratio= c_mul(ratio,runratio);
}

/*:13*/
#line 134 "./mie.w"


while(fabs(c_abs(ratio)-1.0)> 1e-12);

result= c_add(c_sdiv((double)-n,z),runratio);
return result;
}


/*:11*/
#line 36 "./mie.w"

/*20:*/
#line 246 "./mie.w"

/*19:*/
#line 243 "./mie.w"

void Dn_down(struct c_complex z,long nstop,struct c_complex*D)

/*:19*/
#line 247 "./mie.w"

{
long k;
struct c_complex zinv,k_over_z;

D[nstop-1]= Lentz_Dn(z,nstop);
zinv= c_inv(z);

for(k= nstop-1;k>=1;k--){
k_over_z= c_smul((double)k,zinv);
D[k-1]= c_sub(k_over_z,c_inv(c_add(D[k],k_over_z)));
}
}
/*:20*/
#line 37 "./mie.w"

/*17:*/
#line 215 "./mie.w"

/*16:*/
#line 212 "./mie.w"

void Dn_up(struct c_complex z,long nstop,struct c_complex*D)

/*:16*/
#line 216 "./mie.w"

{
struct c_complex zinv,k_over_z;
long k;

D[0]= c_inv(c_tan(z));
zinv= c_inv(z);

for(k= 1;k<nstop;k++){
k_over_z= c_smul((double)k,zinv);
D[k]= c_sub(c_inv(c_sub(k_over_z,D[k-1])),k_over_z);
}
}

/*:17*/
#line 38 "./mie.w"

/*23:*/
#line 284 "./mie.w"

/*22:*/
#line 278 "./mie.w"

void small_Mie(double x,struct c_complex m,double*mu,
long nangles,struct c_complex*s1,
struct c_complex*s2,double*qext,double*qsca,
double*qback,double*g)

/*:22*/
#line 285 "./mie.w"


{
struct c_complex ahat1,ahat2,bhat1;
struct c_complex z0,m2,m4;
double x2,x3,x4;

if((s1==NULL)||(s2==NULL))nangles= 0;

m2= c_sqr(m);
m4= c_sqr(m2);
x2= x*x;
x3= x2*x;
x4= x2*x2;
z0.re= -m2.im;
z0.im= m2.re-1;

/*24:*/
#line 319 "./mie.w"

{struct c_complex z1,z2,z3,z4,D;

z1= c_smul(2.0/3.0,z0);
z2.re= 1.0-0.1*x2+(4.0*m2.re+5.0)*x4/1400.0;
z2.im= 4.0*x4*m2.im/1400.0;
z3= c_mul(z1,z2);

z4= c_smul(x3*(1.0-0.1*x2),z1);
D.re= 2.0+m2.re+(1-0.7*m2.re)*x2-(8.0*m4.re-385.0*m2.re+350.0)/1400.0*x4+z4.re;
D.im= m2.im+(-0.7*m2.im)*x2-(8.0*m4.im-385.0*m2.im)/1400.0*x4+z4.im;

ahat1= c_div(z3,D);

}

/*:24*/
#line 302 "./mie.w"

/*25:*/
#line 339 "./mie.w"

{
struct c_complex z2,z6,z7;
z2= c_smul(x2/45.0,z0);
z6.re= 1.0+(2.0*m2.re-5.0)*x2/70.0;
z6.im= m2.im*x2/35.0;
z7.re= 1.0-(2.0*m2.re-5.0)*x2/30.0;
z7.im= -m2.im*x2/15.0;
bhat1= c_mul(z2,c_div(z6,z7));
}

/*:25*/
#line 303 "./mie.w"

/*26:*/
#line 354 "./mie.w"

{struct c_complex z3,z8;

z3= c_smul((1.0-x2/14.0)*x2/15.0,z0);
z8.re= 2.0*m2.re+3.0-(m2.re/7.0-0.5)*x2;
z8.im= 2.0*m2.im-m2.im/7.0*x2;
ahat2= c_div(z3,z8);

}

/*:26*/
#line 304 "./mie.w"

/*28:*/
#line 396 "./mie.w"

{struct c_complex ss1;
double T;

T= c_norm(ahat1)+c_norm(bhat1)+(5.0/3.0)*c_norm(ahat2);
*qsca= 6.0*x4*T;
*qext= 6.0*x*(ahat1.re+bhat1.re+(5.0/3.0)*ahat2.re);
*g= (ahat1.re*(ahat2.re+bhat1.re)+ahat1.im*(ahat2.im+bhat1.im))/T;
ss1.re= 1.5*x2*(ahat1.re-bhat1.re-(5.0/3.0)*ahat2.re);
ss1.im= 1.5*x2*(ahat1.im-bhat1.im-(5.0/3.0)*ahat2.im);
*qback= 4*c_norm(ss1);
}

/*:28*/
#line 305 "./mie.w"

/*29:*/
#line 418 "./mie.w"

{
double muj,angle;
long j;
x3*= 1.5;
ahat1.re*= x3;
ahat1.im*= x3;
bhat1.re*= x3;
bhat1.im*= x3;
ahat2.re*= x3*(5.0/3.0);
ahat2.im*= x3*(5.0/3.0);
for(j= 0;j<nangles;j++){
muj= mu[j];
angle= 2.0*muj*muj-1.0;
s1[j].re= ahat1.re+(bhat1.re+ahat2.re)*muj;
s1[j].im= ahat1.im+(bhat1.im+ahat2.im)*muj;
s2[j].re= bhat1.re+ahat1.re*muj+ahat2.re*angle;
s2[j].im= bhat1.im+ahat1.im*muj+ahat2.im*angle;
}
}

/*:29*/
#line 306 "./mie.w"

}

/*:23*/
#line 39 "./mie.w"

/*32:*/
#line 447 "./mie.w"

/*31:*/
#line 441 "./mie.w"

void small_conducting_Mie(double x,struct c_complex m,double*mu,
long nangles,struct c_complex*s1,
struct c_complex*s2,double*qext,double*qsca,
double*qback,double*g)

/*:31*/
#line 448 "./mie.w"


{
struct c_complex ahat1,ahat2,bhat1,bhat2;
struct c_complex ss1;
double x2,x3,x4,muj,angle;
long j;

if((s1==NULL)||(s2==NULL))nangles= 0;

m.re+= 0.0;
x2= x*x;
x3= x2*x;
x4= x2*x2;

ahat1= c_div(c_set(0.0,2.0/3.0*(1.0-0.2*x2)),c_set(1.0-0.5*x2,2.0/3.0*x3));
bhat1= c_div(c_set(0.0,(x2-10.0)/30.0),c_set(1+0.5*x2,-x3/3.0));
ahat2= c_set(0.0,x2/30.);
bhat2= c_set(0.0,-x2/45.);

*qsca= 6.0*x4*(c_norm(ahat1)+c_norm(bhat1)+
(5.0/3.0)*(c_norm(ahat2)+c_norm(bhat2)));
*qext= *qsca;
*g= 6.0*x4*(ahat1.im*(ahat2.im+bhat1.im)+
bhat2.im*(5.0/9.0*ahat2.im+bhat1.im)+
ahat1.re*bhat1.re)/(*qsca);

ss1.re= 1.5*x2*(ahat1.re-bhat1.re);
ss1.im= 1.5*x2*(ahat1.im-bhat1.im-(5.0/3.0)*(ahat2.im+bhat2.im));
*qback= 4*c_norm(ss1);

x3*= 1.5;
ahat1.re*= x3;
ahat1.im*= x3;
bhat1.re*= x3;
bhat1.im*= x3;
ahat2.im*= x3*(5.0/3.0);
bhat2.im*= x3*(5.0/3.0);
for(j= 0;j<nangles;j++){
muj= mu[j];
angle= 2.0*muj*muj-1.0;
s1[j].re= ahat1.re+(bhat1.re)*muj;
s1[j].im= ahat1.im+(bhat1.im+ahat2.im)*muj+bhat2.im*angle;;
s2[j].re= bhat1.re+(ahat1.re)*muj;
s2[j].im= bhat1.im+(ahat1.im+bhat2.im)*muj+ahat2.im*angle;
}
}

/*:32*/
#line 40 "./mie.w"

/*35:*/
#line 518 "./mie.w"

/*34:*/
#line 514 "./mie.w"

void Mie(double x,struct c_complex m,double*mu,long nangles,struct c_complex*s1,
struct c_complex*s2,double*qext,double*qsca,double*qback,double*g)

/*:34*/
#line 519 "./mie.w"


{
/*36:*/
#line 546 "./mie.w"

struct c_complex*D;
struct c_complex z1,an,bn,bnm1,anm1,qbcalc;
double*pi0,*pi1,*tau;
struct c_complex xi,xi0,xi1;
double psi,psi0,psi1;
double alpha,beta,factor;
long n,k,nstop,sign;
*qext= -1;
*qsca= -1;
*qback= -1;
*g= -1;

/*:36*/
#line 522 "./mie.w"


/*37:*/
#line 559 "./mie.w"

if(m.im> 0.0){
mie_error("This program requires m.im>=0",1);
return;
}
if(x<=0.0){
mie_error("This program requires positive sphere sizes",2);
return;
}
if(nangles<0){
mie_error("This program requires non-negative angle sizes",3);
return;
}
if(nangles<0){
mie_error("This program requires non-negative angle sizes",4);
return;
}
if((nangles> 0)&&(s1==NULL)){
mie_error("Space must be allocated for s1 if nangles!=0",5);
return;
}
if((nangles> 0)&&(s2==NULL)){
mie_error("Space must be allocated for s2if nangles!=0",6);
return;
}
if(x> 20000){
mie_error("Program not validated for spheres with x>20000",7);
return;
}

/*:37*/
#line 524 "./mie.w"

/*38:*/
#line 589 "./mie.w"

if((m.re==0)&&(x<0.1)){
small_conducting_Mie(x,m,mu,nangles,s1,s2,qext,qsca,qback,g);
return;
}

if((m.re> 0.0)&&(c_abs(m)*x<0.1)){
small_Mie(x,m,mu,nangles,s1,s2,qext,qsca,qback,g);
return;
}

/*:38*/
#line 525 "./mie.w"


/*40:*/
#line 616 "./mie.w"

nstop= floor(x+4.05*pow(x,0.33333)+2.0);

/*:40*/
#line 527 "./mie.w"


/*39:*/
#line 600 "./mie.w"

if(nangles> 0){
set_carray(s1,nangles,c_set(0.0,0.0));
set_carray(s2,nangles,c_set(0.0,0.0));

pi0= new_darray(nangles);
pi1= new_darray(nangles);
tau= new_darray(nangles);

set_darray(pi0,nangles,0.0);
set_darray(tau,nangles,0.0);
set_darray(pi1,nangles,1.0);
}

/*:39*/
#line 529 "./mie.w"

if(m.re> 0)
/*41:*/
#line 634 "./mie.w"

{
struct c_complex z;

z= c_smul(x,m);

D= new_carray(nstop+1);
if(D==NULL){
mie_error("Cannot allocate log array",8);
return;
}

if(fabs(m.im*x)<((13.78*m.re-10.8)*m.re+3.9))
Dn_up(z,nstop,D);
else
Dn_down(z,nstop,D);
}

/*:41*/
#line 531 "./mie.w"


/*42:*/
#line 671 "./mie.w"

psi0= sin(x);
psi1= psi0/x-cos(x);
xi0= c_set(psi0,cos(x));
xi1= c_set(psi1,cos(x)/x+sin(x));
*qsca= 0.0;
*g= 0.0;
*qext= 0.0;
sign= 1;
qbcalc= c_set(0.0,0.0);
anm1= c_set(0.0,0.0);
bnm1= c_set(0.0,0.0);

/*:42*/
#line 533 "./mie.w"


for(n= 1;n<=nstop;n++){
/*43:*/
#line 696 "./mie.w"

if(m.re==0.0){
an= c_sdiv(n*psi1/x-psi0,c_sub(c_smul(n/x,xi1),xi0));
bn= c_sdiv(psi1,xi1);
}else if(m.im==0.0){
z1.re= D[n].re/m.re+n/x;
an= c_sdiv(z1.re*psi1-psi0,c_sub(c_smul(z1.re,xi1),xi0));

z1.re= D[n].re*m.re+n/x;
bn= c_sdiv(z1.re*psi1-psi0,c_sub(c_smul(z1.re,xi1),xi0));
}else{
z1= c_div(D[n],m);
z1.re+= n/x;
an= c_div(c_set(z1.re*psi1-psi0,z1.im*psi1),c_sub(c_mul(z1,xi1),xi0));

z1= c_mul(D[n],m);
z1.re+= n/x;
bn= c_div(c_set(z1.re*psi1-psi0,z1.im*psi1),c_sub(c_mul(z1,xi1),xi0));
}

/*:43*/
#line 536 "./mie.w"

/*44:*/
#line 734 "./mie.w"

for(k= 0;k<nangles;k++){
factor= (2.0*n+1.0)/(n+1.0)/n;
tau[k]= n*mu[k]*pi1[k]-(n+1)*pi0[k];
alpha= factor*pi1[k];
beta= factor*tau[k];
s1[k].re+= alpha*an.re+beta*bn.re;
s1[k].im+= alpha*an.im+beta*bn.im;
s2[k].re+= alpha*bn.re+beta*an.re;
s2[k].im+= alpha*bn.im+beta*an.im;
}

for(k= 0;k<nangles;k++){
factor= pi1[k];
pi1[k]= ((2.0*n+1.0)*mu[k]*pi1[k]-(n+1.0)*pi0[k])/n;
pi0[k]= factor;
}

/*:44*/
#line 537 "./mie.w"

/*45:*/
#line 780 "./mie.w"

factor= 2.0*n+1.0;
*g+= (n-1.0/n)*(anm1.re*an.re+anm1.im*an.im+bnm1.re*bn.re+bnm1.im*bn.im);
*g+= factor/n/(n+1.0)*(an.re*bn.re+an.im*bn.im);
*qsca+= factor*(c_norm(an)+c_norm(bn));
*qext+= factor*(an.re+bn.re);
sign*= -1;
qbcalc.re+= sign*factor*(an.re-bn.re);
qbcalc.im+= sign*factor*(an.im-bn.im);

/*:45*/
#line 538 "./mie.w"

/*46:*/
#line 804 "./mie.w"

factor= (2.0*n+1.0)/x;
xi= c_sub(c_smul(factor,xi1),xi0);
xi0= xi1;
xi1= xi;

psi= factor*psi1-psi0;
psi0= psi1;
psi1= xi1.re;

anm1= an;
bnm1= bn;

/*:46*/
#line 539 "./mie.w"

}

/*47:*/
#line 817 "./mie.w"

*qsca*= 2/(x*x);
*qext*= 2/(x*x);
*g*= 4/(*qsca)/(x*x);
*qback= c_norm(qbcalc)/(x*x);

/*:47*/
#line 542 "./mie.w"

/*48:*/
#line 823 "./mie.w"

if(m.re> 0)free_carray(D);

if(nangles> 0){
free_darray(pi0);
free_darray(pi1);
free_darray(tau);
}

/*:48*/
#line 543 "./mie.w"

}

/*:35*/
#line 41 "./mie.w"

/*51:*/
#line 845 "./mie.w"

/*50:*/
#line 842 "./mie.w"

void ez_Mie(double x,double n,double*qsca,double*g)

/*:50*/
#line 846 "./mie.w"

{
long nangles= 0;
double*mu= NULL;
struct c_complex*s1= NULL;
struct c_complex*s2= NULL;
struct c_complex m;
double qext,qback;

m.re= n;
m.im= 0.0;

Mie(x,m,mu,nangles,s1,s2,&qext,qsca,&qback,g);
}


/*:51*/
#line 42 "./mie.w"

/*53:*/
#line 880 "./mie.w"

/*52:*/
#line 875 "./mie.w"

void ez_Mie_Full(double x,double m_real,double m_imag,long nangles,double*mu,
double*s1_real,double*s1_imag,double*s2_real,double*s2_imag,
double*qext,double*qsca,double*qback,double*g)

/*:52*/
#line 881 "./mie.w"

{
struct c_complex*s1= NULL;
struct c_complex*s2= NULL;
struct c_complex m;
int i;

m.re= m_real;
m.im= m_imag;

s1= new_carray(nangles);
s2= new_carray(nangles);

Mie(x,m,mu,nangles,s1,s2,qext,qsca,qback,g);

for(i= 0;i<nangles;i++){
s1_imag[i]= s1[i].im;
s1_real[i]= s1[i].re;
s2_imag[i]= s2[i].im;
s2_real[i]= s2[i].re;
}

free_carray(s1);
free_carray(s2);
}/*:53*/
#line 43 "./mie.w"


/*:2*/
