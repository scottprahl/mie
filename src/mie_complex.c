/*2:*/
#line 10 "./mie_complex.w"

#include <math.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <float.h> 

#include "mie_complex.h"

/*6:*/
#line 120 "./mie_complex.w"

/*5:*/
#line 117 "./mie_complex.w"

static void complex_error(char*s)

/*:5*/
#line 121 "./mie_complex.w"

{
printf("%s\n",s);
exit(1);
}

/*:6*/
#line 19 "./mie_complex.w"


/*8:*/
#line 133 "./mie_complex.w"

/*7:*/
#line 130 "./mie_complex.w"

struct c_complex c_set(double a,double b)

/*:7*/
#line 134 "./mie_complex.w"

{
struct c_complex c;
c.re= a;
c.im= b;
return c;
}

/*:8*/
#line 21 "./mie_complex.w"

/*10:*/
#line 148 "./mie_complex.w"

/*9:*/
#line 145 "./mie_complex.w"

struct c_complex c_polarset(double r,double theta)

/*:9*/
#line 149 "./mie_complex.w"

{
return c_set(r*cos(theta),r*sin(theta));
}

/*:10*/
#line 22 "./mie_complex.w"

/*12:*/
#line 161 "./mie_complex.w"

/*11:*/
#line 158 "./mie_complex.w"

double c_abs(struct c_complex z)

/*:11*/
#line 162 "./mie_complex.w"

{
double x,y,temp;

x= fabs(z.re);
y= fabs(z.im);
if(x==0.0)return y;
if(y==0.0)return x;

if(x> y){
temp= y/x;
return(x*sqrt(1.0+temp*temp));
}

temp= x/y;
return(y*sqrt(1.0+temp*temp));
}

/*:12*/
#line 23 "./mie_complex.w"

/*16:*/
#line 194 "./mie_complex.w"

/*15:*/
#line 191 "./mie_complex.w"

double c_arg(struct c_complex z)

/*:15*/
#line 195 "./mie_complex.w"

{
return atan2(z.im,z.re);
}

/*:16*/
#line 24 "./mie_complex.w"

/*18:*/
#line 205 "./mie_complex.w"

/*17:*/
#line 202 "./mie_complex.w"

double c_norm(struct c_complex z)

/*:17*/
#line 206 "./mie_complex.w"

{
return(z.re*z.re+z.im*z.im);
}

/*:18*/
#line 25 "./mie_complex.w"

/*20:*/
#line 214 "./mie_complex.w"

/*19:*/
#line 211 "./mie_complex.w"

struct c_complex c_sqrt(struct c_complex z)

/*:19*/
#line 215 "./mie_complex.w"

{
double a,b;

if((z.re==0.0)&&(z.im==0.0))
return c_set(0.0,0.0);

a= sqrt((fabs(z.re)+c_abs(z))*0.5);
if(z.re>=0)
b= z.im/(a+a);
else{
b= z.im<0?-a:a;
a= z.im/(b+b);
}

return c_set(a,b);
}

/*:20*/
#line 26 "./mie_complex.w"

/*22:*/
#line 239 "./mie_complex.w"

/*21:*/
#line 236 "./mie_complex.w"

struct c_complex c_sqr(struct c_complex z)

/*:21*/
#line 240 "./mie_complex.w"

{
return c_mul(z,z);
}

/*:22*/
#line 27 "./mie_complex.w"

/*24:*/
#line 250 "./mie_complex.w"

/*23:*/
#line 247 "./mie_complex.w"

struct c_complex c_inv(struct c_complex w)

/*:23*/
#line 251 "./mie_complex.w"

{
double r,d;

if((w.re==0)&&(w.im==0))complex_error("Attempt to invert 0+0i");

if(fabs(w.re)>=fabs(w.im)){
r= w.im/w.re;
d= 1/(w.re+r*w.im);
return c_set(d,-r*d);
}

r= w.re/w.im;
d= 1/(w.im+r*w.re);
return c_set(r*d,-d);
}

/*:24*/
#line 28 "./mie_complex.w"

/*14:*/
#line 185 "./mie_complex.w"

/*13:*/
#line 182 "./mie_complex.w"

struct c_complex c_conj(struct c_complex z)

/*:13*/
#line 186 "./mie_complex.w"

{
return c_set(z.re,-z.im);
}

/*:14*/
#line 29 "./mie_complex.w"

/*27:*/
#line 275 "./mie_complex.w"

/*26:*/
#line 272 "./mie_complex.w"

struct c_complex c_add(struct c_complex z,struct c_complex w)

/*:26*/
#line 276 "./mie_complex.w"

{
struct c_complex c;

c.im= z.im+w.im;
c.re= z.re+w.re;
return c;
}

/*:27*/
#line 30 "./mie_complex.w"


/*29:*/
#line 290 "./mie_complex.w"

/*28:*/
#line 287 "./mie_complex.w"

struct c_complex c_sub(struct c_complex z,struct c_complex w)

/*:28*/
#line 291 "./mie_complex.w"

{
struct c_complex c;
c.im= z.im-w.im;
c.re= z.re-w.re;
return c;
}

/*:29*/
#line 32 "./mie_complex.w"

/*31:*/
#line 304 "./mie_complex.w"

/*30:*/
#line 301 "./mie_complex.w"

struct c_complex c_mul(struct c_complex z,struct c_complex w)

/*:30*/
#line 305 "./mie_complex.w"

{
struct c_complex c;
c.re= z.re*w.re-z.im*w.im;
c.im= z.im*w.re+z.re*w.im;
return c;
}

/*:31*/
#line 33 "./mie_complex.w"

/*33:*/
#line 319 "./mie_complex.w"

/*32:*/
#line 316 "./mie_complex.w"

struct c_complex c_div(struct c_complex z,struct c_complex w)

/*:32*/
#line 320 "./mie_complex.w"

{
struct c_complex c;
double r,denom;

if((w.re==0)&&(w.im==0))complex_error("Attempt to divide by 0+0i");

if(fabs(w.re)>=fabs(w.im)){
r= w.im/w.re;
denom= w.re+r*w.im;
c.re= (z.re+r*z.im)/denom;
c.im= (z.im-r*z.re)/denom;
}else{
r= w.re/w.im;
denom= w.im+r*w.re;
c.re= (z.re*r+z.im)/denom;
c.im= (z.im*r-z.re)/denom;
}
return c;
}

/*:33*/
#line 34 "./mie_complex.w"

/*35:*/
#line 347 "./mie_complex.w"

/*34:*/
#line 344 "./mie_complex.w"

double c_rdiv(struct c_complex z,struct c_complex w)

/*:34*/
#line 348 "./mie_complex.w"

{
double r,c,denom;

if((w.re==0)&&(w.im==0))complex_error("Attempt to find real part with divisor 0+0i");

if(fabs(w.re)>=fabs(w.im)){
r= w.im/w.re;
denom= w.re+r*w.im;
c= (z.re+r*z.im)/denom;
}else{
r= w.re/w.im;
denom= w.im+r*w.re;
c= (z.re*r+z.im)/denom;
}
return c;
}

/*:35*/
#line 35 "./mie_complex.w"

/*37:*/
#line 371 "./mie_complex.w"

/*36:*/
#line 368 "./mie_complex.w"

double c_rmul(struct c_complex z,struct c_complex w)

/*:36*/
#line 372 "./mie_complex.w"

{
return z.re*w.re-z.im*w.im;
}

/*:37*/
#line 36 "./mie_complex.w"


/*42:*/
#line 398 "./mie_complex.w"

/*41:*/
#line 395 "./mie_complex.w"

struct c_complex c_sadd(double x,struct c_complex z)

/*:41*/
#line 399 "./mie_complex.w"

{
struct c_complex c;
c.re= x+z.re;
c.im= z.im;
return c;
}

/*:42*/
#line 38 "./mie_complex.w"

/*44:*/
#line 413 "./mie_complex.w"

/*43:*/
#line 410 "./mie_complex.w"

struct c_complex c_sdiv(double x,struct c_complex w)

/*:43*/
#line 414 "./mie_complex.w"

{
struct c_complex c;
double r,factor;

if((w.re==0)&&(w.im==0))complex_error("Attempt to divide scalar by 0+0i");

if(fabs(w.re)>=fabs(w.im)){
r= w.im/w.re;
factor= x/(w.re+r*w.im);
c.re= factor;
c.im= -r*factor;
}else{
r= w.re/w.im;
factor= x/(w.im+r*w.re);
c.im= -factor;
c.re= r*factor;
}
return c;
}

/*:44*/
#line 39 "./mie_complex.w"

/*40:*/
#line 384 "./mie_complex.w"

/*39:*/
#line 381 "./mie_complex.w"

struct c_complex c_smul(double x,struct c_complex z)

/*:39*/
#line 385 "./mie_complex.w"

{
struct c_complex c;
c.re= z.re*x;
c.im= z.im*x;
return c;
}

/*:40*/
#line 40 "./mie_complex.w"


/*47:*/
#line 442 "./mie_complex.w"

/*46:*/
#line 439 "./mie_complex.w"

struct c_complex c_sin(struct c_complex z)

/*:46*/
#line 443 "./mie_complex.w"

{
return c_set(sin(z.re)*cosh(z.im),cos(z.re)*sinh(z.im));
}

/*:47*/
#line 42 "./mie_complex.w"

/*49:*/
#line 453 "./mie_complex.w"

/*48:*/
#line 450 "./mie_complex.w"

struct c_complex c_cos(struct c_complex z)

/*:48*/
#line 454 "./mie_complex.w"

{
return c_set(cos(z.re)*cosh(z.im),-(sin(z.re)*sinh(z.im)));
}

/*:49*/
#line 43 "./mie_complex.w"

/*51:*/
#line 484 "./mie_complex.w"

/*50:*/
#line 481 "./mie_complex.w"

struct c_complex c_tan(struct c_complex z)

/*:50*/
#line 485 "./mie_complex.w"

{
double t,x,y;

if(z.im==0)return c_set(tan(z.re),0.0);
if(z.im> DBL_MAX_10_EXP)return c_set(0.0,1.0);
if(z.im<-DBL_MAX_10_EXP)return c_set(0.0,-1.0);

x= 2*z.re;
y= 2*z.im;
t= cos(x)+cosh(y);
if(t==0)complex_error("Complex tangent is infinite");

return c_set(sin(x)/t,sinh(y)/t);
}

/*:51*/
#line 44 "./mie_complex.w"

/*53:*/
#line 506 "./mie_complex.w"

/*52:*/
#line 503 "./mie_complex.w"

struct c_complex c_asin(struct c_complex z)

/*:52*/
#line 507 "./mie_complex.w"

{
struct c_complex x;
x= c_log(c_add(c_set(-z.im,z.re),c_sqrt(c_sub(c_set(1.0,0.0),c_mul(z,z)))));
return c_set(x.im,-x.re);
}

/*:53*/
#line 45 "./mie_complex.w"

/*55:*/
#line 519 "./mie_complex.w"

/*54:*/
#line 516 "./mie_complex.w"

struct c_complex c_acos(struct c_complex z)

/*:54*/
#line 520 "./mie_complex.w"

{
struct c_complex x;
x= c_log(c_add(z,c_mul(c_set(0.0,1.0),c_sqrt(c_sub(c_set(1.0,0.0),c_sqr(z))))));
return c_set(x.im,-x.re);
}

/*:55*/
#line 46 "./mie_complex.w"

/*57:*/
#line 532 "./mie_complex.w"

/*56:*/
#line 529 "./mie_complex.w"

struct c_complex c_atan(struct c_complex z)

/*:56*/
#line 533 "./mie_complex.w"

{
struct c_complex x;
x= c_log(c_div(c_set(z.re,1+z.im),c_set(-z.re,1-z.im)));
return c_set(-x.im/2,x.re/2);
}

/*:57*/
#line 47 "./mie_complex.w"


/*62:*/
#line 554 "./mie_complex.w"

/*61:*/
#line 551 "./mie_complex.w"

struct c_complex c_sinh(struct c_complex z)

/*:61*/
#line 555 "./mie_complex.w"

{
return c_set(sinh(z.re)*cos(z.im),cosh(z.re)*sin(z.im));
}

/*:62*/
#line 49 "./mie_complex.w"

/*60:*/
#line 545 "./mie_complex.w"

/*59:*/
#line 542 "./mie_complex.w"

struct c_complex c_cosh(struct c_complex z)

/*:59*/
#line 546 "./mie_complex.w"

{
return c_set(cosh(z.re)*cos(z.im),sinh(z.re)*sin(z.im));
}

/*:60*/
#line 50 "./mie_complex.w"

/*64:*/
#line 563 "./mie_complex.w"

/*63:*/
#line 560 "./mie_complex.w"

struct c_complex c_tanh(struct c_complex z)

/*:63*/
#line 564 "./mie_complex.w"

{
double x= 2*z.re;
double y= 2*z.im;
double t= 1.0/(cosh(x)+cos(y));

return c_set(t*sinh(x),t*sin(y));
}

/*:64*/
#line 51 "./mie_complex.w"

/*66:*/
#line 576 "./mie_complex.w"

/*65:*/
#line 573 "./mie_complex.w"

struct c_complex c_atanh(struct c_complex z)

/*:65*/
#line 577 "./mie_complex.w"

{
return c_atan(c_set(-z.im,z.re));
}

/*:66*/
#line 52 "./mie_complex.w"

/*68:*/
#line 585 "./mie_complex.w"

/*67:*/
#line 582 "./mie_complex.w"

struct c_complex c_asinh(struct c_complex z)

/*:67*/
#line 586 "./mie_complex.w"

{
return c_asin(c_set(-z.im,z.re));
}

/*:68*/
#line 53 "./mie_complex.w"


/*71:*/
#line 596 "./mie_complex.w"

/*70:*/
#line 593 "./mie_complex.w"

struct c_complex c_exp(struct c_complex z)

/*:70*/
#line 597 "./mie_complex.w"

{
double x= exp(z.re);
return c_set(x*cos(z.im),x*sin(z.im));
}

/*:71*/
#line 55 "./mie_complex.w"

/*73:*/
#line 606 "./mie_complex.w"

/*72:*/
#line 603 "./mie_complex.w"

struct c_complex c_log(struct c_complex z)

/*:72*/
#line 607 "./mie_complex.w"

{
return c_set(log(c_abs(z)),c_arg(z));
}

/*:73*/
#line 56 "./mie_complex.w"

/*75:*/
#line 615 "./mie_complex.w"

/*74:*/
#line 612 "./mie_complex.w"

struct c_complex c_log10(struct c_complex z)

/*:74*/
#line 616 "./mie_complex.w"

{
return c_set(0.2171472409516259*log(c_norm(z)),c_arg(z));
}

/*:75*/
#line 57 "./mie_complex.w"


/*78:*/
#line 628 "./mie_complex.w"

/*77:*/
#line 625 "./mie_complex.w"

struct c_complex*new_carray(long size)

/*:77*/
#line 629 "./mie_complex.w"

{
struct c_complex*a;

if(size<=0)complex_error("Non-positive complex array size chosen");

a= (struct c_complex*)calloc(sizeof(struct c_complex),(unsigned long)size);

if(a==NULL)complex_error("Can't allocate complex array");
return a;
}

/*:78*/
#line 59 "./mie_complex.w"

/*80:*/
#line 644 "./mie_complex.w"

/*79:*/
#line 641 "./mie_complex.w"

void free_carray(struct c_complex*a)

/*:79*/
#line 645 "./mie_complex.w"

{
if(a!=NULL)free(a);
}

/*:80*/
#line 60 "./mie_complex.w"

/*82:*/
#line 656 "./mie_complex.w"

/*81:*/
#line 653 "./mie_complex.w"

struct c_complex*copy_carray(struct c_complex*a,long size)

/*:81*/
#line 657 "./mie_complex.w"

{
struct c_complex*b= NULL;

if(a==NULL)complex_error("Can't duplicate a NULL complex array");

b= new_carray(size);
if(b!=NULL)memcpy(b,a,size*sizeof(struct c_complex));
return b;
}

/*:82*/
#line 61 "./mie_complex.w"

/*84:*/
#line 673 "./mie_complex.w"

/*83:*/
#line 670 "./mie_complex.w"

void set_carray(struct c_complex*a,long size,struct c_complex z)

/*:83*/
#line 674 "./mie_complex.w"

{
long j;

if(a==NULL)complex_error("Can't operate on a NULL complex array");

for(j= 0;j<size;j++)a[j]= z;
}
/*:84*/
#line 62 "./mie_complex.w"


/*:2*/
