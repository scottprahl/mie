/*1:*/
#line 16 "./mie_legendre.w"

#include"mie_legendre.h"

/*4:*/
#line 56 "./mie_legendre.w"

/*3:*/
#line 53 "./mie_legendre.w"

double LegendrePn(long n,double x)

/*:3*/
#line 57 "./mie_legendre.w"

{
double pk,pkp1,pkm1;
long k;

if(n<=0)return 1.0;
if(n==1)return x;

if(x>=1.0)return 1.0;
if(x<=-1.0)return(n%2)?-1.0:1.0;

pk= x;
pkm1= 1.0;
for(k= 1;k<n;k++){
pkp1= ((2*k+1)*x*pk-k*pkm1)/(k+1);
pkm1= pk;
pk= pkp1;
}

return pk;
}

/*:4*/
#line 19 "./mie_legendre.w"

/*6:*/
#line 84 "./mie_legendre.w"

/*5:*/
#line 81 "./mie_legendre.w"

void LegendrePn_and_Pnm1(long n,double x,double*Pnm1Val,double*PnVal)

/*:5*/
#line 85 "./mie_legendre.w"

{
long k;
double Pk,Pkp1;
double Pkm1= 1.0;

*Pnm1Val= 1.0;
*PnVal= 1.0;
if(x>=1.0)x= 1.0;
if(x<=-1.0)x= -1.0;

Pk= x;

for(k= 1;k<n;k++){
Pkp1= ((2*k+1)*x*Pk-k*Pkm1)/(k+1);
Pkm1= Pk;
Pk= Pkp1;
}

*Pnm1Val= Pkm1;
*PnVal= Pk;
}

/*:6*/
#line 20 "./mie_legendre.w"

/*8:*/
#line 124 "./mie_legendre.w"

/*7:*/
#line 121 "./mie_legendre.w"

double LegendrePnd(long n,double x)

/*:7*/
#line 125 "./mie_legendre.w"

{
double p,pminus,pplus;
long i;

if(n<=0)return 0;
if(n==1)return 1;

if(x> 1.0)x= 1.0;
if(x<-1.0)x= -1.0;

pminus= 0;
p= 1;

for(i= 1;i<n;i++){
pplus= ((2*i+1)*x*p-(i+1)*pminus)/i;
pminus= p;
p= pplus;
}
return p;
}


/*:8*/
#line 21 "./mie_legendre.w"

/*10:*/
#line 164 "./mie_legendre.w"

/*9:*/
#line 161 "./mie_legendre.w"

double LegendrePndd(long n,double x)

/*:9*/
#line 165 "./mie_legendre.w"

{
double p,pminus,pplus;
long m;

if(n<=1)return 0;
if(n==2)return 3;

if(x> 1.0)x= 1.0;
if(x<-1.0)x= -1.0;

pminus= 0;
p= 3;

for(m= 2;m<n;m++){
pplus= ((2*m+1)*x*p-(m+2)*pminus)/(m-1);
pminus= p;
p= pplus;
}
return p;

}/*:10*/
#line 22 "./mie_legendre.w"


/*:1*/
