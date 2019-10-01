/*2:*/
#line 38 "mie_main.w"


/*3:*/
#line 68 "mie_main.w"


#include <stdio.h> 
#include <math.h> 
#include <stdlib.h> 
#include <string.h> 
#include "mie_array.h"
#include "mie_complex.h"
#include "mie.h"
#include "mygetopt.h"
#include "version.h"

extern char*optarg;
extern int optind;

/*:3*/
#line 40 "mie_main.w"

/*11:*/
#line 336 "mie_main.w"


static void print_version(void)
{
fprintf(stderr,"mie %s\n\n",Version);
fprintf(stderr,"Copyright (C) 2017 Free Software Foundation, Inc.\n");
fprintf(stderr,"This is free software; see the source for copying conditions.  There is NO\n");
fprintf(stderr,"warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.");
fprintf(stderr,"\n\nWritten by Scott Prahl\n");
exit(0);
}

/*:11*/
#line 41 "mie_main.w"

/*12:*/
#line 348 "mie_main.w"

static void print_usage(void)
{
fprintf(stderr,"mie %s\n\n",Version);
fprintf(stderr,"Calculates spherical Mie scattering phase function\n\n");
fprintf(stderr,"Usage:\n");
fprintf(stderr,"    mie [-l lambda] [-r radius] [-n index] ");
fprintf(stderr,"[-i imag index] [-d density] [-p phase angles]\n\n");
fprintf(stderr,"Options:\n");
fprintf(stderr,"    -o filename           # explicitly specify filename for output\n");
fprintf(stderr,"    -q                    # quiet --- omit output to stderr\n\n");
fprintf(stderr,"    -d  density           # density (spheres/micron^3)  [default=1.000]\n");
fprintf(stderr,"    -i  imag_index        # imag index of refraction    [default=0.000]\n");
fprintf(stderr,"    -l  lambda_vac        # wavelength in vacuum        [default=0.633]\n");
fprintf(stderr,"    -L  last_lambda       # last wavelength in vacuum   [default=0.633]\n");
fprintf(stderr,"    -m  index_of_medium   # refractive index of medium  [default=1.000]\n");
fprintf(stderr,"    -n  real_index        # real index of refraction    [default=1.550]\n");
fprintf(stderr,"    -p  num_of_angles     # number of angles            [default=0    ]\n");
fprintf(stderr,"    -P  num_of_lambda     # number of wavelengths        [default=0    ]\n");
fprintf(stderr,"    -r  radius            # sphere radius [microns]     [default=0.525]\n\n");
fprintf(stderr,"    -h                    # display help\n");
fprintf(stderr,"    -v                    # version information\n\n");
fprintf(stderr,"Examples:\n");
fprintf(stderr,"    mie -p 40         # Bohren & Huffman Appendix A\n");
fprintf(stderr,"    mie -d 1 -i 0 -l 0.6328 -m 1 -n 1.55 -p 40 -r 0.525\n\n");

fprintf(stderr,"Report bugs to scott.prahl@oit.edu\n\n");
exit(0);
}/*:12*/
#line 42 "mie_main.w"


int main(int argc,char**argv){

/*4:*/
#line 85 "mie_main.w"


char*g_out_name= NULL;
double pi= 3.14159265358979;

double area;
long i;
struct c_complex m,mm;

struct c_complex*s1= NULL;
struct c_complex*s2= NULL;
double*parallel= NULL;
double*perpen= NULL;
double*phasefn= NULL;
double*mu= NULL;

double x,qext,qsca,qback,g;
int machine_readable_output= 0;
int quiet= 0;

double radius= 0.525;
double lambda_vac= 0.6328;
double lambda_vac_last= 0.6328;
double lambda= 0.6328;
long nangles= 0;
long nlambda= 0;
double density= 1;
double n_medium= 1.0;
m.re= 1.55;
m.im= 0.00;


/*:4*/
#line 46 "mie_main.w"


/*5:*/
#line 119 "mie_main.w"

{
char c;
double xopt;
while((c= my_getopt(argc,argv,"h?qvm:l:L:n:r:i:o:d:p:P:"))!=EOF){
switch(c){

case'r':
sscanf(optarg,"%lf",&xopt);
if(xopt> 0)radius= xopt;
break;

case'm':
sscanf(optarg,"%lf",&xopt);
if(xopt> 0)n_medium= xopt;
break;

case'n':
sscanf(optarg,"%lf",&xopt);
if(xopt> 0)m.re= xopt;
break;

case'l':
sscanf(optarg,"%lf",&xopt);
if(xopt> 0)lambda_vac= xopt;
break;

case'L':
sscanf(optarg,"%lf",&xopt);
if(xopt> 0)lambda_vac_last= xopt;
break;

case'i':
sscanf(optarg,"%lf",&xopt);
if(xopt<=0)m.im= xopt;
break;

case'p':
sscanf(optarg,"%lf",&xopt);
if(xopt>=0)nangles= (long)xopt;
nlambda= 0;
break;

case'P':
sscanf(optarg,"%lf",&xopt);
if(xopt>=0)nlambda= (long)xopt;
nangles= 0;
break;

case'd':
sscanf(optarg,"%lf",&xopt);
if(xopt>=0)density= xopt;
break;

case'o':
g_out_name= strdup(optarg);
break;

case'q':
machine_readable_output= 1;
quiet= 1;
break;

case'v':
print_version();
break;

default:
case'h':
case'?':
print_usage();
break;
}
}

argc-= optind;
argv+= optind;

if(argc> 0){
fprintf(stderr,"No file support in this version.  Sorry.\n");
exit(1);
}

if(g_out_name!=NULL){
if(freopen(g_out_name,"w",stdout)==NULL){
fprintf(stderr,"Could not open file <%s> for output",g_out_name);
exit(1);
}
}
}

/*:5*/
#line 48 "mie_main.w"


/*6:*/
#line 210 "mie_main.w"

if(nangles> 0){
mu= new_darray(nangles);
for(i= 0;i<nangles;i++)
mu[i]= cos(2*pi/nangles*i);

parallel= new_darray(nangles);
perpen= new_darray(nangles);
phasefn= new_darray(nangles);
s1= new_carray(nangles);
s2= new_carray(nangles);
}

/*:6*/
#line 50 "mie_main.w"

/*8:*/
#line 236 "mie_main.w"

printf("# Mie Scattering                # Version %s\n",Version);
printf("# Oregon Medical Laser Center   # http://omlc.org\n");
printf("# by Scott Prahl                # scott.prahl@oit.edu\n");
printf("#\n");

/*:8*/
#line 51 "mie_main.w"


mm.re= m.re/n_medium;
mm.im= m.im/n_medium;
lambda= lambda_vac/n_medium;
x= 2*3.1415926*radius*n_medium/lambda_vac;

Mie(x,mm,mu,nangles,s1,s2,&qext,&qsca,&qback,&g);

/*9:*/
#line 242 "mie_main.w"

{
double mut,mus,musp;

area= 3.14159265358979*radius*radius;
mut= density*qext*area*1000;
mus= density*qsca*area*1000;
musp= mus*(1.0-g);

printf("# radius     \t%9.5f\t [microns]     (sphere radius)\n",radius);
printf("# n_medium   \t%9.5f\t [---]         (refractive index of medium)\n",n_medium);
printf("# n_real     \t%9.5f\t [---]         (refractive index of sphere)\n",m.re);
printf("# n_imag     \t%9.5f\t [---]         (absorption of sphere)\n",m.im);
printf("# lambda_vac \t%9.5f\t [microns]     (wavelength in vacuum)\n",lambda_vac);
printf("# density    \t%9.5f\t [#/micron^3]  (spheres per cubic micron)\n",density);
printf("#\n");
printf("# lambda     \t%9.5f\t [microns]     (wavelength in medium)\n",lambda);
printf("# X          \t%9.5f\t [---]         (size parameter)\n",x);
printf("# g          \t%9.5g\t [---]         (average cosine of phase function)\n",g);
printf("# Qsca       \t%9.5g\t [---]         (scattering efficiency)\n",qsca);
printf("# Qext       \t%9.5g\t [---]         (extinction efficiency)\n",qext);
printf("# Qback      \t%9.5g\t [---]         (backscattering efficiency)\n",qback);
printf("# Csca       \t%9.5g\t [micron^2]    (scattering cross section)\n",qsca*area);
printf("# Cext       \t%9.5g\t [micron^2]    (extinction cross section)\n",qext*area);
printf("# Cback      \t%9.5g\t [micron^2]    (backscattering cross section)\n",qback*area);
printf("# mu_s       \t%9.5g\t [1/mm]        (scattering coefficient)\n",mus);
printf("# mu_s'      \t%9.5g\t [1/mm]        (reduced scattering coefficient)\n",musp);
printf("# mu_t       \t%9.5g\t [1/mm]        (total attenuation coefficient)\n",mut);
}

/*:9*/
#line 60 "mie_main.w"

/*10:*/
#line 272 "mie_main.w"

if(nangles> 0){
int j;
double max_natural,max_perpen,max_parallel;

for(i= 0;i<nangles;++i){
parallel[i]= c_norm(s2[i])/(x*x*qsca)/3.14159;
perpen[i]= c_norm(s1[i])/(x*x*qsca)/3.14159;
phasefn[i]= (parallel[i]+perpen[i])/2.0;
}

max_natural= phasefn[0];
max_perpen= perpen[0];
max_parallel= parallel[0];

for(i= 0;i<nangles;++i){
if(phasefn[i]> max_natural)max_natural= phasefn[i];
if(parallel[i]> max_parallel)max_parallel= parallel[i];
if(perpen[i]> max_perpen)max_perpen= perpen[i];
}

printf("#\n");
printf("#   The second column is normalized so that the integral of it over \n");
printf("#   4*pi steradians will be unity.  The average of the 3rd & 4th\n");
printf("#   columns is the second.  The next three columns are normalized\n");
printf("#   to the value at 0 degrees. \n");
printf("#\n");
printf("#       natural      = (|S1|^2+|S2|^2)/2*1/(pi X^2 Qsca)\n");
printf("#       perpen       = |S1|^2/(pi X^2 Qsca)\n");
printf("#       parallel     = |S2|^2/(pi X^2 Qsca)\n");
printf("#       polarization = (|S1|^2-|S2|^2)/(|S1|^2+|S2|^2)\n");
printf("#       S33          =  Real(S2 * S1^*)\n");
printf("#       S34          = -Imag(S2 * S1^*)\n");
printf("#\n");

printf("##theta\t natural \t perpen \t parallel");
printf("\t natural \t perpen \t parallel\t polarization \t   S33   \t   S34\n");
for(j= 0;j<nangles;j++){
double angle,d,polar,s33,s34;
struct c_complex t;

i= j+nangles/2+1;
if(i>=nangles)i-= nangles;

if(i<=nangles/2)
angle= 180.0/3.1415926*acos(mu[i]);
else
angle= -180.0/3.1415926*acos(mu[i]);

t= c_mul(s2[i],c_conj(s1[i]));
d= (c_norm(s1[i])+c_norm(s2[i]))/2.0;
polar= (c_norm(s1[i])-c_norm(s2[i]))/2.0/d;
s33= (t.re)/d;
s34= -(t.im)/d;
printf("%6.3f\t%8.5f\t%8.5f\t%8.5f",
angle,phasefn[i],perpen[i],parallel[i]);
printf("\t% 6.5f\t% 6.5f\t% 6.5f",
phasefn[i]/max_natural,
perpen[i]/max_perpen,
parallel[i]/max_parallel);
printf("\t% 8.5f\t% 8.5f\t% 8.5f\n",polar,s33,s34);
}
}

/*:10*/
#line 61 "mie_main.w"

/*7:*/
#line 223 "mie_main.w"

if(nangles> 0){
free_darray(mu);
free_darray(parallel);
free_darray(perpen);
free_darray(phasefn);
free_carray(s1);
free_carray(s2);
}


/*:7*/
#line 62 "mie_main.w"


return 0;
}


/*:2*/
