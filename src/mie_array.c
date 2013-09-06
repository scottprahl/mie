/*2:*/
#line 9 "./mie_array.w"

#include <stdlib.h> 
#include <string.h> 
#include <stdio.h> 
#include <float.h> 
#include "mie_array.h"

/*6:*/
#line 42 "./mie_array.w"

/*5:*/
#line 39 "./mie_array.w"

static void array_error(char*s)

/*:5*/
#line 43 "./mie_array.w"

{
printf("Array -- %s\n",s);
exit(1);
}

/*:6*/
#line 16 "./mie_array.w"

/*8:*/
#line 54 "./mie_array.w"

/*7:*/
#line 51 "./mie_array.w"

double*new_darray(long size)

/*:7*/
#line 55 "./mie_array.w"

{
double*a;

if(size<=0)array_error("Non-positive double array size chosen");

a= (double*)calloc(sizeof(double),(unsigned long)size+2);
if(a==NULL)array_error("Insufficient space to allocate array");

a[0]= DBL_MIN;
a[size+1]= DBL_MAX;
return a+1;
}

/*:8*/
#line 17 "./mie_array.w"

/*10:*/
#line 72 "./mie_array.w"

/*9:*/
#line 69 "./mie_array.w"

void free_darray(double*a)

/*:9*/
#line 73 "./mie_array.w"

{
if(a!=NULL)free(a-1);
}

/*:10*/
#line 18 "./mie_array.w"

/*12:*/
#line 84 "./mie_array.w"

/*11:*/
#line 81 "./mie_array.w"

double*copy_darray(double*a,long size)

/*:11*/
#line 85 "./mie_array.w"

{
double*b= NULL;
if(a==NULL)return b;
b= new_darray(size+2);
if(b==NULL)array_error("Insufficient space to duplicate array");

memcpy(b,a-1,sizeof(double)*(size+2));
return b+1;
}

/*:12*/
#line 19 "./mie_array.w"

/*14:*/
#line 101 "./mie_array.w"

/*13:*/
#line 98 "./mie_array.w"

void set_darray(double*a,long size,double x)

/*:13*/
#line 102 "./mie_array.w"

{
long j;

if(a==NULL)array_error("Attempt to set elements in a NULL array");

for(j= 0;j<size;j++)
a[j]= x;
}

/*:14*/
#line 20 "./mie_array.w"

/*16:*/
#line 117 "./mie_array.w"

/*15:*/
#line 114 "./mie_array.w"

void min_max_darray(double*a,long size,double*min,double*max)

/*:15*/
#line 118 "./mie_array.w"

{
long j;

if(a==NULL)array_error("A NULL array does not have a min or max");

if(size==0)array_error("An array with no elements does not have a min or max");

*min= a[0];
*max= *min;
for(j= 1;j<size;j++){
if(a[j]> *max)*max= a[j];
if(a[j]<*min)*min= a[j];
}

}

/*:16*/
#line 21 "./mie_array.w"

/*19:*/
#line 145 "./mie_array.w"

/*18:*/
#line 142 "./mie_array.w"

void sort_darray(double*a,long size)

/*:18*/
#line 146 "./mie_array.w"

{
long i,ir,j,l;
double aa;

if(a==NULL)array_error("Can't sort a NULL array");

if(size<2)return;
l= (size>>1)+1;
ir= size;
for(;;){
if(l> 1){
aa= a[--l-1];
}else{
aa= a[ir-1];
a[ir-1]= a[0];
if(--ir==1){
a[0]= aa;
break;
}
}
i= l;
j= l+l;
while(j<=ir){
if(j<ir&&a[j-1]<a[j])j++;
if(aa<a[j-1]){
a[i-1]= a[j-1];
i= j;
j<<= 1;
}else j= ir+1;
}
a[i-1]= aa;
}
}

/*:19*/
#line 22 "./mie_array.w"

/*22:*/
#line 189 "./mie_array.w"

/*21:*/
#line 186 "./mie_array.w"

void print_darray(double*a,long size,long ilow,long ihigh)

/*:21*/
#line 190 "./mie_array.w"

{
long j;

if(a==NULL)array_error("Can't print a NULL array");

if(ilow<0)ilow= 0;
if(ihigh> size-1)ihigh= size-1;

for(j= ilow;j<=ihigh;j++)
printf("x[%ld]= %-10.5g \n",j,a[j]);

}

/*:22*/
#line 23 "./mie_array.w"


/*:2*/
