/*22:*/
#line 269 "./mie_lobatto.w"


#include <stddef.h> 
#include <stdio.h> 
#include "mie_array.h"
#include "mie_lobatto.h"

int main()
{
double*x,*w;
double sum;
long i,n;

printf("testing n=10 --- sum should be 2\n");
sum= 0;
n= 10;
x= new_darray(n);
w= new_darray(n);
Lobatto(-1.0,1.0,x,w,n);
printf("The x_i are\n");
print_darray(x,n,0L,n-1);
printf("The w_i are\n");
print_darray(w,n,0L,n-1);
for(i= 0;i<n;i++)sum= sum+w[i];
printf("sum   %20.15f\n",sum);
free_darray(x);
free_darray(w);
printf("\n");

printf("testing n=9 --- sum should be 2\n");
sum= 0;
n= 9;
x= new_darray(n);
w= new_darray(n);
Lobatto(-1.0,1.0,x,w,n);
printf("The x_i are\n");
print_darray(x,n,0L,n-1);
printf("The w_i are\n");
print_darray(w,n,0L,n-1);
for(i= 0;i<n;i++)sum= sum+w[i];
printf("sum   %20.15f\n",sum);
free_darray(x);
free_darray(w);
printf("\n");
return 0;
}

/*:22*/
