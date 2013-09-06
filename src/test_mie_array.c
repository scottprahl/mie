/*24:*/
#line 207 "./mie_array.w"


#include <stdio.h> 
#include "mie_array.h"

void main(){

double*x;
double*y;
long i,size;
double min,max;

size= 10;
printf("starting\n");
fflush(stdout);

x= new_darray(size);

/*25:*/
#line 238 "./mie_array.w"

printf("Testing set_darray\n");
printf("All entries should be 3.0\n");

set_darray(x,size,3.0);
print_darray(x,size,0,size-1);
fflush(stdout);
printf("\n");

/*:25*/
#line 225 "./mie_array.w"


/*26:*/
#line 247 "./mie_array.w"

printf("Testing copy_darray\n");

for(i= 0;i<size;i++)x[i]= size-i;
printf("The original vector was:\n");
print_darray(x,size,0,size-1);
fflush(stdout);

y= copy_darray(x,size);
printf("The copied vector is:\n");
print_darray(y,size,0,size-1);
fflush(stdout);
printf("\n");

/*:26*/
#line 227 "./mie_array.w"


/*27:*/
#line 261 "./mie_array.w"

printf("Testing sort_darray\n");
printf("The original vector is:\n");
print_darray(x,size,0,size-1);
fflush(stdout);

sort_darray(x,size);
printf("The sorted vector is:\n");
print_darray(x,size,0,size-1);
fflush(stdout);
printf("\n");

/*:27*/
#line 229 "./mie_array.w"


/*28:*/
#line 273 "./mie_array.w"


printf("Testing min_max_darray\n");

min_max_darray(x,size,&min,&max);
printf("min=%g  max=%g\n",min,max);
fflush(stdout);
printf("\n");/*:28*/
#line 231 "./mie_array.w"


printf("done\n");
fflush(stdout);
}


/*:24*/
