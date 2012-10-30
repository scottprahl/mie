@** Double Array Routines.

Here are a bunch of routines to deal arrays of doubles.  This file will
create three files when run through {\tt ctangle} --- the usual {\tt .c}
and {\tt .h}, as well as a testing driver. 

@ Here, then, is an overview of document structure

@(mie_array.c@>=
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include "mie_array.h"

@<Definition for |array_error|@>@;
@<Definition for |new_darray|@>@;
@<Definition for |free_darray|@>@;
@<Definition for |copy_darray|@>@;
@<Definition for |set_darray|@>@;
@<Definition for |min_max_darray|@>@;
@<Definition for |sort_darray|@>@;
@<Definition for |print_darray|@>@;

@ Each function has its prototype exported to a header file.

@(mie_array.h@>=
@<Prototype for |new_darray|@>;
@<Prototype for |free_darray|@>;
@<Prototype for |copy_darray|@>;
@<Prototype for |set_darray|@>;
@<Prototype for |min_max_darray|@>;
@<Prototype for |sort_darray|@>;
@<Prototype for |print_darray|@>;

@*1 Allocation.

@ A simple error routine.
@<Prototype for |array_error|@>=
static void array_error(char *s)

@ @<Definition for |array_error|@>=
		@<Prototype for |array_error|@>
{
 printf("Array -- %s\n", s);
 exit(1);
 }

@*1 Double array routines.

@<Prototype for |new_darray|@>=
double * new_darray(long size)

@ @<Definition for |new_darray|@>=
		@<Prototype for |new_darray|@>
{
  double *a;

  if (size<=0) array_error("Non-positive double array size chosen");

  a = (double *) calloc(sizeof(double), (unsigned long) size+2);
  if (a==NULL) array_error("Insufficient space to allocate array");
  
  a[0]=DBL_MIN;
  a[size+1] = DBL_MAX;
  return a+1;
}

@ @<Prototype for |free_darray|@>=
void free_darray(double *a)

@ @<Definition for |free_darray|@>=
		@<Prototype for |free_darray|@>
{
 if (a!=NULL) free(a-1);
}

@ This allocates a new double array data structure and 
copies the contents of |a| into it.

@<Prototype for |copy_darray|@>=
double *copy_darray(double *a, long size)

@ @<Definition for |copy_darray|@>=
		@<Prototype for |copy_darray|@>
{
  double *b=NULL;
  if (a==NULL) return b;
  b = new_darray(size+2);
  if (b==NULL) array_error("Insufficient space to duplicate array");

  memcpy(b, a-1, sizeof(double)*(size+2));
  return b+1;
}

@ This sets all the entries in the array |a[]| to |x|.

@<Prototype for |set_darray|@>=
void set_darray(double * a, long size, double x)

@ @<Definition for |set_darray|@>=
		@<Prototype for |set_darray|@>
{
  long j;

  if (a==NULL) array_error("Attempt to set elements in a NULL array");
    
  for (j=0; j< size; j++)
  	a[j] = x;
}

@ |min_max_darray| finds the minimum and maximum of the array |a|.

@<Prototype for |min_max_darray|@>=
void min_max_darray(double * a, long size, double *min, double *max)

@ @<Definition for |min_max_darray|@>=
		@<Prototype for |min_max_darray|@>
{
  long j;

  if (a==NULL) array_error("A NULL array does not have a min or max");
  
  if (size==0) array_error("An array with no elements does not have a min or max");
   
  *min=a[0];
  *max=*min;
  for (j=1; j< size; j++) {
  	if (a[j] > *max) *max = a[j];
  	if (a[j] < *min) *min = a[j];
  }
  
}

@*1 Sorting.

@ |sort_darray| will sort an array |a| into ascending numerical
order using the Heapsort algorithm.  Adapted to work with zero-based arrays
from {\it Numerical Recipes}.  This could certainly use some sprucing up,
but I can't quite seem to figure out how to do it.  It is kinda tricky.

@<Prototype for |sort_darray|@>=
void sort_darray(double *a, long size)

@ @<Definition for |sort_darray|@>=
		@<Prototype for |sort_darray|@>
{
	long i,ir,j,l;
	double aa;

    if (a==NULL) array_error("Can't sort a NULL array");

	if (size < 2) return;
	l=(size >> 1)+1;
	ir=size;
	for (;;) {
		if (l > 1) {
			aa=a[--l-1];
		} else {
			aa=a[ir-1];
			a[ir-1]=a[0];
			if (--ir == 1) {
				a[0]=aa;
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			if (j < ir && a[j-1] < a[j]) j++;
			if (aa < a[j-1]) {
				a[i-1]=a[j-1];
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		a[i-1]=aa;
	}
}

@*1 Printing.

@ |print_darray| prints the elements of the array |a| from
|ilow| through |ihigh|.

@<Prototype for |print_darray|@>=
void print_darray(double * a, long size, long ilow, long ihigh)

@ @<Definition for |print_darray|@>=
		@<Prototype for |print_darray|@>
{
  long j;

  if (a==NULL)  array_error("Can't print a NULL array");
  
  if (ilow<0) ilow=0;
  if (ihigh>size-1) ihigh=size-1;
  
  for (j=ilow; j<= ihigh; j++)
  	printf("x[%ld]= %-10.5g \n", j, a[j]);
  
}

@*1 Testing.
@ Here are driver routines to test the routines in this file.

@(test_mie_array.c@>=

#include <stdio.h>
#include "mie_array.h"

void main() {

double *x;
double *y;
long   i,size;
double min, max;

size=10;
printf("starting\n");
fflush(stdout);

x = new_darray(size);

@<Test Set Routine@>@;

@<Test Copy Routine@>@;

@<Test Sort Routine@>@;

@<Test Min/Max Routine@>@;

printf("done\n");
fflush(stdout);
}


@ @<Test Set Routine@>=
printf("Testing set_darray\n");
printf("All entries should be 3.0\n");

set_darray(x,size,3.0);
print_darray(x,size,0,size-1);
fflush(stdout);
printf("\n");

@ @<Test Copy Routine@>=
printf("Testing copy_darray\n");

for(i=0; i<size; i++) x[i] = size-i;
printf("The original vector was:\n");
print_darray(x,size,0,size-1);
fflush(stdout);

y = copy_darray(x,size);
printf("The copied vector is:\n");
print_darray(y,size,0,size-1);
fflush(stdout);
printf("\n");

@ @<Test Sort Routine@>=
printf("Testing sort_darray\n");
printf("The original vector is:\n");
print_darray(x,size,0,size-1);
fflush(stdout);

sort_darray(x,size);
printf("The sorted vector is:\n");
print_darray(x,size,0,size-1);
fflush(stdout);
printf("\n");

@ @<Test Min/Max Routine@>=

printf("Testing min_max_darray\n");

min_max_darray(x,size,&min,&max);
printf("min=%g  max=%g\n", min, max);
fflush(stdout);
printf("\n");
