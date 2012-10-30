
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include "mie_array.h"



static void 
array_error(char *s)
{
    printf("Array -- %s\n", s);
    exit(1);
}




double             *
new_darray(long size)
{
    double             *a;

    if (size <= 0)
	array_error("Non-positive double array size chosen");

    a = (double *) calloc(sizeof(double), (unsigned long) size + 2);
    if (a == NULL)
	array_error("Insufficient space to allocate array");

    a[0] = DBL_MIN;
    a[size + 1] = DBL_MAX;
    return a + 1;
}




void 
free_darray(double *a)
{
    if (a != NULL)
	free(a - 1);
}




double             *
copy_darray(double *a, long size)
{
    double             *b = NULL;

    if (a == NULL)
	return b;
    b = new_darray(size + 2);
    if (b == NULL)
	array_error("Insufficient space to duplicate array");

    memcpy(b, a - 1, sizeof(double) * (size + 2));
    return b + 1;
}




void 
set_darray(double *a, long size, double x)
{
    long                j;

    if (a == NULL)
	array_error("Attempt to set elements in a NULL array");

    for (j = 0; j < size; j++)
	a[j] = x;
}




void 
min_max_darray(double *a, long size, double *min, double *max)
{
    long                j;

    if (a == NULL)
	array_error("A NULL array does not have a min or max");

    if (size == 0)
	array_error("An array with no elements does not have a min or max");

    *min = a[0];
    *max = *min;
    for (j = 1; j < size; j++) {
	if (a[j] > *max)
	    *max = a[j];
	if (a[j] < *min)
	    *min = a[j];
    }

}




void 
sort_darray(double *a, long size)
{
    long                i, ir, j, l;
    double              aa;

    if (a == NULL)
	array_error("Can't sort a NULL array");

    if (size < 2)
	return;
    l = (size >> 1) + 1;
    ir = size;
    for (;;) {
	if (l > 1) {
	    aa = a[--l - 1];
	} else {
	    aa = a[ir - 1];
	    a[ir - 1] = a[0];
	    if (--ir == 1) {
		a[0] = aa;
		break;
	    }
	}
	i = l;
	j = l + l;
	while (j <= ir) {
	    if (j < ir && a[j - 1] < a[j])
		j++;
	    if (aa < a[j - 1]) {
		a[i - 1] = a[j - 1];
		i = j;
		j <<= 1;
	    } else
		j = ir + 1;
	}
	a[i - 1] = aa;
    }
}




void 
print_darray(double *a, long size, long ilow, long ihigh)
{
    long                j;

    if (a == NULL)
	array_error("Can't print a NULL array");

    if (ilow < 0)
	ilow = 0;
    if (ihigh > size - 1)
	ihigh = size - 1;

    for (j = ilow; j <= ihigh; j++)
	printf("x[%ld]= %-10.5g \n", j, a[j]);

}
