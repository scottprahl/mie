

void		    bracketroot(double (*fx) (double), double x1, double x2, long n,
       		    double	    xb1[], double xb2[], long *nrequested);

double		    saferoot(void (*funcd) (double, double *, double *), double x1, double x2, double xacc);
