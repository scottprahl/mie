
#define MIE_VERBOSE_ERROR_REPORTING 0

struct c_complex 
Lentz_Dn(struct c_complex z, long n);

void 
Dn_down(struct c_complex z, long nstop, struct c_complex *D);

void 
Dn_up(struct c_complex z, long nstop, struct c_complex *D);

void 
small_Mie(double x, struct c_complex m, double *mu,
	  long nangles, struct c_complex *s1,
	  struct c_complex *s2, double *qext, double *qsca,
	  double *qback, double *g);

void 
small_conducting_Mie(double x, struct c_complex m, double *mu,
		     long nangles, struct c_complex *s1,
		     struct c_complex *s2, double *qext, double *qsca,
		     double *qback, double *g);

void 
Mie(double x, struct c_complex m, double *mu, long nangles, struct c_complex *s1,
struct c_complex *s2, double *qext, double *qsca, double *qback, double *g);

void 
ez_Mie(double x, double n, double *qsca, double *g);

void 
ez_Mie_Full(double x, double m_real, double m_imag, long nangles, double *mu,
	 double *s1_real, double *s1_imag, double *s2_real, double *s2_imag,
	    double *qext, double *qsca, double *qback, double *g);
