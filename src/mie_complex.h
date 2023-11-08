
struct c_complex {
    double		re;
    double		im;
};


struct c_complex 
c_set(double a, double b);

struct c_complex 
c_polarset(double r, double theta);


double 
c_abs(struct c_complex z);

double 
c_arg(struct c_complex z);

struct c_complex 
c_sqr(struct c_complex z);

struct c_complex 
c_conj(struct c_complex z);

double 
c_norm(struct c_complex z);

struct c_complex 
c_sqrt(struct c_complex z);

struct c_complex 
c_inv(struct c_complex w);


struct c_complex 
c_add(struct c_complex z, struct c_complex w);

struct c_complex 
c_sub(struct c_complex z, struct c_complex w);

struct c_complex 
c_mul(struct c_complex z, struct c_complex w);

struct c_complex 
c_div(struct c_complex z, struct c_complex w);


double 
c_rdiv(struct c_complex z, struct c_complex w);

double 
c_rmul(struct c_complex z, struct c_complex w);

struct c_complex 
c_sadd(double x, struct c_complex z);

struct c_complex 
c_sdiv(double x, struct c_complex w);

struct c_complex 
c_smul(double x, struct c_complex z);


struct c_complex 
c_sin(struct c_complex z);

struct c_complex 
c_cos(struct c_complex z);

struct c_complex 
c_tan(struct c_complex z);

struct c_complex 
c_asin(struct c_complex z);

struct c_complex 
c_acos(struct c_complex z);

struct c_complex 
c_atan(struct c_complex z);


struct c_complex 
c_sinh(struct c_complex z);

struct c_complex 
c_cosh(struct c_complex z);

struct c_complex 
c_tanh(struct c_complex z);

struct c_complex 
c_atanh(struct c_complex z);

struct c_complex 
c_asinh(struct c_complex z);


struct c_complex 
c_exp(struct c_complex z);

struct c_complex 
c_log(struct c_complex z);

struct c_complex 
c_log10(struct c_complex z);


struct c_complex   *
new_carray(long size);

void 
free_carray(struct c_complex *a);

struct c_complex   *
copy_carray(struct c_complex *a, long size);

void 
set_carray(struct c_complex *a, long size, struct c_complex z);
