/*3:*/
#line 46 "./mie.w"

#define MIE_VERBOSE_ERROR_REPORTING 0
/*10:*/
#line 122 "./mie.w"

struct c_complex Lentz_Dn(struct c_complex z,long n)

/*:10*/
#line 48 "./mie.w"
;
/*19:*/
#line 243 "./mie.w"

void Dn_down(struct c_complex z,long nstop,struct c_complex*D)

/*:19*/
#line 49 "./mie.w"
;
/*16:*/
#line 212 "./mie.w"

void Dn_up(struct c_complex z,long nstop,struct c_complex*D)

/*:16*/
#line 50 "./mie.w"
;
/*22:*/
#line 278 "./mie.w"

void small_Mie(double x,struct c_complex m,double*mu,
long nangles,struct c_complex*s1,
struct c_complex*s2,double*qext,double*qsca,
double*qback,double*g)

/*:22*/
#line 51 "./mie.w"
;
/*31:*/
#line 441 "./mie.w"

void small_conducting_Mie(double x,struct c_complex m,double*mu,
long nangles,struct c_complex*s1,
struct c_complex*s2,double*qext,double*qsca,
double*qback,double*g)

/*:31*/
#line 52 "./mie.w"
;
/*34:*/
#line 514 "./mie.w"

void Mie(double x,struct c_complex m,double*mu,long nangles,struct c_complex*s1,
struct c_complex*s2,double*qext,double*qsca,double*qback,double*g)

/*:34*/
#line 53 "./mie.w"
;
/*50:*/
#line 842 "./mie.w"

void ez_Mie(double x,double n,double*qsca,double*g)

/*:50*/
#line 54 "./mie.w"
;
/*52:*/
#line 875 "./mie.w"

void ez_Mie_Full(double x,double m_real,double m_imag,long nangles,double*mu,
double*s1_real,double*s1_imag,double*s2_real,double*s2_imag,
double*qext,double*qsca,double*qback,double*g)

/*:52*/
#line 55 "./mie.w"
;

/*:3*/
