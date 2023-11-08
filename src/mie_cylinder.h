/*3:*/
#line 29 "mie_cylinder.w"

/*16:*/
#line 269 "mie_cylinder.w"

void MieCylinderCoefficients(double x,struct c_complex m,double zeta,int n_terms,
struct c_complex an1[],
struct c_complex bn1[],
struct c_complex an2[],
struct c_complex bn2[])

/*:16*/
#line 30 "mie_cylinder.w"
;
/*18:*/
#line 346 "mie_cylinder.w"

void MieCylinder(double x,struct c_complex m,double zeta,const double*theta,int nangles,
struct c_complex*t1,struct c_complex*t2,struct c_complex*t3,
struct c_complex*qexpar,struct c_complex*qexper,
double*qscpar,double*qscper)

/*:18*/
#line 31 "mie_cylinder.w"
;

/*:3*/
