/*3:*/
#line 67 "./mie_complex.w"

struct c_complex{double re;double im;};

/*7:*/
#line 130 "./mie_complex.w"

struct c_complex c_set(double a,double b)

/*:7*/
#line 70 "./mie_complex.w"
;
/*9:*/
#line 145 "./mie_complex.w"

struct c_complex c_polarset(double r,double theta)

/*:9*/
#line 71 "./mie_complex.w"
;

/*11:*/
#line 158 "./mie_complex.w"

double c_abs(struct c_complex z)

/*:11*/
#line 73 "./mie_complex.w"
;
/*15:*/
#line 191 "./mie_complex.w"

double c_arg(struct c_complex z)

/*:15*/
#line 74 "./mie_complex.w"
;
/*21:*/
#line 236 "./mie_complex.w"

struct c_complex c_sqr(struct c_complex z)

/*:21*/
#line 75 "./mie_complex.w"
;
/*13:*/
#line 182 "./mie_complex.w"

struct c_complex c_conj(struct c_complex z)

/*:13*/
#line 76 "./mie_complex.w"
;
/*17:*/
#line 202 "./mie_complex.w"

double c_norm(struct c_complex z)

/*:17*/
#line 77 "./mie_complex.w"
;
/*19:*/
#line 211 "./mie_complex.w"

struct c_complex c_sqrt(struct c_complex z)

/*:19*/
#line 78 "./mie_complex.w"
;
/*23:*/
#line 247 "./mie_complex.w"

struct c_complex c_inv(struct c_complex w)

/*:23*/
#line 79 "./mie_complex.w"
;

/*26:*/
#line 272 "./mie_complex.w"

struct c_complex c_add(struct c_complex z,struct c_complex w)

/*:26*/
#line 81 "./mie_complex.w"
;
/*28:*/
#line 287 "./mie_complex.w"

struct c_complex c_sub(struct c_complex z,struct c_complex w)

/*:28*/
#line 82 "./mie_complex.w"
;
/*30:*/
#line 301 "./mie_complex.w"

struct c_complex c_mul(struct c_complex z,struct c_complex w)

/*:30*/
#line 83 "./mie_complex.w"
;
/*32:*/
#line 316 "./mie_complex.w"

struct c_complex c_div(struct c_complex z,struct c_complex w)

/*:32*/
#line 84 "./mie_complex.w"
;

/*34:*/
#line 344 "./mie_complex.w"

double c_rdiv(struct c_complex z,struct c_complex w)

/*:34*/
#line 86 "./mie_complex.w"
;
/*36:*/
#line 368 "./mie_complex.w"

double c_rmul(struct c_complex z,struct c_complex w)

/*:36*/
#line 87 "./mie_complex.w"
;
/*41:*/
#line 395 "./mie_complex.w"

struct c_complex c_sadd(double x,struct c_complex z)

/*:41*/
#line 88 "./mie_complex.w"
;
/*43:*/
#line 410 "./mie_complex.w"

struct c_complex c_sdiv(double x,struct c_complex w)

/*:43*/
#line 89 "./mie_complex.w"
;
/*39:*/
#line 381 "./mie_complex.w"

struct c_complex c_smul(double x,struct c_complex z)

/*:39*/
#line 90 "./mie_complex.w"
;

/*46:*/
#line 439 "./mie_complex.w"

struct c_complex c_sin(struct c_complex z)

/*:46*/
#line 92 "./mie_complex.w"
;
/*48:*/
#line 450 "./mie_complex.w"

struct c_complex c_cos(struct c_complex z)

/*:48*/
#line 93 "./mie_complex.w"
;
/*50:*/
#line 481 "./mie_complex.w"

struct c_complex c_tan(struct c_complex z)

/*:50*/
#line 94 "./mie_complex.w"
;
/*52:*/
#line 503 "./mie_complex.w"

struct c_complex c_asin(struct c_complex z)

/*:52*/
#line 95 "./mie_complex.w"
;
/*54:*/
#line 516 "./mie_complex.w"

struct c_complex c_acos(struct c_complex z)

/*:54*/
#line 96 "./mie_complex.w"
;
/*56:*/
#line 529 "./mie_complex.w"

struct c_complex c_atan(struct c_complex z)

/*:56*/
#line 97 "./mie_complex.w"
;

/*61:*/
#line 551 "./mie_complex.w"

struct c_complex c_sinh(struct c_complex z)

/*:61*/
#line 99 "./mie_complex.w"
;
/*59:*/
#line 542 "./mie_complex.w"

struct c_complex c_cosh(struct c_complex z)

/*:59*/
#line 100 "./mie_complex.w"
;
/*63:*/
#line 560 "./mie_complex.w"

struct c_complex c_tanh(struct c_complex z)

/*:63*/
#line 101 "./mie_complex.w"
;
/*65:*/
#line 573 "./mie_complex.w"

struct c_complex c_atanh(struct c_complex z)

/*:65*/
#line 102 "./mie_complex.w"
;
/*67:*/
#line 582 "./mie_complex.w"

struct c_complex c_asinh(struct c_complex z)

/*:67*/
#line 103 "./mie_complex.w"
;

/*70:*/
#line 593 "./mie_complex.w"

struct c_complex c_exp(struct c_complex z)

/*:70*/
#line 105 "./mie_complex.w"
;
/*72:*/
#line 603 "./mie_complex.w"

struct c_complex c_log(struct c_complex z)

/*:72*/
#line 106 "./mie_complex.w"
;
/*74:*/
#line 612 "./mie_complex.w"

struct c_complex c_log10(struct c_complex z)

/*:74*/
#line 107 "./mie_complex.w"
;

/*77:*/
#line 625 "./mie_complex.w"

struct c_complex*new_carray(long size)

/*:77*/
#line 109 "./mie_complex.w"
;
/*79:*/
#line 641 "./mie_complex.w"

void free_carray(struct c_complex*a)

/*:79*/
#line 110 "./mie_complex.w"
;
/*81:*/
#line 653 "./mie_complex.w"

struct c_complex*copy_carray(struct c_complex*a,long size)

/*:81*/
#line 111 "./mie_complex.w"
;
/*83:*/
#line 670 "./mie_complex.w"

void set_carray(struct c_complex*a,long size,struct c_complex z)

/*:83*/
#line 112 "./mie_complex.w"
;

/*:3*/
