@** A driver program for Mie scattering.

This should calculate the scattering and anisotropy coefficients
for 360\,nm latex microspheres.

@(test_latex.c@>=

#include <stdio.h>
#include <math.h>
#include "mie_array.h"
#include "mie_complex.h"
#include "mie_legendre.h"
#include "mie_lobatto.h"
#include "mie.h"

int main ()
{
	@<Declare Mie Variables@>@;
	@<Print header@>@;

	for(lambda=300;lambda<=400;lambda+=5){	
		@<Calculate sphere size and number@>@;
        
		ez_Mie(x, n_sphere/n_medium, &qsca, &g);
		@<Print summary@>@;
	}
	return 0;
}


@ Declare and initialize many variables.

@<Declare Mie Variables@>=

double pi=3.14159265358979;

double n_sphere = 1.59;
double n_medium = 1.34;
double diameter= 360;          /* in nm */
double sphere_density = 1.05;    /* in gm/cc */
double medium_density = 1.00;     /* in gm/cc */
double concentration = 0.01;     /* in gm spheres/ gm fluid */

double x, qsca, g;
double sphere_area, sphere_volume;
double number_per_cc;
double lambda;

  
@ Calculating the number of spheres per cubic centimeter is straightforward, if onerous.
The first constraint is that 
$$
1\hbox{ cc} = N_s V_s + V_f
$$
where the subscript $s$ refers to the sphere and $f$ to the surrounding fluid.
$N_s$ is the number of spheres per milliliter, and $V$ is the volume.
The next equation is that the concentration $c$ (in gm spheres/gm fluid) is 
$$
c = {N_s w_s\over w_f} = {N_s \rho_s V_s\over \rho_f V_f}
$$
where $w$ is the weight and $\rho$ is the density.  Putting these two equations
together leads to
$$
\hbox{1 cc} = N_s V_s{\rho_f c +\rho_s \over \rho_f c}
$$
or
$$
N_s = {\hbox{1 cc}  \over V_s}{1\over 1+ \rho_s/(\rho_f c)}
$$

@<Calculate sphere size and number@>=
  x = (diameter * pi) / (lambda / (n_sphere/n_medium));
  sphere_area   = pi*diameter*diameter/4.0/1e14; /* in cm$^2$ */
  sphere_volume = pi*diameter*diameter*diameter/6.0/1e21;  /* in cc */
  number_per_cc = 1/sphere_volume/(1 + sphere_density/medium_density/concentration);
  
@ Print a header then the angles.  Make sure everything lines up.

@<Print header@>=
  printf("Mie Scattering -- Non-absorbing spheres\n");
  printf("Sphere concentration is          %5.2g%% gm spheres/gm medium\n", concentration);
  printf("Medium refractive index is       %7.4g \n", n_medium);
  printf("Sphere refractive index is       %7.4g \n", n_sphere);
  printf("The density of the spheres is    %10g gm/ml\n", sphere_density);
  printf("The density of the medium is     %10g gm/ml\n", medium_density);
  printf("\n");

  printf("  d     lambda    x       N      Qsca      g     mu_s  mu_s(1-g)\n");
  printf(" [nm]    [nm]            #/ml                    [1/cm]  [1/cm] \n");

@ @<Print summary@>=
  printf("%7.1f %7.1f ", diameter, lambda);
  printf("%7.2f %7.2e %7.4f %7.5f ", x, number_per_cc, qsca, g);
  printf("%7.0f ", qsca * sphere_area * number_per_cc);
  printf("%7.1f \n", (1 - g) * qsca * sphere_area * number_per_cc);

