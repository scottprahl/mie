

#include <stdio.h>
#include <math.h>
#include "mie_array.h"
#include "mie_complex.h"
#include "mie_legendre.h"
#include "mie_lobatto.h"
#include "mie.h"

int 
main()
{


    double		pi = 3.14159265358979;

    double		n_sphere = 1.59;
    double		n_medium = 1.34;
    double		diameter = 360;
    double		sphere_density = 1.05;
    double		medium_density = 1.00;
    double		concentration = 0.01;

    double		x          , qsca, g;
    double		sphere_area, sphere_volume;
    double		number_per_cc;
    double		lambda;




    printf("Mie Scattering -- Non-absorbing spheres\n");
    printf("Sphere concentration is          %5.2g%% gm spheres/gm medium\n", concentration);
    printf("Medium refractive index is       %7.4g \n", n_medium);
    printf("Sphere refractive index is       %7.4g \n", n_sphere);
    printf("The density of the spheres is    %10g gm/ml\n", sphere_density);
    printf("The density of the medium is     %10g gm/ml\n", medium_density);
    printf("\n");

    printf("  d     lambda    x       N      Qsca      g     mu_s  mu_s(1-g)\n");
    printf(" [nm]    [nm]            #/ml                    [1/cm]  [1/cm] \n");



    for (lambda = 300; lambda <= 400; lambda += 5) {

	x = (diameter * pi) / (lambda / (n_sphere / n_medium));
	sphere_area = pi * diameter * diameter / 4.0 / 1e14;
	sphere_volume = pi * diameter * diameter * diameter / 6.0 / 1e21;
	number_per_cc = 1 / sphere_volume / (1 + sphere_density / medium_density / concentration);



	ez_Mie(x, n_sphere / n_medium, &qsca, &g);

	printf("%7.1f %7.1f ", diameter, lambda);
	printf("%7.2f %7.2e %7.4f %7.5f ", x, number_per_cc, qsca, g);
	printf("%7.0f ", qsca * sphere_area * number_per_cc);
	printf("%7.1f \n", (1 - g) * qsca * sphere_area * number_per_cc);

    }
    return 0;
}
