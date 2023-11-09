#include <stdio.h>
#include <math.h>
#include "mie.h"

int main() {
    const double n_sphere = 1.59;
    const double n_medium = 1.34;
    const double diameter_nm = 360;
    const double sphere_density = 1.05; // gm/ml
    const double medium_density = 1.00; // gm/ml
    const double concentration = 0.01; // gm spheres/gm medium
    double x, qsca, g, sphere_area_cm2, sphere_volume_ml, number_per_cc, lambda_nm;

    printf("Mie Scattering -- Non-absorbing spheres\n\n");
    printf("Sphere concentration is          %5.2g%% gm spheres/gm\n", concentration * 100);
    printf("Medium refractive index is       %7.4f\n", n_medium);
    printf("Sphere refractive index is       %7.4f\n", n_sphere);
    printf("The density of the spheres is    %7g gm/ml\n", sphere_density);
    printf("The density of the medium is     %7g gm/ml\n", medium_density);
    printf("\n");
    printf("    d    lambda     x      N       Q_sca    g      mu_s    mu_s'\n");
    printf("  [nm]    [nm]     [-]   [#/ml]     [-]    [-]    [1/cm]  [1/cm]\n");

    for (lambda_nm = 300; lambda_nm <= 400; lambda_nm += 5) {
        x = (diameter_nm * M_PI) / (lambda_nm / (n_sphere / n_medium));
        sphere_area_cm2 = M_PI * pow(diameter_nm / 1e7, 2) / 4.0;
        sphere_volume_ml = M_PI * pow(diameter_nm / 1e7, 3) / 6.0;
        number_per_cc = 1 / (sphere_volume_ml * (1 + sphere_density / (medium_density * concentration)));
        ez_Mie(x, n_sphere / n_medium, &qsca, &g);

        printf("%7.1f %7.1f ", diameter_nm, lambda_nm);
        printf("%7.2f %7.2e %7.4f %7.5f ", x, number_per_cc, qsca, g);
        printf("%7.2f ", qsca * sphere_area_cm2 * number_per_cc);
        printf("%7.2f \n", (1 - g) * qsca * sphere_area_cm2 * number_per_cc);
    }
    return 0;
}