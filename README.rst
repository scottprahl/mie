Mie Scattering
==============

.. image:: https://img.shields.io/github/v/tag/scottprahl/mie?label=github&color=68CA66
   :target: https://github.com/scottprahl/mie

.. image:: https://img.shields.io/github/license/scottprahl/mie?color=68CA66)
   :target: https://github.com/scottprahl/mie/blob/master/LICENSE

.. image:: https://img.shields.io/badge/docs-passing-68CA66
   :target: https://github.com/scottprahl/mie/blob/master/doc/mie_doc.pdf

.. image:: https://github.com/scottprahl/mie/actions/workflows/test.yaml/badge.svg
   :target: https://github.com/scottprahl/mie/actions

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.10087653.svg
   :target: https://doi.org/10.5281/zenodo.10087653

Mie Scattering is a computational project for simulating the scattering of
electromagnetic waves by small particles. Utilizing the robust MIEV0 FORTRAN
code by Wiscombe, this software offers researchers and scientists an accurate
and reliable tool for Mie scattering calculations.

Features
--------

- Accurate simulation of scattering for small to very large sphere sizes (𝜋d/λ > 10,000).

- Based on the well-tested MIEV0 FORTRAN code.

- Includes additional scattering code for cylinders.

- Complemented by `miepython <https://github.com/scottprahl/miepython>`_ for Python users.

The source code is written in `CWEB <https://github.com/ascherer/cweb>`_, which
allows excellent documentation of scientific programs. Basically, there is a
program ``ctangle`` that converts cweb code to C. There is another program
``cweave`` that converts cweb code to TeX. This then generates really `nice
documentation <https://github.com/scottprahl/mie/blob/master/doc/mie_doc.pdf>`_.

Getting Started
---------------

Downloading
~~~~~~~~~~~

Clone the repository using git or download the source code as a zip file:

.. code-block:: shell

    git clone https://github.com/scottprahl/mie.git

Installation
~~~~~~~~~~~~

Unix/macOS
^^^^^^^^^^

Run the following command to build and install:

.. code-block:: shell

    make install

The executable (``mie``) will be installed in ``/usr/local/bin/``.

Windows
^^^^^^^

Download the latest executable from the `Releases page <https://github.com/scottprahl/mie/releases>`_.

Python Package
^^^^^^^^^^^^^^

I have written a pure python version of the sphere scattering code that is based on this code.  
It is a bit slower, but because it is python, but it has a number of nice affordances that make doing Mie calculations (and plotting them) less of a hassle.

.. code-block:: shell

    pip install miepython

Usage
-----

.. code-block:: c

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

When compiled (see `src/Makefile`) this produces::

    Sphere concentration is              1% gm spheres/gm
    Medium refractive index is        1.3400
    Sphere refractive index is        1.5900
    The density of the spheres is       1.05 gm/ml
    The density of the medium is           1 gm/ml

        d    lambda     x      N       Q_sca    g      mu_s    mu_s'
      [nm]    [nm]     [-]   [#/ml]     [-]    [-]    [1/cm]  [1/cm]
      360.0   300.0    4.47 3.86e+11  1.2843 0.87699  504.82   62.10 
      360.0   305.0    4.40 3.86e+11  1.2453 0.87409  489.49   61.63 
      360.0   310.0    4.33 3.86e+11  1.2075 0.87131  474.65   61.08 
      360.0   315.0    4.26 3.86e+11  1.1711 0.86861  460.32   60.48 
      360.0   320.0    4.19 3.86e+11  1.1360 0.86600  446.54   59.84 
      360.0   325.0    4.13 3.86e+11  1.1024 0.86347  433.34   59.16 
      360.0   330.0    4.07 3.86e+11  1.0703 0.86102  420.70   58.47 
      360.0   335.0    4.01 3.86e+11  1.0395 0.85866  408.62   57.75 
      360.0   340.0    3.95 3.86e+11  1.0101 0.85640  397.04   57.02 
      360.0   345.0    3.89 3.86e+11  0.9818 0.85422  385.92   56.26 
      360.0   350.0    3.83 3.86e+11  0.9546 0.85212  375.22   55.49 
      360.0   355.0    3.78 3.86e+11  0.9282 0.85006  364.87   54.71 
      360.0   360.0    3.73 3.86e+11  0.9027 0.84803  354.83   53.93 
      360.0   365.0    3.68 3.86e+11  0.8779 0.84597  345.08   53.15 
      360.0   370.0    3.63 3.86e+11  0.8537 0.84384  335.59   52.41 
      360.0   375.0    3.58 3.86e+11  0.8302 0.84159  326.35   51.70 
      360.0   380.0    3.53 3.86e+11  0.8074 0.83919  317.36   51.04 
      360.0   385.0    3.49 3.86e+11  0.7851 0.83658  308.62   50.43 
      360.0   390.0    3.44 3.86e+11  0.7635 0.83375  300.13   49.90 
      360.0   395.0    3.40 3.86e+11  0.7426 0.83065  291.91   49.43 
      360.0   400.0    3.35 3.86e+11  0.7224 0.82729  283.98   49.05 

License
-------

This project is licensed under the BSD 3-clause License.

Citation
--------

If you use this software in your research, please cite it as below:

.. code-block:: bibtex

    @misc{prahl_mie_scattering,
      author = {Scott Prahl},
      title = {Mie Scattering},
      year = {2023},
      doi = {10.5281/zenodo.10087653},
      url = {https://github.com/scottprahl/mie}
    }
