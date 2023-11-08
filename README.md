# Mie Scattering

v2.6.0

February 2019

by Scott Prahl

[![github](https://img.shields.io/github/v/tag/scottprahl/mie?label=github&color=68CA66)](https://github.com/scottprahl/mie)
[![License](https://img.shields.io/github/license/scottprahl/mie?color=68CA66)](https://github.com/scottprahl/mie/blob/master/LICENSE.txt)
[![Docs](https://readthedocs.org/projects/mie/badge?color=68CA66)](https://mie.readthedocs.io)
![Downloads](https://img.shields.io/github/downloads/scottprahl/mie/total?color=68CA66)
[![doi](https://zenodo.org/badge/99259684.svg)](https://zenodo.org/badge/latestdoi/99259684)

___
## OVERVIEW

Yet another Mie scattering program.  This one was written in 1995 and based on Wiscombe's well-tested and documented MIEV0 FORTRAN code.  It works and is the basis for the on-line Mie scattering calculator at [omlc.org](http://omlc.org/calc/mie_calc.html)
I have used Wiscombe's Mie testing data to validate this code and it works for small spheres and very large spheres (sphere sizes with pi*d/lambda>10,000).

The source code is written in [CWEB](http://www-cs-faculty.stanford.edu/~knuth/cweb.html), which allows excellent documentation of scientific programs. Basically, there is a program `ctangle` that converts the cweb code to C. There is another program `cweave` that converts the cweb code to TeX. This then generates really [nice documentation](https://github.com/scottprahl/doc/mie_src.pdf), however if you are using the C code directly it means that you'll see none of my comments. 

There is also a scattering program buried in the `src` directory for cylinders that I translated from some FORTRAN code by Mackowski.  This is much less well tested.

I have written a pure python version of the sphere scattering code [miepython](https://github.com/scottprahl/miepython) that is based on this code, but has a number of nice affordances that make doing Mie calculations less of a hassle.

## Download

Visit [https://github.com/scottprahl/mie](https://github.com/scottprahl/mie) and download a zip file or just use `git`

    git clone https://github.com/scottprahl/mie.git

## INSTALLATION

In principle you should be able to just type

    make

then you can do some basic tests by

    make test
    
To install

    make install

and the binary executable file (`mie`) will be created and placed in a default location of
`/usr/local/bin/`.


### Shared library support.  

Edit the Makefile to select the right type of shared library for your platform

	make install-lib

### Mathematica support.  

If you have Mathematica, then (and only if you have installed the right
tools and edited the Makefile for your platform) and have the libraries installed then you should be able to type

	make mma
	make install mma

and then load the `mie` module to use in Mathematica.  Very cool.

### Python support

You'll just want to go to the [miepython](https://github.com/scottprahl/miepython) and download the pure python package.

## License

This code is licensed under the 3-clause MIT license.

## Author

Scott Prahl

`scott.prahl@oit.edu`

http://omlc.org/~prahl
