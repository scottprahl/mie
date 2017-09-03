# Mie Scattering

** v2.5.0 **

** Scott Prahl **

** September 2017 **

** Oregon Tech **

## OVERVIEW

Yet another Mie scattering program.  This one was based on Wiscombe's well-tested and documented MIEV0 fortran code.  It works and is the basis for the on-line Mie scattering calculator at [omlc.org](http://omlc.org/calc/mie_calc.html)

I recently (2017) wrote a pure python version of the code [miepython](https://github.com/scottprahl/miepython) that is based on this code, but has a number of nice affordances that make doing Mie calculations less of a hassle.

## Download

Visit [https://github.com/scottprahl/mie] and download a zip file or just use `git`

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

## Author

Scott Prahl

`scott.prahl@oit.edu`

http://omlc.org/~prahl
