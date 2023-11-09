# Mie Scattering

by Scott Prahl

[![github](https://img.shields.io/github/v/tag/scottprahl/mie?label=github&color=68CA66)](https://github.com/scottprahl/mie)
[![License](https://img.shields.io/github/license/scottprahl/mie?color=68CA66)](https://github.com/scottprahl/mie/blob/master/LICENSE)
[![Docs](https://img.shields.io/badge/docs-passing-68CA66)](https://github.com/scottprahl/mie/blob/master/doc/mie_doc.pdf)
![Test](https://github.com/scottprahl/mie/actions/workflows/test.yaml/badge.svg)
[![DOI](https://zenodo.org/badge/102227669.svg)](https://zenodo.org/doi/10.5281/zenodo.10087653)

Mie Scattering is a computational project for simulating the scattering of electromagnetic waves by small particles. Utilizing the robust MIEV0 FORTRAN code by Wiscombe, this software offers researchers and scientists an accurate and reliable tool for Mie scattering calculations.

___

## OVERVIEW

Yet another Mie scattering program.  This one was written in 1995 and based on Wiscombe's well-tested and documented MIEV0 FORTRAN code.  It works and is the basis for the on-line Mie scattering calculator at [omlc.org](http://omlc.org/calc/mie_calc.html)
I have used Wiscombe's Mie testing data to validate this code and it works for small spheres and very large spheres (sphere sizes with 𝜋d/λ > 10,000).

The source code is written in [CWEB](https://github.com/ascherer/cweb), which allows excellent documentation of scientific programs. Basically, there is a program `ctangle` that converts the cweb code to C. There is another program `cweave` that converts the cweb code to TeX. This then generates really [nice documentation](https://github.com/scottprahl/mie/blob/master/doc/mie_doc.pdf), however if you are using the C code directly it means that you'll see none of my comments. 

There is also a scattering program buried in the `src` directory for cylinders that I translated from some FORTRAN code by Mackowski.  This is much less well tested.

I have written a pure python version of the sphere scattering code [miepython](https://github.com/scottprahl/miepython) that is based on this code.  It is a bit slower, but because it is python, it has a number of nice affordances that make doing Mie calculations (and plotting them) less of a hassle.

## Download

Visit [https://github.com/scottprahl/mie](https://github.com/scottprahl/mie/releases) and download a zip file or just use `git`

    git clone https://github.com/scottprahl/mie.git

## INSTALLATION

On a Unix/macos machine you should be able to just type

    make install

and the binary executable file (`mie`) will be created and placed in a default location of
`/usr/local/bin/`.

On windows, you can go to [Releases](https://github.com/scottprahl/mie/releases) and download a windows executable.

### Python support

You'll just want to go to the [miepython](https://github.com/scottprahl/miepython) and download the pure python package.

## License

This code is licensed under the 3-clause MIT license.
