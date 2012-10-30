% This program by Scott Prahl
% It is distributed WITHOUT ANY WARRANTY, express or implied.

% Copyright (C) 2012 Scott Prahl

% Permission is granted to make and distribute verbatim copies of this
% document provided that the copyright notice and this permission notice
% are preserved on all copies.

% Permission is granted to copy and distribute modified versions of this
% document under the conditions for verbatim copying, provided that the
% entire resulting derived work is given a different name and distributed
% under the terms of a permission notice identical to this one.

\input version.tex

\def\title{Mie Scattering (\version)}
\def\topofcontents{\null\vfill
  \centerline{\titlefont Mie Scattering}
  \vskip 15pt
  \centerline{\version}
  \vfill}
\def\botofcontents{\vfill
\noindent
Copyright \copyright\ 2012 Scott Prahl
\bigskip\noindent
Permission is granted to make and distribute verbatim copies of this
document provided that the copyright notice and this permission notice
are preserved on all copies.

\smallskip\noindent
Permission is granted to copy and distribute modified versions of this
document under the conditions for verbatim copying, provided that the
entire resulting derived work is given a different name and distributed
under the terms of a permission notice identical to this one.
}

%\pageno=\contentspagenumber \advance\pageno by 1

@** Equations for Mie scattering.

The index of refraction $m$ of the sphere may be complex,
$$
m = n(1-i\kappa)
$$
The imaginary part of the complex refractive index $n\kappa$ is the damping factor while
$\kappa$ is called the index of absorption or the index of attenuation.  Note that
the sign of the imaginary part of the index of refraction is negative.
The complex index of refraction may also be written in terms of the
conductivity $\sigma$, the dielectric constant $\varepsilon$ and
the circular frequency $\omega$ as
$$
m = \sqrt{\varepsilon - i {4\pi\sigma\over\omega}}
$$

The other important parmeter governing scattering by a sphere is the size
parameter $x$ of the sphere, which is given by
$$
x = {2\pi a/\lambda}
$$
Sometimes the value $\rho$ is used to indicate the size of the sphere and
it is defined as
$$
\rho = 2x(m-1)
$$
and really only makes sense if the sphere does not absorb light.

The
absorption coefficient from Beer's law is defined as
$$
I = I_0 \exp(-\mu_a z)
$$
and thus
$$
\mu_a = {4\pi n\kappa\over\lambda_0} = {4\pi \kappa\over \lambda}
$$
where $\lambda_0$ is the wavelength in a vacuum [Kerker, p. 15].

Now to reprise some nomenclature.  The extinction efficiency may be separated into
$$
Q_{\hbox{ext}}=Q_{\hbox{sca}}+Q_{\hbox{abs}}
$$
where $Q_{\hbox{sca}}$ is the scattering efficiency and $Q_{\hbox{abs}}$ is the absorption
efficiency.  Typically $Q_{\hbox{sca}}$ and $Q_{\hbox{ext}}$ are determined by the
Mie scattering program and $Q_{\hbox{abs}}$ is obtained by subtraction.  

The radiation pressure is given by 
$$
Q_{\hbox{pr}}=Q_{\hbox{ext}}-g Q_{\hbox{sca}}
$$
The pressure exerted on the particle of cross-sectional area $\pi r_0^2$
is
$$
P = {F\over\pi r_0^2} = {Q_{\hbox{ext}}\over c}
$$
were $c$ is the velocity of the radiation in the medium
[Kerker, p. 94].

The relation between the efficiency factor for scattering and the cross section for
scattering are obtained by dividing by the actual geometrical cross section
$$
Q_{\hbox{sca}} = {C_{\hbox{sca}}\over \pi r_0^2}
$$
where $r_0$ is the radius of the sphere.

The scattering cross section may be related to the transmission of a beam 
through a dispersion of scatterers of equal size.  For $\rho$ particles per
unit volume, the attenuation due to scattering is
$$
-{dI\over dx} = \rho C_{\hbox{sca}} I
$$
The transmission is
$$
T = I/I_0 = \exp(-\rho C_{\hbox{sca}} x) = \exp(-\mu_s x)
$$
or
$$
\mu_s = \rho C_{\hbox{sca}} = \rho \pi r_0^2 Q_{\hbox{sca}}
$$
[Kerker, p. 38].

@i "mie_array.w"
@i "mie_complex.w"
@i "mie.w"
@i "mie_main.w"
@i "mie_cylinder.w"
@i "mie_cylinder_main.w"

@**Index.
Here is a cross-reference table for the Mie scattering program.
All sections in which an identifier is
used are listed with that identifier, except that reserved words are
indexed only when they appear in format definitions, and the appearances
of identifiers in section names are not indexed. Underlined entries
correspond to where the identifier was declared. Error messages and
a few other things like ``ASCII code dependencies'' are indexed here too.
