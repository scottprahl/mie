(* ::Package:: *)

(* :Name:Optics`Mie*)

(* :Title:Mie Scattering Calculations*)

(* :Author:Scott Prahl*)

(* :Summary:This package finds the scattering parameters for spheres.*)

(* :Context:Optics`Mie`*)

(* :Package Version:1.6*)

(* :Copyright:Copyright 2011, Oregon Medical Laser Center*)

(* :History:Originally written by Scott Prahl, November 2002 *)

(* :Discussion: A bunch of handy Mie scattering routines *)

BeginPackage["Optics`Mie`"]

Mie::usage = "Mie[x,nreal,nimag,{cos_angles}] returns returns {Qsca,g,Qext,Qback,{S1real},{S1imag},{S2real},{S2imag}} for spheres with refractive index n=n-ik and size parameter of x for {Cos[angles]}."

MieMus::usage = "MieMus[f,d,n,lambda] and MieMus[f,d,n,m,lambda] calculate the scattering coefficient (per mm) for a sphere with diameter d (in nanometers), index of refraction n=n-ik (can be complex), in a medium with refractive index m, at the wavelength lambda (in nanometers) for a volume fraction f (f=0.41 means 41% spheres by volume)."

MieMusp::usage = "MieMusp[f,d,n,lambda] and MieMusp[f,d,n,m,lambda] calculate the reduced scattering coefficient (per mm) for a sphere with diameter d (in nanometers), index of refraction n=n-ik (can be complex), in a medium with refractive index m, at the wavelength lambda (in nanometers) for a volume fraction f (f=0.41 means 41% spheres by volume)."

MieG::usage = "MieG[d,n,lambda] and MieG[d,n,m,lambda] calculate the scattering anisotropy for a sphere with diameter d (in nanometers), index of refraction n=n-ik (can be complex), in a medium with real index m, at the wavelength lambda (in nanometers)."

MieQsca::usage = "MieQsca[d,n,lambda] and MieQsca[d,n,m,lambda] calculate the scattering efficiency for a sphere with diameter d (in nanometers), index of refraction n=n-ik (can be complex), in a medium with real index m, at the wavelength lambda (in nanometers)."

MieQscap::usage = "MieQscap[d,n,lambda] and MieQscap[d,n,m,lambda] calculate the reduced scattering efficiency for a sphere with diameter d (in nanometers), index of refraction n=n-ik (can be complex), in a medium with real index m, at the wavelength lambda (in nanometers)."

MieQext::usage = "MieQext[d,n,lambda] and MieQext[d,n,m,lambda] calculate the extinction efficiency for a sphere with diameter d (in nanometers), index of refraction n=n-ik (can be complex), in a medium with real index m, at the wavelength lambda (in nanometers)."

MieQback::usage = "MieQback[d,n,lambda] and MieQback[d,n,m,lambda] calculate the back-scattering efficiency for a sphere with diameter d (in nanometers), index of refraction n=n-ik (can be complex), in a medium with real index m, at the wavelength lambda (in nanometers)."

MiePhaseFunction::usage = "MiePhaseFunction[d,n,lambda,{angles}] and MiePhaseFunction[d,n,m,lambda,{angles}] returns {{S1},{S2}} for a sphere with diameter d (in nanometers), index n, in a medium with real index m, at each of the directions in radians {angles} for light with wavelength lambda (in nanometers in a vacuum)."

MiePhasePar::usage = "MiePhasePar[d,n,lambda,{angles}] and MiePhasePar[d,n,m,lambda,{angles}] returns the phase function for incident light polarized parallel to the plane of incidence on a sphere with diameter d (in nanometers), complex index n, in a medium with real index m, at each of the directions in radians {angles} for light with wavelength lambda (in nanometers in a vacuum)."

MiePhasePer::usage = "MiePhasePer[d,n,lambda,{angles}] and MiePhasePer[d,n,m,lambda,{angles}] returns the phase function for light polarized perpendicular to the plane of incidence on a sphere with diameter d (in nanometers), complex index n, in a medium with real index m, at each of the directions in radians {angles} for light with wavelength lambda (in nanometers in a vacuum)."

MiePhaseUnpolarized::usage = "MiePhaseUnpolarized[d,n,lambda,{angles}] and MiePhaseUnpolarized[d,n,m,{angles}] returns the phase function for unpolarized light incident on a sphere with diameter d (in nanometers), complex index n, in a medium with real index m, at each of the directions in radians {angles} for light with wavelength lambda (in nanometers in a vacuum)."

MiePhasePolarization::usage = "MiePhasePolarization[d,n,lambda,{angles}] and MiePhasePolarization[d,n,m,{angles}] returns the polarization state for unpolarized light incident on a sphere with diameter d (in nanometers), complex index n, in a medium with real index m, at each of the directions in radians {angles} for light with wavelength lambda (in nanometers in a vacuum)."

MieMatrix::usage = "MieMatrix[d,n,m,lambda,theta] and MieMatrix[d,n,m,lambda,theta,phi] return the Mueller matrix for a sphere with diameter d (in nanometers), index n-ik, in a medium with real index m, for light with wavelength lambda (in nanometers in a vacuum) for light scattered at an angle theta (with azimuthal angle phi).\n\nTo find the Stokes vector for an horizontally polarized light and analyzed by a vertical polarizer just do PolarizerMatrixV . MieMatrix[x,n,theta,phi] . StokesVectorH"

StokesVectorH::usage = "StokesVectorH is the Stokes vector for horizontal polarized light"

StokesVectorV::usage = "StokesVectorV is the Stokes vector for vertical polarized light"

StokesVectorPPlus::usage = "StokesVectorPPlus is the Stokes vector for light polarized at +45 degrees"

StokesVectorPMinus::usage = "StokesVectorPMinus is the Stokes vector for light polarized at -45 degrees"

StokesVectorR::usage = "StokesVectorR is the Stokes vector for right circularly polarized light"

StokesVectorL::usage = "StokesVectorL is the Stokes vector for left circularly polarized light"

PolarizerMatrixH::usage = "PolarizerMatrixH is the Mueller matrix for a horizontal polarizer"

PolarizerMatrixV::usage = "PolarizerMatrixV is the Mueller matrix for a vertical polarizer"
      
PolarizerMatrixPPlus::usage = "PolarizerMatrixPPlus is the Mueller matrix for a polarizer oriented at +45 degrees"
      
PolarizerMatrixPMinus::usage = "PolarizerMatrixPMinus is the Mueller matrix for a polarizer oriented at -45 degrees"
      
PolarizerMatrixR::usage ="PolarizerMatrixR is the Mueller matrix for a right circular polarizer"
   
PolarizerMatrixL::usage = "PolarizerMatrixL is the Mueller matrix for a left circular polarizer"

myRotationMatrix::usage = "myRotationMatrix[theta] returns the matrix needed to rotate the Stokes vector by an angle theta. "

Mie::badVolumeFraction = "Volume Fraction must be between zero and one."

Begin["`Private`"]

Install["Optics/External/Mie"]

Mie[x_,n_Complex,mu_List] := Mie[x,Re[n],Im[n],mu]
Mie[x_,n_,mu_List]        := Mie[x,n,0.0,mu]
Mie[x_,n_Complex]         := Mie[x,Re[n],Im[n],{0.5}][[1]]
Mie[x_,n_]                := Mie[x,n,0.0,{0.5}][[1]]

MieQsca[diameter_, nSphere_, lambdaVacuum_] := 
    Mie[Pi diameter/lambdaVacuum, nSphere][[1]]

MieQsca[diameter_, nSphere_, nMedium_, lambdaVacuum_] :=
    MieQsca[diameter, nSphere/nMedium, lambdaVacuum/nMedium]

MieG[diameter_, nSphere_, lambdaVacuum_] := 
    Mie[Pi diameter/lambdaVacuum, nSphere][[2]]

MieG[diameter_, nSphere_, nMedium_, lambdaVacuum_] := 
	MieG[diameter, nSphere/nMedium, lambdaVacuum/nMedium]

MieQscap[diameter_, nSphere_, lambdaVacuum_] := 
	Block[{a},
      a = Mie[Pi diameter/lambdaVacuum, nSphere];
      a[[1]](1 - a[[2]])
	]

MieQscap[diameter_, nSphere_, nMedium_, lambdaVacuum_] :=
    MieQscap[diameter, nSphere/nMedium, lambdaVacuum/nMedium]
    
MieQext[diameter_, nSphere_, lambdaVacuum_] := 
    Mie[Pi diameter/lambdaVacuum, nSphere][[3]]

MieQext[diameter_, nSphere_, nMedium_, lambdaVacuum_] :=
    MieQext[diameter, nSphere/nMedium, lambdaVacuum/nMedium]
    
MieQback[diameter_, nSphere_, lambdaVacuum_] := 
    Mie[Pi diameter/lambdaVacuum, nSphere][[4]]

MieQback[diameter_, nSphere_, nMedium_, lambdaVacuum_] :=
    MieQback[diameter, nSphere/nMedium, lambdaVacuum/nMedium]

(* Properties that apply to bulk *)
  
MieMus[VolumeFraction_, diameter_, nSphere_, lambdaVacuum_] := If[VolumeFraction>1,
	Message[Mie::badVolumeFraction],
    (1500000 VolumeFraction/diameter) MieQsca[diameter,nSphere,lambdaVacuum]]

MieMus[VolumeFraction_, diameter_, nSphere_, nMedium_, lambdaVacuum_] :=
    MieMus[VolumeFraction, diameter, nSphere/nMedium, lambdaVacuum/nMedium]

MieMusp[VolumeFraction_, diameter_, nSphere_, lambdaVacuum_] := 
If[VolumeFraction>1,
	Message[Mie::badVolumeFraction],
    (1500000 VolumeFraction/diameter) MieQscap[diameter,nSphere,lambdaVacuum]]
      
MieMusp[VolumeFraction_, diameter_, nSphere_, nMedium_, lambdaVacuum_]  :=
  MieMusp[VolumeFraction, diameter, nSphere/nMedium, lambdaVacuum/nMedium]
  
(* Properties that vary with angle *)

MiePhaseFunction[diameter_, nSphere_, lambdaVacuum_, angle_] := 
	Block[{a,s1re,s1im,s2re,s2im},
      {a,s1re,s1im,s2re,s2im} = Mie[diameter/lambdaVacuum, nSphere, angle];
      { s1re + I s1im, s2re + I s2im}
	]

MiePhasePar[diameter_, nSphere_, lambdaVacuum_, angle_] := 
 	Block[{a,b,s1re,s1im,s2re,s2im},
      {a,s1re,s1im,s2re,s2im} = Mie[Pi diameter/lambdaVacuum, nSphere, angle];
      b = s2re^2 + s2im^2;
      If[Length[b]>1, b, b[[1]]]
      ]

MiePhasePer[diameter_, nSphere_, lambdaVacuum_, angle_] := 
	Block[{a,b,s1re,s1im,s2re,s2im},
      {a,s1re,s1im,s2re,s2im} = Mie[Pi diameter/lambdaVacuum, nSphere, angle];
      b = s1re^2 + s1im^2;
      If[Length[b]>1, b, b[[1]]]
      ]

MiePhaseUnpolarized[diameter_, nSphere_, lambdaVacuum_, angle_] := 
	Block[{a,b,s1re,s1im,s2re,s2im},
      {a,s1re,s1im,s2re,s2im} = Mie[Pi diameter/lambdaVacuum, nSphere, angle];
      b = (s1re^2 + s1im^2 + s2re^2 + s2im^2)/2.0;
      If[Length[b]>1, b, b[[1]]]
	]

MiePhasePolarization[diameter_, nSphere_, lambdaVacuum_, angle_] := 
	Block[{a,b,s1re,s1im,s2re,s2im},
      {a,s1re,s1im,s2re,s2im}=Mie[Pi diameter/lambdaVacuum, nSphere, angle];
      b = (s1re^2 + s1im^2 - s2re^2 - s2im^2)/(s1re^2 + s1im^2 + s2re^2 + s2im^2);
      If[Length[b]>1, b, b[[1]]]
	]

MiePhaseFunction[ diameter_, nSphere_, nMedium_, lambdaVacuum_, angle_] := 
	MiePhaseFunction[diameter, nSphere/nMedium, lambdaVacuum/nMedium, angle]

MiePhasePar[ diameter_, nSphere_, nMedium_, lambdaVacuum_, angle_] := 
	MiePhasePar[diameter, nSphere/nMedium, lambdaVacuum/nMedium, angle]

MiePhasePer[ diameter_, nSphere_, nMedium_, lambdaVacuum_, angle_] := 
	MiePhasePer[diameter, nSphere/nMedium, lambdaVacuum/nMedium, angle]

MiePhaseUnpolarized[ diameter_, nSphere_, nMedium_, lambdaVacuum_, angle_] := 
	MiePhaseUnpolarized[diameter, nSphere/nMedium, lambdaVacuum/nMedium, angle]

MiePhasePolarization[ diameter_, nSphere_, nMedium_, lambdaVacuum_, angle_] := 
	MiePhasePolarization[diameter, nSphere/nMedium, lambdaVacuum/nMedium, angle]
	
MieMatrix[x_, n_, theta_] := Block[{s1, s2, s11, s12, s33, s34, v},
    v = Mie[x, n, theta];
    s1 = v[[2]] + I v[[3]];
    s2 = v[[4]] + I v[[5]];
    s11 = ((Abs[s2[[1]]]^2 + Abs[s1[[1]]]^2)/2);
    s12 = ((Abs[s2[[1]]]^2 - Abs[s1[[1]]]^2)/2);
    s33 = Re[Conjugate[s1[[1]]] s2[[1]]];
    s34 = Im[Conjugate[s1[[1]]] s2[[1]]];
    {{s11, s12, 0, 0}, {s12, s11, 0, 0}, {0, 0, s33, s34}, {0, 0, s34, s33}}
    ]

MieMatrix[x_, n_, theta_, phi_] := myRotationMatrix[-phi Sign[Pi/2 - theta]] . 
                                    MieMatrix[x, n, theta] . 
                                    myRotationMatrix[phi]


StokesVectorH = {1, 1, 0, 0}

StokesVectorV = {1, -1, 0, 0}

StokesVectorPPlus = {1, 0, 1, 0}

StokesVectorPMinus = {1, 0, -1, 0}

StokesVectorR = {1, 0, 0, 1}

StokesVectorL = {1, 0, 0, -1}

PolarizerMatrixH = 0.5{{1, 1, 0, 0}, {1, 1, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}

PolarizerMatrixV = 0.5{{1, -1, 0, 0}, {-1, 1, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}
      
PolarizerMatrixPPlus = 0.5{{1, 0, 1, 0}, {1, 0, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}
      
PolarizerMatrixPMinus = 0.5{{1, 0, -1, 0}, {-1, 0, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}
      
PolarizerMatrixR = 0.5 {{1, 0, 0, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 0, 0, 1}}
   
PolarizerMatrixL = 0.5 {{1, 0, 0, -1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 0, 0, -1}}

myRotationMatrix[angle_] := {{1,0,0,0},{0,Cos[2 angle],Sin[2 angle],0},{0,-Sin[2 angle],Cos[2 angle],0},{0,0,0,1}}

End [ ]

EndPackage[ ]
