# SSDP

## The Simple Sky Dome Projector

The [Simple/Slick/Stupid/Smart/../Slimy/Sexy] Sky Dome Projector (SSDP)
(have not decided yet which adjective fits best, there is just so many) 
is a library for modeling irradiance taking into accound the local 
topography. It implements a hexagonal mesh for the sky and the Perez 
all weather sky model. For the topography it works with a point cloud 
of x,y,z coordinates or a regular grid. It can compute the horizon as 
seen from a specific location in the topography and can project the 
modelled sky dome onto a surface with a specified tilt and orientation. 
The surface location and orientation may be adapted to the topography 
(e.g. specify height above the surface and adapt the orientation to the 
surface normal).

## Features
* Implements a uniform sky and the Perez all weather sky model [1]
	* variable size hexagonal sky-dome mesh
	* sun taken seperrately as a point source a alleviate sky-mesh 
	resolution requirements
* Computes solar position given the latitude and lonitude and the 
UTC time (specified as Unix Time) using the PCA algorithm [2]
	* also supports a manual specification of azimuth and zenith angles
* Simple and Fast Irradiance computation
	* easy specification of location (x,y,z) and orientation 
	(azimuth,zenith) of the incident plane
	* Traces incident light
	* Crude "one bounce" approximation for albedo (no ray tracing, it 
	is a fully "local" model)
	* only computes irradiance at specified locations
* 2.5D topography
	* Supports both unstructured and regular meshes for the topography

## Examples
In Fig. 1 we show two examples. In Fig 1a the surface irradiance is 
computed ate some location for the 15th of june 2015 at 08:00 for a 
Global Horizontal Irradiance of 500 W/m^2 and a Diffuse Horizontal 
Irradiance of 200 W/m^2. In Fig. 1b we integrated the irradiance for a 
complete summer week
 
a.![Cumputed Irradiance](park_irr.png) b.![Cumputed Irradiance](park_int.png)

_Fig. 1 a. Example irradiance computation at one particular moment. b. 
Example for the integrated irradiance over one summer week_

## Installation
The SSDP program comes with autotools build scripts, hance the 
installation procedure follows the standard:

1. `./configure`
2. `make`
3. `(sudo) make install` 

This installs both the ssdp library and the ssdp cvommand-line program. 

The configure command takes several optional arguments. Specific ssdp 
options are:

* --enable-openmp: Use openmp parallelization (disabled per default) 
* --disable_fastatan2: Disables a fast atan2 approximation

Compiler optimizations:

Do not use unsafe math optimizations. With gcc/clang I would recommend 
'-O3' and '-march=native', e.g.

`./configure CFLAGS='-march=native -O3' --enable-openmp`

## Documentation

The ssdp program comes with man pages. If you have the man program 
installed you can do:

`man ssdp` 

to read the documentation. Otherwise you can produce a readable manual 
using the groff program e.g.:

`groff -man -T pdf src/ssdp.man > ssdp_manual.pdf`

## References
[1]: R. Perez, et al. "All-Weather Model for Sky Luminance Distribution 
-- Preliminary Configuration and Validation." Solar Energy  50.3  
(1993):235-245

[2] Blanco-Muriel, Manuel, et al. "Computing the solar vector." Solar 
energy 70.5 (2001): 431-441

