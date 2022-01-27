# SSDP

## The Simple Sky Dome Projector

The [Simple/Slick/Stupid/Smart/../Slimy/Sexy] Sky Dome Projector (SSDP)
(have not decided yet which adjective fits best, there is just so many) 
is a library for modeling irradiance taking into account the local 
topography. It implements a hexagonal mesh for the sky and the Perez 
all weather sky model. For the topography it works with a point cloud 
of x,y,z coordinates or a regular grid. It can compute the horizon as 
seen from a specific location in the topography and can project the 
modeled sky dome onto a surface with a specified tilt and orientation. 
The surface location and orientation may be adapted to the topography 
(e.g. specify height above the surface and adapt the orientation to the 
surface normal).

## Features
* Implements a uniform sky and the Perez all weather sky model [1]
	* variable size hexagonal sky-dome mesh
	* sun taken separately as a point source a alleviate sky-mesh 
	resolution requirements
* Computes solar position given the latitude and lonitude and the 
UTC time (specified as Unix Time) using the PCA algorithm [2]
	* also supports a manual specification of azimuth and zenith angles
* Simple and fast irradiance computation
	* easy specification of location (x,y,z) and orientation 
	(azimuth,zenith) of the incident plane
	* Traces incident light
	* Crude "one bounce" approximation for albedo (no ray tracing, it 
	is a fully "local" model)
	* only computes irradiance at specified locations
* 2.5D topography
	* Supports both unstructured and regular meshes for the topography

## Examples
As a validation of the irradiance and transposition models in ssdp we 
compare simulated plane of array irradiance values with 
[pvlib](https://github.com/pvlib/pvlib-python). To this end we take
measured GHI, DHI, and POA data from a location somewhere in the south 
of germany, and use the measured GHI and DHI to predict the POA 
irradiance. The POA in this case is a 40 degrees tilt, facing south. We 
assume a ground albedo of 0.25. For this simulation we do not take into 
account any topography (the data is for an open field) so we only look 
at irradiance and transposition modeling. For both ssdp and pvlib we 
use a Perez model (ssdp uses [1] and pvlib a slightly different model 
[3]). The results are shown in Fig. 1 for ssdp (left) and pvlib 
(right). Both models perform well with a coefficient of determination 
of more than R<sup>2</sup>=0.995. Small differences (~0.1% range) are 
observed and expected due to the different methods and models.

![Cumputed Irradiance](POA_ssdp.png) ![Cumputed Irradiance](POA_pvlib.png)

_Fig. 1 Transposition and irradiance modeling with (left) ssdp, and 
(right) pvlib_
  
In Fig. 2 we show two examples. In Fig 2 (left) the surface irradiance 
is computed ate some location for the 15th of june 2015 at 08:00 for a 
Global Horizontal Irradiance of 500 W/m<sup>2</sup> and a Diffuse 
Horizontal Irradiance of 200 W/m<sup>2</sup>. In Fig. 2 (right) we 
integrated the irradiance for a complete summer week
 
![Cumputed Irradiance](park_irr.png) ![Cumputed Irradiance](park_int.png)

_Fig. 2 (left) Example irradiance computation at one particular moment. 
(right) Example for the integrated irradiance over one summer week_

## Installation
The SSDP program comes with autotools build scripts, hence the 
installation procedure follows the standard:

1. `./configure`
2. `make`
3. `(sudo) make install` 

This installs both the ssdp library and the ssdp command-line program. 

The configure command takes several optional arguments. Specific ssdp 
options are:

* --enable-openmp: Use openmp palatalization (disabled per default) 
* --disable_fastatan2: Disables a fast atan2 approximation

Compiler optimizations:

Do not use unsafe math optimizations. With gcc/clang I would recommend 
'-O3' and '-march=native', e.g.

`./configure CFLAGS='-march=native -O3' --enable-openmp`

## Documentation

The ssdp program comes with man pages. If you have the man program and
ssdp installed you can do

`man ssdp` 

to read the documentation. Otherwise you can produce a readable manual 
using the groff program and the ssdp.man file, e.g.:

`groff -man -T pdf src/ssdp.man > ssdp_manual.pdf`

(assuming your working directory is the root source directory of ssdp)

## References
[1] R. Perez, et al. "All-Weather Model for Sky Luminance Distribution 
-- Preliminary Configuration and Validation." Solar Energy  50.3 
(1993):235-245

[2] M. Blanco-Muriel, et al. "Computing the solar vector." Solar 
Energy 70.5 (2001): 431-441

[3] R. Perez, et al. "Modeling daylight availability and irradiance 
components from direct and global irradiance". Solar Energy 44.5 
(1990):271-289.

