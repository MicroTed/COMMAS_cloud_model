## COMMAS: Collaborative Model for Multiscale Atmospheric Simulation

* September 2026: First public release

## Requirements:

Fortran compiler (gfortran works fine). An MPI compiler is recommended. NetCDF4 is needed for I/O, and optionally pnetcdf is recommended for high core counts (more than 80-100). Compression is supported for Netcdf4, which also has a compressed parallel-IO option (experimental: can be unreliable). For electrification with MPI, the BoxMG libraries are required.

BoxMG: https://github.com/MicroTed/boxmg4wrf

Electrification: Initial release does not have 3D branched lightning, but will be added

## Support:

This software is provided 'as is' with no guarantees, etc. Questions/problems can be brought up in the "Issues" section of the github repository:

https://github.com/MicroTed/COMMAS_cloud_model


