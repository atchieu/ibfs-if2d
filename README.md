IBFS-IncompressibleFluids-2D 
====

Immersed boundary fractional step method written in FORTRAN 90

See the following references

Colonius and Taira, 2008, CMAME
Taira and Colonius, 2007, JCP

MGRID VERSION 2.10, Release Notes 

This is a multi-grid, fractional-step, immersed boundary code used to solve the incompressible Navier-Stokes equations in two dimensions. 

Features:

- Solves 2D Navier-Stokes equation with an imbedded immersed boundary/body
- 2nd order accuracy in space away from immersed boundary
- 2nd order accuracy in time
- 1st order accuracy in vicinity of the immersed boundary
- Requires uniform grid for use of multi-grid method
- Discrete streamfunction/vorticity formulation removes divergence to machine precision
- Fast sine transforms used to invert and perform Laplacian operations
- Cholesky factorization used to solve modified Poisson's equaitons for stationary bodies

NOTES: Please make sure that you adjust the source code if you are planning to move the body. This has not been coded to be handled outside the source yet. Please see grid.f90 for more details.

