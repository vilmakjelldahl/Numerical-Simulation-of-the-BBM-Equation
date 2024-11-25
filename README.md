# Numerical-Simulation-of-the-BBM-Equation
Numerical Simulations of the generalized modified Benjamin-Bona-Mahony equation using SBP-SAT in time. 
https://doi.org/10.1016/j.cam.2024.116377

This repository contains Matlab scripts using high-order accurate finite difference approximations for solving the generalized modified Benjamin-Bona-Mahony (BBM) equation, a non-linear soliton model. The spatial discretization uses high-order accurate summation-by-parts (SBP) finite difference operators combined with both weak and strong enforcement of well-posed boundary conditions, using the SAT-method and projection method respectively. For time integration an explicit RK4 method is used, as well as an implicit SBP time integrator with different types of SBP-operators.  It is shown that the implicit SBP time-integrator is more efficient than the explicit RK4 method for non-linear soliton models. 


Authors: Vilma Kjelldahl, Ken Mattsson
