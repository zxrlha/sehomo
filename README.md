[Bootstrapping the Solutions of Scattering Equations](https://arxiv.org/abs/1810.00384)
==================================================

## General Information
Scattering Equations are a system of algebraic equations in quantum field theory.
The algorithm presented in [arXiv:1810.00384](https://arxiv.org/abs/1810.00384) is intended to solve those equations numerically, which is based on numerical algebraic geometry.
Here is the simple implementation in that paper.

## Download:
You can download it directly using git:

```Shell
git clone https://github.com/zxrlha/sehomo.git
```

## Installation and Running

To compile this package, a C++ compiler with C++11 support is required.
We recommend [GCC](https://gcc.gnu.org) version 4.9 or later.

Besides basic compiling tool chain, two external packages are required:

1. [Eigen](http://eigen.tuxfamily.org)
2. [Boost](https://boost.org)

For Linux distributions, usually those packages are available from the your package manager.

With those package installed, you can simple type

```Shell
make
```
to compile those code.
Two executable are generated: one is **prep_ic** and the other is **run_phy**.

You should first run **prep_ic**, then several files with names from **ic-5** to **ic-9** are created.
Those files store the initial conditions from 5-point kinematics to 9-point kinematic.
With those initial conditions, you can run **run_phy** to solve the scattering equations for physical points.


## LICENSE

GNU GPL v3 or later(see "COPYING" for detail)

## Brief explanation for the code

- lorentz.hpp: the class representing four-momentum
- rcm.hpp: generating phasespace for 2->(n-2) scattering
- timer.hpp and timer.cpp: the class for timing
- posimp.hpp and posimp.cpp: the implementation on utilizing inverse soft homotopy in positive kinematics(Step 1.1)
- comimp.hpp and comimp.cpp: the implementation on tracking from positive kinematics to physical kinematics(Step 1.2)
- phyimp.hpp and phyimp.cpp: the implementation on tracking from one physical point to other physical points(Step 2)
- prep_ic.cpp: preparing the initial condition for one physical point(i.e. call posimp and comimp for Step 1.1 and 1.2, correspondingly)
- calc_phys.cpp: calculating points in physical region using the previous generated initial condition(i.e. call phyimp for Step 2)
