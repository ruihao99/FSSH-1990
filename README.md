# FSSH 
A c++ implementation of Tully's fewest switching surface hopping (FSSH).

John C. Tully , "Molecular dynamics with electronic transitions", [J. Chem. Phys. 93, 1061-1071 (1990)](https://doi.org/10.1063/1.459170)

Here is the structure for the code.

```
.
├── DrawPot.cpp           # Draws the potential, from x \in {-10 <= x <= 10}
├── README.md             # this file
└── include
    ├── integrator.hpp    # Header file for the FSSH_integrator
    ├── pot2d.hpp         # Header file for the 2D potentials
    └── utils.hpp         # some utilities (constants, linspace)
```

Two prerequisite libraries for compiling this code are

* Eigen3: For matrix calculation
* boost/odeint: For molecular dynamics (RK4) and wavefunction Integration (RK4).
