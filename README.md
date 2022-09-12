# FSSH 
A c++ implementation of Tully's fewest switching surface hopping (FSSH).

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

Two prerequisite library for compiling this code is

* Eigen3: For matrix calculation
* boost/odeint: For molecular dynamics (RK4) and wavefunction Integration (RK4).
