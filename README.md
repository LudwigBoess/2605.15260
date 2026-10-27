# Entity – Hardware-agnostic Particle-in-Cell Code for Plasma Astrophysics
## III: Higher-order shape functions & generalized field stencils

In this repository you will find all scripts, pgens, and parameter files to reproduce the figures presented in Böss et al. (2026, submitted) and [Böss et al. (2026b)](https://ui.adsabs.harvard.edu/abs/2026arXiv260515260B/abstract).

The scripts require Julia >= 1.10 to be installed. I reccomend using v1.12.

To initialize the dependencies run `julia build.jl`.

This will install all packages.

You can get all plots by running from the main repository directory: `julia scripts/Fig01.jl`.
If you use an interactive session just ignore the warnings about redefinition of global variables.
If you use an interactive session just ignore the warnings about redefinition of global variables.
