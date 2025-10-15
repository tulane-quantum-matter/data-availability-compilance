This is the code for http://arxiv.org/abs/2510.05230, "Boundary criticality in two-dimensional correlated topological superconductors", written by Yang Ge and Shao-Kai Jian. The auxillary-field QMC simulations in it were carried out with the ALF package [1] available at https://git.physik.uni-wuerzburg.de/ALF/ALF. Documentations and tutorials can be found on ALF website.

The files constitute a determinant quantum Monte Carlo submodule to be included in a ALF instance. All the `.F90` and `.inc` files should be included in `${PATH_TO_ALF}/Prog/Hamiltonians/`.
Then, the module name `TSC` should be added to `${PATH_TO_ALF}/Prog/Hamiltonians.list`.

An example `parameters` file is provided as well.

[1] ALF Collaboration, The ALF (Algorithms for Lattice Fermions) project release 2.0. Documentation for the auxiliary-field quantum Monte Carlo code, SciPost Phys. Codebases 1 (2022), URL: https://scipost.org/SciPostPhysCodeb.1.
