## README.txt
## Author -YH 2021-2026

This folder contains a couple versions of the algorithm HBPPDM (Heavy-ball projected primal-dual method) and helper scripts to build and execute ocp problems

# HBPPDM_var_reset.m
Version of the algorithm without infeasibility detection, and with variable stepsize for the heavy-ball acceleration
The variable stepsize is reset once the problem starts to diverge due to high momentum.

# HBPPDM_feas.m
Version of the algorithm with all features of the CDC paper (HBPPDM_var_reset.m) including infeasibility detection.

# HBPPDM_blksqp.m
Version of the algorithm for non-diagonal weighting matrices G. (Version used in Automatica)

# plotting_header.m
Header file for generating plots for double column articles.

# genOCPproblem.m
Legacy version for generating ocp problems. Use fastgenOCPproblem.m instead.

# fastgenOCPproblem.m
Fast generation of (sparse) OCP problems. Sets up all matrices that are required to execute the solver.

# prepareOSQP.m
Reparse the problem such that it runs with osqp (https://osqp.org/)

# prepareSCS.m
Reparse the problem such that it runs with scs (https://github.com/cvxgrp/scs)