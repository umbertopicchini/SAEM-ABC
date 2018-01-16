# SAEM-ABC, maximum likelihood estimation using an appoximate Bayesian computation (ABC) particle filter

This repository contains example code for the paper by Picchini and Samson 2017, "Coupling stochastic EM and Approximate Bayesian Computation for parameter inference in state-space models" https://arxiv.org/abs/1512.04831 and forthcoming on "Computational Statistics".

The code was written by Umberto Picchini (firstname dot secondname at gmail.com)

Notice there are several README files. Please read those before contacting the author. Not everything is commented in great detail, and it is beneficial to read the paper above first.

The core code illustrating the methodology is written in MATLAB. Some comparisons are written in R, see the folders "pomp_IF2", "pomp_PMCMC" and "theophylline_PMCMC".

# Warning:
The function abcsmc_filter in the theophylline folder currently calls a function "XSum". However this one might be not necessary and might be just ok if all XSum() calls are substituted with the regular sum().

XSum will actually call a MEX version, therefore:

(i) if you are a Windows user you are good to go as a MEX file compiled for a Windows x64 architecture is available (XSum.mexw64).

(ii) if you are not a Windows user you will need to compile XSum.c to obtain a MEX file suitable for your architecture.
XSum was retrieved from https://se.mathworks.com/matlabcentral/fileexchange/26800-xsum?s_tid=gn_loc_drop
