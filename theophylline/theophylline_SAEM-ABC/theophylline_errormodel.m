function yobssim = theophylline_errormodel(bigtheta,xhatpredict,numparticles)


X0 = bigtheta(1);
log_Ke = bigtheta(2);
log_Ka = bigtheta(3);
log_Cl = bigtheta(4);
log_sigma = bigtheta(5);
log_sigmaepsilon = bigtheta(6);

sigmaepsilon = exp(log_sigmaepsilon);
yobssim = xhatpredict + sigmaepsilon*randn(1,numparticles);


end

