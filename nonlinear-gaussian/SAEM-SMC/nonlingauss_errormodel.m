function yobssim_next = nonlingauss_errormodel(bigtheta,xhat,numparticles)
% create noisy observations given a noise-free state xhat 

X0 = bigtheta(1);
log_sigmax = bigtheta(2);
log_sigmay = bigtheta(3);

sigmay = exp(log_sigmay);

yobssim_next = xhat + sigmay*randn(1,numparticles);


end

