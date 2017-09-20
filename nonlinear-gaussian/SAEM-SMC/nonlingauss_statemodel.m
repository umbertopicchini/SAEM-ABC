function xhat_next = nonlingauss_statemodel(bigtheta,time_id,xhat_pre,numparticles)


X0 = bigtheta(1);
log_sigmax = bigtheta(2);
log_sigmay = bigtheta(3);

sigmax = exp(log_sigmax);


if time_id == 1
    xhat_next = 2*sin(exp(X0))*ones(1,numparticles) + sigmax*randn(1,numparticles);
else
    xhat_next = 2*sin(exp(xhat_pre)) + sigmax*randn(1,numparticles);
end





end

