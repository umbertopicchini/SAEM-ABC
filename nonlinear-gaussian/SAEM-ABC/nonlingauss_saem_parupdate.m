function theta = nonlingauss_saem_parupdate(bigtheta,sufficientstats,time)
%This is the M-step in SAEM



SS_sigmaxsquared =  sufficientstats(1);
SS_sigmaysquared =  sufficientstats(2);


nobs = length(time);

sigmaxsquared_upd = SS_sigmaxsquared/nobs;
sigmaysquared_upd = SS_sigmaysquared/nobs;

log_sigmax = log(sqrt(sigmaxsquared_upd));
log_sigmay = log(sqrt(sigmaysquared_upd));

theta = [log_sigmax,log_sigmay];


end

