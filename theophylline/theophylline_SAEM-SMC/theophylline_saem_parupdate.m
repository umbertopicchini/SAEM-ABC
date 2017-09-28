function theta = theophylline_saem_parupdate(bigtheta,sufficientstats,time,owntime)
% the M-step in SAEM

SSKeClratio =  sufficientstats(1);
SSKe =  sufficientstats(2);
SSsigmasquared = sufficientstats(3);
SSsigmaepsilonsquared = sufficientstats(4);

nobs = length(time);
N = length(owntime);

Ke_upd = SSKe;
Cl_upd = SSKe/SSKeClratio;
sigmasquared_upd = SSsigmasquared/N;
sigmaepsilonsquared_upd = SSsigmaepsilonsquared/nobs;


log_Ke = log(Ke_upd);
log_Cl = log(Cl_upd);
log_sigma = log(sqrt(sigmasquared_upd));
log_sigmaepsilon = log(sqrt(sigmaepsilonsquared_upd));

%theta = [log_Ke,log_Cl,log_sigma,log_sigmaepsilon];
theta = [log_Ke,log_Cl,log_sigma,log_sigmaepsilon];



end

