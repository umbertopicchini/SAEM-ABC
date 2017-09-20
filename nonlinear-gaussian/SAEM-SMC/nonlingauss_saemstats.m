function sufficientstats = nonlingauss_saemstats(xhatpredict,yobs,bigtheta,time)
%Define calculations for the sufficient statistics for unknow parameters THETA.
%Order of elements in the output should be consistent with the ordering of
%parameters in the THETA vector

xzero     = bigtheta(1);
log_sigmax    = bigtheta(2);
log_sigmay    = bigtheta(3);



%:::::::::::::::::::::::::::::::::::::::::::::::::::

S_sigmaxsquared = (xhatpredict(1) - 2*sin(exp(xzero))).^2 + sum( (xhatpredict(2:end)-2*sin(exp(xhatpredict(1:end-1)))).^2);
S_sigmaysquared = sum((yobs-xhatpredict).^2);


sufficientstats = [S_sigmaxsquared,S_sigmaysquared];


end

