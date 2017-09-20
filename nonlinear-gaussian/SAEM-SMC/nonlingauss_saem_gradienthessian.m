function [gradient,hessian] = nonlingauss_saem_gradienthessian(xhatpredict,yobs,bigtheta,time)
%this is not necessary for SAEM to run

xzero = bigtheta(1);
log_sigmax = bigtheta(2);
log_sigmay = bigtheta(3);

sigmax = exp(log_sigmax);
sigmay = exp(log_sigmay);


nobs = length(time);

gradient_sigmaxsquared = -nobs/(2*sigmax^2) + 1/(2*sigmax^4)* ( (xhatpredict(1) - 2*sin(exp(xzero))).^2 + sum( (xhatpredict(2:end)-2*sin(exp(xhatpredict(1:end-1)))).^2) );
gradient_sigmaysquared = -nobs/(2*sigmay^2) + 1/(2*sigmay^4)*sum((yobs-xhatpredict).^2);

gradient = [gradient_sigmaxsquared; gradient_sigmaysquared];

% diagonal elements of the Hessian matrix
ddsigmaysquared = nobs/(2*sigmay^4) -1/sigmay^6 * sum((yobs-xhatpredict).^2); 
ddsigmaxsquared = nobs/(2*sigmax^4) -1/sigmax^6 * ( (xhatpredict(1) - 2*sin(exp(xzero))).^2 + sum( (xhatpredict(2:end)-2*sin(exp(xhatpredict(1:end-1)))).^2) ); 


hessian = [ddsigmaxsquared          0    ;... 
             0         ddsigmaysquared  ];


end

