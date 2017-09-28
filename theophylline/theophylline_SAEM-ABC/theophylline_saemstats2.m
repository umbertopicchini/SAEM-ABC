function sufficientstats = theophylline_saemstats(xhat_selected_big,xhat_selected_small,yobs,bigtheta,owntime,time)
%Define calculations for the sufficient statistics for unknow parameters THETA.
%Order of elements in the output should be consistent with the ordering of
%parameters in the THETA vector

xzero     = bigtheta(1);
log_Ke    = bigtheta(2);
log_Ka    = bigtheta(3);
log_Cl    = bigtheta(4);
log_sigma = bigtheta(5);
log_sigmaepsilon = bigtheta(6);

Ke = exp(log_Ke);
Ka = exp(log_Ka);
Cl = exp(log_Cl);
sigma = exp(log_sigma);
sigmaepsilon = exp(log_sigmaepsilon);

% xhat
% length(xhat)
% plot(time,xhat,'k-')
% hold on
% plot(time,yobs,'mo-')

Dose = 4;
%delta = time(2)-time(1);
h = owntime(2)-owntime(1);

%:::::::::::::::::::::::::::::::::::::::::::::::::::
N = length(owntime);

S_sigmaepsilonsquared = sum((yobs-xhat_selected_small).^2);

sum_sigmasquared = 0;

for ii = 2:N
        sum_sigmasquared = sum_sigmasquared + (xhat_selected_big(ii) -xhat_selected_big(ii-1)-(Dose*Ka*Ke/Cl*exp(-Ka*owntime(ii-1))-Ke*xhat_selected_big(ii-1))*h)^2/(h*xhat_selected_big(ii-1));
end
S_sigmasquared = sum_sigmasquared;



Yregression = (xhat_selected_big(2:N) - xhat_selected_big(1:N-1))./sqrt(xhat_selected_big(1:N-1));
Zregression = h * [Dose*Ka*exp(-Ka*owntime(1:N-1)')./sqrt(xhat_selected_big(1:N-1)), -sqrt(xhat_selected_big(1:N-1))];
%S_ClKe = inv(Zregression'*Zregression)*Zregression'*Yregression;  % Cl first then Ke ... a 2x1 vector
% here we use a more stable matrix inversion method instead of inv(Zregression'*Zregression)
S_ClKe = (Zregression'*Zregression)\(Zregression'*Yregression);  % Ke/Cl

% Here follows a constrained optimization based version.
%options = optimoptions('lsqlin','MaxIter',10000);
%startparam = [Cl Ke];
%lb = [1e-3;1e-3];
%S_ClKe = lsqlin(Zregression,Yregression,[],[],[],[],lb,[],startparam,options);

S_Ke = S_ClKe(2);
S_KeClratio = S_ClKe(1);

% if S_Ke<=0 || S_KeClratio<=0
%     S_KeClratio
%     S_Ke
% end

%sufficientstats = [S_Ke,S_Cl,S_sigma,S_sigmaepsilon];
%sufficientstats = [S_Ke,S_Cl,S_sigma];
sufficientstats = [S_KeClratio,S_Ke,S_sigmasquared,S_sigmaepsilonsquared];


end

