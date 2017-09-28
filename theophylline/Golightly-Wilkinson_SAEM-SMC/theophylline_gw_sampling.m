function [xhat_next,sum_logtransdensity,sum_logdensitysampler] = theophylline_gw_sampling(bigtheta,subowntime,time_id,xhat_pre,numparticles,startstate,yobs_current)
% particles propagation via the the Golightly-Wilkinson (GW) sampler.
% See Golightly, Andrew, and Darren J. Wilkinson. "Bayesian parameter inference for stochastic biochemical network models using particle Markov chain Monte Carlo." Interface focus 1.6 (2011): 807-820.

% Parameters
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

Dose = 4; % the adminstered drug dose

N= length(subowntime);
xhat_matrix = zeros(N-1,numparticles);

if isempty(startstate) % startstate is the starting state for the used SMC filter, not to be confused with xzero
   init_state = xzero;
else
   init_state = startstate;
end

sum_logdensitysampler = 0;
sum_logtransdensity = 0;

for j=2:N
    
    if time_id==1 && j==2  % here xhat_pre is empty
        xhat_pre = init_state*ones(1,numparticles);
        t = subowntime(j-1);
        if t<0
            error('Write OWNTIME in such a way that ONWTIME(1) is strictly positive.')
        end
    end

    h = subowntime(j)- subowntime(j-1);
    t = subowntime(j-1);

    driftX = (Dose * Ka*Ke)/Cl * exp(-Ka*t) -Ke * xhat_pre ;          % the Ito SDE drift
    diffusionX = sigma*sqrt(abs(xhat_pre));      % the Ito SDE diffusion
    
    % construct the alpha and beta coefficients to conform with the
    % notation in Golightly and Wilkinson 2011, Interface Focus 1, 807-820.
    % see beginning of section 4.2
    alpha = driftX;
    beta = diffusionX.^2;
    
    % sample as from equation 4.8 in Golightly and Wilkinson 2011, Interface Focus 1, 807-820.
    % coefficients a and b are in the appendix A.3
    delta_j = subowntime(end) - t;
    a = alpha + beta./(beta*delta_j + sigmaepsilon^2) .* (yobs_current - (xhat_pre + alpha*delta_j));
    b = beta - beta./(beta*delta_j + sigmaepsilon^2) .* beta * h;
    % now sample as from equation 4.8 in Golightly and Wilkinson 2011
    xhat_new = xhat_pre + a*h + sqrt(b*h).*randn(1,numparticles);
    
    xhat_matrix(j-1,:) = xhat_new;
    
    % here is the sum of log-densities induced by the G-W sampler
    sum_logdensitysampler = sum_logdensitysampler -1/2*log(b*h) - 1./(2*b*h) .* (xhat_new-xhat_pre-a*h).^2;
    % here is the sum of log-transition densities induced by Euler-Maruyama
    sum_logtransdensity = sum_logtransdensity -1/2*log(beta*h) - 1./(2*beta*h) .* (xhat_new-xhat_pre-alpha*h).^2;
    
    xhat_pre = xhat_new;
    
end

if  time_id==1 
   xhat_next = [init_state*ones(1,numparticles);abs(xhat_matrix)];
else
   xhat_next = abs(xhat_matrix);
end

end


