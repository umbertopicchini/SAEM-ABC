function xhat_next = theophylline_statemodel(bigtheta,subowntime,time_id,xhat_pre,numparticles,startstate)


% Parameters
xzero     = bigtheta(1);
log_Ke    = bigtheta(2);
log_Ka    = bigtheta(3);
log_Cl    = bigtheta(4);
log_sigma = bigtheta(5);

Ke = exp(log_Ke);
Ka = exp(log_Ka);
Cl = exp(log_Cl);
sigma = exp(log_sigma);

Dose = 4; % the adminstered drug dose

N= length(subowntime);
xhat_matrix = zeros(N-1,numparticles);

if isempty(startstate) % startstate is the starting state for the used SMC filter, not to be confused with xzero
   init_state = xzero;
else
   init_state = startstate;
end

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
  
    Winc = sqrt(h)*randn(1,numparticles); % the Wiener increment(s) dWj;
    driftX = (Dose * Ka*Ke)/Cl * exp(-Ka*t) -Ke * xhat_pre ;          % the Ito SDE drift
    diffusionX = sigma*sqrt(abs(xhat_pre));      % the Ito SDE diffusion

    xhat_pre = abs(xhat_pre + driftX * h + diffusionX .* Winc);  
    
%     index_negative = xhat_pre<1e-3;
%     if any(index_negative)
%         xhat_pre(index_negative)=1e-3;
%     end
    xhat_matrix(j-1,:) = xhat_pre;
  
end

if  time_id==1 
   xhat_next = [init_state*ones(1,numparticles);abs(xhat_matrix)];
else
   xhat_next = abs(xhat_matrix);
end
%     if any(xhat_next<0) 
%         xhat_next
%         error('x<0')
%     end

end


