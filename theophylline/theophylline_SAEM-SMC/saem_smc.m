function THETAmatrix_saem = saem_smc(model_param,parmask,parbase,yobs,initstate,saem_numit,warmup,fisherestim_iter,numparticles,N_threshold)
%the SAEM-SMC algorithm in Picchini & Samson 2017, Coupling stochastic EM and Approximate Bayesian Computation for parameter inference in state-space models, arXiv:1512.04831.

bigtheta = model_param{1};
problem = model_param{2};
owntime    = model_param{3};
time    = model_param{4};
numdepvars = model_param{5};
vrbl    = model_param{6};

theta = param_mask(bigtheta,parmask);
fprintf('\n SAEM parameters starting values are:');
fprintf('\n %d',theta);
fprintf('\n')


THETAmatrix_saem = zeros(saem_numit+1,length(theta));
THETAmatrix_saem(1,:) = theta;

stats_old = zeros(1,length(theta));


for saem_iter = 2:saem_numit+1
    
     if saem_iter < warmup
         alpha_sequence = 1;
     else
         alpha_sequence = 1/(saem_iter-warmup+2)^0.95;
     end
 

    % now produce a state trajectory using the ABC-SMC sampler
   [xhat_selected_big,xhat_selected_small] = smc_filter(model_param,yobs,initstate,numparticles,N_threshold);
   
   fprintf('\n SAEM iteration %d completed. ',saem_iter);
   % produce a vector of sufficient summary statistics for the unknown
   % parameters
   saem_stats = feval([problem, '_saemstats2'],xhat_selected_big,xhat_selected_small,yobs,bigtheta,owntime,time);

   stats_new = stats_old + alpha_sequence*(saem_stats-stats_old);
   stats_old = stats_new;
   
   % M-STEP of SAEM (only for parameters to be stimated)
   theta = feval([problem, '_saem_parupdate'],bigtheta,stats_new,time,owntime);
   THETAmatrix_saem(saem_iter,:) = theta;
   % recover the full vector of structural model parameters (known +
   % unknown parameters)
   bigtheta = param_unmask(theta,parmask,parbase);
   model_param{1} = bigtheta;
   
end


save('THETAmatrix_saem','THETAmatrix_saem')


fprintf('\n')

end

