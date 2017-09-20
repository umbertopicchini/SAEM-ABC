function [THETAmatrix_saem,Fisher_info,xhat_selected] = saem_abc(model_param,parmask,parbase,yobs,saem_numit,warmup,fisherestim_iter,numparticles,N_threshold,abc_schedule,abc_vector,verbose)
% the SAEM-ABC parameter estimation algorithm
% See Picchini & Samson 2017, Coupling stochastic EM and Approximate Bayesian Computation for parameter inference in state-space models, arXiv:1512.04831.

% retrieve model specific settings
bigtheta = model_param{1};
problem = model_param{2};
time    = model_param{3};
numdepvars = model_param{4};
vrbl    = model_param{5};

% extract FREE parameters (to be estimated)
theta = param_mask(bigtheta,parmask);
fprintf('\n SAEM parameters starting values are:');
fprintf('\n %d',theta);
fprintf('\n')


THETAmatrix_saem = zeros(saem_numit+1,length(theta));
THETAmatrix_saem(1,:) = theta;

stats_old = zeros(1,length(theta));
G_vec = zeros(length(theta),1);
Hessian_matr = zeros(length(theta),length(theta));


for saem_iter = 2:saem_numit+1
     % set the coefficients denoted as gamma_k in Picchini and Samson (equation 3)
     if saem_iter < warmup
         alpha_sequence = 1;
     else
         alpha_sequence = 1/(saem_iter-warmup+1);
     end
 
   [ok,index] = ismember(saem_iter-1,abc_schedule);
   if ok  % update (i.e. decrease) the ABC tolerance
       abc_tolerance = abc_vector(index);
       fprintf('\n*** Iteration #%d, now using ABC tolerance = %d ***',saem_iter,abc_tolerance)
   end
       
   % now produce a state trajectory using the ABC filter
   xhat_selected = abcsmc_filter(model_param,yobs,numparticles,N_threshold,abc_tolerance,verbose);
   if verbose
      fprintf('\n SAEM iteration %d completed. ',saem_iter);
   end
   % produce a vector of sufficient summary statistics for the unknown
   % parameters
   saem_stats = feval([problem, '_saemstats'],xhat_selected,yobs,bigtheta,time);

   stats_new = stats_old + alpha_sequence*(saem_stats-stats_old);
   stats_old = stats_new;
   
   % this is not really needed
   [gradient,secondderivatives] = feval([problem, '_saem_gradienthessian'],xhat_selected,yobs,bigtheta,time);
  
   % M-STEP of SAEM (only for parameters to be stimated)
   theta = feval([problem, '_saem_parupdate'],bigtheta,stats_new,time);
   THETAmatrix_saem(saem_iter,:) = theta;
   % recover the full vector of structural model parameters (known +
   % unknown parameters)
   bigtheta = param_unmask(theta,parmask,parbase);
   model_param{1} = bigtheta;
   
   if saem_iter >=fisherestim_iter
      % THIS IS NOT NECESSARY FOR THE ESTIMATION PROCEDURE
      G_vec = G_vec + alpha_sequence*(gradient - G_vec);
      Hessian_matr = Hessian_matr + alpha_sequence*(secondderivatives + gradient*gradient' - Hessian_matr);   
      Fisher_info = Hessian_matr - G_vec*G_vec';
   end
end

save('THETAmatrix_saem','THETAmatrix_saem')


fprintf('\n')

end

