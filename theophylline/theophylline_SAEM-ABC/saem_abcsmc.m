function THETAmatrix_saem = saem_abcsmc(model_param,parmask,parbase,yobs,initstate,saem_numit,warmup,fisherestim_iter,numparticles,N_threshold,abc_schedule,abc_vector,verbose)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
         alpha_sequence = 1/(saem_iter-warmup+1);
     end
     
     [ok,index] = ismember(saem_iter-1,abc_schedule);
     if ok
       abc_tolerance = abc_vector(index);
       fprintf('\n*** Iteration #%d, now using ABC tolerance = %d ***',saem_iter,abc_tolerance)
    end
 

   [xhat_selected_big,xhat_selected_small]  = abcsmc_filter(model_param,yobs,initstate,numparticles,N_threshold,abc_tolerance,verbose);
   
   if verbose
      fprintf('\n SAEM iteration %d completed. ',saem_iter);
   end
   
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


fprintf('\n')

end

