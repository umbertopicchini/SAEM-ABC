function [xhat_selected_big,xhat_selected_small]  = abcsmc_filter(model_param,yobs,initstate,numparticles,N_threshold,abc_tolerance,verbose)
%the ABC filter for state estimation

        bigtheta = model_param{1};
        problem = model_param{2};
        owntime    = model_param{3};
        time    = model_param{4};
        numdepvars = model_param{5};
        vrbl    = model_param{6};

        stepsize = owntime(2)-owntime(1);
        time_id = 1; 
        nobs = length(time);
        xhat_selected_small = zeros(nobs,1);
        startstate = initstate;  % this is the starting state for the filter, NOT x0 (the first value in bigtheta)
        xhat = feval([problem, '_statemodel'],bigtheta,[owntime(1):stepsize:time(1)],time_id,[],numparticles,startstate);
        
        yobssim = feval([problem, '_errormodel'],bigtheta,xhat(end,:),numparticles);

        logweights =  -log(abc_tolerance) -(yobs(time_id)-yobssim).^2/(2*abc_tolerance^2);
        logweights = (logweights-max(logweights))-log(XSum(exp(logweights-max(logweights)))); % normalize weights; suggestion from page 6 of Cappe et al. "An overview of existing methods and recent advances in sequential Monte Carlo"
        ess = 1/XSum((exp(logweights)).^2);
        resampling_count = 0;
        
        id_selected = stratresample(exp(logweights),1); % sample a single index from the set of probabilities
        xhat_selected_big = xhat(:,id_selected);
        xhat_selected_small(1) = xhat(end,id_selected);
   
        % in this example we assume the initial state
        % to be known and equal to X0 for each particle
        if ess < N_threshold
           resampling_count = resampling_count +1;
           particle_indeces = stratresample(exp(logweights),numparticles);  % Resampling. Notice the function will take care of normalizing weights
           logweights = -log(numparticles) *ones(1,numparticles);
        else
           particle_indeces  = [1:numparticles];
        end
        
        xinit_shuffled = xhat(end,particle_indeces);
        
        for ii=1:nobs-1
            time_id = ii+1;
            xhat = feval([problem, '_statemodel'],bigtheta,[time(ii):stepsize:time(ii+1)],time_id,xinit_shuffled,numparticles,[]);
            yobssim = feval([problem, '_errormodel'],bigtheta,xhat(end,:),numparticles);
            logweights = logweights - log(abc_tolerance) -(yobs(time_id)-yobssim).^2/(2*abc_tolerance^2);
            logweights = (logweights-max(logweights))-log(XSum(exp(logweights-max(logweights)))); % normalize weights; suggestion from page 6 of Cappe et al. "An overview of existing methods and recent advances in sequential Monte Carlo"
            
            id_selected = stratresample(exp(logweights),1); % sample a single index from the set of probabilities
            xhat_selected_big = [xhat_selected_big; xhat(:,id_selected)];
            xhat_selected_small(ii+1) = xhat(end,id_selected);
            ess = 1/XSum((exp(logweights)).^2);  % the Effective Sample Size
            if ess < N_threshold && time_id<nobs %(we do not want to reset weights nor resample if time_id==nobs) 
               resampling_count = resampling_count +1;
               particle_indeces = stratresample(exp(logweights),numparticles);  % Resampling. Notice the function will take care of normalizing weights
               logweights = -log(numparticles) *ones(1,numparticles);
            else
               particle_indeces  = [1:numparticles];
            end
            xinit_shuffled = xhat(end,particle_indeces);
        end
        if verbose
           fprintf('\nResampled %d times out of %d. Final ESS %d',resampling_count,nobs,ess);
        end
   
end

