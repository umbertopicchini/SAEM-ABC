function [xhat_selected_big,xhat_selected_small] = smc_filter_gw(model_param,yobs,initstate,numparticles,N_threshold)
% SMC filtering using the Golightly-Wilkinson (GW) sampler.
% See Golightly, Andrew, and Darren J. Wilkinson. "Bayesian parameter inference for stochastic biochemical network models using particle Markov chain Monte Carlo." Interface focus 1.6 (2011): 807-820.

        bigtheta = model_param{1};
        problem = model_param{2};
        owntime    = model_param{3};
        time    = model_param{4};
        numdepvars = model_param{5};
        vrbl    = model_param{6};
        
        log_sigmay = bigtheta(6);
        sigmay = exp(log_sigmay);

        stepsize = owntime(2)-owntime(1);
        time_id = 1;
        nobs = length(time);
        xhat_selected_small = zeros(nobs,1);

        startstate = initstate;  % this is the starting state for the filter, NOT x0 (the first value in bigtheta)
        %use the Golightly-Wilkinson sampler
        [xhat,sum_logtransdensity,sum_logdensitysampler] = feval([problem, '_gw_sampling'],bigtheta,[owntime(1):stepsize:time(1)],time_id,[],numparticles,startstate,yobs(time_id));
        logweights =  -log(sigmay) -(yobs(time_id)-xhat(end,:)).^2/(2*sigmay^2) + sum_logtransdensity - sum_logdensitysampler;
        
        logweights = (logweights-max(logweights))-log(sum(exp(logweights-max(logweights)))); % normalize weights; suggestion from page 6 of Cappe et al. "An overview of existing methods and recent advances in sequential Monte Carlo"
        ess = 1/sum((exp(logweights)).^2);  % the Effective Sample Size
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
        num_diff_particle(1) = length(unique(particle_indeces));
        
        xinit_shuffled = xhat(end,particle_indeces);
       
        for ii=1:nobs-1
            time_id = ii+1;
            %use the Golightly-Wilkinson sampler
            [xhat,sum_logtransdensity,sum_logdensitysampler] = feval([problem, '_gw_sampling'],bigtheta,[time(ii):stepsize:time(ii+1)],time_id,xinit_shuffled,numparticles,[],yobs(time_id));
            logweights = logweights - log(sigmay) -(yobs(time_id)-xhat(end,:)).^2/(2*sigmay^2) + sum_logtransdensity - sum_logdensitysampler;            
            logweights = (logweights-max(logweights))-log(sum(exp(logweights-max(logweights)))); % normalize weights; suggestion from page 6 of Cappe et al. "An overview of existing methods and recent advances in sequential Monte Carlo"
            id_selected = stratresample(exp(logweights),1); % sample a single index from the set of probabilities
            xhat_selected_big = [xhat_selected_big; xhat(:,id_selected)];
            xhat_selected_small(ii+1) = xhat(end,id_selected);
            
            
            ess = 1/sum((exp(logweights)).^2);  % the Effective Sample Size
            if ess < N_threshold && ii<nobs %(we do not want to reset weights nor resample if ii==nobs) 
               resampling_count = resampling_count +1;
               particle_indeces = stratresample(exp(logweights),numparticles);  % Resampling. Notice the function will take care of normalizing weights
               logweights = -log(numparticles) *ones(1,numparticles);
            else
               particle_indeces  = [1:numparticles];
            end
            xinit_shuffled = xhat(end,particle_indeces);
        end
        fprintf('\nResampled %d times out of %d. Final ESS %d',resampling_count,nobs,ess);
        

end

