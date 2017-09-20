function xhat_selected = smc_filter(model_param,yobs,numparticles,N_threshold)
%the SMC filter for state estimation: uses a boostrap filter

        bigtheta = model_param{1};
        problem = model_param{2};
        time    = model_param{3};
        numdepvars = model_param{4};
        vrbl    = model_param{5};
        
        log_sigmay = bigtheta(3);
        sigmay = exp(log_sigmay);

        
        time_id = 1;
        nobs = length(time);

        xhat = feval([problem, '_statemodel'],bigtheta,time_id,[],numparticles);
        logweights =  -log(sigmay) -(yobs(time_id)-xhat).^2/(2*sigmay^2);
        logweights = (logweights-max(logweights))-log(sum(exp(logweights-max(logweights)))); % normalize weights; suggestion from page 6 of Cappe et al. "An overview of existing methods and recent advances in sequential Monte Carlo"
        ess = 1/sum((exp(logweights)).^2);  % the Effective Sample Size
        resampling_count = 0;
   
        % in this example we assume the initial state
        % to be known and equal to X0 for each particle
        if ess < N_threshold
           resampling_count = resampling_count +1;
           particle_indeces = stratresample(exp(logweights),numparticles);  % Resampling. Notice the function will take care of normalizing weights
           logweights = -log(numparticles) *ones(1,numparticles);
        else
           particle_indeces  = [1:numparticles];
        end

        xinit_shuffled = xhat(particle_indeces);
        
        genealogy_indeces = zeros(nobs-1,numparticles);
        genealogy_indeces(1,:) = particle_indeces;
        
        genealogy_states = zeros(nobs,numparticles);
        genealogy_states(1,:)=xinit_shuffled;
        
        for ii=2:nobs
            time_id = ii;
            xhat = feval([problem, '_statemodel'],bigtheta,time_id,xinit_shuffled,numparticles);
            logweights = logweights - log(sigmay) -(yobs(time_id)-xhat).^2/(2*sigmay^2);
            logweights = (logweights-max(logweights))-log(sum(exp(logweights-max(logweights)))); % normalize weights; suggestion from page 6 of Cappe et al. "An overview of existing methods and recent advances in sequential Monte Carlo"
            ess = 1/sum((exp(logweights)).^2);  % the Effective Sample Size
            if ess < N_threshold && ii<nobs %(we do not want to reset weights nor resample if ii==nobs) 
               resampling_count = resampling_count +1;
               particle_indeces = stratresample(exp(logweights),numparticles);  % Resampling. Notice the function will take care of normalizing weights
               logweights = -log(numparticles) *ones(1,numparticles);
            else
               particle_indeces  = [1:numparticles];
            end
            xinit_shuffled = xhat(particle_indeces);
            genealogy_states(ii,:)=xinit_shuffled;
            if ii<nobs
               genealogy_indeces(ii,:)=particle_indeces;
            end
        end
        fprintf('\nResampled %d times out of %d. Final ESS %d',resampling_count,nobs,ess);
        
        xhat_selected = zeros(length(time),1);        

        % reconstruct the genealogy of particles ancestors
        % Start by selecting an index from {1,2,...,numparticles} as in
        % Particle Marginal Metropolis-Hastings (Andrieu et al 2010)
        id_selected = stratresample(exp(logweights),1); % sample a single index from the last set of probabilities
        % now start filling  the selected vector backwards
        xhat_selected(nobs) = genealogy_states(nobs,id_selected);
        for ii= nobs-1:-1:1
           index = genealogy_indeces(ii,id_selected);
           xhat_selected(ii) = genealogy_states(ii,index);
           id_selected = index;
        end

end

