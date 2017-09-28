% runs 50 estimation procedures using SAEM-GW for the example in section 4.2
% in Picchini & Samson 2017, Coupling stochastic EM and Approximate Bayesian Computation for parameter inference in state-space models, arXiv:1512.04831.
% in particular, the settings below reproduce the SAEM-GW results in table
% when using M=200, \bar{M}=10.

numattempts = 50; % number of attempted estimation procedures
all_theta_estimated = zeros(numattempts,4);
all_mean_ess = zeros(numattempts,1);
allmean_num_diff_particle = zeros(numattempts,1);

% for each trial we change the seed of the random numbers generator, hence
% for each trail we simulate a different dataset, using the same parameter
% values. This is the only difference between the trials.

for trial=1:numattempts

rng(trial*100)

saem_numit = 300;  % K, the total number of SAEM iterations
warmup = 250; % K1, the number of SAEM iterations before using the (descreasing) "gamma" sequence into SAEM 
numparticles = 200; % M, the number of particles
N_threshold = 10; % the resampling threshold (\bar{M} in the paper)
fisherestim_iter = 200; % ignore, this is not used in the paper.

% Model Setup
problem = 'theophylline';            % a string identifying the given problem at hand. 
numdepvars = 1;                      % the number of modellized "state variables"
time = [1:1:100]';                   % the observational times, MUST BE A COLUMN VECTOR    
stepsize = (time(2)-time(1))/20;     % the SDE integration stepsize
owntime = [0:stepsize:time(end)];    % the discretised time interval for the numerical solution of the SDE

fprintf('\n\nCOMPUTING...\n\n');

%:::::::::::: HERE GO TRUE PARAMETER VALUES FOR DATA GENERATION :::::::::::::::
X0 = 8;
log_Ke = -3;  % Ke = 0.0498
log_Ka = 0.4;   % Ka = 1.49
log_Cl = -3.22;  % Cl = 0.04
sigma = 0.1;  % log_sigma = -2.3
sigmaepsilon = 0.1; % log_sigmaepsilon = -2.3


log_sigma = log(sigma);
log_sigmaepsilon = log(sigmaepsilon);


vrbl = [1:numdepvars];
% Now replicate the values in vrbl depending on the time value 
vrbl = repmat(vrbl,1,length(unique(time)));


%::::::::: GENERATE DATA ::::::::::::::::::::::::::::::::::::::::

bigtheta_true = [X0,log_Ke,log_Ka,log_Cl,log_sigma,log_sigmaepsilon];   % store here all parameters needed for SDE simulation

% simulate the latent state
xhat_true_full = [];
yobs = zeros(length(time),1);
time_id = 1; 
xhat_true = feval([problem, '_statemodel'],bigtheta_true,[owntime(1):stepsize:time(1)],time_id,[],1,[]);
xhat_true_full = [xhat_true_full;xhat_true];
% keep simulating the latent state
for ii=1:length(time)-1
    time_id = ii+1;
    xhat_true = feval([problem, '_statemodel'],bigtheta_true,[time(ii):stepsize:time(ii+1)],time_id,xhat_true(end,:),1,[]);
    xhat_true_full = cat(1,xhat_true_full,xhat_true);
end

% take the sequence of simulated states of the fine time-grid "owntime" the produce linearly interpolated states at observational time points
xhat_interpolated = interp1(owntime,xhat_true_full,time);

% add observational error
for ii = 1:length(time)
   yobs(ii) = feval([problem, '_errormodel'],bigtheta_true,xhat_interpolated(ii),1);
end
% 
% plot(time,yobs,'o-')
%return
% %::::::::: END OF DATA GENERATION ::::::::::::::::::::::::::::::::::::


%                 X0    log_Ke   log_Ka   log_Cl      log_sigma    log_sigmaepsilon
bigtheta_start = [8,      -0.22,      0.4,   2.3,             -2,          0   ];
parmask        = [0 ,      1,        0,       1,               1,            1  ];
parbase = bigtheta_start;

initstate = X0; 

   model_param = {bigtheta_start,problem,owntime,time,numdepvars,vrbl};  
   % produces inference via SAEM-GW
   THETAmatrix_saem = saem_gw(model_param,parmask,parbase,yobs,initstate,saem_numit,warmup,fisherestim_iter,numparticles,N_threshold);
   all_theta_estimated(trial,:) = THETAmatrix_saem(end,:);
 
   if trial==1
        figure
   end
   % exponentiate so we transform log-parameters to their natural scale
   subplot(2,2,1)
   plot(exp(real(THETAmatrix_saem(:,1))))
   xlabel('Ke')
   hline(exp(log_Ke))  % line showing the true parameter value
   hold on
   subplot(2,2,2)
   plot(exp(real(THETAmatrix_saem(:,2))))
   xlabel('Cl')
   hline(exp(log_Cl))  % line showing the true parameter value
   hold on
   subplot(2,2,3)
   plot(exp(real(THETAmatrix_saem(:,3))))
   xlabel('\sigma')
   hline(exp(log_sigma))  % line showing the true parameter value
   hold on
   subplot(2,2,4)
   plot(exp(real(THETAmatrix_saem(:,4))))
   xlabel('\sigma_{\epsilon}')
   hline(exp(log_sigmaepsilon))  % line showing the true parameter value
   hold on
   saveas(gcf,'saem_trajectories')
%   end
end
save('all_theta_estimated','all_theta_estimated')
hold off
% create boxplots for parameters on log-scale
% boxplot((all_theta_estimated))
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'XTickLabel',{'$\log Ke$','$\log Cl$','$\log \sigma$','$\log \sigma_\epsilon$'});


