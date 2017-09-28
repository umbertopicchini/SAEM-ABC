% runs 30 estimation procedures using SAEM-ABC for the example in section 4.1
% in Picchini & Samson 2017, Coupling stochastic EM and Approximate Bayesian Computation for parameter inference in state-space models, arXiv:1512.04831.
% in particular, the settings below reproduce the SAEM-ABC results in the second column of Table 1

rng(10)  % set a random numbers seed for reproducibility


warmup = 300; % K1, the number of SAEM iterations before using the (descreasing) "gamma" sequence into SAEM 
saem_numit = 400; % K, the total number of SAEM iterations
numparticles = 1000;  % M, the number of particles
N_threshold = 200; % the resampling threshold (\bar{M} in the paper)
fisherestim_iter = 300; % ignore, this is not used in the paper.
abc_schedule = [1 80 150 200]; % interpret this together with abc_vector below: from SAEM iterations between 1 and 80 we use delta = 2,
                               % from iteration 80 to iteration 150 we use delta=1.7,
                               % from iteration 150 to iteration 200 we use delta=1.3,
                               % from iteration 200 to iteration saem_numit we use delta=1.
abc_vector = [2 1.7 1.3 1];  % the DECREASING sequence of ABC thresholds delta
verbose = 1;

% Model Setup
problem = 'nonlingauss';      % a string identifying the given problem at hand. 
numdepvars = 1;               % the number of modelled "state variables"
integrator = [];              % keep this empty, we do not perform numerical discretization of the state process here
sampletime = [1:1:50];        % the observational times
time = sampletime';  % MUST BE A COLUMN VECTOR              
owntime = time;  % NOT NECESSARY, WILL NOT BE USED



fprintf('\n\nCOMPUTING...\n\n');



%:::::::::::: HERE GO TRUE PARAMETER VALUES FOR DATA GENERATION :::::::::::::::
X0 = 0;
sigmax = sqrt(5);  % sigmax = 2.24
sigmay = sqrt(5);  % sigmay = 2.24
log_sigmax = log(sigmax); % log_sigmax = 0.8
log_sigmay = log(sigmay); % log_sigmay = 0.8


vrbl = [1:numdepvars];
% Now replicate the values in vrbl depending on the time value 
vrbl = repmat(vrbl,1,length(unique(time)));


%::::::::: GENERATE DATA ::::::::::::::::::::::::::::::::::::::::

bigtheta_true = [X0,log_sigmax,log_sigmay];   % store here all parameters 
model_param = {bigtheta_true,problem,time,numdepvars,vrbl};  

% create 'observations' by adding measurement noise to the noise-free trajectory
xhat_true = zeros(length(time),1);
yobs = zeros(length(time),1);
time_id = 1;
xhat_true(1) = feval([problem, '_statemodel'],bigtheta_true,time_id,[],1);
yobs(1) =  feval([problem, '_errormodel'],bigtheta_true,xhat_true(1),1);

for ii=2:length(time)
    time_id = ii;
    xhat_true(ii) = feval([problem, '_statemodel'],bigtheta_true,time_id,xhat_true(ii-1),1);
    yobs(ii) = feval([problem, '_errormodel'],bigtheta_true,xhat_true(ii),1);
end

% % plot data
% plot(time,yobs,'ko-')
% return

%:::


%::::::::: END OF DATA GENERATION ::::::::::::::::::::::::::::::::::::

% starting values for SAEM-ABC
%                 X0   log_sigmax    log_sigmay
bigtheta_start = [0,      0.8,           0.8];    % mean parameter values, see below
parmask        = [0 ,       1,            1 ];    % 0 for fixed parameter, 1 for parameters to be estimated 
parbase = bigtheta_start;

numattempts = 30;  % number of attempted estimation procedures
% starting parameter values (log-scale), randomly sampled from a Gaussian with mean
% bigtheta_mean and diagonal variance matrix with diagonal [2 2]
theta_start_matrix = mvnrnd(bigtheta_start(2:3),[2 2],numattempts);
bigtheta_start_matrix = [zeros(numattempts,1),theta_start_matrix];
all_theta_estimated = zeros(numattempts,sum(parmask));

tic
for ii=1:numattempts
    
    model_param{1} = bigtheta_start_matrix(ii,:);
    parbase = bigtheta_start_matrix(ii,:);

   [THETAmatrix_saem,Fisher_info,xhat_selected] = saem_abc(model_param,parmask,parbase,yobs,saem_numit,warmup,fisherestim_iter,numparticles,N_threshold,abc_schedule,abc_vector,verbose);

   all_theta_estimated(ii,:) = THETAmatrix_saem(end,:);
 
    if ii==1
        figure
    end
   % exponentiate so we transform log-parameters to their natural scale
   subplot(1,2,1)
   plot(exp(THETAmatrix_saem(:,1)))
   xlabel('\sigma_x')
   hline(sigmax)  % line showing the true parameter value
   hold on
   subplot(1,2,2)
   plot(exp(THETAmatrix_saem(:,2)))
   xlabel('\sigma_y')
   hline(sigmay)  % line showing the true parameter value
   hold on
   saveas(gcf,'saem_trajectories')
end
eval = toc
save('all_theta_estimated','all_theta_estimated')
hold off
% exponentiate parameters to go from log-scale to natural scale
figure
subplot(1,2,1)
hist(exp(all_theta_estimated(:,1)))
subplot(1,2,2)
hist(exp(all_theta_estimated(:,2)))
fprintf('\nfirst, second and third quartile for sigma_x')
quantile(exp(all_theta_estimated(:,1)),[0.25,0.5,0.75])
fprintf('\nfirst, second and third quartile for sigma_y')
quantile(exp(all_theta_estimated(:,2)),[0.25,0.5,0.75])



