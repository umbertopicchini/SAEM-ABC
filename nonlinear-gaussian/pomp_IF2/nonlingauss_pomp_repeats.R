# runs the (improved) iterated filtering (IF2) procedure for the experiment in section 4.1.2 of the reference below:
# Picchini & Samson 2017, Coupling stochastic EM and Approximate Bayesian Computation for parameter inference in state-space models, arXiv:1512.04831.
# notice results were obtained using v.1.4.1.1 of the pomp package

# WARNING: THE IF2 PROCEDURE BELOW IS SLOW! It runs 30 estimations procedures just as in the cited paper. This can be 
# changed by setting a different value for "numattempts", see below.

# Notice this script uses the same data and the same parameter starting values as in the MATLAB scripts for SAEM-ABC and SAEM-SMC


library(pomp)  # we used v.1.4.1.1 to produce results
stopifnot(packageVersion("pomp")>="1.0.0.0")
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8",
  scipen=5
)
set.seed(594709947L)

## ----prelims2,echo=FALSE,cache=FALSE-------------------------------------
library(ggplot2)
library(plyr)
library(reshape2)
library(magrittr)
theme_set(theme_bw())


# here follows the same data that is simulated in the SAEM-ABC and SAEM-SMC folders
data.dat <- read.csv(text="
                      time,y
                 1.0000000e+00,   -1.4910201e-01	
                 2.0000000e+00,    1.9068885e+00	
                 3.0000000e+00,    4.6734116e-01	
                 4.0000000e+00,   -6.7779706e-01	
                 5.0000000e+00,   -3.0348450e+00	
                 6.0000000e+00,	   6.5399287e+00	
                 7.0000000e+00,	  -6.6112666e-01	
                 8.0000000e+00,	   2.5084616e+00	
                 9.0000000e+00,	   1.7287442e-01	
                 1.0000000e+01,	   5.7009581e+00	
                 1.1000000e+01,	  -1.3298102e+00	
                 1.2000000e+01,	   1.7562693e+00	
                 1.3000000e+01,	  -2.2904698e-02	
                 1.4000000e+01,	   2.5330131e+00	
                 1.5000000e+01,	   4.8754377e-01	
                 1.6000000e+01,	   1.8761293e+00	
                 1.7000000e+01,	  -4.4095573e-01	
                 1.8000000e+01,	  -2.2305301e+00	
                 1.9000000e+01,	   2.9528561e+00	
                 2.0000000e+01,	  -1.2645165e+00	
                 2.1000000e+01,	   2.0093663e+00	
                 2.2000000e+01,	   2.4434073e+00	
                 2.3000000e+01,	   1.7706251e+00	
                 2.4000000e+01,	   4.8372965e+00	
                 2.5000000e+01,	  -3.8144987e+00	
                 2.6000000e+01,	  -6.4620759e-01	
                 2.7000000e+01,	   1.7273159e+00	
                 2.8000000e+01,	   1.1849207e+00	
                 2.9000000e+01,	   2.8706871e+00	
                 3.0000000e+01,	  -7.8792890e-01	
                 3.1000000e+01,	   5.7410034e+00	
                 3.2000000e+01,	  -1.4632518e-01	
                 3.3000000e+01,	  -4.2673067e-01	
                 3.4000000e+01,	  -1.0232973e+00	
                 3.5000000e+01,	   1.1854848e+00	
                 3.6000000e+01,	  -8.8354174e-02	
                 3.7000000e+01,	   9.6616186e-01	
                 3.8000000e+01,	   2.3321992e+00	
                 3.9000000e+01,	   5.2167439e+00	
                 4.0000000e+01,	   1.8078763e+00	
                 4.1000000e+01,	   2.4969157e+00	
                 4.2000000e+01,	  -3.6976843e+00	
                 4.3000000e+01,	  -2.0314177e+00	
                 4.4000000e+01,	   9.1071866e-01	
                 4.5000000e+01,	   9.6122893e-01	
                 4.6000000e+01,	   1.5510186e+00	
                 4.7000000e+01,	  -2.0222503e+00	
                 4.8000000e+01,	   1.6072386e+00	
                 4.9000000e+01,	  -4.2592852e-01	
                 5.0000000e+01,	   5.0207822e+00"
)


# plot data
ggplot(data=data.dat,mapping=aes(x=time,y=y))+
  geom_line()+geom_point()+
  expand_limits(y=0)+
  theme_classic()

# state equation
nlingauss.step <- Csnippet("
  double eps = rnorm(0,sigmax);
  x = 2*sin(exp(x)) + eps;
  ")

# the deterministic skeleton (not sure if this is necessary!)
#nlingauss.skel <- Csnippet("
#                    Dx = 2*sin(exp(x));
#                    ")

# that's how to simulate from the observation model
#rmeas <- Csnippet("
#  y = rnorm(x,sigmay);
#                  ")

# the observatiosal density
dmeas <- Csnippet("
  lik = dnorm(y,x,sigmay,give_log);
                  ")

logtrans <- Csnippet("
            Tsigmax = log(sigmax);
            Tsigmay = log(sigmay);
                     ")

exptrans <- Csnippet("
                     Tsigmax = exp(sigmax);
                     Tsigmay = exp(sigmay);
                     ")

# construct the state-space model object via POMP
nlingauss <- pomp(data.dat,time="time",t0=0,rprocess=discrete.time.sim(nlingauss.step,delta.t=1),
                  dmeasure=dmeas,statenames="x",paramnames=c("sigmax","sigmay"),
                  toEstimationScale=logtrans,fromEstimationScale=exptrans)


# just produce a test of the above, to see if it works
#sim <- simulate(nlingauss,nsim=10,params=c(sigmax=sqrt(5),sigmay=sqrt(5),x.0=0),states=TRUE)
#traj <- trajectory(nlingauss,params=c(sigmax=sqrt(5),sigmay=sqrt(5),x.0=0))
#pf <- pfilter(nlingauss,Np=1000,params=c(sigmax=sqrt(5),sigmay=sqrt(5),x.0=0))

# iterated filtering
theta.true<-c(sigmax=sqrt(5),sigmay=sqrt(5),x.0=0)
estpars <- c("sigmax","sigmay")  # we only make inference for sigmax and sigmay
theta.guess <- theta.true

numattempts = 30

# these are the same starting values for the parameters as in the SAEM-ABC and SAEM-SMC folders.
# the (apparent difference) is that this R code optimizes sigmax and sigmay, while the MATLAB code 
# uses parameters on their log-scale, log_sigmax and log_sigmay. This means that by taking the natural logarithms
# of the values in theta_start_matrix we get the same starting values as used in the MATLAB scripts.
theta_start_matrix = matrix(c(    0.4778,   10.1808,
                                  1.8706,    4.8454,
                                  4.1900,    2.7098,
                                  6.0509,    2.4152,
                                  3.5876,    0.1399,
                                  0.1004,   12.4406,
                                  1.5329,    3.1448,
                                  2.5547,    4.8199,
                                  0.9399,    2.4406,
                                  2.2637,   23.5831,
                                  5.1497,    3.7482,
                                  0.8396,    3.9676,
                                  0.6215,    0.5988,
                                  0.7757,    1.4798,
                                  1.1044,    3.7894,
                                  0.5519,    0.4715,
                                  7.3608,    1.8551,
                                  1.3885,    0.3983,
                                  0.0555,    2.3181,
                                  4.7019,   17.0991,
                                  1.8776,    0.6765,
                                  5.0303,   26.8553,
                                  0.2483,    2.7346,
                                  0.4064,    4.0470,
                                  3.3812,    1.0913,
                                  0.3001,    1.7552,
                                  0.2990,    0.7597,
                                  2.7989,    0.0646,
                                  0.3150,    0.2505,
                                  1.2140,    1.7385),nrow=numattempts,ncol=2,byrow=T)

all_theta_estimated <- matrix(nrow=numattempts,ncol=2)

for (ii in 1:numattempts){
theta.guess[estpars]  <- theta_start_matrix[ii,]
m1 <- mif2(
             nlingauss,
             Nmif=500,
             start=theta.guess,
             transform=TRUE,
             rw.sd=rw.sd(sigmax=0.2,sigmay=0.2),
             cooling.fraction.50=0.9,
             Np=1000
             )
m1 <- continue(m1,Nmif=100,cooling.fraction=0.7)
m1 <- continue(m1,Nmif=100,cooling.fraction=0.4)
m1 <- continue(m1,Nmif=100,cooling.fraction=0.3)
m1 <- continue(m1,Nmif=100,cooling.fraction=0.2)

all_theta_estimated[ii,] <- coef(m1)[c(1, 2)]


}

save(file="all_theta_estimated.Rdata","all_theta_estimated")
