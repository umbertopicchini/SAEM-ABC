# runs the pseudo-marginal Bayesian procedure for the experiment in section 4.1.2 of the reference below:
# Picchini & Samson 2017, Coupling stochastic EM and Approximate Bayesian Computation for parameter inference in state-space models, arXiv:1512.04831.
# notice results were obtained using v.1.4.1.1 of the pomp package

# Notice this script uses the same data and the same parameter starting values as in the MATLAB scripts for SAEM-ABC and SAEM-SMC


library(pomp)
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

# the observational density
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

densprior=Csnippet("
    lik = dunif(sigmax,0.1,15,1)+dunif(sigmay,0.1,15,1);
                lik = (give_log) ? lik : exp(lik);
                ")

nlingauss <- pomp(data.dat,time="time",t0=0,rprocess=discrete.time.sim(nlingauss.step,delta.t=1),
                  dmeasure=dmeas,dprior=densprior,statenames="x",paramnames=c("sigmax","sigmay"),
                  toEstimationScale=logtrans,fromEstimationScale=exptrans)

theta.true<-c(sigmax=sqrt(5),sigmay=sqrt(5),x.0=0)
estpars <- c("sigmax","sigmay")  # we only make inference for sigmax and sigmay
theta.guess <- theta.true


# try to achieve a 7% acceptance rate, which is supposed to be optimal for pmcmc algorithms with 
# perturbations on parameters driven by Gaussian random walk, see Sherlock et al. 2015, Annals of Statistics
pmcmc1 <- pmcmc(nlingauss,Nmcmc=4000,Np=2000,start=theta.guess,
                proposal=mvn.rw.adaptive(rw.sd=c(sigmax=0.1,sigmay=0.1),target = 0.07,scale.start=100,shape.start=100)) 
library(coda)
pmcmc1 %>% conv.rec() -> traces
plot(traces)
mean(traces[2000:4000,4])  # sigmax
quantile(traces[2000:4000,4],c(0.025,0.975)) # sigmax
mean(traces[2000:4000,5])  # sigmay
quantile(traces[2000:4000,5],c(0.025,0.975)) # sigmay
# from the above is clear that at about iteration 1000 we approach convergence
# lets perform more iterations using the parameters covariance obtained from above
#pmcmc1 %>% pmcmc(Nmcmc=2000,proposal=mvn.rw(covmat(pmcmc1)))
#pmcmc1 %>% conv.rec() -> traces
