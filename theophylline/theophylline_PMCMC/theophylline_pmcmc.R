# runs the pseudo-marginal Bayesian procedure for the experiment in section 4.3 of the reference below:
# Picchini & Samson 2017, Coupling stochastic EM and Approximate Bayesian Computation for parameter inference in state-space models, arXiv:1512.04831.
# notice results were obtained using v.1.4.1.1 of the pomp package

# Notice results from SAEM-ABC and SAEM-SMC have been obtained on 50 different datasets. The Bayesian procedure below
# works on a single data set. Hence the result of the Bayesian procedure is not directly comparable with results from SAEM-ABC and SAEM-SMC.


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

#data <-read.table('data.dat',header=T)
#attach(data)
data.dat <- read.csv(text="
                      time,y
                     1.0000,   11.6661
                     2.0000,   11.8163
                     3.0000,   11.0020
                     4.0000,   10.4717
                     5.0000,    9.8763
                     6.0000,    9.7701
                     7.0000 ,   9.0750
                     8.0000,    8.3915
                     9.0000,    7.9834
                     10.0000,    7.7654
                     11.0000,    7.6254
                     12.0000,    7.1055
                     13.0000,    7.6130
                     14.0000,    7.1441
                     15.0000,    6.4846
                     16.0000,    5.5977
                     17.0000,    5.7072
                     18.0000,    5.5185
                     19.0000,    4.9659
                     20.0000,    4.6396
                     21.0000,    4.8931
                     22.0000,    4.4881
                     23.0000,    4.2300
                     24.0000,    4.1558
                     25.0000,    3.4457
                     26.0000,    3.1985
                     27.0000,    2.8334
                     28.0000,    2.4181
                     29.0000,    2.4155
                     30.0000 ,   2.4749
                     31.0000,    2.2192
                     32.0000,    2.3229
                     33.0000,    2.5917
                     34.0000,    2.3822
                     35.0000,    2.2993
                     36.0000,    2.5573
                     37.0000,    2.2876
                     38.0000,    1.9056
                     39.0000,    1.9657
                     40.0000 ,   1.8270
                     41.0000,    1.9023
                     42.0000,    1.5322
                     43.0000,    1.3206
                     44.0000,    1.1870
                     45.0000,    1.1325
                     46.0000 ,   1.0559
                     47.0000,    0.9806
                     48.0000,    1.0880
                     49.0000,    0.9145
                     50.0000,    0.8110
                     51.0000,    0.5750
                     52.0000,    0.7885
                     53.0000 ,   1.0197
                     54.0000,    0.8659
                     55.0000,    0.5663
                     56.0000,    0.5942
                     57.0000,    0.6952
                     58.0000 ,   0.4462
                     59.0000,    0.5617
                     60.0000,    0.4589
                     61.0000,    0.2535
                     62.0000,    0.2706
                     63.0000,    0.3573
                     64.0000,    0.2909
                     65.0000,    0.2830
                     66.0000 ,   0.1156
                     67.0000,    0.0602
                     68.0000,    0.0264
                     69.0000,    0.0122
                     70.0000,    0.0629
                     71.0000,    0.0322
                     72.0000,    0.0953
                     73.0000,    0.1107
                     74.0000,    0.0718
                     75.0000,    0.1430
                     76.0000 ,  -0.0605
                     77.0000,    0.0359
                     78.0000,    0.0807
                     79.0000 ,   0.1118
                     80.0000,    0.0525
                     81.0000,    0.0272
                     82.0000 ,   0.0823
                     83.0000,    0.0724
                     84.0000,    0.0664
                     85.0000 ,  -0.1317
                     86.0000,   -0.1023
                     87.0000,   -0.0623
                     88.0000 ,  -0.0198
                     89.0000  ,  0.0639
                     90.0000,   -0.1204
                     91.0000,   -0.1131
                     92.0000,    0.0300
                     93.0000,   -0.1724
                     94.0000 ,   0.1341
                     95.0000,    0.0246
                     96.0000,   -0.0417
                     97.0000,    0.0148
                     98.0000,   -0.0473
                     99.0000,    0.0282
                     100.0000,   -0.0633"
)


# plot data
ggplot(data=data.dat,mapping=aes(x=time,y=y))+
  geom_line()+geom_point()+
  expand_limits(y=0)+
  theme_classic()

# state equation
step.fun <- Csnippet("
  double dW = rnorm(0,sqrt(dt));
  x += (Dose*Ka*Ke/Cl*exp(-Ka*t)-Ke*x)*dt+sigma*sqrt(abs(x))*dW;
  ")


# the observatiosal density
dmeas <- Csnippet("
  lik = dnorm(y,x,sigmaepsilon,give_log);
                  ")

logtrans <- Csnippet("
            TKe = log(Ke);
            TCl = log(Cl);
            Tsigma = log(sigma);
            Tsigmaepsilon = log(sigmaepsilon);
                     ")

exptrans <- Csnippet("
                     TKe = exp(Ke);
                     TCl = exp(Cl);
                     Tsigma = exp(sigma);
                     Tsigmaepsilon = exp(sigmaepsilon);
                     ")

densprior=Csnippet("
    lik = dunif(Ke,0.01,0.2,1)+dunif(Cl,0.01,0.2,1)+dunif(sigma,0.01,0.3,1)+dunif(sigmaepsilon,0.05,0.5,1);
                lik = (give_log) ? lik : exp(lik);
                ")

theop <- pomp(data.dat,time="time",t0=0,rprocess=euler.sim(step.fun=step.fun,delta.t=0.05),
              dmeasure=dmeas,dprior=densprior,statenames="x",paramnames=c("Dose","Ka","Ke","Cl","sigma","sigmaepsilon"),
              toEstimationScale=logtrans,fromEstimationScale=exptrans)

# iterated filtering
theta.true<-c(Dose=4,Ka=exp(0.4),Ke=exp(-3),Cl=exp(-3.22),sigma=0.1,sigmaepsilon=0.1,x.0=8)
estpars <- c("Ke","Cl","sigma","sigmaepsilon")  
theta.guess <- theta.true
theta.guess[estpars]  <- c(Ke=0.05,Cl=0.04,sigma=0.2,sigmaepsilon=0.3)

# try to achieve a 7% acceptance rate, which is supposed to be optimal for pmcmc algorithms with 
# perturbations on parameters driven by MRW, Sherlock et al. 2015, Annals of Statistics
pmcmc1 <- pmcmc(theop,Nmcmc=2000,Np=1000,start=theta.guess,
                proposal=mvn.rw.adaptive(rw.sd=c(Ke=0.005,Cl=0.005,sigma=0.01,sigmaepsilon=0.01),target = 0.07)) 
library(coda)
pmcmc1 %>% conv.rec() -> traces
plot(traces[,6])
plot(traces[,7])
plot(traces[,8])
plot(traces[,9])
mean(traces[1000:2000,6])  # Ke
quantile(traces[1000:2000,6],c(0.025,0.975)) # Ke
mean(traces[1000:2000,7])  # Cl
quantile(traces[1000:2000,7],c(0.025,0.975)) # Cl
mean(traces[1000:2000,8])  # sigma
quantile(traces[1000:2000,8],c(0.025,0.975)) # sigma
mean(traces[1000:2000,9])  # sigma_epsilon
quantile(traces[1000:2000,9],c(0.025,0.975)) # sigma_epsilon
