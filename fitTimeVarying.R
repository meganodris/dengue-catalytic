#----- Code to simulate & fit dengue time-varying FOI model -----#
rm(list=ls())
library(ggplot2)
library(tidyr)
library(rstan)
library(cmdstanr)

#--- simulate some data
simfoi <- function(nT, nA, lamH, lam, rho, gamma, pop, amin, amax){
  
  # immune profiles at beginning of time period
  age <- seq(0,99)
  susc0 <- exp(-4*lamH*age)
  mono0 <- 4*exp(-3*lamH*age)*(1-exp(-lamH*age))
  multi0 <- 1 - susc0 - mono0
  
  # immune profiles
  susc <- mono <- multi <- inc1 <- inc2 <- matrix(NA, nrow=nT, ncol=100)
  susc[,1] <- 1
  mono[,1] <- 0
  susc[1,2:100] <- susc0[1:99] - 4*lam[1]*susc0[1:99]
  mono[1,2:100] <- mono0[1:99] + 4*lam[1]*susc0[1:99] - 3*lam[1]*mono0[1:99]
  multi[1,] <- 1 - susc0 - mono0 
  inc1[1,] <- 4*lam[1]*susc0
  inc2[1,] <- 3*lam[1]*mono0
  
  # loop through time
  for(t in 2:nT){
    susc[t,2:100] <- susc[t-1,1:99] - 4*lam[t]*susc[t-1,1:99]
    mono[t,2:100] <- mono[t-1,1:99] + 4*lam[t]*susc[t-1,1:99] - 3*lam[t]*mono[t-1,1:99]
    multi[t,] <- 1 - susc[t,] - mono[t,]
    inc1[t,] <- 4*lam[t]*susc[t-1,]
    inc2[t,] <- 3*lam[t]*mono[t-1,]
  }
  
  # reported cases
  cases <- matrix(NA, nrow=nT, ncol=nA)
  for(t in 1:nT) for(a in 1:nA){
    cases[t,a] = rho*(mean(inc2[t,amin[a]:amax[a]]) + gamma*mean(inc1[t,amin[a]:amax[a]]))*pop[t,a]
  }
  
  return(list(cases=cases, susc=susc, mono=mono, multi=multi))
}

# choose parameter values for simulation
lam <- runif(10,0.001,0.05) # draw lambda t values
lamH <- 0.03
rho <- 0.12
gamma <- 0.1
plot(lam*4, type='l')
amin <- c(1,seq(10,90,10)) # lower bound of age groups
amax <- seq(9,100,10) # upper bound of age groups
pop <- matrix(100000, nrow=10, ncol=10) # population
sim <- simfoi(nT=10, nA=10, lamH=lamH, lam=lam, rho=rho, gamma=gamma,
                pop=pop, amin=amin, amax=amax)


#--- Fit to simulated data

# data for model fitting
cases <- round(sim$cases)
data <- list(nA=10, nT=10, cases=cases, pop=pop, age=seq(0,99), ageLims=rbind(amin,amax))

# fit model
check_cmdstan_toolchain(fix=T)
set_cmdstan_path('C:/Users/Megan/Documents/.cmdstanr/cmdstan-2.26.0')
setwd('C:/Users/Megan/Documents/GitHub/dengue-catalytic/StanCode')
mod <- cmdstan_model('timeVarying.stan', pedantic=T)
fit <- mod$sample(data=data, chains=3, parallel_chains=3, iter_sampling=5000, refresh=100, iter_warmup=1000)
stanfit <- rstan::read_stan_csv(fit$output_files())


# check convergence
chains <- rstan::extract(stanfit)
traceplot(stanfit, pars=c('rho','gamma','lam_H','lam_t'))
quantile(chains$rho, c(0.5,0.025,0.975))
quantile(chains$gamma, c(0.5,0.025,0.975))
quantile(chains$lam_H, c(0.5,0.025,0.975))

# extract lambda estimates
lam <- data.frame(t=seq(1,data$nT),true=lam, med=NA, ciL=NA, ciU=NA)
for(t in 1:data$nT) lam[t,3:5] <- quantile(chains$lam_t[,t],c(0.5,0.025,0.975))
ggplot(lam, aes(t,true))+ geom_point()+ theme_minimal()+ ylab('lambda')+
  geom_point(aes(t+0.2,med), col='orange')+ ylab('per serotype FOI')+ 
  geom_linerange(aes(t+0.2,ymin=ciL,ymax=ciU),col='orange')+ xlab('year')

# model fit
ageG <- c('0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80-89','90-99')
fit <- data.frame(cases)
colnames(fit) <- ageG
fit$year <- seq(1,10)
fit <- tidyr::gather(fit, key='age',value='cases',1:10)
fit[,c('pred','ciL','ciU')] <- NA
for(t in 1:data$nT) for(a in 1:data$nA){
  fit[fit$year==t & fit$age==ageG[a],4:6] <- quantile(chains$Ecases[,t,a], c(0.5,0.025,0.975))
}
ggplot(fit, aes(age, cases))+ geom_point()+ theme_minimal()+
  geom_line(aes(age,pred,group=year),col='dodgerblue')+ facet_wrap(~year)+
  geom_ribbon(aes(age,ymin=ciL, ymax=ciU,group=year),fill='dodgerblue',alpha=0.2)+
  theme(axis.text.x=element_text(angle=60,hjust=1), text=element_text(size=16))


