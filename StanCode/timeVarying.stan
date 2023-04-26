//--- Time-varying dengue catalytic model ---//
// assumes constant endemic FOI prior to data
// assumes complete immunity after 2nd infection
// assumes equal transmissability of 4 serotypes

data {
  
  int nA; // N age groups
  int nT; // N time points
  int cases[nT,nA]; // reported case data
  matrix[nT,nA] pop; // population data
  int ageLims[2,nA]; // lower & upper bounds of age groups
  row_vector[100] age;
  
}


parameters {
  
  real<lower=0,upper=0.15> lam_H; // historic average FOI
  real<lower=0,upper=0.15> lam_t[nT]; // time varying FOI
  real<lower=0,upper=1> rho; // reporting rate of 2nd infections
  real<lower=0,upper=1> gamma; // relative reporting rate of 1st infections

}

transformed parameters {
  
  row_vector<lower=0,upper=1>[100] susc0;
  row_vector<lower=0,upper=1>[100] mono0;
  row_vector<lower=0,upper=1>[100] multi0;
  matrix<lower=0,upper=1>[nT,100] susc; // proportion susceptible
  matrix<lower=0,upper=1>[nT,100] mono; // proportion monotypic
  matrix<lower=0,upper=1>[nT,100] multi; // proportion multitypic
  matrix<lower=0,upper=1>[nT,100] inc1; // incidence of primary infections
  matrix<lower=0,upper=1>[nT,100] inc2; // incidence of secondary infections
  vector<lower=0>[nA] Ecases[nT]; // expected reported cases
  
  
  //--- immune profiles at beginning of time series
  susc0 = exp(-4*lam_H*age);
  mono0 = 4*exp(-3*lam_H*age).*(1-exp(-lam_H*age));
  multi0 = 1 - (susc0 + mono0);
  
  //--- infants (assumes no infections in <1 year olds)
  for(t in 1:nT) susc[t,1] = 1;
  for(t in 1:nT) mono[t,1] = 0;
  for(t in 1:nT) multi[t,1] = 0;
  
  
  //--- subsequent time steps
  susc[1,2:100] = susc0[1:99] - 4*lam_t[1]*susc0[1:99];
  mono[1,2:100] = mono0[1:99] + 4*lam_t[1]*susc0[1:99] - 3*lam_t[1]*mono0[1:99];
  multi[1,2:100] = multi0[1:99] + 3*lam_t[1]*mono0[1:99];
  inc1[1,] = 4*lam_t[1]*susc0;
  inc2[1,] = 3*lam_t[1]*mono0;
  
  for(t in 2:nT){
    
    susc[t,2:100] = susc[t-1,1:99] - 4*lam_t[t]*susc[t-1,1:99];
    mono[t,2:100] = mono[t-1,1:99] + 4*lam_t[t]*susc[t-1,1:99] - 3*lam_t[t]*mono[t-1,1:99];
    multi[t,2:100] = multi[t-1,1:99] + 3*lam_t[t]*mono[t-1,1:99]; 
    inc1[t,] = 4*lam_t[t]*susc[t-1,];
    inc2[t,] = 3*lam_t[t]*mono[t-1,];
    
  }
  
  //--- reported cases
  for(t in 1:nT) for(a in 1:nA){
    Ecases[t,a] = rho*(mean(inc2[t,ageLims[1,a]:ageLims[2,a]]) + gamma*mean(inc1[t,ageLims[1,a]:ageLims[2,a]]))*pop[t,a];
  }

  
}


model {
  
  //--- priors
  lam_H ~ normal(0,0.05);
  lam_t ~ normal(0,0.05);
  rho ~ normal(0,0.25);
  gamma ~ normal(0,0.25);
  
  
  
  //--- likelihood 
  for(t in 1:nT) cases[t,] ~ poisson(Ecases[t,]); // poisson likelihood
  //for(t in 1:nT) cases[t,] ~ neg_binomial_2(Ecases[t,], phi); // negative binomial likelihood with estimated phi
  
}



