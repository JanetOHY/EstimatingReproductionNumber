########################################################################
## Estimate reproduction numbers with modified EpiEstim and EpiFilter ##
########################################################################

rm(list=ls())

source("~/function/epiFilter.R")
source("~/function/epiSmoother.R")

# load case data 
case.data = read.csv("~/Daily dengue cases (2010-2020).csv")

# Parameters for generation interval distribution 
min.mu_m = 0
max.mu_m = 0.2

min.theta_m = 1/15
max.theta_m = 1/2

min.mu_h = 1/83.6/365
max.mu_h = 1/81.7/365

min.theta_h = 1/10
max.theta_h = 1/3

min.alpha_h = 1/7
max.alpha_h = 1/2

min.c_m = 0
max.c_m = 1

# hyperparameter
alpha = 1
beta = 0.5
tau = 7

#####################################################
## Estimate reproduction numbers at national level ##
#####################################################

# Estimate Rt
Iday = case.data$count
dates = case.data$date

Rt_EpiEstim = Rt_EpiFilter = matrix(rep(NA,length(y)*5000),nrow=length(y))

for(i in 1:5000){
  print(i)
  g = rep(NA,length(y))
  
  # Generation interval distribution 
  theta_m = runif(1,min.theta_m,max.theta_m)
  mu_m = runif(1,min.mu_m,max.mu_m)
  c_m = runif(1,min.c_m,max.c_m)
  theta_h = runif(1,min.theta_h,max.theta_h)
  mu_h = runif(1,min.mu_h,max.mu_h)
  alpha_h = runif(1,min.alpha_h,max.alpha_h)
  
  s1 = theta_m+mu_m+c_m
  s2 = mu_m+c_m
  s3 = theta_h+mu_h
  s4 = alpha_h+mu_h
  
  for(t in 1:length(Iday)){
    g[t] = (s1*s2*s3*s4)*((exp(-s1*t)/((s2-s1)*(s3-s1)*(s4-s1)))+
                          (exp(-s2*t)/((s1-s2)*(s3-s2)*(s4-s2)))+
                          (exp(-s3*t)/((s1-s3)*(s2-s3)*(s4-s3)))+
                          (exp(-s4*t)/((s1-s4)*(s2-s4)*(s3-s4))))
  }
  wdist = g/sum(g)
  
  # EpiEstim method
  for(day in 32:length(Iday)){
    It = Iday[(day-tau+1):day]
    
    # Total infectiousness
    gs=c()
    for(s in (day-tau+1):day){
      Is = Iday[1:(s-1)]
      ws = wdist[1:(s-1)]
      gs[s] = ws%*%rev(Is)
    }
    Rt_EpiEstim[day,i] = rgamma(1,alpha+sum(It),beta+sum(gs,na.rm=T))
  }
  
  # EpiFilter method
  nday = length(dates);tday = 1:nday
  
  # Total infectiousness
  Lday = rep(0, nday) 
  for(i in 2:nday){
    Lday[i] = sum(Iday[seq(i-1, 1, -1)]*wdist[1:(i-1)])    
  }
  
  # Setup grid and noise parameters
  Rmin = 0.01; Rmax = 10; eta = 0.1
  
  # Uniform prior over grid of size m
  m = 200; pR0 = (1/m)*rep(1, m)
  
  # Delimited grid defining space of R
  Rgrid = seq(Rmin, Rmax, length.out = m)
  
  # Filtered estimates
  Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Lday[tday], Iday[tday], 0.025)
 
  # Smoothed estimates
  Rsmooth = epiSmoother(Rgrid, m, Rfilt[[4]], Rfilt[[5]], nday, Rfilt[[6]], 0.025)
  Rt_EpiFilter[2:nday,i]=Rsmooth[[3]][2:nday]
}


#########################################################
## Estimate reproduction numbers at fine spatial scale ##
#########################################################

rm(list=ls())

# load case data 
case.data = read.csv("~/Spatial dengue cases (2010-2020).csv")

# Estimate Rt
Iday = case.data$count
dates = unique(case.data$date)

Rt_EpiEstim=c()
for(i in 1:length(unique(case.data$SpatialUnitID))){
  Rt_EpiEstim[[i]] = matrix(rep(NA,5000*(length(dates))),nrow=(length(dates)))
}

Rt_EpiEstim = Rt_EpiFilter

for(i in 1:5000){
  print(i)
  # Generation interval distribution 
  theta_m = runif(1,min.theta_m,max.theta_m)
  mu_m = runif(1,min.mu_m,max.mu_m)
  c_m = runif(1,min.c_m,max.c_m)
  theta_h = runif(1,min.theta_h,max.theta_h)
  mu_h = runif(1,min.mu_h,max.mu_h)
  alpha_h = runif(1,min.alpha_h,max.alpha_h)
  
  s1 = theta_m+mu_m+c_m
  s2 = mu_m+c_m
  s3 = theta_h+mu_h
  s4 = alpha_h+mu_h
  
  g = rep(NA,length(dates))
  for(t in 1:length(dates)){
    g[t] = (s1*s2*s3*s4)*((exp(-s1*t)/((s2-s1)*(s3-s1)*(s4-s1)))+
                          (exp(-s2*t)/((s1-s2)*(s3-s2)*(s4-s2)))+
                          (exp(-s3*t)/((s1-s3)*(s2-s3)*(s4-s3)))+
                          (exp(-s4*t)/((s1-s4)*(s2-s4)*(s3-s4))))
  }
  wdist = g/sum(g)
  
  # EpiEstim method
  count = 1
  for(ID in sort(unique(case.data$SpatialUnitID))){
    Iday = case.data$count[case.data$SpatialUnitID==ID]
    
    for(day in 32:length(y)){
      It = y[(day-tau+1):day]
      gs=c()
      
      # Total infectiousness
      for(s in (day-tau+1):day){
        Is=Iday[1:(s-1)]
        ws=wdist[1:(s-1)]
        gs[s]=ws%*%rev(Is)
      }
      Rt_EpiEstim[[count]][day,i]=rgamma(1,alpha+sum(It),beta+sum(gs,na.rm=T))
    }
    
    # EpiFilter method
    nday = length(dates);tday = 1:nday
    
    # Total infectiousness
    Lday = rep(0, nday) 
    for(i in 2:nday){
      Lday[i] = sum(Iday[seq(i-1, 1, -1)]*wdist[1:(i-1)])    
    }
    
    # Setup grid and noise parameters
    Rmin = 0.01; Rmax = 10; eta = 0.1
    
    # Uniform prior over grid of size m
    m = 200; pR0 = (1/m)*rep(1, m)
    
    # Delimited grid defining space of R
    Rgrid = seq(Rmin, Rmax, length.out = m)
    
    # Filtered estimates
    Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Lday[tday], Iday[tday], 0.025)
    
    # Smoothed estimates
    Rsmooth = epiSmoother(Rgrid, m, Rfilt[[4]], Rfilt[[5]], nday, Rfilt[[6]], 0.025)
    Rt_EpiFilter[2:nday,i]=Rsmooth[[3]][2:nday]
    
    count=count+1
  }
}
