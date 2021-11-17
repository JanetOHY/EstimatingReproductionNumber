rm(list=ls())

source("~/function/ModelAssessment.R")

case.data = read.csv("~/Daily dengue cases (2010-2020).csv")

Iday = case.data$count
dates = case.data$date

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

#######################################################
## Model assessment of national reproduction numbers ##
#######################################################

# Load output from 01_EstimateReproductionNUmber first 

## EpiEstim Estimates

## Estimated case counts
It=c()
for(day in 32:length(Iday)){
  Is=Iday[1:(day-1)]
  ws=wdist[1:(day-1)]
  gs=sum(ws*rev(Is))
  It[day]=Rt[day]*sum(gs)
}

ModelAssessment(Iday[32:length(dates)],It[32:length(dates)])

## EpiFilter Estimates

## Estimated case counts
Lday = rep(0, nday) 
for(i in 2:nday){
  Lday[i] = sum(Iday[seq(i-1, 1, -1)]*wdist[1:(i-1)])    
}

pred = Lday*Rmean

ModelAssessment(Iday[32:length(dates)],pred[32:length(dates)])

####################################################
## Model assessment of local reproduction numbers ##
####################################################

rm(list=ls())

library(spdep)
library(maptools)

spatialunit = readShapeSpatial("~/Spatial unit layer.shp")

case.data = read.csv("~/Spatial dengue cases (2010-2020).csv")

dates = case.data$date

# Load output from 01_EstimateReproductionNUmber first

## EpiEstim Estimates

## Process output from 01_EstimateReproductionNUmber 
dates = unique(case.data$date)
ID = sort(unique(case.data$SpatialUnitID))

Rt_est=NULL
for(i in 1:length(ID)){
  Rt_mean = rowMeans(Rt_EpiEstim[[i]],na.rm=T)
  LB = apply(Rt_EpiEstim[[i]],1,quantile,prob=0.025,na.rm=T)
  UB = apply(Rt_EpiEstim[[i]],1,quantile,prob=0.975,na.rm=T)
  tmp = data.frame(date=dates,ID=ID[i],Rt=Rt_mean,LowerCI=LB,UpperCI=UB)
  Rt_est = rbind(Rt_est,tmp)
}

EvaluateEpiEstim=NULL
for(i in 1:length(ID)){
  Iday = case.data$count[case.data$SpatialUnitID==ID[i]]
  Rt = Rt_est$Rt[Rt_est$ID==ID[i]]
  
  ## Estimated case counts
  It=c()
  for(day in 32:length(Iday)){
    Is = Iday[1:(day-1)]
    ws = wdist[1:(day-1)]
    gs = sum(ws*rev(Is))
    It[day] = Rt[day]*sum(gs)
  }
  
  metrics = ModelAssessment(Iday[32:length(dates)],pred[32:length(dates)])
  tmp = data.frame(ID=ID,Rsq=metrics$Rsq,MSE=metrics$MSE,MASE=metrics$MASE)
  EvaluateEpiEstim = rbind(EvaluateEpiEstim,tmp)
}

## Moran's I test
Rt_est$ind = ifelse(Rt_est$LowerCI>1,1,0)
y1 = aggregate(ind~spatialunit,data=Rt_est,sum)
y2 = aggregate(Rt~SpatialUnitID,data=Rt_est,mean)
y1$proportion = y1$ind/length(dates)

combined.data = merge(spatialunit,y1,by="SpatialUnitID",all.x=T)
combined.data = merge(combined.data,y2,by="SpatialUnitID",all.x=T)
neighbour = poly2nb(combined.data, queen=TRUE)

# Mean reproduction number 
moran.test(combined.data$Rt, nb2listw(neighbour))

## Percentage of time reproduction number > 1.0
moran.test(combined.data$proportion, nb2listw(neighbour))


## EpiFilter Estimates

## Process output from 01_EstimateReproductionNUmber 
Rt_est=NULL
for(i in 1:length(ID)){
  Rt_mean = rowMeans(Rt_EpiFilter[[i]],na.rm=T)
  LB = apply(Rt_EpiFilter[[i]],1,quantile,prob=0.025,na.rm=T)
  UB = apply(Rt_EpiFilter[[i]],1,quantile,prob=0.975,na.rm=T)
  tmp = data.frame(date=dates,ID=ID[i],Rt=Rt_mean,LowerCI=LB,UpperCI=UB)
  Rt_est = rbind(Rt_est,tmp)
}

EvaluateEpiFilter=NULL
for(i in 1:length(ID)){
  Iday = case.data$count[case.data$SpatialUnitID==ID[i]]
  Rt = Rt_est$Rt[Rt_est$ID==ID[i]]
  
  ## Estimated case counts
  Lday = rep(0, length(dates)) 
  for(i in 2:nday){
    Lday[i] = sum(Iday[seq(i-1, 1, -1)]*wdist[1:(i-1)])    
  }
  
  pred = Lday*Rt
  
  metrics = ModelAssessment(Iday[32:length(dates)],pred[32:length(dates)])
  tmp = data.frame(ID=ID,Rsq=metrics$Rsq,MSE=metrics$MSE,MASE=metrics$MASE)
  EvaluateEpiFilter = rbind(EvaluateEpiFilter,tmp)
}

## Moran's I test
Rt_est$ind = ifelse(Rt_est$LowerCI>1,1,0)
y1 = aggregate(ind~spatialunit,data=Rt_est,sum)
y2 = aggregate(Rt~SpatialUnitID,data=Rt_est,mean)
y1$proportion = y1$ind/length(dates)

combined.data = merge(spatialunit,y1,by="SpatialUnitID",all.x=T)
combined.data = merge(combined.data,y2,by="SpatialUnitID",all.x=T)
neighbour = poly2nb(combined.data, queen=TRUE)

# Mean reproduction number 
moran.test(combined.data$Rt, nb2listw(neighbour))

## Percentage of time reproduction number > 1.0
moran.test(combined.data$proportion, nb2listw(neighbour))