rm(list=ls())

library(glmnet)

###################################################################
## Determining temporal drivers of national reproduction numbers ##
###################################################################

# Load output from 01_EstimateReproductionNUmber first

var = read.csv("~/Temporal Covariates.csv")

var$date = as.Date(as.character(var$date),"%Y-%m-%d")

reg.est = matrix(rep(NA,(length(names(var))-1)*5000),nrow=(length(var)-1))

# LASSO model 
for(i in 1:5000){
  y = Rt_EpiEstim[32:nrow(Rt_EpiEstim),i]
  x = as.matrix(var[32:nrow(var),2:ncol(var)])
  lasso.fit = cv.glmnet(x=x, y=outcome, nfolds=10, alpha=1, family="gaussian", keep= T, standardize=TRUE)
  reg.est[,i] = coef(lasso.fit,s="lambda.min")[2:ncol(var),1]
}


###############################################################
## Determining spatial drivers of local reproduction numbers ##
###############################################################

rm(list=ls())

# Load output from 01_EstimateReproductionNUmber first

# Process output from 01_EstimateReproductionNUmber 
dates = unique(case.data$date)
ID = sort(unique(case.data$SpatialUnitID))

Rt_est=NULL
for(i in 1:length(ID)){
  Rt_mean = rowMeans(Rt_EpiFilter[[i]],na.rm=T)
  LB = apply(Rt_EpiFilter[[i]],1,quantile,prob=0.025,na.rm=T)
  UB = apply(Rt_EpiFilter[[i]],1,quantile,prob=0.975,na.rm=T)
  tmp = data.frame(date=dates,ID=ID[i],Rt=Rt_mean,LowerCI=LB,UpperCI=UB)
  Rt_est = rbind(Rt_est,tmp)
}

var = read.csv("~/Spatial Covariates.csv")

# Mean reproduction number as dependent variable
y = aggregate(Rt~SpatialUnitID,data=Rt_est,mean)

combined.data = merge(var,y,by="SpatialUnitID",all.y=T)
summary(step(glm(y ~ ., data=combined.data[,-c(1)])))

# Percentage of time reproduction number > 1.0 as dependent variable
Rt_est$ind = ifelse(Rt_est$LowerCI>1,1,0)
y = aggregate(ind~SpatialUnitID,data=Rt_est,sum)
y$proportion = y$ind/length(dates)

combined.data = merge(var,y[,c("SpatialUnitID","proportion")],by="SpatialUnitID",all.y=T)
summary(step(glm(proportion ~ ., data=combined.data[,-c(1)])))
