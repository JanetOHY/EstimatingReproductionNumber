ModelAssessment = function(obs,pred){
  rsq = (cor.test(pred,obs)$estimate)**2
  mse = sum((pred-obs)**2)/length(obs])
  mase = mean(abs(obs-pred))/sum(abs(diff(obs)))/(length(obs)-1)
  return(list(Rsq=rsq, MSE=mse, MASE=MASE))
}

