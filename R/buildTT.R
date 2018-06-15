
buildTT <- function(TrainingSet){
  #Use TrainingSet to get equations for all four methods
  TrainEqSVR <<- e1071::svm(Y~., data = TrainingSet, cost = 150, epsilon = 0.05, gamma = 0.00014)
  TrainEqPLS <<- pls::plsr(Y~., data = TrainingSet, ncomp = 3, scale=TRUE)
  TrainEqRF <<- randomForest::randomForest(Y~., data = TrainingSet, ntree = 500, nodesize = 5, mtry = 200, importance = TRUE, na.action = na.omit)
  TrainEqMLR <<- lm(Y~., data = TrainingSet)

  #Add equations to package data

}




