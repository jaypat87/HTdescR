
buildTT <- function(TrainingSet, TestSet){
  #Use TrainingSet to get equations for all four methods
  TrainEqSVR <<- getEQ(TrainingSet, method = "SVR")
  TrainEqPLS <<- getEQ(TrainingSet, method = "PLS")
  TrainEqRF <<- getEQ(TrainingSet, method = "RF")
  TrainEqMLR <<- getEQ(TrainingSet, method = "MLR")

}




