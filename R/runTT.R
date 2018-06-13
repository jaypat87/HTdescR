
runTT <- function(TrainingSet, TestSet) {
  MTrainingSet <<- TrainingSet[]
  MTestSet <<- TestSet[]

  #Predict TrainingSet Y for all methods
  TrainYSVR <<- predY(TrainEqSVR, TrainingSet)
  TrainYPLS <<- predY(TrainEqPLS, TrainingSet)
  TrainYRF <<- predY(TrainEqRF, TrainingSet)
  TrainYMLR <<- predY(TrainEqMLR, TrainingSet)

  #Get TrainingSet r2 for all methods
  Trainr2SVR <<- getr2(TrainYSVR, MTrainingSet)
  Trainr2PLS <<- getr2(TrainYPLS, MTrainingSet)
  Trainr2RF <<- getr2(TrainYRF, MTrainingSet)
  Trainr2MLR <<- getr2(TrainYMLR, MTrainingSet)


  #Predict TestSet Y for all methods
  TestYSVR <<- predY(TrainEqSVR, TestSet)
  TestYPLS <<- predY(TrainEqPLS, TestSet)
  TestYRF <<- predY(TrainEqRF, TestSet)
  TestYMLR <<- predY(TrainEqMLR, TestSet)

  #Get TestSet r2 for all methods
  Testr2SVR <<- getr2(TestYSVR, MTestSet)
  Testr2PLS <<- getr2(TestYPLS, MTestSet)
  Testr2RF <<- getr2(TestYRF, MTestSet)
  Testr2MLR <<- getr2(TestYMLR, MTestSet)




}



