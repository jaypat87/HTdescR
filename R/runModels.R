
runModel <- function(data, method){

  n <- nrow(data)
  PREDset <- data[c(2:6)]   #Prep data for use - only descriptors
  ADset <- dplyr::mutate(PREDset, index = seq_len(n))
  SetAppDom <- data.frame(Index=numeric(), AD=character(), ADLevel=numeric())
  predictionsDF <- data.frame()

  if (method == "SVR"){
    SVRpredictions <- stats::predict(TrainEqSVR, newdata = PREDset)
    PredictedRateSVR <- as.data.frame(SVRpredictions)
    SVRout <- dplyr::bind_cols(data, PredictedRateSVR)
    predictionsDF <- as.data.frame(SVRout[c(1,7,2:6)])


  } else if (method == "MLR"){
    MLRpredictions <- stats::predict(TrainEqMLR, newdata = PREDset)
    PredictedRateMLR <- as.data.frame(MLRpredictions)
    MLRout <- bind_cols(data, PredictedRateMLR)
    predictionsDF <- as.data.frame(MLRout[c(1,7,2:6)])


  } else if (method == "RF"){
    RFpredictions <- stats::predict(TrainEqRF, newdata = PREDset)
    PredictedRateRF <- as.data.frame(RFpredictions)
    RFout <- bind_cols(data, PredictedRateRF)
    predictionsDF <- as.data.frame(RFout[c(1,7,2:6)])


  } else if (method == "PLS"){
    PLSpredictions <- stats::predict(TrainEqPLS, newdata = PREDset)
    PLSpredictions <- as.data.frame(PLSpredictions)
    PredictedRatePLS <- PLSpredictions$`log.rate.exp.3 comps`
    PLSout <- bind_cols(data, as.data.frame(PredictedRatePLS))
    predictionsDF <- as.data.frame(PLSout[c(1,7,2:6)])


  } else {
    stop("Error: Please specify a valid method.")
  }


  for(i in 1:n){

    TestChem <- dplyr::filter(ADset, index == i)
    CompareSet <- dplyr::bind_rows(ADTrain, TestChem)
    CompareSet <- CompareSet[c(1:5)]
    CompareAD1 <- FNN::knn.dist(CompareSet, k=5, algorithm=c("kd_tree", "cover_tree", "CR", "brute"))
    CompareAD1 <- as.data.frame(CompareAD1)
    CompareAD2 <- dplyr::mutate(CompareAD1, avgDist = rowMeans(CompareAD1))
    TestChemAD <- CompareAD2$avgDist[57]

    if(TestChemAD <= P95){
      SetAppDom <- dplyr::add_row(SetAppDom, Index = i, AD = "IN", ADLevel = TestChemAD)
    }else {
      SetAppDom <- dplyr::add_row(SetAppDom, Index = i, AD = "OUt", ADLevel = TestChemAD)
    }
  }

  AppDomInfo <- SetAppDom[c(2:3)]
  Indeces <- SetAppDom[1]
  predicted <- predictionsDF[c(1:2)]
  desc <- predictionsDF[c(3:7)]
  Predictions <- bind_cols(Indeces, predicted, AppDomInfo, desc)
  return(Predictions)
}



