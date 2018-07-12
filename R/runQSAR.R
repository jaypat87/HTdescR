
runQSAR <- function(data, method){

  #Get models from package data

  allPredictions <- extractPrediction(extractedModels, unkX = data, unkOnly = TRUE, verbose = FALSE)
  SVMpredictions <- predict(SVMmodels, newdata = data, type = "raw", na.action = "na.omit")
  SVMpredictions <- data.frame(SVMpredictions)
  SVMpredictions <- tibble::add_column(SVMpredictions, index = 1:nrow(SVMpredictions))

  if (method =="SVMR") {
    SVMRpredictions <- SVMpredictions[c("SVMRmodel", "index")]
    return (SVMRpredictions)

  } else if (method == "SVMP") {
    SVMPpredictions <- SVMpredictions[c("SVMPmodel", "index")]
    return(SVMPpredictions)

  } else if (method == "PLS") {
    PLSpredictions <- filter(allPredictions, model == "pls")
    PLSpredictions <- tibble::add_column(PLSpredictions, index = 1:nrow(PLSpredictions))
    return(PLSpredictions)

  } else if (method == "RF") {
    RFpredictions <- filter(allPredictions, model == "rf")
    RFpredictions <- tibble::add_column(RFpredictions, index = 1:nrow(RFpredictions))
    return(RFpredictions)

  } else if (method == "MLR") {
    MLRpredictions <- filter(allPredictions, model == "lm")
    MLRpredictions <- tibble::add_column(MLRpredictions, index = 1:nrow(MLRpredictions))
    return(MLRpredictions)

  } else {
    stop ("Specify valid method")
  }

}


