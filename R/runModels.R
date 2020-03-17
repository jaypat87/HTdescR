#' Chemical prediction value generator
#'
#' R function that generates predictions for chemicals based on the machine learning method selected.
#'
#' @param data chemical descriptor data frame
#' @param method The model from which the prediction should be generated. Valid inputs are "SVR", "MLR", "RF", and "PLS".
#'
#' @return a data frame containing indices, rate constant prediction, applicability domain info, and descriptors
#'
#' @importFrom stats predict
#' @importFrom dplyr mutate filter bind_rows bind_cols add_row
#' @importFrom FNN knn.dist
#'
#' @export
#'
#' @examples ##NOT RUN:runModel (dataacidester, method = "SVR")
#'
runModel <- function(data, method){

  n <- nrow(data)
  preds <- data.frame(row.names = 1:n) #Initialize dataframe for PLS predictions
  PREDset <- data[c(2:3,5:6)]   #Prep data for use - only descriptors
  ADset <- dplyr::mutate(PREDset, index = seq_len(n)) #Prep data for use in AD calculation - descriptors and index numbers
  SetAppDom <- data.frame(Index = numeric(), AD = character(), ADLevel = numeric()) #Initialize applicability domain dataframe
  predictionsDF <- data.frame() #Initialize dataframe for predictions -- do we need this?

  if (method == "SVR") {
    i <- 1
    for (m in SVMmods) {
      preds[i] <- stats::predict(m, newdata = PREDset, type = "response") #Predict rates over all 100 saved models
      i = i + 1
    }

    #preds <- stats::predict(SVMmods, newdata = PREDset, type = "response") #Predict rates
    preds <- as.data.frame(preds)
    l <- ncol(preds)
    names(preds) <- 1:l

    #Average predictions
    avgpreds <- as.data.frame(rowMeans(preds))

    #Rotate predictions
    rotpreds <- (avgpreds - coeffics[[4]][1])/coeffics[[4]][2]
    logpreds <- rotpreds
    colnames(logpreds) <- "logpreds"

    #Antilog predictions
    ratepreds <- 10^rotpreds
    colnames(ratepreds) <- "ratepreds"

    SVRout <- data.frame(data, logpreds, ratepreds) #Bind original data with predicted rates
    predictionsDF <- as.data.frame(SVRout[c(1,8:9,2:3,5:6)]) #Order columns for predictions dataframe


  } else if (method == "MLR") {
    i <- 1
    for (m in MLRmods) {
      preds[i] <- stats::predict(m, newdata = PREDset, type = "response") #Predict rates over all 100 saved models
      i = i + 1
    }

    #preds <- stats::predict(MLRmods, newdata = PREDset, type = "response") #Predict rates
    preds <- as.data.frame(preds)
    l <- ncol(preds)
    names(preds) <- 1:l

    #Average predictions
    avgpreds <- as.data.frame(rowMeans(preds))

    #Rotate predictions
    rotpreds <- (avgpreds - coeffics[[1]][1])/coeffics[[1]][2]
    logpreds <- rotpreds
    colnames(logpreds) <- "logpreds"

    #Antilog predictions
    ratepreds <- 10^rotpreds
    colnames(ratepreds) <- "ratepreds"

    MLRout <- data.frame(data, logpreds, ratepreds) #Bind original data with predicted rates
    predictionsDF <- as.data.frame(MLRout[c(1,8:9,2:3,5:6)]) #Order columns for predictions dataframe


  } else if (method == "RF") {
    i <- 1
    for (m in RFmods) {
      preds[i] <- stats::predict(m, newdata = PREDset, type = "response") #Predict rates over all 100 saved models
      i = i + 1
    }

    #preds <- stats::predict(RFmods, newdata=PREDset, type="response")
    preds <- as.data.frame(preds)
    l <- ncol(preds)
    names(preds) <- 1:l

    #Average predictions
    avgpreds <- as.data.frame(rowMeans(preds))

    #Rotate predictions
    rotpreds <- (avgpreds - coeffics[[3]][1])/coeffics[[3]][2]
    logpreds <- rotpreds
    colnames(logpreds) <- "logpreds"

    #Antilog predictions
    ratepreds <- 10^rotpreds
    colnames(ratepreds) <- "ratepreds"

    RFout <- data.frame(data, logpreds, ratepreds) #Bind original data with predicted rates
    predictionsDF <- as.data.frame(RFout[c(1,8:9,2:3,5:6)]) #Order columns for predictions dataframe


  } else if (method == "PLS") {
    i <- 1
    for (m in PLSmods) {
      pred_pls <- stats::predict(m, newdata = PREDset, type = "response") #Predict rates over all 100 saved models, output includes predictions for 1, 2, 3, and 4 components
      pred_pls <- data.frame(pred_pls) #Convert to dataframe to pull out the "4 components" predictions
      preds[i] <- data.frame(preds = pred_pls$log.rate.exp.4.comps) #Pull out "4 components" prediction column
      i = i + 1
    }

    #preds <- stats::predict(PLSmods, newdata = PREDset, type = "response") #Predict rates
    preds <- as.data.frame(preds)
    l <- ncol(preds)
    names(preds) <- 1:l

    #Average predictions
    avgpreds <- as.data.frame(rowMeans(preds))

    #Rotate predictions
    rotpreds <- (avgpreds - coeffics[[2]][1])/coeffics[[2]][2]
    logpreds <- rotpreds
    colnames(logpreds) <- "logpreds"

    #Antilog predictions
    ratepreds <- 10^rotpreds
    colnames(ratepreds) <- "ratepreds"

    PLSout <- data.frame(data, logpreds, ratepreds) #Bind original data with predicted rates
    predictionsDF <- as.data.frame(PLSout[c(1,8:9,2:3,5:6)]) #Order columns for predictions dataframe


  } else {
    stop("Error: Please specify a valid method.")
  }



  #AD calculation - for each chemical in new data, calculate AD and add to dataframe
  for (i in 1:n) {

    TestChem <- dplyr::filter(ADset, index == i) #Filter the new data to pick out each chemical individually
    CompareSet <- dplyr::bind_rows(ADTrain, TestChem) #Add individual chemical to model data
    CompareSet1 <- CompareSet[c(1:4)] #Pull out just descriptors, no loopnum or index columns
    CompareAD1 <- FNN::knn.dist(CompareSet1, k = 5, algorithm = c("kd_tree", "cover_tree", "CR", "brute")) #Find distances of 5 nearest neighbors
    CompareAD1 <- as.data.frame(CompareAD1)
    CompareAD2 <- dplyr::mutate(CompareAD1, avgDist = rowMeans(CompareAD1)) #Calculate average distance
    TestChemAD <- CompareAD2$avgDist[71] #Store average distance for new chemical

    if (TestChemAD <= P95) {
      SetAppDom <- dplyr::add_row(SetAppDom, Index = i, AD = "IN", ADLevel = TestChemAD) #If avg distance for new chemical is less than 95% threshold, chemical is in AD
    }else {
      SetAppDom <- dplyr::add_row(SetAppDom, Index = i, AD = "OUT", ADLevel = TestChemAD) #If not, chemical is out of AD
    }
  }

  AppDomInfo <- SetAppDom[c(2:3)] #Pull out AD in/out and AD level
  Indices <- SetAppDom[1] #Pull out indices
  predicted <- predictionsDF[c(1:3)] #Pull out experimental values, logratepreds, ratepreds
  desc <- predictionsDF[c(4:7)] #Pull out descriptors
  Predictions <- dplyr::bind_cols(Indices, predicted, AppDomInfo, desc) #Bind all together to form final dataframe
  return(Predictions)
}



