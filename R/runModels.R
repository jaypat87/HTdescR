#' Chemical prediction value generator
#'
#' R function that generates predictions for chemicals based on the machine learning method selected.
#'
#' @param data chemical descriptor data frame
#' @param method The model from which the prediction should be generated.
#'
#' @return a data frame containing indeces, rate constant prediction, applicability domain info, and descriptors
#'
#' @import stats
#' @import dplyr
#' @import FNN
#'
#' @export
#'
#' @examples ##NOT RUN:runModel (dataacidester, method = "SVR")
#'
runModel <- function(data, method){

  n <- nrow(data)
  PREDset <- data[c(2:3,5:7)]   #Prep data for use - only descriptors
  ADset <- dplyr::mutate(PREDset, index = seq_len(n)) #Prep data for use in AD calculation
  SetAppDom <- data.frame(Index=numeric(), AD=character(), ADLevel=numeric())
  predictionsDF <- data.frame() #Initialize dataframe for predictions -- do we need this?

  if (method == "SVR"){
    SVRpredictions <- stats::predict(ModelSVR, newdata = PREDset) #Predict rates
    PredictedRateSVR <- as.data.frame(SVRpredictions)
    SVRout <- dplyr::bind_cols(data, PredictedRateSVR) #Bind original data with predicted rates
    predictionsDF <- as.data.frame(SVRout[c(1,8,2:3,5:7)])


  } else if (method == "MLR"){
    MLRpredictions <- stats::predict(ModelMLR, newdata = PREDset)
    PredictedRateMLR <- as.data.frame(MLRpredictions)
    MLRout <- dplyr::bind_cols(data, PredictedRateMLR)
    predictionsDF <- as.data.frame(MLRout[c(1,8,2:3,5:7)])


  } else if (method == "RF"){
    RFpredictions <- stats::predict(ModelRF, newdata = PREDset)
    PredictedRateRF <- as.data.frame(RFpredictions)
    RFout <- dplyr::bind_cols(data, PredictedRateRF)
    predictionsDF <- as.data.frame(RFout[c(1,8,2:3,5:7)])


  } else if (method == "PLS"){
    PLSpredictions <- stats::predict(ModelPLS, newdata = PREDset)
    PLSpredictions <- as.data.frame(PLSpredictions)
    PredictedRatePLS <- PLSpredictions$`log.rate.exp.3 comps`
    PLSout <- dplyr::bind_cols(data, as.data.frame(PredictedRatePLS))
    predictionsDF <- as.data.frame(PLSout[c(1,8,2:3,5:7)])


  } else {
    stop("Error: Please specify a valid method.")
  }


  #AD calculation
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
  Indices <- SetAppDom[1]
  predicted <- predictionsDF[c(1:2)]
  desc <- predictionsDF[c(3:7)]
  Predictions <- dplyr::bind_cols(Indices, predicted, AppDomInfo, desc)
  return(Predictions)
}



