
runQSAR <- function(newdata, method){
  #Get equations from package data


  if (method =="SVR") {
    SVRpredY <- predY(TrainEqSVR, newdata)
    return (SVMpredY)
  } else if (method == "PLS") {
    PLSpredY <- predY(TrainEqPLS, newdata)
    return(PLSpredY)
  } else if (method == "RF") {
    RFpredY <- predY(TrainEqRF, newdata)
    return(RFpredY)
  } else if (method == "MLR") {
    MLRpredY <- predY(TrainEqMLR, newdata)
    return(MLSpredY)
  } else {
    stop ("Specify valid method")
  }

}


