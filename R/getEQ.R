
getEQ <- function(dataset, method){
  if (method =="SVR") {
      # support vector regression using e1071 library
      SVRqsar <- e1071::svm(Y~., data = dataset, cost = 150, epsilon = 0.05, gamma = 0.00014)
      return (SVMqsar)
  } else if (method == "PLS") {
      # partial least squares using pls library
      PLSqsar<-pls::plsr(Y~., data = dataset, ncomp = 3, scale=TRUE)
      return(PLSqsar)
  } else if (method == "RF") {
      # random forest regression using randomforest library
      RFqsar <- randomForest::randomForest(Y~., data = dataset, ntree = 500, nodesize = 5, mtry = 200, importance = TRUE, na.action = na.omit)
      return(RFqsar)
  } else if (method == "MLR") {
      # multiple linear regression using base R linear model functions
      MLRqsar <- lm(Y~., data = dataset)
      return(MLSqsar)
  } else {
      stop ("Specify valid method")
  }

}














