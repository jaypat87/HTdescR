#Do we even need this equation? Will people ever need to make their own equations?

getEQ <- function(dataset, method){
  #Need to change defaults in SVR, PLS, and RF

  if (method =="SVR") {
      # support vector regression using e1071 library
      SVReq <- e1071::svm(Y~., data = dataset, cost = 150, epsilon = 0.05, gamma = 0.00014)
      return (SVMeq)
  } else if (method == "PLS") {
      # partial least squares using pls library
      PLSeq <- pls::plsr(Y~., data = dataset, ncomp = 3, scale=TRUE)
      return(PLSeq)
  } else if (method == "RF") {
      # random forest regression using randomforest library
      RFeq <- randomForest::randomForest(Y~., data = dataset, ntree = 500, nodesize = 5, mtry = 200, importance = TRUE, na.action = na.omit)
      return(RFeq)
  } else if (method == "MLR") {
      # multiple linear regression using base R linear model functions
      MLReq <- lm(Y~., data = dataset)
      return(MLSeq)
  } else {
      stop ("Specify valid method")
  }

}














