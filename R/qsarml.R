#' Predicting hydrolysis rates based on Quantitative Structure-Activity Relationships (QSAR) Models using Machine Learning Regressions (ML)
#'
#' Predict hydrolysis rates based on QSAR model using a variety of machine learning regressions such as random forest, support vector machine, partial least squares, and multiple linear regression.
#'
#' @param smile SMILES string for a chemical fragment in character, factor, or SMIset datatype.
#' @param HT.type The type of Hammett-Taft (HT) descriptor; valid inputs include "taft", "meta", "para", "ortho", "induction", "es" and "user" for a user defined sigma descriptor.
#' @param ... arguments from the fmcsR::fmcbatch such as a and b
#' @return A list containing tanimoto coefficient for the closest matching MCS, SMILES string of the MCS, and index nuber of the matched fragment from the library.
#'
#' @examples
#' htdesc (smile = "CCC", HT.type = "taft", sigma.selection = "A")

#create seperate functions for qsar based on random forests, pls, SVM, etc but put it in this file.

#building a qsar model

plsqsar <- function (TrainingSet, TestSet, ...) {
PLSqsar<-pls::plsr(LogRateExp~.,data=TrainingSet)
return(PLSqsar)
}

# calculating the r^2 of the model



predictqsar <- function (model, TestSet) {
PredictedTest<-predict(PLSqsar, newdata=TestSet)
PredictedTest
}

# get qsardataframe from htbatch
# clean it using cleanmerge
# subset it to keep only descriptors (X1, X2..) and property (Y) to be modeled
# split it into test and training set
# use training set for buildmlmodel

buildMLmodel <- function (qsardataframe, regression.method = "SVR", ...) {

  if (regression.method =="SVR") {
    # support vector regression using e1071 library
      SVMqsar <- e1071::svm(Y~., data= qsardataframe, cost = 150, epsilon = 0.05, gamma = 0.00014)
      return (SVMqsar)
    # call AD function
  } else if (regression.method == "pls") {
      # partial least squares using pls library
      PLSqsar<-pls::plsr(Y~.,data = qsardataframe, ncomp=2,scale=TRUE)
      return(PLSqsar)
    # call AD function
  } else if (regression.method == "RF") {
    # random forest regression using randomforest library
      RFqsar <- randomForest::randomForest(Y~., datab = qsardataframe, ntree=500, nodesize=5, mtry=200, importance=TRUE, na.action = na.omit)
      return(RFqsar)
    # call AD function
  } else if (regression.method == "MLR") {
      # multiple linear regression using base R linear model functions
      MLRqsar <- lm(Y~., data= qsardataframe)
    # call AD function
  } else {
    stop ("Specify valid regression.method")
  }
}

# PLSqsar <- buildMLmodel(qsardataframe= acidesterfilled, regression.method = "pls")
# using model above, use predict function to predict Y for training set
# PredPLSR <- predict(PLSqsar, newdata = qsardataframe)
# calculate the r^2
# CorrYtrainingPLSR<-lm(PredPLSR[,,2] ~ MeasuredYTraining)
# summary(CorrYtrainingPLSR)
# using model above, use predict function to predict Y for test set
# PredPLSR <- predict(PLSqsar, newdata = testset)
# calculate the r^2
# CorrYtrainingPLSR<-lm(PredPLSR[,,2] ~ MeasuredYTest)
# summary(CorrYtrainingPLSR)
