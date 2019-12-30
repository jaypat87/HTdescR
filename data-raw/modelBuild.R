#Building Models with "Bootstrapping" Train/Test Split Scheme

# Libraries ---------------------------------------------------------------


library(caTools)
library(stats)
library(randomForest)
library(pls)
library(e1071)
library(Matrix)
library(xgboost)
library(gbm)
library(plyr)


# Data Preparation --------------------------------------------------------


#Get data
dataset <- read.csv("modeldata.csv")

#Remove low variance feature
modeldata <- dataset[2:6]

#Set seed and shuffle data
set.seed(40)
n <- nrow(modeldata)
modeldata <- modeldata[sample(n),]

#Divide dataset into features and targets
features <- modeldata[-1]
targets <- modeldata[,1]

#Generate random numbers for bootstrap seeds
set.seed(40)
randoms <- floor(runif(100, min=1, max=1000))

#Init results dataframes
MLRallstat <- data.frame()
PLSallstat <- data.frame()
RFallstat <- data.frame()
SVMallstat <- data.frame()

MLRpreds <- data.frame(targets)
PLSpreds <- data.frame(targets)
RFpreds <- data.frame(targets)
SVMpreds <- data.frame(targets)

MLRmods <- list(100)
PLSmods <- list(100)
RFmods <- list(100)
SVMmods <- list(100)


# Models ------------------------------------------------------------------



#"Bootstrap" loop including train/test split
for(i in 1:100){
  #Set randomly generated seed
  seed <- randoms[i]
  set.seed(seed)

  #Split train/test
  sample <- caTools::sample.split(targets, SplitRatio = 0.8) #Random splitting into Training and Test Sets
  Training <- subset(modeldata, sample == TRUE) #Training set with features and targets
  Testing <- subset(modeldata, sample == FALSE) #Testing set with features and targets
  TrainF <- subset(features, sample == TRUE) #Training set features
  TrainT <- subset(targets, sample == TRUE) #Training set targets
  TestF <- subset(features, sample == FALSE) #Testing set features
  TestT <- subset(targets, sample == FALSE) #Testing set targets



  #Run models

  #MLR model

  model <- stats::lm(log.rate.exp~., data=Training) #Train model
  MLRmods[[i]] <- model
  predIN <- stats::predict(model, newdata=TrainF, type="response") #Predict on training set to check internal model accuracy
  res <- predIN - TrainT #Training set residuals
  RMSEin <- sqrt(mean(res^2)) #Training set RMSE
  R2 <- 1 - (sum(res^2)) / (sum((TrainT - mean(TrainT))^2)) #Training set R2
  predOUT <- stats::predict(model, newdata=TestF, type="response") #Predict on testing set to check predictive accuracy
  res <- predOUT - TestT #Testing set residuals
  RMSEout <- sqrt(mean(res^2)) #Testing set RMSE
  Q2 <- 1 - (sum(res^2)) / (sum((TrainT - mean(TrainT))^2)) #Testing set R2 (called Q2 here)
  MLRstat <- data.frame(model="MLR",RMSEin, R2, RMSEout, Q2)
  MLRallstat <- rbind(MLRallstat, MLRstat)
  predALL <- stats::predict(model, newdata=features, type="response")
  MLRpreds <- data.frame(MLRpreds, predALL, check.names=TRUE)


  #PLS model

  model <- pls::plsr(log.rate.exp~., data = Training, ncomp = 4)
  PLSmods[[i]] <- model
  predIN <- stats::predict(model, newdata=TrainF, type="response") #Predict on training set to check internal model accuracy
  predIN <- data.frame(predIN)
  predIN <- predIN$log.rate.exp.4.comps
  res <- predIN - TrainT #Training set residuals
  RMSEin <- sqrt(mean(res^2)) #Training set RMSE
  R2 <- 1 - (sum(res^2)) / (sum((TrainT - mean(TrainT))^2)) #Training set R2
  predOUT <- stats::predict(model, newdata=TestF, type="response") #Predict on testing set to check predictive accuracy
  predOUT <- data.frame(predOUT)
  predOUT <- predOUT$log.rate.exp.4.comps
  res <- predOUT - TestT #Testing set residuals
  RMSEout <- sqrt(mean(res^2)) #Testing set RMSE
  Q2 <- 1 - (sum(res^2)) / (sum((TrainT - mean(TrainT))^2)) #Testing set R2 (called Q2 here)
  PLSstat <- data.frame(model="PLS",RMSEin, R2, RMSEout, Q2)
  PLSallstat <- rbind(PLSallstat, PLSstat)
  predALL <- stats::predict(model, newdata=features, type="response")
  predALL <- data.frame(predALL)
  predALL <- data.frame(predALL = predALL$log.rate.exp.4.comps)
  PLSpreds <- data.frame(PLSpreds, predALL, check.names=TRUE)


  #RF model

  model <- randomForest::randomForest(log.rate.exp~., data = Training, ntree = 500, nodesize = 5, mtry = 1, importance = TRUE, na.action = na.omit)
  RFmods[[i]] <- model
  predIN <- stats::predict(model, newdata=TrainF, type="response")
  res <- predIN - TrainT
  RMSEin <- sqrt(mean(res^2))
  R2 <- 1 - (sum(res^2)) / (sum((TrainT - mean(TrainT))^2))
  predOUT <- stats::predict(model, newdata=TestF, type="response")
  res <- predOUT - TestT
  RMSEout <- sqrt(mean(res^2))
  Q2 <- 1 - (sum(res^2)) / (sum((TrainT - mean(TrainT))^2))
  RFstat <- data.frame(model="RF",RMSEin, R2, RMSEout, Q2)
  RFallstat <- rbind(RFallstat, RFstat)
  predALL <- stats::predict(model, newdata=features, type="response")
  RFpreds <- data.frame(RFpreds, predALL, check.names=TRUE)


  #SVM model

  model <- e1071::svm(log.rate.exp~., data = Training, cost = 75.5, epsilon = 0.033, gamma = 0.9970)
  SVMmods[[i]] <- model
  predIN <- stats::predict(model, newdata=TrainF, type="response")
  res <- predIN - TrainT
  RMSEin <- sqrt(mean(res^2))
  R2 <- 1 - (sum(res^2)) / (sum((TrainT - mean(TrainT))^2))
  predOUT <- stats::predict(model, newdata=TestF, type="response")
  res <- predOUT - TestT
  RMSEout <- sqrt(mean(res^2))
  Q2 <- 1 - (sum(res^2)) / (sum((TrainT - mean(TrainT))^2))
  SVMstat <- data.frame(model="SVM",RMSEin, R2, RMSEout, Q2)
  SVMallstat <- rbind(SVMallstat, SVMstat)
  predALL <- stats::predict(model, newdata=features, type="response")
  SVMpreds <- data.frame(SVMpreds, predALL, check.names=TRUE)

}


#Mean responses

MLRpreds <- mutate(MLRpreds, meanpred = rowMeans(MLRpreds[,-1]))
PLSpreds <- mutate(PLSpreds, meanpred = rowMeans(PLSpreds[,-1]))
RFpreds <- mutate(RFpreds, meanpred = rowMeans(RFpreds[,-1]))
SVMpreds <- mutate(SVMpreds, meanpred = rowMeans(SVMpreds[,-1]))


#Rotating responses

mlrfit <- lm(MLRpreds$meanpred ~ MLRpreds$targets)
coMLRfit <- coef(mlrfit)
MLRpreds <- mutate(MLRpreds, rot = (MLRpreds$meanpred-coMLRfit[1])/coMLRfit[2])
mlrfitROT <- lm(MLRpreds$rot ~ MLRpreds$targets)
coMLRfitROT <- coef(mlrfitROT)

plsfit <- lm(PLSpreds$meanpred ~ PLSpreds$targets)
coPLSfit <- coef(plsfit)
PLSpreds <- mutate(PLSpreds, rot = (PLSpreds$meanpred-coPLSfit[1])/coPLSfit[2])
plsfitROT <- lm(PLSpreds$rot ~ PLSpreds$targets)
coPLSfitROT <- coef(plsfitROT)

rffit <- lm(RFpreds$meanpred ~ RFpreds$targets)
coRFfit <- coef(rffit)
RFpreds <- mutate(RFpreds, rot = (RFpreds$meanpred-coRFfit[1])/coRFfit[2])
rffitROT <- lm(RFpreds$rot ~ RFpreds$targets)
coRFfitROT <- coef(rffitROT)

svmfit <- lm(SVMpreds$meanpred ~ SVMpreds$targets)
coSVMfit <- coef(svmfit)
SVMpreds <- mutate(SVMpreds, rot = (SVMpreds$meanpred-coSVMfit[1])/coSVMfit[2])
svmfitROT <- lm(SVMpreds$rot ~ SVMpreds$targets)
coSVMfitROT <- coef(svmfitROT)


coeffics <- data.frame(MLR=coMLRfit, PLS=coPLSfit, RF=coRFfit, SVM=coSVMfit)




# Summary Statistics ------------------------------------------------------


#Summary stats

summary(MLRallstat)
summary(PLSallstat)
summary(RFallstat)
summary(SVMallstat)



# Calculate Applicability Domain ------------------------------------------

#Calculate AD using kNN
AD1 <- FNN::knn.dist(modeldata, k=5, algorithm=c("kd_tree", "cover_tree", "CR", "brute")) #Calculating the five nearest neighbors' distances for the Training Set
AD2 <- as.data.frame(AD1)
AD3 <- dplyr::mutate(AD2, Index = seq_len(70)) #Labeling rows with index numbers for sorting purposes
AD4 <- dplyr::mutate(AD3, avgDist = rowMeans(AD2)) #Finding the average KNN distance for each molecule (row)
AD4 <- dplyr::arrange(AD4, avgDist) #Sorting by the average distance

P95 <- stats::quantile(AD4$avgDist, prob = 0.95) #Finding the 95 percentile of the average KNN distance - this is the Applicability Domain Threshold

ADTrainSet <- tibble::rownames_to_column(modeldata)
ADTrainSet <- tibble::add_column(ADTrainSet, loopnum = 0)
ADTrain <- ADTrainSet[c(3:7)]


# Store Internal Data -----------------------------------------------------

list = c("es", "esSDF", "indSDF", "indsigma", "metaSDF", "metasigma", "orthoSDF", "orthosigma", "paraSDF", "parasigma", "taftSDF", "taftsigma", "userSDF", "usersigma")

for (a in list) {
  file = paste(c("./data/",a,".rda"), collapse = "")
  load(file)
}

htdescHelper <- read.csv("./inst/extdata/htdescHelper.csv")

usethis::use_data(MLRmods, PLSmods, RFmods, SVMmods, coeffics,
                  htdescHelper, P95, ADTrain,
                  esSDF, es, indSDF, indsigma, metaSDF, metasigma,
                  orthoSDF, orthosigma, paraSDF, parasigma, taftSDF, taftsigma,
                  userSDF, usersigma,
                  internal = TRUE, overwrite = TRUE)
