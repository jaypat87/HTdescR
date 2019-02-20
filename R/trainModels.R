#Random Split for Training and Test Sets

dataacidester <- read.csv("./data/dataacidester.csv") #Import data file
dataacidester <- dataacidester[c(2:8)] #Trim index numbers from data frame
moddata <- dataacidester
colnames(moddata)[1] <- "rate" #Renaming rate column
moddata <- dplyr::mutate(moddata, log.rate.exp = log10(rate)) #Making a column of log rates
moddata <- moddata[c(8,2:3,5:7)] #Sorting columns of rate and descriptors to make sense visually
moddata <- dplyr::slice(moddata, c(1:22,28:70)) #Cutting out outliers for testing


set.seed(529) #Setting seed for Train/Test splitting
sample <- caTools::sample.split(moddata$log.rate.exp, SplitRatio = 0.8) #Random splitting into Training and Test Sets
moddataTrain <- subset(moddata, sample == TRUE) #Saving the Training Set
moddataTest <- subset(moddata, sample == FALSE) #Saving the Test Set

TrainingSet <- moddataTrain #Saving sets with new names
TestSet <- moddataTest
TrainingDesc <- TrainingSet[c(2:6)] #Making data frames with only the descriptors for the Training and Test Sets (no rates)
TestDesc <- TestSet[c(2:6)]
MTrain <- moddataTrain$log.rate.exp #Saving data frames with only the experimental rates for the Training and Test Sets
MTest <- moddataTest$log.rate.exp



#SVMradial model

SVMRmod <- caret::train(log.rate.exp ~ ., data = TrainingSet, method = "svmRadial", metric = "RMSE",
                        trControl = trainControl(method = "cv", number = 5, summaryFunction = defaultSummary,
                                                 returnData = TRUE, returnResamp = "all", savePredictions = "final"))
SVMRresults <- SVMRmod$results
SVMRresample <- SVMRmod$resample
SVMRbestTune <- SVMRmod$bestTune
SVMRfinalModel <- SVMRmod$finalModel
SVMRpred <- SVMRmod$pred
SVMRpred2 <- arrange(SVMRpred, rowIndex)
SVMRtrain <- SVMRmod$trainingData
SVMRobs <- SVMRtrain[c(-1)]
SVMRpred3 <- select(SVMRpred2, pred)
SVMRobs <- bind_cols(SVMRobs, SVMRpred3)
SVMRobs <- SVMRobs[c(6,1:5)]
SVMRtrain <- rename(SVMRtrain, meas = .outcome)
SVMRdata <- full_join(SVMRtrain, SVMRobs)
SVMRdata <- SVMRdata[c(2:6,1,7)]

SVMRrmseTrain <- mutate(SVMRpred, residuals = obs - pred)
SVMRrmseTrain <- mutate(SVMRrmseTrain, residuals2 = residuals^2)
SVMRrmseTrain <- mutate(SVMRrmseTrain, mean = mean(residuals2))
SVMRrmseTrain <- mutate(SVMRrmseTrain, RMSE = sqrt(mean))


#PLS model

PLSmod <- caret::train(log.rate.exp ~ ., data = TrainingSet, method = "pls", metric = "RMSE",
                       trControl = trainControl(method = "cv", number = 5, summaryFunction = defaultSummary,
                                                returnData = TRUE, returnResamp = "all", savePredictions = "final"))
PLSresults <- PLSmod$results
PLSresample <- PLSmod$resample
PLSbestTune <- PLSmod$bestTune
PLSfinalModel <- PLSmod$finalModel
PLSpred <- PLSmod$pred
PLSpred2 <- arrange(PLSpred, rowIndex)
PLStrain <- PLSmod$trainingData
PLSobs <- PLStrain[c(-1)]
PLSpred3 <- select(PLSpred2, pred)
PLSobs <- bind_cols(PLSobs, PLSpred3)
PLSobs <- PLSobs[c(6,1:5)]
PLStrain <- rename(PLStrain, meas = .outcome)
PLSdata <- full_join(PLStrain, PLSobs)
PLSdata <- PLSdata[c(2:6,1,7)]




#MLR model

MLRmod <- caret::train(log.rate.exp ~ ., data = TrainingSet, method = "lm", metric = "RMSE",
                       trControl = trainControl(method = "cv", number = 5, summaryFunction = defaultSummary,
                                                returnData = TRUE, returnResamp = "all", savePredictions = "final"))
MLRresults <- MLRmod$results
MLRresample <- MLRmod$resample
MLRbestTune <- MLRmod$bestTune
MLRfinalModel <- MLRmod$finalModel
MLRpred <- MLRmod$pred
MLRpred2 <- arrange(MLRpred, rowIndex)
MLRtrain <- MLRmod$trainingData
MLRobs <- MLRtrain[c(-1)]
MLRpred3 <- select(MLRpred2, pred)
MLRobs <- bind_cols(MLRobs, MLRpred3)
MLRobs <- MLRobs[c(6,1:5)]
MLRtrain <- rename(MLRtrain, meas = .outcome)
MLRdata <- full_join(MLRtrain, MLRobs)
MLRdata <- MLRdata[c(2:6,1,7)]

MLRrmseTrain <- mutate(MLRpred, residuals = obs - pred)
MLRrmseTrain <- mutate(MLRrmseTrain, residuals2 = residuals^2)
MLRrmseTrain <- mutate(MLRrmseTrain, mean = mean(residuals2))
MLRrmseTrain <- mutate(MLRrmseTrain, RMSE = sqrt(mean))



#RF model

RFmod <- caret::train(log.rate.exp ~ ., data = TrainingSet, method = "rf", metric = "RMSE",
                      trControl = trainControl(method = "cv", number = 5, summaryFunction = defaultSummary,
                                               returnData = TRUE, returnResamp = "all", savePredictions = "final"))
RFresults <- RFmod$results
RFresample <- RFmod$resample
RFbestTune <- RFmod$bestTune
RFfinalModel <- RFmod$finalModel
RFpred <- RFmod$pred
RFpred2 <- arrange(RFpred, rowIndex)
RFtrain <- RFmod$trainingData
RFobs <- RFtrain[c(-1)]
RFpred3 <- select(RFpred2, pred)
RFobs <- bind_cols(RFobs, RFpred3)
RFobs <- RFobs[c(6,1:5)]
RFtrain <- rename(RFtrain, meas = .outcome)
RFdata <- full_join(RFtrain, RFobs)
RFdata <- RFdata[c(2:6,1,7)]
