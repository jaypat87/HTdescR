#MLR comparisons


#Caret CV TrainingSet MLR model

set.seed(300)
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
CVtrainMLRrmse <- MLRrmseTrain$RMSE[[1]]


#Caret CV moddata MLR model

set.seed(300)
MLRmod <- caret::train(log.rate.exp ~ ., data = moddata, method = "lm", metric = "RMSE",
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
CVallMLRrmse <- MLRrmseTrain$RMSE[[1]]








#Caret LOOCV TrainingSet MLR model

set.seed(300)
MLRmod <- caret::train(log.rate.exp ~ ., data = TrainingSet, method = "lm", metric = "RMSE",
                       trControl = trainControl(method = "LOOCV", summaryFunction = defaultSummary,
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
LOOCVtrainMLRrmse <- MLRrmseTrain$RMSE[[1]]


#Caret LOOCV moddata MLR model

set.seed(300)
MLRmod <- caret::train(log.rate.exp ~ ., data = moddata, method = "lm", metric = "RMSE",
                       trControl = trainControl(method = "LOOCV", summaryFunction = defaultSummary,
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
LOOCVallMLRrmse <- MLRrmseTrain$RMSE[[1]]





#Caret reg train MLR model

set.seed(300)
MLRmod <- caret::train(log.rate.exp ~ ., data = TrainingSet, method = "lm", metric = "RMSE",
                       trControl = trainControl(method = "none", summaryFunction = defaultSummary,
                       returnData = TRUE, savePredictions = "final"))

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
CaretRegTrainMLRrmse <- MLRrmseTrain$RMSE[[1]]




#Caret reg moddata MLR model

set.seed(300)
MLRmod <- caret::train(log.rate.exp ~ ., data = moddata, method = "lm", metric = "RMSE",
                       trControl = trainControl(summaryFunction = defaultSummary,
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
CaretRegAllMLRrmse <- MLRrmseTrain$RMSE[[1]]


#Regular MLR model

set.seed(300)

TrainEqMLR <- stats::lm(log.rate.exp~., data = TrainingSet) #Training the MLR model on the Training Set
TrainYMLR <- stats::predict(TrainEqMLR, newdata = TrainingDesc) #Using the model to predict the Training Set
Trainr2MLR <- stats::lm(TrainYMLR ~ MTrain) #Making a simple regression of the experimental vs predicted rates for the Training Set
sumTrainMLR <- summary(Trainr2MLR)
multR2MLR <- as.numeric(sumTrainMLR[8]) #Grabbing mult. R2 for this model
adjR2MLR <- as.numeric(sumTrainMLR[9]) #Grabbing adj. R2 for this model
ResStErrMLR <- as.numeric(sumTrainMLR[6]) #Grabbing residual standard error for this model

TestYMLR <- stats::predict(TrainEqMLR, newdata = TestDesc) #Using the model to predict the Test Set
TestMLRreg <- stats::lm(TestYMLR ~ MTest)
sumTestMLR <- summary(TestMLRreg)


MLRresiduals <- as.data.frame(sumTrainMLR$residuals)
MLRresiduals2 <- data.frame(MLRresiduals^2)
meanMLRresiduals2 <- mutate(MLRresiduals2, mean(MLRresiduals2$sumTrainMLR.residuals))
MLRtrainRMSE <- sqrt(meanMLRresiduals2$`mean(MLRresiduals2$sumTrainMLR.residuals)`[1])

MLRresidualsTest <- as.data.frame(sumTestMLR$residuals)
MLRresidualsTest2 <- data.frame(MLRresidualsTest^2)
meanMLRresidualsTest2 <- mutate(MLRresidualsTest2, mean(MLRresidualsTest2$sumTestMLR.residuals))
MLRtestRMSE <- sqrt(meanMLRresidualsTest2$`mean(MLRresidualsTest2$sumTestMLR.residuals`[1])




#Aggregate

MLRcompareCVrmse <- data.frame(CVtrainMLRrmse, LOOCVtrainMLRrmse, CVallMLRrmse, LOOCVallMLRrmse)
MLRcompareRMSE <- data.frame(CVtrainMLRrmse, CVallMLRrmse, MLRtrainRMSE, MLRtestRMSE)
MLRcompareRMSE2 <- data.frame(CVtrainMLRrmse, CVallMLRrmse, CaretRegTrainMLRrmse, CaretRegAllMLRrmse,
                             MLRtrainRMSE, MLRtestRMSE)
