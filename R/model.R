##Data Prep

library(tidyr)
library(caret)
library(cvq2)
library(dplyr)
library(matrixpls)
library(kernlab)

#aesample <- aefilled[c(3,53:70)]
#aeready <- replace_na(aesample, replace = list(r1.taft.value = 0, r1.es.value = 0, r1.ind.value = 0, r1.meta1.value = 0, r1.meta2.value = 0, r1.ortho1.value = 0, r1.ortho2.value = 0, r1.para1.value = 0, r1.hammett.value = 0, r2.taft.value = 0, r2.es.value = 0, r2.ind.value = 0, r2.meta1.value = 0, r2.meta2.value = 0, r2.ortho1.value = 0, r2.ortho2.value = 0, r2.para1.value = 0, r2.hammett.value = 0))

devtools::use_data(dataacidester, internal = TRUE)

moddata <- dataacidester
moddata <- rename(moddata, rate = log.rate.exp)
moddata <- mutate(moddata, log.rate.exp = log10(rate))
moddata <- moddata[c(8,2:7)]



##Model Building and Cross Validation

options(warn=-1)
seeds <- c(1, 2)

##do svm here and send it to cvq2?
#q2CV <- cvq2::cvq2(modelData = moddata, formula = log.rate.exp ~ ., nFold = 5, extOut = TRUE)



#SVMR

set.seed(2)
SVMRmod <- caret::train(log.rate.exp ~ ., data = moddata, method = "svmRadial", metric = "Rsquared",
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
SVMRobs <- SVMRobs[c(7,1:6)]
SVMRtrain <- rename(SVMRtrain, meas = .outcome)
SVMRdata <- full_join(SVMRtrain, SVMRobs)
SVMRdata <- SVMRdata[c(2:7,1,8)]

SVMRdatm <- mutate(SVMRdata, SE = (meas - pred)^2)
SVMRdatm <- mutate(SVMRdatm, meanmeas = mean(meas))
SVMRdatm <- mutate(SVMRdatm, SR = (pred - meanmeas)^2)
SVMRdatm <- mutate(SVMRdatm, ST = (SR + SE))
SVMRdatm <- mutate(SVMRdatm, SSE = sum(SE))
SVMRdatm <- mutate(SVMRdatm, SSR = sum(SR))
SVMRdatm <- mutate(SVMRdatm, SST = (SSE + SSR))
SVMRdatm <- mutate(SVMRdatm, STcalc = (meas - meanmeas)^2)
SVMRdatm <- mutate(SVMRdatm, SSTcalc = sum(STcalc))
SVMRdatm <- mutate(SVMRdatm, R2div = (SSR/SSTcalc))
SVMRdatm <- mutate(SVMRdatm, R2sub = (1-(SSE/SSTcalc)))
SVMRdatm <- mutate(SVMRdatm, Q2eq = (1-(SSE/SSTcalc)))
SVMRdatm <- mutate(SVMRdatm, Q2other = (1-(SSE/SST)))




#SVMP

set.seed(2)
SVMPmod <- caret::train(log.rate.exp ~ ., data = moddata, method = "svmPoly", metric = "Rsquared",
                        trControl = trainControl(method = "cv", number = 5, summaryFunction = defaultSummary,
                                                 returnData = TRUE, returnResamp = "all", savePredictions = "final"))
SVMPresults <- SVMPmod$results
SVMPresample <- SVMPmod$resample
SVMPbestTune <- SVMPmod$bestTune
SVMPfinalModel <- SVMPmod$finalModel
SVMPpred <- SVMPmod$pred
SVMPpred2 <- arrange(SVMRpred, rowIndex)
SVMPtrain <- SVMPmod$trainingData

SVMPobs <- SVMPtrain[c(-1)]
SVMPpred3 <- select(SVMPpred2, pred)
SVMPobs <- bind_cols(SVMPobs, SVMPpred3)
#SVMPobs <- rename(SVMPobs, y = pred)
SVMPobs <- SVMPobs[c(7,1:6)]
SVMPtrain <- rename(SVMPtrain, meas = .outcome)
SVMPdata <- full_join(SVMPtrain, SVMPobs)
SVMPdata <- SVMPdata[c(2:7,1,8)]

SVMPdatm <- mutate(SVMPdata, SE = (meas - pred)^2)
SVMPdatm <- mutate(SVMPdatm, meanmeas = mean(meas))
SVMPdatm <- mutate(SVMPdatm, SR = (pred - meanmeas)^2)
SVMPdatm <- mutate(SVMPdatm, ST = (SR + SE))
SVMPdatm <- mutate(SVMPdatm, SSE = sum(SE))
SVMPdatm <- mutate(SVMPdatm, SSR = sum(SR))
SVMPdatm <- mutate(SVMPdatm, SST = (SSE + SSR))
SVMPdatm <- mutate(SVMPdatm, STcalc = (meas - meanmeas)^2)
SVMPdatm <- mutate(SVMPdatm, SSTcalc = sum(STcalc))
SVMPdatm <- mutate(SVMPdatm, R2div = (SSR/SST))
SVMPdatm <- mutate(SVMPdatm, R2sub = (1-(SSE/SST)))
SVMPdatm <- mutate(SVMPdatm, Q2eq = (1-(SSE/SSTcalc)))


#PLS

set.seed(2)
PLSmod <- caret::train(log.rate.exp ~ ., data = moddata, method = "pls", metric = "Rsquared",
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
PLSobs <- PLSobs[c(7,1:6)]
PLStrain <- rename(PLStrain, meas = .outcome)
PLSdata <- full_join(PLStrain, PLSobs)
PLSdata <- PLSdata[c(2:7,1,8)]

PLSdatm <- mutate(PLSdata, SE = (meas - pred)^2)
PLSdatm <- mutate(PLSdatm, meanmeas = mean(meas))
PLSdatm <- mutate(PLSdatm, SR = (pred - meanmeas)^2)
PLSdatm <- mutate(PLSdatm, ST = (SR + SE))
PLSdatm <- mutate(PLSdatm, SSE = sum(SE))
PLSdatm <- mutate(PLSdatm, SSR = sum(SR))
PLSdatm <- mutate(PLSdatm, SST = (SSE + SSR))
PLSdatm <- mutate(PLSdatm, STcalc = (meas - meanmeas)^2)
PLSdatm <- mutate(PLSdatm, SSTcalc = sum(STcalc))
PLSdatm <- mutate(PLSdatm, R2div = (SSR/SST))
PLSdatm <- mutate(PLSdatm, R2sub = (1-(SSE/SST)))
PLSdatm <- mutate(PLSdatm, Q2eq = (1-(SSE/SSTcalc)))



#RF

set.seed(2)
RFmod <- caret::train(log.rate.exp ~ ., data = moddata, method = "rf", metric = "Rsquared",
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
RFobs <- RFobs[c(7,1:6)]
RFtrain <- rename(RFtrain, meas = .outcome)
RFdata <- full_join(RFtrain, RFobs)
RFdata <- RFdata[c(2:7,1,8)]

RFdatm <- mutate(RFdata, SE = (meas - pred)^2)
RFdatm <- mutate(RFdatm, meanmeas = mean(meas))
RFdatm <- mutate(RFdatm, SR = (pred - meanmeas)^2)
RFdatm <- mutate(RFdatm, ST = (SR + SE))
RFdatm <- mutate(RFdatm, SSE = sum(SE))
RFdatm <- mutate(RFdatm, SSR = sum(SR))
RFdatm <- mutate(RFdatm, SST = (SSE + SSR))
RFdatm <- mutate(RFdatm, STcalc = (meas - meanmeas)^2)
RFdatm <- mutate(RFdatm, SSTcalc = sum(STcalc))
RFdatm <- mutate(RFdatm, R2div = (SSR/SST))
RFdatm <- mutate(RFdatm, R2sub = (1-(SSE/SST)))
RFdatm <- mutate(RFdatm, Q2eq = (1-(SSE/SSTcalc)))




#MLR

set.seed(2)
MLRmod <- caret::train(log.rate.exp ~ ., data = moddata, method = "lm", metric = "Rsquared",
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
MLRobs <- MLRobs[c(7,1:6)]
MLRtrain <- rename(MLRtrain, meas = .outcome)
MLRdata <- full_join(MLRtrain, MLRobs)
MLRdata <- MLRdata[c(2:7,1,8)]

MLRdatm <- mutate(MLRdata, SE = (meas - pred)^2)
MLRdatm <- mutate(MLRdatm, meanmeas = mean(meas))
MLRdatm <- mutate(MLRdatm, SR = (pred - meanmeas)^2)
MLRdatm <- mutate(MLRdatm, ST = (SR + SE))
MLRdatm <- mutate(MLRdatm, SSE = sum(SE))
MLRdatm <- mutate(MLRdatm, SSR = sum(SR))
MLRdatm <- mutate(MLRdatm, SST = (SSE + SSR))
MLRdatm <- mutate(MLRdatm, STcalc = (meas - meanmeas)^2)
MLRdatm <- mutate(MLRdatm, SSTcalc = sum(STcalc))
MLRdatm <- mutate(MLRdatm, R2div = (SSR/SST))
MLRdatm <- mutate(MLRdatm, R2sub = (1-(SSE/SST)))
MLRdatm <- mutate(MLRdatm, Q2eq = (1-(SSE/SSTcalc)))




#Output

MLRstats <- c("MLR", MLRdatm$R2div[1], MLRdatm$Q2eq[1])
PLSstats <- c("PLS", PLSdatm$R2div[1], PLSdatm$Q2eq[1])
RFstats <- c("RF", RFdatm$R2div[1], RFdatm$Q2eq[1])
SVMPstats <- c("SVMP", SVMPdatm$R2div[1], SVMPdatm$Q2eq[1])
SVMRstats <- c("SVMR", SVMRdatm$R2div[1], SVMRdatm$Q2eq[1])

type <- c("MLR", "PLS", "RF", "SVMP", "SVMR")
allR2 <- c(MLRdatm$R2div[1], PLSdatm$R2div[1], RFdatm$R2div[1], SVMPdatm$R2div[1], SVMRdatm$R2div[1])
allQ2 <- c(MLRdatm$Q2eq[1], PLSdatm$Q2eq[1], RFdatm$Q2eq[1], SVMPdatm$Q2eq[1], SVMRdatm$Q2eq[1])
allStats <- data.frame(type, allR2, allQ2)
allStats


##Final Models

allFinalModels <- list("SVMRmodel"=SVMRmod, "SVMPmodel"=SVMPmod,
                    "PLSmodel"=PLSmod, "RFmodel"=RFmod, "MLRmodel"=MLRmod)

extractedModels <- list("PLSmodel"=PLSmod, "RFmodel"=RFmod, "MLRmodel"=MLRmod)
SVMmodels <- list("SVMRmodel"=SVMRmod, "SVMPmodel"=SVMPmod)


#Save models here

devtools::use_data(dataacidester, allFinalModels, extractedModels, SVMmodels,
                   SVMPmod, SVMRmod, PLSmod, RFmod, MLRmod, internal = TRUE, overwrite = TRUE)

