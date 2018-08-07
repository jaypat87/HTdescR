#Random Split for Training and Test Sets

library(caTools)
library(svMisc)
library(tidyr)
library(tibble)
library(dplyr)
library(FNN)
library(pls)
library(randomForest)
library(stats)
library(e1701)

devtools::use_data(dataacidester, internal = TRUE) #Keep this?

moddata <- dataacidester
moddata <- dplyr::rename(moddata, rate = log.rate.exp)
moddata <- dplyr::mutate(moddata, log.rate.exp = log10(rate))
moddata <- moddata[c(8,2:3,5:7)]


set.seed(529)
sample <- caTools::sample.split(moddata$log.rate.exp, SplitRatio = 0.8)
moddataTrain <- subset(moddata, sample == TRUE)
moddataTest <- subset(moddata, sample == FALSE)

TrainingSet <- moddataTrain
TestSet <- moddataTest
TrainingDesc <- TrainingSet[c(2:6)]
TestDesc <- TestSet[c(2:6)]
MTrain <- moddataTrain$log.rate.exp
MTest <- moddataTest$log.rate.exp


#MLR model

TrainEqMLR <- stats::lm(log.rate.exp~., data = TrainingSet)
TrainYMLR <- stats::predict(TrainEqMLR, newdata = TrainingDesc)
Trainr2MLR <- stats::lm(TrainYMLR ~ MTrain)
sumTrainMLR <- summary(Trainr2MLR)
multR2MLR <- as.numeric(sumTrainMLR[8])
adjR2MLR <- as.numeric(sumTrainMLR[9])
ResStErrMLR <- as.numeric(sumTrainMLR[6])

#TotalMLR <- dplyr::bind_cols(as.data.frame(TrainYMLR), TrainingSet)
#TotalMLR <- dplyr::mutate(TotalMLR, SE = (log.rate.exp - TrainYMLR)^2)
#TotalMLR <- dplyr::mutate(TotalMLR, SSE = sum(SE))
#TotalMLR <- dplyr::mutate(TotalMLR, meanmeas = mean(log.rate.exp))
#TotalMLR <- dplyr::mutate(TotalMLR, ST = (log.rate.exp - meanmeas)^2)
#TotalMLR <- dplyr::mutate(TotalMLR, SST = sum(ST))

TestYMLR <- stats::predict(TrainEqMLR, newdata = TestDesc)
TestMLR <- dplyr::bind_cols(as.data.frame(TestYMLR), TestSet)
TestMLR <- dplyr::mutate(TestMLR, SE = (log.rate.exp - TestYMLR)^2)
TestMLR <- dplyr::mutate(TestMLR, PRESS = sum(SE))
TestMLR <- dplyr::mutate(TestMLR, meanmeas = mean(log.rate.exp))
TestMLR <- dplyr::mutate(TestMLR, ST = (log.rate.exp - meanmeas)^2)
TestMLR <- dplyr::mutate(TestMLR, SST = sum(ST))
TestMLR <- dplyr::mutate(TestMLR, Q2 = (1-(PRESS/SST)))
MLRtestQ2 <- TestMLR$Q2[1]



#PLS model

TrainEqPLS <- pls::plsr(log.rate.exp~., data = TrainingSet, ncomp = 3, scale=TRUE)
PLSpredictions <- stats::predict(TrainEqPLS, newdata = TrainingDesc)
PLSpredDF <- as.data.frame(PLSpredictions)
TrainYPLS <- PLSpredDF$`log.rate.exp.3 comps`
Trainr2PLS <- stats::lm(TrainYPLS ~ MTrain)
sumTrainPLS <- summary(Trainr2PLS)
multR2PLS <- as.numeric(sumTrainPLS[8])
adjR2PLS <- as.numeric(sumTrainPLS[9])
ResStErrPLS <- as.numeric(sumTrainPLS[6])

#TotalPLS <- dplyr::bind_cols(as.data.frame(TrainYPLS), TrainingSet)
#TotalPLS <- dplyr::mutate(TotalPLS, SE = (log.rate.exp - TrainYPLS)^2)
#TotalPLS <- dplyr::mutate(TotalPLS, SSE = sum(SE))
#TotalPLS <- dplyr::mutate(TotalPLS, meanmeas = mean(log.rate.exp))
#TotalPLS <- dplyr::mutate(TotalPLS, ST = (log.rate.exp - meanmeas)^2)
#TotalPLS <- dplyr::mutate(TotalPLS, SST = sum(ST))

PLSTestPred <- stats::predict(TrainEqPLS, newdata = TestDesc)
PLSTestPredDF <- as.data.frame(PLSTestPred)
TestYPLS <- PLSTestPredDF$`log.rate.exp.3 comps`
TestPLS <- dplyr::bind_cols(as.data.frame(TestYPLS), TestSet)
TestPLS <- dplyr::mutate(TestPLS, SE = (log.rate.exp - TestYPLS)^2)
TestPLS <- dplyr::mutate(TestPLS, PRESS = sum(SE))
TestPLS <- dplyr::mutate(TestPLS, meanmeas = mean(log.rate.exp))
TestPLS <- dplyr::mutate(TestPLS, ST = (log.rate.exp - meanmeas)^2)
TestPLS <- dplyr::mutate(TestPLS, SST = sum(ST))
TestPLS <- dplyr::mutate(TestPLS, Q2 = (1-(PRESS/SST)))
PLStestQ2 <- TestPLS$Q2[1]



#RF model

TrainEqRF <- randomForest::randomForest(log.rate.exp~., data = TrainingSet, ntree = 500, nodesize = 5,
                                         mtry = 5, importance = TRUE, na.action = na.omit)
TrainYRF <- stats::predict(TrainEqRF, newdata = TrainingDesc)
Trainr2RF <- stats::lm(TrainYRF ~ MTrain)
sumTrainRF <- summary(Trainr2RF)
multR2RF <- as.numeric(sumTrainRF[8])
adjR2RF <- as.numeric(sumTrainRF[9])
ResStErrRF <- as.numeric(sumTrainRF[6])

#TotalRF <- dplyr::bind_cols(as.data.frame(TrainYRF), TrainingSet)
#TotalRF <- dplyr::mutate(TotalRF, SE = (log.rate.exp - TrainYRF)^2)
#TotalRF <- dplyr::mutate(TotalRF, SSE = sum(SE))
#TotalRF <- dplyr::mutate(TotalRF, meanmeas = mean(log.rate.exp))
#TotalRF <- dplyr::mutate(TotalRF, ST = (log.rate.exp - meanmeas)^2)
#TotalRF <- dplyr::mutate(TotalRF, SST = sum(ST))

TestYRF <- stats::predict(TrainEqRF, newdata = TestDesc)
TestRF <- dplyr::bind_cols(as.data.frame(TestYRF), TestSet)
TestRF <- dplyr::mutate(TestRF, SE = (log.rate.exp - TestYRF)^2)
TestRF <- dplyr::mutate(TestRF, PRESS = sum(SE))
TestRF <- dplyr::mutate(TestRF, meanmeas = mean(log.rate.exp))
TestRF <- dplyr::mutate(TestRF, ST = (log.rate.exp - meanmeas)^2)
TestRF <- dplyr::mutate(TestRF, SST = sum(ST))
TestRF <- dplyr::mutate(TestRF, Q2 = (1-(PRESS/SST)))
RFtestQ2 <- TestRF$Q2[1]



#SVR model

TrainEqSVR <- e1071::svm(log.rate.exp~., data = TrainingSet, cost = 150, epsilon = 0.05, gamma = 0.00014)
TrainYSVR <- stats::predict(TrainEqSVR, newdata = TrainingDesc)
Trainr2SVR <- stats::lm(TrainYSVR ~ MTrain)
sumTrainSVR <- summary(Trainr2SVR)
multR2SVR <- as.numeric(sumTrainSVR[8])
adjR2SVR <- as.numeric(sumTrainSVR[9])
ResStErrSVR <- as.numeric(sumTrainSVR[6])

#TotalSVR <- dplyr::bind_cols(as.data.frame(TrainYSVR), TrainingSet)
#TotalSVR <- dplyr::mutate(TotalSVR, SE = (log.rate.exp - TrainYSVR)^2)
#TotalSVR <- dplyr::mutate(TotalSVR, SSE = sum(SE))
#TotalSVR <- dplyr::mutate(TotalSVR, meanmeas = mean(log.rate.exp))
#TotalSVR <- dplyr::mutate(TotalSVR, ST = (log.rate.exp - meanmeas)^2)
#TotalSVR <- dplyr::mutate(TotalSVR, SST = sum(ST))

TestYSVR <- stats::predict(TrainEqSVR, newdata = TestDesc)
TestSVR <- dplyr::bind_cols(as.data.frame(TestYSVR), TestSet)
TestSVR <- dplyr::mutate(TestSVR, SE = (log.rate.exp - TestYSVR)^2)
TestSVR <- dplyr::mutate(TestSVR, PRESS = sum(SE))
TestSVR <- dplyr::mutate(TestSVR, meanmeas = mean(log.rate.exp))
TestSVR <- dplyr::mutate(TestSVR, ST = (log.rate.exp - meanmeas)^2)
TestSVR <- dplyr::mutate(TestSVR, SST = sum(ST))
TestSVR <- dplyr::mutate(TestSVR, Q2 = (1-(PRESS/SST)))
SVRtestQ2 <- TestSVR$Q2[1]




#Compare all Q2

allTrainR2 <- data.frame(multR2MLR, multR2PLS, multR2RF, multR2SVR)
allTestQ2 <- data.frame(MLRtestQ2, PLStestQ2, RFtestQ2, SVRtestQ2)


#Calculate AD using kNN

AD1 <- FNN::knn.dist(TrainingSet, k=5, algorithm=c("kd_tree", "cover_tree", "CR", "brute"))
AD2 <- as.data.frame(AD1)
AD3 <- dplyr::mutate(AD2, Index = seq_len(56))
AD4 <- dplyr::mutate(AD3, avgDist = rowMeans(AD2))
AD4 <- dplyr::arrange(AD4, avgDist)

P95 <- stats::quantile(AD4$avgDist, prob = 0.95)



#Compare Test Set with AD

ADTestSet <- dplyr::rownames_to_column(TestSet)
ADTestSet <- dplyr::mutate(ADTestSet, loopnum = seq_len(14))
ADTrainSet <- dplyr::rownames_to_column(TrainingSet)
ADTrainSet <- dplyr::add_column(ADTrainSet, loopnum = 0)
TestSetAppDom <- data.frame(TestSetIndex=numeric(), ADLevel= numeric(), AD=character())


for(i in 1:14){

  TestChem <- dplyr::filter(ADTestSet, loopnum == i)
  TestChemInd <- TestChem$rowname
  CompareSet <- dplyr::bind_rows(ADTrainSet, TestChem)
  CompareSet <- CompareSet[c(2:7)]
  CompareAD1 <- FNN:knn.dist(CompareSet, k=5, algorithm=c("kd_tree", "cover_tree", "CR", "brute"))
  CompareAD1 <- as.data.frame(CompareAD1)
  CompareAD2 <- dplyr::mutate(CompareAD1, avgDist = rowMeans(CompareAD1))
  TestChemAD <- CompareAD2$avgDist[57]

  if(TestChemAD <= P95){
    TestSetAppDom <- dplyr::add_row(TestSetAppDom, TestSetIndex = TestChemInd, ADLevel = TestChemAD, AD = "IN")
  }
  else {
    TestSetAppDom <- dplyr::add_row(TestSetAppDom, TestSetIndex = TestchemInd, ADLevel = TestchemAD, AD = "OUt")
  }
}

#

ADTrain <- ADTrainSet[c(3:8)]


#Store allmodels

devtools::use_data(dataacidester, TrainEqMLR, TrainEqPLS, TrainEqRF, TrainEqSVR,
                   htdescHelper, TestSet, TrainingSet, P95, ADTestSet, ADTrain,
                   internal = TRUE, overwrite = TRUE)






