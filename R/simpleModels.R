#Random Split for Training and Test Sets

dataacidester <- read.csv("./data/dataacidester.csv") #Import data file
dataacidester <- dataacidester[c(2:8)] #Trim index numbers from data frame
moddata <- dataacidester
colnames(moddata)[1] <- "rate" #Renaming rate column
moddata <- dplyr::mutate(moddata, log.rate.exp = log10(rate)) #Making a column of log rates
moddata <- moddata[c(8,2:3,5:7)] #Sorting columns of rate and descriptors to make sense visually
moddata <- dplyr::slice(moddata, c(1:22,28:70)) #Cutting out outliers for testing

# EFSA conversionns

EFSACAE2 <- EFSACAE[c(5, 55, 56, 64, 65, 72)]
EFSACAE3 <- slice(EFSACAE2, 5:9)
EFSACAE4 <- mutate_at(EFSACAE3, vars(log.rate.exp), funs(log(., 10)))

moddata2 <- bind_rows(moddata, EFSACAE4)
moddata3 <- tidyr::replace_na(moddata2, list(r2.hammett.value = 0))

moddata <- moddata3




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


#MLR model

TrainEqMLR <- stats::lm(log.rate.exp~., data = TrainingSet) #Training the MLR model on the Training Set
TrainYMLR <- stats::predict(TrainEqMLR, newdata = TrainingDesc) #Using the model to predict the Training Set
Trainr2MLR <- stats::lm(TrainYMLR ~ MTrain) #Making a simple regression of the experimental vs predicted rates for the Training Set
sumTrainMLR <- summary(Trainr2MLR)
multR2MLR <- as.numeric(sumTrainMLR[8]) #Grabbing mult. R2 for this model
adjR2MLR <- as.numeric(sumTrainMLR[9]) #Grabbing adj. R2 for this model
ResStErrMLR <- as.numeric(sumTrainMLR[6]) #Grabbing residual standard error for this model

#Old code for manually calculating R2
#TotalMLR <- dplyr::bind_cols(as.data.frame(TrainYMLR), TrainingSet)
#TotalMLR <- dplyr::mutate(TotalMLR, SE = (log.rate.exp - TrainYMLR)^2)
#TotalMLR <- dplyr::mutate(TotalMLR, SSE = sum(SE))
#TotalMLR <- dplyr::mutate(TotalMLR, meanmeas = mean(log.rate.exp))
#TotalMLR <- dplyr::mutate(TotalMLR, ST = (log.rate.exp - meanmeas)^2)
#TotalMLR <- dplyr::mutate(TotalMLR, SST = sum(ST))

#Getting Test R2 through manual calculation
TestYMLR <- stats::predict(TrainEqMLR, newdata = TestDesc) #Using the model to predict the Test Set
TestMLRreg <- stats::lm(TestYMLR ~ MTest)
sumTestMLR <- summary(TestMLRreg)
TestMLR <- dplyr::bind_cols(as.data.frame(TestYMLR), TestSet) #Pasting together the Test Set descriptors and predicted values into one dataframe
TestMLR <- dplyr::mutate(TestMLR, SE = (log.rate.exp - TestYMLR)^2)
TestMLR <- dplyr::mutate(TestMLR, PRESS = sum(SE))
TestMLR <- dplyr::mutate(TestMLR, meanmeas = mean(log.rate.exp))
TestMLR <- dplyr::mutate(TestMLR, ST = (log.rate.exp - meanmeas)^2)
TestMLR <- dplyr::mutate(TestMLR, SST = sum(ST))
TestMLR <- dplyr::mutate(TestMLR, Q2 = (1-(PRESS/SST)))
MLRtestR2 <- TestMLR$Q2[1]

#Getting Test R2 by regression (same way as Training R2) for comparison purposes
#Testr2MLR <- stats::lm(TestYMLR ~ MTest)
#sumTestMLR <- summary(Testr2MLR)
#TestMultR2MLR <- as.numeric(sumTestMLR[8])
#TestAdjR2MLR <- as.numeric(sumTestMLR[9])
#TestRSEMLR <- as.numeric(sumTestMLR[6])



#PLS model

TrainEqPLS <- pls::plsr(log.rate.exp~., data = TrainingSet, ncomp = 3, scale=TRUE) #Training the PLS model on the Training Set
PLSpredictions <- stats::predict(TrainEqPLS, newdata = TrainingDesc) #Using the model to predict the Training Set
PLSpredDF <- as.data.frame(PLSpredictions)
TrainYPLS <- PLSpredDF$`log.rate.exp.3 comps` #Pulling out the predicted values for the Training Set
Trainr2PLS <- stats::lm(TrainYPLS ~ MTrain) #Making a simple regression of the experimental vs predicted values for the Training Set
sumTrainPLS <- summary(Trainr2PLS)
multR2PLS <- as.numeric(sumTrainPLS[8]) #Grabbing the mult. R2 for this model
adjR2PLS <- as.numeric(sumTrainPLS[9]) #Grabbing the adj. R2 for this model
ResStErrPLS <- as.numeric(sumTrainPLS[6]) #Grabbing the residual standard error for this model

#Old code for manually calculating R2
#TotalPLS <- dplyr::bind_cols(as.data.frame(TrainYPLS), TrainingSet)
#TotalPLS <- dplyr::mutate(TotalPLS, SE = (log.rate.exp - TrainYPLS)^2)
#TotalPLS <- dplyr::mutate(TotalPLS, SSE = sum(SE))
#TotalPLS <- dplyr::mutate(TotalPLS, meanmeas = mean(log.rate.exp))
#TotalPLS <- dplyr::mutate(TotalPLS, ST = (log.rate.exp - meanmeas)^2)
#TotalPLS <- dplyr::mutate(TotalPLS, SST = sum(ST))

#Getting Test R2 through manual calculation
PLSTestPred <- stats::predict(TrainEqPLS, newdata = TestDesc) #Using the model to predict the Test Set
PLSTestPredDF <- as.data.frame(PLSTestPred)
TestYPLS <- PLSTestPredDF$`log.rate.exp.3 comps`
TestPLSreg <- stats::lm(TestYPLS ~ MTest)
sumTestPLS <- summary(TestPLSreg)
TestPLS <- dplyr::bind_cols(as.data.frame(TestYPLS), TestSet)
TestPLS <- dplyr::mutate(TestPLS, SE = (log.rate.exp - TestYPLS)^2)
TestPLS <- dplyr::mutate(TestPLS, PRESS = sum(SE))
TestPLS <- dplyr::mutate(TestPLS, meanmeas = mean(log.rate.exp))
TestPLS <- dplyr::mutate(TestPLS, ST = (log.rate.exp - meanmeas)^2)
TestPLS <- dplyr::mutate(TestPLS, SST = sum(ST))
TestPLS <- dplyr::mutate(TestPLS, Q2 = (1-(PRESS/SST)))
PLStestR2 <- TestPLS$Q2[1]



#RF model

TrainEqRF <- randomForest::randomForest(log.rate.exp~., data = TrainingSet, ntree = 500, nodesize = 5,
                                         mtry = 5, importance = TRUE, na.action = na.omit) #Training the RF model on the Training Set
TrainYRF <- stats::predict(TrainEqRF, newdata = TrainingDesc) #Using the model to predict the Training Set
Trainr2RF <- stats::lm(TrainYRF ~ MTrain) #Making a simple regression of the experimental vs predicted values for the Training Set
sumTrainRF <- summary(Trainr2RF)
multR2RF <- as.numeric(sumTrainRF[8]) #Grabbing the mult. R2 for this model
adjR2RF <- as.numeric(sumTrainRF[9]) #Grabbing the adj. R2 for this model
ResStErrRF <- as.numeric(sumTrainRF[6]) #Grabbing the residual standard error for this model

#Old code for manually calculating R2
#TotalRF <- dplyr::bind_cols(as.data.frame(TrainYRF), TrainingSet)
#TotalRF <- dplyr::mutate(TotalRF, SE = (log.rate.exp - TrainYRF)^2)
#TotalRF <- dplyr::mutate(TotalRF, SSE = sum(SE))
#TotalRF <- dplyr::mutate(TotalRF, meanmeas = mean(log.rate.exp))
#TotalRF <- dplyr::mutate(TotalRF, ST = (log.rate.exp - meanmeas)^2)
#TotalRF <- dplyr::mutate(TotalRF, SST = sum(ST))

#Getting Test R2 through manual calculation
TestYRF <- stats::predict(TrainEqRF, newdata = TestDesc) #Using the model to predict the Test Set
TestRFreg <- stats::lm(TestYRF ~ MTest)
sumTestRF <- summary(TestRFreg)
TestRF <- dplyr::bind_cols(as.data.frame(TestYRF), TestSet)
TestRF <- dplyr::mutate(TestRF, SE = (log.rate.exp - TestYRF)^2)
TestRF <- dplyr::mutate(TestRF, PRESS = sum(SE))
TestRF <- dplyr::mutate(TestRF, meanmeas = mean(log.rate.exp))
TestRF <- dplyr::mutate(TestRF, ST = (log.rate.exp - meanmeas)^2)
TestRF <- dplyr::mutate(TestRF, SST = sum(ST))
TestRF <- dplyr::mutate(TestRF, Q2 = (1-(PRESS/SST)))
RFtestR2 <- TestRF$Q2[1]



#SVR model

TrainEqSVR <- e1071::svm(log.rate.exp~., data = TrainingSet, cost = 75.5, epsilon = 0.033, gamma = 0.9970) #Training the SVR model on the Training Set
TrainYSVR <- stats::predict(TrainEqSVR, newdata = TrainingDesc) #Using the model to predict the Training Set
Trainr2SVR <- stats::lm(TrainYSVR ~ MTrain) #Making a simple regression of the experimental vs predicted values for the Training Set
sumTrainSVR <- summary(Trainr2SVR)
multR2SVR <- as.numeric(sumTrainSVR[8]) #Grabbing the mult. R2 for this model
adjR2SVR <- as.numeric(sumTrainSVR[9]) #Grabbing the adj. R2 for this model
ResStErrSVR <- as.numeric(sumTrainSVR[6]) #Grabbing the residual standard error for this model

#Old code for manually calculating R2
#TotalSVR <- dplyr::bind_cols(as.data.frame(TrainYSVR), TrainingSet)
#TotalSVR <- dplyr::mutate(TotalSVR, SE = (log.rate.exp - TrainYSVR)^2)
#TotalSVR <- dplyr::mutate(TotalSVR, SSE = sum(SE))
#TotalSVR <- dplyr::mutate(TotalSVR, meanmeas = mean(log.rate.exp))
#TotalSVR <- dplyr::mutate(TotalSVR, ST = (log.rate.exp - meanmeas)^2)
#TotalSVR <- dplyr::mutate(TotalSVR, SST = sum(ST))

#Getting Test R2 through manual calculation
TestYSVR <- stats::predict(TrainEqSVR, newdata = TestDesc) #Using the model to predict the Test Set
TestSVRreg <- stats::lm(TestYSVR ~ MTest)
sumTestSVR <- summary(TestSVRreg)
TestSVR <- dplyr::bind_cols(as.data.frame(TestYSVR), TestSet)
TestSVR <- dplyr::mutate(TestSVR, SE = (log.rate.exp - TestYSVR)^2)
TestSVR <- dplyr::mutate(TestSVR, PRESS = sum(SE))
TestSVR <- dplyr::mutate(TestSVR, meanmeas = mean(log.rate.exp))
TestSVR <- dplyr::mutate(TestSVR, ST = (log.rate.exp - meanmeas)^2)
TestSVR <- dplyr::mutate(TestSVR, SST = sum(ST))
TestSVR <- dplyr::mutate(TestSVR, Q2 = (1-(PRESS/SST)))
SVRtestR2 <- TestSVR$Q2[1]




#Compare all R2

allTrainMultR2 <- data.frame(multR2MLR, multR2PLS, multR2RF, multR2SVR) #Pasting together the different models' R2s for comparison
allTrainAdjR2 <- data.frame(adjR2MLR, adjR2PLS, adjR2RF, adjR2SVR)
allTrainRSE <- data.frame(ResStErrMLR, ResStErrPLS, ResStErrRF, ResStErrSVR)
allTestR2 <- data.frame(MLRtestR2, PLStestR2, RFtestR2, SVRtestR2)


#Calculate AD using kNN

AD1 <- FNN::knn.dist(TrainingSet, k=5, algorithm=c("kd_tree", "cover_tree", "CR", "brute")) #Calculating the five nearest neighbors' distances for the Training Set
AD2 <- as.data.frame(AD1)
AD3 <- dplyr::mutate(AD2, Index = seq_len(56)) #Labeling rows with index numbers for sorting purposes
AD4 <- dplyr::mutate(AD3, avgDist = rowMeans(AD2)) #Finding the average KNN distance for each molecule (row)
AD4 <- dplyr::arrange(AD4, avgDist) #Sorting by the average distance

P95 <- stats::quantile(AD4$avgDist, prob = 0.95) #Finding the 95 percentile of the average KNN distance - this is the Applicability Domain Threshold



#Compare Test Set with AD

ADTestSet <- tibble::rownames_to_column(TestSet) #Preparing a dataframe of the Test Set for AD comparison
ADTestSet <- dplyr::mutate(ADTestSet, loopnum = seq_len(14))
ADTrainSet <- tibble::rownames_to_column(TrainingSet)
ADTrainSet <- tibble::add_column(ADTrainSet, loopnum = 0)
TestSetAppDom <- data.frame(TestSetIndex=numeric(), ADLevel= numeric(), AD=character()) #Initializing an empty dataframe for the Test Set AD loop

#Loop through Test Set to find which Test Set molecules are in/out of AD
for(i in 1:14){

  TestChem <- dplyr::filter(ADTestSet, loopnum == i)
  TestChemInd <- TestChem$rowname
  CompareSet <- dplyr::bind_rows(ADTrainSet, TestChem) #Adding the Test Set chemical in question to the training set to be able to find KNN
  CompareSet <- CompareSet[c(2:7)]
  CompareAD1 <- FNN::knn.dist(CompareSet, k=5, algorithm=c("kd_tree", "cover_tree", "CR", "brute")) #Finding the KNN distances of training set + test set chemical in question
  CompareAD1 <- as.data.frame(CompareAD1)
  CompareAD2 <- dplyr::mutate(CompareAD1, avgDist = rowMeans(CompareAD1)) #Averaging the KNN distances
  TestChemAD <- CompareAD2$avgDist[57] #Pulling out the average KNN distance for the Test Set chemical in question

  #Labeling Test Set chemical in question as In or Out of AD based on if above or below AD threshold
  if(TestChemAD <= P95){
    TestSetAppDom <- dplyr::add_row(TestSetAppDom, TestSetIndex = TestChemInd, ADLevel = TestChemAD, AD = "IN")
  }
  else {
    TestSetAppDom <- dplyr::add_row(TestSetAppDom, TestSetIndex = TestChemInd, ADLevel = TestChemAD, AD = "OUt")
  }
}

#

ADTrain <- ADTrainSet[c(3:8)]

#Store allmodels
htdescHelper <- read.csv("./data/htdescHelper.csv")
devtools::use_data(dataacidester, TrainEqMLR, TrainEqPLS, TrainEqRF, TrainEqSVR,
                   htdescHelper, TestSet, TrainingSet, P95, ADTestSet, ADTrain,
                   internal = TRUE, overwrite = TRUE)






