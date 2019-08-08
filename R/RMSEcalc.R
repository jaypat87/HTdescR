#MLR RMSE

MLRresiduals <- as.data.frame(sumTrainMLR$residuals)
MLRresiduals2 <- data.frame(MLRresiduals^2)
meanMLRresiduals2 <- mutate(MLRresiduals2, mean(MLRresiduals2$sumTrainMLR.residuals))
MLRtrainRMSE <- sqrt(meanMLRresiduals2$`mean(MLRresiduals2$sumTrainMLR.residuals)`[1])

MLRresidualsTest <- as.data.frame(sumTestMLR$residuals)
MLRresidualsTest2 <- data.frame(MLRresidualsTest^2)
meanMLRresidualsTest2 <- mutate(MLRresidualsTest2, mean(MLRresidualsTest2$sumTestMLR.residuals))
MLRtestRMSE <- sqrt(meanMLRresidualsTest2$`mean(MLRresidualsTest2$sumTestMLR.residuals`[1])



#PLS RMSE

PLSresiduals <- as.data.frame(sumTrainPLS$residuals)
PLSresiduals2 <- data.frame(PLSresiduals^2)
meanPLSresiduals2 <- mutate(PLSresiduals2, mean(PLSresiduals2$sumTrainPLS.residuals))
PLStrainRMSE <- sqrt(meanPLSresiduals2$`mean(PLSresiduals2$sumTrainPLS.residuals)`[1])

PLSresidualsTest <- as.data.frame(sumTestPLS$residuals)
PLSresidualsTest2 <- data.frame(PLSresidualsTest^2)
meanPLSresidualsTest2 <- mutate(PLSresidualsTest2, mean(PLSresidualsTest2$sumTestPLS.residuals))
PLStestRMSE <- sqrt(meanPLSresidualsTest2$`mean(PLSresidualsTest2$sumTestPLS.residuals`[1])



#RF RMSE

RFresiduals <- as.data.frame(sumTrainRF$residuals)
RFresiduals2 <- data.frame(RFresiduals^2)
meanRFresiduals2 <- mutate(RFresiduals2, mean(RFresiduals2$sumTrainRF.residuals))
RFtrainRMSE <- sqrt(meanRFresiduals2$`mean(RFresiduals2$sumTrainRF.residuals)`[1])

RFresidualsTest <- as.data.frame(sumTestRF$residuals)
RFresidualsTest2 <- data.frame(RFresidualsTest^2)
meanRFresidualsTest2 <- mutate(RFresidualsTest2, mean(RFresidualsTest2$sumTestRF.residuals))
RFtestRMSE <- sqrt(meanRFresidualsTest2$`mean(RFresidualsTest2$sumTestRF.residuals`[1])



#SVR RMSE

SVRresiduals <- as.data.frame(sumTrainSVR$residuals)
SVRresiduals2 <- data.frame(SVRresiduals^2)
meanSVRresiduals2 <- mutate(SVRresiduals2, mean(SVRresiduals2$sumTrainSVR.residuals))
SVRtrainRMSE <- sqrt(meanSVRresiduals2$`mean(SVRresiduals2$sumTrainSVR.residuals)`[1])

SVRresidualsTest <- as.data.frame(sumTestSVR$residuals)
SVRresidualsTest2 <- data.frame(SVRresidualsTest^2)
meanSVRresidualsTest2 <- mutate(SVRresidualsTest2, mean(SVRresidualsTest2$sumTestSVR.residuals))
SVRtestRMSE <- sqrt(meanSVRresidualsTest2$`mean(SVRresidualsTest2$sumTestSVR.residuals`[1])



#All together

dataRange <- range(moddata3$log.rate.exp)
Range <- abs(dataRange[1]) + dataRange[2]
RMSEthreshold <- Range / 10

TrainRMSEs <- data.frame(MLRtrainRMSE, PLStrainRMSE, RFtrainRMSE, SVRtrainRMSE)
TestRMSEs <- data.frame(MLRtestRMSE, PLStestRMSE, RFtestRMSE, SVRtestRMSE)

TrainRMSEs
TestRMSEs
RMSEthreshold

