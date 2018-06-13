
getr2 <- function(predictedY, dataset) {
  r2<-lm(predictedY[,,2] ~ dataset = MeasuredYTest)
  summary(r2)

}
