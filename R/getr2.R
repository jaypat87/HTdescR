
getr2 <- function(predictedY, mdataset) {
  r2<-lm(predictedY[,,2] ~ mdataset = MeasuredYTest)
  summary(r2)
  return(r2)
}
