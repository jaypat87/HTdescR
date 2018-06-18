
getr2 <- function(predictedY, MeasuredY) {
  r2<-lm(predictedY[,,2] ~ MeasuredY)
  summary(r2) #Is this needed? Should summary be returned instead?
  return(r2)
}
