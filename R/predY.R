
predY <- function(eq, newdata){
  predictedY <- predict(eq, newdata = newdata)
  return(predictedY)
}

