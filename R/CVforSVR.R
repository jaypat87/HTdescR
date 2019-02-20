#CV for SVR

library(purrr)
library(e1071)
library(dplyr)

#Initializing the function to tune parameters for
mod <- function(cost, gamma, epsilon){
  set.seed(202)
  eq <- e1071::svm(log.rate.exp~., data = TrainingSet, cross = 10, cost = cost, gamma = gamma, epsilon = epsilon)
  eqPred <- stats::predict(eq, newdata = TrainingDesc)
  eqTrainr2 <- stats::lm(eqPred ~ MTrain)
  eqSumTrain <- summary(eqTrainr2)
  eqMultR2 <- as.numeric(eqSumTrain[8])
  eqAdjR2 <- as.numeric(eqSumTrain[9])
  eqResStErr <- as.numeric(eqSumTrain[6])
  return(eqMultR2)
}

costRange <- seq(from = 0, to = 1000, by = 0.5) #Creating sequences of parameter values over range
gammaRange <- seq(from = 0, to = 1, by = 0.0001)
epsilonRange <- seq(from = 0, to = 1, by = 0.001)

costSample <- sample(costRange, 50) #Randomly sampling from the sequences created so that we can make a grid that isn't too big
gammaSample <- sample(gammaRange, 50)
epsilonSample <- sample(epsilonRange, 50)

gs <- list(cost = costSample,
           gamma = gammaSample,
           epsilon = epsilonSample) %>%
  cross_df() #Creating a grid of all the combinations of randomly sampled parameter values


options(warn = -1)
fit <- apply(gs, 1, function(x) mod(x['cost'],x['gamma'],x['epsilon'])) #Running the SVM function using all of the combinations of parameter values and returning the R2 for each run
fitted <- as.data.frame(fit) #Converting the R2 output to a data.frame for easier use
gs <- dplyr::bind_cols(gs, fitted) #Binding the R2 with the grid of combinations of parameter values used to find those R2s
gs <- dplyr::arrange(gs, desc(fit)) #Arranging the output by highest R2
View(gs)


##############################################


svmTest <- best.svm(x = log.rate.exp~., data = TrainingSet, tunecontrol = tune.control(sampling = "cross", cross = 10, best.model = TRUE))

