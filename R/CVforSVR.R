#CV for SVR

library(purrr)
library(e1071)
library(dplyr)

mod <- function(cost, gamma, epsilon){
  set.seed(202)
  eq <- e1071::svm(log.rate.exp~., data = TrainingSet, cross = 10, cost = cost, gamma = gamma, epsilon = epsilon)
  r2 <- eq$scorrcoeff
  return(r2)
}


gs <- list(cost = c(1.0,50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0, 550.0, 600.0),
           gamma = c(1.0,5.0,0.1,0.5,0.05,0.01,0.005,0.001, 0.0005, 0.0001),
           epsilon = c(0.1,0.2,0.3,0.4,0.05,0.01,0.001,0.005)) %>%
  cross_df()


options(warn = -1)
fit <- apply(gs, 1, function(x) mod(x['cost'],x['gamma'],x['epsilon']))
r2 <- as.data.frame(fit)
gs <- dplyr::bind_cols(gs, r2)
gs <- dplyr::arrange(gs, desc(fit))
View(gs)



