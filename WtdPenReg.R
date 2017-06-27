######################################################
#
#              Weighted Penalized Regression
#
######################################################

WtdPenReg <- function(y.train, X.train, y.test, X.test, 
                      cv4alpha = FALSE, alpha = 0.2, A = 0, k = 1,
                      wt.method = c("equal", "zscore", "sigmoid", 
                                    "expsigmoid", "curtailed"),
                      tailToWeight = c("left","right","both"),
                      left.cut.train = quantile(y.train, 1/4),
                      right.cut.train = quantile(y.train, 3/4),
                      left.cut.test = quantile(y.test, 1/4),
                      right.cut.test = quantile(y.test, 3/4),
                      nfolds4alpha = ifelse(cv4alpha, 5, NA),
                      nfolds4lambda = 10,
                      print.out = TRUE) {

  #############################################
  trnsize <- length(y.train)
  tstsize <- length(y.test)
  
  
  ##########      Creating Weights     #########
  
  wt.y.train <- ObsWeights(y = y.train, method = wt.method, 
                         tail = tailToWeight,
                         left.cut = left.cut.train,
                         right.cut = right.cut.train,
                         A = A, k = k)
  
  #############################################
  
  if(cv4alpha) {
    MyTrainControl = trainControl(method = "cv",
                                  number = nfolds4alpha, 
                                  returnResamp = "none")
    
    t <- proc.time()
    EN.mod.cv <- train(x = X.train, y = y.train, 
                       method='glmnet',
                       metric = "RMSE",
                       weights = wt.y.train,
                       tuneGrid = expand.grid(.alpha = seq(0, 1, by = 0.1),
                                              .lambda = seq(0.01, 2, by=0.01)),
                       trControl = MyTrainControl)
    
    timeTaken <- proc.time()-t
    cat("*Cross Validation for alpha & lambda complete.\n")
    opt.EN.mod <- glmnet(x = X.train, y = y.train,
                         family="gaussian",
                         lambda = EN.mod.cv$bestTune["lambda"],
                         weights = wt.y.train,
                         alpha = EN.mod.cv$bestTune["alpha"],
                         intercept = TRUE,
                         standardize = TRUE)
    
    cat("*Final model fit complete.\n")
    alpha <- unlist(EN.mod.cv$bestTune["alpha"])
    lambda <- EN.mod.cv$bestTune["lambda"]
  } else {
    t <- proc.time()
    EN.mod.cv <- cv.glmnet(x = X.train, y = y.train,
                           nfolds = nfolds4lambda,
                           family="gaussian",
                           weights = wt.y.train,
                           type.measure = "mse",
                           alpha = alpha,
                           intercept = TRUE,
                           standardize = TRUE)
    timeTaken <- proc.time()-t
    
    cat("*Cross Validation for lambda complete.\n")
    opt.EN.mod <- glmnet(x = X.train, y = y.train,
                         family="gaussian",
                         weights = wt.y.train,
                         lambda = EN.mod.cv$lambda.min,
                         alpha = alpha,
                         intercept = TRUE,
                         standardize = TRUE)
    cat("*Final model fit complete.\n")
    lambda <- EN.mod.cv$lambda.min
  }
  
  optBeta <- opt.EN.mod$beta
  sparsity <- length(which(optBeta!=0))
  optIntcept <- opt.EN.mod$a0
  
  ####### Prediction for Test Data ##########
  
  yhat.test <- predict(opt.EN.mod, newx = X.test)
  colnames(yhat.test) <- NULL
  
  # RMSE results
  
  rmse.all <- sqrt(sum((yhat.test-y.test)^2))/sqrt(length(y.test))
  
  rmse.left <- rmse(y = y.test, yhat = yhat.test, 
                    direction = "left",
                    left.cut = left.cut.test,
                    right.cut = right.cut.test)
  
  rmse.right <- rmse(y = y.test, yhat = yhat.test, 
                     direction = "right", 
                    left.cut = left.cut.test,
                    right.cut = right.cut.test)
  
  rmse.both <- rmse(y = y.test, yhat = yhat.test, 
                     direction = "both", 
                     left.cut = left.cut.test,
                     right.cut = right.cut.test)
  
  if(print.out){
    cat("\n====================================================",
        "\nWeighting Method        :", paste(wt.method, ", k=",k, sep=""),
        "\nWhich tail is weighted  :", ifelse(wt.method != "equal", tailToWeight, "None"),
        "\nCV for alpha done?      :", cv4alpha,
        "\n----------------------------------------------------",
        "\nRMSE - Total                :", round(rmse.all, 4), 
        "\nRMSE - Left Tail            :", round(rmse.left, 4),
        "\nRMSE - Right Tail           :", round(rmse.right, 4),
        "\nRMSE - Both Tails           :", round(rmse.both, 4),
        "\n----------------------------------------------------",
#         "\nLeft Cut off (Train)        :", round(left.cut.train, 4),
#         "\nRight Cut off (Train)       :", round(right.cut.train, 4),
#         "\nLeft Cut off (Test)         :", round(left.cut.test, 4),
#         "\nRight Cut off (Test)        :", round(right.cut.test, 4),
#         "\n----------------------------------------------------",
        "\nSparsity                    :", sparsity, 
        "\nalpha                       :", round(alpha, 4), 
        "\nlambda                      :", round(lambda, 4),
        "\nTime Taken                  :", timeTaken["elapsed"], " seconds",
        "\n====================================================\n")
  }
 
  
  res <- list(yhat.test = yhat.test,
              finalENModel = opt.EN.mod,
              sparsity = sparsity,
              timeTaken = timeTaken["elapsed"],
              rmse.all = rmse.all,
              rmse.left = rmse.left,
              rmse.right = rmse.right,
              rmse.both = rmse.both
              )
  return(res)
}

