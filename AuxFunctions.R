######################################################
#
#               Observation Weights
#
######################################################

ObsWeights <- function(y, method = c("equal",
                                     "zscore",
                                     "sigmoid", 
                                     "curtailed"),
                       tail = c("left", "right", "both"),
                       left.cut = quantile(y, 1/4),
                       right.cut = quantile(y, 3/4),
                       A = 0, k = 1) {
  n <- length(y)
  #==================================
  #   Error Checking
  #==================================
  if((tail != "left") & (tail != "right") & (tail != "both")) {
    stop("tail must be one of - left, right or both")
  }

  
  #==================================
  #   Equal weights
  #==================================
  
  if(method == "equal") {
    w <- rep(1,n)
  }
  
  #==================================
  #   Zscore weights
  #==================================
  
  if(method == "zscore" & tail == "both") {
    w <- abs(y-mean(y))/sd(y)
    
    w <- w^k
    w <- (w/sum(w))*n
  }
  
  if(method == "zscore" & tail == "left") {
    w <- abs(y-mean(y))/sd(y)

    ylarge.ind <- which(y >  left.cut)
    w[ylarge.ind] <- abs(left.cut-mean(y))/sd(y)
    
    w <- w^k
    w <- (w/sum(w))*n
  }
  
  if(method == "zscore" & tail == "right") {
    w <- abs(y-mean(y))/sd(y)
    
    ysmall.ind <- which(y <  right.cut)
    w[ysmall.ind] <- abs(right.cut-mean(y))/sd(y)
    
    w <- w^k
    w <- (w/sum(w))*n
  }
  
  #==================================
  #   Sigmoid weights
  #==================================
  
  if(method == "sigmoid" & tail == "both") {
    w <- 1/(1+ exp(-abs(y-mean(y))))
    
    w <- w^k
    w <- (w/sum(w))*n
  }
  
  if(method == "sigmoid" & tail == "left") {
    w <- 1/(1+ exp(-abs(y-mean(y))))
 
    ylarge.ind <- which(y >  left.cut)
    w[ylarge.ind] <- 1/(1+ exp(-abs(left.cut-mean(y))))
    
    w <- w^k
    w <- (w/sum(w))*n
  }
  
  if(method == "sigmoid" & tail == "right") {
    w <- 1/(1 + exp(-abs(y-mean(y))))
    
    ysmall.ind <- which(y < right.cut)
    w[ysmall.ind] <- 1/(1 + exp(-abs(right.cut-mean(y))))
    
    w <- w^k
    w <- (w/sum(w))*n
  }
  
  #==================================
  #   Curtailed weights
  #==================================
  
  if(method == "curtailed" & tail == "both" ) {

    w <- rep(1, n)
    ysmall.ind <- which(y <= left.cut)
    ylarge.ind <- which(y >= right.cut)
    
    for(j in 1:length(ysmall.ind)) {
      w[ysmall.ind[j]] <- A + exp(1 + abs(y[ysmall.ind[j]]-left.cut))
    }
    for(j in 1:length(ylarge.ind)) {
      w[ylarge.ind[j]] <- A + exp(1 + abs(y[ylarge.ind[j]]-right.cut))
    }
    
    w <- w^k
    w <- (w/sum(w))*n
  }
  
  if(method == "curtailed" & tail == "left" ) {
    w <- rep(1, n)                               # 0 for drug 6 doesn't improve by much
    ysmall.ind <- which(y <= left.cut)
    
    for(j in 1:length(ysmall.ind)) {
      w[ysmall.ind[j]] <- A + (exp(1 + abs(y[ysmall.ind[j]]-left.cut)))
    }
    
    w <- w^k
    w <- (w/sum(w))*n
  }
  
  if(method == "curtailed" & tail == "right" ) {
    w <- rep(1, n)
    ylarge.ind <- which(y >= right.cut)
    
    for(j in 1:length(ylarge.ind)) {
      w[ylarge.ind[j]] <- A + exp(1 + abs(y[ylarge.ind[j]]-right.cut))
    }
    
    w <- w^k
    w <- (w/sum(w))*n
  }
  
  return(w)
}


######################################################
#
#                RMSE Measure
#
######################################################

rmse <- function(y, yhat, direction = c("left", "right", "both"), 
                 left.cut = quantile(y, 1/4),
                 right.cut = quantile(y, 3/4)) {
  if(direction == "left") {
    lefttailed.ind <- which((y <= left.cut))
    lefttailed.n <- length(lefttailed.ind )
    SS <- sum((y[lefttailed.ind] - yhat[lefttailed.ind])^2)
    rmse <- sqrt(SS/lefttailed.n)
  }
  
  if(direction == "right") {
    righttailed.ind <- which((y >= right.cut))
    righttailed.n <- length(righttailed.ind )
    SS <- sum((y[righttailed.ind] - yhat[righttailed.ind])^2)
    rmse <- sqrt(SS/righttailed.n)
  }
  
  if(direction == "both") {
    twotailed.ind <- which((y <= left.cut) | (y >= right.cut))
    twotailed.n <- length(twotailed.ind )
    SS <- sum((y[twotailed.ind] - yhat[twotailed.ind])^2)
    rmse <- sqrt(SS/twotailed.n)
  }
  
  return(rmse)
}


######################################################
#
#               Histogram Plots
#
######################################################


HistPlot <- function(drugind, data ) {
  dvcl <- as.numeric(data[drugind,-which(is.na(data[drugind,]))])
  df <- data.frame(dvcl = dvcl)
  
  p1 <- ggplot(df, aes(x=dvcl)) + 
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=.2,
                   colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") +
    ggtitle("Histogram + Density Curve") +
    labs(x=paste("Dose-Response AUC for Drug: ", drugind),y="Density")
  
  
  p2 <- ggplot(df, aes(x=dvcl)) + 
    geom_histogram(aes(y=..count..),      # Histogram with density instead of count on y-axis
                   binwidth=.2,
                   colour="black", fill="skyblue")+
    ggtitle("Histogram of Count") +
    labs(x=paste("Dose-Response AUC for Drug: ", drugind), y="Count")
  grid.arrange(p2, p1, ncol=2)
}

######################################################
#
#               Density Plots
#
######################################################

DensityPlot <- function(y, yhat, drugind) {
  x <- data.frame(Original.TestY = y, Predicted.TestY = yhat)
  data<- melt(x)
  p1 <- ggplot(data, aes(x=value, fill=variable)) + geom_density(alpha=0.55) +
        ylim(0, 1) +
        ggtitle(paste("Density Plots for Drug ", drugind, sep=""))
  
  p2 <- ggplot(data,aes(x=variable, y=value, fill=variable)) + geom_boxplot() +
        ggtitle(paste("Box Plots for Drug ", drugind, sep=""))
  grid.arrange(p2, p1, ncol=2)
}