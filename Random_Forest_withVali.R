###########################
#     RandomForest        #
#     Parallel            #
###########################

##### Packages
# Machine Learning
library(caret)
library(randomForest)
##### Foreach and multicore for running the script in parallel
library(foreach)
library(doMC)
##### Statistics
library(DescTools)
library(quantchem)
##### Plotting
library(ggplot2)

##### Load the preprocessed and preselected data here
##### Type the path to datafile or use

setwd("~/wz/ZWAge/Age_Sex/Data/Buffy_Coat/BuffyCoat_IntegratedData")

##### Set the amount of cores to uses
registerDoMC(cores = 4)

##### SET PARAMETERS
s1 <- Sys.time()
iters <- 50
ntrees <- 501
set.seed(100)

##### structures to save results of parallel threads
RMSE <- matrix(ncol = 5)
colnames(RMSE) <- c('a','b','c','d','e')
RMSE.outer <- c()
models <- list()
plots <- list()
feat.sel <- c()
RMSE.nr.feats <- c()

##### Make three folds
threeFCV <- createFolds(Buffy_Coat_Y_Transposed$Age, k = 3, list = TRUE)

##### Foreach loop with dopar, the iterations are done parallel
vali <- foreach (it =1:iters, .multicombine = F ) %dopar% {
  
  for(j in 1:3){
    
  ##### Assign Validation & Traintest set
  validation <- as.data.frame(Buffy_Coat_Y_Transposed [threeFCV[[j]],])
  trainTest <- as.data.frame(Buffy_Coat_Y_Transposed [-threeFCV[[j]],])
  
  
  ##### Create ten folds
  tenFCV <-  createFolds(y=trainTest$Age, k = 10, list = TRUE) 
  
  ##### Vectors for save the results
  feats <-c(rep(x = 0, times = dim(Buffy_Coat_Y_Transposed )[2]))
  names(feats) <- colnames(Buffy_Coat_Y_Transposed [,-1])
  plot.tune <- list()
  RMSE.inner <- matrix(nrow = 10, ncol = 2)
  n.feats <- c()
  
  for (inner in 1:10){
    print(inner)
    test <- as.data.frame(trainTest[tenFCV[[inner]],])
    train <- as.data.frame(trainTest[-tenFCV[[inner]],])
    # Repeat RF 10 times to rank the variables
    # Give points to the 100 most important features of each model
    # 100 points for #1 and 1 point for #100
    
    rf <- randomForest(x = train[,-1],
                       y = train[,1],
                       do.trace = T,
                       ntree = ntrees,
                       replace = TRUE,
                       importance = TRUE,
                       localImp = FALSE,
                       nPerm = 3,
                       corr.bias = FALSE,
                       keep.inbag = FALSE)
    ImpCpGs <- names(sort(rf$importance[,2], decreasing = T))
    
    # Count the points
    for(p in 1:100){
      i<- rev(ImpCpGs[1:100])[p]
      feats[i] <- feats[i] +p
    }
    
    feat.sel <- cbind(feat.sel, feats)
    
    # Order the important features vector and select the 100 probes with highest score
    Imp <- names(feats[order(feats, decreasing = T)])[1:150]
    
    res <- matrix(nrow =47, ncol = 3)
    res[,1] <- seq(4,50)
    colnames(res) <- c("features", "RMSE.in", "cor")
    for (t in seq(4,50)) {
      
      # RF train
      rf <- randomForest(x = train[,colnames(train) %in% Imp[1:t]],
                         y = train[,1],
                         do.trace = FALSE,
                         xtest = NULL,
                         ytest = NULL,
                         ntree = ntrees,
                         replace = TRUE,
                         importance = TRUE,
                         localImp = FALSE,
                         nPerm = 3,
                         norm.votes = TRUE,
                         corr.bias = FALSE,
                         keep.inbag = FALSE)
      # Save OOB in matrix and mtry
      
      # Predict the samples in folds that are not in training
      pred <- predict(rf, test[colnames(test) %in% Imp[1:t]])
      cor <- cor(x = pred, test$Age)
      RMSE.in <- round(sqrt(mean((pred - test$Age)^2)), 2)
      
      # Save RMSE and CC
      res[res[,1] == t , 2:3] <- c(RMSE.in, cor)
    }
    res <- as.data.frame(res)
    
    # Bind it to the number of features that was used in the model into one matrix
    n.feats <- rbind(n.feats, res)
    
    # Save the plots of RMSE over no.variables
    p1 <- ggplot(res, aes(x = features, y = RMSE.in)) + 
      geom_point() +
      xlab("# features") +
      ylab("RMSE(Years)") +
      geom_smooth(method = "lm", formula = y ~ poly(x,2), colour = "green")
    
    # Make model, deduce amount of features with lowest RMSE
    lm.RMSE <- lm(RMSE.in ~ features + I(features^2), data = res)
    deriv <- derivative(lm.RMSE, 1:150)
    min.RMSE <- which(deriv == deriv[deriv >0][1], arr.ind = T)
    to.RMSE.inner <- lm.RMSE$coefficients[1] + lm.RMSE$coefficients[2] * min.RMSE
    RMSE.inner[inner,] <- c (to.RMSE.inner, min.RMSE)
    plot.tune[[inner]] <- p1
  }
  
  RMSE.inner <- RMSE.inner[RMSE.inner[,2] !=1,]
  
  # Select CpG for the trainTest set model
  ImpCpGs.sel <- Imp[1:median(RMSE.inner[,2], na.rm=TRUE)]
  
  # RF on trainTest with optimal amount of features
  rfFinal <- randomForest(x = trainTest[,colnames(trainTest) %in% ImpCpGs.sel],
                          y = trainTest[,1],
                          do.trace = FALSE,
                          ntree = ntrees,
                          replace = TRUE,
                          importance = TRUE,
                          localImp = FALSE,
                          nPerm = 3,
                          norm.votes = TRUE,
                          corr.bias = FALSE,
                          keep.inbag = FALSE)
  
  groups <- data.frame(a = c(1:20), b = c(21:40), c = c(41:60), d = c(61:80), e = c(81:100))
  
  # Predict validation set ages be test set model
  for (g in colnames(groups)) {
    selvec <- validation[,1] %in% as.numeric(groups[,g])
    if(all(selvec == FALSE) == TRUE) {
      next}
    j.pred <- predict(rfFinal, validation[selvec, colnames(validation) %in% ImpCpGs.sel])
    
    # Calculate RMSE
    j.error <- validation$Age[selvec] - j.pred
    j.RMSE <- sqrt(mean(j.error^2))
    j.MAD <- MeanAD(j.error)
    RMSE[,g] <- j.RMSE
    pred_j <- predict(rfFinal, validation)
  }
  
  # Save models, plots, RMSE, #features
  models[[j]] <- rfFinal
  plots[[j]] <- plot.tune
  RMSE.outer <- rbind(RMSE.outer, RMSE)
  RMSE.nr.feats <- rbind(RMSE.nr.feats, n.feats)
}
  result <- list(RMSE.outer, models, plots, feat.sel, RMSE.nr.feats)
}

cat("Time :\t")
cat(Sys.time() - s1)
cat("\n")
