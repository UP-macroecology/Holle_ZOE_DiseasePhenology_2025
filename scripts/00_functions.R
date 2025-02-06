# ZOE project 


#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                            00. Functions                               #
# ---------------------------------------------------------------------- #


#-------------------------------------------------------------------------------

# 1. Function for spatial thinning ---------------------------------------------


thin <- function(sf, thin_dist = 3000, runs = 1, ncores = 1){
  
  require(sf, quietly = TRUE)
  require(purrr, quietly = TRUE)
  require(furrr, quietly = TRUE)
  
  sample.vec <- function(x, ...) x[sample(length(x), ...)]
  
  sf_buffer <- st_buffer(sf, thin_dist)
  buff_int <- st_intersects(sf, sf_buffer) 
  buff_int <- setNames(buff_int, 1:length(buff_int))
  
  n_int <- map_dbl(buff_int, length)
  
  plan(multisession, workers = ncores)
  
  seeds <- sample.int(n = runs)
  results_runs <- future_map(seeds, function(i){
    
    set.seed(i)
    while (max(n_int) > 1) {
      max_neighbors <- names(which(n_int == max(n_int)))
      
      # remove point with max neighbors
      sampled_id <- sample.vec(max_neighbors, 1)
      
      pluck(buff_int, sampled_id) <- NULL
      buff_int <- map(buff_int, function(x) setdiff(x, as.numeric(sampled_id)))
      n_int <- map_dbl(buff_int, length)
    }
    
    unlist(buff_int) %>% unique()
    
  })
  
  lengths <- map_dbl(results_runs, length)
  
  selected_run <- results_runs[[sample.vec(which(lengths == max(lengths)), 1)]]
  
  out <- sf[selected_run,]
  
  out <- sf_to_df(out)[,3:4] %>%
    rename("lon" = "x", "lat" = "y")
  
  out
}



#-------------------------------------------------------------------------------

# 2. Function for variable selection -------------------------------------------

select07_cv <- function(X, y, kfold=5, family="binomial",univar="glm2", threshold=0.7, method="spearman", sequence=NULL, weights=rep(1, length(y)))
{
  require(mgcv)
  
  # X .. Matrix or data.frame containing the predictor variables
  # y .. vector of response variable
  # kfold .. number of folds for random cross-validation or vector with group assignments (indexing which data point belongs to which fold)
  # family .. a description of the error distribution and link function to be used in the model.
  # univar .. a character string indicating the regression method to be used for estimating univariate importance. Must be one of the strings "glm1", "glm2" (default), or "gam". "glm1" will estimate a generalised linear model (GLM) with a linear predictor, "glm2" a GLM with a second order polynomial, and "gam" a generalised additive model (GAM) with smooting splines
  # threshold .. a numeric value indicating the absolute value of the correlation coefficient above which the paired correlation are judged as problematic
  # method .. a character string indicating which correlation coefficient (or covariance) is to be computed. One of "spearman" (default), "kendall", or "pearson"
  # sequence .. an optional character vector providing the order of importance of the predictors. This overrides the univar method
  # weights .. an optional vector of prior weights to be used in univariate GLMs or GAMs
  
  # The function is based on Dormann et al. (2013, https://doi.org/10.1111/j.1600-0587.2012.07348.x) and Zurell et al. (2020, https://doi.org/10.1111/jbi.13608).
  # It selects variables based on removing correlations > 0.7, retaining those
  # variables more important with respect to y
  # Order of importance can be provided by the character vector 'sequence'
  
  # 1. step: cor-matrix
  # 2. step: importance vector
  # 3. step: identify correlated pairs
  # 4. step: in order of importance: remove collinear less important variable,
  #           recalculate correlation matrix a.s.f.
  
  # Make k-fold data partitions
  if (length(kfold)==1) {
    ks <- dismo::kfold(y, k = kfold)
  } else {
    ks <- kfold
  }
  
  compute.univar.cv <- function(variable, response, family,univar,ks,weights){
    preds <- numeric(length(response))
    
    for (n in unique(ks)) {
      df <- data.frame(occ=response,env=variable)
      train_df <- df[!ks ==n,]
      test_df <-  df[ks==n, ]
      
      m1 <- switch(univar,
                   glm1 = glm(occ ~ env, data=train_df, family=family, weights=weights[!ks ==n]),
                   glm2 = glm(occ ~ poly(env,2), data=train_df, family=family, weights=weights[!ks ==n]),
                   gam = mgcv::gam(occ ~ s(env,k=4), data=train_df, family=family, weights=weights[!ks ==n]))
      
      preds[ks==n] <- predict(m1,newdata=test_df,type='response')
    }
    d2 <- expl_deviance(response,preds, weights=weights)
    ifelse(d2<0,0,d2)
  }
  
  imp <- apply(X, 2, compute.univar.cv, response=y, family=family, univar=univar, ks=ks,weights=weights)
  
  cm <- cor(X, method=method)
  
  if (is.null(sequence)) {
    sort.imp <- colnames(X)[order(imp,decreasing=T)]
  } else { 
    sort.imp <- sequence 
  }
  
  pairs <- which(abs(cm) >= threshold, arr.ind=T) # identifies correlated variable pairs
  index <- which(pairs[,1]==pairs[,2])           # removes entry on diagonal
  pairs <- pairs[-index,]                        # -"-
  
  exclude <- NULL
  for (i in 1:length(sort.imp))
  {
    if ((sort.imp[i] %in% row.names(pairs))&
        ((sort.imp[i] %in% exclude)==F)) {
      cv<-cm[setdiff(row.names(cm),exclude),sort.imp[i]]
      cv<-cv[setdiff(names(cv),sort.imp[1:i])]
      exclude<-c(exclude,names(which((abs(cv)>=threshold)))) }
  }
  
  pred_sel <- sort.imp[!(sort.imp %in% unique(exclude)),drop=F]
  return(list(D2=sort(imp, decreasing = T), cor_mat=cm, pred_sel=pred_sel))
}




#-------------------------------------------------------------------------------

# 3. Function for explained deviance -------------------------------------------

expl_deviance <- function(obs, pred, family='binomial',weights=rep(1, length(obs))){
  require(dismo)
  
  if (family=='binomial') {pred <- ifelse(pred<.00001,.00001,ifelse(pred>.9999,.9999,pred))}
  
  null_pred <- rep(mean(obs), length(obs))
  
  1 - (dismo::calc.deviance(obs, pred, family=family, weights=weights) / 
         dismo::calc.deviance(obs, null_pred, family=family, weights=weights))
}



#-------------------------------------------------------------------------------

# 4. Function for SDM evaluation -----------------------------------------------


#' evalSDM
#'
#' Evaluate SDM perfomance
#' @param observation vector containing the observed response 
#' @param predictions vector containing the predictions
#' @param thresh threshold to use for calculating threshold dependent performance measures. If NULL (the default) then threshold is optimised based on the observed presence/absence data provided and the thresh.method
#' @param thresh.method a string indicating which method to use for optimising the binarising threshold (see ?PresenceAbsence::optimal.thresholds. Defaults to "MaxSens+Spec" (the maximum of sensitivity+specificity). Will be ignored if thresh is provided.
#' @param req.sens additional argument to PresenceAbsence::optimal.thresholds(). Will be ignored if thresh is provided.
#' @param req.spec additional argument to PresenceAbsence::optimal.thresholds(). Will be ignored if thresh is provided.
#' @param FPC additional argument to PresenceAbsence::optimal.thresholds(). Will be ignored if thresh is provided.
#' @param FNC additional argument to PresenceAbsence::optimal.thresholds(). Will be ignored if thresh is provided.
#' @param weights an optional vector of prior weights used in the model
#' @return A dataframe with performance statistics.
#' @examples 
#' data(Anguilla_train)
#' m1 <- glm(Angaus ~ poly(SegSumT,2), data=Anguilla_train, family='binomial')
#' preds_cv <- crossvalSDM(m1, kfold=5, traindat=Anguilla_train, colname_species = 'Angaus', colname_pred = 'SegSumT')
#' evalSDM(Anguilla_train$Angaus, preds_cv)
#' @export
evalSDM <- function(observation, predictions, thresh=NULL, thresh.method='MaxSens+Spec', req.sens=0.85, req.spec = 0.85, FPC=1, FNC=1, weights=rep(1, length(observation))){
  thresh.dat <- data.frame(ID=seq_len(length(observation)), 
                           obs = observation,
                           pred = predictions)
  
  if (is.null(thresh)) {
    thresh.mat <- PresenceAbsence::optimal.thresholds(DATA= thresh.dat, req.sens=req.sens, req.spec = req.spec, FPC=FPC, FNC=FNC)
    thresh <- thresh.mat[thresh.mat$Method==thresh.method,2]
  }
  
  cmx.opt <- PresenceAbsence::cmx(DATA= thresh.dat, threshold=thresh)
  
  data.frame(AUC = PresenceAbsence::auc(thresh.dat, st.dev=F),
             TSS = TSS(cmx.opt), 
             Kappa = PresenceAbsence::Kappa(cmx.opt, st.dev=F),
             Sens = PresenceAbsence::sensitivity(cmx.opt, st.dev=F),
             Spec = PresenceAbsence::specificity(cmx.opt, st.dev=F),
             PCC = PresenceAbsence::pcc(cmx.opt, st.dev=F),
             D2 = expl_deviance(observation, predictions, weights=weights),
             thresh = thresh)
}




#-------------------------------------------------------------------------------

# 5. Function for TSS calculation ----------------------------------------------

#' TSS
#'
#' Calculates the true skill statistic (sensitivity+specificity-1) \insertCite{allouche2006}{mecofun}.
#' 
#' @importFrom Rdpack reprompt
#' 
#' @param cmx a confusion matrix
#' 
#' @return A numeric value.
#' 
#' @examples  TSS()
#' 
#' @references
#' \insertAllCited{}
#' 
#' @export
TSS = function(cmx){
  PresenceAbsence::sensitivity(cmx, st.dev=F) + 
    PresenceAbsence::specificity(cmx, st.dev=F) - 1
}





#-------------------------------------------------------------------------------

# 6. Function for Boyce index calculation using smoothing methods --------------

# Code from Liu et al (2024)

sfbi <- function(prd1, prd0, ktry=10) {
  p <- c(prd1, prd0)
  n1 <- length(prd1)
  n0 <- length(prd0)
  prd <- seq(min(p), max(p), length=n0)
  oc <- c(rep(1, n1), rep(0, n0))
  
  md_tp = mgcv::gam(oc ~ s(p,bs="tp",k=min(ktry,length(unique(p)))), family=binomial)
  prd_tp = predict(md_tp,newdata=data.frame(p=prd),type='response')
  md_cr = mgcv::gam(oc ~ s(p,bs="cr",k=min(ktry,length(unique(p)))), family=binomial)
  prd_cr = predict(md_cr,newdata=data.frame(p=prd),type='response')
  md_bs = mgcv::gam(oc ~ s(p,bs="bs",k=min(ktry,length(unique(p)))), family=binomial)
  prd_bs = predict(md_bs,newdata=data.frame(p=prd),type='response')
  md_ps = mgcv::gam(oc ~ s(p,bs="ps",k=min(ktry,length(unique(p)))), family=binomial)
  prd_ps = predict(md_ps,newdata=data.frame(p=prd),type='response')
  md_ad = mgcv::gam(oc ~ s(p, bs = "ad",k=min(ktry,length(unique(p)))), family=binomial)
  prd_ad = predict(md_ad,newdata=data.frame(p=prd),type='response')
  prd_m = (prd_tp + prd_cr + prd_bs + prd_ps + prd_ad)/5
  SBI_tp <- cor(prd,prd_tp,method="spearman")
  SBI_cr <- cor(prd,prd_cr,method="spearman")
  SBI_bs <- cor(prd,prd_bs,method="spearman")
  SBI_ps <- cor(prd,prd_ps,method="spearman")
  SBI_ad <- cor(prd,prd_ad,method="spearman")
  SBI_m <- cor(prd,prd_m,method="spearman")
  
  return(c(SBI_tp, SBI_cr, SBI_bs, SBI_ps, SBI_ad, SBI_m))
}



#-------------------------------------------------------------------------------

# 7. Function for predictions to new data to generate ensemble response plots --


#' predictSDM
#'
#' Make SDM predictions 
#' @param model model object
#' @param newdata a data frame in which to look for variables for which predictions should be made. 
#' @return A numeric vector with predictions.
#' @examples 
#' data(Anguilla_train)
#' data(Anguilla_test)
#' m1 <- glm(Angaus ~ poly(SegSumT,2), data=Anguilla_train, family='binomial')
#' predictSDM(m1, Anguilla_test)
#' @export
predictSDM <- function(model, newdata) {
  switch(class(model)[1],
         Bioclim = predict(model, newdata),
         Domain = predict(model, newdata),
         glm = predict(model, newdata, type='response'),
         Gam = predict(model, newdata, type='response'),
         gam = predict(model, newdata, type='response'),
         negbin = predict(model, newdata, type='response'),
         rpart = predict(model, newdata),
         randomForest.formula = switch(model$type,
                                       regression = predict(model, newdata, type='response'),
                                       classification = predict(model, newdata, type='prob')[,2]),
         randomForest = switch(model$type,
                               regression = predict(model, newdata, type='response'),
                               classification = predict(model, newdata, type='prob')[,2]),
         gbm = switch(ifelse(is.null(model$gbm.call),"GBM","GBM.STEP"), 
                      GBM.STEP =  predict.gbm(model, newdata, 
                                              n.trees=model$gbm.call$best.trees, type="response"),
                      GBM = predict.gbm(model, newdata, 
                                        n.trees=model$n.trees, type="response")),
         maxnet = predict(model, newdata, type="logistic"))
}
