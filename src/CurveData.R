CurveData <- function(data, upper_asym=100){
  
  #' assumes a 4-parameter sigmoidal curve fit with the following parameterization:
  #' f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}
  #' where (e) is the inflection point of the curve, (b) is the slope, (c) is the lower asymptote, 
  #' and (d) is the upper asymptote
  
  require(plyr)
  require(drc)
  
  cat("\nFitting dose response...")
  models <- data %>%
              plyr::dlply(.(ccle_name, compound), function(x){
                
                model <- NULL
                try(model <- drc::drm(viability ~ dose, 
                                      data=x, 
                                      fct=drc::LL.4(fixed=c(NA, NA, upper_asym, NA)),
                                      robust="median"))
                if(class(model) == "try-error"){
                  model <- NULL
                }
                  
                return(list(mod=model, dat=x))
                
              })
  
  # for calculating auc
  integral.logistic <- function(x, scal, xmid, B, A){
    # if lower asymptote is negative, set it to be zero
    B <- max(B, 0)
    # if lower asymptote is greater than upper asymptote, switch A, B
    if(A < B){
      temp <- A
      A <- B
      B <- temp
      scal <- -scal
    }
    # prevent infinite values from arising
    if((x-xmid)/scal > 700){
      return((A-B)*(x-xmid) + B*x)
    } else {
      return(scal*(A-B)*log(exp((x-xmid)/scal)+1) + B*x)
    }
  }
  
  # for calculating ic50
  calculate.ic50 <- function(slope, ec50, B, A){
    ic50 <- exp((log((((A-B)/((A/2) - B)) - 1))/slope + log(ec50)))
  }
  
  cat("\nExtracting model parameters...")
  model_parameters <- lapply(names(models), function(curve){
    
    mod <- models[[curve]][["mod"]]
    dat <- models[[curve]][["dat"]]
    
    if(is.null(mod)){
      return(data.frame(curve_id=curve, converged=FALSE, r2=NA, robust_r2=NA,
                        ec50=NA, ic50=NA, slope=NA, lower_asym=NA, upper_asym=upper_asym, 
                        auc=median(dat$viability, na.rm=T)/upper_asym))
    }
    
    mod_params <- data.frame(mod$parmMat) %>%
                    t() %>%
                    as.data.frame() %>%
                    magrittr::set_colnames(mod$parNames[[2]])
    
    # extract convergence
    converged <- mod$fit$convergence
    # extract R2 (non-robust)
    pred <- predict(mod, newdata=data.frame(mod$dataList$dose))
    true <- mod$dataList$origResp
    r2 <- 1 - (sum((true - pred)^2) / sum((true - mean(true))^2))
    rr2 <- 1 - ((mad(true - pred)^2) / mad(true)^2)
    # extract ec50
    ec50 <- mod_params$e
    # extract slope
    slope <- mod_params$b
    # extract lower asymptote
    lower_asym <- mod_params$c
    # extract upper asymptote
    if(is.na(upper_asym)){
      upper_asym <- mod_params$d
    }
    # calculate auc
    auc <- integral.logistic(max(log(mod$dataList$dose)), -1/slope, log(ec50), lower_asym, upper_asym) %>%
      magrittr::subtract(integral.logistic(min(log(mod$dataList$dose)), -1/slope, log(ec50), lower_asym, upper_asym)) %>%
      magrittr::divide_by(upper_asym*(max(log(mod$dataList$dose)) - min(log(mod$dataList$dose))))
    # calculate ic50
    ic50 <- calculate.ic50(slope, ec50, lower_asym, upper_asym)
    
    
    return(list(curve_id=curve, converged=converged, 
                      r2=r2, robust_r2=rr2, ec50=ec50, ic50=ic50, slope=slope, lower_asym=lower_asym,
                      upper_asym=upper_asym, auc=auc))
    
  }) %>%
    data.table::rbindlist()
  
  model_parameters %<>% tidyr::separate(curve_id, into=c("ccle_name", "compound"),
                                        sep="[.]")
  
  return(model_parameters)
  
}