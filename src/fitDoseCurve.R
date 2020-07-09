#' @title Fit dose-response curves to a data frame
#' 
#' @description 
#' \code{CurvePlate} fits sigmoidal curves to each compound/cell_line pair in the input plate data.
#' Assumes a 4-parameter fit with the following parameterization:
#' f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}
#' where (e) is the inflection point of the curve, (b) is the slope, (c) is the lower asymptote, 
#' and (d) is the upper asymptote
#' 
#' @return A data.frame with fields for the IC50 ('ic50'), potency ('slope'), 
#' auc, and R2 of curve-fit for each compound/cell_line pair


# helper functions -----------------------------------------------
# for calculating auc
compute_AUC <- function(max_logdose, min_logdose, slope, ec50, lower_asym, upper_asym){
  integra_max <- integral.logistic(max_logdose, slope, log(ec50), lower_asym, upper_asym)
  integra_min <- integral.logistic(min_logdose, slope, log(ec50), lower_asym, upper_asym)
  
  (integra_max - integra_min)/(upper_asym * (max_logdose - min_logdose))
}
  
integral.logistic <- function(log_x, slope, log_ec50, lower_asym, upper_asym){
  # the analytical solution for integral of sigmoid function
  # integral of c/(1 + exp(-a*(x + b))) = (c/a)*log(1 + exp(a*(x + b)))
  
  # for the dose curve f(x) = lower_asym + (upper_asym-lower_asym) / (1 + exp(slope * (log_dose(x)-log_ec50)))
  # the integral at dose x is
  # lower_asym*log_dose(x) + (upper_asym-lower_asym)/(-slope)* log(1 + exp((-slope) * (log_dose(x)-log_ec50)))

  
  # if lower asymptote is negative, set it to be zero
  lower_asym <- max(lower_asym, 0)
  # if lower asymptote is greater than upper asymptote, switch A, B
  if(upper_asym < lower_asym){
    temp <- upper_asym
    upper_asym <- lower_asym
    lower_asym <- temp
    slope <- -slope # are you sure?
  }
  # prevent infinite values from arising
  exponent <- (-slope) * (log_x -log_ec50)
  if(exponent > 700){
    (upper_asym - lower_asym) * (log_x -log_ec50) + lower_asym*log_x
  } else {
    (upper_asym - lower_asym)/(-slope) * log(1 + exp(exponent)) + lower_asym*log_x
  }
}

# for calculating ic50
calculate.ic50 <- function(slope, ec50, B, A){
  ic50 <- exp((log((((A-B)/((A/2) - B)) - 1))/slope + log(ec50)))
  ic50
}

# parse fitted model
parse_mod <- function(mod, upper_asym, lower_asym){
  if(is.null(mod)){
    data.frame(converged=FALSE, r2=NA, robust_r2=NA,
               ec50=NA, slope=NA, lower_asym=NA, upper_asym=NA, 
               ic50=NA, auc=NA)
  }else{
    mod_params <- mod$fit$par %>% set_names(mod$parNames[[2]])

    # extract convergence
    converged <- mod$fit$convergence
    
    # extract R2
    pred <- predict(mod, newdata=data.frame(mod$dataList$dose))
    true <- mod$dataList$origResp
    r2 <- 1 - (sum((true - pred)^2) / sum((true - mean(true))^2))
    rr2 <- 1 - ((mad(true - pred)^2) / mad(true)^2)
    
    # extract ec50
    ec50 <- mod_params[['e']]
    # extract slope
    slope <- mod_params[['b']]
    
    # extract lower asymptote
    if(is.na(lower_asym)){
      lower_asym <- mod_params[['c']]
    }
    
    # extract upper asymptote
    if(is.na(upper_asym)){
      upper_asym <- mod_params$d
    }
    # calculate auc
    auc <- compute_AUC(max(log(mod$dataList$dose)), min(log(mod$dataList$dose)),
                       slope, ec50, lower_asym, upper_asym)
    # calculate ic50
    ic50 <- calculate.ic50(slope, ec50, lower_asym, upper_asym)
    
    data.frame(converged=converged, 
               r2=r2, robust_r2=rr2, ec50=ec50, slope=slope, lower_asym=lower_asym,
               upper_asym=upper_asym, ic50=ic50, auc=auc)
  }

}

# main functions -----------------------------------------------
fitDoseCurve <- function(curve.data, upper_asym=100, lower_asym=NA, fit.top5 = F){
  
  #' @title Fit dose-response curves to each PRISM plate
  #' 
  #' @description 
  #' \code{CurvePlate} fits sigmoidal curves to each compound/cell_line pair in the input plate data.
  #' Assumes a 4-parameter fit with the following parameterization:
  #' f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}
  #' where (e) is the inflection point of the curve, (b) is the slope, (c) is the lower asymptote, 
  #' and (d) is the upper asymptote
  #' 
  #' @return A data.frame with fields for the IC50 ('ic50'), potency ('slope'), 
  if(fit.top5) curve.data <- filter(curve.data, pert_effective_dose > 0.03)
  mod <- NULL
  try(suppressWarnings(mod <- drc::drm(viability ~ pert_effective_dose, 
                                       curveid = grouping,
                                       data=curve.data, 
                                       fct=drc::LL.4(fixed=c(NA, lower_asym, upper_asym, NA)),
                                       robust="median")), silent = T)
  
  if(class(mod) == "try-error") mod <- NULL
  res <- parse_mod(mod, upper_asym=upper_asym, lower_asym=lower_asym)
  res
}


