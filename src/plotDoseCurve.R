sigmoid <- function(x, high = 1, slope = 1, ec50 = 0, low = 0) {
  a <- high - low
  scaler <- slope* (x - log(ec50))
  a/(1 + exp(scaler)) + low
}

plot_sigmoid <- function(data_dose_cpd, data_curved_cpd){
  
  #data_dose_cpd %<>% group_by(grouping) 
  #              %<>% summarize(dose_range = seq(min(log(pert_effective_dose)), log(25), length.out = 100))
  
  dose_range <- seq(min(log(data_dose_cpd$pert_effective_dose)), log(25), length.out = 100)
  #df_sigmoid1 <- group_by(data_dose_cpd, grouping) %<>%
  #               summarize(x = dose_range y = )
  df_sigmoid <- data.frame(x = dose_range, 
                           y = sigmoid(dose_range, 
                                       high = data_curved_cpd$upper_asym,
                                       low = data_curved_cpd$lower_asym,
                                       slope = data_curved_cpd$slope,
                                       ec50 = data_curved_cpd$ec50))
  df_mean_se <- data_dose_cpd %>% 
    dplyr::group_by(logdose) %>% 
    dplyr::summarise(mean = mean(viability), 
              se = sd(viability)/sqrt(n()),
              dose = mean(pert_effective_dose)) 
  
  # title <- paste0('erlotinib in ', data_curved_cpd$ccle_name,
  #                 '\n AUC = ', round(data_curved_cpd$auc, 2))
  
  
  ggplot(df_mean_se, aes(log(dose), mean)) + 
    geom_point() +
    geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0.2) + 
    geom_path(data = df_sigmoid, aes(x, y), color = 'darkred') +
    geom_hline(yintercept = c(0, 100), linetype = 3, size = 0.5) +
    scale_x_continuous(breaks = log(c(0.001, 0.01, 0.1, 1, 10)),
                       labels = c(0.001, 0.01, 0.1, 1, 10)) +
    theme_bw() + xlab('Dose (uM)') + ylab('Viability') +
    #ggtitle(title) + 
    theme(plot.title = element_text(size = 12))
}

get_sigmoid <- function(data_dose_cpd, data_curved_cpd) {
  dose_range <- seq(min(log(data_dose_cpd$pert_effective_dose)), log(25), length.out = 100)
  df_sigmoid <- data.frame(x = dose_range, 
                           y = sigmoid(dose_range, 
                                       high = data_curved_cpd$upper_asym,
                                       low = data_curved_cpd$lower_asym,
                                       slope = data_curved_cpd$slope,
                                       ec50 = data_curved_cpd$ec50))
}

get_points <- function(data_dose_cpd){
  df_mean_se <- data_dose_cpd %>% 
    dplyr::group_by(logdose) %>% 
    dplyr::summarise(mean = mean(viability), 
              se = sd(viability)/sqrt(n()),
              dose = mean(pert_effective_dose)) 
}

plot_sigmoid_part <- function(data_dose_cpd, data_curved_cpd){
  
  #data_dose_cpd %<>% group_by(grouping) 
  #              %<>% summarize(dose_range = seq(min(log(pert_effective_dose)), log(25), length.out = 100))
  
  
  dose_range <- seq(min(log(data_dose_cpd$pert_effective_dose)), log(25), length.out = 300)
  #df_sigmoid1 <- group_by(data_dose_cpd, grouping) %<>%
  #               summarize(x = dose_range y = )
  df_sigmoid <- data.frame(x = dose_range, 
                           y = sigmoid(dose_range, 
                                       high = data_curved_cpd$upper_asym,
                                       low = data_curved_cpd$lower_asym,
                                       slope = data_curved_cpd$slope,
                                       ec50 = data_curved_cpd$ec50))
  df_mean_se <- data_dose_cpd %>% 
    group_by(logdose) %>% 
    summarise(mean = mean(viability), 
              se = sd(viability)/sqrt(n()),
              dose = mean(pert_effective_dose)) 
}