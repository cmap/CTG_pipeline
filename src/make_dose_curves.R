source('fitDoseCurve.R', chdir = T)
source('plotDoseCurve.R', chdir = T)
library(plyr)

# requires a df with "dose", "brd_id" and "ccle_name" columns
make_dose_response_curves <- function(df) {
    sigmoids <- df %>%
        dplyr::mutate(pert_effective_dose = dose, grouping = 1) %>%
        plyr::ddply(.(brd_id, ccle_name), function(x) get_sigmoid(x, fitDoseCurve(x))) %>%
        dplyr::rename(viability = y, logdose = x) %>%
        dplyr::mutate(cl_compound = paste0(brd_id, '-', ccle_name)) 
        
    
    p <- df %>% 
        dplyr::mutate(cl_compound = paste0(brd_id, '-', ccle_name)) %>%
        ggplot() +
        geom_point(aes(log(dose), viability)) +
        geom_path(data = sigmoids, mapping = aes(logdose, viability)) + 
        scale_x_continuous(breaks = log(c(0.001, 0.01, 0.1, 1, 10)),
                           labels = -3:1) +
        facet_wrap('cl_compound') + 
        geom_hline(yintercept = c(0, 100), linetype = 3, size = 0.5) +
        scale_y_continuous(limits = c(-25, 200))
    
    return(p)
}

make_dose_response_curves_overlap_cell_lines <- function(df) {
    sigmoids <- df %>%
        dplyr::mutate(pert_effective_dose = dose, grouping = 1) %>%
        plyr::ddply(.(brd_id, ccle_name), function(x) get_sigmoid(x, fitDoseCurve(x))) %>%
        dplyr::rename(viability = y, logdose = x) %>%
        dplyr::mutate(cl_compound = paste0(brd_id, '-', ccle_name)) 
    
    
    p <- df %>% 
        dplyr::mutate(cl_compound = paste0(brd_id, '-', ccle_name)) %>%
        ggplot() +
        geom_point(aes(log(dose), viability, color = ccle_name)) +
        geom_path(data = sigmoids, mapping = aes(logdose, viability, color = ccle_name)) + 
        scale_x_continuous(breaks = log(c(0.001, 0.01, 0.1, 1, 10)),
                           labels = -3:1) +
        geom_hline(yintercept = c(0, 100), linetype = 3, size = 0.5) +
        scale_y_continuous(limits = c(-25, 200)) + 
        theme(legend.position="bottom",legend.direction="vertical")
    
    return(p)
}

make_dose_response_curves_overlap_cell_lines_collapse_replicates <- function(df) {
    sigmoids <- df %>%
        dplyr::mutate(pert_effective_dose = dose, grouping = 1) %>%
        plyr::ddply(.(brd_id, ccle_name), function(x) get_sigmoid(x, fitDoseCurve(x))) %>%
        dplyr::rename(viability = y, logdose = x) %>%
        dplyr::mutate(cl_compound = paste0(brd_id, '-', ccle_name)) 
    
    
    p <- df %>% 
        dplyr::mutate(cl_compound = paste0(brd_id, '-', ccle_name)) %>%
        dplyr::group_by(dose, ccle_name, brd_id) %>%
        dplyr::summarise(viability = median(viability)) %>%
        dplyr::ungroup() %>%
        ggplot() +
        geom_point(aes(log(dose), viability, color = ccle_name)) +
        geom_path(data = sigmoids, mapping = aes(logdose, viability, color = ccle_name)) + 
        scale_x_continuous(breaks = log(c(0.001, 0.01, 0.1, 1, 10)),
                           labels = -3:1) +
        geom_hline(yintercept = c(0, 100), linetype = 3, size = 0.5) +
        scale_y_continuous(limits = c(-25, 200)) + 
        theme(legend.position="bottom",legend.direction="vertical") + 
        ggtitle(df$brd_id[1])
    
    return(p)
}



get_curve_stats <- function(df) {
    df %>%
        dplyr::mutate(pert_effective_dose = dose, grouping = 1) %>%
        plyr::ddply(.(brd_id, ccle_name), fitDoseCurve)
}


make_dose_response_curves_for_CTG_qc <- function(main_data, main_data_name, 
                                                 background_data, background_name) {
    
    shared_compounds <- intersect(background_data$brd_id, main_data$brd_id)
    shared_cell_lines <- intersect(background_data$ccle_name, main_data$ccle_name)
    
    filtered_background_data <- background_data %>%
        dplyr::filter(brd_id %in% shared_compounds) %>%
        dplyr::filter(ccle_name %in% shared_cell_lines)
    
    filtered_main_data <- main_data %>%
        dplyr::filter(brd_id %in% shared_compounds) %>%
        dplyr::filter(ccle_name %in% shared_cell_lines)
    
    main_sigmoids <- filtered_main_data %>%
        dplyr::mutate(pert_effective_dose = as.numeric(dose), grouping = 1) %>%
        plyr::ddply(.(brd_id, ccle_name), function(x) get_sigmoid(x, fitDoseCurve(x))) %>%
        dplyr::rename(viability = y, logdose = x) %>%
        dplyr::mutate(cl_compound = paste0(brd_id, '-', ccle_name)) %>%
        dplyr::mutate(source = main_data_name)
    
    background_sigmoids <- filtered_background_data %>%
        dplyr::mutate(pert_effective_dose = dose, grouping = 1) %>%
        plyr::ddply(.(brd_id, ccle_name), function(x) get_sigmoid(x, fitDoseCurve(x))) %>%
        dplyr::rename(viability = y, logdose = x) %>%
        dplyr::mutate(cl_compound = paste0(brd_id, '-', ccle_name)) %>%
        dplyr::mutate(source = background_name)
    
    full_sigmoids <- rbind(main_sigmoids, background_sigmoids)
    
    shared_data_columns <- intersect(colnames(filtered_main_data), colnames(filtered_background_data))
    
    
    combined_data <- rbind(dplyr::select(filtered_main_data, shared_data_columns) %>% dplyr::mutate(source = main_data_name),
                           dplyr::select(filtered_background_data, shared_data_columns) %>% dplyr::mutate(source = background_name))
    
    p <- combined_data %>%
            ggplot() +
            geom_point(aes(log(dose), viability, color = source)) +
            geom_path(data = full_sigmoids, mapping = aes(logdose, viability, color = source)) + 
            scale_x_continuous(breaks = log(c(0.001, 0.01, 0.1, 1, 10)),
                               labels = -3:1) +
            facet_grid(brd_id ~ ccle_name) + 
            geom_hline(yintercept = c(0, 100), linetype = 3, size = 0.5) +
            scale_y_continuous(limits = c(-25, 200))
    
    return(p)
}

make_dose_response_curves_for_CTG_qc_single_cell_line<- function(main_data, main_data_name, 
                                                 background_data, background_name) {
    
    shared_compounds <- intersect(background_data$brd_id, main_data$brd_id)
    shared_cell_lines <- intersect(background_data$ccle_name, main_data$ccle_name)
    
    filtered_background_data <- background_data %>%
        dplyr::filter(brd_id %in% shared_compounds) %>%
        dplyr::filter(ccle_name %in% shared_cell_lines)
    
    filtered_main_data <- main_data %>%
        dplyr::filter(brd_id %in% shared_compounds) %>%
        dplyr::filter(ccle_name %in% shared_cell_lines)
    
    main_sigmoids <- filtered_main_data %>%
        dplyr::mutate(pert_effective_dose = as.numeric(dose), grouping = 1) %>%
        plyr::ddply(.(brd_id, ccle_name), function(x) get_sigmoid(x, fitDoseCurve(x))) %>%
        dplyr::rename(viability = y, logdose = x) %>%
        dplyr::mutate(cl_compound = paste0(brd_id, '-', ccle_name)) %>%
        dplyr::mutate(source = main_data_name)
    
    background_sigmoids <- filtered_background_data %>%
        dplyr::mutate(pert_effective_dose = dose, grouping = 1) %>%
        plyr::ddply(.(brd_id, ccle_name), function(x) get_sigmoid(x, fitDoseCurve(x))) %>%
        dplyr::rename(viability = y, logdose = x) %>%
        dplyr::mutate(cl_compound = paste0(brd_id, '-', ccle_name)) %>%
        dplyr::mutate(source = background_name)
    
    full_sigmoids <- rbind(main_sigmoids, background_sigmoids)
    
    shared_data_columns <- intersect(colnames(filtered_main_data), colnames(filtered_background_data))
    
    
    combined_data <- rbind(dplyr::select(filtered_main_data, shared_data_columns) %>% dplyr::mutate(source = main_data_name),
                           dplyr::select(filtered_background_data, shared_data_columns) %>% dplyr::mutate(source = background_name))
    
    p <- combined_data %>%
        ggplot() +
        geom_point(aes(log(dose), viability, color = source)) +
        geom_path(data = full_sigmoids, mapping = aes(logdose, viability, color = source)) + 
        scale_x_continuous(breaks = log(c(0.001, 0.01, 0.1, 1, 10)),
                           labels = -3:1) +
        facet_wrap('brd_id') + 
        geom_hline(yintercept = c(0, 100), linetype = 3, size = 0.5) +
        scale_y_continuous(limits = c(-25, 200))
    
    return(p)
}
