library(tidyverse)
library(magrittr)

# function to read data from EnSpire output to long table
read_enspire <- function(data_path) {
  # read in raw file (fill and skip to allow for format)
  raw_data <- read.csv(data_path, fill = T, header = T,
                       skip = 1, blank.lines.skip = T, check.names = F)
  raw_data <- raw_data[apply(raw_data, 1, function(x) !all(is.na(x) | x=="")), ]
  
  # get row where meta and data are separated
  split_row <- which(raw_data$Plate=="No background information available.")
  
  # get meta
  meta_data <- raw_data[1:(split_row-1), 1:4]
  
  # get data
  n_row <- max(as.integer(meta_data$Plate)) * 384 + (max(as.integer(meta_data$Plate))-2)
  ctg_data <- raw_data[(split_row+2):(split_row+2+n_row), 1:3]
  colnames(ctg_data) <- raw_data[(split_row+1), 1:3]
  
  # join meta and data, convert readout to integer
  output_data <- dplyr::left_join(meta_data, ctg_data,
                                  by=c("Plate"="PlateNumber")) %>%
    dplyr::rename(lum = `MeasA:Result`)
  output_data$lum<- as.integer(output_data$lum)
  
  return(output_data)
}

# calculate robust SSMDs
calculate_ssmd_robust <- function(pop1, pop2) {
  return(median(pop1, na.rm = T) - median(pop2, na.rm=T)) /
    sqrt(mad(pop1, na.rm=T)^2 + mad(pop2, na.rm=T)^2)
}

# map broad IDs to structure IDs
brd_id_to_structure_id <- function(brd_ids) {
  id_mapping <- data.frame(brd_ids = unique(brd_ids)) %>%
    dplyr::mutate(structure_id = purrr::map_chr(brd_ids, ~paste(str_split(.x, fixed('-'))[[1]][1:2], collapse = '-')))
  
  data.frame(brd_ids = brd_ids) %>%
    dplyr::left_join(id_mapping) %>%
    dplyr::select(structure_id) %>%
    unlist() %>%
    return()
}
