# script to process CTG data for PRISM
suppressMessages(source("src/helperFunctions.R"))
suppressMessages(source("src/make_dose_curves.R", chdir = T))

# for commandline running
script_args <- commandArgs(trailingOnly = TRUE)
if (length(script_args) != 1) {
  script_args <- c(readline(prompt = "Please supply path to data directory: "))
}

base_dir <- script_args[1]
project <- word(basename(base_dir), 1, -2, sep = fixed("_"))

#---- Load the data ----
# read in treatment meta
meta_treat_path <- paste0(base_dir, "/mergedfile.txt")
treatment_meta <- data.table::fread(meta_treat_path, data.table = F) %>%
  dplyr::rename(color = Color, 
                replicate = Replicate,
                broad_id = broad_sample,
                dose = mmoles_per_liter) %>%
  dplyr::mutate(perturbation_type = dplyr::case_when(
    str_detect(broad_id, fixed("BRD-K88510285")) & dose > 10 ~ "poscon",
    broad_id == "" ~ "negcon",
    T ~ "treatment")) %>%
  dplyr::mutate(source = "project") %>%
  dplyr::rename(pert_type = perturbation_type)

# read in plate meta
meta_plate_path <- paste(base_dir, "mapping.csv")
plate_meta <- data.table::fread(meta_plate_path, data.table = F) %>%
  dplyr::rename(plate_map_name = PLATE_MAP_NAME, 
                color = Mapping, 
                ccle_name = `Cell Line`,
                replicate = Replicate,
                arp_barcode = `Assay Plate Barcode`) %>%
  dplyr::mutate(plate_map_name = str_replace_all(plate_map_name, fixed("_"), "-"))

combined_meta <- dplyr::left_join(plate_meta, treatment_meta)

# read in raw data
ctg_data_file_paths <- list.files(paste0(base_dir, "/data"))
raw_data <- ctg_data_file_paths %>%
  purrr::set_names(str_sub(., 1, 10)) %>%
  purrr::map_dfr(~read_enspire(file.path(base_dir, "data", .x)),
                 .id = "arp_barcode") %>%
  dplyr::mutate(`Dest Well` = paste0(row, col)) %>%
  dplyr::full_join(combined_meta)

if (nrow(raw_data) != nrow(combined_meta)) {
  stop("Something is weird with the meta mapping. The number of rows in the
       combined meta should be the same as the number of data points")
}

#---- Normalize and write output ----
# normalize
normalized_data <- raw_data %>%
  dplyr::group_by(arp_barcode, color) %>%
  dplyr::mutate(viability = 100 * (lum - median(lum[pert_type == 'poscon'])) / 
                  (median(lum[pert_type == 'negcon']) - median(lum[pert_type == 'poscon'])),
                log_fc = log2(lum / median(lum[pert_type == 'negcon']))) %>%
  dplyr::ungroup()

readr::write_rds(normalized_data,
                 paste0(base_dir, "/", project, "ctg_viability_data.Rds"))