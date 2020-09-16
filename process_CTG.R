# script to process CTG data for PRISM
suppressMessages(source("src/helperFunctions.R"))
suppressMessages(source("src/make_dose_curves.R", chdir = T))

# for commandline running
script_args <- commandArgs(trailingOnly = TRUE)
if (length(script_args) != 3) {
  script_args <- c(readline(prompt = "Please supply path to data directory: "))
}

base_dir <- script_args[1]
output_dir <- script_args[2]
project <- script_args[3]

#---- Load the data ----
# read in treatment meta
meta_treat_path <- list.files(base_dir, pattern = "mergedfile*", full.names = T)
treatment_meta <- data.table::fread(meta_treat_path, data.table = F) %>%
  dplyr::rename(broad_id = broad_sample,
                dose = mmoles_per_liter,
                Mapping = Color,
                Well_Position = "Dest Well",
                replicate = Replicate) %>%
  dplyr::mutate(perturbation_type = dplyr::case_when(
    str_detect(broad_id, fixed("BRD-K88510285")) & dose > 10 ~ "poscon",
    broad_id == "" ~ "negcon",
    T ~ "treatment"),
    Well_Position = paste0(str_sub(Well_Position, 1, 1),
                           str_pad(str_sub(Well_Position, 2, -1), width = 2, pad = "0"))) %>%
  dplyr::mutate(source = "project") %>%
  dplyr::rename(pert_type = perturbation_type)

# read in plate meta
meta_plate_path <- list.files(base_dir, pattern = "mapping*", full.names = T)
plate_meta <- data.table::fread(meta_plate_path, data.table = F) %>%
  dplyr::rename(plate_map_name = PLATE_MAP_NAME,
                ccle_name = `Cell Line`,
                replicate = Replicate,
                arp_barcode = `Assay Plate Barcode`,) %>%
  dplyr::mutate(plate_map_name = str_replace_all(plate_map_name, fixed("_"), "-"))

combined_meta <- dplyr::left_join(plate_meta, treatment_meta,
                                  by = c("plate_map_name", "Mapping", "replicate"))

# read in raw data
data_path <- list.files(base_dir, pattern = "*CTG_raw_data*", full.names = T)
raw_data <- read_enspire(data_path) %>%
  dplyr::mutate(Barcode = str_sub(Barcode, 2, -1)) %>%
  dplyr::rename(arp_barcode = Barcode,
                Well_Position = Well) %>%
  dplyr::full_join(combined_meta, by = c("arp_barcode", "Well_Position"))

if (nrow(raw_data) != nrow(combined_meta)) {
  stop("Something is weird with the meta mapping. The number of rows in the
       combined meta should be the same as the number of data points")
}

#---- Normalize and write output ----
# normalize
normalized_data <- raw_data %>%
  dplyr::group_by(arp_barcode, Mapping) %>%
  dplyr::mutate(viability = 100 * (lum - median(lum[pert_type == 'poscon'])) /
                  (median(lum[pert_type == 'negcon']) - median(lum[pert_type == 'poscon'])),
                log_fc = log2(lum / median(lum[pert_type == 'negcon']))) %>%
  dplyr::ungroup()

readr::write_rds(normalized_data,
                 paste0(output_dir, "/", project, "_ctg_viability_data.Rds"))
