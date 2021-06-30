# Script to process CTG data for PRISM
# Andrew Boghossian
# Requires 5 arguments to run (in the following order)
# 1) path to treatment meta (mergedfile)
# 2) path to plate meta (mapping)
# 3) path to raw data (CTG...)
# 4) path to directory to store output (a single .Rds file)
# 5) a string name for the project

# LOAD FUNCTIONS AND LIBRARIES ----
suppressMessages(source("src/helperFunctions.R"))
suppressMessages(source("src/make_dose_curves.R", chdir = T))

# ARGUMENTS ----
script_args <- commandArgs(trailingOnly = TRUE)
if (length(script_args) != 5) {
  stop("Please supply necessary arguments", call. = FALSE)
}

meta_treat_path <- script_args[1]
meta_plate_path <- script_args[2]
data_path <- script_args[3]
output_dir <- script_args[4]
project <- script_args[5]

if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = T)}
# LOAD DATA ----
print("Loading the meta data")
# read in treatment meta (mergedfile)
treatment_meta <- data.table::fread(meta_treat_path, data.table = F) %>%
  dplyr::rename(broad_id = broad_sample,
                dose = mmoles_per_liter,
                Mapping = Color,
                Well_Position = "Dest Well",
                replicate = Replicate) %>%
  dplyr::mutate(pert_type = dplyr::case_when(
    str_detect(broad_id, fixed("BRD-K88510285")) & dose > 9.5 ~ "poscon",
    broad_id == "" | broad_id == "DMSO" ~ "negcon",
    T ~ "treatment"),
    Well_Position = paste0(str_sub(Well_Position, 1, 1),
                           str_pad(str_sub(Well_Position, 2, -1), width = 2, pad = "0"))) %>%
  dplyr::mutate(source = project, Replicate = toString(Replicate),
                plate_map_name = str_replace_all(plate_map_name, fixed("_"), "-"))

if ("broad_sample_2" %in% colnames(treatment_meta)) {
  treatment_meta %<>%
    dplyr::rename(broad_id_2 = broad_sample_2, dose_2 = mmoles_per_liter_2)
}

# read in plate meta (mapping)
plate_meta <- data.table::fread(meta_plate_path, data.table = F) %>%
  dplyr::rename(plate_map_name = PLATE_MAP_NAME,
                ccle_name = `Cell Line`,
                replicate = Replicate,
                arp_barcode = `Assay Plate Barcode`,) %>%
  dplyr::mutate(plate_map_name = str_replace_all(plate_map_name, fixed("_"), "-"),
                Replicate = toString(Replicate))
if(is.integer(plate_meta$replicate)) {
  plate_meta$replicate <- paste0("X", plate_meta$replicate)
}

# join treatment and plate meta
combined_meta <- dplyr::left_join(plate_meta, treatment_meta,
                                  by = c("plate_map_name", "Mapping", "replicate"))

print("Reading in raw data")
# read in raw data (CTG...)
raw_data <- read_enspire(data_path) %>%
  dplyr::mutate(Barcode = str_replace_all(Barcode, "[^[:alnum:]]", "")) %>%
  dplyr::rename(arp_barcode = Barcode,
                Well_Position = Well) %>%
  dplyr::full_join(combined_meta, by = c("arp_barcode", "Well_Position"))

# test to make sure additional entries not created
if (nrow(raw_data) != nrow(combined_meta)) {
  stop("Something is weird with the meta mapping. The number of rows in the
       combined meta should be the same as the number of data points")
}

# CALCULATE VIABILITY and NORMALIZE ----
print("Calculating viabilities")
normalized_data <- raw_data %>%
  dplyr::group_by(arp_barcode, Mapping) %>%
  dplyr::mutate(viability = 100 * (lum - median(lum[pert_type == 'poscon'])) /
                  (median(lum[pert_type == 'negcon']) - median(lum[pert_type == 'poscon'])),
                log_fc = log2(lum / median(lum[pert_type == 'negcon']))) %>%
  dplyr::ungroup()

# if unable to calculate viability (no poscons usually)
if (all(is.na(normalized_data$viability))) {
  print("Unable to calculate normalized viability, calculating 2^lfc instead")
  normalized_data$viability <- 2^normalized_data$log_fc
}

# WRITE ----
readr::write_rds(normalized_data,
                 paste0(output_dir, "/", project, "_ctg_viability_data.Rds"))
readr::write_csv(normalized_data,
                 paste0(output_dir, "/", project, "_ctg_viability_data.csv"))
print(paste("Results written to:",
            paste0(output_dir, "/", project, "_ctg_viability_data.Rds")))
