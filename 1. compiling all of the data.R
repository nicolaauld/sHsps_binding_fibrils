# =============================================================================
# Compile all data gathered on a single experiment day (all treatments) into
# combined CSVs for downstream analysis.
# =============================================================================

library(tidyverse)
library(tidylog)

# -----------------------------------------------------------------------------
# 0. Paths
# -----------------------------------------------------------------------------

# images and ImageJ output csvs live on the network drive ("vorg")
parent_folder_vorg <- "Z:/Chaperone_subgroup/NicolaA/Fibrils with sHsps/Testing code/"

# combined csvs get written here
parent_folder_local <- "C:/Nicola/UNIVERSITY/PhD/SM/Fibrils and sHsps/Testing code/full_data/"

# creates this local folder if it doesn't exist already
if (!dir.exists(parent_folder_local)) {
  dir.create(parent_folder_local, recursive = TRUE)
}

experiment_date <- basename(sub("/$", "", parent_folder_vorg))

dir.create(parent_folder_local, showWarnings = FALSE, recursive = TRUE)

file_name_normal_fibrils    <- paste0(parent_folder_local, experiment_date, "_all_fibrils_normal.csv")
file_name_scrambled_fibrils <- paste0(parent_folder_local, experiment_date, "_all_fibrils_scrambled.csv")
file_name_hsps              <- paste0(parent_folder_local, experiment_date, "_all_Hsps.csv")

# -----------------------------------------------------------------------------
# 1. Choose which treatments to include
# -----------------------------------------------------------------------------

all_items       <- list.files(parent_folder_vorg, full.names = TRUE)
directories     <- all_items[file.info(all_items)$isdir]
directory_names <- basename(directories)
print(directory_names)

treatments <- rep(NA, length(directory_names))

for (a in seq_along(directory_names)) {
  current_treatment_name <- directory_names[a]
  console_prompt <- paste("Is", current_treatment_name, "a treatment you want to include? Y or N (case sensitive)  ")
  include <- readline(console_prompt)
  if (include == "Y") {
    treatments[a] <- current_treatment_name
  }
}

treatments <- na.omit(treatments)
print(treatments)

# -----------------------------------------------------------------------------
# 2. Normal fibril data
# -----------------------------------------------------------------------------

all_fibrils_list <- list()

for (b in seq_along(treatments)) {
  current_treat_folder <- paste0(parent_folder_vorg, treatments[b], "/")
  current_treatment    <- treatments[b]
  print(current_treat_folder)
  
  all_items           <- list.files(current_treat_folder, full.names = TRUE)
  directories         <- all_items[file.info(all_items)$isdir]
  treat_directory_names <- basename(directories)
  print(treat_directory_names)
  
  fibril_files_list <- list()
  for (c in seq_along(treat_directory_names)) {
    data_path <- paste0(current_treat_folder, treat_directory_names[c], "/AF647 photobleach/Colocalisation_analysis/")
    data_files <- list.files(path = data_path, pattern = "\\.csv$", full.names = TRUE)
    # exclude the "_colocalized" files - those are just the colocalised chaperones,
    # we want the ones titled with the image name containing the fibril info
    selected_files <- data_files[!grepl("_colocalized", data_files)]
    print(selected_files)
    
    fibril_files_list[[c]] <- selected_files |>
      rlang::set_names() |>
      purrr::map(readr::read_delim) |>
      purrr::list_rbind(names_to = "file")
  }
  fibril_data_files_all <- bind_rows(fibril_files_list)
  print(unique(fibril_data_files_all$file))
  
  # clean up file name, add tracking columns: treatment, experiment date,
  # image number, analysis method, Fibril_ID
  fibril_data_files_all <- fibril_data_files_all %>%
    mutate(file = sub(".*Channel488", "Channel488", file),
           just_image_number = str_extract(file, "#\\d+"),
           experiment_date = experiment_date,
           Fibril_ID = `Contour ID` * Length,
           Treatment_name = current_treatment,
           analysis_method = "normal")
  
  all_fibrils_list[[b]] <- fibril_data_files_all
}

all_fibrils_df <- bind_rows(all_fibrils_list) %>%
  mutate(Image_treatment = paste0(just_image_number, "-", Treatment_name),
         Protein = str_extract(Treatment_name, "(?<=_).*"))

print(unique(all_fibrils_df$Image_treatment)) # confirm every expected file loaded before saving

write.csv(all_fibrils_df, file = file_name_normal_fibrils, row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Scrambled fibril data
# -----------------------------------------------------------------------------

all_fibrils_list <- list()

for (b in seq_along(treatments)) {
  current_treat_folder <- paste0(parent_folder_vorg, treatments[b], "/")
  current_treatment    <- treatments[b]
  
  all_items           <- list.files(current_treat_folder, full.names = TRUE)
  directories         <- all_items[file.info(all_items)$isdir]
  treat_directory_names <- basename(directories)
  print(treat_directory_names)
  
  fibril_files_list <- list()
  for (c in seq_along(treat_directory_names)) {
    data_path <- paste0(current_treat_folder, treat_directory_names[c], "/AF647 photobleach/Scrambled_fibril_colocalisation_analysis/")
    selected_files <- list.files(path = data_path, pattern = "\\Fibril_scramble_Hsp_results.csv$", full.names = TRUE)
    print(selected_files)
    
    data_files_all <- selected_files |>
      rlang::set_names() |>
      purrr::map(readr::read_delim) |>
      purrr::list_rbind(names_to = "file")
    
    if (ncol(data_files_all) == 12) {
      data_files_all <- data_files_all %>% select(-12)
    }
    
    fibril_files_list[[c]] <- data_files_all
  }
  fibril_data_files_all <- bind_rows(fibril_files_list)
  
  fibril_data_files_all <- fibril_data_files_all %>%
    mutate(file = sub(".*Channel488", "Channel488", file),
           just_image_number = str_extract(file, "#\\d+"),
           experiment_date = experiment_date,
           Fibril_ID = `Contour ID` * Length,
           Treatment_name = current_treatment,
           analysis_method = "scrambled")
  
  all_fibrils_list[[b]] <- fibril_data_files_all
}

all_fibrils_df_scrambled <- bind_rows(all_fibrils_list) %>%
  mutate(Image_treatment = paste0(just_image_number, "-", Treatment_name),
         Protein = str_extract(Treatment_name, "(?<=_).*"))

print(unique(all_fibrils_df_scrambled$Image_treatment))

write.csv(all_fibrils_df_scrambled, file = file_name_scrambled_fibrils, row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Hsp data (normal + scrambled analysis methods)
# -----------------------------------------------------------------------------

all_Hsps_list <- list()

for (b in seq_along(treatments)) {
  current_treat_folder <- paste0(parent_folder_vorg, treatments[b], "/")
  current_treatment    <- treatments[b]
  print(current_treat_folder)
  
  all_items           <- list.files(current_treat_folder, full.names = TRUE)
  directories         <- all_items[file.info(all_items)$isdir]
  treat_directory_names <- basename(directories)
  print(treat_directory_names)
  
  Hsp_files_list_normal   <- list()
  Hsp_files_list_scramble <- list()
  
  for (c in seq_along(treat_directory_names)) {
    data_path <- paste0(current_treat_folder, treat_directory_names[c], "/AF647 photobleach/Scrambled_fibril_colocalisation_analysis/")
    
    data_files <- list.files(path = data_path, pattern = "\\Hsp_results.csv$", full.names = TRUE)
    normal_analysed_files <- data_files[!grepl("_Fibril_scramble", data_files)]
    
    Hsp_files_list_normal[[c]] <- normal_analysed_files |>
      rlang::set_names() |>
      purrr::map(readr::read_delim) |>
      purrr::list_rbind(names_to = "file")
    
    scramble_analysed_files <- list.files(path = data_path, pattern = "\\Hsp_scrambled_results.csv$", full.names = TRUE)
    
    Hsp_files_list_scramble[[c]] <- scramble_analysed_files |>
      rlang::set_names() |>
      purrr::map(readr::read_delim) |>
      purrr::list_rbind(names_to = "file")
  }
  
  normal_Hsp_files_all   <- bind_rows(Hsp_files_list_normal)
  scramble_Hsp_files_all <- bind_rows(Hsp_files_list_scramble)
  
  normal_Hsp_files_all$analysis_method   <- "normal"
  scramble_Hsp_files_all$analysis_method <- "scrambled"
  
  full_Hsp_files <- bind_rows(normal_Hsp_files_all, scramble_Hsp_files_all) %>%
    mutate(file = sub(".*Channel488", "Channel488", file),
           just_image_number = str_extract(file, "#\\d+"),
           experiment_date = experiment_date,
           Treatment_name = current_treatment)
  
  all_Hsps_list[[b]] <- full_Hsp_files
}

all_Hsps_df <- bind_rows(all_Hsps_list) %>%
  mutate(Image_treatment_type = paste0(just_image_number, "-", Treatment_name, "-", analysis_method),
         Protein = str_extract(Treatment_name, "(?<=_).*"))

print(unique(all_Hsps_df$Image_treatment_type))

write.csv(all_Hsps_df, file = file_name_hsps, row.names = FALSE)

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------

rm(list = ls())
#dev.off()
cat("\014")
