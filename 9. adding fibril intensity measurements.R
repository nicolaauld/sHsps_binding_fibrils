# =============================================================================
# Compile all data for each day for the extra fibril analysis
# =============================================================================

library(tidyverse)
library(tidylog)

# -----------------------------------------------------------------------------
# 0. Paths
# -----------------------------------------------------------------------------

# images and ImageJ output csvs live on the network drive ("vorg")
vorg_folder <- "Z:/Chaperone_subgroup/NicolaA/Fibrils with sHsps/Testing code/"

experiment_date <- basename(sub("/$", "", vorg_folder))

# combined csvs get written here
local_folder <- paste0("C:/Nicola/UNIVERSITY/PhD/SM/Fibrils and sHsps/", experiment_date, "/full_data/")
dir.create(local_folder, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Choose which treatments to include
# -----------------------------------------------------------------------------

all_items       <- list.files(vorg_folder, full.names = TRUE)
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
# 2. Fibril intensity data
# -----------------------------------------------------------------------------

all_intensity_list <- list()

for (b in seq_along(treatments)) {
  current_treat_folder <- paste0(vorg_folder, treatments[b], "/")
  current_treatment    <- treatments[b]
  print(current_treat_folder)
  
  all_items             <- list.files(current_treat_folder, full.names = TRUE)
  directories           <- all_items[file.info(all_items)$isdir]
  treat_directory_names <- basename(directories)
  print(treat_directory_names)
  
  fibril_files_list <- list()
  for (c in seq_along(treat_directory_names)) {
    data_path <- paste0(current_treat_folder, treat_directory_names[c], "/488 fibrils/Extra_fibril_analysis/")
    data_files <- list.files(path = data_path, pattern = "\\.csv$", full.names = TRUE)
    # keep only the intensity files - exclude the width ones
    selected_files <- data_files[!grepl("_width", data_files)]
    print(selected_files)
    
    fibril_files_list[[c]] <- selected_files |>
      rlang::set_names() |>
      purrr::map(readr::read_delim) |>
      purrr::list_rbind(names_to = "file")
  }
  fibril_data_files_all <- bind_rows(fibril_files_list)
  print(unique(fibril_data_files_all$file))
  
  # clean up file name and add tracking columns: treatment, experiment date,
  # image number, analysis method
  fibril_data_files_all <- fibril_data_files_all %>%
    mutate(file = sub(".*Channel488", "Channel488", file),
           just_image_number = str_extract(file, "#\\d+"),
           experiment_date = experiment_date,
           Treatment_name = current_treatment,
           analysis_method = "extra")
  
  all_intensity_list[[b]] <- fibril_data_files_all
}

all_intensity_df <- bind_rows(all_intensity_list) %>%
  mutate(Image_treatment = paste0(just_image_number, "-", Treatment_name),
         Protein = str_extract(Treatment_name, "(?<=_).*"))

print(unique(all_intensity_df$Image_treatment)) # confirm every expected file loaded before saving

write.csv(all_intensity_df,
          file = paste0(local_folder, experiment_date, "_all_fibril_intensity_data.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 3. Adding intensity measurements to the data
# -----------------------------------------------------------------------------

all_fibrils_normal <- read.csv(paste0("C:/Nicola/UNIVERSITY/PhD/SM/Fibrils and sHsps/", experiment_date, "/full_data/", experiment_date, "_all_fibrils_normal.csv"))

# This contains the fibrils in the same order that the intensity measurements were taken in 
# we are going to use that to match the intensity data back to the fibrils

individual_fibrils_in_order <- all_fibrils_normal %>%
  reframe(mean_x = mean(X),
          mean_y = mean(Y),
          .by = c(Treatment_name, just_image_number, Contour.ID, Length))

# 1. same number of fibrils?
nrow(individual_fibrils_in_order) == nrow(all_intensity_df)

# 2. does every row sit in the same treatment and image in both dataframes?
all(individual_fibrils_in_order$Treatment_name == all_intensity_df$Treatment_name &
      individual_fibrils_in_order$just_image_number == all_intensity_df$just_image_number)

# 3. if both are TRUE, attach by position
individual_fibrils_in_order <- individual_fibrils_in_order %>%
  mutate(mean_intensity = all_intensity_df$Mean)

# -----------------------------------------------------------------------------
# 4. Saving the intensity data
# -----------------------------------------------------------------------------

write.csv(all_intensity_df,
          file = paste0("C:/Nicola/UNIVERSITY/PhD/SM/Fibrils and sHsps/", experiment_date, "/analysis_outputs/", "_individual_fibril_info_intensity_data.csv"),
          row.names = FALSE)
