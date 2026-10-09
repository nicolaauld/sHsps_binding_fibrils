# =============================================================================
# Rename each trajectory column to "<fibril contour ID>_<X>_<Y>" (or 0 if not
# colocalised with a fibril), e.g. "Mean_56_19813.58_24.5_63.5". The first
# column is renamed to blank, since it's the (0,0) coordinate and doesn't
# correspond to a molecule.
# =============================================================================

library(tidyverse)
library(tidylog)

# -----------------------------------------------------------------------------
# 0. Paths
# -----------------------------------------------------------------------------

wdpath <- "C:/Nicola/UNIVERSITY/PhD/SM/Fibrils and sHsps/20260203/"
setwd(wdpath)

# -----------------------------------------------------------------------------
# 1. Load colocalisation info and build per-chaperone unique IDs
# -----------------------------------------------------------------------------

all_Hsps_coloc_info <- read.csv("analysis_outputs/all_Hsps_only_w_coloc_info.csv")

all_Hsps_coloc_info <- all_Hsps_coloc_info %>%
  mutate(coloc_Fibril_ContourID = replace_na(coloc_Fibril_ContourID, 0),
         Lauren_ID = paste(coloc_Fibril_ContourID, X, Y, sep = "_"),
         Image_treatment_type = paste(just_image_number, Treatment_name, analysis_method, sep = "-")) %>%
  filter(analysis_method == "normal")

vector_of_unique_image_treatments <- unique(all_Hsps_coloc_info$Image_treatment_type)
print(vector_of_unique_image_treatments)

# -----------------------------------------------------------------------------
# 2. Parse each Image_treatment_type into image number / concentration / protein
# -----------------------------------------------------------------------------

pattern <- "^#(\\d+)-([^-_]+(?:-[^-_]+)*)_([^_\\-]+)-.*$"

image_number  <- vector()
concentration <- vector()
protein_name  <- vector()

for (entry in vector_of_unique_image_treatments) {
  match <- regmatches(entry, regexec(pattern, entry))
  if (length(match[[1]]) > 1) {
    image_number  <- c(image_number, match[[1]][2])
    concentration <- c(concentration, match[[1]][3])
    protein_name  <- c(protein_name, match[[1]][4])
  }
}

extracted_info <- data.frame(image_number = image_number,
                             concentration = concentration,
                             protein_name = protein_name,
                             original_entry = vector_of_unique_image_treatments)

print(extracted_info)

# -----------------------------------------------------------------------------
# 3. Rename trajectory columns for each treatment/image
# -----------------------------------------------------------------------------

for (i in seq_len(nrow(extracted_info))) {
  
  current_path <- paste0(wdpath, "imagejresults/", extracted_info$protein_name[i], "/", extracted_info$concentration[i], "/")
  
  all_Hsps_current_treatment <- all_Hsps_coloc_info %>%
    filter(Image_treatment_type == extracted_info$original_entry[i])
  
  # --- Colocalised trajectories ---------------------------------------------
  
  current_path_coloc <- paste0(current_path, "coloc/")
  csv_name <- paste0(extracted_info$protein_name[i], "_", extracted_info$concentration[i], "_coloc", extracted_info$image_number[i], ".csv")
  full_file_path <- paste0(current_path_coloc, csv_name)
  print(full_file_path)
  
  coloc_trajectories <- read.csv(full_file_path)
  coloc_Hsps <- all_Hsps_current_treatment %>% filter(coloc_w_fibril == "Yes")
  
  original_names <- colnames(coloc_trajectories)
  extra_info <- as.data.frame(coloc_Hsps$Lauren_ID)
  
  new_names <- original_names
  new_names[1] <- ""
  new_names[2] <- ""
  for (j in 3:length(original_names)) {
    new_names[j] <- paste0(original_names[j], "_", extra_info[j - 2, 1])
  }
  colnames(coloc_trajectories) <- new_names
  print(colnames(coloc_trajectories))
  
  write.csv(coloc_trajectories, file = full_file_path, row.names = FALSE)
  
  # --- Non-colocalised trajectories ------------------------------------------
  # These files can end up containing ALL trajectories (not just the true
  # non-coloc ones), so after renaming we also strip out any columns that
  # actually match a colocalised Hsp ID.
  
  current_path_noncoloc <- paste0(current_path, "non-coloc/")
  csv_name_noncoloc <- paste0(extracted_info$protein_name[i], "_", extracted_info$concentration[i], "_non-coloc", extracted_info$image_number[i], ".csv")
  full_file_path_noncoloc <- paste0(current_path_noncoloc, csv_name_noncoloc)
  print(full_file_path_noncoloc)
  
  noncoloc_trajectories <- read.csv(full_file_path_noncoloc)
  noncoloc_Hsps <- all_Hsps_current_treatment %>% filter(coloc_w_fibril == "No")
  
  original_names <- colnames(noncoloc_trajectories)
  extra_info <- as.data.frame(noncoloc_Hsps$Lauren_ID)
  
  new_names <- original_names
  new_names[1] <- ""
  new_names[2] <- ""
  if (length(original_names) > 2) {
    for (j in 3:length(original_names)) {
      new_names[j] <- paste0(original_names[j], "_", extra_info[j - 2, 1])
    }
  }
  colnames(noncoloc_trajectories) <- new_names
  print(colnames(noncoloc_trajectories))
  
  # remove any noncoloc trajectory columns that actually match a colocalised Hsp ID
  coloc_Hsps_LRID <- coloc_Hsps$Lauren_ID
  
  columns_to_remove <- sapply(colnames(noncoloc_trajectories), function(col) {
    any(sapply(coloc_Hsps_LRID, function(id) grepl(id, col, fixed = TRUE)))
  })
  
  actual_noncoloc_trajectories <- noncoloc_trajectories[, !columns_to_remove]
  colnames(actual_noncoloc_trajectories)[2] <- ""
  
  write.csv(actual_noncoloc_trajectories, file = full_file_path_noncoloc, row.names = FALSE)
}

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------

rm(list = ls())
dev.off()
cat("\014")
