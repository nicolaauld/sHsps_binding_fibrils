# =============================================================================
# Relate py4bleaching results (max fluorescence, molecule counts) back to the
# fibril colocalisation / end-middle analysis already done.
# =============================================================================

library(tidyverse)
library(tidylog)

# -----------------------------------------------------------------------------
# 0. Paths
# -----------------------------------------------------------------------------

experiment_date <- "20260203"
parent_folder   <- paste0("C:/Nicola/UNIVERSITY/PhD/SM/Fibrils and sHsps/", experiment_date, "/")
setwd(parent_folder)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------

molecule_counts      <- read.csv("python_results/calculate_molecule_size/molecule_counts.csv")
all_Hsps_coloc_info  <- read.csv("analysis_outputs/all_Hsps_only_w_coloc_info.csv")

# =============================================================================
# Coloc
# =============================================================================

coloc_Hsps_only <- all_Hsps_coloc_info %>% filter(coloc_w_fibril == "Yes" & analysis_method == "normal")

# building matching IDs
molecule_counts <- molecule_counts %>%
  mutate(Lauren_ID = paste(Contour_ID, coordsX, coordsY, sep = "_"))

coloc_Hsps_only <- coloc_Hsps_only %>%
  mutate(Lauren_ID = paste(coloc_Fibril_ContourID, X, Y, sep = "_"),
         max_fluorescence = NA,
         mol_count_all_small = NA,
         mol_count_last_step = NA,
         mol_count_single_step = NA,
         conc = NA)

# for every colocalised Hsp, pull in its py4bleaching results if it has one.
# some Hsps will have been filtered out during py4bleaching's own cleaning
# step, so those will stay NA here.
for (i in seq_len(nrow(coloc_Hsps_only))) {
  
  current_molecule_ID <- coloc_Hsps_only$Lauren_ID[i]
  current_molecule_in_py4bleaching_result <- molecule_counts %>% filter(Lauren_ID == current_molecule_ID)
  
  if (nrow(current_molecule_in_py4bleaching_result) > 0) {
    coloc_Hsps_only$max_fluorescence[i]      <- unique(current_molecule_in_py4bleaching_result$max_fluorescence)
    coloc_Hsps_only$mol_count_all_small[i]   <- unique(current_molecule_in_py4bleaching_result$all_small_mol_count)
    coloc_Hsps_only$mol_count_last_step[i]   <- unique(current_molecule_in_py4bleaching_result$last_step_mol_count)
    coloc_Hsps_only$mol_count_single_step[i] <- unique(current_molecule_in_py4bleaching_result$single_step_mol_count)
    coloc_Hsps_only$conc[i]                  <- unique(current_molecule_in_py4bleaching_result$protein1) # concentration ended up stored in the "protein1" column
  }
}

# -----------------------------------------------------------------------------
# 2. Push the same info onto the 12-pixel-min END/MIDDLE dataframe
# -----------------------------------------------------------------------------

all_Hsps_coloc_12pixels_ENDMIDDLE <- read.csv("analysis_outputs/all_Hsps_coloc_info_12pixelmin_ENDMIDDLEinfo_edge_info.csv")

coloc_Hsps_only <- coloc_Hsps_only %>%
  mutate(Image_treatment = paste(just_image_number, Treatment_name, sep = "-"),
         unique_hsp_id = paste(experiment_date, Image_treatment, X_and_Y, sep = "-"))

all_Hsps_coloc_12pixels_ENDMIDDLE <- all_Hsps_coloc_12pixels_ENDMIDDLE %>%
  mutate(max_fluorescence = NA,
         mol_count_all_small = NA,
         mol_count_last_step = NA,
         mol_count_single_step = NA,
         conc = NA)

for (j in seq_len(nrow(all_Hsps_coloc_12pixels_ENDMIDDLE))) {
  
  current_molecule_ID <- all_Hsps_coloc_12pixels_ENDMIDDLE$unique_hsp_id[j]
  current_molecule_in_coloc_Hsps_result <- coloc_Hsps_only %>% filter(unique_hsp_id == current_molecule_ID)
  
  if (nrow(current_molecule_in_coloc_Hsps_result) > 0) {
    all_Hsps_coloc_12pixels_ENDMIDDLE$max_fluorescence[j]      <- unique(current_molecule_in_coloc_Hsps_result$max_fluorescence)
    all_Hsps_coloc_12pixels_ENDMIDDLE$mol_count_all_small[j]   <- unique(current_molecule_in_coloc_Hsps_result$mol_count_all_small)
    all_Hsps_coloc_12pixels_ENDMIDDLE$mol_count_last_step[j]   <- unique(current_molecule_in_coloc_Hsps_result$mol_count_last_step)
    all_Hsps_coloc_12pixels_ENDMIDDLE$mol_count_single_step[j] <- unique(current_molecule_in_coloc_Hsps_result$mol_count_single_step)
    all_Hsps_coloc_12pixels_ENDMIDDLE$conc[j]                  <- unique(current_molecule_in_coloc_Hsps_result$conc)
  }
}

write.csv(coloc_Hsps_only, file = "analysis_outputs/coloc_Hsps_only_w_trajectory_data.csv", row.names = FALSE)
write.csv(all_Hsps_coloc_12pixels_ENDMIDDLE, file = "analysis_outputs/all_Hsps_coloc_12pixels_ENDMIDDLE_w_trajectory_data.csv", row.names = FALSE)

# =============================================================================
# Non-coloc
# =============================================================================
# Note: a lot more non-coloc molecules seem to get filtered out during
# py4bleaching's cleaning step than coloc ones.

molecule_counts <- read.csv("python_results/calculate_molecule_size/molecule_counts.csv") %>%
  filter(variable1 == "Non-coloc") %>%
  mutate(Lauren_ID = paste(Contour_ID, coordsX, coordsY, protein1, sep = "_"))

non_coloc_Hsps_only <- all_Hsps_coloc_info %>% filter(coloc_w_fibril == "No" & analysis_method == "normal")

non_coloc_Hsps_only <- non_coloc_Hsps_only %>%
  mutate(conc = purrr::map_chr(strsplit(Treatment_name, "_"), 1),
         Lauren_ID = paste("0", X, Y, conc, sep = "_"),
         max_fluorescence = NA,
         mol_count_all_small = NA,
         mol_count_last_step = NA,
         mol_count_single_step = NA)

# NOTE: the ID Lauren's approach generates isn't fully unique here, even with
# concentration added - molecule_counts doesn't retain enough per-image info
# to disambiguate further, so some collisions are unavoidable with this
# matching approach.
for (i in seq_len(nrow(non_coloc_Hsps_only))) {
  
  current_molecule_ID <- non_coloc_Hsps_only$Lauren_ID[i]
  current_molecule_in_py4bleaching_result <- molecule_counts %>% filter(Lauren_ID == current_molecule_ID)
  
  if (nrow(current_molecule_in_py4bleaching_result) > 0) {
    non_coloc_Hsps_only$max_fluorescence[i]      <- unique(current_molecule_in_py4bleaching_result$max_fluorescence)
    non_coloc_Hsps_only$mol_count_all_small[i]   <- unique(current_molecule_in_py4bleaching_result$all_small_mol_count)
    non_coloc_Hsps_only$mol_count_last_step[i]   <- unique(current_molecule_in_py4bleaching_result$last_step_mol_count)
    non_coloc_Hsps_only$mol_count_single_step[i] <- unique(current_molecule_in_py4bleaching_result$single_step_mol_count)
  }
}

write.csv(non_coloc_Hsps_only, file = "analysis_outputs/non_coloc_Hsps_only_w_trajectory_data.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------

rm(list = ls())
dev.off()
cat("\014")
