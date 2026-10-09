# =============================================================================
# End vs middle sHsp foci positions on fibrils (2024 edition)
# =============================================================================

library(tidyverse)
library(tidylog)

# -----------------------------------------------------------------------------
# 0. Paths + constants
# -----------------------------------------------------------------------------

experiment_date   <- "20260203"
parent_folder_C   <- paste0("C:/Nicola/UNIVERSITY/PhD/SM/Fibrils and sHsps/", experiment_date, "/")
analysis_folder   <- paste0(parent_folder_C, "analysis_outputs/")
full_data_folder  <- paste0(parent_folder_C, "full_data/")

# minimum fibril length (pixels) to include in this analysis
min_fibril_length_px <- 12

# how many pixels from a fibril tip count as the "end" region (also used as the
# margin for flagging fibrils that cross the image edge)
number_of_pixels_per_end <- 3

px_to_um <- 0.16

desired_order <- c("50nM_aBc3D", "100nM_aBc3D", "250nM_aBc3D", "500nM_aBc3D", "1uM_aBc3D")
custom_labels <- str_wrap(c("50 nM aBc3D", "100 nM aBc3D", "250 nM aBc3D", "500 nM aBc3D", "1 uM aBc3D"), 10)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------

all_Hsps_coloc_info            <- read.csv(paste0(analysis_folder, "all_Hsps_only_w_coloc_info.csv"))
individual_fibril_info         <- read.csv(paste0(analysis_folder, "individual_fibril_info.csv"))
truly_all_fibrils_clean_normal <- read.csv(paste0(analysis_folder, "truly_all_fibrils_clean_normal.csv"))

# -----------------------------------------------------------------------------
# 2. Filter to long fibrils only and build shared IDs
# -----------------------------------------------------------------------------

truly_all_fibrils_clean_normal_long <- truly_all_fibrils_clean_normal %>%
  filter(Length > min_fibril_length_px)

# unique IDs so we can match Hsp foci to fibrils without looping per-treatment
all_Hsps_coloc_info <- all_Hsps_coloc_info %>%
  mutate(unique_hsp_id = paste(experiment_date, just_image_number, Treatment_name, X_and_Y, sep = "-"))

truly_all_fibrils_clean_normal_long <- truly_all_fibrils_clean_normal_long %>%
  mutate(unique_hsp_id = paste(experiment_date, Image_treatment, X_and_Y, sep = "-"))

# keep only the Hsps colocalised with a long fibril
XY_long <- truly_all_fibrils_clean_normal_long$unique_hsp_id
all_Hsps_coloc_info_long <- all_Hsps_coloc_info %>% filter(unique_hsp_id %in% XY_long)

all_Hsps_coloc_info_long <- all_Hsps_coloc_info_long %>%
  mutate(unique_fibril_id = paste(experiment_date, just_image_number, Treatment_name, coloc_Fibril_ID, sep = "-"))

truly_all_fibrils_clean_normal_long <- truly_all_fibrils_clean_normal_long %>%
  mutate(unique_fibril_id = paste(experiment_date, Image_treatment, Fibril_ID, sep = "-"))

# -----------------------------------------------------------------------------
# 3. Get fibril end-point coordinates for each colocalised Hsp
# -----------------------------------------------------------------------------

all_Hsps_coloc_info_long <- all_Hsps_coloc_info_long %>%
  mutate(fibril_EndX1 = NA, fibril_EndX2 = NA, fibril_EndY1 = NA, fibril_EndY2 = NA)

for (i in seq_len(nrow(all_Hsps_coloc_info_long))) {
  coloc_fibril_ID <- all_Hsps_coloc_info_long$unique_fibril_id[i]
  this_fibril <- truly_all_fibrils_clean_normal_long %>% filter(unique_fibril_id == coloc_fibril_ID)
  
  j <- nrow(this_fibril) # last row = other end of the fibril
  all_Hsps_coloc_info_long$fibril_EndX1[i] <- this_fibril$X[1]
  all_Hsps_coloc_info_long$fibril_EndX2[i] <- this_fibril$X[j]
  all_Hsps_coloc_info_long$fibril_EndY1[i] <- this_fibril$Y[1]
  all_Hsps_coloc_info_long$fibril_EndY2[i] <- this_fibril$Y[j]
}

# -----------------------------------------------------------------------------
# 4. Classify each Hsp focus as END or MIDDLE
# -----------------------------------------------------------------------------

all_Hsps_coloc_info_long <- all_Hsps_coloc_info_long %>%
  mutate(dist_1 = sqrt((fibril_EndX1 - X)^2 + (fibril_EndY1 - Y)^2),
         dist_2 = sqrt((fibril_EndX2 - X)^2 + (fibril_EndY2 - Y)^2),
         radius = number_of_pixels_per_end,
         Where  = NA)

for (k in seq_len(nrow(all_Hsps_coloc_info_long))) {
  if (all_Hsps_coloc_info_long$dist_1[k] < all_Hsps_coloc_info_long$radius[k] ||
      all_Hsps_coloc_info_long$dist_2[k] < all_Hsps_coloc_info_long$radius[k]) {
    all_Hsps_coloc_info_long$Where[k] <- "END"
  } else {
    all_Hsps_coloc_info_long$Where[k] <- "MIDDLE"
  }
}

# -----------------------------------------------------------------------------
# 5. Flag fibrils crossing the edge of the image
# -----------------------------------------------------------------------------
# Uses the raw (unfiltered) full_data file, since no fibrils have been removed
# from it - so it gives the true max/min x and y coordinates per image.

all_fibrils_raw <- read.csv(paste0(full_data_folder, experiment_date, "_all_fibrils_normal.csv"))

edge_coords <- all_fibrils_raw %>%
  group_by(Treatment_name, just_image_number) %>%
  reframe(min_x = min(X), min_y = min(Y), max_x = max(X), max_y = max(Y))

edge_crossing_fibrils <- all_Hsps_coloc_info_long %>%
  filter(fibril_EndX1 > (max(edge_coords$max_x) - number_of_pixels_per_end) |
           fibril_EndX2 > (max(edge_coords$max_x) - number_of_pixels_per_end) |
           fibril_EndX1 < number_of_pixels_per_end | fibril_EndX2 < number_of_pixels_per_end |
           fibril_EndY1 > (max(edge_coords$max_y) - number_of_pixels_per_end) |
           fibril_EndY2 > (max(edge_coords$max_y) - number_of_pixels_per_end) |
           fibril_EndY1 < number_of_pixels_per_end | fibril_EndY2 < number_of_pixels_per_end)

edge_crossing_IDs <- edge_crossing_fibrils$unique_fibril_id

all_Hsps_coloc_info_long <- all_Hsps_coloc_info_long %>%
  mutate(edge_crossing = case_when(unique_fibril_id %in% edge_crossing_IDs ~ "Y", TRUE ~ "N"))

write.csv(all_Hsps_coloc_info_long,
          paste0(analysis_folder, "all_Hsps_coloc_info_12pixelmin_ENDMIDDLEinfo_edge_info.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. Calculate END/MIDDLE on the per-fibril level
# -----------------------------------------------------------------------------

individual_fibril_info_long <- individual_fibril_info %>%
  filter(fibril_length > min_fibril_length_px) %>%
  mutate(experiment_date = experiment_date,
         Image_treatment = paste(just_image_number, Treatment_name, sep = "-"),
         unique_fibril_id = paste(experiment_date, just_image_number, Treatment_name, Fibril_ID, sep = "-"),
         number_of_end_foci = NA,
         number_of_middle_foci = NA,
         end_foci_names = NA,
         middle_foci_names = NA,
         edge_crossing = NA)

# for each fibril, count and name its END and MIDDLE colocalised foci
for (l in seq_len(nrow(individual_fibril_info_long))) {
  current_fibril <- individual_fibril_info_long$unique_fibril_id[l]
  
  end_foci <- all_Hsps_coloc_info_long %>% filter(unique_fibril_id == current_fibril & Where == "END")
  individual_fibril_info_long$number_of_end_foci[l] <- nrow(end_foci)
  
  middle_foci <- all_Hsps_coloc_info_long %>% filter(unique_fibril_id == current_fibril & Where == "MIDDLE")
  individual_fibril_info_long$number_of_middle_foci[l] <- nrow(middle_foci)
  
  if (individual_fibril_info_long$number_of_end_foci[l] > 0) {
    individual_fibril_info_long$end_foci_names[l] <- paste(end_foci$X_and_Y, collapse = ",")
  }
  if (individual_fibril_info_long$number_of_middle_foci[l] > 0) {
    individual_fibril_info_long$middle_foci_names[l] <- paste(middle_foci$X_and_Y, collapse = ",")
  }
}

individual_fibril_info_long <- individual_fibril_info_long %>%
  mutate(end_region_length_pixels = 2 * number_of_pixels_per_end,
         end_region_length_um = end_region_length_pixels * px_to_um,
         middle_region_length_pixels = fibril_length - end_region_length_pixels,
         middle_region_length_um = fibril_length_um - end_region_length_um,
         hsps_per_end_length = number_of_end_foci / end_region_length_um,
         hsps_per_middle_length = number_of_middle_foci / middle_region_length_um)

# NOTE: edge_crossing_IDs only covers colocalised fibrils - any fibril with no
# colocalised sHsps that happens to touch the edge won't be flagged here.
individual_fibril_info_long <- individual_fibril_info_long %>%
  mutate(edge_crossing_COLOC_ONLY = case_when(unique_fibril_id %in% edge_crossing_IDs ~ "Y", TRUE ~ "N"))

fibrils_removed_for_edge_crossing <- individual_fibril_info_long %>%
  filter(edge_crossing_COLOC_ONLY == "Y")

write.csv(individual_fibril_info_long,
          paste0(analysis_folder, "individual_fibril_info_12pixelmin_ENDMIDDLEinfo_edge_info.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. Summary tables (edge-crossing fibrils excluded)
# -----------------------------------------------------------------------------

individual_fibril_info_long <- individual_fibril_info_long %>%
  filter(edge_crossing_COLOC_ONLY == "N")

end_middle_summary_by_treatment <- individual_fibril_info_long %>%
  group_by(Treatment_name) %>%
  summarise(
    mean_number_of_end_foci = mean(number_of_end_foci, na.rm = TRUE),
    mean_number_of_middle_foci = mean(number_of_middle_foci, na.rm = TRUE),
    mean_hsps_per_end_length = mean(hsps_per_end_length, na.rm = TRUE),
    mean_hsps_per_middle_length = mean(hsps_per_middle_length, na.rm = TRUE),
    .groups = "drop"
  )

individual_fibril_info_long_non_zero <- individual_fibril_info_long %>%
  filter(how_many_coloc_hsp_on_this_fibril > 0)

standard_error <- function(x) {
  sd(x, na.rm = TRUE) / sqrt(length(x))
}

end_middle_summary_by_treatment_non_zero <- individual_fibril_info_long_non_zero %>%
  group_by(Treatment_name) %>%
  summarise(
    mean_hsps_per_end_length = mean(hsps_per_end_length, na.rm = TRUE),
    sem_hsps_per_end_length = standard_error(hsps_per_end_length),
    mean_hsps_per_middle_length = mean(hsps_per_middle_length, na.rm = TRUE),
    sem_hsps_per_middle_length = standard_error(hsps_per_middle_length),
    .groups = "drop"
  )

# reshape means and standard errors to long format, then join back together
long_means <- end_middle_summary_by_treatment_non_zero %>%
  pivot_longer(cols = c(mean_hsps_per_end_length, mean_hsps_per_middle_length),
               names_to = "position_on_fibril", values_to = "mean_value") %>%
  mutate(position_on_fibril = ifelse(position_on_fibril == "mean_hsps_per_end_length", "End", "Middle")) %>%
  select(Treatment_name, position_on_fibril, mean_value)

long_sems <- end_middle_summary_by_treatment_non_zero %>%
  pivot_longer(cols = c(sem_hsps_per_end_length, sem_hsps_per_middle_length),
               names_to = "position_on_fibril", values_to = "sem_value") %>%
  mutate(position_on_fibril = ifelse(position_on_fibril == "sem_hsps_per_end_length", "End", "Middle")) %>%
  select(Treatment_name, position_on_fibril, sem_value)

combined_summary <- left_join(long_means, long_sems, by = c("Treatment_name", "position_on_fibril")) %>%
  mutate(Treatment_name = factor(Treatment_name, levels = desired_order)) %>%
  arrange(Treatment_name)

# -----------------------------------------------------------------------------
# 8. Plot
# -----------------------------------------------------------------------------

column_plot <- ggplot(combined_summary, aes(x = Treatment_name, y = mean_value, fill = position_on_fibril)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_value - sem_value, ymax = mean_value + sem_value),
                position = position_dodge(width = 0.9), width = 0.25) +
  labs(x = "Treatment Name", y = "Mean Hsp foci per Length", fill = "Foci Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = custom_labels)
print(column_plot)

ggsave(paste0(analysis_folder, experiment_date, "_fibril_positions_of_foci.png"),
       plot = column_plot, units = "px", width = 2500, height = 1600, dpi = 300)

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------

rm(list = ls())
dev.off()
cat("\014")
