# =============================================================================
# Testing for specific chaperone-fibril binding using the scrambled-data method
# =============================================================================

library(tidyverse)
library(tidylog)
library(ggpubr)

# -----------------------------------------------------------------------------
# 0. Paths
# -----------------------------------------------------------------------------

experiment_date   <- "Testing code"
parent_folder_C   <- paste0("C:/Nicola/UNIVERSITY/PhD/SM/Fibrils and sHsps/", experiment_date, "/")
full_data_folder  <- paste0(parent_folder_C, "full_data/")
coloc_out_folder  <- paste0(parent_folder_C, "analysis_outputs/")
binding_folder    <- paste0(parent_folder_C, "testing_for_specific_binding/")

# pixel -> micron conversion factor for fibril length
px_to_um <- 0.16

# size exclusion threshold (this is the size of "dirt"/debris, not real fibrils)
mean_seed_size_um <- 1.207
stdev_seed_um     <- 0.648
seed_cutoff_um    <- mean_seed_size_um + stdev_seed_um

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------

all_Hsps_only         <- read.csv(paste0(full_data_folder, experiment_date, "_all_Hsps.csv"))
all_fibrils_normal    <- read.csv(paste0(full_data_folder, experiment_date, "_all_fibrils_normal.csv"))
all_fibrils_scrambled <- read.csv(paste0(full_data_folder, experiment_date, "_all_fibrils_scrambled.csv"))

# -----------------------------------------------------------------------------
# 2. Build coordinate keys + combine fibril data
# -----------------------------------------------------------------------------

# X/Y coordinate key lets us check whether an Hsp focus falls on a fibril's coordinates
all_Hsps_only <- all_Hsps_only %>%
  mutate(X_and_Y = paste0(X, "&", Y))

all_fibrils_normal <- all_fibrils_normal %>%
  mutate(X_and_Y   = paste0(X2, "&", Y2),
         Length_um = Length * px_to_um)

all_fibrils_scrambled <- all_fibrils_scrambled %>%
  mutate(X_and_Y   = paste0(X2, "&", Y2),
         Length_um = Length * px_to_um)

truly_all_fibrils <- bind_rows(all_fibrils_normal, all_fibrils_scrambled) %>%
  mutate(Image_treatment_type = paste0(just_image_number, "-", Treatment_name, "-", analysis_method))

print(unique(truly_all_fibrils$Image_treatment_type))

# "clean" dataset: excludes anything at/below seed size, i.e. genuine fibrils only
all_fibrils_normal_clean    <- all_fibrils_normal    %>% filter(Length_um > seed_cutoff_um)
all_fibrils_scrambled_clean <- all_fibrils_scrambled %>% filter(Length_um > seed_cutoff_um)

truly_all_fibrils_clean <- bind_rows(all_fibrils_normal_clean, all_fibrils_scrambled_clean) %>%
  mutate(Image_treatment_type = paste0(just_image_number, "-", Treatment_name, "-", analysis_method))

# -----------------------------------------------------------------------------
# 3. Colocalisation matching: does each Hsp focus sit on a fibril?
# -----------------------------------------------------------------------------

all_Hsps_only <- all_Hsps_only %>%
  mutate(coloc_w_fibril      = NA_character_,
         coloc_Fibril_ID     = NA_character_,
         coloc_Fibril_ContourID = NA_character_)

for (i in seq_len(nrow(all_Hsps_only))) {
  
  current_treat_image_type <- all_Hsps_only$Image_treatment_type[i]
  current_XY_of_Hsp        <- all_Hsps_only$X_and_Y[i]
  
  only_these_fibrils <- truly_all_fibrils %>%
    filter(Image_treatment_type == current_treat_image_type)
  
  if (current_XY_of_Hsp %in% only_these_fibrils$X_and_Y) {
    
    all_Hsps_only$coloc_w_fibril[i] <- "Yes"
    
    matched_fibril <- only_these_fibrils %>% filter(X_and_Y == current_XY_of_Hsp)
    
    all_Hsps_only$coloc_Fibril_ID[i] <-
      paste(unique(matched_fibril$Fibril_ID), collapse = "&")
    all_Hsps_only$coloc_Fibril_ContourID[i] <-
      paste(unique(matched_fibril$Contour.ID), collapse = "&")
    
  } else {
    all_Hsps_only$coloc_w_fibril[i] <- "No"
  }
}

# -----------------------------------------------------------------------------
# 4. Sanity check: do the Hsp and fibril datasets cover the same images?
# -----------------------------------------------------------------------------

all_Hsps_identif_codes   <- unique(all_Hsps_only$Image_treatment_type)
all_fibril_identif_codes <- unique(truly_all_fibrils$Image_treatment_type)

extra_elements  <- setdiff(all_Hsps_identif_codes, all_fibril_identif_codes)
sets_are_equal  <- setequal(all_Hsps_identif_codes, all_fibril_identif_codes)

cat("Hsp image codes:   ", length(all_Hsps_identif_codes), "\n")
cat("Fibril image codes:", length(all_fibril_identif_codes), "\n")
cat("Hsp codes with no matching fibril entry:\n")
print(extra_elements)
cat("Datasets cover the same images: ", sets_are_equal, "\n")

if (!sets_are_equal) {
  warning("Hsp and fibril datasets do not cover the same set of images - check extra_elements.")
}

# -----------------------------------------------------------------------------
# 5. Save colocalisation results
# -----------------------------------------------------------------------------

dir.create(coloc_out_folder, showWarnings = FALSE, recursive = TRUE)

write.csv(all_Hsps_only,
          file = paste0(coloc_out_folder, "all_Hsps_only_w_coloc_info.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. Summary tables
# -----------------------------------------------------------------------------

all_Hsps_only <- all_Hsps_only %>%
  mutate(Image_treatment = paste0(just_image_number, "-", Treatment_name))

# per-image summary
coloc_summary <- all_Hsps_only %>%
  group_by(Image_treatment, analysis_method) %>%
  summarize(count_yes   = sum(coloc_w_fibril == "Yes", na.rm = TRUE),
            total_count = n(),
            .groups = "drop")

# per-treatment summary (used for plotting) - full dataset, seeds included
coloc_summary_v2 <- all_Hsps_only %>%
  group_by(Treatment_name, analysis_method) %>%
  summarize(number_of_images     = length(unique(Image_treatment_type)),
            number_of_coloc_Hsps = sum(coloc_w_fibril == "Yes", na.rm = TRUE),
            total_number_of_Hsps = n(),
            .groups = "drop") %>%
  mutate(percent_coloc_w_fibril = (number_of_coloc_Hsps / total_number_of_Hsps) * 100,
         dataset = "full_seeds_inc")

# fibril-centric summary (clean dataset - seeds excluded), per image
fibril_summary <- truly_all_fibrils_clean %>%
  group_by(Treatment_name, just_image_number, analysis_method) %>%
  summarize(
    number_of_fibrils                = length(unique(Fibril_ID)),
    total_number_coloc_chaps         = sum(distance != -1),
    number_of_fibrils_with_coloc_chap = length(unique(Fibril_ID[distance != -1])),
    total_number_of_fibril_pixels    = sum(Length[!duplicated(Fibril_ID)]),
    .groups = "drop"
  ) %>%
  mutate(percent_of_fibrils_w_coloc_chap = (number_of_fibrils_with_coloc_chap / number_of_fibrils) * 100,
         dataset = "clean_no_seeds")

# fibril-centric summary (clean dataset), per treatment - used for plotting
fibril_summary2 <- truly_all_fibrils_clean %>%
  group_by(Treatment_name, analysis_method) %>%
  summarize(
    number_of_fibrils                = length(unique(Fibril_ID)),
    total_number_coloc_chaps         = sum(distance != -1),
    number_of_fibrils_with_coloc_chap = length(unique(Fibril_ID[distance != -1])),
    total_number_of_fibril_pixels    = sum(Length[!duplicated(Fibril_ID)]),
    number_of_images                 = length(unique(just_image_number)),
    .groups = "drop"
  ) %>%
  mutate(percent_of_fibrils_w_coloc_chap = (number_of_fibrils_with_coloc_chap / number_of_fibrils) * 100,
         dataset = "clean_no_seeds")

# -----------------------------------------------------------------------------
# 7. Save summary tables
# -----------------------------------------------------------------------------

dir.create(binding_folder, showWarnings = FALSE, recursive = TRUE)

write.csv(coloc_summary_v2,
          file = paste0(binding_folder, "all_Hsps_coloc_info_summary.csv"),
          row.names = FALSE)

write.csv(fibril_summary,
          file = paste0(binding_folder, "clean_fibrils_coloc_info_summary_images.csv"),
          row.names = FALSE)

write.csv(fibril_summary2,
          file = paste0(binding_folder, "clean_fibrils_coloc_info_summary_treatment.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 8. Optional: detailed per-fibril colocalisation breakdown
#    (not run by default - uncomment if you need per-Fibril_ID coloc counts)
# -----------------------------------------------------------------------------

# fibril_summary_detailed <- truly_all_fibrils_clean %>%
#   group_by(Treatment_name, just_image_number, analysis_method) %>%
#   summarize(
#     number_of_fibrils         = length(unique(Fibril_ID)),
#     total_number_coloc_chaps  = sum(distance != -1),
#     number_of_fibrils_with_coloc = length(unique(Fibril_ID[distance != -1])),
#     coloc_per_fibril = list(
#       data.frame(
#         Fibril_ID   = unique(Fibril_ID),
#         coloc_count = sapply(unique(Fibril_ID), function(id) sum(distance[Fibril_ID == id] != -1))
#       )
#     )
#   ) %>%
#   unnest(cols = c(coloc_per_fibril))

# -----------------------------------------------------------------------------
# 9. Plotting setup
# -----------------------------------------------------------------------------

desired_order <- c("50nM_aBc3D", "100nM_aBc3D", "250nM_aBc3D", "500nM_aBc3D", "1uM_aBc3D")
custom_labels <- str_wrap(c("50 nM aBc3D", "100 nM aBc3D", "250 nM aBc3D", "500 nM aBc3D", "1 uM aBc3D"), 10)

coloc_summary_v2 <- coloc_summary_v2 %>%
  mutate(Treatment_name = factor(Treatment_name, levels = desired_order)) %>%
  arrange(Treatment_name)

fibril_summary2 <- fibril_summary2 %>%
  mutate(Treatment_name = factor(Treatment_name, levels = desired_order)) %>%
  arrange(Treatment_name)

# shared theme object so the plots below look consistent - just a variable, reused with `+`
plot_theme <- theme(
  axis.text        = element_text(colour = "black", size = 12),
  axis.title.x     = element_text(colour = "black", size = 14, face = "bold"),
  axis.title.y     = element_text(colour = "black", size = 14, face = "bold"),
  legend.text      = element_text(colour = "black", size = 12),
  legend.title     = element_text(colour = "black", size = 14, face = "bold")
)

# -----------------------------------------------------------------------------
# 10. Generate plots
# -----------------------------------------------------------------------------

fibril_percent <- ggplot(fibril_summary2, aes(x = Treatment_name, y = percent_of_fibrils_w_coloc_chap, fill = analysis_method)) +
  geom_col(position = "dodge") +
  labs(x = "Concentration of sHsp",
       y = str_wrap("Percentage of fibrils with at least 1 colocalised foci", 30),
       fill = "Analysis method") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  plot_theme
print(fibril_percent)

hsp_percent <- ggplot(coloc_summary_v2, aes(x = Treatment_name, y = percent_coloc_w_fibril, fill = analysis_method)) +
  geom_col(position = "dodge") +
  labs(x = "Concentration of sHsp",
       y = str_wrap("Percentage of Hsp foci colocalised with a fibril", 25),
       fill = "Analysis method") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  plot_theme
print(hsp_percent)

hsp_number <- ggplot(coloc_summary_v2, aes(x = Treatment_name, y = number_of_coloc_Hsps, fill = analysis_method)) +
  geom_col(position = "dodge") +
  labs(x = "Concentration of sHsp",
       y = str_wrap("Number of Hsp foci colocalised with a fibril", 25),
       fill = "Analysis method") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  plot_theme
print(hsp_number)

number_images <- ggplot(coloc_summary_v2, aes(x = Treatment_name, y = number_of_images, fill = analysis_method)) +
  geom_col(position = "dodge") +
  labs(x = "Concentration of sHsp",
       y = str_wrap("Number of images in each treatment", 25),
       fill = "Analysis method") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  plot_theme
print(number_images)

# -----------------------------------------------------------------------------
# 11. Combine + save main figure
# -----------------------------------------------------------------------------

facet_plot <- ggarrange(hsp_number, hsp_percent, fibril_percent, number_images,
                        ncol = 2, nrow = 2, common.legend = TRUE)
print(facet_plot)

ggsave(paste0(binding_folder, experiment_date, "_Figure.png"),
       plot = facet_plot, units = "px", width = 2500, height = 1600, dpi = 300)

# -----------------------------------------------------------------------------
# 12. Optional variant: number_of_fibrils instead of number_of_images,
#     with custom x-axis labels
# -----------------------------------------------------------------------------

number_fibrils <- ggplot(fibril_summary2, aes(x = Treatment_name, y = number_of_fibrils, fill = analysis_method)) +
  geom_col(position = "dodge") +
  scale_x_discrete(labels = custom_labels) +
  labs(x = "Concentration of Hsp27 WT",
       y = str_wrap("Number of fibrils in each treatment", 25),
       fill = "Analysis method") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  plot_theme
print(number_fibrils)

facet_plot_v2 <- ggarrange(hsp_number, hsp_percent, fibril_percent, number_fibrils,
                           ncol = 2, nrow = 2, common.legend = TRUE)
print(facet_plot_v2)

ggsave(paste0(binding_folder, experiment_date, "_Figure_v2.png"),
       plot = facet_plot_v2, units = "px", width = 2500, height = 1600, dpi = 300)

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------

rm(list = ls())
dev.off()
cat("\014")
