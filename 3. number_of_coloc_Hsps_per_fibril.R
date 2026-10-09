# =============================================================================
# Graphing number of sHsp foci per fibril
# =============================================================================

library(tidyverse)
library(tidylog)
library(viridisLite)
library(scales)

# -----------------------------------------------------------------------------
# 0. Paths
# -----------------------------------------------------------------------------

experiment_date  <- "Testing code"
parent_folder_C  <- paste0("C:/Nicola/UNIVERSITY/PhD/SM/Fibrils and sHsps/", experiment_date, "/")
setwd(parent_folder_C)

dir.create("analysis_outputs", showWarnings = FALSE, recursive = TRUE)

# pixel -> micron conversion factor for fibril length
px_to_um <- 0.16

# seed-size exclusion threshold (this is the size of "dirt"/seed debris, not real fibrils)
mean_seed_size_um <- 1.207
stdev_seed_um     <- 0.648
seed_cutoff_um    <- mean_seed_size_um + stdev_seed_um

desired_order <- c("50nM_aBc3D", "100nM_aBc3D", "250nM_aBc3D", "500nM_aBc3D", "1uM_aBc3D")
custom_labels <- str_wrap(c("50 nM aBc3D", "100 nM aBc3D", "250 nM aBc3D", "500 nM aBc3D", "1 uM aBc3D"), 10)

# -----------------------------------------------------------------------------
# 1. Load + clean fibril data
# -----------------------------------------------------------------------------

all_fibrils_normal    <- read.csv(paste0("full_data/", experiment_date, "_all_fibrils_normal.csv"))
all_fibrils_scrambled <- read.csv(paste0("full_data/", experiment_date, "_all_fibrils_scrambled.csv"))

all_fibrils_normal <- all_fibrils_normal %>%
  mutate(X_and_Y = paste0(X2, "&", Y2), Length_um = Length * px_to_um)

all_fibrils_scrambled <- all_fibrils_scrambled %>%
  mutate(X_and_Y = paste0(X2, "&", Y2), Length_um = Length * px_to_um)

truly_all_fibrils <- bind_rows(all_fibrils_normal, all_fibrils_scrambled) %>%
  mutate(Image_treatment_type = paste0(just_image_number, "-", Treatment_name, "-", analysis_method))

print(unique(truly_all_fibrils$Image_treatment_type))

# "clean" dataset: excludes anything at/below seed size, i.e. genuine fibrils only
all_fibrils_normal_clean    <- all_fibrils_normal    %>% filter(Length_um > seed_cutoff_um)
all_fibrils_scrambled_clean <- all_fibrils_scrambled %>% filter(Length_um > seed_cutoff_um)

truly_all_fibrils_clean <- bind_rows(all_fibrils_normal_clean, all_fibrils_scrambled_clean) %>%
  mutate(Image_treatment_type = paste0(just_image_number, "-", Treatment_name, "-", analysis_method))

truly_all_fibrils_scrambled <- truly_all_fibrils_clean %>% filter(analysis_method == "scrambled")
truly_all_fibrils_clean     <- truly_all_fibrils_clean %>% filter(analysis_method == "normal")

write.csv(truly_all_fibrils_clean, file = "analysis_outputs/truly_all_fibrils_clean_normal.csv", row.names = FALSE)
write.csv(truly_all_fibrils_scrambled, file = "analysis_outputs/truly_all_fibrils_clean_scrambled.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. Per-fibril summary: foci count + density (normal and scrambled)
# -----------------------------------------------------------------------------

truly_all_fibrils_clean <- truly_all_fibrils_clean %>%
  mutate(Fibril_ID_detailed = paste(Fibril_ID, Image_treatment, sep = "-"))

individual_fibrils_v2 <- truly_all_fibrils_clean %>%
  group_by(Fibril_ID, Treatment_name, just_image_number) %>%
  summarise(fibril_length = unique(Length),
            how_many_coloc_hsp_on_this_fibril = sum(distance != -1),
            hsp_names_all = paste(X_and_Y[distance != -1], collapse = ","),
            fibril_length_um = unique(Length_um),
            hsp_density = how_many_coloc_hsp_on_this_fibril / fibril_length_um,
            Contour_ID = unique(Contour.ID),
            analysis_method = unique(analysis_method),
            .groups = "drop")

write.csv(individual_fibrils_v2, file = "analysis_outputs/individual_fibril_info.csv", row.names = FALSE)

truly_all_fibrils_scrambled <- truly_all_fibrils_scrambled %>%
  mutate(Fibril_ID_detailed = paste(Fibril_ID, Image_treatment, sep = "-"))

individual_fibrils_v2_scrambled <- truly_all_fibrils_scrambled %>%
  group_by(Fibril_ID, Treatment_name, just_image_number) %>%
  summarise(fibril_length = unique(Length),
            how_many_coloc_hsp_on_this_fibril = sum(distance != -1),
            hsp_names_all = paste(X_and_Y[distance != -1], collapse = ","),
            fibril_length_um = unique(Length_um),
            hsp_density = how_many_coloc_hsp_on_this_fibril / fibril_length_um,
            Contour_ID = unique(Contour.ID),
            analysis_method = unique(analysis_method),
            .groups = "drop")

write.csv(individual_fibrils_v2_scrambled, file = "analysis_outputs/individual_fibril_info_scrambled.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Foci-count categories per treatment (normal data only)
# -----------------------------------------------------------------------------

fibrils_w_coloc_hsps_only <- individual_fibrils_v2 %>% filter(how_many_coloc_hsp_on_this_fibril != 0)

all_combinations <- expand.grid(
  Treatment_name  = unique(individual_fibrils_v2$Treatment_name),
  fibril_category = c("0", "1", "2", "3", "4", "5", "5+")
)

fibril_counts <- individual_fibrils_v2 %>%
  group_by(Treatment_name) %>%
  mutate(fibril_category = case_when(
    how_many_coloc_hsp_on_this_fibril == 0 ~ "0",
    how_many_coloc_hsp_on_this_fibril == 1 ~ "1",
    how_many_coloc_hsp_on_this_fibril == 2 ~ "2",
    how_many_coloc_hsp_on_this_fibril == 3 ~ "3",
    how_many_coloc_hsp_on_this_fibril == 4 ~ "4",
    how_many_coloc_hsp_on_this_fibril == 5 ~ "5",
    how_many_coloc_hsp_on_this_fibril > 5  ~ "5+"
  )) %>%
  count(Treatment_name, fibril_category, name = "fibril_count") %>%
  right_join(all_combinations, by = c("Treatment_name", "fibril_category")) %>%
  replace_na(list(fibril_count = 0))

total_fibrils <- fibril_counts %>%
  group_by(Treatment_name) %>%
  summarise(total_count = sum(fibril_count))

fibril_counts_with_percentage <- fibril_counts %>%
  left_join(total_fibrils, by = "Treatment_name") %>%
  mutate(percentage = (fibril_count / total_count) * 100) %>%
  arrange(Treatment_name, fibril_category)

fibril_counts_non_zero <- fibril_counts_with_percentage %>%
  filter(fibril_category != "0") %>%
  mutate(Treatment_name = factor(Treatment_name, levels = desired_order)) %>%
  arrange(Treatment_name)

# -----------------------------------------------------------------------------
# 4. Foci-count plots
# -----------------------------------------------------------------------------

fibril_counts_plot <- ggplot(fibril_counts_non_zero, aes(x = fibril_category, y = fibril_count, fill = Treatment_name)) +
  geom_col(position = "dodge") +
  facet_wrap(~Treatment_name) +
  scale_fill_manual(labels = custom_labels, values = rev(viridis(6))) +
  theme_minimal()
print(fibril_counts_plot)

fibril_percentages_plot <- ggplot(fibril_counts_non_zero, aes(x = fibril_category, y = percentage, fill = Treatment_name)) +
  geom_col(position = "dodge") +
  facet_wrap(~Treatment_name) +
  scale_fill_manual(labels = custom_labels, values = rev(viridis(6))) +
  theme_minimal()
print(fibril_percentages_plot)

fibril_percentages_plot_v2 <- ggplot(fibril_counts_non_zero, aes(x = fibril_category, y = percentage, fill = Treatment_name)) +
  geom_col() +
  scale_fill_manual(labels = custom_labels, values = rev(viridis(6))) +
  theme_minimal()
print(fibril_percentages_plot_v2)

fibril_percentages_plot_v3 <- ggplot(fibril_counts_non_zero, aes(x = Treatment_name, y = percentage, fill = fibril_category)) +
  geom_col(width = 0.75) +
  theme_minimal()
print(fibril_percentages_plot_v3)

ggsave("analysis_outputs/20260203_fibril_how_many_coloc_percent.png",
       plot = fibril_percentages_plot, units = "px", width = 2500, height = 1600, dpi = 300)

# -----------------------------------------------------------------------------
# 5. Foci density plots
# -----------------------------------------------------------------------------

individual_fibrils_non_zero <- individual_fibrils_v2 %>%
  filter(how_many_coloc_hsp_on_this_fibril != 0) %>%
  mutate(Treatment_name = factor(Treatment_name, levels = desired_order)) %>%
  arrange(Treatment_name)

density_boxplots <- ggplot(individual_fibrils_non_zero, aes(x = Treatment_name, y = hsp_density, fill = Treatment_name)) +
  geom_boxplot() +
  labs(x = "Concentration", y = "Hsp foci density") +
  theme_minimal()
print(density_boxplots)

ggsave("analysis_outputs/20260203_fibril_density_boxplots.png",
       plot = density_boxplots, units = "px", width = 2500, height = 1600, dpi = 300)

density_histograms <- ggplot(individual_fibrils_non_zero, aes(x = hsp_density, fill = Treatment_name)) +
  geom_histogram() +
  facet_wrap(~Treatment_name) +
  labs(x = "Hsp foci density on fibrils", y = "Frequency") +
  theme_minimal()
print(density_histograms)

density_histograms_percent <- ggplot(individual_fibrils_non_zero, aes(x = hsp_density, fill = Treatment_name)) +
  geom_histogram(aes(y = after_stat(count) / sum(after_stat(count)))) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  facet_wrap(~Treatment_name) +
  labs(x = "Hsp foci density on fibrils", y = "Percentage") +
  theme_minimal()
print(density_histograms_percent)

ggsave("analysis_outputs/20260203_fibril_density_histograms.png",
       plot = density_histograms, units = "px", width = 2500, height = 1600, dpi = 300)

# -----------------------------------------------------------------------------
# 6. Fibril length vs. foci count
# -----------------------------------------------------------------------------

individual_fibrils_non_zero_classified <- individual_fibrils_non_zero %>%
  mutate(fibril_category = case_when(
    how_many_coloc_hsp_on_this_fibril == 0 ~ "0",
    how_many_coloc_hsp_on_this_fibril == 1 ~ "1",
    how_many_coloc_hsp_on_this_fibril == 2 ~ "2",
    how_many_coloc_hsp_on_this_fibril == 3 ~ "3",
    how_many_coloc_hsp_on_this_fibril == 4 ~ "4",
    how_many_coloc_hsp_on_this_fibril == 5 ~ "5",
    how_many_coloc_hsp_on_this_fibril > 5  ~ "5+"
  ))

write.csv(individual_fibrils_non_zero_classified,
          file = "analysis_outputs/individual_fibril_info_w_classification.csv", row.names = FALSE)

lengths_vs_how_many <- ggplot(individual_fibrils_non_zero_classified, aes(x = fibril_category, y = fibril_length_um, fill = Treatment_name)) +
  geom_boxplot() +
  facet_wrap(~Treatment_name) +
  theme_minimal()
print(lengths_vs_how_many)

lengths_vs_how_many_log <- ggplot(individual_fibrils_non_zero_classified, aes(x = fibril_category, y = log10(fibril_length_um), fill = Treatment_name)) +
  geom_boxplot() +
  facet_wrap(~Treatment_name) +
  theme_minimal()
print(lengths_vs_how_many_log)

ggsave("analysis_outputs/20260203_fibril_lengths_vs_foci_count.png",
       plot = lengths_vs_how_many_log, units = "px", width = 2500, height = 1600, dpi = 300)

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------

rm(list = ls())
dev.off()
cat("\014")