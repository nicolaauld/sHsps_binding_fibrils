# Fibril–sHsp colocalisation and single-molecule analysis

Analysis pipeline for single-molecule TIRF movies of AF647-labelled sHsps (or control proteins) binding to fibrils. It covers image analysis in ImageJ, colocalisation and specific-binding tests in R, photobleaching-step analysis of individual foci with py4bleaching, and linking those results back to the fibril data.

These experiments are designed for AF647-labelled chaperones bound to ASCP-bound fibrils that does not require alignment of the two channels (as the two emissions are in the same channel). 

## Pipeline at a glance

| Step | File | What it does | Needs | Produces |
|---|---|---|---|---|
| 0 | `Fibrils_AF647_macro_w_scramble_width_intensity_annotated.ijm` | Per-image ImageJ analysis: background correction, fibril ridge detection (with width), fibril intensity, Hsp peak finding, colocalisation, trajectories, scrambled control | Raw fibril image + raw Hsp (AF647) movie | Per-image CSVs, trajectories and tifs (see below) |
| 1 | `1__compiling_all_of_the_data.R` | Compiles every image's CSVs for one experiment day into three combined files | Step 0 output on vorg | `full_data/` combined CSVs |
| 2 | `2__testing_for_specific_binding.R` | Tests each Hsp focus for colocalisation with a fibril, real vs scrambled coordinates; summary tables and figure | Step 1 | Colocalisation tables, summaries, figure |
| 3 | `3__number_of_coloc_Hsps_per_fibril.R` | Foci count and density per fibril; length vs foci count | Step 1 | Per-fibril tables and plots |
| 4 | `4__coloc_Hsp_end_or_middle_edge_excluded.R` | Classifies each focus as fibril END or MIDDLE; excludes fibrils crossing the image edge | Steps 1–3 | END/MIDDLE tables and plot |
| 5 | `5__moving_the_trajectory_files_to_local_folder.R` | Copies trajectory CSVs from vorg to a local folder structure | Step 0 trajectories | `imagejresults/` |
| 6 | `6__rename_trajectories.R` | Renames each trajectory column to a unique per-molecule ID | Steps 2 and 5 | Renamed trajectory CSVs (overwritten in place) |
| 7 | py4bleaching (Python) | Photobleaching-step analysis of each trajectory: molecule counts and max fluorescence | Step 6 output | `python_results/calculate_molecule_size/molecule_counts.csv` |
| 8 | `8__matching_R_results_to_py4bleaching_results.R` | Matches py4bleaching results back to each Hsp focus and the END/MIDDLE table | Steps 2, 4 and 7 | `*_w_trajectory_data.csv` files |
| 9 | `9__adding_fibril_intensity_measurements.R` | Compiles fibril intensity measurements and attaches them to individual fibrils | Step 1 and step 0 | Fibril intensity tables |

Steps 2 to 4 only need step 1, so they can be run in any order after it, apart from step 4, which also needs the outputs of steps 2 and 3. Step 9 is independent of steps 5 to 8.

## Folder layout

**On vorg** (network drive), one folder per experiment day:

```
<experiment_date>/
└── <concentration>_<protein>/              e.g. 250nM_aBc3D   (one per treatment)
    └── <timestamp folder>/                 e.g. 20260203_102907_567
        ├── 488 fibrils/
        │   ├── <image number>/             raw fibril image + per-image outputs
        │   └── Extra_fibril_analysis/      #N_Fibril_results_w_intensity.csv
        └── AF647 photobleach/
            ├── <image number>/             raw Hsp movie + Hsp_results.csv
            ├── Colocalisation_analysis/
            ├── Scrambled_fibril_colocalisation_analysis/
            ├── coloc_trajectories/
            └── non-coloc_trajectories/647/
```

**Local** (`.../Fibrils and sHsps/<experiment_date>/`):

```
<experiment_date>/
├── full_data/                  step 1 and step 9 combined CSVs
├── analysis_outputs/           steps 2, 3, 4, 8, 9 outputs
├── testing_for_specific_binding/   step 2 summaries and figures
├── imagejresults/              step 5 output, renamed in step 6
│   └── <protein>/<concentration>/{coloc,non-coloc}/
└── python_results/calculate_molecule_size/molecule_counts.csv   step 7 output
```

## Requirements

- **R packages:** `tidyverse`, `tidylog`, `ggpubr`, `viridisLite`, `scales`, `fs`
- **ImageJ/Fiji** with the Ridge Detection, Peak Finder and Beam Profile Correction plugins, and the `s_m_b.jar` (Single Molecule Biophysics) plugin
- **Extra ImageJ files the macro calls:** `count_colocalized_peaks_Fibril.js` and the "count colocalized peaks Nicola" plugin, `Filter_for_colocal.ijm`, `Filter_for_noncolocal.ijm`, `v2_results_to_ROI_X2.txt`, `results_to_ROI.txt`, `peak_intensity.txt`, `scramble_xy.js`. The macro currently points at `C:/Nicola/UNIVERSITY/PhD/SM/ImageJ macros/` and `.../LRplugins/plugins/Macros/`
- **py4bleaching** (Python)
- Also needed next to the image folders for the macro: `ROI.zip` and the two background images (`fibril_background.tif`, `AF647_background.tif`)

## Step details

### 0. ImageJ macro

Run once per image. It asks you to pick the fibril image, then the Hsp AF647 image, and then:

1. Background-corrects both channels and makes a max projection of the first 20 frames.
2. Detects fibril ridges (Ridge Detection, line width 3.5, sigma 1.2, width estimated at the same time). You delete junctions from the ROI manager before continuing.
3. Measures the intensity of each fibril ROI and saves it to `Extra_fibril_analysis/#N_Fibril_results_w_intensity.csv`.
4. Finds Hsp peaks (Peak Finder, checked by eye).
5. Colocalises Hsp peaks with fibril coordinates (max distance 3 px) and saves the results to `Colocalisation_analysis/`.
6. Extracts photobleaching trajectories for colocalised peaks (`coloc_trajectories/`) and non-colocalised peaks (`non-coloc_trajectories/647/`).
7. Scrambles the Hsp peak coordinates within the x/y range of the real peaks and repeats the colocalisation. This estimates how much colocalisation happens by chance. Results go to `Scrambled_fibril_colocalisation_analysis/`.

### 1. Compile data (`1__compiling_all_of_the_data.R`)

Prompts you (Y/N in the console) for each treatment folder to include, then combines the per-image CSVs. It adds tracking columns (treatment, experiment date, image number, analysis method, `Fibril_ID` = Contour ID × Length, protein).

- **In:** `Colocalisation_analysis/*.csv` (excluding `_colocalized` files) and the `*_Fibril_scramble_Hsp_results.csv`, `*_Hsp_results.csv` and `*_Hsp_scrambled_results.csv` files in `Scrambled_fibril_colocalisation_analysis/`
- **Out** (`full_data/`): `<date>_all_fibrils_normal.csv`, `<date>_all_fibrils_scrambled.csv`, `<date>_all_Hsps.csv`

### 2. Test for specific binding (`2__testing_for_specific_binding.R`)

Matches each Hsp focus to fibril coordinates image by image, for both the real and scrambled analyses. Reports the percentage of Hsp foci on fibrils and of fibrils with at least one focus. The fibril summaries use only fibrils longer than the debris cutoff (1.207 + 0.648 µm); the Hsp summaries use all fibrils.

- **In:** the three step 1 files
- **Out:** `analysis_outputs/all_Hsps_only_w_coloc_info.csv`; in `testing_for_specific_binding/`: `all_Hsps_coloc_info_summary.csv`, `clean_fibrils_coloc_info_summary_images.csv`, `clean_fibrils_coloc_info_summary_treatment.csv`, `<date>_Figure.png`, `<date>_Figure_v2.png`

### 3. Foci per fibril (`3__number_of_coloc_Hsps_per_fibril.R`)

Counts colocalised Hsp foci on each fibril (normal and scrambled), calculates foci density (foci per µm), bins fibrils by foci count (0, 1 to 5, 5+), and plots density and fibril length against foci count.

- **In:** `full_data/<date>_all_fibrils_normal.csv`, `..._scrambled.csv`
- **Out** (`analysis_outputs/`): `truly_all_fibrils_clean_normal.csv`, `truly_all_fibrils_clean_scrambled.csv`, `individual_fibril_info.csv`, `individual_fibril_info_scrambled.csv`, `individual_fibril_info_w_classification.csv`, plus PNGs for foci-count percentages, density boxplots, density histograms and length vs foci count

### 4. END vs MIDDLE (`4__coloc_Hsp_end_or_middle_edge_excluded.R`)

For fibrils at least 12 px long, classifies each colocalised focus as END (within 3 px of either fibril tip) or MIDDLE, counts END and MIDDLE foci per fibril, and calculates foci per µm for each region. Fibrils crossing the image edge are flagged and excluded from the summary. Only fibrils with at least one colocalised focus can be flagged as edge-crossing.

- **In:** `all_Hsps_only_w_coloc_info.csv`, `individual_fibril_info.csv`, `truly_all_fibrils_clean_normal.csv`, and the raw `full_data/<date>_all_fibrils_normal.csv`
- **Out** (`analysis_outputs/`): `all_Hsps_coloc_info_12pixelmin_ENDMIDDLEinfo_edge_info.csv`, `individual_fibril_info_12pixelmin_ENDMIDDLEinfo_edge_info.csv`, `<date>_fibril_positions_of_foci.png`

### 5. Move trajectory files (`5__moving_the_trajectory_files_to_local_folder.R`)

Copies the trajectory CSVs from vorg into a local folder structure, one folder per concentration with `coloc/` and `non-coloc/` subfolders. Files are renamed `<protein>_<concentration>_<coloc|non-coloc><image number>.csv`. Run separately for each fluorophore, because py4bleaching has to be run separately on each differently labelled molecule set.

- **In:** `coloc_trajectories/` and `non-coloc_trajectories/647/` for each treatment on vorg
- **Out:** `imagejresults/<protein>/<concentration>/{coloc,non-coloc}/`

### 6. Rename trajectories (`6__rename_trajectories.R`)

Renames each trajectory column to `<name>_<fibril contour ID>_<X>_<Y>` (contour ID is 0 for non-colocalised foci), so each molecule keeps its identity through py4bleaching. The first two column names are blanked. Non-coloc files are also stripped of any column that matches a colocalised focus.

**The files are overwritten in place.** Run it once per fresh copy of the trajectory files from step 5, otherwise the IDs get appended a second time.

- **In:** `analysis_outputs/all_Hsps_only_w_coloc_info.csv` and the step 5 trajectory CSVs
- **Out:** the same trajectory CSVs, renamed

### 7. py4bleaching

This code is available at 10.5281/zenodo.10616736.

Run on the renamed trajectory CSVs in `imagejresults/`, separately for each fluorophore. The R scripts that follow expect a single results file, `python_results/calculate_molecule_size/molecule_counts.csv`, with these columns: `Contour_ID`, `coordsX`, `coordsY`, `max_fluorescence`, `all_small_mol_count`, `last_step_mol_count`, `single_step_mol_count`, `protein1` (holds the concentration), and `variable1` (`"Non-coloc"` marks non-colocalised molecules).

Molecules removed by py4bleaching's own cleaning step will have no result downstream. A lot more non-coloc molecules are removed than coloc ones.

### 8. Link py4bleaching results (`8__matching_R_results_to_py4bleaching_results.R`)

Matches each molecule in `molecule_counts.csv` back to its Hsp focus using the ID `<contour ID>_<X>_<Y>` (with the concentration appended for non-coloc molecules). Adds max fluorescence and the three molecule counts (all-small, last-step, single-step), and passes the same information onto the END/MIDDLE table.

The non-coloc ID is not fully unique even with the concentration added, so some collisions are unavoidable with this matching approach.

- **In:** `python_results/calculate_molecule_size/molecule_counts.csv`, `all_Hsps_only_w_coloc_info.csv`, `all_Hsps_coloc_info_12pixelmin_ENDMIDDLEinfo_edge_info.csv`
- **Out** (`analysis_outputs/`): `coloc_Hsps_only_w_trajectory_data.csv`, `all_Hsps_coloc_12pixels_ENDMIDDLE_w_trajectory_data.csv`, `non_coloc_Hsps_only_w_trajectory_data.csv`

### 9. Fibril intensity (`9__adding_fibril_intensity_measurements.R`)

Compiles the `Extra_fibril_analysis` intensity CSVs (same Y/N treatment prompt as step 1), then attaches each fibril's mean intensity by position. The intensity measurements are in the same order as the fibril contours, so the script checks that the fibril counts and treatment/image labels line up before attaching them.

- **In:** `488 fibrils/Extra_fibril_analysis/*.csv` on vorg (excluding `_width` files), and `full_data/<date>_all_fibrils_normal.csv`
- **Out:** `full_data/<date>_all_fibril_intensity_data.csv`, and `analysis_outputs/<date>_individual_fibril_info_intensity_data.csv`

## Editing for each new experiment

| What | Where |
|---|---|
| `experiment_date` and folder paths | Hard-coded at the top of every R script |
| The list of treatment folders (including the timestamp folder names) | `folders` list in step 5 |
| Treatment names and plotting order (`desired_order`, `custom_labels`) | Steps 2, 3, 4 |
| Date-stamped figure file names | `ggsave()` calls in step 3 |
| Macro paths to the extra ImageJ scripts | Macro, `runMacro(...)` lines |
| py4bleaching settings | Step 7 |

Thresholds used across scripts: pixel size 0.16 µm/px; colocalisation distance 3 px; END region 3 px; minimum fibril length 12 px (step 4); debris cutoff 1.207 + 0.648 µm (steps 2 and 3).
