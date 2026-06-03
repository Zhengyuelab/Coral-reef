#!/usr/bin/env Rscript

# ============================================================
# Step 04. Prepare input table for KEGG pathway bubble plot
# ============================================================

packages <- c("dplyr", "readr", "stringr")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

occurrence_file <- file.path(
  "KO_occurrence_statistics",
  "KO_occurrence_count_and_prevalence_by_state.csv"
)

annotation_file <- file.path(
  "KEGG_annotation",
  "KO_with_KEGG_pathway_annotation.csv"
)

out_dir <- "bubble_plot"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(occurrence_file)) {
  stop("Cannot find occurrence file: ", occurrence_file)
}
if (!file.exists(annotation_file)) {
  stop("Cannot find annotation file: ", annotation_file)
}

occ <- readr::read_csv(occurrence_file, show_col_types = FALSE)
anno <- readr::read_csv(annotation_file, show_col_types = FALSE)

required_occ_cols <- c("State", "KEGG_KO", "MAG_occurrence_count")
required_anno_cols <- c("KEGG_KO", "Pathway_class_level2")

missing_occ <- setdiff(required_occ_cols, colnames(occ))
missing_anno <- setdiff(required_anno_cols, colnames(anno))

if (length(missing_occ) > 0) {
  stop("Occurrence file missing columns: ", paste(missing_occ, collapse = ", "))
}
if (length(missing_anno) > 0) {
  stop("Annotation file missing columns: ", paste(missing_anno, collapse = ", "))
}

anno_level2 <- anno %>%
  dplyr::filter(!is.na(Pathway_class_level2), Pathway_class_level2 != "") %>%
  dplyr::select(KEGG_KO, Pathway_class_level2) %>%
  dplyr::distinct()

bubble_input <- occ %>%
  dplyr::select(
    State,
    KEGG_KO,
    `Frequency of occurrence` = MAG_occurrence_count
  ) %>%
  dplyr::left_join(anno_level2, by = "KEGG_KO") %>%
  dplyr::filter(!is.na(Pathway_class_level2), Pathway_class_level2 != "") %>%
  dplyr::mutate(
    State = factor(State, levels = c("Good", "Moderate", "Poor")),
    KEGG_KO = trimws(as.character(KEGG_KO)),
    Pathway_class_level2 = trimws(as.character(Pathway_class_level2))
  ) %>%
  dplyr::arrange(State, Pathway_class_level2, KEGG_KO)

readr::write_csv(
  bubble_input,
  file.path(out_dir, "KO_pathway_bubble_input.csv")
)

cat("\nStep 04 completed.\n")
cat("Output:", file.path(out_dir, "KO_pathway_bubble_input.csv"), "\n")
