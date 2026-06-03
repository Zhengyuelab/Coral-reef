#!/usr/bin/env Rscript

# ============================================================
# Step 02. Count occurrence of characteristic KOs in MAGs
# ============================================================

packages <- c("dplyr", "readr", "tidyr")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

input_file <- file.path(
  "MAG_KO_presence_absence_results",
  "11_MAG_characteristic_KO_detail.tsv"
)

out_dir <- "KO_occurrence_statistics"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n--- Step 02: Reading MAG-KO detail table ---\n")

if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file)
}

df <- readr::read_tsv(input_file, show_col_types = FALSE)

required_cols <- c("State", "MAG_ID", "KEGG_KO", "KO_presence")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df <- df %>%
  dplyr::mutate(
    State = trimws(as.character(State)),
    MAG_ID = trimws(as.character(MAG_ID)),
    KEGG_KO = trimws(as.character(KEGG_KO)),
    KO_presence = as.numeric(KO_presence),
    KO_copy_number_in_MAG = if ("KO_copy_number_in_MAG" %in% colnames(.)) {
      as.numeric(KO_copy_number_in_MAG)
    } else {
      NA_real_
    }
  ) %>%
  dplyr::filter(
    !is.na(State), State != "",
    !is.na(MAG_ID), MAG_ID != "",
    !is.na(KEGG_KO), KEGG_KO != ""
  )

cat("Rows after cleaning:", nrow(df), "\n")
cat("State counts:\n")
print(table(df$State))

mag_count_by_state <- df %>%
  dplyr::distinct(State, MAG_ID) %>%
  dplyr::group_by(State) %>%
  dplyr::summarise(
    Total_MAGs_in_state = dplyr::n_distinct(MAG_ID),
    .groups = "drop"
  )

ko_occurrence <- df %>%
  dplyr::filter(KO_presence == 1) %>%
  dplyr::group_by(State, KEGG_KO) %>%
  dplyr::summarise(
    MAG_occurrence_count = dplyr::n_distinct(MAG_ID),
    MAG_list = paste(sort(unique(MAG_ID)), collapse = ";"),
    Genus_list = if ("Genus" %in% colnames(df)) {
      paste(sort(unique(Genus)), collapse = ";")
    } else {
      NA_character_
    },
    mean_KO_copy_number_in_MAG = mean(KO_copy_number_in_MAG, na.rm = TRUE),
    total_KO_copy_number = sum(KO_copy_number_in_MAG, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(mag_count_by_state, by = "State") %>%
  dplyr::mutate(
    KO_prevalence_in_state = MAG_occurrence_count / Total_MAGs_in_state
  ) %>%
  dplyr::arrange(
    State,
    dplyr::desc(KO_prevalence_in_state),
    dplyr::desc(MAG_occurrence_count),
    KEGG_KO
  )

wide_count <- ko_occurrence %>%
  dplyr::select(KEGG_KO, State, MAG_occurrence_count) %>%
  tidyr::pivot_wider(
    names_from = State,
    values_from = MAG_occurrence_count,
    values_fill = 0
  )

wide_prevalence <- ko_occurrence %>%
  dplyr::select(KEGG_KO, State, KO_prevalence_in_state) %>%
  tidyr::pivot_wider(
    names_from = State,
    values_from = KO_prevalence_in_state,
    values_fill = 0
  )

readr::write_csv(
  ko_occurrence,
  file.path(out_dir, "KO_occurrence_count_and_prevalence_by_state.csv")
)

readr::write_csv(
  wide_count,
  file.path(out_dir, "KO_occurrence_count_wide_by_state.csv")
)

readr::write_csv(
  wide_prevalence,
  file.path(out_dir, "KO_prevalence_wide_by_state.csv")
)

cat("\nStep 02 completed.\n")
cat("Main output:", file.path(out_dir, "KO_occurrence_count_and_prevalence_by_state.csv"), "\n")
