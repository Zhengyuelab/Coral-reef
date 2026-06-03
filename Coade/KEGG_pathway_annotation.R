#!/usr/bin/env Rscript

# ============================================================
# Step 03. Annotate characteristic KOs with KEGG pathway information
# ============================================================
# This script queries KEGG REST API and therefore requires internet access.
# ============================================================

packages <- c("dplyr", "tidyr", "readr", "stringr")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

input_file <- file.path(
  "MAG_KO_presence_absence_results",
  "13_characteristic_KO_prevalence_heatmap_matrix.tsv"
)

out_dir <- "KEGG_annotation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n--- Step 03: Reading characteristic KO table ---\n")

if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file)
}

ko_df <- readr::read_tsv(input_file, show_col_types = FALSE)

if (!"KEGG_KO" %in% colnames(ko_df)) {
  stop("Input file must contain a column named KEGG_KO.")
}

ko_df$KEGG_KO <- trimws(as.character(ko_df$KEGG_KO))

ko_list <- unique(ko_df$KEGG_KO)
ko_list <- ko_list[!is.na(ko_list) & ko_list != ""]

cat("Number of unique KOs:", length(ko_list), "\n")

safe_readLines <- function(url) {
  cat("Querying:", url, "\n")
  x <- tryCatch(
    readLines(url, warn = FALSE),
    error = function(e) {
      warning("Failed to query: ", url)
      character(0)
    }
  )
  Sys.sleep(0.2)
  x
}

cat("\n--- Downloading KO descriptions from KEGG ---\n")

ko_desc_lines <- safe_readLines("https://rest.kegg.jp/list/ko")

ko_desc <- data.frame()

if (length(ko_desc_lines) > 0) {
  ko_desc <- data.frame(raw = ko_desc_lines, stringsAsFactors = FALSE) %>%
    tidyr::separate(
      raw,
      into = c("KEGG_KO_raw", "KO_description"),
      sep = "\t",
      remove = TRUE,
      fill = "right"
    ) %>%
    dplyr::mutate(
      KEGG_KO = gsub("^ko:", "", KEGG_KO_raw),
      KO_description = trimws(KO_description)
    ) %>%
    dplyr::select(KEGG_KO, KO_description)
}

ko_desc <- ko_desc %>%
  dplyr::filter(KEGG_KO %in% ko_list)

cat("KO descriptions matched:", nrow(ko_desc), "\n")

cat("\n--- Mapping KOs to KEGG pathways ---\n")

ko_pathway_list <- list()

for (ko in ko_list) {
  url <- paste0("https://rest.kegg.jp/link/pathway/ko:", ko)
  lines <- safe_readLines(url)
  if (length(lines) == 0) next

  tmp <- data.frame(raw = lines, stringsAsFactors = FALSE) %>%
    tidyr::separate(
      raw,
      into = c("KEGG_KO_raw", "Pathway_ID_raw"),
      sep = "\t",
      remove = TRUE,
      fill = "right"
    ) %>%
    dplyr::mutate(
      KEGG_KO = gsub("^ko:", "", KEGG_KO_raw),
      Pathway_ID = gsub("^path:", "", Pathway_ID_raw)
    ) %>%
    dplyr::select(KEGG_KO, Pathway_ID)

  ko_pathway_list[[ko]] <- tmp
}

ko_pathway <- dplyr::bind_rows(ko_pathway_list)

cat("KO-pathway links:", nrow(ko_pathway), "\n")
cat("KOs mapped to at least one pathway:", length(unique(ko_pathway$KEGG_KO)), "\n")

cat("\n--- Downloading KEGG pathway names ---\n")

pathway_lines <- safe_readLines("https://rest.kegg.jp/list/pathway/ko")

pathway_names <- data.frame()

if (length(pathway_lines) > 0) {
  pathway_names <- data.frame(raw = pathway_lines, stringsAsFactors = FALSE) %>%
    tidyr::separate(
      raw,
      into = c("Pathway_ID_raw", "Pathway_name"),
      sep = "\t",
      remove = TRUE,
      fill = "right"
    ) %>%
    dplyr::mutate(
      Pathway_ID = gsub("^path:", "", Pathway_ID_raw),
      Pathway_name = trimws(Pathway_name)
    ) %>%
    dplyr::select(Pathway_ID, Pathway_name)
}

cat("\n--- Downloading KEGG pathway classes ---\n")

unique_pathways <- unique(ko_pathway$Pathway_ID)
pathway_class_list <- list()

for (pid in unique_pathways) {
  url <- paste0("https://rest.kegg.jp/get/", pid)
  lines <- safe_readLines(url)
  if (length(lines) == 0) next

  class_line <- lines[grepl("^CLASS", lines)]
  class_text <- if (length(class_line) == 0) NA_character_ else gsub("^CLASS\\s+", "", class_line[1])

  pathway_class_list[[pid]] <- data.frame(
    Pathway_ID = pid,
    Pathway_class = class_text,
    stringsAsFactors = FALSE
  )
}

pathway_class <- dplyr::bind_rows(pathway_class_list) %>%
  tidyr::separate(
    Pathway_class,
    into = c("Pathway_class_level1", "Pathway_class_level2"),
    sep = "; ",
    remove = FALSE,
    fill = "right"
  )

cat("\n--- Building final KO annotation table ---\n")

ko_annotation <- ko_pathway %>%
  dplyr::left_join(ko_desc, by = "KEGG_KO") %>%
  dplyr::left_join(pathway_names, by = "Pathway_ID") %>%
  dplyr::left_join(pathway_class, by = "Pathway_ID") %>%
  dplyr::arrange(KEGG_KO, Pathway_ID)

readr::write_csv(
  ko_annotation,
  file.path(out_dir, "KO_KEGG_pathway_annotation.csv")
)

ko_df_annotated <- ko_df %>%
  dplyr::left_join(ko_annotation, by = "KEGG_KO")

readr::write_csv(
  ko_df_annotated,
  file.path(out_dir, "KO_with_KEGG_pathway_annotation.csv")
)

if ("State" %in% colnames(ko_df_annotated)) {
  pathway_summary <- ko_df_annotated %>%
    dplyr::filter(!is.na(Pathway_ID), !is.na(Pathway_name)) %>%
    dplyr::group_by(
      State,
      Pathway_ID,
      Pathway_name,
      Pathway_class_level1,
      Pathway_class_level2
    ) %>%
    dplyr::summarise(
      KO_count = dplyr::n_distinct(KEGG_KO),
      KO_list = paste(sort(unique(KEGG_KO)), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::group_by(State) %>%
    dplyr::mutate(KO_ratio = KO_count / sum(KO_count)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(State, dplyr::desc(KO_count))

  readr::write_csv(
    pathway_summary,
    file.path(out_dir, "KO_KEGG_pathway_summary_by_state.csv")
  )
}

mapped_kos <- unique(ko_annotation$KEGG_KO)
unmapped_kos <- setdiff(ko_list, mapped_kos)

readr::write_csv(
  data.frame(KEGG_KO = unmapped_kos),
  file.path(out_dir, "KO_unmapped_KOs.csv")
)

cat("\nStep 03 completed.\n")
cat("Main output:", file.path(out_dir, "KO_with_KEGG_pathway_annotation.csv"), "\n")
