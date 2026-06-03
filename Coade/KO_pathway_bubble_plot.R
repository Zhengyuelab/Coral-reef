#!/usr/bin/env Rscript

# ============================================================
# Step 05. Draw KEGG pathway class level-2 bubble plot
# ============================================================

packages <- c("dplyr", "readr", "ggplot2", "stringr")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
})

input_file <- file.path("bubble_plot", "KO_pathway_bubble_input.csv")
output_pdf <- file.path(
  "bubble_plot",
  "KO_Pathway_class_level2_bubble_size_KOcount_color_OccurrenceRatio.pdf"
)

cat("\n--- Step 05: Reading bubble plot input file ---\n")

if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file)
}

df <- readr::read_csv(
  input_file,
  show_col_types = FALSE,
  locale = readr::locale(encoding = "UTF-8")
)

required_cols <- c(
  "State",
  "KEGG_KO",
  "Frequency of occurrence",
  "Pathway_class_level2"
)

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df_clean <- df %>%
  dplyr::mutate(
    State = trimws(as.character(State)),
    KEGG_KO = trimws(as.character(KEGG_KO)),
    Pathway_class_level2 = trimws(as.character(Pathway_class_level2)),
    Frequency_of_occurrence = as.numeric(`Frequency of occurrence`)
  ) %>%
  dplyr::filter(
    !is.na(State), State != "",
    !is.na(KEGG_KO), KEGG_KO != "",
    !is.na(Pathway_class_level2), Pathway_class_level2 != "",
    !is.na(Frequency_of_occurrence)
  )

df_clean$State <- factor(
  df_clean$State,
  levels = c("Good", "Moderate", "Poor")
)

bubble_df <- df_clean %>%
  dplyr::group_by(State, Pathway_class_level2) %>%
  dplyr::summarise(
    KO_count = dplyr::n_distinct(KEGG_KO),
    Total_occurrence = sum(Frequency_of_occurrence, na.rm = TRUE),
    Mean_occurrence = mean(Frequency_of_occurrence, na.rm = TRUE),
    KO_list = paste(sort(unique(KEGG_KO)), collapse = ";"),
    .groups = "drop"
  ) %>%
  dplyr::group_by(State) %>%
  dplyr::mutate(
    KO_ratio = KO_count / sum(KO_count, na.rm = TRUE),
    Occurrence_ratio = Total_occurrence / sum(Total_occurrence, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

pathway_order <- bubble_df %>%
  dplyr::group_by(Pathway_class_level2) %>%
  dplyr::summarise(
    Total_KO_count_all_status = sum(KO_count, na.rm = TRUE),
    Total_occurrence_all_status = sum(Total_occurrence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Total_KO_count_all_status, Total_occurrence_all_status) %>%
  dplyr::pull(Pathway_class_level2)

bubble_df$Pathway_class_level2 <- factor(
  bubble_df$Pathway_class_level2,
  levels = pathway_order
)

p_bubble <- ggplot(
  bubble_df,
  aes(
    x = State,
    y = Pathway_class_level2
  )
) +
  geom_point(
    aes(
      size = KO_count,
      color = Occurrence_ratio
    ),
    alpha = 0.85
  ) +
  scale_size_continuous(
    range = c(3, 10),
    name = "Number of KOs"
  ) +
  scale_color_gradient(
    low = "#FDEBD0",
    high = "#A04000",
    name = "Occurrence ratio\nwithin state"
  ) +
  labs(
    x = "State",
    y = "KEGG pathway class level 2",
    title = "Relative functional composition of state-associated KOs"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_bubble)

ggsave(
  filename = output_pdf,
  plot = p_bubble,
  width = 8,
  height = 6
)

cat("\nStep 05 completed.\n")
cat("Only output PDF:", output_pdf, "\n")
