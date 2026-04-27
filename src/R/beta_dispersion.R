# beta_dispersion.R
# Pairwise beta dispersion (betadisper) analysis for the beta diversity pipeline.
# Outputs CSVs only — no plotting code.

suppressPackageStartupMessages({
  library(optparse)
  library(phyloseq)
  library(vegan)
  library(ape)
  library(phangorn)
  library(stringr)
})

# ---- Source shared helper ---------------------------------------------------
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
if (is.na(script_dir) || script_dir == "") script_dir <- "."
source(file.path(script_dir, "load_feature_table.R"))

# ---- Option parsing ---------------------------------------------------------
option_list <- list(
  # Primary input (--biom_file kept as alias for backward compat)
  make_option("--feature_table",   type = "character", default = NULL,
              help = "Path to feature table (BIOM / TSV / GTDB) [required]"),
  make_option("--biom_file",       type = "character", default = NULL,
              help = "Alias for --feature_table (backward compat)"),
  make_option("--input_format",    type = "character", default = "biom",
              help = "Input format: biom | tsv | gtdb [default: biom]"),
  make_option("--taxonomy_table",  type = "character", default = "",
              help = "Path to taxonomy TSV (required for tsv/gtdb)"),
  make_option("--meta_table",      type = "character", default = NULL,
              help = "Path to metadata CSV [required]"),
  make_option("--tree_file",       type = "character", default = "",
              help = "Path to phylogenetic tree (.nwk); required for unifrac/wunifrac"),
  make_option("--output_dir",      type = "character", default = ".",
              help = "Output directory [default: .]"),
  # Distance metric
  make_option("--distance_metric", type = "character", default = "bray",
              help = "Distance metric: bray | jaccard | unifrac | wunifrac | aitchison [default: bray]"),
  make_option("--output_metric",  type = "character", default = "",
              help = "Metric tag used in output filenames (set by Nextflow; defaults to distance_metric)"),
  # Label
  make_option("--label",           type = "character", default = "analysis",
              help = "Analysis label used in output filenames [default: analysis]"),
  # Filtering
  make_option("--min_library_size", type = "integer",  default = 5000,
              help = "Minimum reads per sample [default: 5000]"),
  make_option("--exclude_column",  type = "character", default = "",
              help = "Metadata column for sample exclusion"),
  make_option("--exclude_values",  type = "character", default = "",
              help = "Comma-separated values to exclude"),
  # Grouping
  make_option("--group", type = "character", default = "",
              help = "One metadata column, or comma-separated columns to paste as the group label"),
  # p-value adjustment
  make_option("--p_adjust_method",  type = "character", default = "BH",
              help = "p-value adjustment method: BH | bonferroni | holm | none [default: BH]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- Backward compat: --biom_file alias ------------------------------------
if (is.null(opt$feature_table) && !is.null(opt$biom_file)) {
  opt$feature_table <- opt$biom_file
}
if (is.null(opt$feature_table)) stop("--feature_table (or --biom_file) is required.")
if (is.null(opt$meta_table))    stop("--meta_table is required.")

# ---- Validate enums --------------------------------------------------------
valid_formats   <- c("biom", "tsv", "gtdb")
valid_distances <- c("bray", "jaccard", "unifrac", "wunifrac", "aitchison")
valid_padj      <- c("BH", "bonferroni", "holm", "none")

if (!opt$input_format    %in% valid_formats)   stop("--input_format must be one of: ", paste(valid_formats, collapse=", "))
if (!opt$distance_metric %in% valid_distances) stop("--distance_metric must be one of: ", paste(valid_distances, collapse=", "))
if (!opt$p_adjust_method %in% valid_padj)      stop("--p_adjust_method must be one of: ", paste(valid_padj, collapse=", "))

if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)

# ---- Tree / unifrac validation ---------------------------------------------
distance_metric <- opt$distance_metric
tree_needed     <- distance_metric %in% c("unifrac", "wunifrac")
tree_available  <- opt$tree_file != "" && file.exists(opt$tree_file)

if (tree_needed && !tree_available) {
  message("SKIP: '", opt$distance_metric, "' requires a phylogenetic tree but none was provided.",
          " Re-run with --tree_file to enable this metric.")
  quit(save = "no", status = 0)
}

output_metric <- if (nchar(opt$output_metric) > 0) opt$output_metric else distance_metric

# ---- Load data --------------------------------------------------------------
message("Loading feature table...")
ft_data <- load_feature_table(
  feature_table  = opt$feature_table,
  input_format   = opt$input_format,
  taxonomy_table = if (opt$taxonomy_table != "") opt$taxonomy_table else NULL
)
abund_table  <- ft_data$abund_table
feature_taxonomy <- ft_data$feature_taxonomy

message("Loading metadata: ", opt$meta_table)
meta_table <- local({
  sep <- if (grepl("\t", readLines(opt$meta_table, n = 1, warn = FALSE))) "\t" else ","
  read.table(opt$meta_table, header = TRUE, sep = sep, row.names = 1,
             check.names = FALSE, stringsAsFactors = FALSE)
})

# ---- Validate metadata columns ---------------------------------------------
check_col <- function(col, arg) {
  if (col != "" && !col %in% colnames(meta_table))
    stop("Column '", col, "' specified by ", arg, " not found in metadata.")
}
check_col(opt$exclude_column, "--exclude_column")
if (opt$group != "") {
  for (col in trimws(strsplit(opt$group, ",")[[1]])) check_col(col, "--group")
}

# ---- Library size filter ---------------------------------------------------
abund_table <- abund_table[rowSums(abund_table) >= opt$min_library_size, , drop = FALSE]
if (nrow(abund_table) == 0)
  stop("No samples remain after minimum library size filter (", opt$min_library_size, ").")

abund_table  <- abund_table[, colSums(abund_table) > 0, drop = FALSE]

# ---- Align samples ---------------------------------------------------------
common_samples <- intersect(rownames(abund_table), rownames(meta_table))
if (length(common_samples) == 0)
  stop("No samples are shared between the feature table and metadata.")
abund_table  <- abund_table[common_samples, , drop = FALSE]
meta_table   <- meta_table[common_samples, , drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

# ---- Exclusion filter ------------------------------------------------------
if (opt$exclude_column != "" && opt$exclude_values != "") {
  exc_vals    <- trimws(strsplit(opt$exclude_values, ",")[[1]])
  keep_rows   <- !meta_table[[opt$exclude_column]] %in% exc_vals
  meta_table  <- meta_table[keep_rows, , drop = FALSE]
  abund_table <- abund_table[rownames(meta_table), , drop = FALSE]
}

# ---- Build Groups factor ---------------------------------------------------
if (opt$group != "") {
  cols <- trimws(strsplit(opt$group, ",")[[1]])
  meta_table$Groups <- if (length(cols) == 1) {
    as.factor(as.character(meta_table[[cols]]))
  } else {
    as.factor(do.call(paste, c(meta_table[, cols, drop = FALSE], sep = " ")))
  }
} else {
  stop("--group is required for beta dispersion.")
}

# ---- Re-align -------------------------------------------------------------
abund_table  <- abund_table[rownames(meta_table), , drop = FALSE]
abund_table  <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

# ---- Minimum sample count check -------------------------------------------
if (nrow(abund_table) < 3)
  stop("Fewer than 3 samples remain after filtering. Cannot run betadisper.")

n_groups <- nlevels(meta_table$Groups)
if (n_groups < 2)
  stop("Fewer than 2 groups found. Betadisper requires at least 2 groups.")

group_counts <- table(meta_table$Groups)
small_groups <- names(group_counts[group_counts < 2])
if (length(small_groups) > 0)
  message("WARNING: Groups with fewer than 2 samples (will be skipped in pairwise comparisons): ",
          paste(small_groups, collapse = ", "))

# ---- Build phyloseq --------------------------------------------------------
OTU <- otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(feature_taxonomy))
SAM <- sample_data(meta_table)

if (tree_available) {
  feature_tree           <- read.tree(opt$tree_file)
  feature_tree$tip.label <- gsub("'", "", feature_tree$tip.label)
  physeq <- merge_phyloseq(phyloseq(OTU, TAX), SAM, feature_tree)
} else {
  physeq <- merge_phyloseq(phyloseq(OTU, TAX), SAM)
}

# ---- CLR transform helper for aitchison ------------------------------------
compute_dist_sub <- function(physeq_sub, distance_metric) {
  if (distance_metric == "aitchison") {
    mat <- as.matrix(otu_table(physeq_sub))
    if (taxa_are_rows(physeq_sub)) mat <- t(mat)
    mat_clr <- mat + 0.5
    clr_mat <- t(apply(mat_clr, 1, function(x) log(x) - mean(log(x))))
    return(vegdist(clr_mat, method = "euclidean"))
  }
  phyloseq::distance(physeq_sub, method = distance_metric)
}

# ---- Pairwise betadisper ---------------------------------------------------
group_levels <- levels(meta_table$Groups)
if (length(group_levels) < 2)
  stop("Need at least 2 groups for pairwise betadisper.")

pairs <- combn(group_levels, 2)
df    <- NULL

for (i in seq_len(ncol(pairs))) {
  g1 <- pairs[1, i]
  g2 <- pairs[2, i]
  comparison_label <- paste(g1, "-", g2)
  message("Processing comparison: ", comparison_label)

  # Subset physeq to the two groups
  physeq_sub <- tryCatch(
    prune_samples(sample_data(physeq)$Groups %in% c(g1, g2), physeq),
    error = function(e) NULL
  )
  if (is.null(physeq_sub)) {
    message("WARNING: Could not subset physeq for comparison '", comparison_label, "'. Skipping.")
    next
  }
  physeq_sub <- prune_taxa(taxa_sums(physeq_sub) > 0, physeq_sub)

  # Check per-group sample counts
  sub_groups  <- as.character(sample_data(physeq_sub)$Groups)
  count_g1    <- sum(sub_groups == g1)
  count_g2    <- sum(sub_groups == g2)

  if (count_g1 < 2 || count_g2 < 2) {
    message("WARNING: Skipping comparison '", comparison_label,
            "' — fewer than 2 samples in one or both groups (",
            g1, "=", count_g1, ", ", g2, "=", count_g2, ").")
    next
  }

  dist_sub <- tryCatch(
    compute_dist_sub(physeq_sub, distance_metric),
    error = function(e) {
      message("WARNING: Distance computation failed for '", comparison_label, "': ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(dist_sub)) next

  mod <- tryCatch(
    betadisper(dist_sub, sample_data(physeq_sub)$Groups, type = "centroid"),
    error = function(e) {
      message("WARNING: betadisper failed for '", comparison_label, "': ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(mod)) next

  dist_to_centroid <- mod[["distances"]]
  sample_ids       <- names(dist_to_centroid)

  df2 <- data.frame(
    sample               = sample_ids,
    distance_to_centroid = as.numeric(dist_to_centroid),
    Groups               = as.character(meta_table[sample_ids, "Groups"]),
    comparison           = comparison_label,
    stringsAsFactors     = FALSE
  )

  # ANOVA for significance
  aov_res <- tryCatch(
    summary(aov(distance_to_centroid ~ Groups, data = df2))[[1]][["Pr(>F)"]][1],
    error = function(e) NA_real_
  )
  df2$pvalue <- aov_res

  df <- if (is.null(df)) df2 else rbind(df, df2)
}

if (is.null(df) || nrow(df) == 0) {
  message("No valid pairwise comparisons were completed. No output written.")
  quit(status = 0)
}

# ---- p-value adjustment across comparisons ---------------------------------
# One padj per unique comparison pvalue
padj_method     <- opt$p_adjust_method
unique_pvalues  <- tapply(df$pvalue, df$comparison, function(x) x[1])
pvals_vec       <- as.numeric(unique_pvalues)
comp_names      <- names(unique_pvalues)

if (padj_method == "none") {
  padj_vec <- pvals_vec
} else {
  padj_vec <- p.adjust(pvals_vec, method = padj_method)
}

padj_map <- setNames(padj_vec, comp_names)
df$padj  <- padj_map[df$comparison]

# Significant flag
df$significant <- !is.na(df$padj) & df$padj <= 0.05

# ---- Write CSV -------------------------------------------------------------
out_file <- file.path(
  opt$output_dir,
  paste0("Betadisper_", output_metric, "_", opt$label, ".csv")
)
write.csv(df, out_file, row.names = FALSE)
message("Written: ", out_file)

message("beta_dispersion.R complete.")
