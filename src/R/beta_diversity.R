# beta_diversity.R
# PCoA / NMDS ordination + PERMANOVA (adonis2) for the beta diversity pipeline.
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
  # Taxonomy level
  make_option("--taxon_rank",     type = "character", default = "Feature",
              help = "Taxonomy level: Feature | Genus | Family | Order | Class | Phylum [default: Feature]"),
  # Ordination
  make_option("--ordination_method", type = "character", default = "pcoa",
              help = "Ordination method: pcoa | nmds [default: pcoa]"),
  make_option("--distance_metric", type = "character", default = "bray",
              help = "Distance metric: bray | jaccard | unifrac | wunifrac | aitchison [default: bray]"),
  # Ellipses
  make_option("--ellipse_kind",    type = "character", default = "se",
              help = "Ellipse type: sd | se [default: se]"),
  # Label / naming
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
  make_option("--type",  type = "character", default = "",
              help = "One metadata column, or comma-separated columns to paste as the point-style label (optional)"),
  # PERMANOVA
  make_option("--permanova_variables",    type = "character", default = "",
              help = "Comma-separated metadata columns for PERMANOVA (adonis2)"),
  make_option("--permanova_permutations", type = "integer",   default = 999,
              help = "Number of PERMANOVA permutations [default: 999]"),
  make_option("--p_adjust_method",        type = "character", default = "BH",
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
valid_formats    <- c("biom", "tsv", "gtdb")
valid_ordination <- c("pcoa", "nmds")
valid_distances  <- c("bray", "jaccard", "unifrac", "wunifrac", "aitchison")
valid_kinds      <- c("sd", "se")
valid_levels     <- c("Feature", "Genus", "Family", "Order", "Class", "Phylum")
valid_padj       <- c("BH", "bonferroni", "holm", "none")

if (!opt$input_format    %in% valid_formats)    stop("--input_format must be one of: ", paste(valid_formats, collapse=", "))
if (!opt$ordination_method %in% valid_ordination) stop("--ordination_method must be one of: ", paste(valid_ordination, collapse=", "))
if (!opt$distance_metric %in% valid_distances)  stop("--distance_metric must be one of: ", paste(valid_distances, collapse=", "))
if (!opt$ellipse_kind    %in% valid_kinds)      stop("--ellipse_kind must be one of: sd, se")
if (!opt$taxon_rank     %in% valid_levels)     stop("--taxon_rank must be one of: ", paste(valid_levels, collapse=", "))
if (!opt$p_adjust_method %in% valid_padj)       stop("--p_adjust_method must be one of: ", paste(valid_padj, collapse=", "))

if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)

# ---- Tree / unifrac validation ---------------------------------------------
distance_metric <- opt$distance_metric
tree_needed     <- distance_metric %in% c("unifrac", "wunifrac")
tree_available  <- opt$tree_file != "" && file.exists(opt$tree_file)

if (tree_needed && !tree_available) {
  message("WARNING: Falling back from '", distance_metric, "' to 'bray' (no tree file).")
  distance_metric <- "bray"
  tree_needed     <- FALSE
}

if (tree_needed && opt$taxon_rank != "Feature") {
  message("WARNING: Falling back from '", distance_metric,
          "' to 'bray' (phylogenetic distances require --taxon_rank Feature).")
  distance_metric <- "bray"
  tree_needed     <- FALSE
}

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
check_col(opt$exclude_column,        "--exclude_column")
if (opt$group != "") {
  for (col in trimws(strsplit(opt$group, ",")[[1]])) check_col(col, "--group")
}
if (opt$type != "") {
  for (col in trimws(strsplit(opt$type, ",")[[1]])) check_col(col, "--type")
}
if (opt$permanova_variables != "") {
  perm_vars <- trimws(strsplit(opt$permanova_variables, ",")[[1]])
  for (col in perm_vars) check_col(col, "--permanova_variables")
}

# ---- Library size filter ---------------------------------------------------
abund_table <- abund_table[rowSums(abund_table) >= opt$min_library_size, , drop = FALSE]
if (nrow(abund_table) == 0)
  stop("No samples remain after minimum library size filter (", opt$min_library_size, ").")

# ---- Remove zero-sum features after size filter ----------------------------
abund_table  <- abund_table[, colSums(abund_table) > 0, drop = FALSE]

# ---- Align samples between table and metadata ------------------------------
common_samples <- intersect(rownames(abund_table), rownames(meta_table))
if (length(common_samples) == 0)
  stop("No samples are shared between the feature table and metadata.")
abund_table  <- abund_table[common_samples, , drop = FALSE]
meta_table   <- meta_table[common_samples, , drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

# ---- Exclusion filter ------------------------------------------------------
if (opt$exclude_column != "" && opt$exclude_values != "") {
  exc_vals   <- trimws(strsplit(opt$exclude_values, ",")[[1]])
  keep_rows  <- !meta_table[[opt$exclude_column]] %in% exc_vals
  meta_table <- meta_table[keep_rows, , drop = FALSE]
  abund_table <- abund_table[rownames(meta_table), , drop = FALSE]
}

# ---- Resolve grouping columns ----------------------------------------------
resolve_columns <- function(param_val, param_name, meta) {
  cols <- trimws(strsplit(param_val, ",")[[1]])
  missing <- setdiff(cols, colnames(meta))
  if (length(missing) > 0)
    stop(param_name, " references columns not in metadata: ", paste(missing, collapse = ", "))
  if (length(cols) == 1) {
    as.factor(as.character(meta[[cols]]))
  } else {
    as.factor(do.call(paste, c(meta[, cols, drop = FALSE], sep = " ")))
  }
}

if (opt$group != "") {
  meta_table$Groups <- resolve_columns(opt$group, "--group", meta_table)
} else {
  meta_table$Groups <- as.factor(rep("All", nrow(meta_table)))
  message("No --group specified — all samples assigned to group 'All'.")
}

if (opt$type != "") {
  meta_table$Type <- resolve_columns(opt$type, "--type", meta_table)
} else {
  meta_table$Type <- NULL
}

# ---- Re-align after all filtering ------------------------------------------
abund_table  <- abund_table[rownames(meta_table), , drop = FALSE]
abund_table  <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

# ---- Minimum sample count check --------------------------------------------
if (nrow(abund_table) < 3)
  stop("Fewer than 3 samples remain after filtering. Cannot run ordination.")

# ---- Collate at taxonomic level --------------------------------------------
taxon_rank <- opt$taxon_rank
if (taxon_rank == "Feature") {
  new_abund_table <- abund_table
} else {
  lvl_list        <- unique(feature_taxonomy[, taxon_rank])
  new_abund_table <- NULL
  for (i in lvl_list) {
    feature_idx <- rownames(feature_taxonomy)[feature_taxonomy[, taxon_rank] == i]
    tmp <- data.frame(rowSums(abund_table[, feature_idx, drop = FALSE]))
    colnames(tmp) <- if (i == "") "__Unknowns__" else i
    new_abund_table <- if (is.null(new_abund_table)) tmp else cbind(new_abund_table, tmp)
  }
}
abund_table <- as.data.frame(as.matrix(new_abund_table))

# ---- Build phyloseq object -------------------------------------------------
OTU <- otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)

# When collapsed to a non-feature rank, feature_taxonomy still has ASV-level
# row names and no longer matches the collapsed abund_table columns.
# Rebuild a minimal taxonomy aligned to the current column names.
if (taxon_rank != "Feature") {
  rank_names_std <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Feature")
  feature_taxonomy <- data.frame(
    matrix("", nrow = ncol(abund_table), ncol = 7,
           dimnames = list(colnames(abund_table), rank_names_std)),
    stringsAsFactors = FALSE
  )
  feature_taxonomy[[taxon_rank]] <- colnames(abund_table)
  feature_taxonomy$Feature       <- colnames(abund_table)
}
TAX <- tax_table(as.matrix(feature_taxonomy))
SAM <- sample_data(meta_table)

if (tree_needed && tree_available) {
  feature_tree           <- read.tree(opt$tree_file)
  feature_tree$tip.label <- gsub("'", "", feature_tree$tip.label)
  physeq <- merge_phyloseq(phyloseq(OTU, TAX), SAM, feature_tree)
} else {
  physeq <- merge_phyloseq(phyloseq(OTU, TAX), SAM)
}

# ---- Compute distance matrix -----------------------------------------------
message("Computing distance matrix (", distance_metric, ")...")

compute_dist <- function(physeq, distance_metric, abund_mat) {
  if (distance_metric == "aitchison") {
    # CLR transform then Euclidean
    mat  <- as.matrix(otu_table(physeq))
    if (taxa_are_rows(physeq)) mat <- t(mat)
    # Add pseudo-count to avoid log(0)
    mat_clr <- mat + 0.5
    clr_mat  <- t(apply(mat_clr, 1, function(x) log(x) - mean(log(x))))
    return(vegdist(clr_mat, method = "euclidean"))
  }
  phyloseq::distance(physeq, method = distance_metric)
}

dist_mat <- tryCatch(
  compute_dist(physeq, distance_metric, abund_table),
  error = function(e) stop("Failed to compute distance matrix: ", conditionMessage(e))
)

# ---- Ordination ------------------------------------------------------------
ordination_method <- tolower(opt$ordination_method)
message("Running ", toupper(ordination_method), " ordination...")

# Ellipse helper
veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta  <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

sol        <- NULL
stress_val <- NA_real_

if (ordination_method == "pcoa") {
  sol <- tryCatch(
    cmdscale(dist_mat, eig = TRUE),
    error = function(e) {
      message("ERROR: PCoA (cmdscale) failed: ", conditionMessage(e))
      NULL
    }
  )
} else if (ordination_method == "nmds") {
  sol <- tryCatch(
    metaMDS(dist_mat, k = 2, trymax = 100, trace = FALSE),
    error = function(e) {
      message("ERROR: NMDS (metaMDS) failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(sol)) stress_val <- sol$stress
}

# ---- Build coords data.frame -----------------------------------------------
if (!is.null(sol)) {

  if (ordination_method == "pcoa") {
    coords        <- sol$points[, 1:2, drop = FALSE]
    eig_total     <- sum(abs(sol$eig))
    variance_dim1 <- (abs(sol$eig[1]) / eig_total) * 100
    variance_dim2 <- (abs(sol$eig[2]) / eig_total) * 100
  } else {
    coords        <- sol$points[, 1:2, drop = FALSE]
    variance_dim1 <- NA_real_
    variance_dim2 <- NA_real_
  }

  PCOA <- data.frame(
    x      = coords[, 1],
    y      = coords[, 2],
    meta_table,
    stringsAsFactors = FALSE
  )

  # ---- CSV 1: ordination coordinates -----------------------------------------
  coords_out <- data.frame(
    sample            = rownames(PCOA),
    x                 = PCOA$x,
    y                 = PCOA$y,
    Groups            = as.character(PCOA$Groups),
    variance_dim1     = variance_dim1,
    variance_dim2     = variance_dim2,
    ordination_method = ordination_method,
    stress            = stress_val,
    stringsAsFactors  = FALSE
  )
  if (!is.null(meta_table$Type))
    coords_out$Type <- as.character(PCOA$Type)

  coords_file <- file.path(
    opt$output_dir,
    paste0("PCOA_coords_", distance_metric, "_", taxon_rank, "_", opt$label, ".csv")
  )
  write.csv(coords_out, coords_file, row.names = FALSE)
  message("Written: ", coords_file)

  # ---- CSV 2: ellipse path points --------------------------------------------
  groups_levels <- levels(factor(PCOA$Groups))

  # Generate ellipses via ordiellipse (suppress graphics device noise)
  pdf(nullfile())
  ord_ell <- tryCatch(
    ordiellipse(
      ord      = if (ordination_method == "pcoa") sol else sol,
      groups   = PCOA$Groups,
      display  = "sites",
      kind     = opt$ellipse_kind,
      conf     = 0.95,
      label    = FALSE
    ),
    error = function(e) {
      message("WARNING: ordiellipse failed: ", conditionMessage(e))
      NULL
    }
  )
  dev.off()

  df_ell <- data.frame()
  if (!is.null(ord_ell)) {
    for (h in groups_levels) {
      if (h %in% names(ord_ell)) {
        tryCatch({
          sub_df <- PCOA[PCOA$Groups == h, ]
          if (nrow(sub_df) >= 2) {
            ell_pts <- veganCovEllipse(
              ord_ell[[h]]$cov,
              ord_ell[[h]]$center,
              ord_ell[[h]]$scale
            )
            df_ell <- rbind(df_ell,
              data.frame(
                x      = ell_pts[, 1],
                y      = ell_pts[, 2],
                Groups = h,
                stringsAsFactors = FALSE
              )
            )
          }
        }, error = function(e) {
          message("WARNING: ellipse skipped for group '", h, "': ", conditionMessage(e))
        })
      }
    }
  }

  ell_file <- file.path(
    opt$output_dir,
    paste0("PCOA_ellipses_", distance_metric, "_", taxon_rank, "_", opt$label, ".csv")
  )
  write.csv(df_ell, ell_file, row.names = FALSE)
  message("Written: ", ell_file)

  # ---- CSV 3: ADONIS (PERMANOVA) ---------------------------------------------
  if (opt$permanova_variables != "") {
    perm_vars <- trimws(strsplit(opt$permanova_variables, ",")[[1]])
    message("Running PERMANOVA for: ", paste(perm_vars, collapse = ", "))

    # Align metadata to samples in distance matrix
    meta_perm <- meta_table[rownames(as.matrix(dist_mat)), , drop = FALSE]

    adonis_result <- tryCatch(
      adonis2(
        as.formula(paste("dist_mat ~", paste(perm_vars, collapse = "+"))),
        data         = meta_perm,
        permutations = opt$permanova_permutations
      ),
      error = function(e) {
        message("ERROR: PERMANOVA (adonis2) failed: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(adonis_result)) {
      adonis_df           <- as.data.frame(adonis_result)
      adonis_df           <- cbind(variable = rownames(adonis_df), adonis_df)
      rownames(adonis_df) <- NULL
      colnames(adonis_df) <- gsub("Pr\\(>F\\)", "Pr.F", colnames(adonis_df))

      # p-value adjustment (only for real term rows, not Total/Residual)
      real_rows  <- !is.na(adonis_df$Pr.F)
      padj_vals  <- rep(NA_real_, nrow(adonis_df))
      padj_method <- opt$p_adjust_method
      if (padj_method == "none") {
        padj_vals[real_rows] <- adonis_df$Pr.F[real_rows]
      } else {
        padj_vals[real_rows] <- p.adjust(adonis_df$Pr.F[real_rows], method = padj_method)
      }
      adonis_df$padj <- padj_vals

      adonis_file <- file.path(
        opt$output_dir,
        paste0("ADONIS_", distance_metric, "_", taxon_rank, "_", opt$label, ".csv")
      )
      write.csv(adonis_df, adonis_file, row.names = FALSE)
      message("Written: ", adonis_file)
    }
  }

} else {
  message("Ordination returned NULL — no output files written.")
}

message("beta_diversity.R complete.")
