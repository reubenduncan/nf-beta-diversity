# load_feature_table.R
# Shared ingestion helper for beta diversity pipeline.
# Returns a list with:
#   $abund_table  — samples x features numeric matrix
#   $OTU_taxonomy — data.frame with cols Kingdom,Phylum,Class,Order,Family,Genus,Otus

library(phyloseq)
library(stringr)
library(data.table)

# ---------------------------------------------------------------------------
# Internal helper: parse a semicolon-separated taxonomy string into a named
# vector with seven standard ranks. Handles both SILVA (D_0__, D_1__, ...)
# and GTDB (d__, p__, c__, o__, f__, g__, s__) prefixes.
# ---------------------------------------------------------------------------
.parse_tax_string <- function(tax_string) {
  ranks  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Otus")
  result <- setNames(rep("", 7), ranks)
  if (is.na(tax_string) || tax_string == "") return(result)

  parts <- trimws(strsplit(tax_string, ";")[[1]])

  # Strip SILVA / GTDB prefixes
  strip_prefix <- function(x) {
    x <- gsub("^D_0__|^d__", "", x)
    x <- gsub("^D_1__|^p__", "", x)
    x <- gsub("^D_2__|^c__", "", x)
    x <- gsub("^D_3__|^o__", "", x)
    x <- gsub("^D_4__|^f__", "", x)
    x <- gsub("^D_5__|^g__", "", x)
    x <- gsub("^D_6__|^s__", "", x)
    trimws(x)
  }

  parts <- sapply(parts, strip_prefix, USE.NAMES = FALSE)
  # Fill as many ranks as we have, leave the rest ""
  n <- min(length(parts), 7)
  result[seq_len(n)] <- parts[seq_len(n)]
  result
}

# ---------------------------------------------------------------------------
# Internal helper: strip rank prefixes from a taxonomy data.frame that already
# has columns Kingdom … Otus (used after import_biom).
# ---------------------------------------------------------------------------
.strip_tax_df_prefixes <- function(tax_df) {
  tax_df[] <- lapply(tax_df, as.character)
  tax_df[is.na(tax_df)] <- ""

  tax_df$Kingdom <- gsub("D_0__|d__", "", tax_df$Kingdom)
  tax_df$Phylum  <- gsub("D_1__|p__", "", tax_df$Phylum)
  tax_df$Class   <- gsub("D_2__|c__", "", tax_df$Class)
  tax_df$Order   <- gsub("D_3__|o__", "", tax_df$Order)
  tax_df$Family  <- gsub("D_4__|f__", "", tax_df$Family)
  tax_df$Genus   <- gsub("D_5__|g__", "", tax_df$Genus)
  tax_df$Otus    <- gsub("D_6__|s__", "", tax_df$Otus)
  tax_df[] <- lapply(tax_df, trimws)
  tax_df
}

# ---------------------------------------------------------------------------
# Internal helper: prune contaminant / unclassified features.
# ---------------------------------------------------------------------------
.prune_taxonomy <- function(abund_table, tax_df) {
  keep <- !(
    tax_df$Kingdom %in% c("Unassigned") |
    tax_df$Phylum  == ""               |
    tax_df$Order   %in% c("Chloroplast") |
    tax_df$Family  %in% c("Mitochondria")
  )
  if (sum(keep) == 0) stop("No features remain after taxonomy pruning.")
  abund_table <- abund_table[, keep, drop = FALSE]
  tax_df      <- tax_df[keep, , drop = FALSE]
  list(abund_table = abund_table, OTU_taxonomy = tax_df)
}

# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------
load_feature_table <- function(feature_table, input_format, taxonomy_table = NULL) {

  input_format <- tolower(trimws(input_format))
  if (!input_format %in% c("biom", "tsv", "gtdb")) {
    stop("input_format must be one of: biom, tsv, gtdb. Got: ", input_format)
  }

  # ---- BIOM ---------------------------------------------------------------
  if (input_format == "biom") {
    message("Loading BIOM file: ", feature_table)
    if (!file.exists(feature_table))
      stop("BIOM file not found: ", feature_table)

    physeq <- tryCatch(
      import_biom(feature_table),
      error = function(e) stop("Failed to import BIOM file: ", conditionMessage(e))
    )

    abund_mat <- as.matrix(otu_table(physeq))
    # Ensure samples are rows
    if (taxa_are_rows(physeq)) abund_mat <- t(abund_mat)

    tax_raw         <- as.data.frame(tax_table(physeq))
    n_cols          <- ncol(tax_raw)
    rank_names_std  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Otus")
    if (n_cols >= 7) {
      colnames(tax_raw)[1:7] <- rank_names_std
    } else {
      colnames(tax_raw) <- rank_names_std[seq_len(n_cols)]
      for (r in rank_names_std[seq(n_cols + 1, 7)]) tax_raw[[r]] <- ""
    }
    tax_raw <- .strip_tax_df_prefixes(tax_raw)
    # Align features
    common_features <- intersect(colnames(abund_mat), rownames(tax_raw))
    abund_mat <- abund_mat[, common_features, drop = FALSE]
    tax_raw   <- tax_raw[common_features, , drop = FALSE]
    pruned    <- .prune_taxonomy(abund_mat, tax_raw)
    return(pruned)
  }

  # ---- TSV / GTDB ---------------------------------------------------------
  if (input_format %in% c("tsv", "gtdb")) {
    message("Loading TSV feature table: ", feature_table)
    if (!file.exists(feature_table))
      stop("Feature table file not found: ", feature_table)

    ft <- tryCatch(
      as.data.frame(fread(feature_table, sep = "\t", header = TRUE, check.names = FALSE)),
      error = function(e) stop("Failed to read feature table: ", conditionMessage(e))
    )

    # First column is feature IDs
    rownames(ft) <- ft[[1]]
    ft           <- ft[, -1, drop = FALSE]
    ft[]         <- lapply(ft, function(x) suppressWarnings(as.numeric(x)))
    ft[is.na(ft)] <- 0

    # Auto-detect orientation: features as rows if nrow > ncol → transpose
    if (nrow(ft) > ncol(ft)) {
      message("Feature table has more rows than columns — assuming features-as-rows, transposing to samples-as-rows.")
      ft <- as.data.frame(t(ft))
    }

    abund_mat <- as.matrix(ft)
    storage.mode(abund_mat) <- "numeric"

    # ---- Taxonomy table ---------------------------------------------------
    ranks_std <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Otus")

    if (!is.null(taxonomy_table) && taxonomy_table != "" && file.exists(taxonomy_table)) {
      message("Loading taxonomy table: ", taxonomy_table)
      tax_raw <- tryCatch(
        as.data.frame(fread(taxonomy_table, sep = "\t", header = FALSE,
                            check.names = FALSE, fill = TRUE)),
        error = function(e) stop("Failed to read taxonomy table: ", conditionMessage(e))
      )
      # Col 1 = feature ID, col 2 = taxonomy string (semicolon-separated)
      feat_ids   <- as.character(tax_raw[[1]])
      tax_strs   <- as.character(tax_raw[[2]])
      tax_parsed <- as.data.frame(
        t(mapply(.parse_tax_string, tax_strs)),
        stringsAsFactors = FALSE
      )
      rownames(tax_parsed) <- feat_ids
      colnames(tax_parsed) <- ranks_std

      # Align features between abundance table and taxonomy
      common_features <- intersect(colnames(abund_mat), rownames(tax_parsed))
      if (length(common_features) == 0)
        stop("No feature IDs overlap between feature table and taxonomy table.")
      abund_mat  <- abund_mat[, common_features, drop = FALSE]
      tax_parsed <- tax_parsed[common_features, , drop = FALSE]

    } else {
      message("No taxonomy table provided — creating empty taxonomy for all features.")
      tax_parsed <- as.data.frame(
        matrix("", nrow = ncol(abund_mat), ncol = 7,
               dimnames = list(colnames(abund_mat), ranks_std)),
        stringsAsFactors = FALSE
      )
    }

    pruned <- .prune_taxonomy(abund_mat, tax_parsed)
    return(pruned)
  }
}
