# merge_parquet.R
# Reads all beta-diversity CSVs staged into the Nextflow work directory,
# classifies each by type and extracts the distance metric from the filename,
# wide-merges with NA-fill, and writes a single Parquet file.

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
})

option_list <- list(
  make_option("--label",
              type = "character", default = "analysis",
              help = "Label used in the output filename [default: analysis]"),
  make_option("--output_dir",
              type = "character", default = ".",
              help = "Directory for the output Parquet file [default: .]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!dir.exists(opt$output_dir))
  dir.create(opt$output_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Classify CSV and extract distance metric from basename.
# All output filenames embed the metric immediately after their prefix, so we
# use a prefix-specific regex for each type to avoid off-by-one on underscores.
#   PCOA_coords_{metric}_{level}_{label}.csv
#   PCOA_ellipses_{metric}_{level}_{label}.csv
#   ADONIS_{metric}_{level}_{label}.csv
#   Betadisper_{metric}_{label}.csv
# ---------------------------------------------------------------------------
classify_csv <- function(fname) {
  b <- basename(fname)
  if      (grepl("^PCOA_coords_",   b)) list(analysis = "ordination", table = "pcoa_coords",   metric = sub("^PCOA_coords_([^_]+)_.*\\.csv$",   "\\1", b))
  else if (grepl("^PCOA_ellipses_", b)) list(analysis = "ordination", table = "pcoa_ellipses", metric = sub("^PCOA_ellipses_([^_]+)_.*\\.csv$", "\\1", b))
  else if (grepl("^ADONIS_",        b)) list(analysis = "permanova",  table = "adonis",        metric = sub("^ADONIS_([^_]+)_.*\\.csv$",        "\\1", b))
  else if (grepl("^Betadisper_",    b)) list(analysis = "dispersion", table = "betadisper",    metric = sub("^Betadisper_([^_]+)_.*\\.csv$",    "\\1", b))
  else                                   list(analysis = "unknown",    table = "unknown",       metric = NA_character_)
}

# ---------------------------------------------------------------------------
# Read all CSVs in the work directory
# ---------------------------------------------------------------------------
csv_files <- list.files(".", pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
if (length(csv_files) == 0)
  stop("No CSV files found in working directory.")

message("Found ", length(csv_files), " CSV file(s) to merge.")

frames <- lapply(csv_files, function(f) {
  cls <- classify_csv(f)
  df  <- tryCatch(
    read.csv(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) { warning("Failed to read ", basename(f), ": ", conditionMessage(e)); NULL }
  )
  if (is.null(df)) return(NULL)

  # Strip row-index artifacts
  if (ncol(df) > 0 && colnames(df)[1] == "") df <- df[, -1, drop = FALSE]
  idx_cols <- grepl("^X$|^X\\.\\d+$", colnames(df))
  if (any(idx_cols)) df <- df[, !idx_cols, drop = FALSE]

  df$analysis       <- cls$analysis
  df$table          <- cls$table
  df$distance_metric <- cls$metric
  df$source_file    <- basename(f)

  message("  ", basename(f), " → analysis=", cls$analysis,
          ", table=", cls$table, ", metric=", cls$metric,
          " (", nrow(df), " rows)")
  df
})

frames <- Filter(Negate(is.null), frames)
if (length(frames) == 0)
  stop("All CSV files failed to read — cannot write Parquet.")

# ---------------------------------------------------------------------------
# Wide merge with NA-fill
# ---------------------------------------------------------------------------
all_cols       <- unique(unlist(lapply(frames, colnames)))
frames_aligned <- lapply(frames, function(df) {
  df[setdiff(all_cols, colnames(df))] <- NA
  df[, all_cols, drop = FALSE]
})

merged <- do.call(rbind, frames_aligned)
rownames(merged) <- NULL

# ---------------------------------------------------------------------------
# Write Parquet
# ---------------------------------------------------------------------------
out_path <- file.path(opt$output_dir,
                      paste0("beta_diversity_", opt$label, ".parquet"))
arrow::write_parquet(merged, out_path)
message("Written: ", out_path,
        " (", nrow(merged), " rows x ", ncol(merged), " columns)")
