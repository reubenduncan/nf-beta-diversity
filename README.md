# beta-diversity

A Nextflow pipeline for beta diversity analysis of microbiome feature tables, supporting multiple distance metrics run in parallel, PCoA/NMDS ordination, PERMANOVA, and beta dispersion.

## Introduction

This pipeline accepts a feature table (BIOM, TSV, or GTDB format) and a sample metadata CSV. For each specified distance metric it runs PCoA or NMDS ordination, computes PERMANOVA (adonis2) if variables are specified, and calculates beta dispersion (betadisper) between groups. Multiple metrics are processed in parallel. An optional `--merge_parquet` flag consolidates all outputs into a single Parquet file with a `distance_metric` column identifying each row's origin.

## Quick start

```bash
nextflow run main.nf \
  --feature_table   /path/to/feature_table.biom \
  --meta_table      /path/to/meta_table.csv \
  --groups_column   Treatment \
  --distance_metric bray,jaccard \
  --label           my_analysis
```

## Parameters

### Input / Output

| Parameter | Default | Description |
|---|---|---|
| `--feature_table` | *(required)* | Path to feature table (BIOM, TSV, or GTDB format) |
| `--meta_table` | *(required)* | Path to sample metadata CSV (first column = sample IDs) |
| `--taxonomy_table` | `""` | Taxonomy TSV (required for `tsv`/`gtdb` input formats) |
| `--tree_file` | `""` | Newick phylogenetic tree (required for `unifrac`/`wunifrac`) |
| `--input_format` | `biom` | `biom` \| `tsv` \| `gtdb` |
| `--output_dir` | `results/` | Directory for output files |

### Filtering

| Parameter | Default | Description |
|---|---|---|
| `--min_library_size` | `5000` | Minimum per-sample read depth; samples below this are dropped |
| `--exclude_column` | `""` | Metadata column used to identify samples for exclusion |
| `--exclude_values` | `""` | Comma-separated values in `exclude_column` to remove |

### Grouping

| Parameter | Default | Description |
|---|---|---|
| `--groups_column` | `""` | Metadata column for the primary grouping variable |
| `--groups_paste_columns` | `""` | Comma-separated columns pasted together to form groups |
| `--type_column` | `""` | Secondary metadata column for shape mapping |
| `--type2_column` | `""` | Tertiary metadata column for ellipse grouping |
| `--connections_column` | `""` | Column identifying sample connections |
| `--subconnections_column` | `""` | Column for secondary connection level |

### Ordination

| Parameter | Default | Description |
|---|---|---|
| `--taxon_rank` | `Feature` | Taxonomic level for feature collation (`Feature` \| `Genus` \| `Family` \| `Order` \| `Class` \| `Phylum`) |
| `--ordination_method` | `pcoa` | `pcoa` (PCoA via cmdscale) \| `nmds` (metaMDS) |
| `--distance_metric` | `bray` | Comma-separated list of metrics: `bray`, `jaccard`, `unifrac`, `wunifrac`, `aitchison` |
| `--ellipse_kind` | `se` | Ellipse method passed to `vegan::ordiellipse`: `se` \| `sd` |
| `--label` | `analysis` | Label appended to output file names |

### PERMANOVA

| Parameter | Default | Description |
|---|---|---|
| `--permanova_variables` | `""` | Comma-separated metadata columns for the PERMANOVA model (leave empty to skip) |
| `--permanova_permutations` | `999` | Number of permutations |
| `--p_adjust_method` | `BH` | P-value adjustment: `BH` \| `bonferroni` \| `holm` \| `none` |

### Output options

| Parameter | Default | Description |
|---|---|---|
| `--merge_parquet` | `false` | Merge all output CSVs into a single Parquet file |

## Outputs

All files are written to `--output_dir`. Metric names are embedded in each filename.

| File | Description |
|---|---|
| `PCOA_coords_{metric}_{level}_{label}.csv` | Ordination coordinates per sample with group metadata |
| `PCOA_ellipses_{metric}_{level}_{label}.csv` | Ellipse polygon coordinates per group |
| `ADONIS_{metric}_{level}_{label}.csv` | PERMANOVA results (R², F, p-value, adjusted p-value per variable) |
| `Betadisper_{metric}_{label}.csv` | Beta dispersion pairwise group comparisons |
| `beta_diversity_{label}.parquet` | All CSVs merged with `analysis`, `table`, and `distance_metric` columns (`--merge_parquet` only) |

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/) (default executor — environment built automatically from `environment.yml`)
- **or** Docker with `-profile docker`
- **or** Singularity with `-profile singularity`
- **or** a local R installation with: `optparse`, `vegan`, `ape`, `phangorn`, `stringr`, `data.table`, `phyloseq`, `arrow`

## Running with a local R installation

Add `-profile` to select your execution environment (conda is used by default if no profile is specified):

```bash
nextflow run main.nf \
  -c nextflow.config \
  --feature_table   /path/to/table.biom \
  --meta_table      /path/to/meta.csv \
  --groups_column   Treatment \
  --distance_metric bray,jaccard \
  --label           my_analysis
```

Available profiles: `conda` (default), `docker`, `singularity`.