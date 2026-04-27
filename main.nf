#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// Beta Diversity Pipeline
// PCoA/NMDS ordination + PERMANOVA (adonis2) + beta dispersion
// Runs each distance metric in parallel; optionally merges outputs to Parquet.
// ============================================================================

process BETA_DIVERSITY {
    tag "${distance_metric}:${params.label}"
    publishDir params.output_dir, mode: 'copy'

    input:
    val  distance_metric
    path feature_table
    path meta_table
    path tree_file

    output:
    path "PCOA_coords_*.csv",    optional: true, emit: pcoa_coords
    path "PCOA_ellipses_*.csv",  optional: true, emit: pcoa_ellipses
    path "ADONIS_*.csv",         optional: true, emit: adonis

    script:
    def tree_arg = (tree_file.name != 'NO_TREE') ? "--tree_file ${tree_file}" : ""
    def tax_arg  = params.taxonomy_table ? "--taxonomy_table ${params.taxonomy_table}" : ""
    """
    Rscript ${projectDir}/src/R/beta_diversity.R \\
        --feature_table              ${feature_table} \\
        --input_format               ${params.input_format} \\
        ${tax_arg} \\
        --meta_table                 ${meta_table} \\
        ${tree_arg} \\
        --output_dir                 . \\
        --taxon_rank                ${params.taxon_rank} \\
        --ordination_method          ${params.ordination_method} \\
        --distance_metric            ${distance_metric} \\
        --ellipse_kind               ${params.ellipse_kind} \\
        --label                      ${params.label} \\
        --min_library_size           ${params.min_library_size} \\
        --exclude_column             "${params.exclude_column}" \\
        --exclude_values             "${params.exclude_values}" \\
        --groups_column              "${params.groups_column}" \\
        --groups_paste_columns       "${params.groups_paste_columns}" \\
        --type_column                "${params.type_column}" \\
        --type2_column               "${params.type2_column}" \\
        --connections_column         "${params.connections_column}" \\
        --subconnections_column      "${params.subconnections_column}" \\
        --permanova_variables        "${params.permanova_variables}" \\
        --permanova_permutations     ${params.permanova_permutations} \\
        --p_adjust_method            ${params.p_adjust_method}
    """
}

process BETA_DISPERSION {
    tag "${distance_metric}:${params.label}"
    publishDir params.output_dir, mode: 'copy'

    input:
    val  distance_metric
    path feature_table
    path meta_table
    path tree_file

    output:
    path "Betadisper_*.csv", optional: true, emit: betadisper

    script:
    def tree_arg = (tree_file.name != 'NO_TREE') ? "--tree_file ${tree_file}" : ""
    def tax_arg  = params.taxonomy_table ? "--taxonomy_table ${params.taxonomy_table}" : ""
    """
    Rscript ${projectDir}/src/R/beta_dispersion.R \\
        --feature_table              ${feature_table} \\
        --input_format               ${params.input_format} \\
        ${tax_arg} \\
        --meta_table                 ${meta_table} \\
        ${tree_arg} \\
        --output_dir                 . \\
        --distance_metric            ${distance_metric} \\
        --label                      ${params.label} \\
        --min_library_size           ${params.min_library_size} \\
        --exclude_column             "${params.exclude_column}" \\
        --exclude_values             "${params.exclude_values}" \\
        --groups_column              "${params.groups_column}" \\
        --groups_paste_columns       "${params.groups_paste_columns}" \\
        --p_adjust_method            ${params.p_adjust_method}
    """
}

// ============================================================================
// Process: MERGE_PARQUET — merge all output CSVs into one Parquet file
// ============================================================================
process MERGE_PARQUET {
    tag "${params.label}"
    publishDir params.output_dir, mode: 'copy'

    input:
    path csvs

    output:
    path "beta_diversity_${params.label}.parquet"

    script:
    """
    Rscript ${projectDir}/src/R/merge_parquet.R \\
        --label      '${params.label}' \\
        --output_dir '.'
    """
}

// ============================================================================
// Workflow
// ============================================================================
workflow {
    // Value channels — reused once per metric without being consumed
    feat_ch = Channel.fromPath(params.feature_table, checkIfExists: true).first()
    meta_ch = Channel.fromPath(params.meta_table,    checkIfExists: true).first()
    tree_ch = params.tree_file
        ? Channel.fromPath(params.tree_file, checkIfExists: true).first()
        : Channel.value(file('NO_TREE'))

    // One invocation per metric, all running in parallel
    metrics_ch = Channel.of(params.distance_metric.tokenize(',').collect { it.trim() }).flatten()

    bd_out   = BETA_DIVERSITY(metrics_ch, feat_ch, meta_ch, tree_ch)
    disp_out = BETA_DISPERSION(metrics_ch, feat_ch, meta_ch, tree_ch)

    if (params.merge_parquet) {
        all_csvs = bd_out.pcoa_coords
            .mix(bd_out.pcoa_ellipses)
            .mix(bd_out.adonis)
            .mix(disp_out.betadisper)
            .collect()
        MERGE_PARQUET(all_csvs)
    }
}
