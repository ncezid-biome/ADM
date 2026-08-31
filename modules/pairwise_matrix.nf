
/*
 * Module: pairwise_matrix
 *
 * Generate a pairwise difference matrix for all amplicon FASTA files
 *
 * Inputs:
 *   Channel of FASTA files from primersearch workflow
 * Outputs:
 *   CSV matrix with pairwise differences
 *
 */

nextflow.enable.dsl = 2

workflow pairwise_matrix {

    take:
        all_amplicon_fastas   // tuple(sample_id, fasta)

    main:

        all_fasta_ch = all_amplicon_fastas.map { t -> t[1] }.collect()

        // Feed into runPairwiseRow
        row_files = runPairwiseRow(all_fasta_ch.flatten())

        // Feed into mergeRows
        matrix_csv = mergeRows(row_files.collect())

    emit:
        matrix_csv
}

process runPairwiseRow {

    tag "${fasta_file.baseName}"

    input:
        path(fasta_file)

    output:
        path("${fasta_file.baseName}.txt")

    script:
    """
    pairwise_compare.py \
        --query ${fasta_file} \
        --all ${params.amplicon_folder} \
        --primers ${params.primers} \
        --output ${fasta_file.baseName}.txt \
        --diff_only
    """
}

process mergeRows {

    tag "merge"

    input:
        path row_files

    output:
        path "${params.output_csv}"

    script:
    """
    merge_rows.py \
        -i ${row_files.join(' ')} \
        -o ${params.output_csv}
    """
}
