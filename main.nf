nextflow.enable.dsl = 2

include { download_reads } from './modules/download_reads'
include { shovill_assembly } from './modules/shovill_assembly'
include { primersearch_workflow } from './modules/primersearch_workflow'  // new module
include { pairwise_matrix }   from './modules/pairwise_matrix'

workflow {

    /*
     * Build SRR channel ONCE
     */
    srr_ch = params.srr_list ? Channel
        .fromPath(params.srr_list)
        .splitText()
        .map { it.trim() }
        .filter { it && !it.startsWith('#') } : null

    /*
     * Optionally run download_reads
     */
    reads_ch = srr_ch ? download_reads(srr_ch) : null

    /*
     * Optionally run shovill_assembly
     * If reads_ch is null, user must provide --reads directly for primersearch
     */
    assemblies_ch = reads_ch ? shovill_assembly(reads_ch) : null

    /*
     * Run primersearch pipeline
     * Input can be from shovill output or directly from folder
     */
    primersearch_input_ch = assemblies_ch ?: Channel.fromPath("${params.reads}/**/*.{fasta,fa}")
                                                .map { file -> tuple(file.getBaseName(), file) }

    amplicon_fastas = primersearch_workflow(primersearch_input_ch)

    /*
     * Pairwise difference matrix (FINAL STEP)
     */
    pairwise_matrix(amplicon_fastas.amplicon_ch)
}
