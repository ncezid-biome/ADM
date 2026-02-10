nextflow.enable.dsl = 2

include { download_reads } from './modules/download_reads'
include { shovill_assembly } from './modules/shovill_assembly'
include { primersearch_workflow } from './modules/primersearch_workflow'  // new module
include { pairwise_matrix }   from './modules/pairwise_matrix'

/*
    The Full workflow
*/
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



/*
    nextflow run main.nf -entry download_only \
    --srr_list srr_ids.txt \
    --outdir results
*/
workflow download_only {
    
    srr_ch = params.srr_list ? Channel.fromPath(params.srr_list)
        .splitText()
        .map { it.trim() }
        .filter { it && !it.startsWith('#') } : null

    if(!srr_ch) error "--srr_list is required for download_only"
    download_reads(srr_ch)
}

/*
    nextflow run main.nf -entry assembly_only \
    --reads /path/to/fastqs \
    --outdir results

*/
workflow assembly_only {

    // Take FASTQ files directly from folder (instead of downloading)
    reads_ch = Channel.fromFilePairs(["${params.reads}/*_{1,2}.fastq.gz", "${params.reads}/**/*_{1,2}.fastq.gz"], flat: true)

    shovill_assembly(reads_ch)
}


/*
    download + assembly

    nextflow run main.nf -entry da \
        --srr_list srr_ids.txt \
        --outdir results
*/
workflow da {

    srr_ch = params.srr_list ? Channel.fromPath(params.srr_list)
        .splitText()
        .map { it.trim() }
        .filter { it && !it.startsWith('#') } : null

    if(!srr_ch) error "--srr_list is required for download_only"

    reads_ch = download_reads(srr_ch)
    shovill_assembly(reads_ch)
}


/*
    nextflow run main.nf -entry primersearch_only \
    --reads assemblies \
    --primers primers.txt \
    --outdir results

*/
workflow primersearch_only {

    // No inputs needed — module will fall back to --reads
    primersearch_workflow(null)
}



/*
 * nextflow run main.nf -entry pairwise_only \
    --amplicon_folder $PWD/amplicons (need to use absolute path, or prefix it with $PWD \
    --primers primers.txt \
    --output_csv pairwise.csv \
    --outdir output_folder
 *
 */
workflow pairwise_only {

    Channel
        .fromPath(["${params.amplicon_folder}/*.{fa,fasta}", "${params.amplicon_folder}/**/*.{fa,fasta}"])
        .map { file -> tuple(file.getBaseName(), file) }
        | pairwise_matrix
}



/*
    primersearch + pairwise_matrix

    nextflow run main.nf -entry pp \
    --reads assemblies \
    --primers primers.txt \
    --outdir $PWD/results  (need to use absolute path, or prefix it with $PWD\
    --output_csv pairwise.csv
*/
workflow pp {

    ps = primersearch_workflow(null)
    pairwise_matrix(ps.amplicon_ch)
}
