nextflow.enable.dsl = 2

/*
 * Module: download_reads
 *
 * Download raw reads from NCBI SRA using prefetch + fasterq-dump and compress with pigz.
 *
 * Inputs:
 *   Channel of SRR accessions (optional, fallback to params.srr_list)
 * Outputs:
 *   Emits paired FASTQ files into: ${params.outdir}/fastqs/
 *   Channel emits: tuple( sample_id, [R1.fastq.gz, R2.fastq.gz] )
 */

workflow download_reads {

    take:
        srr_ch

    main:
        downloaded_reads = FASTQ_SRR(srr_ch)

    emit:
        downloaded_reads
}

process FASTQ_SRR {

    tag "$srr_id"

    input:
        val srr_id

    output:
        tuple val(srr_id), path("${srr_id}_1.fastq.gz"), path("${srr_id}_2.fastq.gz")

    script:
    """
    set -euo pipefail

    # Download .sra file
    prefetch ${srr_id} --output-directory .

    # Convert SRA to paired FASTQ
    fasterq-dump ${srr_id} --split-files --threads ${task.cpus} --outdir .

    # Compress FASTQ files
    gzip ${srr_id}_*.fastq

    # Optional: remove .sra to save space
    rm -rf ${srr_id}
    """
}
