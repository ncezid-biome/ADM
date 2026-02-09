nextflow.enable.dsl = 2

/*
 * Module: shovill_assembly
 *
 * Assemble paired-end FASTQ reads using Shovill.
 *
 * Inputs:
 *   Channel of tuples: (sample_id, R1.fastq, R2.fastq)
 * Outputs:
 *   Emits assembled contigs: tuple(sample_id, assembled_fasta)
 *   Assembled FASTA files are published in: ${params.outdir}/assemblies/
 */

workflow shovill_assembly {

    take:
        paired_reads  // tuple(sample_id, R1, R2)

    main:
        assembled_ch = SHOVILL_ASSEMBLY(paired_reads)

    emit:
        assembled_ch
}


process SHOVILL_ASSEMBLY {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(reads1), path(reads2)

    output:
        tuple val(sample_id), path("${sample_id}_assembled.fasta")

    script:
    """
    set -euo pipefail

    shovill -R1 ${reads1} -R2 ${reads2} \
            --outdir ${sample_id} --force --trim ON --cpu ${task.cpus} > /dev/null 2>&1

    mv ${sample_id}/contigs.fa ${sample_id}_assembled.fasta
    """
}
