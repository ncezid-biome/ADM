nextflow.enable.dsl = 2
import java.nio.file.Paths

/*
 * Workflow: primersearch_workflow
 *
 * Combines three processes:
 * 1) Run primersearch (EMBOSS)
 * 2) Parse primersearch results (Python)
 * 3) Convert extracted FASTA to JSON (Python)
 *
 * Inputs:
 *   Channel of tuples: (sample_id, assembly.fasta)
 * Outputs:
 *   Emits parsed amplicons and JSON files
 * Published in: ${params.outdir}/primersearch/
 */

workflow primersearch_workflow {

    take:
        assemblies_ch// optional: true

    main:
        // fallback to folder if channel not provided
        fasta_ch = assemblies_ch ?: Channel.fromPath(["${params.reads}/*_{assembled,contigs}.fasta", "${params.reads}/**/*_{assembled,contigs}.fasta"])
            .map { file ->
                def sample = file.baseName.replaceAll(/_assembled|_contigs/, '')
                tuple(sample, file)
            }

        // 1) Run primersearch
        primersearch_ch = RUN_PRIMERSEARCH(fasta_ch)

        // 2) Parse primersearch output
        amplicon = PARSE_PRIMERSEARCH(primersearch_ch)

        // 3) Convert FASTA to JSON
        json_ch = FASTA_TO_JSON(amplicon.amplicons)

    emit:
        amplicon_ch = amplicon.amplicons
        // amplicon_folder_ch = amplicon.amplicon_folder
        json_ch

}

////////////////////////////////////////////
// PROCESS: RUN_PRIMERSEARCH
////////////////////////////////////////////
process RUN_PRIMERSEARCH {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(fasta_file)

    output:
        tuple val(sample_id), path("${sample_id}.ps"), path(fasta_file)

    script:
    """
    set -euo pipefail

    ${params.EMBOSS_PRIMERSEARCH_6_4_0} \
        -infile ${params.primers} \
                    -seqall ${fasta_file} \
                    -outfile ${sample_id}.ps \
                    -mismatchpercent ${params.mismatchpercent}
    """
}

////////////////////////////////////////////
// PROCESS: PARSE_PRIMERSEARCH
////////////////////////////////////////////
process PARSE_PRIMERSEARCH {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(ps_file), path(fasta_file)

    output:
        tuple val(sample_id), path("*_extractedAmplicons.fasta"), emit: amplicons
        path("*.txt"), emit: empty_primers, optional: true
        // path("${params.outdir}/primersearch/amplicon/*"), emit: amplicon_folder

    script:
    """
    set -euo pipefail

    parse_primersearch.py --results ${ps_file} \
                          --sequence ${fasta_file} \
                          --amp_len ${params.max_amplicon_len}
    """
}

////////////////////////////////////////////
// PROCESS: FASTA_TO_JSON
////////////////////////////////////////////
process FASTA_TO_JSON {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(amplicon_fasta)

    output:
        path("*.json"), emit: json_file

    script:
    """
    set -euo pipefail

    fasta_to_json.py --fasta_file ${amplicon_fasta} \
                     --sample_id ${sample_id} \
                     --primers ${params.primers}
    """
}
