#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ======================================================
// PARAMETERS
// ======================================================

params.read1    = null
params.read2    = null
params.fasta    = null
params.sample_id = "sample"
params.outdir   = "results"

// ======================================================
// PROCESS: ILLUMINA FASTQ → EL_GATO
// ======================================================

process EL_GATO_FASTQ {

    container 'staphb/elgato:1.22.0'

    publishDir "${params.outdir}/illumina_fastq", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "${sample_id}"

    script:
    """
    el_gato.py \
      --read1 ${read1} \
      --read2 ${read2} \
      --out ${sample_id}
    """
}

// ======================================================
// PROCESS: ILLUMINA FASTA → EL_GATO
// ======================================================

process EL_GATO_FASTA {

    container 'staphb/elgato:1.22.0'

    publishDir "${params.outdir}/illumina_fasta", mode: 'copy', overwrite: true

    input:
    path fasta

    output:
    path "${fasta.baseName}"

    script:
    """
    el_gato.py \
      --read1 ${fasta} \
      --out ${fasta.baseName}
    """
}

// ======================================================
// ONT PIPELINE
// ======================================================

process FLYE {

    publishDir "${params.outdir}/ont/flye", mode: 'copy', overwrite: true

    input:
    path reads

    output:
    path "assembly.fasta"

    script:
    """
    flye --nano-raw ${reads} --out-dir flye_out --genome-size 5m --threads ${task.cpus} --asm-coverage 60
    cp flye_out/assembly.fasta .
    """
}

process MEDAKA {

    publishDir "${params.outdir}/ont/medaka", mode: 'copy', overwrite: true

    input:
    path assembly

    output:
    path "consensus.fasta"

    script:
    """
    medaka_consensus \
      -i ${assembly} \
      -d ${assembly} \
      -o medaka_out

    cp medaka_out/consensus.fasta .
    """
}

process EL_GATO_ONT {

    container 'staphb/elgato:1.22.0'

    publishDir "${params.outdir}/ont/el_gato", mode: 'copy', overwrite: true

    input:
    path fasta

    output:
    path "ont_result"

    script:
    """
    el_gato.py \
      --read1 ${fasta} \
      --out ont_result
    """
}

// ======================================================
// WORKFLOW (FIXED ROUTER)
// ======================================================

workflow {

    // ==================================================
    // MODE 1: FASTQ
    // ==================================================
    if( params.read1 && params.read2 ) {

        log.info "RUNNING MODE 1: ILLUMINA FASTQ"

        Channel.of(
            tuple(
                params.sample_id,
                file(params.read1),
                file(params.read2)
            )
        ).set { reads_ch }

        EL_GATO_FASTQ(reads_ch)
    }

    // ==================================================
    // MODE 2: FASTA
    // ==================================================
    else if( params.fasta ) {

        log.info "RUNNING MODE 2: ILLUMINA FASTA"

        Channel.of(file(params.fasta))
            .set { fasta_ch }

        EL_GATO_FASTA(fasta_ch)
    }

    // ==================================================
    // MODE 3: ONT
    // ==================================================
     else if( params.read1 ) {

        log.info "RUNNING MODE 3: ONT PIPELINE (SINGLE FASTQ)"

        Channel.of(file(params.read1))
            .set { ont_ch }

        FLYE(ont_ch)
        MEDAKA(FLYE.out)
        EL_GATO_ONT(MEDAKA.out)
    }
}