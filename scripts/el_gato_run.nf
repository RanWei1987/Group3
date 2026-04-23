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

    cpus 2
    memory '6 GB'

    publishDir "${params.outdir}/ont/flye", mode: 'copy', overwrite: true

    input:
    path reads

   output:
   tuple val(reads), path("assembly.fasta"), path("flye_out")

    script:
    """
    flye \
      --nano-raw ${reads} \
      --out-dir flye_out \
      --genome-size 3m \
      --threads ${task.cpus} \
      --asm-coverage 10

    cp flye_out/assembly.fasta .
    """
}

process SUBSAMPLE_ONT {

        cpus 1
    memory '1 GB'

    input:
    path reads

    output:
    path "subsampled.fastq.gz"

    script:
    """
    ls -la ${reads}
    seqtk sample -s100 ${reads} 5000 | gzip > subsampled.fastq.gz
    """
}

process RAVEN {

    cpus 2
    memory '4 GB'   // MUCH lower than Flye

    input:
    path reads

    output:
    path "assembly.fasta"

    container 'staphb/raven:latest'

    script:
    """
    ls -lh
    raven -t 1 $reads > assembly.fasta
    """
}

process MEDAKA {

    cpus 2
    memory '6 GB'

    input:
    tuple path(reads), path(assembly)

    output:
    path "consensus.fasta"

    script:
    """
    medaka_consensus \
      -i ${reads} \
      -d ${assembly} \
      -o medaka_out \
      -t ${task.cpus}

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
      --assembly ${fasta} \
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

    // Single source channel
    raw_reads = Channel.of(file(params.read1))

    // STEP 1: subsample
    subsampled_reads = SUBSAMPLE_ONT(raw_reads)

    // STEP 2: assemble
    raven_out = RAVEN(subsampled_reads)

    // STEP 3: pair reads + assembly
    medaka_input = raw_reads.combine(raven_out)

    // STEP 4: polish
    medaka_out = MEDAKA(medaka_input)

    // STEP 5: final tool
    EL_GATO_ONT(medaka_out)
}
}