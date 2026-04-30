#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ======================================================
// PARAMETERS
// ======================================================

params.read1     = null
params.read2     = null
params.fasta     = null
params.sample_id = "sample"
params.outdir    = "results"

// ======================================================
// ILLUMINA MODE
// ======================================================

process EL_GATO_FASTQ {

    container 'staphb/elgato:1.22.0'
    publishDir "${params.outdir}/illumina_fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "${sample_id}"

    script:
    """
    el_gato.py --read1 ${read1} --read2 ${read2} --out ${sample_id}
    """
}

process EL_GATO_FASTA {

    container 'staphb/elgato:1.22.0'
    publishDir "${params.outdir}/illumina_fasta", mode: 'copy'

    input:
    path fasta

    output:
    path "${fasta.baseName}"

    script:
    """
    el_gato.py --read1 ${fasta} --out ${fasta.baseName}
    """
}

// ======================================================
// ONT PIPELINE (CLEAN LINEAR DAG)
// ======================================================

process SUBSAMPLE_ONT {

    container 'staphb/seqtk:latest'
    publishDir "${params.outdir}/ont/subsample", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("subsampled.fastq.gz")

    script:
    """
    seqtk sample -s100 ${reads} 5000 | gzip > subsampled.fastq.gz
    """
}

process RAVEN {

    container 'staphb/raven:latest'
    publishDir "${params.outdir}/ont/assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("assembly.fasta")

    script:
    """
    raven -t 1 ${reads} > assembly.fasta
    """
}

process MEDAKA {

    container 'staphb/medaka:latest'
    publishDir "${params.outdir}/ont/medaka", mode: 'copy'
    cpus 2

    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id), path("consensus.fasta")

    script:
    """
    medaka_consensus \
        -i ${reads} \
        -d ${assembly} \
        -o medaka_out \
        -t ${task.cpus} \
	-b 40

    cp medaka_out/consensus.fasta .
    """
}

// ======================================================
// ONT FAN-OUT TOOLS
// ======================================================

process NANOPLOT {
    container 'staphb/nanoplot:latest'
    publishDir "${params.outdir}/ont/nanoplot", mode: 'copy'
    cpus 2

    input:
    tuple val(sample_id), path(reads)

    output:
    path "nanoplot/NanoPlot-report.html"

    script:
    """
    NanoPlot --fastq ${reads} -o nanoplot --threads ${task.cpus}
    """
}

process CHECKM {

    container 'staphb/checkm:latest'
    publishDir "${params.outdir}/ont/checkm", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    path "checkm.tsv"

    script:
    """
    mkdir -p in
    cp ${fasta} in/

    checkm taxon_set genus Legionella legionella || true
    checkm analyze legionella in out -x fasta
    checkm qa legionella out -f checkm.tsv -o 2
    """
}

process PROKKA {
    container 'staphb/prokka:latest'
    publishDir "${params.outdir}/ont/prokka", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    val sample_id, emit: id
    // 给路径套上 sample_id，形成 tuple
    tuple val(sample_id), path("prokka/${sample_id}.gff"), emit: gff
    tuple val(sample_id), path("prokka/${sample_id}.faa"), emit: faa

    script:
    """
    prokka \
        --outdir prokka \
        --prefix ${sample_id} \
        --cpus ${task.cpus} \
        ${fasta}
    """
}

process CHECKM2 {
    container 'staphb/checkm2:latest'
    publishDir "${params.outdir}/ont/checkm2", mode: 'copy'
    cpus 4
    memory '4 GB'

    input:
    tuple val(sample_id), path(fasta)

    output:
    path "checkm2_out/quality_report.tsv"

    script:
    """
    checkm2 predict \
        --input ${fasta} \
        --output-directory checkm2_out \
        --threads ${task.cpus} \
        --database_path /mnt/checkm2_db/uniref100.KO.1.dmnd \
        --force
    """
}
process AMRFINDER {

    container 'staphb/ncbi-amrfinderplus:3.12.8' 
    publishDir "${params.outdir}/ont/amrfinder", mode: 'copy'

    input:
    tuple val(sample_id), path(gff), path(faa)

    output:
    path "${sample_id}_amr.tsv"

    script:
    """
    amrfinder \
        -g ${gff} \
        -p ${faa} \
        -a prokka \
        --plus \
        -o ${sample_id}_amr.tsv
    """
}
process EL_GATO_ONT {

    container 'staphb/elgato:1.22.0'
    publishDir "${params.outdir}/ont/el_gato", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    path "report.json"

    script:
    """
    el_gato.py --assembly ${fasta} --out report.json
    """
}

// ======================================================
// WORKFLOW ROUTER (SAFE)
// ======================================================

workflow {

    // =========================
    // ILLUMINA FASTQ
    // =========================
    if( params.read1 && params.read2 ) {

        Channel.of(tuple(params.sample_id, file(params.read1), file(params.read2)))
            .set { illumina_ch }

        EL_GATO_FASTQ(illumina_ch)
    }

    // =========================
    // ILLUMINA FASTA
    // =========================
    else if( params.fasta ) {

        Channel.of(file(params.fasta))
            .set { fasta_ch }

        EL_GATO_FASTA(fasta_ch)
    }

    // =========================
    // ONT PIPELINE (FIXED & SAFE)
    // =========================
    else if( params.read1 ) {

        log.info "RUNNING ONT PIPELINE (FINAL STABLE VERSION)"

        raw = Channel.of(tuple(params.sample_id, file(params.read1)))

        subsampled = SUBSAMPLE_ONT(raw)
        raven_out  = RAVEN(subsampled)

        medaka_in = raven_out.join(subsampled)

        polished = MEDAKA(medaka_in)
        NANOPLOT(subsampled)
        //CHECKM(polished)
        CHECKM2(polished)
        PROKKA(polished)

        amr_input_ch = PROKKA.out.gff.join(PROKKA.out.faa)

    // Now amr_input_ch is [sample_id, gff, faa]
        AMRFINDER(amr_input_ch)

        EL_GATO_ONT(polished)
    }
}