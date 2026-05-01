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
params.checkm2_db  = "/mnt/checkm2_db/uniref100.KO.1.dmnd"
params.dashboard   = "${projectDir}/dashboard.py"

log.info "projectDir = ${projectDir}"

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

    container 'quay.io/biocontainers/raven-assembler:1.8.3--h43eeafb_0'
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

    container 'nanozoo/checkm:latest'
    publishDir "${params.outdir}/ont/checkm", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    path "checkm.tsv"

    script:
    """
    mkdir -p in
    cp ${fasta} in/

    checkm taxon_set genus Legionella legionella
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
    container 'nanozoo/checkm2:1.1.0--f5e828a'
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
        -database_path ${params.checkm2_db} \
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
// DASHBOARD BUILD AND LAUNCH
// ======================================================

process LAUNCH_DASHBOARD {
    cache    false   // always re-launch, never skip

    input:
    val elgato_path
    val amrfinder_path
    val checkm_path
    val nanoplot_path
    val consensus_path

    output:
    stdout

    script:
    def elgato_out    = "${launchDir}/${params.outdir}/ont/el_gato/report.json/report.json"
    def amrfinder_out = "${launchDir}/${params.outdir}/ont/amrfinder/${params.sample_id}_amr.tsv"
    def checkm_out    = "${launchDir}/${params.outdir}/ont/checkm/checkm.tsv"
    def nanoplot_out  = "${launchDir}/${params.outdir}/ont/nanoplot/nanoplot/NanoPlot-report.html"
    def consensus_out = "${launchDir}/${params.outdir}/ont/medaka/consensus.fasta"
    """
    ARGS=""
    [ -f ${elgato_out}    ] && ARGS="\$ARGS --elgato    ${elgato_out}"
    [ -f ${amrfinder_out} ] && ARGS="\$ARGS --amrfinder ${amrfinder_out}"
    [ -f ${checkm_out}    ] && ARGS="\$ARGS --checkm    ${checkm_out}"
    [ -f ${nanoplot_out}  ] && ARGS="\$ARGS --nanoplot  ${nanoplot_out}"
    [ -f ${consensus_out} ] && ARGS="\$ARGS --consensus ${consensus_out}"

    echo "ARGS built: \$ARGS"
    echo "Dashboard: ${params.dashboard}"

    # Run in background so Nextflow process can complete
    nohup streamlit run ${params.dashboard} -- \$ARGS \
        > ${launchDir}/streamlit.log 2>&1 &

    echo "Streamlit PID: \$!"
    echo "Dashboard launching at http://localhost:8501"
    echo "Log: ${launchDir}/streamlit.log"

    # Give it a moment to start then exit cleanly
    sleep 5
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
        CHECKM(polished)
        //CHECKM2(polished)
        PROKKA(polished)

        amr_input_ch = PROKKA.out.gff.join(PROKKA.out.faa)

    // Now amr_input_ch is [sample_id, gff, faa]
        AMRFINDER(amr_input_ch)

        EL_GATO_ONT(polished)

    // ifEmpty("MISSING") is a trick to make sure that each input produces one single value that goes to the dashboard - so even if onoe of these processes fails, all five inputs still exist
        // elgato_path    = EL_GATO_ONT.out
        //                     .map  { f -> f.toString() }
        //                     .ifEmpty( "MISSING" )

        // amrfinder_path = AMRFINDER.out
        //                     .map  { f -> f.toString() }
        //                     .ifEmpty( "MISSING" )

        // checkm_path    = CHECKM.out
        //                     .map  { f -> f.toString() }
        //                     .ifEmpty( "MISSING" )

        // nanoplot_path  = NANOPLOT.out
        //                     .map  { f -> f.toString() }
        //                     .ifEmpty( "MISSING" )

        // consensus_path = polished
        //                     .map  { sid, f -> f.toString() }
        //                     .ifEmpty( "MISSING" )
        // Replace the individual channel mappings with this
    Channel.empty()
        .mix(
            EL_GATO_ONT.out.map { f -> f.toString() },
            AMRFINDER.out.map { f -> f.toString() },
            CHECKM.out.map { f -> f.toString() },
            NANOPLOT.out.map { f -> f.toString() },
            polished.map { sid, f -> f.toString() }
        )
        .collect()
        .map { paths ->
            def findPath = { name -> paths.find { it.contains(name) } ?: "MISSING" }
            [
                findPath("report.json"),
                findPath("_amr.tsv"),
                findPath("checkm.tsv"),
                findPath("NanoPlot-report.html"),
                findPath("consensus.fasta")
            ]
        }
        .set { dashboard_inputs }

        LAUNCH_DASHBOARD(
            dashboard_inputs.map { it[0] },
            dashboard_inputs.map { it[1] },
            dashboard_inputs.map { it[2] },
            dashboard_inputs.map { it[3] },
            dashboard_inputs.map { it[4] }
        )
        // LAUNCH_DASHBOARD(
        //     elgato_path,
        //     amrfinder_path,
        //     checkm_path,
        //     nanoplot_path,
        //     consensus_path
        // )
    }

    // error message for if user messes up inputs
    else {
        error "No valid input provided. Use --read1 (ONT), --read1 + --read2 (Illumina PE), or --fasta."
    }
}