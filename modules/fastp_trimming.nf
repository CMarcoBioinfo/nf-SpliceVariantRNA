process FASTP_TRIMMING {
    tag "$id ($group_id)"

    input:
    tuple val(id), path(r1), path(r2), val(group_id), val(tech)

    output:
    tuple val(id), val(group_id), path("${id}.${params.min_length}bp.1.fastq.gz"), path("${id}.${params.min_length}bp.2.fastq.gz"), path("${id}.${params.min_length}bp.html"), path("${id}.${params.min_length}bp.json"), path("${id}.${params.min_length}bp.out"), path("${id}.${params.min_length}bp.err"), emit: fastp_raw

    script:
    def outdir_abs = file(params.outdir).toAbsolutePath()
    """
    ${params.fastp} \\
        --thread ${task.cpus} \\
        --length_required ${params.min_length} \\
        --qualified_quality_phred ${params.mean_quality} \\
        --detect_adapter_for_pe \\
        --in1 ${r1} \\
        --in2 ${r2} \\
        --out1 ${id}.${params.min_length}bp.1.fastq.gz \\
        --out2 ${id}.${params.min_length}bp.2.fastq.gz \\
        --html ${id}.${params.min_length}bp.html \\
        --json ${id}.${params.min_length}bp.json \\
        1> ${id}.${params.min_length}bp.out \\
        2> ${id}.${params.min_length}bp.err

    # Copie robuste des FASTQ trimés vers results/fastq_trimmed/<groupe>/
    mkdir -p ${outdir_abs}/fastq_trimmed/${group_id}
    rsync -ac ${id}.${params.min_length}bp.1.fastq.gz ${id}.${params.min_length}bp.2.fastq.gz ${outdir_abs}/fastq_trimmed/${group_id}/

    # Copie robuste des rapports fastp vers results/qc/fastp/<groupe>/
    mkdir -p ${outdir_abs}/qc/fastp/${group_id}
    rsync -ac ${id}.${params.min_length}bp.html ${id}.${params.min_length}bp.json ${outdir_abs}/qc/fastp/${group_id}/
    """
}