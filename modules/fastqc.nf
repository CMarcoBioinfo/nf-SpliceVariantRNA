process FASTQC {
    tag "$id ($qc_type)"

    input:
    tuple val(id), val(group_id), path(reads), val(qc_type)

    output:
    tuple val(id), val(group_id), path("*_fastqc.html"), path("*_fastqc.zip"), val(qc_type), emit: qc_files
    path "${id}.${qc_type}.fastqc.out"                                                      , emit: log_out
    path "${id}.${qc_type}.fastqc.err"                                                      , emit: log_err

    script:
    def outdir_abs = file(params.outdir).toAbsolutePath()
    """
    ${params.fastqc} -t ${task.cpus} ${reads} \\
        1> ${id}.${qc_type}.fastqc.out \\
        2> ${id}.${qc_type}.fastqc.err

    # Copie directe et robuste avec rsync vers results/qc/fastqc_<raw|trimmed>/<groupe>/
    mkdir -p ${outdir_abs}/qc/fastqc_${qc_type}/${group_id}
    rsync -ac *_fastqc.html *_fastqc.zip ${outdir_abs}/qc/fastqc_${qc_type}/${group_id}/
    """
}