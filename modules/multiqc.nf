process MULTIQC {
    tag "$type"

    input:
    path(qc_files)
    val(type)
    val(report_name)

    output:
    path("${report_name}.html")     , emit: html
    path("${report_name}_data")     , emit: data
    tuple val(type), path("${report_name}.html"), path("${report_name}_data"), emit: multiqc_raw

    script:
    def report_title = type == 'raw' ? 'Quality Control of raw fastq files' : 'Quality Control of trimmed fastq files'
    """
    ${params.multiqc} \\
        . \\
        --filename ${report_name} \\
        --title "${report_title}" \\
        --dirs --dirs-depth 1
    """
}