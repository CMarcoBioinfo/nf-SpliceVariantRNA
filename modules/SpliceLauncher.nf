process SPLICELAUNCHER_INSTALL {
    tag "$genome_name"

    input:
    path gff3
    path mane
    path fasta
    val genome_name

    output:
    path "${genome_name}"                              , emit: ref_dir
    path "${genome_name}/BEDannotation.bed"            , emit: bed
    path "${genome_name}/SJDBannotation.sjdb"          , emit: sjdb
    path "${genome_name}/SpliceLauncherAnnot.txt"      , emit: annot
    path "${genome_name}/STARgenome/*.tab"             , optional: true, emit: tabs
    path "${genome_name}/STARgenome/Genome"            , optional: true, emit: star_genome
    path "${genome_name}.splicelauncher_install.out"   , emit: log_out
    path "${genome_name}.splicelauncher_install.err"   , emit: log_err

    script:
    def outdir_abs = file(params.outdir).toAbsolutePath()
    """
    bash ${params.splicelauncher} --runMode INSTALL \\
        --output ${genome_name} \\
        --gff ${gff3} \\
        --fasta ${fasta} \\
        --mane ${mane} \\
        --STAR ${params.star} \\
        --samtools ${params.samtools} \\
        --bedtools ${params.bedtools} \\
        -t ${task.cpus} \\
        1> ${genome_name}.splicelauncher_install.out \\
        2> ${genome_name}.splicelauncher_install.err

    mkdir -p ${outdir_abs}/references
    rsync -ac ${genome_name} ${outdir_abs}/references/

    CONFIG_FILE="\$(dirname "${params.splicelauncher}")/config.cfg"
    if [ -f "\$CONFIG_FILE" ]; then
        sed -i "s|^genome=.*|genome=\"${outdir_abs}/references/${genome_name}/STARgenome\"|" "\$CONFIG_FILE"
        sed -i "s|^SJDBannot=.*|SJDBannot=\"${outdir_abs}/references/${genome_name}/SJDBannotation.sjdb\"|" "\$CONFIG_FILE"
        sed -i "s|^spliceLaucherAnnot=.*|spliceLaucherAnnot=\"${outdir_abs}/references/${genome_name}/SpliceLauncherAnnot.txt\"|" "\$CONFIG_FILE"
        sed -i "s|^BEDrefPath=.*|BEDrefPath=\"${outdir_abs}/references/${genome_name}/BEDannotation.bed\"|" "\$CONFIG_FILE"
    fi
    """
}

process SPLICELAUNCHER_ALIGN {
    tag "${id}.${params.min_length}bp ($group_id)"
    label 'process_high'

    input:
    tuple val(id), val(group_id), path(reads)
    path genome_dir

    output:
    tuple val(id), val(group_id), path("align_out/Bam/*.bam")             , emit: bam
    tuple val(id), val(group_id), path("align_out/Bam/*.bam.bai")         , optional: true, emit: bai
    tuple val(id), val(group_id), path("align_out/Bam/*.bam.csi")         , optional: true, emit: csi
    tuple val(id), val(group_id), path("align_out/Bam/*.SJ.out.tab")      , optional: true, emit: sj_tab
    tuple val(id), val(group_id), path("align_out/Bam/*.Log.final.out")   , optional: true, emit: log_final
    path "${id}.${params.min_length}bp.splicelauncher_align.out"          , emit: log_out
    path "${id}.${params.min_length}bp.splicelauncher_align.err"          , emit: log_err

    script:
    def outdir_abs = file(params.outdir).toAbsolutePath()
    """
    mkdir -p fastq_input align_out align_tmp
    ln -s \$(readlink -f ${reads[0]}) fastq_input/${id}.${params.min_length}bp_R1_001.fastq.gz
    ln -s \$(readlink -f ${reads[1]}) fastq_input/${id}.${params.min_length}bp_R2_001.fastq.gz
    
    STAR_GENOME_PATH=\$(readlink -f ${genome_dir}/STARgenome)

    bash ${params.splicelauncher} --runMode Align \\
        --fastq fastq_input \\
        --output align_out \\
        -p \\
        --threads ${task.cpus} \\
        --tmpDir align_tmp \\
        --genome \$STAR_GENOME_PATH \\
        --STAR ${params.star} \\
        --samtools ${params.samtools} \\
        1> ${id}.${params.min_length}bp.splicelauncher_align.out \\
        2> ${id}.${params.min_length}bp.splicelauncher_align.err

    for bam in align_out/Bam/*.bam; do
        if [ -f "\$bam" ]; then
            if [ ! -f "\${bam}.bai" ]; then
                ${params.samtools} index -@ ${task.cpus} -b "\$bam"
            fi
        fi
    done

    mkdir -p ${outdir_abs}/splicelauncher/mapping/${group_id}
    rsync -ac align_out/Bam/ ${outdir_abs}/splicelauncher/mapping/${group_id}/

    rm -rf align_tmp fastq_input
    """
}

process SPLICELAUNCHER_COUNT {
    tag "$run_id"

    input:
    path bams
    path bed_annot
    val run_id

    output:
    tuple val(run_id), path("count_out/${run_id}.txt")            , emit: count_matrix
    tuple val(run_id), path("count_out/getClosestExons/*.count")  , optional: true, emit: counts
    path "${run_id}.splicelauncher_count.out"                     , emit: log_out
    path "${run_id}.splicelauncher_count.err"                     , emit: log_err

    script:
    def outdir_abs = file(params.outdir).toAbsolutePath()
    def paired_end = params.paired_end ? "-p" : ""

    """
    mkdir -p bam_input count_out
    for b in ${bams}; do
        ln -s "\$(readlink -f "\$b")" bam_input/
    done

    bash ${params.splicelauncher} --runMode Count \\
        --bam bam_input \\
        --output count_out \\
        --BEDannot ${bed_annot} \\
        ${paired_end} \\
        --bedtools ${params.bedtools} \\
        --samtools ${params.samtools} \\
        -p
        1> ${run_id}.splicelauncher_count.out \\
        2> ${run_id}.splicelauncher_count.err

    mv count_out/*.txt count_out/${run_id}.txt 2>/dev/null || true

    mkdir -p ${outdir_abs}/splicelauncher/${run_id}/sample_counts
    mkdir -p ${outdir_abs}/splicelauncher/${run_id}/count_matrix

    rsync -ac count_out/getClosestExons/ ${outdir_abs}/splicelauncher/${run_id}/sample_counts/ 2>/dev/null || true
    rsync -ac count_out/${run_id}.txt ${outdir_abs}/splicelauncher/${run_id}/count_matrix/ 2>/dev/null || true

    rm -rf bam_input
    """
}

process SPLICELAUNCHER_ANALYSIS {
    tag "$run_id"

    input:
        tuple val(run_id), path(count_matrix)
        path annot_txt

    output:
        tuple val(run_id), path("${run_id}_results/*") , emit: results_files
        path "${run_id}.splicelauncher_analysis.out"   , emit: log_out
        path "${run_id}.splicelauncher_analysis.err"   , emit: log_err

    script:
    def outdir_abs   = file(params.outdir).toAbsolutePath()
    def graphics_opt = params.graphics ? "--Graphics" : ""
    def txt_opt      = params.txt_out ? "--txtOut" : ""
    def bed_opt      = params.bed_out ? "--bedOut" : ""

    """
    mkdir -p ${run_id}_results

    bash << 'EOF'

    header=\$(head -n1 ${count_matrix})
    samples=\$(echo "\$header" | cut -f6-)

    sample_names=""
    for col in \$samples; do
        clean=\${col%.Aligned.sortedByCoord.out.count}

        if [[ \$clean =~ ^[0-9] ]]; then
            clean="X\$clean"
        fi

        if [[ -z "\$sample_names" ]]; then
            sample_names="\$clean"
        else
            sample_names="\$sample_names|\$clean"
        fi
    done

    echo "SampleNames utilisés : \$sample_names"

    bash ${params.splicelauncher} --runMode SpliceLauncher \
        -I ${count_matrix} \
        -O ${run_id}_results \
        -R ${annot_txt} \
        --SampleNames "\$sample_names" \
        --min_cov ${params.min_cov} \
        --threshold ${params.threshold} \
        ${graphics_opt} \
        ${txt_opt} \
        ${bed_opt} \
        1> ${run_id}.splicelauncher_analysis.out \
        2> ${run_id}.splicelauncher_analysis.err

    EOF

    mkdir -p ${outdir_abs}/splicelauncher/${run_id}/analysis_results
    rsync -ac ${run_id}_results/ ${outdir_abs}/splicelauncher/${run_id}/analysis_results/
    """
}

