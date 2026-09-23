nextflow.enable.dsl = 2

include { FASTP_TRIMMING }                                                         from './modules/fastp_trimming.nf'
include { FASTQC as FASTQC_RAW; FASTQC as FASTQC_TRIMMED }                         from './modules/fastqc.nf'
include { SPLICELAUNCHER_INSTALL; SPLICELAUNCHER_ALIGN; SPLICELAUNCHER_COUNT; SPLICELAUNCHER_ANALYSIS } from './modules/SpliceLauncher.nf'

workflow {

    main:
    // =============================================================
    // 0. IDENTIFIANT UNIQUE DE RUN
    // =============================================================
    def run_date         = new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())
    def samplesheet_name = file(params.samplesheet).baseName
    def unique_run_id    = params.run_id ?: "${samplesheet_name}_${run_date}"

    println "--> [INFO] Analyse de la série : ${unique_run_id}"

    // =============================================================
    // 1. RÉFÉRENCE
    // =============================================================
    def outdir_abs   = file(params.outdir).toAbsolutePath()
    def ref_dir_path = file("${outdir_abs}/references/${params.genome_name}")
    def ref_bed      = file("${outdir_abs}/references/${params.genome_name}/BEDannotation.bed")
    def ref_sjdb     = file("${outdir_abs}/references/${params.genome_name}/SJDBannotation.sjdb")
    def ref_annot    = file("${outdir_abs}/references/${params.genome_name}/SpliceLauncherAnnot.txt")

    if (ref_bed.exists() && ref_bed.size() > 0 && ref_sjdb.exists() && ref_annot.exists()) {
        println "--> [SKIP] Référence SpliceLauncher '${params.genome_name}' déjà existante."
        db_ref_dir   = Channel.value(ref_dir_path)
        db_ref_bed   = Channel.value(ref_bed)
        db_ref_sjdb  = Channel.value(ref_sjdb)
        db_ref_annot = Channel.value(ref_annot)
    } else {
        println "--> [RUN]  Installation SpliceLauncher..."
        gff3_file  = file(params.gff3, checkIfExists: true)
        mane_file  = file(params.mane, checkIfExists: true)
        fasta_file = file(params.fasta, checkIfExists: true)

        SPLICELAUNCHER_INSTALL(gff3_file, mane_file, fasta_file, params.genome_name)

        db_ref_dir   = SPLICELAUNCHER_INSTALL.out.ref_dir
        db_ref_bed   = SPLICELAUNCHER_INSTALL.out.bed
        db_ref_sjdb  = SPLICELAUNCHER_INSTALL.out.sjdb
        db_ref_annot = SPLICELAUNCHER_INSTALL.out.annot
    }

    // =============================================================
    // 2. INPUT
    // =============================================================
    raw_samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.id,
                file(row.path_read1),
                file(row.path_read2),
                row.group_id,
                row.technologie
            )
        }

    // =============================================================
    // 3. FASTQC RAW
    // =============================================================
    raw_samples_ch
        .branch { id, r1, r2, group, tech ->
            def base_r1 = r1.name.replaceAll(/(\.gz|\.fastq|\.fq)+$/, '')
            def base_r2 = r2.name.replaceAll(/(\.gz|\.fastq|\.fq)+$/, '')
            def out_zip1 = file("${outdir_abs}/qc/fastqc_raw/${group}/${base_r1}_fastqc.zip")
            def out_zip2 = file("${outdir_abs}/qc/fastqc_raw/${group}/${base_r2}_fastqc.zip")

            already_done: out_zip1.exists() && out_zip1.size() > 0 && out_zip2.exists() && out_zip2.size() > 0
            to_process:   true
                return tuple(id, group, [r1, r2], 'raw')
        }
        .set { fastqc_raw_branch }

    FASTQC_RAW(fastqc_raw_branch.to_process)

    // =============================================================
    // 4. FASTP TRIMMING
    // =============================================================
    raw_samples_ch
        .branch { id, r1, r2, group, tech ->
            def target_r1 = file("${outdir_abs}/fastq_trimmed/${group}/${id}.${params.min_length}bp.1.fastq.gz")
            def target_r2 = file("${outdir_abs}/fastq_trimmed/${group}/${id}.${params.min_length}bp.2.fastq.gz")

            already_done: target_r1.exists() && target_r1.size() > 0 && target_r2.exists() && target_r2.size() > 0
                return tuple(id, group, target_r1, target_r2)
            to_process:   true
                return tuple(id, r1, r2, group, tech)
        }
        .set { fastp_branch }

    FASTP_TRIMMING(fastp_branch.to_process)

    // =============================================================
    // 5. FASTQC TRIMMED
    // =============================================================
    all_trimmed_reads_ch = FASTP_TRIMMING.out.fastp_raw
        .map { id, group, r1, r2, html, json, out, err -> tuple(id, group, r1, r2) }
        .mix(fastp_branch.already_done)

    all_trimmed_reads_ch
        .branch { id, group, r1, r2 ->
            def base_r1 = r1.name.replaceAll(/(\.gz|\.fastq|\.fq)+$/, '')
            def base_r2 = r2.name.replaceAll(/(\.gz|\.fastq|\.fq)+$/, '')
            def out_zip1 = file("${outdir_abs}/qc/fastqc_trimmed/${group}/${base_r1}_fastqc.zip")
            def out_zip2 = file("${outdir_abs}/qc/fastqc_trimmed/${group}/${base_r2}_fastqc.zip")

            already_done: out_zip1.exists() && out_zip1.size() > 0 && out_zip2.exists() && out_zip2.size() > 0
            to_process:   true
                return tuple(id, group, [r1, r2], 'trimmed')
        }
        .set { fastqc_trimmed_branch }

    FASTQC_TRIMMED(fastqc_trimmed_branch.to_process)

    // =============================================================
    // 6. ALIGNEMENT
    // =============================================================
    all_trimmed_reads_ch
        .branch { id, group, r1, r2 ->
            def target_bam = file("${outdir_abs}/splicelauncher/mapping/${group}/${id}.${params.min_length}bp.Aligned.sortedByCoord.out.bam")

            already_done: target_bam.exists() && target_bam.size() > 0
                return tuple(id, group, target_bam)
            to_process:   true
                return tuple(id, group, [r1, r2])
        }
        .set { align_branch }

    SPLICELAUNCHER_ALIGN(align_branch.to_process, db_ref_dir.collect())

    ch_final_bams = SPLICELAUNCHER_ALIGN.out.bam
        .mix(align_branch.already_done)

    // =============================================================
    // 7. COUNT
    // =============================================================
    all_bams_collected = ch_final_bams
        .map { id, group, bam -> bam }
        .collect()

    SPLICELAUNCHER_COUNT(all_bams_collected, db_ref_bed, unique_run_id)

    // =============================================================
    // 8. ANALYSIS (AVEC CRÉATION DE SampleNames.txt)
    // =============================================================
    SPLICELAUNCHER_ANALYSIS(
        SPLICELAUNCHER_COUNT.out.count_matrix,
        db_ref_annot
    )

}


// output {
//     pipeline_outputs {
//         mode 'copy'
//         overwrite true

//         path { record ->
//             record.trimmed_r1        >> "fastq_trimmed/${record.group_id}/${record.trimmed_r1.name}"
//             record.trimmed_r2        >> "fastq_trimmed/${record.group_id}/${record.trimmed_r2.name}"
//             record.fastp_html        >> "qc/fastp/${record.group_id}/${record.fastp_html.name}"
//             record.fastp_json        >> "qc/fastp/${record.group_id}/${record.fastp_json.name}"
//             record.fastqc_raw_html1  >> "qc/fastqc_raw/${record.group_id}/${record.fastqc_raw_html1.name}"
//             record.fastqc_raw_html2  >> "qc/fastqc_raw/${record.group_id}/${record.fastqc_raw_html2.name}"
//             record.fastqc_raw_zip1   >> "qc/fastqc_raw/${record.group_id}/${record.fastqc_raw_zip1.name}"
//             record.fastqc_raw_zip2   >> "qc/fastqc_raw/${record.group_id}/${record.fastqc_raw_zip2.name}"
//             record.fastqc_trim_html1 >> "qc/fastqc_trimmed/${record.group_id}/${record.fastqc_trim_html1.name}"
//             record.fastqc_trim_html2 >> "qc/fastqc_trimmed/${record.group_id}/${record.fastqc_trim_html2.name}"
//             record.fastqc_trim_zip1  >> "qc/fastqc_trimmed/${record.group_id}/${record.fastqc_trim_zip1.name}"
//             record.fastqc_trim_zip2  >> "qc/fastqc_trimmed/${record.group_id}/${record.fastqc_trim_zip2.name}"
//         }

//         index {
//             path 'samplesheet_processed.csv'
//             header true
//         }
//     }
// }
