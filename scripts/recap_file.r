#!/usr/bin/env Rscript

# Installs and loads a package if not already available.
check_install_load_package <- function(package_name, source = "CRAN") {
    if (!requireNamespace(package_name, quietly = TRUE)) {
        if (source == "Bioconductor") {
            BiocManager::install(package_name)
        } else {
            install.packages(package_name)
        }
    }
    library(package_name, character.only = TRUE)
}

# Ensures BiocManager is installed for Bioconductor dependencies.
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}

# Loads required packages.
check_install_load_package("stringr")
check_install_load_package("Rsamtools")

# Verifies NM transcript identifiers; updates them using MANE table if necessary.
NM_check <- function(NM, mane) {
    original <- NM
    NM_clean <- sub("\\..*$", "", original)
    if (grepl("^NM_", NM_clean)) {
        return(original)
    }

    mane_other_clean <- sub("\\..*$", "", mane$Other)
    match_index <- which(mane_other_clean == NM_clean)
    if (length(match_index) > 0) {
        nm_candidate <- mane$NM_refseq[match_index[1]]
        if (grepl("^NM_", nm_candidate)) {
            return(paste0(original, "; ", nm_candidate))
        }
    }
    return(original)
}

# Detects if the junction annotation represents a deletion (∆ but not ▼).
is_del <- function(AnnotJuncs){
    DeltaSymb = intToUtf8(0x2206)  # Symbole "∆"
    BlackInvTriangle = intToUtf8(0x25BC)  # Symbole "▼"
    
    return(grepl(DeltaSymb, AnnotJuncs) && !grepl(BlackInvTriangle, AnnotJuncs))
}

# Detects if the junction annotation represents an insertion (▼ but not ∆).
is_ins <- function(AnnotJuncs){
    BlackInvTriangle = intToUtf8(0x25BC)  # Symbole "▼"
    DeltaSymb = intToUtf8(0x2206)  # Symbole "∆"
    
    return(grepl(BlackInvTriangle, AnnotJuncs) && !grepl(DeltaSymb, AnnotJuncs))
}

# Builds an HGVS-compliant annotation string based on input coordinates and type of event (insertion/deletion).
generate_hgvs <- function(AnnotJuncs, NM, chr, start, end, strand, cStart, cEnd, fasta){
    parts <- unlist(strsplit(NM, ";\\s*"))
    NM_clean <- if (any(grepl("^NM_", parts))) {
    parts[grepl("^NM_", parts)][1]
    } else {
        parts[1]
    }

    rStart_clean <- sub("^c\\.", "", cStart)
    rEnd_clean   <- sub("^c\\.", "", cEnd)

    if (is_del(AnnotJuncs)) {
        return(paste0(NM_clean, ":r.", rStart_clean, "_", rEnd_clean, "del"))
    }   

    if (is_ins(AnnotJuncs)) {
        coords <- GRanges(chr, IRanges(start, end))
        seq <- getSeq(fasta, coords)

        if (strand == "-") {
            seq <- reverseComplement(seq)
        }

        seq_rna <- gsub("T", "U", as.character(seq))
        return(paste0(NM_clean, ":r.", rStart_clean, "_", rEnd_clean, "ins", as.character(seq_rna)))
    }

    return("Complex annotation")
}

# Reads input data from a file, supporting text and Excel formats.
read_file <- function(file, xlsx) {
    if(xlsx) {
        data_genes <- read_excel(file)
        data_genes <- as.data.frame(data_genes)
        return(data_genes)
    } else {
        data_genes <- read.csv(file, sep = "\t" , header = TRUE, check.names = FALSE)
        data_genes <- as.data.frame(data_genes)
        return(data_genes)
    }
}

# Writes two output data frames (statistical/non_statistical junctions) into text or Excel format.
write_file <- function(df_statistical, df_no_model, df_non_statistical, df_threshold, df_too_complex, output, xlsx) {
    if (xlsx) {
        wb <- createWorkbook()
        addWorksheet(wb, "Statistical Junctions")
        addWorksheet(wb, "Unique Junctions")
        addWorksheet(wb, "Threshold Exceeded Junctions")
        addWorksheet(wb, "No Model Junctions")
        addWorksheet(wb, "Event too complex")


        writeData(wb, "Statistical Junctions", df_statistical, rowNames = FALSE)
        writeData(wb, "Unique Junctions", df_non_statistical, rowNames = FALSE)
        writeData(wb, "Threshold Exceeded Junctions", df_threshold, rowNames = FALSE)
        writeData(wb, "No Model Junctions", df_no_model, rowNames = FALSE)
        writeData(wb, "Event too complex", df_too_complex, rowNames = FALSE)

        saveWorkbook(wb, file = output, overwrite = TRUE)

    } else {
        writeLines("##### Statistical Junctions #####", output)
        write.table(df_statistical, output, sep = ";", quote = FALSE, row.names = FALSE, append = TRUE)

        writeLines("\n##### Non Statistical Junctions - Unique #####", output, append = TRUE)
        write.table(df_non_statistical, output, sep = ";", quote = FALSE, row.names = FALSE, append = TRUE)

        writeLines("\n##### Non Statistical Junctions - Threshold exceeded #####", output, append = TRUE)
        write.table(df_threshold, output, sep = ";", quote = FALSE, row.names = FALSE, append = TRUE)

        writeLines("\n##### Non Statistical Junctions - No model #####", output, append = TRUE)
        write.table(df_no_model, output, sep = ";", quote = FALSE, row.names = FALSE, append = TRUE)

        writeLines("\n##### Non Statistical Junctions - Event too complex #####", output, append = TRUE)
        write.table(df_too_complex, output, sep = ";", quote = FALSE, row.names = FALSE, append = TRUE)
    }
}


# Parses command-line arguments for dynamic script configuration.
args <- commandArgs(trailingOnly = TRUE)
extension <- if (any(args %in% "--text")) ".txt" else ".xlsx"
bool_xlsx <- if (any(args %in% "--text")) FALSE else TRUE
sample <- NULL
statistical_file <- NULL
non_statistical_file <- NULL
output_file <- NULL
fasta <- NULL
mane <- NULL


for (i in seq_along(args)) {
    if (args[i] == "--statisticalFile"){
        statistical_file <- args[i + 1]
    } else if (args[i] == "--nonStatisticalFile") {
        non_statistical_file <- args[i + 1]
    } else if (args[i] == "--reference") {
        fasta <- args[i + 1]
    } else if (args[i] == "--output") {
        output_file<- args[i + 1]
    } else if (args[i] == "--sample") {
        sample <- args[i + 1]
    } else if (args[i] == "--mane") {
        mane <- args[i + 1]
    }
}


# Stops execution if the input file is missing or invalid.
if (is.null(statistical_file) || !file.exists(statistical_file)) {
    stop("Statistical file is missing or invalid.")
}

# Stops execution if the input file is missing or invalid.
if (is.null(non_statistical_file) || !file.exists(non_statistical_file)) {
    stop("Non Statistical file is missing or invalid.")
}

if (bool_xlsx) {
    check_install_load_package("readxl")
    check_install_load_package("openxlsx")
}


fasta_file <- FaFile(fasta)
open(fasta_file)

mane_file <- mane
mane <- read.delim(mane_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(mane) <- c("Gene", "NM_refseq", "Other", "Strand")


statistical_genes <- NULL
statistical_genes <- read_file(statistical_file, bool_xlsx)
P_sample <- paste0("P_",sample)
col_statistical <- c("Conca", "chr", "start", "end", "strand", "NM", "Gene", sample, P_sample, "event_type", "AnnotJuncs", "cStart", "cEnd", "DistribAjust", "Significative", "filterInterpretation", "nbSignificantSamples", "p_value", "SignificanceLevel")
statistical_filtered <- statistical_genes[, col_statistical]
statistical_filtered$NM <- sapply(statistical_filtered$NM, NM_check, mane = mane)
statistical_filtered$HGVS <- mapply(generate_hgvs,
    AnnotJuncs = statistical_filtered$AnnotJuncs,
    NM = statistical_filtered$NM,
    chr = statistical_filtered$chr,
    start = statistical_filtered$start,
    end = statistical_filtered$end,
    strand = statistical_filtered$strand,
    cStart = statistical_filtered$cStart,
    cEnd = statistical_filtered$cEnd,
    MoreArgs = list(fasta = fasta_file)
)
cols <- colnames(statistical_filtered)
i <- which(cols == "HGVS")
j <- which(cols == "SignificanceLevel")
cols[c(i, j)] <- cols[c(j, i)]
statistical_filtered <- statistical_filtered[, cols]


non_statistical_genes <- NULL
non_statistical_genes <- read_file(non_statistical_file, bool_xlsx)
col_non_statistical <- c("Conca", "chr", "start", "end", "strand", "NM", "Gene", sample, P_sample, "event_type", "AnnotJuncs", "cStart", "cEnd", "SampleReads", "nbSampFilter","filterInterpretation")
non_statistical_filtered <- non_statistical_genes[, col_non_statistical]
non_statistical_filtered$NM <- sapply(non_statistical_filtered$NM, NM_check, mane = mane)
non_statistical_filtered$HGVS <- mapply(generate_hgvs,
    AnnotJuncs = non_statistical_filtered$AnnotJuncs,
    NM = non_statistical_filtered$NM,
    chr = non_statistical_filtered$chr,
    start = non_statistical_filtered$start,
    end = non_statistical_filtered$end,
    strand = non_statistical_filtered$strand,
    cStart = non_statistical_filtered$cStart,
    cEnd = non_statistical_filtered$cEnd,
    MoreArgs = list(fasta = fasta_file)
)
cols <- colnames(non_statistical_filtered)
i <- which(cols == "HGVS")
j <- which(cols == "filterInterpretation")
cols[c(i, j)] <- cols[c(j, i)]
non_statistical_filtered <- non_statistical_filtered[, cols]
non_statistical_no_model <- non_statistical_filtered[non_statistical_filtered$filterInterpretation == "No model", ]
non_statistical_non_statistical <- non_statistical_filtered[non_statistical_filtered$filterInterpretation == "Unique junction", ]
non_statistical_threshold <- non_statistical_filtered[non_statistical_filtered$filterInterpretation == "Percentage threshold execeeded", ]
non_statistical_too_complex <-non_statistical_filtered[non_statistical_filtered$filterInterpretation == "Event too complex", ]

write_file(statistical_filtered, non_statistical_no_model, non_statistical_non_statistical, non_statistical_threshold, non_statistical_too_complex, output_file, bool_xlsx)
