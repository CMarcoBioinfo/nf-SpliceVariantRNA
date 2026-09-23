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

check_install_load_package("stringr")

# Parses command-line arguments for dynamic script configuration.
args <- commandArgs(trailingOnly = TRUE)
extension <- if (any(args %in% "--text")) ".txt" else ".xlsx"
bool_xlsx <- if (any(args %in% "--text")) FALSE else TRUE
input_file <- NULL
lenght_file <- NULL
graphics <- FALSE
StatAnalysis <- "YES"
thr=1

for (i in seq_along(args)) {
    if (args[i] == "--input"){
        input_file <- args[i + 1]
    } else if (args[i] == "--outputStatistical") {
        output_file_significant <- args[i + 1]
    } else if (args[i] == "--outputNonStatistical") {
        output_file_unique <- args[i + 1]
    } else if (args[i] == "--Graphics") {
        graphics <- TRUE
    } else if (args[i] == "--directory") {
        chem_dessin <- args[i + 1]
    } else if (args[i] == "--length") {
        length <- args[i + 1]
        length <- as.integer(length)
    } else if (args[i] == "--minNonStatisticalReads") {
        minUniqueReads <- args[i + 1]
        minUniqueReads <- as.integer(minUniqueReads)
    } else if (args[i] == "--maxNonStatisticalSamples") {
        maxUniqueSamples <- args[i + 1]
        maxUniqueSamples <- as.integer(maxUniqueSamples)
    } else if (args[i] == "--maxStatisticalSamples") {
        maxSignificantSamples <- args[i + 1]
        maxSignificantSamples <- as.integer(maxSignificantSamples)
    } else if (args[i] == "--thresholdSignificanceLevel") {
        thresholdSignificanceLevel <- args[i + 1]
        thresholdSignificanceLevel <- as.integer(thresholdSignificanceLevel)
    }
}

# Stops execution if the input file is missing or invalid.
if (is.null(input_file) || !file.exists(input_file)) {
    stop("Input file is missing or invalid.")
}

# Writes data to a file in text or Excel format, with optional highlighting.
write_file <- function(df, output, xlsx, append_bool = FALSE, colNames = TRUE, startRow=1, cols_to_highlight=NULL) {
    if(xlsx) {
        if (append_bool && file.exists(output)) {
            wb <- loadWorkbook(output)
            writeData(wb, sheet = "Results", df, startCol = 1, startRow = startRow, rowNames = FALSE, colNames = colNames)
            saveWorkbook(wb, file = output, overwrite = TRUE)        
        } else {
             write.xlsx(df, file = output, rowNames = FALSE, colNames = colNames, sheetName = "Results", encoding = "UTF-8", na = "", startRow = startRow)
        }

        if (!is.null(cols_to_highlight) && file.exists(output)) {
            wb <- loadWorkbook(output)
            for (col in cols_to_highlight) {
                addStyle(wb, "Results", style = createStyle(fgFill = "#FFFF00"), cols = col, rows = 1:(nrow(df)+1), gridExpand = TRUE)
            }
            saveWorkbook(wb, file = output, overwrite = TRUE)
        }
    } else {
        write.table(df, file = output, row.names=FALSE, col.names=colNames, dec=".",sep="\t",quote=FALSE , append=append_bool)
    }
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

data_genes <- NULL
if (bool_xlsx) {
    check_install_load_package("readxl")
    check_install_load_package("openxlsx")
}

# Formats non-zero read counts with sample names for easy interpretation.
getSignificantReads <- function(row, names, minUniqueReads) {
  reads <- row[row != 0 & row >= minUniqueReads] # Filter non-zero reads 
  indivs <- names[row != 0 & row >= minUniqueReads]  # Filter names of matching individuals
  if (length(indivs) == 0) {
    return("No reads")
  } else {
    return(paste(indivs, "reads = ", reads, collapse = "; "))
  }
}

# Counts the number of samples with read values that pass the MinNonStatisticalReads threshold.
# Filters out artifacts and ensures accurate sample representation in the analysis.
# Used to generate nbSampFilter for refined junction classification.
getNbSampFilter <- function(sampleReads) {
    if (sampleReads == "No reads") {
        return(0)
    }
    entries <- strsplit(sampleReads, "; ")[[1]]
    return(length(entries))
}


# Reads input data into a dataframe and filters it to extract Statistical and Non Statistical splice junctions.
# Defines input/output sample columns dynamically based on the data structure.
# Processes read counts and formats filtered results for output files.

data_genes <- read_file(input_file, bool_xlsx)
length_input <- 8 + length
EchName <- names(data_genes)[c(9:length_input)]
SampleOutput <- names(data_genes)[c((length_input+1):(length_input+length))] 
SampleInput <- EchName

if (as.integer(maxSignificantSamples) == as.integer(-1)) {
    maxSignificantSamples = as.integer(length)
}

if (as.integer(maxUniqueSamples) == as.integer(-1)) {
    maxUniqueSamples = as.integer(length)
}

significant_genes <- data_genes[!is.na(data_genes$Significative) & data_genes$Significative != "No" & !is.na(data_genes$filterInterpretation) & data_genes$event_type != "NoData" & data_genes$event_type != "Physio", ]
write_file(significant_genes, output_file_significant, bool_xlsx)
data_genes$filterInterpretation[data_genes$event_type == "NoData" & data_genes$AnnotJuncs == "Event too complex"] <- "Event too complex"
unique_genes <- data_genes[!is.na(data_genes$filterInterpretation) & (data_genes$event_type != "NoData" | data_genes$AnnotJuncs == "Event too complex" ) & data_genes$event_type != "Physio" & (data_genes$filterInterpretation == "Unique junction" | data_genes$filterInterpretation == "No model" | data_genes$filterInterpretation == "Percentage threshold execeeded" | data_genes$filterInterpretation == "Event too complex"), ]
unique_genes <- unique_genes[, colSums(is.na(unique_genes)) != nrow(unique_genes)]
unique_reads <- unique_genes[,EchName]
SampleReads <- apply(unique_reads, 1, getSignificantReads, names(unique_reads), minUniqueReads = minUniqueReads)
unique_genes <- cbind(unique_genes[, 1:( ncol(unique_genes) - 1)], SampleReads, unique_genes[,  ncol(unique_genes)])
nbSampFilter <- sapply(unique_genes$SampleReads, getNbSampFilter)
unique_genes <- cbind(unique_genes[, 1:( ncol(unique_genes) - 1)], nbSampFilter, unique_genes[,  ncol(unique_genes)])
colnames(unique_genes)[ncol(unique_genes)] <- "filterInterpretation"
write_file(unique_genes, output_file_unique, bool_xlsx)

# Ensure thresholdSignificanceLevel is correctly set 
# Default is 0 (no filtering) if the parameter is missing or invalid
if (is.null(thresholdSignificanceLevel) || !(thresholdSignificanceLevel %in% 0:3)) {
    thresholdSignificanceLevel <- 0  
}

# Define significance levels for filtering
significance_levels <- c("*", "**", "***")

# Filters Statistical splice junctions for each sample based on p-value thresholds.
# Writes results to separate files for Statistical junctions and filtered junctions.
# Highlights key columns in Excel outputs for visualization (if enabled).
message("######################")
message("#Filter Statistical genes...")
message("######################")

new_columns <- data.frame(nbSignificantSamples = NA, p_value = NA, SignificanceLevel = "")    
significant_genes <- cbind(significant_genes, new_columns)
significant_genes$nbSignificantSamples <- str_count(significant_genes$Significative,"p-value")

for(j in 1:length(EchName)){
    message(paste("   Filter Statistical genes for:",EchName[j]))

    dir.create(file.path(chem_dessin, EchName[j]), recursive = TRUE)

    output_EchName = paste(chem_dessin,EchName[j], "/", EchName[j],".statistical_junctions",extension,sep="")
    output_EchName_filter = paste(chem_dessin,EchName[j], "/", EchName[j],".statistical_junctions.filter",extension,sep="")
    P_EchName <- paste("P",EchName[j],sep="_") 
    cols_to_highlight <- which(names(significant_genes) %in% c(P_EchName, EchName[j],"nbSignificantSamples"))
    write_file(significant_genes[0,], output_EchName, bool_xlsx)
    write_file(significant_genes[0,], output_EchName_filter, bool_xlsx)

    k <- 2
    l <- 2
    for(i in 1:nrow(significant_genes))
        if(str_detect(significant_genes$Significative[i], EchName[j])){
            pattern <- paste0(EchName[j], ", p-value = ([0-9.eE+-]+)")
            p_value_sample <- as.numeric(str_match(significant_genes$Significative[i], pattern)[, 2])
            SignificanceLevel <- ''
            SignificanceLevel[p_value_sample<0.05 & p_value_sample>0.01] <- '*'
            SignificanceLevel[p_value_sample<0.01 & p_value_sample>0.001] <- '**'
            SignificanceLevel[p_value_sample<0.001] <- '***' 
            significant_genes$p_value[i] <- p_value_sample
            significant_genes$SignificanceLevel <- SignificanceLevel
            write_file(significant_genes[i,], output_EchName, bool_xlsx, append_bool=TRUE, colNames=FALSE, startRow=k)

            if (thresholdSignificanceLevel == 0 || which(significance_levels == SignificanceLevel) >= thresholdSignificanceLevel) {
                if (significant_genes$nbSignificantSamples[i] <= maxSignificantSamples) {
                    write_file(significant_genes[i,], output_EchName_filter, bool_xlsx, append_bool=TRUE, colNames=FALSE, startRow=l)
                        l <- l + 1
    }
}

            k <- k + 1
    }
    if (bool_xlsx) {
        output <- read_file(output_EchName,bool_xlsx)
        write_file(output, output_EchName, bool_xlsx, append_bool=TRUE, cols_to_highlight=cols_to_highlight)

        output_filter <- read_file(output_EchName_filter,bool_xlsx)
        write_file(output_filter, output_EchName_filter, bool_xlsx, append_bool=TRUE, cols_to_highlight=cols_to_highlight)
    }
}


# Filters unique splice junctions for each sample based on read count thresholds.
# Outputs full and filtered results to separate files, with optional column highlighting.
message("######################")
message("#Filter Non Statistical genes...")
message("######################")

unique_genes <- read_file(output_file_unique, bool_xlsx)
for(j in 1:length(EchName)){
    message(paste("   Filter Non Statistical genes for:",EchName[j]))
    P_EchName <- paste("P",EchName[j],sep="_") 
    output_EchName = paste(chem_dessin,EchName[j], "/", EchName[j],".non_statistical_junctions",extension,sep="")
    output_EchName_filter = paste(chem_dessin,EchName[j], "/", EchName[j],".non_statistical_junctions.filter",extension,sep="")
    
    unique_junctions <- unique_genes[(unique_genes[[P_EchName]] != 0 & !is.na(unique_genes[[P_EchName]]) | unique_genes[["filterInterpretation"]] == "Event too complex") & unique_genes[[EchName[j]]] >= minUniqueReads,]
    cols_to_highlight <- which(names(unique_junctions) %in% c(P_EchName, EchName[j], "nbSamp"))
    write_file(unique_junctions, output_EchName, bool_xlsx, append_bool=TRUE, cols_to_highlight=cols_to_highlight)

    unique_junctions_filter <- unique_genes[(unique_genes[[P_EchName]] != 0 & !is.na(unique_genes[[P_EchName]]) | unique_genes[["filterInterpretation"]] == "Event too complex") & unique_genes[[EchName[j]]] >= minUniqueReads & unique_genes$nbSampFilter <= maxUniqueSamples,]
    cols_to_highlight <- which(names(unique_junctions_filter) %in% c(P_EchName, EchName[j], "nbSamp"))
    write_file(unique_junctions_filter, output_EchName_filter, bool_xlsx, append_bool=TRUE, cols_to_highlight=cols_to_highlight)
}