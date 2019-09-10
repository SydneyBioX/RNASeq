processRNAreads <- function(baseFolder, sampleNameDelimiter = '_', sampleNameParts = 1, sampleNameReadOrder = 5,
                            adapterFirst = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", adapterSecond = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA",
                            trimOptions = "-q 20,20 -m 50", organismGenome = "hg38", STARgenomePath = "/savona/nobackup/biostat/indexes/STAR/hg38maskedWithGENCODE31/",
                            RSEMreference = "/savona/nobackup/biostat/databases/GRCh38/RSEMgencode31masked",
                            bowtie2polymorphicIndex = "/savona/nobackup/biostat/indexes/bowtie2/polymorphic/polyUniqH337K280",
                            stranded = c("reverse", "forward", "none"), genesGFF3 = "/savona/nobackup/biostat/databases/GRCh38/gencode.v31.annotation.gff3",
                            cores = 8, doCombine = TRUE, doTrim = TRUE, multipleDatasetsInStudy = FALSE,
                            datasetID = strsplit(basename(baseFolder), '_')[[1]][2],
                            datasetsCountsDirectory = "/savona/nobackup/biostat/datasets/HeadNeckCancer/RNAseq")
{
  library(parallel)
  library(BiocParallel)
  library(GenomicFeatures)
  library(GenomicAlignments)
  library(rtracklayer)
  library(edgeR)
  library(limma)
  library(RUVnormalize)

  stranded <- match.arg(stranded)

  # Concatenate samples of patient sequenced in different lanes.
  if(doCombine == TRUE)
  {
    originalFiles <- list.files(baseFolder, ".*L00.*gz$", full.names = TRUE, recursive = TRUE)
    filenameComponents <- strsplit(basename(originalFiles), sampleNameDelimiter)
    sampleAndReadSideID <- sapply(filenameComponents, function(components) paste(components[c(sampleNameParts, sampleNameReadOrder)], sep = '', collapse = '_'))
    fileGroupings <- split(1:length(originalFiles), sampleAndReadSideID)
    invisible(mapply(function(grouping, combinedName)
	  {
	    if(!grepl("fastq.gz", combinedName))
	      combinedName <- paste(combinedName, "fastq.gz", sep = '.')
		  system(paste("cat ", paste('\'', originalFiles[grouping], '\'', collapse = ' ', sep = ''), " > ", '\'', file.path("/savona/nobackup/biostat/datasets/HeadNeckCancer/RNAseq/merged", combinedName), '\'', sep = ''))
	  }, fileGroupings, names(fileGroupings)))
  }
  # Output files: 1 pair (R1, R2) per patient for paired-end data, or one file per single-end data.
  # Calculate reads per sample before any preprocessing is done.
  mergedReads <- list.files(baseFolder, ".fastq.gz$", full.names = TRUE, recursive = TRUE)
  if(doCombine == TRUE) 
    ignorePattern <- "L00|trimmed"
  else
    ignorePattern <- "trimmed"
  ignoreFiles <- grepl(ignorePattern, mergedReads)
  mergedReads <- mergedReads[!ignoreFiles] # Not separate lanes and not already trimmed.
  if(length(mergedReads) > 0)
  {
    if(is.numeric(sampleNameReadOrder)) # Has components such as R1, R2.
    {
      filenameComponents <- strsplit(basename(mergedReads), sampleNameDelimiter)
      sampleID <- sapply(filenameComponents, function(components) paste(components[1:length(sampleNameParts)], sep = '', collapse = '_'))
      mergedReads <- split(mergedReads, sampleID)
    }
    readsOriginal <- unlist(mclapply(mergedReads, function(readFiles)
    {
      as.integer(system(paste("echo $(zcat '", readFiles[1], "' | wc -l) / 4 | bc", sep = ''), intern = TRUE))
    }, mc.cores = 5))
  }
  
  # Do adapter trimming to remove sequences due to kit or lots of low-quality bases.
  # cutadapt's containing folder needs to be in PATH variable.
  if(doTrim == TRUE)
  {
    if(doCombine == FALSE)
    {
      filenameComponents <- strsplit(basename(mergedReads), sampleNameDelimiter)
      sampleAndReadSideID <- sapply(filenameComponents, function(components)
                             paste(components[c(sampleNameParts, sampleNameReadOrder)], sep = '', collapse = '_'))
      fileGroupings <- split(1:length(mergedReads), sampleAndReadSideID)
      mergedReadsRenamed <- file.path(baseFolder, paste(sampleAndReadSideID, ".fastq.gz", sep = ''))
      sampleID <- sapply(filenameComponents, function(components) paste(components[sampleNameParts], sep = '', collapse = '_'))
      mergedReadsRenamed <- split(mergedReadsRenamed, sampleID)
    }
    else{
      mergedReadsRenamed <- mergedReads
    }
    
    invisible(mcmapply(function(oldNames, newNames)
    {
      cutCommand <- paste("cutadapt -a", adapterFirst)
      if(is.numeric(sampleNameReadOrder)) # Paired-end
        cutCommand <- paste(cutCommand, "-A", adapterSecond)
      cutCommand <- paste(cutCommand, trimOptions)
      cutCommand <- paste(cutCommand, " -o ", '\'', gsub("fastq.gz", "trimmed.fastq.gz", newNames[1]), '\'', sep = '')
      if(is.numeric(sampleNameReadOrder)) # Paired-end
        cutCommand <- paste(cutCommand, " -p ", '\'', gsub("fastq.gz", "trimmed.fastq.gz", newNames[2]), '\'', sep = '')
      cutCommand <- paste(cutCommand, ' \'', oldNames[1], '\'', sep = '')
      if(is.numeric(sampleNameReadOrder)) # Paired-end
        cutCommand <- paste(cutCommand, ' \'', oldNames[2], '\'', sep = '')
      system(cutCommand)
    }, mergedReads, mergedReadsRenamed, mc.cores = cores))
  }
  
  if(doTrim == TRUE)
  {
    mappableReads <- list.files(baseFolder, "trimmed.*fastq.gz", full.names = TRUE, recursive = TRUE)
    filenameComponents <- strsplit(basename(mappableReads), sampleNameDelimiter)
    sampleID <- sapply(filenameComponents, function(components) paste(components[1:length(sampleNameParts)], sep = '', collapse = '_'))
    mappableReads <- split(mappableReads, sampleID)
  } else { # Adapter trimming not done.
    mappableReads <- mergedReads
  }
  
  # Count the number of reads remaining after trimming.
  readsToMap <- unlist(mclapply(mappableReads, function(readFiles)
  {
    as.integer(system(paste("echo $(zcat '", readFiles[1], "' | wc -l) / 4 | bc", sep = ''), intern = TRUE))
  }, mc.cores = 5))
  
  # Map the FASTQ files to the genome with STAR.
  
  # First, load genome index into RAM.
  system(paste("STARnew --genomeDir", STARgenomePath, "--genomeLoad LoadAndExit"))
  
  # Second, do the mapping to create BAM files.
  invisible(mapply(function(sampleReads, sampleName)
  {
    mapCommand <- paste("STARnew --runThreadN 16 --readFilesIn ", '\'', sampleReads[1], '\' ', sep = '')
    if(is.numeric(sampleNameReadOrder))
    {
      for(index in 2:length(sampleReads))
        mapCommand <- paste(mapCommand, '\'', sampleReads[index], '\'', sep = '')
    }
    mapCommand <- paste(mapCommand, "--readFilesCommand gunzip -c --genomeDir", STARgenomePath, "--genomeLoad LoadAndKeep --outFilterMultimapNmax 999 --outFilterMismatchNoverReadLmax 0.1")
    if(is.numeric(sampleNameReadOrder))
      mapCommand <- paste(mapCommand, "--chimSegmentMin 20 --alignMatesGapMax 1000000 --alignIntronMax 1000000 --alignEndsProtrude 10 ConcordantPair --chimSegmentReadGapMax 3 --peOverlapNbasesMin 10 --peOverlapMMp 0.1")
    mapCommand <- paste(mapCommand, "--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --limitBAMsortRAM 16000000000 --chimOutJunctionFormat 1 --chimMultimapNmax 1")
    mapCommand <- paste(mapCommand, " --outReadsUnmapped Fastx --outFileNamePrefix ", '\'', paste(file.path(baseFolder, "mapped/"), sampleName, organismGenome, sep = ''), '\'', sep = '')
    system(mapCommand)
  }, mappableReads, names(mappableReads)))
  
  # Create indexes for BAM files (i.e. Filenames with .bai suffix).
  alignedGenomeFiles <- list.files(baseFolder, paste(organismGenome, "Aligned.sortedByCoord.out.bam$", sep = ''), full.names = TRUE, recursive = TRUE)
  sampleIDs <- gsub(paste(organismGenome, "Aligned.*", sep = ''), '', basename(alignedGenomeFiles))
  names(alignedGenomeFiles) <- sampleIDs
  
  invisible(mclapply(alignedGenomeFiles, function(alignedSample)
  {
    system(paste("samtools index '", alignedSample, '\'', sep = ''))
  }, mc.cores = 4)) # Don't change to any more than 4, since lots of disk access is done.
  
  alignmentReports <- list.files(file.path(baseFolder, "mapped"), "Log.final.out", full.names = TRUE)
  reportsSample <- gsub(paste(organismGenome, "Log.final.out", sep = ''), '', basename(alignmentReports))
  readsAlignedGenome <- sapply(alignmentReports, function(reportFile)
  {
    alignSummary <- read.delim(reportFile, sep = '\t', stringsAsFactors = FALSE, header = FALSE)
    as.numeric(alignSummary[5, 2]) - (as.numeric(alignSummary[28, 2]) + as.numeric(alignSummary[30, 2]) + as.numeric(alignSummary[32, 2])) # Total - unmapped.
  })
  names(readsAlignedGenome) <- reportsSample
  readsAlignedGenome <- readsAlignedGenome[match(names(readsOriginal), names(readsAlignedGenome))]
  
  # Estimate gene counts using RSEM.
  
  alignedTranscriptomeFiles <- list.files(baseFolder, paste(organismGenome, "Aligned.toTranscriptome.out.bam", sep = ''), full.names = TRUE, recursive = TRUE)
  sampleIDs <- gsub(paste(organismGenome, "Aligned.toTranscriptome.out.bam", sep = ''), '', basename(alignedTranscriptomeFiles))
  invisible(mapply(function(alignFile, sample)
  {
    system(paste("rsem-calculate-expression --strandedness", stranded, "--alignments", if(length(mappableReads[[1]]) == 2) "--paired-end", "-p", cores, alignFile, RSEMreference, file.path(baseFolder, "genomicAbundances", sample)))
  }, alignedTranscriptomeFiles, sampleIDs))
  
  # Import the estimated counts from RSEM and combine into a table.
  sampleEstimatesFiles <- list.files(file.path(baseFolder, "genomicAbundances"), "genes.results", full.names = TRUE)
  sampleIDs <- gsub(".genes.results", '', basename(sampleEstimatesFiles))
  allGeneData <- lapply(sampleEstimatesFiles, function(estimatesFile) 
  {
    sampleEstimates <- read.delim(estimatesFile, stringsAsFactors = FALSE)
    sampleEstimates[, c("gene_id", "effective_length", "expected_count")]
  })
  
  countsTable <- do.call(cbind, lapply(allGeneData, '[[', "expected_count"))
  countsTable <- round(countsTable) # Make integers.
  colnames(countsTable) <- sampleIDs
  
  # Also explort FPKM table for tasks such as cell deconvolution or clustering.
  
  geneInfo <- readGFF(genesGFF3)
  RSEMgeneIDs <- read.delim(sampleEstimatesFiles[1], stringsAsFactors = FALSE)[, "gene_id"]
  rowSymbols <- geneInfo[match(RSEMgeneIDs, geneInfo[, "gene_id"]), "gene_name"]
  rownames(countsTable) <- rowSymbols
  duplicates <- names(table(rowSymbols)[table(rowSymbols) > 1])
  
  mergedDuplicatesTable <- do.call(rbind, lapply(duplicates, function(duplicateGene)
  {
    colSums(countsTable[which(rowSymbols == duplicateGene), , drop = FALSE])
  }))
  rownames(mergedDuplicatesTable) <- duplicates
  dropRows <- which(rowSymbols %in% duplicates)
  countsTable <- countsTable[-dropRows, , drop = FALSE]
  countsTable <- rbind(countsTable, mergedDuplicatesTable)
  
  # Estimate HLA and KIR genes and their alleles using IMGT HLA and KIR databases and incorporate their values into the gene expression table.
  
  unalignedGenomeFilesR1 <- list.files(file.path(baseFolder, "mapped"), paste(organismGenome, "Unmapped.out.mate1", sep = ''), full.names = TRUE)
  unalignedGenomeFilesR2 <- list.files(file.path(baseFolder, "mapped"), paste(organismGenome, "Unmapped.out.mate2", sep = ''), full.names = TRUE)
  sampleNames <- gsub(paste(organismGenome, ".*", sep = ''), '', basename(unalignedGenomeFilesR1))
  mapAndEstimatePolymorphic(unalignedGenomeFilesR1, unalignedGenomeFilesR2, sampleNames = sampleNames,
                                    bowtie2polymorphicIndex = bowtie2polymorphicIndex)
  
  # Incorporate the estimated polymorphic gene counts into the all-gene count table and correctly calculate FPKMs.
  
  polyGenesQuantFiles <- list.files(file.path(baseFolder, "polymorphicAbundances"), "genes.results", full.names = TRUE)
  polyGenesQuantFiles <- polyGenesQuantFiles[match(names(readsOriginal), gsub(".genes.results", '', basename(polyGenesQuantFiles)))]

  polyGeneData <- lapply(polyGenesQuantFiles, function(estimatesFile) 
  {
    sampleEstimates <- read.delim(estimatesFile, stringsAsFactors = FALSE)
    sampleEstimates[, c("gene_id", "effective_length", "expected_count")]
  })
  
  # IMGT database and GENCODE use different names for some genes. These are the GENCODE versions.
  
  polyGeneIDs <- c("HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPA2", "HLA-DPB1", "HLA-DPB2", "HLA-DQA1", "HLA-DQB1",
                   "HLA-DRA", "HLA-DRB1", "HLA-DRB2", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRB6", "HLA-DRB7", "HLA-DRB8", "HLA-DRB9",
                   "HFE", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-H", "HLA-J", "HLA-K", "HLA-L", "HLA-N", "HLA-S", "HLA-T",
                   "HLA-U", "HLA-V", "HLA-W", "HLA-Y", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DP1",
                   "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DP1", "KIR3DS1",
                   "MICA", "MICB", "TAP1", "TAP2")
  IMGTpolyIDs <- polyGeneData[[1]][, "gene_id"]
  
  polyGeneCounts <- round(do.call(cbind, lapply(polyGeneData, "[[", "expected_count")))
  polyInAllRows <- match(polyGeneIDs, rownames(countsTable))
  polyGeneCounts <- polyGeneCounts[!is.na(polyInAllRows), , drop = FALSE]
  polyNotInGenesDB <- !is.na(polyInAllRows)
  polyGeneIDs <- polyGeneIDs[polyNotInGenesDB] 
  IMGTpolyIDs <- IMGTpolyIDs[polyNotInGenesDB]
  
  countsTable[match(polyGeneIDs, rownames(countsTable)), ] <- polyGeneCounts
  
  countsFile <- "unnormalisedGeneCounts.txt"
  countsPath <- file.path(baseFolder, countsFile)
  write.table(countsTable, countsPath, sep = '\t', quote = FALSE)
  
  # Calculate FPKM based on effective gene length, read count per gene and total read count.
  
  allRows <- match(rownames(countsTable), geneInfo[match(allGeneData[[1]][, "gene_id"], geneInfo[, "gene_id"]), "gene_name"])
  polyRows <- match(IMGTpolyIDs, polyGeneData[[1]][, "gene_id"])
  polyRowsInAll <- match(polyGeneIDs, rownames(countsTable))
  mapply(function(sampleCounts, allData, polyData)
  {
    allLengths <- allData[allRows, "effective_length"]
    polyLengths <- polyData[polyRows, "effective_length"]
    FPKMs <- sampleCounts / (allLengths / 1000) / (sum(sampleCounts) / 1000000)
    FPKMs[is.nan(FPKMs)] <- 0 # If count and length are 0, this causes NaN.
    FPKMs[polyRowsInAll] <- sampleCounts[polyRowsInAll] / (polyLengths / 1000) / (sum(sampleCounts) / 1000000)
  }, as.data.frame(countsTable), allGeneData, polyGeneData)
  
  FPKMfile <- "geneFPKMs.txt"
  FPKMpath <- file.path(baseFolder, FPKMfile)
  write.table(FPKMsTable, FPKMpath, sep = '\t', quote = FALSE)
  
  if(multipleDatasetsInStudy == TRUE)
  {
    datasetCountsFile <- paste(datasetID, countsFile, sep = '_')
    allCountsPath <- file.path(datasetsCountsDirectory, datasetCountsFile)
    write.table(countsTable, allCountsPath, sep = '\t', quote = FALSE)
  }
  
  # Estimate which alleles are expressed.
  
  predictQuantifyPolymorphicAlleles(file.path(baseFolder, "mapped"), file.path(baseFolder, "polymorphicAbundances"), colSums(countsTable))
  
  # TMM normalise the counts, then normalise with voom, which is essentially just log2-CPM.
  DGEListObject <- DGEList(countsTable)
  DGEListObject <- calcNormFactors(DGEListObject) # TMM.
  
  # Filter out genes with large CPM values in only one or less samples and voom transform.
  keepGenes <- rowSums(cpm(DGEListObject) >= 3) >= 2
  DGEListObject <- DGEListObject[keepGenes, , keep.lib.sizes = FALSE]
  normalisedValues <- voom(DGEListObject)
  colnames(normalisedValues[["weights"]]) <- colnames(normalisedValues[['E']])
  rownames(normalisedValues[["weights"]]) <- rownames(normalisedValues[['E']])
  
  valuesPath <- file.path(baseFolder, "normalisedValues.txt")
  weightsPath <- file.path(baseFolder, "valuesWeights.txt")
  write.table(normalisedValues[['E']], valuesPath, sep = '\t', quote = FALSE)
  write.table(normalisedValues[["weights"]], weightsPath, sep = '\t', quote = FALSE)
  
  if(multipleDatasetsInStudy == TRUE) # Do the same for the dataset tables combined.
  {
    allCountsTables <- lapply(list.files(datasetsCountsDirectory, ".*_unnormalisedGeneCounts.txt*", full.names = TRUE), read.delim, check.names = FALSE, stringsAsFactors = FALSE)
    allCountsTable <- do.call(cbind, allCountsTables)
    write.table(allCountsTable, file.path(datasetsCountsDirectory, "unnormalisedGeneCounts.txt"), sep = '\t', quote = FALSE)
    allDGEListObject <- DGEList(allCountsTable)
    allDGEListObject <- calcNormFactors(allDGEListObject) # TMM.
    
    # Filter out genes with large CPM values in only one or less samples and voom transform.
    keepGenes <- rowSums(cpm(allDGEListObject) >= 3) >= 2
    allDGEListObject <- allDGEListObject[keepGenes, , keep.lib.sizes = FALSE]
    allDatasetsNormalisedValues <- voom(allDGEListObject)
    colnames(allDatasetsNormalisedValues[["weights"]]) <- colnames(allDatasetsNormalisedValues[['E']])
    rownames(allDatasetsNormalisedValues[["weights"]]) <- rownames(allDatasetsNormalisedValues[['E']])
    
    valuesPath <- file.path(datasetsCountsDirectory, "normalisedValues.txt")
    weightsPath <- file.path(datasetsCountsDirectory, "valuesWeights.txt")
    allValues <- allDatasetsNormalisedValues[['E']]
    write.table(allValues, valuesPath, sep = '\t', quote = FALSE)
    write.table(allDatasetsNormalisedValues[["weights"]], weightsPath, sep = '\t', quote = FALSE)
    
    samplesDatasetsFilePath <- file.path(datasetsCountsDirectory, "samplesDatasets.txt")
    samplesDatasets <- read.delim(samplesDatasetsFilePath, stringsAsFactors = FALSE)
    samplesDatasets <- rbind(samplesDatasets, data.frame(Dataset = datasetID, Sample = sampleIDs))
    write.table(samplesDatasets, samplesDatasetsFilePath, row.names = FALSE, sep = '\t', quote = FALSE)
    
    # Provide RUV normalised table. 100 least variable genes as controls.
    leastVariable <- rownames(allValues)[order(apply(allValues, 1, sd))[1:100]]
    batchCorrectedAll <- t(naiveRandRUV(t(allValues), leastVariable, nu.coeff = 1, k = 2))
    valuesPath <- file.path(datasetsCountsDirectory, "normalisedValuesRUVtransformed.txt")
    write.table(batchCorrectedAll, valuesPath, sep = '\t', quote = FALSE)
  }
  
  # Create a read total summary each sample.
  if(length(mergedReads) > 0)
  {
    readsSummary <- data.frame(Sample = names(readsToMap),
                               `Sequenced Reads` = readsOriginal,
                               `Trimmed and Quality Filtered Reads` = readsToMap,
                               `Mapped to Genome` = readsAlignedGenome, check.names = FALSE)
  } else {
    readsSummary <- data.frame(Sample = names(readsToMap),
                               `Trimmed and Quality Filtered Reads` = readsToMap,
                               `Mapped to Genome` = readsAlignedGenome, check.names = FALSE)    
  }
  write.table(readsSummary, file.path(baseFolder, "samplesTotalReads.txt"), row.names = FALSE, sep = '\t', quote = FALSE)
}

# Genome index
# STARnew --runThreadN 8 --runMode genomeGenerate --genomeDir /savona/nobackup/biostat/indexes/STAR/hg38withGENCODE29 --genomeFastaFiles /savona/nobackup/biostat/sequence/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --sjdbGTFfile /savona/nobackup/biostat/databases/GRCh38/gencode.v29.annotation.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 149

# RSEM index
# rsem-prepare-reference --gtf /savona/nobackup/biostat/databases/GRCh38/gencode.v29.annotation.gtf /savona/nobackup/biostat/sequence/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna /savona/nobackup/biostat/databases/GRCh38/RSEMgencode29
# rsem-prepare-reference --gtf /savona/nobackup/biostat/databases/GRCh38/gencode.v31.annotation.gtf /savona/nobackup/biostat/sequence/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.masked.fna /savona/nobackup/biostat/databases/GRCh38/RSEMgencode31masked

# processRNAreads("/savona/nobackup/biostat/datasets/HeadNeckCancer/RNAseq")