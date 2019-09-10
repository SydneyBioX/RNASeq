predictQuantifyPolymorphicAlleles <- function(mappingsFolder, predictionsFolder, samplesTotalReads, minCoverage = 2, maxFlankingBasesUncovered = 100, verbose = TRUE)
{
  require(GenomicAlignments)
  require(BiocParallel)
  
  allelesQuantFiles <- list.files(predictionsFolder, "isoforms.results", full.names = TRUE)
  sampleIDs <- basename(gsub(".isoforms.results", '', allelesQuantFiles))
  alleleTables <- lapply(allelesQuantFiles, function(file) read.delim(file, stringsAsFactors = FALSE)[, c("transcript_id", "gene_id", "expected_count", "IsoPct")])
  
  datasetAllelesTable <- do.call(cbind, mcmapply(function(alleleTable, sampleName, sampleTotalReads)
  {
    if(verbose == TRUE)
      message("Estimating alleles and their expression for ", sampleName)
    polymorphicAlignments <- list.files(mappingsFolder, paste(sampleName, "toPolymorphic.bam$", sep = ".*"), full.names = TRUE)
    polymorphicGenomicSequenceMapping  <- readGAlignmentPairs(polymorphicAlignments, param = ScanBamParam(mapqFilter = 1), strandMode = 0, use.names = TRUE)
    readsToPolymorphic <- length(unique(names(polymorphicGenomicSequenceMapping)))
    allelesCoverage <- coverage(polymorphicGenomicSequenceMapping)

    do.call(rbind, by(alleleTable, alleleTable[, "gene_id"], function(allelesByGeneTable)
    {
      topAlleles <- allelesByGeneTable[order(allelesByGeneTable[, "expected_count"], decreasing = TRUE)[1:min(nrow(allelesByGeneTable), 2)], ]
      alleleStatus <- sapply(topAlleles[, "transcript_id"], function(alleleID)
      {
        aCoverage <- allelesCoverage[[alleleID]]
        whichUncovered <- aCoverage == 0
        if(all(whichUncovered[1:(maxFlankingBasesUncovered/2)]) || all(whichUncovered[(length(whichUncovered) - (maxFlankingBasesUncovered/2) + 1):length(whichUncovered)]))
          return("Not in Database")
        if(any(whichUncovered[(maxFlankingBasesUncovered/2 + 1):(length(whichUncovered) - (maxFlankingBasesUncovered/2))]))
          return("Not in Database")
        return("Detected")
      })
      
      whichFailed <- which(alleleStatus == "Not in Database")
      if(length(whichFailed) > 0)
      {
        for(failed in whichFailed)
        {
          if(all(allelesCoverage[[failed]] == 0) || mean(allelesCoverage[[failed]][allelesCoverage[[failed]] != 0]) < minCoverage)
            alleleStatus[whichFailed] <- "Not Detected"
        }
      }

      alleleStatus[alleleStatus == "Detected"] <- topAlleles[alleleStatus == "Detected", "transcript_id"]
      if(length(alleleStatus) == 1) alleleStatus <- c(alleleStatus, NA)
      alleleStatus <- gsub(".*\\*", '', alleleStatus)
      counts <- topAlleles[, "expected_count"]
      counts[alleleStatus %in% c("Not in Database", "Not Detected", NA)] <- NA
      counts <- round(counts)
      lengths <- sapply(allelesCoverage[topAlleles[, "transcript_id"]], length)
      lengths[alleleStatus %in% c("Not in Database", "Not Detected", NA)] <- NA
      FPKMs <- round((counts / (1000 / lengths) / (sampleTotalReads / 1000000)), 2)
      
      alleles <- c(alleleStatus, lengths, counts, FPKMs)
      alleles <- t(data.frame(alleles))
      colnames(alleles) <- c(paste(sampleName, "Allele", 1:2),
                             paste(sampleName, "Allele", 1:2, "Length (nt)"),
                             paste(sampleName, "Allele", 1:2, "Read Count"),
                             paste(sampleName, "Allele", 1:2, "FPKM"))
      alleles
    }))
  }, alleleTables, sampleIDs, samplesTotalReads, mc.cores = threads, SIMPLIFY = FALSE))
  rownames(datasetAllelesTable) <- unique(alleleTables[[1]][, "gene_id"])
  write.table(datasetAllelesTable, file.path(predictionsFolder, file = "polymorphicAllelesValues.txt"), sep = '\t', quote = FALSE)
  
  system(paste("rm -r", file.path(mappingsFolder, paste('*', maskedGenomeName, '*', sep = ''))))
  
  if(verbose == TRUE)
    message("Allele predictions saved to disk.")
}