percentGroupsReads <- function(baseFolder, genesGFF3 = "/savona/nobackup/biostat/databases/GRCh38/gencode.v29.annotation.gff3")
{
  library(rtracklayer)
  
  abundanceFiles <- list.files(baseFolder, "genes.results", full.names = TRUE, recursive = TRUE)
  sampleIDs <- gsub("\\..*", '', basename(abundanceFiles))
  genesInfo <- readGFF(genesGFF3)
  genesInfo <- subset(genesInfo, type == "gene") # Reduce memory use.
  rowOrdering <- match(read.delim(abundanceFiles[1])[, "gene_id"], genesInfo[, "ID"])
  genesInfo <- genesInfo[rowOrdering, ]
  
  geneTypeSummary <- do.call(rbind, mapply(function(file, ID)
  {
    geneCounts <- read.delim(file)
    countsPerGroup <- sapply(split(geneCounts[, "expected_count"], genesInfo[, "gene_type"]), sum)
    percentsPerGroup <- round(countsPerGroup / sum(countsPerGroup) * 100)
    data.frame(sample = ID, `mRNA Percent` = percentsPerGroup["protein_coding"],
               `lincRNA Percent` = percentsPerGroup["lincRNA"],
               `rRNA Percent` = percentsPerGroup["rRNA"], check.names = FALSE)
  }, abundanceFiles, sampleIDs, SIMPLIFY = FALSE))
  
  write.table(geneTypeSummary, file.path(baseFolder, "biotypeSummary.txt"), row.names = FALSE, quote = FALSE, sep = '\t')
}