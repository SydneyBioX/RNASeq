prepareAndIndexSequences <- function(genesGFF3, genomeFASTA, STARindexFolder, bowtie2IndexPrefix, polymorphicFASTAfile, polymorphicGTFfile, readLength = 100)
{
  require(ShortRead)
  require(rtracklayer)
  polymorphicGenes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DOA", "HLA-DRA", 
                        "HLA-DMB", "HLA-DOB", "HLA-DRB1", "HLA-DRB2", "HLA-DRB3", "HLA-DRB4", 
                        "HLA-DRB5", "HLA-DRB6", "HLA-DRB7", "HLA-DRB8", "HLA-DRB9", "HLA-DPA1", 
                        "HLA-DPA2", "HLA-DPB1", "HLA-DPB2", "HLA-DQA1", "HLA-DQB1", "HLA-E", 
                        "HLA-F", "HLA-G", "HFE", "HLA-H", "HLA-J", "HLA-K", "HLA-L", 
                        "HLA-V", "TAP1", "TAP2", "MICA", "MICB", "KIR3DL1", "KIR3DL3", 
                        "KIR2DL4", "KIR2DS5", "KIR3DL2", "KIR2DS1", "KIR2DP1", "KIR2DS3", 
                        "KIR2DL5A", "KIR3DP1", "KIR2DS4", "KIR2DS2", "KIR2DL5B", "KIR2DL2", 
                        "KIR2DL3", "KIR3DS1", "KIR2DL1") # Constant.
  
  genesLocations <- import.gff3(genesGFF3)
  genesLocations <- subset(genesLocations, type == "gene")
  polyGenesLocations <- subset(genesLocations, gene_name %in% polymorphicGenes)
  genesOverlapPoly <- subjectHits(findOverlaps(polyGenesLocations, genesLocations, ignore.strand = TRUE))
  mcols(polyGenesLocations) <- NULL
  overlappingPolyGenesLocations <- GRanges(seqnames(genesLocations)[genesOverlapPoly],
                                           IRanges(start(genesLocations)[genesOverlapPoly], end(genesLocations)[genesOverlapPoly]))
  polyGenesLocations <- c(polyGenesLocations, overlappingPolyGenesLocations)
  polyGenesLocations <- as(polyGenesLocations, "data.frame")
  
  chromosomesToMask <- as.character(unique(polyGenesLocations[, "seqnames"]))
  genomeSequence <- readFasta(genomeFASTA)
  genomeChrs <- sapply(strsplit(as.character(id(genomeSequence)), ' '), '[', 1)
  lapply(chromosomesToMask, function(chromosome)
  {
    whichChr <- which(genomeChrs == chromosome)
    chromosomeSequence <- sread(genomeSequence)[[whichChr]]
    maskCoordinates <- subset(polyGenesLocations, seqnames == chromosome)[, c("start", "end")]
    mapply(function(start, end)
    {
      chromosomeSequence[start:end] <<- 'N'
    }, maskCoordinates[, "start"], maskCoordinates[, "end"])
    genomeSequence@sread[[whichChr]] <<- chromosomeSequence
  })
  
  maskedGenomeFile <- gsub("^(.*)\\.", "\\1.masked.", genomeFASTA)
  writeFasta(genomeSequence, file = maskedGenomeFile)
  indexFa(maskedGenomeFile)
    
  # Generate a genome index for STAR
  if(!dir.exists(STARindexFolder))
    dir.create(STARindexFolder)
  STARindexCommand <- paste("STARnew --runMode genomeGenerate --genomeDir", STARindexFolder, "--genomeFastaFiles", maskedGenomeFile,
                            "--sjdbGTFfile", genesGFF3, "--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang", readLength - 1)
  system(STARindexCommand)
  
  # Generate a polymorphic genes index for Bowtie 2
  Bowtie2indexCommand <- paste("bowtie2-build --threads ", 8, " \'", polymorphicFASTAfile, "\' ", "\'", bowtie2IndexPrefix, "\'", sep = '')
  system(Bowtie2indexCommand)
  
  # Generate a reference for RSEM
  RSEMrefCommand <- paste("rsem-prepare-reference --gtf \'", polymorphicGTFfile, "\' \'", polymorphicFASTAfile, "\' '", bowtie2IndexPrefix, "\'", sep = '')
  system(RSEMrefCommand)
}