mapAndEstimatePolymorphic <- function(R1unmappedReadsPaths, R2unmappedReadsPaths = NULL, strandedness = c("reverse", "forward", "none"),
                                      sampleNames, threads = 8, mappingsFolder = dirname(R1unmappedReadsPaths)[1],
                                      bowtie2polymorphicIndex, predictionsFolder = file.path(dirname(mappingsFolder), "polymorphicAbundances"),
                                      verbose = TRUE)
{

  if(missing(sampleNames))
    stop("Error: No sample names provided. Please provide a vector of sample names the same length as the read file vector(s).")
  strandedness <- match.arg(strandedness)
  BOWTIE2_INDEXES <- system("echo $BOWTIE2_INDEXES", intern = TRUE)
  strandMode <- switch(strandedness, none = 0, forward = 1, reverse = 2)
  
  # Mapping to masked genome.
  invisible(lapply(1:(length(R1readsPaths)), function(sampleNumber)
  {
    currentSample <- sampleNames[sampleNumber]
    bowtie2Command <- paste("bowtie2 -p", threads, "--local --no-mixed --no-unal -a --ma 0 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min C,0 -x", bowtie2polymorphicIndex)
    if(is.null(R2unmappedReadsPaths))
    {
      bowtie2Command <- paste(bowtie2Command, "-U", R1unmappedReadsPaths[sampleNumber])
    } else { # Paired files.
      bowtie2Command <- paste(bowtie2Command, "-1", R1unmappedReadsPaths[sampleNumber], "-2", R2unmappedReadsPaths[sampleNumber])
    }
    polymorphicSAM <- file.path(mappingsFolder, paste(currentSample, "ToPolymorphic.sam", sep = ''))
    bowtie2Command <- paste(bowtie2Command, "-S", polymorphicSAM)
    
    # Run Bowtie 2 mapping to polymorphic genes with the specified options.
    if(verbose == TRUE)
      message("Mapping reads not mapped to genome to polymorphic genes with Bowtie 2")
    system(bowtie2Command)
    
    RSEMcommand <- paste("rsem-calculate-expression -p", threads, "--sam")
    if(!is.null(R2unmappedReadsPaths))
    {
      RSEMcommand <- paste(RSEMcommand, "--paired-end")
    }
    if(grepl(BOWTIE2_INDEXES, bowtie2polymorphicIndex))
      bowtie2indexPath <- bowtie2polymorphicIndex
    else
      bowtie2indexPath <- file.path(BOWTIE2_INDEXES, bowtie2polymorphicIndex)
    RSEMcommand <- paste(RSEMcommand, "--strandedness", strandedness, polymorphicSAM, bowtie2indexPath, file.path(predictionsFolder, currentSample))
    
    # Run RSEM allele and gene expression estimation with the specified options.
    if(verbose == TRUE)
      message("Estimating allele expression with RSEM.")
    system(RSEMcommand)
    
    if(verbose == TRUE)
      message("Sorting and indexing reads mapped to alleles.")
    polymorphicAlignments <- list.files(predictionsFolder, paste(currentSample, "toPolymorphic.bam$", sep = ".*"), full.names = TRUE)
    system(paste("samtools sort -@ ", threads, " \'", polymorphicAlignments, '\'', " -o \'", gsub("bam", "sorted.bam", polymorphicAlignments), '\'', sep = ''))
    system(paste("mv '", gsub("bam", "sorted.bam", polymorphicAlignments), '\'', " \'", polymorphicAlignments, "\' ", sep = ''))
    system(paste("samtools index '", polymorphicAlignments, '\'', sep = ''))
    system(paste("rm", polymorphicSAM))
  }))
}