HLAfile <- "/savona/nobackup/biostat/databases/IMGT/hla_nuc337.fasta"
KIRfile <- "/savona/nobackup/biostat/databases/IMGT/KIR_nuc280.fasta"
outputFile <- "/savona/nobackup/biostat/databases/IMGT/poly_nuc_uniqueH337K280.fasta"
preparePolymorphicSequences(HLAfile, KIRfile, outputFile)

allAllelesFASTA <- "/savona/nobackup/biostat/databases/IMGT/poly_nuc_uniqueH337K280.fasta"
outputGTFfile <- "/savona/nobackup/biostat/databases/IMGT/polymorphicGenesH337K280.gtf"
createPolymorphicGTF(allAllelesFASTA, outputGTFfile)

genesGFF3 <- "/savona/nobackup/biostat/databases/GRCh38/gencode.v31.annotation.gff3"
hg38sequenceFile <- "/savona/nobackup/biostat/sequence/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
STARmaskedFolder <- "/savona/nobackup/biostat/indexes/STAR/hg38maskedWithGENCODE31"
bowtie2IndexPrefix <- "/savona/nobackup/biostat/indexes/bowtie2/polymorphic/polyUniqH337K280"
prepareAndIndexSequences(genesGFF3, hg38sequenceFile, STARmaskedFolder, bowtie2IndexPrefix, allAllelesFASTA, outputGTFfile,
                         readLength = 150)

R1s <- list.files("/savona/nobackup/biostat/datasets/HelenRizos/AGRF_CAGRF10479_C79UVANXX",
                  ".*R1.*trimmed.*gz", full.names = TRUE) 
R2s <- list.files("/savona/nobackup/biostat/datasets/HelenRizos/AGRF_CAGRF10479_C79UVANXX",
                  ".*R2.*trimmed.*gz", full.names = TRUE)
sampleNames <- sapply(lapply(strsplit(basename(R1s), '_'), '[', 1:2), paste, collapse = '_')
predictQuantifyPolymorphicAlleles(R1s, R2s, sampleNames = sampleNames, maskedGenomePath = "/savona/nobackup/biostat/indexes/STAR/hg38maskedWithGENCODE26/",
                          maskedGenomeName = "hg38masked", threads = 8, bowtie2polymorphicIndex = "polymorphic/polyUniq", predictionsFolder = "/savona/nobackup/dario/TestRun/"
                          )

RNAfolders <- c("/savona/nobackup/biostat/datasets/HelenRizos/AGRF_CAGRF10479_C79UVANXX",
                "/savona/nobackup/biostat/datasets/HelenRizos/AGRF_CAGRF11240_C8768ANXX/",
                "/savona/nobackup/biostat/datasets/HelenRizos/AGRF_CAGRF11395_C8768ANXX/",
                "/savona/nobackup/biostat/datasets/HelenRizos/AGRF_CAGRF11830_C9CWMANXX/",
                "/savona/nobackup/biostat/datasets/HelenRizos/AGRF_CAGRF13546_CAFANANXX/",
                "/savona/nobackup/biostat/datasets/HelenRizos/AGRF_CAGRF14197_CAT5WANXX/",
                "/savona/nobackup/biostat/datasets/HelenRizos/AGRF_CAGRF14665_CAUJPANXX/"
                )

invisible(lapply(RNAfolders, function(folder)
{
  R1s <- list.files(folder, ".*R1.*trimmed.*gz", full.names = TRUE) 
  R2s <- list.files(folder, ".*R2.*trimmed.*gz", full.names = TRUE)
  sampleNames <- sapply(lapply(strsplit(basename(R1s), '_'), '[', 1:2), paste, collapse = '_')
  predictQuantifyPolymorphicAlleles(R1s, R2s, sampleNames = sampleNames, maskedGenomePath = "/savona/nobackup/biostat/indexes/STAR/hg38maskedWithGENCODE26/",
                                    maskedGenomeName = "hg38masked", threads = 8, bowtie2polymorphicIndex = "polymorphic/polyUniq", predictionsFolder = folder
                                    )
}))