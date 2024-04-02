#### Overview ####
# Test to check for overlaps

#### Preparation ####
# loading libraries
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(GenomicAlignments)
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(ggbio) # https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
library(dplyr)
library(tidyr)
library(Rsamtools)
library(ggplot2)
library(BSgenome)
library(karyoploteR)

# BiocManager::install("karyoploteR")
library(msa)

#require(BiocManager)

#### loading data ####
setwd("/powerplant/workspace/cfngle")

# defining objects 
save_path <- "/workspace/cfngle/results-data/02_conserved_seq/"
suffix <- ".fasta"


# Path to your .sam files
bam_files <- c("results-data/bowtie2/ZF_AC_CpG_1000bp_bt2_.bam", "results-data/bowtie2/ZF_AC_CpG_1000bp_bt2_.bam", "results-data/bowtie2/ZF_EH_CpG_1000bp_bt2_.bam")

# Function to read sequences from a .sam file and extract identifiers
read_bam_sequences <- function(bam_file) {
  header <- scanBamHeader(bam_file)
  param <- ScanBamParam(what="rname")  # qname for query names (identifiers)
  pos <-ScanBamParam(what="pos") 
  seqs <- scanBam(bam_file, param=param)
  unlist(lapply(seqs, function(x) x$rname))
}

# Extracting identifiers from each .sam file
identifiers_list <- lapply(bam_files, read_bam_sequences)

# Find overlapping identifiers between the first two files
overlaps <- Reduce(intersect, identifiers_list)

# Print overlapping identifiers
print(unique(overlaps))


#### testing MSA in R ####

mySequences <- readDNAStringSet("results-data/02_conserved_seq/AC_AS_EH_1000_conserved.fasta", format = "fasta")
msaResult <- msa(mySequences)

# You can print or plot the alignment result
msaResult

# msaPrettyPrint(msaResult) # needs more parameters set and an output for the PDF file

#### testing plotting data ####

# BiocManager::install("BSgenome.Drerio.UCSC.danRer11")
library(BSgenome.Drerio.UCSC.danRer11)
Drerio
genome_rerio <- readDNAStringSet("raw-data/ZF/rgenome/GCF_000002035.6_GRCz11_genomic.fna")
names(genome_rerio) <- gsub(" .*", "", names(genome_rerio))
genome_rerio_gr <- GRanges(
  seqnames = Rle(genome_rerio@ranges@NAMES),
  ranges = IRanges(c(start = genome_rerio@ranges@start), end = c(genome_rerio@ranges@start + genome_rerio@ranges@width), names = genome_rerio@ranges@NAMES),
)

plot_test <- granges(overlap_ZF_EH_AC_AS)
plot_test_24 <- plot_test[seqnames(plot_test) %in% genome_rerio_gr@ranges@NAMES[1:24]]
# plot_test@strand <- Rle(strand(rep("*", length(plot_test))))


# chromosomeNames = (paste0("chr", as.character(seq(1,24,1))))

plot_1 <- plotKaryotype(genome_rerio_gr[genome_rerio_gr@ranges@NAMES[1:24]])

# kpAddBaseNumbers(plot_1, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=0.5, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "grey")

kpPlotMarkers(plot_1, chr=seqnames(plot_test), x=plot_test_24@ranges@start, labels = as.character(1:length(plot_test_24)), cex=0.7, text.orientation = "horizontal", line.color = "red")

# create costum df for plotting 
df_plot <- data.frame(chr=seqnames(plot_test_24), start=plot_test_24@ranges@start, end=plot_test_24@ranges@start+plot_test_24@ranges@width)

kpPlotRegions(plot_1, data=df_plot, col = "red", r0=0.3, layer.margin = 0.6)
kpPlotDensity(plot_1, df_plot)

# kpLines(plot_1, chr=seqnames(plot_test), x=plot_test@ranges@start, y=1, col="#440000", lwd=1.5)
kpPoints(plot_1, chr=seqnames(plot_test), x=plot_test@ranges@start, y=0, cex=3)


test_regions <- createRandomRegions(nregions=10, length.mean = 1e6, mask=NA, non.overlapping = FALSE)




print(autoplot("danRer11"))

autoplot(genome_rerio_gr[genome_rerio_gr@ranges@NAMES == "NC_007113.7 Danio rerio strain Tuebingen chromosome 2, GRCz11 Primary Assembly"], which = granges(overlap_ZF_EH_AC_AS), layout= "karyogram")
      