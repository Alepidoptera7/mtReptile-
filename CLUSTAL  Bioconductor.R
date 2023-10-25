#ensure that the msa bioconductor package is installed. 

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

system.file("tex", "texshade.sty", package = "msa")

BiocManager::install("msa")
BiocManager::install("Biostrings")

library(msa)
library(Biostrings)

#now, to unpack the fasta file one sequence at a time, performing a CLUSTAL
#alignment on each seq against the reference genome. 

#fa_file <- system.file(package="Biostrings", "SRR8351024.fastq.gz", "extdata")

seqs <- readDNAStringSet("M2_seqs.fa", "extdata", format="fasta")

seqs

lines <- readLines(seqs, 1)

alignment <- msa(lines)

print(alignment, show="complete")



