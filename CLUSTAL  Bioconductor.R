#ensure that the msa bioconductor package is installed. 
if (!requireNamespace("BiocManager", quietly=TRUE)){
  install.packages("BiocManager")
  library(Biostrings)
}

if(!require('msa')){
  install.packages('msa')
  library('msa')
}

system.file("tex", "texshade.sty", package="msa")

BiocManager::install("msa")
BiocManager::install("Biostrings")

#check that requsite files extis
#for m1 use 30, for m2 use 29
file.exists("C:/Users/Quin The Conquoror!/Desktop/MN864229.fa")
file.exists("C:/Users/Quin The Conquoror!/Desktop/M2_seqs.fa")

#access the file paths

ref_path <- file.path("C:/Users/Quin The Conquoror!/Desktop/MN864229.fa")#"C:/Users/Quin The Conquoror!/PycharmProjects/Tuatara_mTs/SRR8351024.fastq.gz")
ref_path

fa_path <- file.path("C:/Users/Quin The Conquoror!/Desktop/M2_seqs.fa")
fa_path

#develop DNAStringSets for the ref and seqs
seqs <- readDNAStringSet(fa_path, format="fasta")
seqs

ref <- readDNAStringSet(ref_path, format="fasta")
ref

#now, to unpack the dnastringset one read at a time, developing a new fasta file. 
for(i in 1:length((seqs))) names(seqs) <- paste0("M2_seqs_", seq(i))

#write a fasta file for each seq, paired with the reference genome for CLUSTAL
for(s in names(seqs)) writeXStringSet(seqs[s], paste0(s, ".fa"))

#performing CLUSTAL on each seq against the reference genome.
for(s in names(seqs)){
alignment <- msa(seqs[s], ref, method="ClustalOmega")
alignment 
}

print(alignment, show="complete")
msaPrettyPrint(alignment, output="tex", showNames = "none", verbose=FALSE, askForOverwrite = FALSE, showLogo="none")
tools::texi2pdf("alignment.tex", texinputs=system.file("tex", package="msa"),clean=TRUE)


