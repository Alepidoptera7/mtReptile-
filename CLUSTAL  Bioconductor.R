#ensure that the msa bioconductor package is installed. 

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

system.file("tex", "texshade.sty", package = "msa")

BiocManager::install("msa")
BiocManager::install("Biostrings")

library(msa)
library(Biostrings)

#check that requsite files extis
file.exists("C:/Users/Quin The Conquoror!/PycharmProjects/Tuatara_mTs/SRR8351024.fastq.gz")
file.exists("C:/Users/Quin The Conquoror!/Desktop/M2_seqs.fa")

#access the file paths

ref_path <- file.path("C:/Users/Quin The Conquoror!/Desktop/M1_seqs.fa")#"C:/Users/Quin The Conquoror!/PycharmProjects/Tuatara_mTs/SRR8351024.fastq.gz")
ref_path

fa_path <- file.path("C:/Users/Quin The Conquoror!/Desktop/M2_seqs.fa")
fa_path

#develop DNAStringSets for the ref and seqs
seqs <- readDNAStringSet(fa_path)
seqs

ref <- readDNAStringSet(ref_path)
ref

#now, to unpack the dnastringset one read at a time, developing a new fasta file. 

for(i in 1:length((seqs))) names(seqs) <- paste0("M2_seqs_", seq(i))

#write a fasta file for each seq, paired with the reference genome for CLUSTAL

for(s in names(seqs)) writeXStringSet(seqs[s], paste0(s,".fa"))


#performing CLUSTAL on each seq against the reference genome.

alignment <- msa(names(seqs))

alignment 

print(alignment, show="complete")

msaPrettyPrint(alignment, output="pdf", showNames = "True")
}

