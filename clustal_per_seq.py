import Bio
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import pairwise2

class clustal_per_seq:

    def __init__(self):
        self.read_dict = {}
        self.file_dict = {}

    def fasta_splitter(self):
        """Splits fasta file into one fasta file for each read sequence, containing only that read sequence and the reference genome."""

        read_list, reference_list = [], []
        with open("M2_seqs.fa") as fasta:
            for line in fasta:
                read_list.append(line)

        with open("MN864229.fa") as reference_fasta:
            for line in reference_fasta:
                reference_list.append(line)
                #print(len(line))
        print(len(read_list), len(reference_list))

        for i in range(len(read_list)):
            if ">" in read_list[i]:
                self.read_dict[read_list[i]] = read_list[i+1]

        for key in self.read_dict.keys():
            new_file_name = key[2:-14] +".fa"
            new_file = open(new_file_name, "w")
            new_file.write(key)
            new_file.write(self.read_dict[key])
            new_file.write("\n")

            for line in reference_list:
                new_file.write(line)

            self.biopython_clustalw(new_file_name)

    def biopython_clustalw(self, infile):
        clustalOmega_exe = r"C:/Users/Quin The Conquoror!/Desktop/clustal-omega-1.2.2-win64/clustalo"
        outfile = infile + "_aligned"
        cline = ClustalOmegaCommandline(clustalOmega_exe, infile=infile, verbose=True, outfile=outfile, outfmt="fasta")
        print("./"+str(cline))


    def driver(self):
        #self.fasta_reader()
        self.fasta_splitter()

def main():
    class_access = clustal_per_seq()
    class_access.driver()
    print("OK")

if __name__ == "__main__":
    main()