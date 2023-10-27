from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

class clustal_per_seq:

    def __init__(self):
        self.read_dict = {}

    def fasta_splitter(self):
        """Splits fasta file into one fasta file for each read sequence, containing only that read sequence and the reference genome."""

        read_list = []
        with open("M2_seqs.fa") as fasta:
            for line in fasta:
                read_list = line.split("\n")

        print(len(read_list))

        for i in range(len(read_list)):

            if ">" in read_list[i]:

                self.read_dict[read_list[i]] = read_list[i+1]

                new_file_name = read_list[i].replace(">", "").replace("@", "") + ".fa"
                print(new_file_name)
                new_file = open(new_file_name, "w")
                new_file.write(read_list[i])
                new_file.write(read_list[i+1])

        #infile = "M2_seqs.fasta"


        #outfile = "clustalw"+ seq_header + ".aln"



        #self.aligner(infile, outfile)

    def aligner(self, infile, outfile):
        cline = ClustalwCommandline("clustalw2", infile=infile, outfile=outfile)
        print(cline)

        format = "clustal"
        align = AlignIO.read(outfile, format)

        print(align)

    def driver(self):
        #self.fasta_reader()
        self.fasta_splitter()
        print("OK")

def main():
    class_access = clustal_per_seq()
    class_access.driver()
    print("OK")

if __name__ == "__main__":
    main()