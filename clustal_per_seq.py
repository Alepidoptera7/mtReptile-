
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

class clustal_per_seq:

    def __init__(self):
        self.read_dict = {}

    def fasta_reader(self):

        read_list = []
        with open("M2_seqs.fasta") as fasta:
            for line in fasta:
                read_list = line.split("\\n")

        for i in range(len(read_list)):

            if ">" in read_list[i]:
                self.read_dict[read_list[i]] = read_list[i+1]

    def bio_python(self):

        infile = "M2_seqs.fasta"
        outfile = "clustalw.aln"
        cline = ClustalwCommandline("clustalw2", infile=infile, outfile=outfile)
        print(cline)

        format = "clustal"
        align = AlignIO.read(outfile, format)

        print(align)

    def driver(self):
        self.fasta_reader()
        self.bio_python()

def main():
    class_access = clustal_per_seq()
    class_access.driver()

if __name__ == "__main__":
    main()