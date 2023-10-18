"""
The purpose of this script is to generate fastq files containing only the sequences which are referenced in .txt files.

The .txt files represent collections of reads proposed to align with the M1 and M2 tuatara genomes.

"""

import gzip
from Bio import AlignIO

class fastq_generation:
    """Opens a given fastq.gz file and packs the contents into a dictionary.

    Input: .txt file, source .fastq.gz file
    Output: .fastq file containing only sequences found in a given mt genome
    """

    def __init__(self):
        self.read_dict = {}
        self.read_composition_dict = {}
        self.light_strand_dict = {}

    def fastq_file_parse(self):
        """In parsing the fasta files, the bait sequences must first be assessed.

        Afterwards, the read files are unpacked sequentially.
        """

        read_file_name = "SRR8351024.fastq.gz"

        print("File Name: ", read_file_name)
        print(" ")

        reads = gzip.open(read_file_name, "r")

        read_list = []
        for read in reads:
            read = str(read).replace("\n", "")
            read = read.strip("b").replace("\n", "")
            read_list.append(read)


        for i in range(len(read_list)-1):
            if "@" in read_list[i]:
                self.read_dict[read_list[i]] = read_list[i+1]

                print(read_list[i])

    def text_file_parse(self, infile):
        """Opens text file and unpacks items"""


       #new_file = open("M2_seqs.txt", "w")

        file_list = ["M1.txt", "M2.txt"]
        new_file_list = ["M1_seqs.fasta", "M2_seqs.fasta"]


        for i in range(len(file_list)):

            new_file = open(new_file_list[i], "w")
            with open(file_list[i]) as text:
                for line in text:
                    trim_index = line.find("M")
                    read_header = line[:trim_index-1]

                    for key in self.read_dict.keys():
                        if read_header in key:

                            new_file.write('\n')
                            new_file.write(">"+key.replace("'",""))
                            new_file.write('\n')
                            new_file.write(self.read_dict[key].replace("'",""))



    def driver(self):
        """"""

        self.fastq_file_parse()
        self.text_file_parse("M1.txt")


def main():
    class_access = fastq_generation()
    class_access.driver()

if __name__ == '__main__':
    main()