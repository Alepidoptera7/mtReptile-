import sys
import os
import gzip
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

class strand_inspection:
    """"""

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


    def composition_measurer(self):

        #print("G, A, T, C composition")
        for key in self.read_dict.keys():

            read = self.read_dict[key]
            read_len = len(read)
            self.read_composition_dict[key] = (read.count("G")/read_len * 100,  read.count("A")/read_len * 100, read.count("T")/read_len * 100, read.count("C")/read_len * 100)

            #print(key, self.read_composition_dict[key])

    def reverse_complement(self, seq):
        """Develops the reverse complement of a parameter sequence."""

        return seq.lower().replace('a', 'T').replace('t', 'A').replace('g', 'C').replace('c', 'G').replace(' ', '')[::-1]

    def heavy_to_light_conversion(self):
        key_count = 0
        for key in self.read_composition_dict.keys():
            read = self.read_dict[key]

            if self.read_composition_dict[key][0] > 14:
                self.light_strand_dict[key] = self.reverse_complement(read)

            else:
                self.light_strand_dict[key] = read

            key_count += 1

            #if key_count > 745:
            print("Read Header: >" + key)
            #    print("Sequence: ", self.light_strand_dict[key].strip("\n"))
            #    self.BLAST(self.light_strand_dict[key])

    def BLAST(self, seq):

        """Develops a BLAST request for each clustered sequence via BioPython.

           Design taken from Biopython Manual.
           """
        E_VALUE_THRESH = 0.0000000000000000001

        result_handle = NCBIWWW.qblast("blastn", "nt", sequence=seq, expect=1, alignments=3, perc_ident=95)

        blast_record = NCBIXML.read(result_handle)

        print(" ")
        print("****Alignments****")

        count = 0
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    count += 1
                    print("Sequence title:", alignment.title)
                    print("Alignment length: ", alignment.length)
                    print('Percent Identity: %', round(hsp.identities / hsp.align_length * 100, 2))
                    #print(hsp.query[0:75] + "...")
                    #print(hsp.match[0:75] + "...")
                    #print(hsp.sbjct[0:75] + "...")
                    #print(f"e-value: {hsp.expect}")

        print(" ")
        print("Number of alignments: ", count)
        print("=========================================================== ")
        print(" ")

    def driver(self):

        self.fastq_file_parse()
        self.composition_measurer()
        self.heavy_to_light_conversion()

def main():
    class_access = strand_inspection()
    class_access.driver()

if __name__ == '__main__':
    main()