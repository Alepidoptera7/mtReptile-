from Bio.Align.Applications import ClustalOmegaCommandline
import subprocess
import argparse

class clustal_per_seq:

    """"This class is designed to develop pairwise CLUSTAL Omega alignments
     between any number of sequences in a given FASTA file and a reference genome in FASTA format.

    A file will be created to hold each CLUSTAL Omega alignment.

    CLI example, "clustal_per_seq.py -f M2_seqs.fa -r MN864229.fa -o MN864229-M2_seqs"

    """

    def __init__(self, fasta, reference, path, outfile):
        self.read_dict = {}
        self.file_dict = {}
        self.alignment_fasta = fasta
        self.reference_genome = reference
        self.path = path
        self.outfile = outfile
        self.aln_outfile_list = []
        self.percent_id_dict = {}
        self.percent_id_list = []

    def fasta_splitter(self):
        """Splits fasta file into one fasta file for each read sequence, containing only that read sequence and the reference genome."""

        read_list, reference_list = [], []
        with open(self.alignment_fasta) as fasta:
            for line in fasta:
                read_list.append(line)

        with open(self.reference_genome) as reference_fasta:
            for line in reference_fasta:
                reference_list.append(line)
                #print(len(line))
        print(len(read_list), len(reference_list))

        count_seqs = 0
        for i in range(len(read_list)):
            if ">" in read_list[i]:
                count_seqs += 1
                self.read_dict[read_list[i]] = read_list[i+1]
        print("number of sequences: ", count_seqs)
        alignment_count = 0

        for key in self.read_dict.keys():
            new_file_name = key[2:-14] +".fa"
            new_file = open(new_file_name, "w")
            new_file.write(key)
            new_file.write(self.read_dict[key])
            new_file.write("\n")

            for line in reference_list:
                new_file.write(line)

            cline = self.biopython_clustalw(new_file_name)
            self.sub_process(cline)
            alignment_count += 1
            print(cline)
        print("number of alignments: ", alignment_count)

    def biopython_clustalw(self, infile):
        # r"C:/Users/Quin The Conquoror!/Desktop/clustal-omega-1.2.2-win64/clustalo"
        clustalOmega_exe = r"C:/Users/Quin The Conquoror!/Desktop/clustal-omega-1.2.2-win64/clustalo"
        # "M1-MN864230.1_clustalOmega_aligned"
        cline_outfile = infile + self.outfile
        print(cline_outfile)
        self.aln_outfile_list.append(cline_outfile)

        cline = ClustalOmegaCommandline(clustalOmega_exe, infile=infile, verbose=True, outfile=cline_outfile, outfmt="fasta", percentid=True)

        #powershell_cline = "./"+str(cline)
        return str(cline)

    def sub_process(self, cline):
        subprocess.run(cline)

    def percentid_calculator(self):
        """The purpose of this def is to caculate percent identity between
         a given sequence and the reference genome from previously generated alignments.

         Input: previously generated .fasta alignment file
         Output: percent identy per alignment
         """

        base_list = ["N", "A", "G", "C", "T", "n"]

        for aln_outfile in self.aln_outfile_list:
            alignment_string = ""
            base_index_list = []
            alignment_string_list = []
            gap_count = 0
            with open(aln_outfile) as aln:
                for line in aln:
                    alignment_string_list.append(line)

            for line in alignment_string_list[1:]:
                alignment_string += line

            for base in base_list:
                base_index_list.append(alignment_string.index(base))

            base_index_list.sort()

            for base in alignment_string[base_index_list[0]:base_index_list[-1]]:
                if base == "-":
                    gap_count += 1

            aln_len = base_index_list[-1] - base_index_list[0]
            match_count = aln_len-gap_count

            self.percent_id_dict[aln_outfile] = ((match_count * 100)/aln_len)

            percent_identity = ((match_count * 100)/aln_len)

            self.percent_id_list.append((aln_outfile, percent_identity))

            #print(aln_outfile)
            #print("base index list", base_index_list)
            #print("gap count", gap_count)
            #print("alignment length", aln_len)
            #print("match count", match_count)
            #print("Percent Identity: ", ((match_count * 100)/aln_len))

    def driver(self):
        #self.fasta_reader()
        self.fasta_splitter()
        self.percentid_calculator()

        self.percent_id_list.sort(key=lambda a: a[1])
        print(self.percent_id_list)

class CommandLine:
    """
    Allows use of command line for program inputs.
    Input: command line rendered variables.
    Output: values assigned to variables.
    """

    def __init__(self, inOpts = None):

        self.parser = argparse.ArgumentParser()

        self.parser.add_argument('-f', '--fasta', type=str,
                                 action='store', default="", help='alignment fasta')
        self.parser.add_argument('-r', '--reference', type=str,
                                 action='store', default="", help='reference genome fasta')
        self.parser.add_argument('-p', '--path', type=str,
                                 action='store', default="", help='path of CLUSTAL variant')
        self.parser.add_argument('-o', '--outfile', type=str,
                                 action='store', default="", help='name of outfile')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


class Usage(Exception):
    """
    Used to signal a Usage error, evoking a usage
    statement and eventual exit when raised.
    """

    def __init__(self, msg):
        self.msg = msg


def main(myCommandLine=None):

    if myCommandLine == None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(myCommandLine)
    try:
        print(myCommandLine.args)

    except Usage as err:
        print(err.msg)

    class_access = clustal_per_seq(myCommandLine.args.fasta, myCommandLine.args.reference,
                                   myCommandLine.args.path, myCommandLine.args.outfile)
    class_access.driver()
    print("OK")

if __name__ == "__main__":
    main()