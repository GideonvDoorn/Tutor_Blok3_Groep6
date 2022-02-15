# Tutor weektaak 3 - groep 6
# Codon bias

from cds import cds
from translation import codon_dict

def get_sequences_from_file(file_name):
    """
    Reads content of a FASTA file (containing multiple sequence),
    to return a list of cds objects

    Input:
    file_name : string
        path the file is located

    Output:
    sequences : [cds]
        List of cds objects
    """

    # Open and read file, catch FileNotFoundError
    try:
        file = open(file_name, "r")
    except FileNotFoundError:
        print("file not found!")
        return "", ""

    sequences = []

    sequence = ""
    header = ""

    for line in file:

        if line.startswith(">"):
            header = line
            if sequence != "":
                sequences.append(cds(header, sequence))
                sequence = ""
        else:
            sequence += line[:-1]


    sequences.append(cds(header, sequence))

    file.close()

    return sequences

def print_codon_frequency(cds_objects):

    for seq in cds_objects:
        print("-" * 150)
        print(seq.get_gene_name())
        print()

        for codon_freq in seq.get_codon_frequencies():
            print(codon_freq)

        print()


if __name__ == '__main__':

    # Paths to sequences
    hiv1_path = "CDS sequenties/HIV-1 coding sequence nucleotide.txt"
    hiv1_coding_sequences = get_sequences_from_file(hiv1_path)

    print_codon_frequency(hiv1_coding_sequences)
