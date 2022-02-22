# Weektaak 4 - Groep 6
# Aminozuur


from protein import *

def get_sequences_from_file(file_name):
    """
    Reads content of a FASTA file (containing multiple protein sequence),
    to return a list of ProteinSequence objects

    Input:
    file_name : string
        path the file is located

    Output:
    sequences : [ProteinSequence]
        List of ProteinSequence objects
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
                sequences.append(ProteinSequence(header, sequence.upper()))
                sequence = ""
        else:
            sequence += line[:-1]


    sequences.append(ProteinSequence(header, sequence.upper()))

    file.close()

    return sequences

def print_protein_count(protein_sequence):

    for aa in aa_data:
        letter = aa.get_letter()
        print(f"{letter} : {protein_sequence.percentage_by_letter(letter)}")

if __name__ == '__main__':

    fumh_path = "sequences/FUMH_HUMAN_protein.fasta"

    fumh_sequence = get_sequences_from_file(fumh_path)[0]

    print_protein_count(fumh_sequence)
    print(fumh_sequence.percentage_by_class(AAClass.inbetween))


