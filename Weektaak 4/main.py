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

    print("-"*80)
    print(protein_sequence.get_header())

    for aa in aa_data:
        letter = aa.get_letter()
        print(f"{letter} : {protein_sequence.percentage_by_letter(letter)}")

    print(f"Hydrophobic: {protein_sequence.percentage_by_class(AAClass.hydrophobic)}%")
    print("")

if __name__ == '__main__':

    fumh_path = "sequences/FUMH_HUMAN_protein.fasta"
    gtr8_path = "sequences/GTR8_HUMAN_protein.fasta"
    hrh3_path = "sequences/HRH3_HUMAN_protein.fasta"
    ntrk1_path = "sequences/NTRK1_HUMAN_protein.fasta"

    hiv1_path = "sequences/HIV-1 coding sequence protein.txt"
    hiv2_path = "sequences/HIV-2 coding sequence protein.txt"
    siv1_path = "sequences/SIV-1 coding sequence protein.txt"
    siv2_path = "sequences/SIV-2 coding sequence protein.txt"

    fumh_sequence = get_sequences_from_file(fumh_path)[0]
    gtr8_sequence = get_sequences_from_file(gtr8_path)[0]
    hrh3_sequence = get_sequences_from_file(hrh3_path)[0]
    ntrk1_sequence = get_sequences_from_file(ntrk1_path)[0]

    print_protein_count(fumh_sequence)
    print_protein_count(gtr8_sequence)
    print_protein_count(hrh3_sequence)
    print_protein_count(ntrk1_sequence)



