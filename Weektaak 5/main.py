import random

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
            if sequence != "":
                sequences.append(
                    sequence.upper())
                sequence = ""
            header = line
        else:
            sequence += line[:-1]

    sequences.append(
        sequence.upper())

    file.close()

    return sequences


if __name__ == '__main__':
    # HIV/SIV genome paths
    hiv1_path = "sequences/HIV-1 coding sequence protein.txt"
    hiv2_path = "sequences/HIV-2 coding sequence protein.txt"
    siv1_path = "sequences/SIV-1 coding sequence protein.txt"
    siv2_path = "sequences/SIV-2 coding sequence protein.txt"

    hiv1_sequence = get_sequences_from_file(hiv1_path)[0]
    hiv2_sequence = get_sequences_from_file(hiv1_path)[0]
    siv1_sequence = get_sequences_from_file(hiv1_path)[0]
    siv2_sequence = get_sequences_from_file(hiv1_path)[0]

    print(hiv1_sequence)

    random.sample(hiv1_sequence)


