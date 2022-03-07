# Weektaak 4 - Groep 6
# Aminozuur

from protein import *
from plotter import Plotter


def get_sequences_from_file(file_name, organism, gene):
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
                    ProteinSequence(header, sequence.upper(), organism, gene))
                sequence = ""
            header = line
        else:
            sequence += line[:-1]

    sequences.append(ProteinSequence(header, sequence.upper(), organism, gene))

    file.close()

    return sequences


def get_env_ps(protein_sequences):
    for ps in protein_sequences:
        if ps.is_env_protein():
            return ps


def combine_ps_sequences(protein_sequences, organism):
    new_sequence = ""
    for ps in protein_sequences:
        if not ps.is_env_protein():
            new_sequence += ps.get_sequence()

    return ProteinSequence("[gene=viral]", new_sequence, organism, "Viral")


def print_protein_count(protein_sequence, organism):
    print("-" * 80)
    print(organism + "\n")
    print("AA:  %")
    print("--------")
    for aa in aa_data:
        letter = aa.get_letter()
        print(f"{letter} : {protein_sequence.percentage_by_letter(letter)}")

    print(f"\nMass : {protein_sequence.calculate_mass()}u")
    print(f"Charge : {protein_sequence.predicted_charge().name}")
    print("\n-Hydropathy-")
    print(
        f"Hydrophobic: {protein_sequence.percentage_by_class(AAClass.hydrophobic)}%")
    print(
        f"Hydrophilic: {protein_sequence.percentage_by_class(AAClass.hydrophilic)}%")
    print(
        f"In between: {protein_sequence.percentage_by_class(AAClass.in_between)}%")
    print("")


if __name__ == '__main__':
    # Human gene paths
    fumh_path = "sequences/FUMH_HUMAN_protein.fasta"
    gtr8_path = "sequences/GTR8_HUMAN_protein.fasta"
    hrh3_path = "sequences/HRH3_HUMAN_protein.fasta"
    ntrk1_path = "sequences/NTRK1_HUMAN_protein.fasta"

    # HIV/SIV genome paths
    hiv1_path = "sequences/HIV-1 coding sequence protein.txt"
    hiv2_path = "sequences/HIV-2 coding sequence protein.txt"
    siv1_path = "sequences/SIV-1 coding sequence protein.txt"
    siv2_path = "sequences/SIV-2 coding sequence protein.txt"

    # Retrieve sequences as objects
    fumh_sequence = get_sequences_from_file(fumh_path, "Human", "FUMH")[0]
    gtr8_sequence = get_sequences_from_file(gtr8_path, "Human", "GTR8")[0]
    hrh3_sequence = get_sequences_from_file(hrh3_path, "Human", "HRH3")[0]
    ntrk1_sequence = get_sequences_from_file(ntrk1_path, "Human", "NTRK1")[0]

    # Add domains to human proteins
    fumh_sequence.add_domain("Transit Peptide", 1, 44)
    fumh_sequence.add_domain("Rest Of Chain", 44, 510)

    gtr8_sequence.add_domain("Cytoplasmic", 1, 25)
    gtr8_sequence.add_domain("Transmembrane", 26, 46)
    gtr8_sequence.add_domain("Extracellular", 47, 70)
    gtr8_sequence.add_domain("Transmembrane", 71, 91)
    gtr8_sequence.add_domain("Cytoplasmic", 92, 96)
    gtr8_sequence.add_domain("Transmembrane", 97, 117)
    gtr8_sequence.add_domain("Extracellular", 118, 127)
    gtr8_sequence.add_domain("Transmembrane", 128, 148)
    gtr8_sequence.add_domain("Cytoplasmic", 149, 156)
    gtr8_sequence.add_domain("Transmembrane", 157, 177)
    gtr8_sequence.add_domain("Extracellular", 178, 182)
    gtr8_sequence.add_domain("Transmembrane", 183, 203)
    gtr8_sequence.add_domain("Cytoplasmic", 204, 256)
    gtr8_sequence.add_domain("Transmembrane", 257, 277)
    gtr8_sequence.add_domain("Extracellular", 278, 292)
    gtr8_sequence.add_domain("Transmembrane", 293, 313)
    gtr8_sequence.add_domain("Cytoplasmic", 314, 319)
    gtr8_sequence.add_domain("Transmembrane", 320, 340)
    gtr8_sequence.add_domain("Extracellular", 341, 367)
    gtr8_sequence.add_domain("Transmembrane", 368, 388)
    gtr8_sequence.add_domain("Cytoplasmic", 389, 404)
    gtr8_sequence.add_domain("Transmembrane", 405, 425)
    gtr8_sequence.add_domain("Extracellular", 426, 438)
    gtr8_sequence.add_domain("Transmembrane", 439, 459)
    gtr8_sequence.add_domain("Cytoplasmic", 460, 477)

    hrh3_sequence.add_domain("Cytoplasmic", 1, 39)
    hrh3_sequence.add_domain("Transmembrane", 40, 60)
    hrh3_sequence.add_domain("Extracellular", 61, 70)
    hrh3_sequence.add_domain("Transmembrane", 71, 91)
    hrh3_sequence.add_domain("Cytoplasmic", 92, 108)
    hrh3_sequence.add_domain("Transmembrane", 109, 129)
    hrh3_sequence.add_domain("Extracellular", 130, 156)
    hrh3_sequence.add_domain("Transmembrane", 157, 177)
    hrh3_sequence.add_domain("Cytoplasmic", 178, 196)
    hrh3_sequence.add_domain("Transmembrane", 197, 217)
    hrh3_sequence.add_domain("Extracellular", 218, 359)
    hrh3_sequence.add_domain("Transmembrane", 360, 380)
    hrh3_sequence.add_domain("Cytoplasmic", 381, 395)
    hrh3_sequence.add_domain("Transmembrane", 396, 416)
    hrh3_sequence.add_domain("Extracellular", 417, 445)

    ntrk1_sequence.add_domain("Signal Peptide", 1, 32)
    ntrk1_sequence.add_domain("Extracellular ", 33, 423)
    ntrk1_sequence.add_domain("Transmembrane", 424, 439)
    ntrk1_sequence.add_domain("Cytoplasmic ", 440, 796)


    hiv1_coding_sequences = get_sequences_from_file(hiv1_path, "hiv-1", "env")
    hiv2_coding_sequences = get_sequences_from_file(hiv2_path, "hiv-2", "env")
    siv1_coding_sequences = get_sequences_from_file(siv1_path, "siv-1", "env")
    siv2_coding_sequences = get_sequences_from_file(siv2_path, "siv-2", "env")

    # Separate env from viral genes
    env_hiv1 = get_env_ps(hiv1_coding_sequences)
    viral_hiv1 = combine_ps_sequences(hiv1_coding_sequences, "hiv-1")
    env_hiv2 = get_env_ps(hiv2_coding_sequences)
    viral_hiv2 = combine_ps_sequences(hiv2_coding_sequences, "hiv-2")
    env_siv1 = get_env_ps(siv1_coding_sequences)
    viral_siv1 = combine_ps_sequences(siv1_coding_sequences, "siv-1")
    env_siv2 = get_env_ps(siv2_coding_sequences)
    viral_siv2 = combine_ps_sequences(siv2_coding_sequences, "siv-2")

    # Print sequence data
    print_protein_count(fumh_sequence, "Human FUMH")
    print_protein_count(gtr8_sequence, "Human GTR8")
    print_protein_count(hrh3_sequence, "Human HRH3")
    print_protein_count(ntrk1_sequence, "Human NTRK1")

    print_protein_count(env_hiv1, "HIV-1 Env")
    print_protein_count(viral_hiv1, "HIV-1 Viral")
    print_protein_count(env_hiv2, "HIV-2 Env")
    print_protein_count(viral_hiv2, "HIV-2 Viral")
    print_protein_count(env_siv1, "SIV-1 Env")
    print_protein_count(viral_siv1, "SIV-1 Viral")
    print_protein_count(env_siv2, "SIV-2 Env")
    print_protein_count(viral_siv2, "SIV-2 Viral")

    all_sequences = [fumh_sequence, gtr8_sequence, hrh3_sequence,
                     ntrk1_sequence, env_hiv1, viral_hiv1, env_hiv2,
                     viral_hiv2, env_siv1, viral_siv1, env_siv2, viral_siv2]

    hiv_sequences = [env_hiv1, viral_hiv1, env_hiv2,
                     viral_hiv2, env_siv1, viral_siv1, env_siv2, viral_siv2]

    plt = Plotter()

    # plt.plot_percentage_by_domain_by_class(fumh_sequence, AAClass.hydrophilic)
    # plt.plot_percentage_by_domain_by_class(gtr8_sequence, AAClass.hydrophilic)
    # plt.plot_percentage_by_domain_by_class(hrh3_sequence, AAClass.hydrophilic)
    # plt.plot_percentage_by_domain_by_class(ntrk1_sequence, AAClass.hydrophilic)

    # plt.plot_hydropathy_percentage_by_domain(fumh_sequence)
    # plt.plot_hydropathy_percentage_by_domain(gtr8_sequence)
    # plt.plot_hydropathy_percentage_by_domain(hrh3_sequence)
    # plt.plot_hydropathy_percentage_by_domain(ntrk1_sequence)
    #
    for hiv_seq in hiv_sequences:
        plt.plot_three(hiv_seq)
        plt.plot_three(hiv_seq, False)

    # plt.plot_percentage_single_aa(hiv_sequences, "C")
    # plt.plot_percentage_single_aa(hiv_sequences, "W")
    # plt.plot_percentage_hydro(hiv_sequences)




    # for aa in aa_data:
    #
    #     plt.plot_percentage_single_aa(
    #         all_sequences, aa.get_letter())


    # plt.plot_percentage_class(
    #     all_sequences,
    #     AAClass.hydrophobic)
