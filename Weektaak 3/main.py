# Tutor weektaak 3 - groep 6
# Codon bias

from cds import CDS
from plotter import Plotter
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
                sequences.append(CDS(header, sequence.upper()))
                sequence = ""
        else:
            sequence += line[:-1]


    sequences.append(CDS(header, sequence.upper()))

    file.close()

    return sequences

def get_env_cds(cds_objects):

    for cds in cds_objects:
        if cds.is_env_protein():
            return cds

def combine_cds_sequences(cds_objects):

    new_sequence = ""
    for cds in cds_objects:
        new_sequence += cds.get_sequence()

    return CDS("[gene=combined]", new_sequence)


def print_codon_frequency(cds_objects):

    """
    print codon frequencies from a list of CDS objects
    """

    for seq in cds_objects:
        print("-" * 150)
        print(f"Gene: {seq.get_gene_name()}")
        print()

        for codon_freq in seq.get_codon_frequencies():
            print(codon_freq)

        print()


def plot_all_charts_for_cds(plt, cds, organism, gene):

    #plt.plot_stacked_bar(cds.get_codon_frequencies(), organism, gene)
    #plt.plot_stacked_bar(cds.get_codon_frequencies(), organism, gene, True)
    plt.plot_nested_pie_from_dict(cds.get_codon_frequencies(), organism, gene)


if __name__ == '__main__':

    # Paths to sequences
    human_path = "CDS sequenties/Human FH.fasta"
    plant_path = "CDS sequenties/Arabidopsis thaliana FUM1.fasta"
    bact1_path = "CDS sequenties/Treponema caldarium FUM.fasta"
    bact2_path = "CDS sequenties/Yersinia pestis FumA.fasta"

    hiv1_path = "CDS sequenties/HIV-1 coding sequence nucleotide.txt"
    hiv2_path = "CDS sequenties/HIV-2 coding sequence nucleotide.txt"
    SIV1_path = "CDS sequenties/SIV-1 coding sequence nucleotide.txt"
    SIV2_path = "CDS sequenties/SIV-2 coding sequence nucleotide.txt"

    # Get sequences from files
    human_coding_sequences = get_sequences_from_file(human_path)
    plant_coding_sequences = get_sequences_from_file(plant_path)
    bact1_coding_sequences = get_sequences_from_file(bact1_path)
    bact2_coding_sequences = get_sequences_from_file(bact2_path)

    hiv1_coding_sequences = get_sequences_from_file(hiv1_path)
    hiv2_coding_sequences = get_sequences_from_file(hiv2_path)
    siv1_coding_sequences = get_sequences_from_file(SIV1_path)
    siv2_coding_sequences = get_sequences_from_file(SIV2_path)

    plt = Plotter()

    # plt.plot_stacked_bar(human_coding_sequences[0].get_codon_frequencies(), "Human", "FUM")
    # plt.plot_stacked_bar(human_coding_sequences[0].get_codon_frequencies(), "Human", "FUM", True)
    # plt.plot_nested_pie_from_dict(human_coding_sequences[0].get_codon_frequencies(), "Human", "FUM")

    # Plot human
    for cds in human_coding_sequences:
        plot_all_charts_for_cds(plt, cds, "Human", cds.get_gene_name())

    # Plot plant
    for cds in plant_coding_sequences:
        plot_all_charts_for_cds(plt, cds, "Arabidopsis thaliana", cds.get_gene_name())

    # Plot bact 1
    for cds in bact1_coding_sequences:
        plot_all_charts_for_cds(plt, cds, "Treponema caldarium", cds.get_gene_name())

    # Plot bact 2
    for cds in bact2_coding_sequences:
        plot_all_charts_for_cds(plt, cds, "Yersinia pestis", cds.get_gene_name())

    # Plot HIV-1
    env_hiv1 = get_env_cds(hiv1_coding_sequences)
    viral_hiv1 = combine_cds_sequences(hiv1_coding_sequences)
    plot_all_charts_for_cds(plt, env_hiv1, "HIV-1", env_hiv1.get_gene_name())
    plot_all_charts_for_cds(plt, viral_hiv1, "HIV-1", "Viral")

    # Plot HIV-2
    env_hiv2 = get_env_cds(hiv2_coding_sequences)
    viral_hiv2 = combine_cds_sequences(hiv2_coding_sequences)
    plot_all_charts_for_cds(plt, env_hiv2, "HIV-2", env_hiv2.get_gene_name())
    plot_all_charts_for_cds(plt, viral_hiv2, "HIV-2", "Viral")

    # Plot SIV-1
    env_siv1 = get_env_cds(siv1_coding_sequences)
    viral_siv1 = combine_cds_sequences(siv1_coding_sequences)
    plot_all_charts_for_cds(plt, env_siv1, "SIV-1", env_siv1.get_gene_name())
    plot_all_charts_for_cds(plt, viral_siv1, "SIV-1", "Viral")

    # Plot SIV-2
    env_siv2 = get_env_cds(siv2_coding_sequences)
    viral_siv2 = combine_cds_sequences(siv2_coding_sequences)
    plot_all_charts_for_cds(plt, env_siv2, "SIV-2", env_siv2.get_gene_name())
    plot_all_charts_for_cds(plt, viral_siv2, "SIV-2", "Viral")


