# Tutor weektaak 2 - groep 6
# GC% analyseren

import matplotlib.pyplot as plt
import numpy as np

def get_sequence_from_file(file_name):
    """
    Reads content of a FASTA file (containing one sequence),
    to return a header and the isolated sequence

    Input:
    file_name : string
        path the file is located

    Output:
    header : string
        isolated header
    sequence : string
        sequence without headers
    """

    # Open and read file, catch FileNotFoundError
    try:
        file = open(file_name, "r")
    except FileNotFoundError:
        print("file not found!")
        return "", ""

    header = ""
    sequence = ""

    # Separate headers , from sequence
    for line in file:
        if line.startswith(">"):
            header = line
        else:
            sequence += line
    file.close()

    return header, sequence


def get_gc_per_step(sequence, step):

    gc_step_percentages = []
    count = 0
    gc_count = 0


    for char in sequence:

        if char == "G" or char == "C":
            gc_count += 1

        if count == step:
            gc_step_percentages.append(gc_count / step)
            gc_count = 0
            count = 0
        else:
            count += 1

    return gc_step_percentages

def plot_charts(gc_step_percentages, step):

    # Create mean line values
    mean_line = [np.mean(gc_step_percentages)] * len(gc_step_percentages)
    # Create median line values
    median_line = [np.median(gc_step_percentages)] * len(gc_step_percentages)


    # Plot line chart
    plt.plot(range(len(gc_step_percentages)), gc_step_percentages)
    plt.title("GC Percentage by steps of: " + str(step) + " bp")
    plt.xlabel("Step")
    plt.ylabel("GC Percentage")
    plt.ylim([min(gc_step_percentages) - 0.015, max(gc_step_percentages) + 0.015])
    # Plot mean line
    plt.plot(range(len(gc_step_percentages)), mean_line)
    # Plot median line
    plt.plot(range(len(gc_step_percentages)), median_line)

    plt.show()

    # Plot bar chart
    fig = plt.figure(figsize=(10, 5))
    plt.bar(range(len(gc_step_percentages)), gc_step_percentages)
    plt.xlabel("Step")
    plt.ylabel("GC Percentage")
    plt.title("GC Percentage by steps of: " + str(step) + " bp")
    plt.ylim([min(gc_step_percentages) - 0.015, max(gc_step_percentages) + 0.015])
    # Plot mean line
    plt.plot(range(len(gc_step_percentages)), mean_line, color="orange")

    plt.show()

if __name__ == '__main__':
    human_chromosome = "Genomes/chr17_human.fasta"
    plant_chromosome = "Genomes/Chr1_ChloropiconPrimus.fasta"
    syphilis_genome = "Genomes/Treponema pallidum full sequence.fasta"
    plague_genome = "Genomes/Yersinia pestis full sequence.fasta"

    step = 50000

    print(get_sequence_from_file(syphilis_genome)[0])
    gc_step_percentages = get_gc_per_step(get_sequence_from_file(syphilis_genome)[1], step)

    for count in gc_step_percentages:
        print(count)

    plot_charts(gc_step_percentages, step)


