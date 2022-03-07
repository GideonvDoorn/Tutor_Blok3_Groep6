import matplotlib.pyplot as plt
import numpy as np

from amino_acid import *

# Define amino acid label colors
aa_colors = {"Ala": "indianred", "Arg": "lavender", "Asn": "royalblue",
             "Asp": "blueviolet", "Cys": "darkmagenta",
             "Gln": "goldenrod", "Glu": "deeppink", "Gly": "salmon",
             "His": "olivedrab", "Ile": "peru",
             "Leu": "darkorchid", "Lys": "thistle", "Met": "darkseagreen",
             "Phe": "slateblue",
             "Pro": "aquamarine", "Ser": "dimgray", "Stop": "darkgreen",
             "Thr": "lightyellow", "Trp": "rosybrown",
             "Tyr": "cadetblue", "Val": "darkorange"}

class Plotter:

    def __init__(self):
        pass

    def __get_color_by_abr(self, abr):

        try:
            return aa_colors[abr]
        except KeyError:
            print("No such amino acid : Plotter.get_color_by_abr")


    def plot_percentage_single_aa(self, sequences, letter):

        categories = []
        letter_values = []

        for seq in sequences:
            categories.append(f"{seq.get_organism()}\n{seq.get_gene()}")
            letter_values.append(seq.percentage_by_letter(letter))

        # define figure
        fig, ax = plt.subplots(1)

        # plot bars
        x_values = np.arange(len(categories))
        plt.bar(categories, letter_values, label=categories, color=self.__get_color_by_abr(get_abr_by_letter(letter)))

        # Scale fontsize with amount of categories
        fontsize = 8 / (len(categories) /10)
        # bar labels
        for x in x_values:
            if letter_values[x] < 1:
                y_offset = 0
            else:
                y_offset = 0.5

            plt.text(x-0.22, letter_values[x] - y_offset, f"{round(letter_values[x], 1)}%", size=fontsize)

        # remove spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # grid lines
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray', linestyle='dashed', alpha=0.2)

        plt.ylabel(f"% {get_full_name_by_letter(letter)}")
        plt.title(f"Percentage of amino acid: {letter}")
        plt.tight_layout()
        plt.show()
        fig.savefig(f"Single_all {letter}", dpi=1000,
                    transparent=True)

    def plot_three(self, sequence, highest_occurrence=True):

        aa_frequencies = {}

        for aa in aa_data:

            aa_count = round(sequence.count_by_letter(aa.get_letter()) / sequence.get_length() * 100, 2)

            aa_frequencies[aa.get_abbreviation()] = aa_count

        aa_frequencies = dict(sorted(aa_frequencies.items(), key=lambda item: item[1]))

        if highest_occurrence:

            categories = list(aa_frequencies.keys())[-3:]
            values = list(aa_frequencies.values())[-3:]
        else:
            categories = list(aa_frequencies.keys())[:3]
            values = list(aa_frequencies.values())[:3]

        colors = []
        for i in range(len(categories)):
            colors.append(self.__get_color_by_abr(categories[i]))
            categories[i] += f" {values[i]}%"



        plt.pie(values, labels=categories, colors=colors, wedgeprops=dict(width=0.3,
                                                        edgecolor='w'))

        if highest_occurrence:
            plt.title(f"Three Most occurring amino acids in {sequence.get_gene()} {sequence.get_organism()}")
        else:
            plt.title(f"Three Least occurring amino acids in {sequence.get_gene()} {sequence.get_organism()}")
        centre_circle = plt.Circle((0, 0), 0.4, color="white", linewidth=0)
        fig = plt.gcf()
        fig.gca().add_artist(centre_circle)
        plt.axis('equal')
        fig.savefig(f"Three_{sequence.get_organism()} {sequence.get_gene()} {highest_occurrence}", dpi=1000, transparent=True)




        plt.show()


    def plot_percentage_class(self, sequences, class_):

        categories = []
        class_values = []

        for seq in sequences:
            categories.append(f"{seq.get_organism()}\n{seq.get_gene()}")
            class_values.append(seq.percentage_by_class(class_))

        # define figure
        fig, ax = plt.subplots(1)

        # plot bars
        x_values = np.arange(len(categories))

        plt.bar(categories, class_values, label=categories, color="coral")

        # Scale fontsize with amount of categories
        fontsize = 10 / (len(categories) /10)
        # bar labels
        for x in x_values:
            if class_values[x] < 4:
                y_offset = 0
            else:
                y_offset = 5

            plt.text(x-0.19, class_values[x] - y_offset, f"{round(class_values[x])}%", size=fontsize)

        # remove spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # grid lines
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray', linestyle='dashed', alpha=0.2)
        plt.ylabel(f"% {class_.name}")
        plt.title(f"Percentage of amino acid: {class_.name}")
        plt.tight_layout()
        plt.show()

        fig.savefig(f"Hydro_{sequences[0].get_gene()}", dpi=1000, transparent=True)

    def plot_percentage_hydro(self, sequences):

        categories = []
        philic_values = []
        phobic_values = []

        for seq in sequences:
            categories.append(f"{seq.get_organism()}\n{seq.get_gene()}")
            philic_values.append(seq.percentage_by_class(AAClass.hydrophilic))
            phobic_values.append(seq.percentage_by_class(AAClass.hydrophobic))

        # define figure
        fig, ax = plt.subplots(1)

        # plot bars
        x_values = np.arange(len(categories))
        plt.bar(x_values - 0.2, philic_values, 0.4, color='aquamarine')
        plt.bar(x_values + 0.2, phobic_values, 0.4, color='coral')

        # Scale fontsize with amount of categories

        fontsize = 6 / (len(categories) / 10)
        # bar labels
        for x in x_values:
            if philic_values[x] < 4:
                y_offset = 0
            else:
                y_offset = 3

            plt.text(x - 0.38, philic_values[x] - y_offset,
                     f"{round(philic_values[x])}%", size=fontsize)

            if phobic_values[x] < 4:
                y_offset = 0
            else:
                y_offset = 3

            plt.text(x + 0.03, phobic_values[x] - y_offset,
                     f"{round(phobic_values[x])}%", size=fontsize)
        # remove spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # grid lines
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray', linestyle='dashed', alpha=0.2)

        plt.xticks(x_values, categories)
        plt.ylabel(f"%")
        plt.title(f"Percentage of amino acid: Hydropathy")
        plt.tight_layout()
        plt.show()

        fig.savefig(f"Hydro_{sequences[0].get_gene()}", dpi=1000, transparent=True)

    def plot_percentage_by_domain_by_class(self, sequence, class_):

        class_percentages = {}
        percentages = sequence.domain_percentages_by_class(class_)

        domains = sequence.get_domains()
        for i in range(len(domains)):

            if domains[i].get_name() not in list(class_percentages.keys()):
                class_percentages[domains[i].get_name()] = [percentages[i]]
            else:
                class_percentages[domains[i].get_name()].append(percentages[i])

        categories = list(class_percentages.keys())
        class_values = []

        for category in categories:
            total = sum(class_percentages[category])
            class_values.append(total / len(class_percentages[category]))

        # define figure
        fig, ax = plt.subplots(1)
        plt.bar(categories, class_values, label=categories, color="coral")
        plt.ylabel(f"% {class_.name}")
        plt.title(f"Percentage {class_.name} by domain for {sequence.get_gene()}")
        plt.tight_layout()
        plt.show()
        fig.savefig(f"Hydro_{sequence.get_gene()}", dpi=1000,
                    transparent=True)

    def plot_hydropathy_percentage_by_domain(self, sequence):
        # Get hydrophilic values
        philic_dict = {}
        philic_percentages = sequence.domain_percentages_by_class(AAClass.hydrophilic)
        domains = sequence.get_domains()
        for i in range(len(domains)):
            if domains[i].get_name() not in list(philic_dict.keys()):
                philic_dict[domains[i].get_name()] = [philic_percentages[i]]
            else:
                philic_dict[domains[i].get_name()].append(
                    philic_percentages[i])
        categories = list(philic_dict.keys())
        philic_values = []
        for category in categories:
            total = sum(philic_dict[category])
            philic_values.append(total / len(philic_dict[category]))

        # Get hydrophobic values
        phobic_dict = {}
        phobic_percentages = sequence.domain_percentages_by_class(AAClass.hydrophobic)
        domains = sequence.get_domains()
        for i in range(len(domains)):
            if domains[i].get_name() not in list(phobic_dict.keys()):
                phobic_dict[domains[i].get_name()] = [phobic_percentages[i]]
            else:
                phobic_dict[domains[i].get_name()].append(
                    phobic_percentages[i])
        phobic_values = []
        for category in categories:
            total = sum(phobic_dict[category])
            phobic_values.append(total / len(phobic_dict[category]))

        # define figure
        fig, ax = plt.subplots(1)

        # plot bars
        x_values = np.arange(len(categories))
        plt.bar(x_values - 0.1, philic_values, 0.2, color='aquamarine')
        plt.bar(x_values + 0.1, phobic_values, 0.2, color='coral')

        # Plot average lines
        plt.axhline(y = sequence.percentage_by_class(AAClass.hydrophilic),
                    color = 'turquoise', linestyle = '-', alpha=0.5,label=None)
        plt.axhline(y=sequence.percentage_by_class(AAClass.hydrophobic),
                    color='lightcoral', linestyle='-', alpha=0.5,label=None)

        # Scale fontsize with amount of categories

        fontsize = 3 / (len(categories) /10)
        # bar labels
        for x in x_values:
            if philic_values[x] < 4:
                y_offset = 0
            else:
                y_offset = 3

            plt.text(x-0.17, philic_values[x] - y_offset, f"{round(philic_values[x])}%", size=fontsize)

            if phobic_values[x] < 4:
                y_offset = 0
            else:
                y_offset = 3

            plt.text(x+0.03, phobic_values[x] - y_offset, f"{round(phobic_values[x])}%", size=fontsize)


        # remove spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # grid lines
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray', linestyle='dashed', alpha=0.2)

        # Info
        plt.xticks(x_values, categories)
        plt.ylabel(f"%")
        plt.title(f"Percentage class by domain for {sequence.get_gene()}")

        plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True, labels=["Mean Hydrophilic", "Mean Hydrophobic", "Hydrophilic", "Hydrophobic"])
        plt.tight_layout()
        plt.show()
        fig.savefig(f"Hydro_{sequence.get_gene()}", dpi=1000,
                    transparent=True)

    def plot_class_by_domain(self, sequence, class_):

        class_percentages = {}
        percentages = sequence.domain_percentages_by_class(class_, True)

        domains = sequence.get_domains()
        for i in range(len(domains)):

            if domains[i].get_name() not in list(class_percentages.keys()):
                class_percentages[domains[i].get_name()] = [percentages[i]]
            else:
                class_percentages[domains[i].get_name()].append(percentages[i])

        categories = list(class_percentages.keys())
        class_values = []

        for category in categories:
            total = sum(class_percentages[category])
            class_values.append(total)

        plt.bar(categories, class_values, label=categories, color="coral")
        plt.ylabel(f"# {class_.name}")
        plt.title(f"{class_.name} by domain for {sequence.get_gene()}")
        plt.tight_layout()
        plt.show()