import math

import matplotlib.pyplot as plt
import pandas as pd

#
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import translation


class Plotter:

    # Define amino acid label colors
    aa_colors = {"Ala": "indianred", "Arg": "lavender", "Asn": "royalblue", "Asp": "blueviolet", "Cys": "darkmagenta",
                 "Gln": "goldenrod", "Glu": "deeppink", "Gly": "salmon", "His": "olivedrab", "Ile": "peru",
                 "Leu": "darkorchid", "Lys": "thistle", "Met": "darkseagreen", "Phe": "slateblue",
                 "Pro": "aquamarine", "Ser": "dimgray", "Stop": "darkgreen", "Thr": "lightyellow", "Trp": "rosybrown",
                 "Tyr": "cadetblue", "Val": "darkorange"}

    def __init__(self):
        pd.set_option("display.max_rows", 200, "display.max_columns",
                      200)

    def convert_aa_data_to_df(self, codon_data, by_proportion=False):

        """
        Converts AACodonFrequency object list into a Pandas dataframe

        Input:
        codon_data : [AACodonFrequency]
            List of AACodonFrequency objects holding codon frequency data
        by_proportion : bool (default=False)
            Fills codon data by proportion instead of absolute values

        """

        # Dictionary to hold codon data
        _dict = {}
        # List to hold totals
        totals = []
        for obj in codon_data:

            # Get copy of the codon frequencies
            freq_dict = obj.get_freq_dict().copy()

            # Convert to proportions
            if by_proportion:

                total_codons = obj.get_total_codons()

                for key in freq_dict.keys():
                    freq_dict[key] = freq_dict[key] / total_codons * 100

            # Add frequency to corresponding amino acid
            _dict[obj.get_abr()] = freq_dict
            totals.append(sum(freq_dict.values()))




        # Make pandas dataframe using created _dict
        df = pd.DataFrame(_dict)
        df = df.fillna(0)
        df = df.transpose()
        df["Total"] = totals

        # # Add empty amino acids not translated by sequence
        # all_aa = list(translation.aa_translating_codons.keys())
        # for aa in all_aa:
        #
        #     if aa not in list(df.index):
        #         df.loc[len(df)] = [0.00001] * len(df.columns)
        #         df = df.rename(index = lambda x: aa)


        #(df)
        return df

    def plot_stacked_bar(self, codon_data, organism, gene, by_proportion=False):

        # Convert codon data to dataframe
        df = self.convert_aa_data_to_df(codon_data, by_proportion)
        row_names = list(df.index)
        columns_names = df.columns

        # figure and axis
        fig, ax = plt.subplots(1, figsize=(12, 10))

        # plot bars
        left = len(row_names) * [0]
        print(len(left))
        for idx, name in enumerate(columns_names):

            col_vals = df.iloc[:, idx].tolist()
            ln = len(row_names)
            ln2 = len(col_vals)

            df = df

            plt.barh(row_names, col_vals, left=left)

            for i in range(len(left)):
                left[i] += col_vals[i]

        # title, legend, labels
        plt.title(f'Gene: {gene}\n', loc='left')
        plt.suptitle(f'Codon Frequency in {organism}\n')
        plt.legend(columns_names, bbox_to_anchor=([0.55, 1, 0, 0]), ncol=4,
                   frameon=False)
        plt.tight_layout()
        plt.xlabel('Frequency')

        # remove spines
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        # adjust limits and draw grid lines
        plt.ylim(-0.5, ax.get_yticks()[-1] + 0.5)
        ax.set_axisbelow(True)
        ax.xaxis.grid(color='gray', linestyle='dashed')

        plt.show()

    def generate_color_list(self, size):

        """ Generate list of alternating colours (first and second)"""

        first = [1 / 255 * 255, 1 / 255 * 102, 1 / 255 * 102]
        second = [1 / 255 * 153, 1 / 255 * 51, 1 / 255 * 255]

        colors = []

        for i in range(size):

            if i % 2 == 0:
                colors.append(first)
            else:
                colors.append(second)

        return colors

    def generate_aa_colours(self, aa_list):

        colors = []
        for i in range(len(aa_list)):
            colors.append(self.aa_colors[aa_list[i]])

        return colors

    def plot_nested_pie_from_dict(self, codon_data, organism, gene):

        # Convert data into dataframe
        df = self.convert_aa_data_to_df(codon_data, True)
        df = df.sort_index()
        df = df.drop('Total', 1)
        print(df)
        row_names = list(df.index)
        column_names = df.columns

        # Create inner dictionary
        inner_dict = {}
        for i in range(len(row_names)):

            row_vals = df.iloc[i, :].tolist()

            total = 0
            for val in row_vals:
                total += val

            inner_dict[row_names[i]] = total

        # Create outer dictionary
        codon_frequencies = {}

        for i in range(len(row_names)):

            row_vals = df.iloc[i, :].tolist()


            inner_freqs = {}
            for j in range(len(row_vals)):
                if row_vals[j] > 0:
                    inner_freqs[column_names[j]] = row_vals[j]

            inner_total = sum(list(inner_freqs.values()))
            temp_dict = {}
            for key in inner_freqs.keys():
                percentage = round(inner_freqs[key] / inner_total * 100)
                new_key = f"{key} {percentage}%"
                temp_dict[new_key] = inner_freqs[key]

            inner_freqs = temp_dict



            inner_freqs = dict(sorted(inner_freqs.items(), key=lambda item: item[1]))
            codon_frequencies.update(inner_freqs)

        # Generate colors for aa wedges
        aa_colors = self.generate_aa_colours(list(inner_dict.keys()))
        # Generate alternating colors for codon labels
        codon_colors = self.generate_color_list(len(column_names))

        # Plot pie circles
        codon_patches, codon_texts = plt.pie(codon_frequencies.values(), labels=codon_frequencies.keys(),
                                             startangle=90, frame=True,
                                             colors=codon_colors, labeldistance=1.1,
                                             wedgeprops=dict(width=0.3,
                                                             edgecolor='w'),
                                             rotatelabels=270)
        aa_patches, aa_texts = plt.pie(inner_dict.values(),
                                       labels=inner_dict.keys()
                                       , radius=0.7, startangle=90,
                                       colors=aa_colors,
                                       labeldistance=0.7,
                                       wedgeprops=dict(width=0.3,
                                                       edgecolor='w'),
                                       rotatelabels=270)



        # Set label properties by looping over the aa text objects
        for i in range(len(aa_texts)):



            # If aa frequency is small move text outside pie
            aa_sum = sum(list(inner_dict.values()))
            aa_value = list(inner_dict.values())[i]

            if aa_value / aa_sum < 0.01:
                # Set AA font size
                aa_texts[i].set_fontsize(5)

                # Text coordinates
                x = aa_texts[i]._x
                y = aa_texts[i]._y

                # Set new coords
                aa_texts[i]._x = x - x /1.8
                aa_texts[i]._y = y - y /1.8
            elif aa_value / aa_sum > 0.2:
                # Set AA font size
                aa_texts[i].set_fontsize(7)
            else:
                # Set AA font size
                aa_texts[i].set_fontsize(5)





        # Set label properties by looping over the codon text objects
        angle = 0
        frequency_values = list(codon_frequencies.values())
        for i in range(len(codon_texts)):

            # Hide codons with frequency of 0
            if frequency_values[i] == 0:
                codon_texts[i].set(visible=False)
                continue

            # Set font size
            codon_texts[i].set_fontsize(6)

            # Text coordinates
            x = codon_texts[i]._x
            y = codon_texts[i]._y

            # Calc x and y where labels should move towards
            radial_x, radial_y = 2 * math.sin(
                math.radians(-angle)), 2 * math.cos(math.radians(-angle))
            # plt.plot(radial_x, radial_y, 'ro')
            # plt.plot([x,radial_x],[y,radial_y])

            # increment angle based on codon position
            angle += 360 / sum(frequency_values) * frequency_values[i]

            # Even indexed codons should have different appearance from uneven indexed codons
            if i % 2 == 0:

                # Set color and alpha
                codon_texts[i]._color = (
                1 / 255 * 255, 1 / 255 * 102, 1 / 255 * 102)
                codon_texts[i]._alpha = 0.7

                # Calc new coords
                new_x = x + (radial_x - x) / 1.7
                new_y = y + (radial_y - y) / 1.7

                # Set new coords
                codon_texts[i]._y = new_y
                codon_texts[i]._x = new_x

                # Plot label line
                plt.plot([codon_texts[i]._x,
                          codon_texts[i]._x - codon_texts[i]._x / 2.7],
                         [codon_texts[i]._y,
                          codon_texts[i]._y - codon_texts[i]._y / 2.7],
                         color=(
                         1 / 255 * 255, 1 / 255 * 102, 1 / 255 * 102, 0.2))


            else:

                # Set color and alpha
                codon_texts[i]._color = (
                1 / 255 * 153, 1 / 255 * 51, 1 / 255 * 255)
                codon_texts[i]._alpha = 0.5

                # Calc new coords
                new_x = x + (radial_x - x) / 100
                new_y = y + (radial_y - y) / 100

                # Set new coords
                codon_texts[i]._y = new_y
                codon_texts[i]._x = new_x

                # Plot label line
                plt.plot([codon_texts[i]._x,
                          codon_texts[i]._x - codon_texts[i]._x / 13],
                         [codon_texts[i]._y,
                          codon_texts[i]._y - codon_texts[i]._y / 13],
                         color=(
                         1 / 255 * 153, 1 / 255 * 51, 1 / 255 * 255, 0.2))


        # Make it a donut
        centre_circle = plt.Circle((0, 0), 0.4, color="beige", linewidth=0)
        fig = plt.gcf()
        fig.gca().add_artist(centre_circle)
        plt.axis('equal')

        # Other properties
        plt.tight_layout()
        # fig.suptitle(f'Codon Frequency\nOrg: {organism[:8]}... '
        #              f'\nGene: {gene[:8]}...', ha="left", x=0)

        # plt.legend(row_names, bbox_to_anchor=([0.2, 1, 0, 0]), ncol=2, loc="right")
        # plt.tight_layout()


        plt.show()
        fig.savefig(f"Pie_{organism}_{gene}", dpi=1000, transparent=True)

    def randrange(self, n, vmin, vmax):
        '''
        Helper function to make an array of random numbers having shape (n, )
        with each number distributed Uniform(vmin, vmax).
        '''
        return (vmax - vmin) * np.random.rand(n) + vmin
