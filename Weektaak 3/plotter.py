import math

import matplotlib.pyplot as plt
import pandas as pd


class Plotter:

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

        # Make pandas dataframe using created _dict
        df = pd.DataFrame(_dict)
        df = df.fillna(0)
        df = df.transpose()

        #print(df)
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

    def plot_nested_pie_from_dict(self, codon_data, organism, gene):

        # Convert data into dataframe
        df = self.convert_aa_data_to_df(codon_data)
        row_names = list(df.index)
        columns_names = df.columns

        # Create inner dictionary
        inner_dict = {}
        for i in range(len(row_names)):

            row_vals = df.iloc[i, :].tolist()

            total = 0
            for val in row_vals:
                total += val

            inner_dict[row_names[i]] = total

        colors = []

        # Create outer lists [codons] [frequencies]
        codons = df.columns
        codon_frequencies = []

        for i in range(len(columns_names)):

            row_vals = df.iloc[:, i].tolist()

            total = 0
            for val in row_vals:
                total += val

            codon_frequencies.append(total)

        # Generate alternating colors for codon labels
        colors = self.generate_color_list(len(codons))

        # Plot pie circles
        codon_patches, codon_texts = plt.pie(codon_frequencies, labels=codons,
                                             startangle=90, frame=True,
                                             colors=colors, labeldistance=1.1,
                                             wedgeprops=dict(width=0.3,
                                                             edgecolor='w'))
        aa_patches, aa_texts = plt.pie(inner_dict.values(),
                                       labels=inner_dict.keys()
                                       , radius=0.7, startangle=90,
                                       labeldistance=0.7,
                                       wedgeprops=dict(width=0.3,
                                                       edgecolor='w'))

        for text in aa_texts:
            text.set_fontsize(8)

        # Move special text
        aa_texts[-1]._y -= 0.1
        aa_texts[-1]._x -= 0.05
        aa_texts[-2]._x -= 0.05
        aa_texts[-3]._x -= 0.05
        aa_texts[0]._x += 0.05

        # Set label properties by looping over the codon text objects
        angle = 0
        for i in range(len(codon_texts)):

            # Hide codons with frequency of 0
            if codon_frequencies[i] == 0:
                codon_texts[i].set(visible=False)
                continue

            # Set font size
            codon_texts[i].set_fontsize(8)

            # Text coordinates
            x = codon_texts[i]._x
            y = codon_texts[i]._y

            # Calc x and y where labels should move towards
            radial_x, radial_y = 2 * math.sin(
                math.radians(-angle)), 2 * math.cos(math.radians(-angle))
            # plt.plot(radial_x, radial_y, 'ro')
            # plt.plot([x,radial_x],[y,radial_y])

            # increment angle based on codon position
            angle += 360 / sum(codon_frequencies) * codon_frequencies[i]

            # Even indexed codons should have different appearance from uneven indexed codons
            if i % 2 == 0:

                # Set color and alpha
                codon_texts[i]._color = (
                1 / 255 * 255, 1 / 255 * 102, 1 / 255 * 102)
                codon_texts[i]._alpha = 0.7

                # Calc new coords
                new_x = x + (radial_x - x) / 2
                new_y = y + (radial_y - y) / 2

                # Set new coords
                codon_texts[i]._y = new_y
                codon_texts[i]._x = new_x

                # Plot label line
                plt.plot([codon_texts[i]._x,
                          codon_texts[i]._x - codon_texts[i]._x / 2.9],
                         [codon_texts[i]._y,
                          codon_texts[i]._y - codon_texts[i]._y / 2.9],
                         color=(
                         1 / 255 * 255, 1 / 255 * 102, 1 / 255 * 102, 0.2))


            else:

                # Set color and alpha
                codon_texts[i]._color = (
                1 / 255 * 153, 1 / 255 * 51, 1 / 255 * 255)
                codon_texts[i]._alpha = 0.5

                # Calc new coords
                new_x = x + (radial_x - x) / 4
                new_y = y + (radial_y - y) / 4

                # Set new coords
                codon_texts[i]._y = new_y
                codon_texts[i]._x = new_x

                # Plot label line
                plt.plot([codon_texts[i]._x,
                          codon_texts[i]._x - codon_texts[i]._x / 4.2],
                         [codon_texts[i]._y,
                          codon_texts[i]._y - codon_texts[i]._y / 4.2],
                         color=(
                         1 / 255 * 153, 1 / 255 * 51, 1 / 255 * 255, 0.2))

        # Set special text properties
        codon_texts[0]._y += 0.1

        # Make it a donut
        centre_circle = plt.Circle((0, 0), 0.4, color='white', linewidth=0)
        fig = plt.gcf()
        fig.gca().add_artist(centre_circle)
        plt.axis('equal')

        # Other properties
        plt.tight_layout()
        fig.suptitle(f'Codon Frequency\nOrg: {organism[:8]}... '
                     f'\nGene: {gene}', ha="left", x=0)


        plt.show()
