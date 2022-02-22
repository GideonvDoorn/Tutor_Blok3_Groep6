from translation import aa_translating_codons, codon_dict

class CDS:

    def __init__(self, header, sequence):
        self.__header = header
        self.__sequence = sequence
        self.__codon_frequencies = self.__calc_codon_frequencies()

    def get_header(self):
        return self.__header

    def get_sequence(self):
        return self.__sequence

    def get_codon_frequencies(self):
        return self.__codon_frequencies

    def get_gene_name(self):

        """ Extracts gene name from header"""

        end_index = 0
        try:
            gene_index = self.__header.index("[gene=")
            end_index = len("[gene=")
        except ValueError:

            try:
                gene_index = self.__header.index("[protein=")
                end_index = len("[protein=")
            except ValueError:
                return "Fumerate Hydratase"
        header_cut = self.__header[gene_index :]

        return header_cut[end_index : header_cut.index("]")]


    def is_env_protein(self):
        return "env" in self.__header

    def get_codons(self):

        """ Get codons from cds sequence """

        codons = []

        for i in range(0, len(self.__sequence), 3):
            codons.append(self.__sequence[i: i + 3])

        return codons

    def count_codon(self, search_codon):

        """ Count specific codon {search_codon} in sequence"""

        count = 0
        for codon in self.get_codons():
            if search_codon == codon:
                count += 1
        return count

    def __calc_codon_frequencies(self):

        """ Returns list of AACodonFrequency objects,
        that hold amino acid and corresponding codon + number of occurrences
         in sequence """

        frequencies = []

        for codon in self.get_codons():

            aa_abr = codon_dict[codon.lower()]

            aa_registered = False
            for aa in frequencies:
                if aa.get_abr() == aa_abr:
                    aa_registered = True
                    aa.incr_codon_freq(codon)


            if not aa_registered:
                new_entry = AACodonFrequency(aa_abr, aa_translating_codons[aa_abr])
                new_entry.incr_codon_freq(codon)
                frequencies.append(new_entry)

        return frequencies

    def __str__(self):
        return f"{self.__header}{self.__sequence}\n"


class AACodonFrequency():

    """ Holds amino acid name and codon frequencies"""

    def __init__(self, aa_abr, codons):
        self.__aa_abr = aa_abr
        self.__codons = codons
        self.__codon_frequencies = {}

        for codon in self.__codons:
            self.__codon_frequencies[codon] = 0

    def get_abr(self):
        return self.__aa_abr

    def codon_translates(self, codon):
        return codon in self.__codons

    def get_freq_dict(self):
        return self.__codon_frequencies

    def get_total_codons(self):

        total = 0
        for freq in self.__codon_frequencies.values():
            total += freq

        return total

    def incr_codon_freq(self, codon, incr = 1):

        try:
            self.__codon_frequencies[codon] += 1
        except KeyError:
            print(f"Codon {codon} not found!")


    def __str__(self):
        _str = self.__aa_abr

        for key in self.__codon_frequencies.keys():

            _str += f"\n {key} {self.__codon_frequencies[key]}"

        return _str





