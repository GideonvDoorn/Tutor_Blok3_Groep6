from amino_acid import *

class ProteinSequence:

    def __init__(self, header, sequence, organism, gene):
        self.__header = header
        self.__sequence = sequence
        self.__organism = organism
        self.__gene = gene
        self.__domains = []

    def get_header(self):
        return self.__header

    def get_sequence(self):
        return self.__sequence

    def get_length(self):
        return len(self.__sequence)

    def get_organism(self):
        return self.__organism

    def get_gene(self):
        return self.__gene

    def get_domains(self):
        return self.__domains

    def add_domain(self, name, start, end, description=""):
        self.__domains.append(ProteinDomain(name, start, end, description))

    def count_by_letter(self,letter):
        return self.__sequence.count(letter)

    def percentage_by_letter(self, letter, decimals=2):
        return round(self.__sequence.count(letter) / len(self.__sequence) * 100, decimals)

    def percentage_by_class(self, aa_class, absolutes=False, decimals=2):

        """ Returns percentage of given class of amino acid"""

        if type(aa_class) != AAClass:
            print("No valid AAClass given : ProteinSequence.percentage_by_class")
            return 0

        total_of_class = 0

        # Calculate for hydrophobic subclasses
        if aa_class == AAClass.hydrophobic or aa_class == AAClass.hydrophilic or aa_class == AAClass.in_between:

            for aa in aa_data:
                if aa.get_hydropathy_class() == aa_class:
                    total_of_class += self.count_by_letter(aa.get_letter())

        # If absolutes was specified, return value in absolutes
        if absolutes:
            return total_of_class
        else:
            return round(total_of_class / len(self.__sequence) * 100, decimals)

    def domain_percentages_by_class(self, aa_class, absolutes=False, decimals=2):

        """ Return all domain and the percentages of a given class"""

        if type(aa_class) != AAClass:
            print("No valid AAClass given : ProteinSequence.percentage_by_class")
            return 0

        domain_percentages = []

        for domain in self.__domains:

            start = domain.get_start()
            end = domain.get_end()
            sub_sequence = self.__sequence[start : end]

            domain_sequence = ProteinSequence(domain.get_name(), sub_sequence, "", "")

            domain_percentages.append(domain_sequence.percentage_by_class(aa_class,absolutes))

        return domain_percentages




    def calculate_hydropathy(self, as_class=False, decimals=2):

        total = 0
        for aa in aa_data:
            total += aa.get_hydropathy() * self.count_by_letter(aa.get_letter())

        if as_class:
            return get_hydropathy_class(total)
        else:
            return round(total, decimals)

    def calculate_mass(self, decimals=2):

        total_mass = 0
        for aa in aa_data:
            total_mass += aa.get_size() * self.count_by_letter(aa.get_letter())

        return round(total_mass, decimals)

    def predicted_charge(self, decimals=2):

        total_charge = 0
        for char in self.__sequence:
            if get_charge_by_letter(char) == AAClass.positive:
                total_charge += 1
            elif get_charge_by_letter(char) == AAClass.negative:
                total_charge -= 1

        if total_charge > 0:
            return AAClass.positive
        elif total_charge < 0:
            return  AAClass.negative
        else:
            return  AAClass.neutral


    def get_gene_name(self):

        """ Extracts gene name from HIV/SIV header"""

        end_index = 0
        try:
            gene_index = self.__header.index("[gene=")
            end_index = len("[gene=")
        except ValueError:

            try:
                gene_index = self.__header.index("[protein=")
                end_index = len("[protein=")
            except ValueError:
                return "No gene name found!"
        header_cut = self.__header[gene_index :]

        return header_cut[end_index : header_cut.index("]")]


    def is_env_protein(self):
        return "env" in self.__header


class ProteinDomain:

    def __init__(self, name, start, end, description=""):
        self.__name = name
        self.__start = start
        self.__end = end
        self.__description = description

    def get_name(self):
        return self.__name

    def get_start(self):
        return self.__start

    def get_end(self):
        return self.__end

    def get_description(self):
        return self.__description

    def get_length(self):
        return self.__end - self.__start