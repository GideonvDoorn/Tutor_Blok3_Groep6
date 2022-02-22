from amino_acid import aa_data, AAClass, get_hydropathy_class

class ProteinSequence:

    def __init__(self, header, sequence):
        self.__header = header
        self.__sequence = sequence

    def get_header(self):
        return self.__header

    def get_sequence(self):
        return self.__sequence

    def count_by_letter(self,letter):
        return self.__sequence.count(letter)

    def percentage_by_letter(self, letter, decimals=2):
        return round(self.__sequence.count(letter) / len(self.__sequence) * 100, decimals)

    def percentage_by_class(self, aa_class, decimals=2):

        if type(aa_class) != AAClass:
            print("No valid AAClass given : ProteinSequence.percentage_by_class")
            return 0

        total_of_class = 0
        if aa_class == AAClass.hydrophobic or aa_class == AAClass.hydrophilic or aa_class == AAClass.inbetween:

            for aa in aa_data:
                if aa.get_hydropathy_class() == aa_class:
                    total_of_class += self.count_by_letter(aa.get_letter())


        return round(total_of_class / len(self.__sequence) * 100, decimals)

    def calculate_hydropathy(self, as_class=False):

        total = 0
        for aa in aa_data:
            total += aa.get_hydropathy() * self.count_by_letter(aa.get_letter())

        if as_class:
            return get_hydropathy_class(total)
        else:
            return total
