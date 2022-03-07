from enum import Enum

class AAClass(Enum):
    hydrophilic = 0
    hydrophobic = 1
    in_between = 2
    positive = 3
    negative = 4
    neutral = 5


class AminoAcid:

    def __init__(self, letter, abbreviation, full_name, size, hydropathy, charge, special_properties=""):
        self.__letter = letter
        self.__abbreviation = abbreviation
        self.__full_name = full_name
        self.__size = size
        self.__hydropathy = hydropathy
        self.__charge = charge
        self.__special_properties = special_properties

    def get_letter(self):
        return self.__letter
    def get_abbreviation(self):
        return self.__abbreviation
    def get_full_name(self):
        return self.__full_name
    def get_size(self):
        return self.__size
    def get_hydropathy(self):
        return self.__hydropathy
    def get_charge(self):
        return self.__charge
    def get_other_properties(self):
        return self.__special_properties

    def get_hydropathy_class(self):
        return get_hydropathy_class(self.__hydropathy)

    def get_size_class(self):
        if self.__size > 140:
            return 'large'
        elif self.__size < 100:
            return 'small'
        else:
            return 'medium'

def get_hydropathy_class(hydropathy):
    if hydropathy > 0:
        return AAClass.hydrophobic
    elif hydropathy < 0:
        return AAClass.hydrophilic
    else:
        return AAClass.in_between

def get_charge_by_letter(letter):
    for aa in aa_data:
        if aa.get_letter() == letter:
            return aa.get_charge()

def get_full_name_by_letter(letter):
    for aa in aa_data:
        if aa.get_letter() == letter:
            return aa.get_full_name()

def get_abr_by_letter(letter):
    for aa in aa_data:
        if aa.get_letter() == letter:
            return aa.get_abbreviation()

aa_data = [AminoAcid('A', 'Ala', 'Alanine', 89.094, 1.8, AAClass.neutral),
           AminoAcid('R', 'Arg', 'Arginine', 174.203, -4.5, AAClass.positive),
           AminoAcid('N', 'Asn', 'Asparagine', 132.119, -3.5, AAClass.neutral),
           AminoAcid('D', 'Asp', 'Aspartic Acid', 133.104, -3.5,
                     AAClass.negative),
           AminoAcid('C', 'Cys', 'Cysteine', 121.154, 2.5, AAClass.neutral),
           AminoAcid('Q', 'Gln', 'Glutamine', 146.146, -3.5, AAClass.neutral),
           AminoAcid('E', 'Glu', 'Glutamatic Acid', 147.131, -3.5,
                     AAClass.negative),
           AminoAcid('G', 'Gly', 'Glycine', 75.067, -0.4, AAClass.neutral),
           AminoAcid('H', 'His', 'Histidine', 155.156, -3.2, AAClass.neutral),
           AminoAcid('I', 'Ile', 'Isoleucine', 131.175, 4.5, AAClass.neutral),
           AminoAcid('L', 'Leu', 'Leucine', 131.175, 3.8, AAClass.neutral),
           AminoAcid('K', 'Lys', 'Lysine', 146.189, -3.9, AAClass.positive),
           AminoAcid('M', 'Met', 'Methionine', 149.208, 1.9, AAClass.neutral),
           AminoAcid('F', 'Phe', 'Phenylalanine', 165.192, 2.8,
                     AAClass.neutral),
           AminoAcid('P', 'Pro', 'Proline', 115.132, -1.6, AAClass.neutral),
           AminoAcid('S', 'Ser', 'Serine', 105.093, -0.8, AAClass.neutral),
           AminoAcid('T', 'Thr', 'Threonine', 119.119, -0.7, AAClass.neutral),
           AminoAcid('W', 'Trp', 'Tryptophan', 204.228, -0.9, AAClass.neutral),
           AminoAcid('Y', 'Tyr', 'Tyrosine', 181.191, -1.3, AAClass.neutral),
           AminoAcid('V', 'Val', 'Valine', 117.148, 4.2, AAClass.neutral),
           ]