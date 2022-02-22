from enum import Enum

class AAClass(Enum):
    hydrophilic = 0
    hydrophobic = 1
    inbetween = 2
    positive = 3
    negative = 4
    neutral = 5


class AminoAcid:

    def __init__(self, letter, abbreviation, full_name, size, hydropathy, charge, other_properties=""):
        self.__letter = letter
        self.__abbreviation = abbreviation
        self.__full_name = full_name
        self.__size = size
        self.__hydropathy = hydropathy
        self.__charge = charge
        self.__other_properties = other_properties

    def get_letter(self):
        return self.__letter
    def get_abbreviation(self):
        return self.__abbreviation
    def get___full_name(self):
        return self.__full_name
    def get_size(self):
        return self.__size
    def get_hydropathy(self):
        return self.__hydropathy
    def get_charge(self):
        return self.__charge
    def get_other_properties(self):
        return self.__other_properties

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
    elif hydropathy < -2:
        return AAClass.hydrophilic
    else:
        return AAClass.inbetween



aa_data = [AminoAcid('A', 'Ala', 'Alanine', 89.094,1.8,'neutral'),
           AminoAcid('R', 'Arg', 'Arginine', 174.203,-4.5,'positive'),
           AminoAcid('N', 'Asn', 'Asparagine', 132.119,-3.5,'neutral'),
           AminoAcid('D', 'Asp', 'Aspartic Acid', 133.104,-3.5,'negative'),
           AminoAcid('C', 'Cys', 'Cysteine', 121.154,2.5,'neutral'),
           AminoAcid('Q', 'Gln', 'Glutamine', 146.146,-3.5,'neutral'),
           AminoAcid('E', 'Glu', 'Glutamatic Acid', 147.131,-3.5, 'negative'),
           AminoAcid('G', 'Gly', 'Glycine', 75.067,-0.4,'neutral'),
           AminoAcid('H', 'His', 'Histidine', 155.156,-3.2,'neutral'),
           AminoAcid('I', 'Ile', 'Isoleucine', 131.175,4.5,'neutral'),
           AminoAcid('L', 'Leu', 'Leucine', 131.175,3.8,'neutral'),
           AminoAcid('K', 'Lys', 'Lysine', 146.189,-3.9,'positive'),
           AminoAcid('M', 'Met', 'Methionine', 149.208,1.9,'neutral'),
           AminoAcid('F', 'Phe', 'Phenylalanine', 165.192,2.8,'neutral'),
           AminoAcid('P', 'Pro', 'Proline', 115.132,-1.6,'neutral'),
           AminoAcid('S', 'Ser', 'Serine', 105.093,-0.8,'neutral'),
           AminoAcid('T', 'Thr', 'Threonine', 119.119,-0.7,'neutral'),
           AminoAcid('W', 'Trp', 'Tryptophan', 204.228,-0.9,'neutral'),
           AminoAcid('Y', 'Tyr', 'Tyrosine', 181.191,-1.3,'neutral'),
           AminoAcid('V', 'Val', 'Valine', 117.148,4.2,'neutral'),
           ]