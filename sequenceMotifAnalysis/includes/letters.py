'''
Created on 07.11.2021

@author: grustentier
'''

AMINO_ACID_LETTERS_ALPHABETICAL = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
AMINO_ACID_LETTERS_HYDROPHOB = ["A", "C", "F", "I", "L", "M", "P", "T", "V", "W", "Y", "D", "E", "G", "H", "K", "N", "Q", "R", "S"]
AMINO_ACID_LETTERS_SMALL = ["A", "C", "D", "G", "N", "P", "V", "S", "T", "E", "F", "H", "I", "K", "L", "M", "Q", "R", "W", "Y"]
AMINO_ACID_LETTERS_POLAR = ["C", "D", "E", "H", "K", "N", "Q", "R", "S", "T", "W", "Y", "A", "F", "G", "I", "L", "M", "P", "V"]


def getAminoAcidLetters(sort_by):
    if str(sort_by).lower() == "hydrophob":
        return AMINO_ACID_LETTERS_HYDROPHOB
    elif str(sort_by).lower() == "small":
        return AMINO_ACID_LETTERS_SMALL
    elif str(sort_by).lower() == "polar":
        return AMINO_ACID_LETTERS_POLAR
    else:return AMINO_ACID_LETTERS_ALPHABETICAL