from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo

# Alignment of two DNA sequences, using match score of 1, mismatch and gap score are 0
def global_Align_Python(adn1, adn2):
    alignments = pairwise2.align.globalxx(adn1, adn2)

# Alignment of two DNA sequences using BLOSUM62 substitution matrix
# gap score are -4 and extension penalty of -1
def Align_using_BLOSUM62(prot1, prot2):
    matrix = MatrixInfo.blosum62
    for a in pairwise2.align.globalds(prot1, prot2, matrix, -4, -1):
        print(format_alignment(*a))

#perform local alignments: in the first case of DNA
#sequences using a match score of 3, mismatch score of −2 and
# constant gap penalty g of −3.
def local_Align_Python(adn1, adn2, prot1, prot2):
    matrix = MatrixInfo.blosum62
    # Align DNA sequence
    local_dna = pairwise2.align.localms(adn1, adn2, 3, -2, -3, -3)
    # Align protein sequence
    local_prot = pairwise2.align.localds(prot1, prot2, matrix, -4, -1)
    return local_dna, local_prot

def test():
    seqAdn_1 = input("Input DNA sequence 1:")
    seqAdn_2 = input("Input DNA enzyme 2:")
    seqPr_1 = input("Input protein sequence 1:")
    seqPr_2 = input("Input protein enzyme 2:")
    alignment = local_Align_Python(seqAdn_1, seqAdn_2, seqPr_1, seqPr_2)
    print("ADN alignment result is ", alignment[0])
    print("ADN alignment result is ", alignment[1])

test()