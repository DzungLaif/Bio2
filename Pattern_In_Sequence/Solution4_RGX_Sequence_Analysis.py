import re

class Regex_Sequence_Analysis():
    @classmethod
    def __init__(self, seq, pattern):
        self.sequence = seq
        self.pattern = pattern

    # Find if a pattern  on a sequence
    def validate_dna_re (self):
        from re import search
        if search("[^pattern]",self.seq) != None:
            return False
        else:
            return True

    # translate a DNA sequence to protein
    def translate_codon_re(self):
        table = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
        }
        protein =""
        if len(self.sequence) % 3 == 0:
            for i in range(0, len(self.sequence), 3):
                codon = self.sequence[i:i+3]
                protein += table[codon]
        return protein

    #  Find largest protein start with "M" and end with "-"
    def find_largest_protein(self):
        mos = re.finditer("M[^_]*_", self.sequence)
        sizem = 0
        largest_prot = ""
        for x in mos:
            init = x.span()[0]
            fin = x.span()[1]
            s = fin - init + 1
            if s > sizem:
                largest_prot = x.group()
                sizem = s
        return largest_prot

def test():
    seq = input("Input sequence:")
    pat = input("Input pattern (as a regular expression):")
    new = Regex_Sequence_Analysis(seq,pat)

    res1 = new.translate_codon_re(seq, pat)
    res2 = new.find_largest_protein(seq,pat)
    print("Protein translated from the DNA sequence: ", res1)
    print("Largest protein is: ", res2)

test()

