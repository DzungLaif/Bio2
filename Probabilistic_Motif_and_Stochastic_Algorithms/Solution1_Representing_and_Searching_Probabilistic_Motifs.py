# -*- coding: utf-8 -*-
def create_matrix_zeros (nrows, ncols):
    res = [ ]
    for i in range(0, nrows):
        res.append([0]*ncols)
    return res

def print_matrix(mat):
    for i in range(0, len(mat)): print(mat[i])

def translate_codon(cod):
    """Translates a codon into an aminoacid using an internal dictionary with the standard genetic code."""
    tc = {"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
          "TGT": "C", "TGC": "C",
          "GAT": "D", "GAC": "D",
          "GAA": "E", "GAG": "E",
          "TTT": "F", "TTC": "F",
          "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
          "CAT": "H", "CAC": "H",
          "ATA": "I", "ATT": "I", "ATC": "I",
          "AAA": "K", "AAG": "K",
          "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
          "ATG": "M", "AAT": "N", "AAC": "N",
          "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
          "CAA": "Q", "CAG": "Q",
          "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
          "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
          "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
          "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
          "TGG": "W",
          "TAT": "Y", "TAC": "Y",
          "TAA": "_", "TAG": "_", "TGA": "_"}
    if cod in tc:
        return tc[cod]
    else:
        return None


class MySeq:
    """ Class for biological sequences. """

    def __init__(self, seq, seq_type="DNA"):
        self.seq = seq.upper()
        self.seq_type = seq_type

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, n):
        return self.seq[n]

    def __getslice__(self, i, j):
        return self.seq[i:j]

    def __str__(self):
        return self.seq

    def get_seq_biotype(self):
        return self.seq_type

    def show_info_seq(self):
        print("Sequence: " + self.seq + " biotype: " + self.seq_type)

    def alphabet(self):
        if (self.seq_type == "DNA" or self.seq_type == "dna"):
            return "ACGT"
        elif (self.seq_type == "RNA" or self.seq_type == "rna"):
            return "ACGU"
        elif (self.seq_type == "PROTEIN" or self.seq_type == "protein"):
            return "ACDEFGHIKLMNPQRSTVWY"
        else:
            return None

    def validate(self):
        alp = self.alphabet()
        res = True
        i = 0
        while i < len(self.seq) and res:
            if self.seq[i] not in alp:
                res = False
            else:
                i += 1
        return res

    def transcription(self):
        if (self.seq_type == "DNA"):
            return MySeq(self.seq.replace("T", "U"), "RNA")
        else:
            return None

    def reverse_comp(self):
        if (self.seq_type != "DNA"): return None
        comp = ""
        for c in self.seq:
            if (c == 'A'):
                comp = "T" + comp
            elif (c == "T"):
                comp = "A" + comp
            elif (c == "G"):
                comp = "C" + comp
            elif (c == "C"):
                comp = "G" + comp
        return MySeq(comp, "DNA")

    def translate(self, iniPos=0):
        if (self.seq_type != "DNA"): return None
        seq_aa = ""
        for pos in range(iniPos, len(self.seq) - 2, 3):
            cod = self.seq[pos:pos + 3]
            seq_aa += translate_codon(cod)
        return MySeq(seq_aa, "PROTEIN")

class MyMotifs:

    def __init__(self, seqs = [], pwm = [], alphabet = None):
        if seqs:
            self.size = len(seqs[0])
            self.seqs = seqs # objet from class MySeq
            self.alphabet = seqs[0].alphabet()
            self.do_counts()
            self.create_pwm()
        else:
            self.pwm = pwm
            self.size = len(pwm[0])
            self.alphabet = alphabet

    def __len__ (self):
        return self.size

    def do_counts(self):
        self.counts = create_matrix_zeros(len(self.alphabet), self.size)
        for s in self.seqs:
            for i in range(self.size):
                lin = self.alphabet.index(s[i])
                self.counts[lin][i] += 1

    def create_pwm(self):
        if self.counts == None: self.do_counts()
        self.pwm = create_matrix_zeros(len(self.alphabet), self.size)
        for i in range(len(self.alphabet)):
            for j in range(self.size):
                self.pwm[i][j] = float(self.counts[i][j]) / len(self.seqs)

    def consensus(self):
        """ returns the sequence motif obtained with the most frequent symbol at each position of the motif"""
        res = ""
        for j in range(self.size):
            maxcol = self.counts[0][j]
            maxcoli = 0
            for i in range(1, len(self.alphabet) ):
                if self.counts[i][j] > maxcol:
                    maxcol = self.counts[i][j]
                    maxcoli = i
            res += self.alphabet[maxcoli]
        return res

    def masked_consensus(self):
        """ returns the sequence motif obtained with the symbol that occurrs in at least 50% of the input sequences"""
        res = ""
        for j in range(self.size):
            maxcol = self.counts[0][j]
            maxcoli = 0
            for i in range(1, len(self.alphabet) ):
                if self.counts[i][j] > maxcol:
                    maxcol = self.counts[i][j]
                    maxcoli = i
            if maxcol > len(self.seqs) / 2:
                res += self.alphabet[maxcoli]
            else:
                res += "-"
        return res

    def probability_sequence(self, seq):
        res = 1.0
        for i in range(self.size):
            lin = self.alphabet.index(seq[i])
            res *= self.pwm[lin][i]
        return res

    def probability_all_positions(self, seq):
        res = []
        for k in range(len(seq)-self.size+1):
            res.append(self.probability_sequence(seq))
        return res

    def most_probable_sequence(self, seq):
        """ Returns the index of the most probable sub-sequence of the input sequence"""
        maximum = -1.0
        maxind = -1
        for k in range(len(seq)-self.size):
            p = self.probability_sequence(seq[k:k+self.size])
            if(p > maximum):
                maximum = p
                maxind = k
        return maxind

    def create_motif(self, seqs):
        l = []
        for s in seqs:
            ind = self.most_probable_sequence(s.seq)
            subseq = MySeq(s[ind:(ind+self.size)], s.get_seq_biotype() )
            l.append(subseq)

        return MyMotifs(l)


def test():
    seq1 = MySeq("AAAGTT")
    seq2 = MySeq("CACGTG")
    seq3 = MySeq("TTGGGT")
    seq4 = MySeq("GACCGT")
    seq5 = MySeq("AACCAT")
    seq6 = MySeq("AACCCT")
    seq7 = MySeq("AAACCT")
    seq8 = MySeq("GAACCT")
    lseqs = [seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8]
    motifs = MyMotifs(lseqs)
    print ("Counts matrix")
    print_matrix (motifs.counts)
    print ("PWM")
    print_matrix (motifs.pwm)
    print ("Sequence alphabet")
    print(motifs.alphabet)

    print(motifs.probability_sequence("AAACCT"))
    print(motifs.probability_sequence("ATACAG"))
    print(motifs.most_probable_sequence("CTATAAACCTTACATC"))

    for s in lseqs:
        print (s)
    print ("Consensus sequence")
    print(motifs.consensus())
    print ("Masked Consensus sequence")
    print(motifs.masked_consensus())

    s1 = MySeq("TAAAGTTATGA")
    s2 = MySeq("ATGACACGTG")
    s3 = MySeq("TTTGGGTAT")

    newmotif = motifs.create_motif([s1,s2,s3])
    print_matrix(newmotif.counts)
    print(newmotif.most_probable_sequence("AAAACT"))


if __name__ == '__main__':
    test()