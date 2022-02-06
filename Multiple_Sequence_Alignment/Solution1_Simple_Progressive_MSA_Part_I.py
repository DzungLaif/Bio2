#Step 1: Create class for biological sequences.
class MySeq:
    def __init__(self, seq, seq_type = "DNA"):
        self .seq = seq.upper()
        self .seq_type = seq_type
    def __len__(self):
        return len(self.seq)
    def __getitem__(self, n):
        return self.seq[n]
    def __getslice__(self, i, j):
        return self.seq[i:j]
    def __str__(self):
        return self.seq
    def get_seq_biotype (self):
        return self.seq_type
    def show_info_seq ( self):
        print ("Sequence: " + self.seq + " biotype: " + self.seq_type)
    def alphabet(self):
        if(self.seq_type == "DNA"): return "ACGT"
        elif (self.seq_type == "RNA"): return "ACGU"
        elif (self.seq_type == "PROTEIN"): return "ACDEFGHIKLMNPQRSTVWY"
        else: return None
    def validate(self):
        alp = self.alphabet()
        res = True
        i = 0
        while i < len(self.seq) and res:
            if self.seq[i] not in alp: res = False
            else: i += 1
        return res
    def transcription(self):
        if(self.seq_type == "DNA"):
            return MySeq(self.seq.replace("T", "U"), "RNA")
        else:
            return None
    def reverse_comp(self):
        if(self.seq_type != "DNA"): return None
        comp = ""
        for c in self.seq:
            if(c == "A"): comp = "T" + comp
            elif (c == "T"): comp = "A" + comp
            elif (c == "G"): comp = "C" + comp
            elif (c == "C"): comp = "G" + comp
        return MySeq(comp, "DNA")
    def translate (self, iniPos = 0):
        if (self.seq_type != "DNA"): return None
        seq_aa = ""
        for pos in range (iniPos, len(self.seq)-2,3):
            cod = self.seq[pos:pos+3]
            seq_aa += self.translate_codon(cod)
        return MySeq(seq_aa, "PROTEIN")

    def translate_codon(self, cod):
        """Translates a codon into an aminoacid using an internal
        dictionary with the standard genetic code."""
        tc = {"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TGT": "C", "TGC": "C",
              "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TTT": "F", "TTC": "F",
              "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", "CAT": "H", "CAC": "H",
              "ATA": "I", "ATT": "I", "ATC": "I", "AAA": "K", "AAG": "K","TTA": "L",
              "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "ATG":"M",
              "AAT": "N", "AAC": "N", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
              "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
              "AGA": "R", "AGG": "R", "TCT":"S", "TCC": "S", "TCA": "S", "TCG": "S",
              "AGT": "S", "AGC": "S", "ACT":"T", "ACC": "T", "ACA": "T", "ACG": "T",
              "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TGG": "W", "TAT": "Y",
              "TAC": "Y", "TAA": "_", "TAG": "_", "TGA": "_"}
        if cod in tc: return tc[cod]
        else: return None

# Alignments using two variables: the alignment type (DNA,RNA, protein)
# and the list of sequences included in the alignment that will be defined
# as strings, including the symbol “-” to represent gaps
class MyAlign:
    def __init__(self, listseqs, al_types="protein"):
        self.listOfSeqs = listseqs
        self.allTypes = al_types

    def __len__(self):
        return len(self.listOfSeqs[0])

    def __getitem__(self, n):
        if type(n) is tuple and len(n) == 2:
            i, j = n
            return self.listOfSeqs[i][j]
        elif type(n) is int:
            return self.listOfSeqs[n]
        return None

    def __str__(self):
        res = ""
        for seq in self.listOfSeqs:
            res += "\n" + seq
        return res

    def num_seqs(self):
        return len(self.listOfSeqs)

    def column(self, indice):
        res = []
        for k in range(len(self.listOfSeqs)):
            res.append(self.listOfSeqs[k][indice])
        return res

    # Find consensus which is defined as the sequence
    # composed of the most frequent character
    # in each column of the alignment, ignoring gaps.
    def concecus(self):
        cons = ""
        for i in range(len(self)):
            cont = {}
            for k in range(len(self.listseqs)):
                c = self.listseqs[k][i]
            if c in cont:
                cont[c] = cont[c] + 1
            else:
                cont[c] = 1
            maximum = 0
            cmax = None
            for ke in cont.keys():
                if ke != "−" and cont[ke] > maximum:
                    maximum = cont[ke]
                    cmax = ke
            cons = cons + cmax
        return cons

#Step2: Create submatrix
class SubstMatrix:
    def __init__( self ):
        self .alphabet = ""
        self .sm = {}

    def __getitem__( self , ij):
        i, j = ij
        return self .score_pair(i, j)

    def score_pair( self , c1, c2):
        if c1 not in self.alphabet or c2 not in self.alphabet: return None
        return self.sm[c1+c2]

    #?? sep
    def read_submat_file(self, filename, sep):
        self.sm = {}
        f = open(filename, "r")
        line = f.readline()
        tokens = line.split("\t")
        ns = len(tokens)
        alphabet = []
        for i in range(0, ns):
            alphabet.append(tokens[i][0])
        for i in range(0, ns):
            line = f.readline();
            tokens = line.split("\t")
            for j in range(0, len(tokens)):
                k = alphabet[i] + alphabet[j]
                self.sm[k] = int(tokens[j])
        return self.sm

    def create_submat(self, match, mismatch, alphabet):
        sm = {}
        for c1 in alphabet:
            for c2 in alphabet:
                if c1 == c2: sm[c1 + c2] = match
                else: sm[c1 + c2] = mismatch
        return sm

# Step3: Pairwase Alignment
class PairwiseAlignment:
    def __init__(self, sm, g):
        self.g = g
        self.sm = sm
        self.S = None
        self.T = None
        self.seq1 = None
        self.seq2 = None

    def score_pos (self, c1, c2):
        if c1 == "−" or c2=="−":
            return self.g
        else:
            return self.sm[c1+c2]
    #?? align
    def score_align(self, align):
        res = 0;
        for i in range(len(self.seq1)):
            res += self.score_pos(self.seq1[i], self.seq2[i], self.sm, self.g)
        return res

    def max3traceback(self, v1, v2, v3):
        if v1 > v2:
            if v1 > v3: return 1
            else: return 3
        else:
            if v2 > v3: return 2
            else: return 3

    def needleman_Wunsh_Algorithsm(self, seq1, seq2):
        if (seq1.seq_type != seq2.seq_type): return None
        else:
            self.S = [[0]]
            self.T = [[0]]
            # initialize gap's row
            for j in range(1, len(seq2) + 1):
                (self.S[0][j]).append(self.g * j)
                (self.T[0][j]).append(3)
            # initialize gap's column
            for i in range(1, len(seq1) + 1):
                (self.S[i][0]).append([self.g * i])
                (self.T[i][0]).append(2)
            # use a recurrence relation to fill the remaining of the matrix
            for i in range(0, len(seq1)):
                for j in range(len(seq2)):
                    s1 = self.S[i][j] + self.score_pos(seq1[i], seq2[j], self.sm, self.g)
                    s2 = self.S[i][j + 1] + self.g
                    s3 = self.S[i + 1][j] + self.g
                    self.S[i + 1][j].append(max(s1, s2, s3))
                    self.T[i + 1][j].append(self.max3traceback(s1, s2, s3))
        return self.S[len(seq1)][len(seq2)]

    def recover_align(self):
        res = ["", ""]
        i = len(self.seq1)
        j = len(self.seq2)
        while i > 0 or j > 0:
            if self.T[i][j] == 1:
                res[0] = self.seq1[i - 1] + res[0]
                res[1] = self.seq2[j - 1] + res[1]
                i -= 1
                j -= 1
            elif self.T[i][j] == 3:
                res[0] = "−" + res[0]
                res[1] = self.seq2[j - 1] + res[1]
                j -= 1
            else:
                res[0] = self.seq1[i - 1] + res[0]
                res[1] = "−" + res[1]
                i -= 1
        return MyAlign(res, self.seq1.seq_type)

    def smith_Waterman(self, seq1, seq2):
        if(seq1.seq_type != seq2.seq_type):return None
        else:
            self.S = [[0]]
            self.T = [[0]]
            maxscore = 0
            # initialize gap's row
            for j in range(1, len(seq2) + 1):
                (self.S[0][j]).append(0)
                (self.T[0][j]).append(0)
            # initialize gap's column
            for i in range(1, len(seq1) + 1):
                (self.S[i][0]).append(0)
                (self.T[i][0]).append(0)
            # use a recurrence relation to fill the remaining of the matrix
            for i in range(0, len(seq1)):
                for j in range(len(seq2)):
                    s1 = self.S[i][j] + self.score_pos(seq1[i], seq2[j], self.sm, self.g)
                    s2 = self.S[i][j + 1] + self.g
                    s3 = self.S[i + 1][j] + self.g
                    b = max(s1, s2, s3)
                    if b <= 0:
                        self.S[i + 1].append(0)
                        self.T[i + 1].append(0)
                    else:
                        self.S[i + 1].append(b)
                        self.T[i + 1].append(self.max3traceback(s1, s2, s3))
                        if b > maxscore:
                            maxscore = b
        return maxscore

    def recover_align_local(self):
        res = ["", ""]
        i = len(self.seq1)
        j = len(self.seq2)
        while i > 0 or j > 0:
            if self.T[i][j] == 1:
                res[0] = self.seq1[i - 1] + res[0]
                res[1] = self.seq2[j - 1] + res[1]
                i -= 1
                j -= 1
            elif self.T[i][j] == 3:
                res[0] = "−" + res[0]
                res[1] = self.seq2[j - 1] + res[1]
                j -= 1
            else:
                res[0] = self.seq1[i - 1] + res[0]
                res[1] = "−" + res[1]
                i -= 1
        return MyAlign(res, self.seq1.seq_type)