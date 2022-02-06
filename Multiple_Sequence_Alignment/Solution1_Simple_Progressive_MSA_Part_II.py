from Solution1_Simple_Progressive_MSA_Part_I import PairwiseAlignment
from Solution1_Simple_Progressive_MSA_Part_I import MyAlign
from Solution1_Simple_Progressive_MSA_Part_I import MySeq
from Solution1_Simple_Progressive_MSA_Part_I import SubstMatrix

class MultipleAlignment():
    def __init__(self, seqs, alignseq: PairwiseAlignment):
        self.seqs = seqs
        self.alignpars = alignseq

    def add_sequence_alignment(self, alignment: MyAlign, sequence):
        res = []
        for i in range(len(alignment.listOfSeqs)+1):
            res.append("")
        cons = MySeq(alignment.concecus(), alignment.al_type)
        self.alignpars.needleman_Wunsh_Algorithsm(cons, sequence)
        align2 = self.alignpars.recover_align()
        orig = 0
        for i in range(len(align2)):
            if align2[0, i] == '−':
                for k in range(len(alignment.listOfSeqs)):
                    res[k] += "−"
            else:
                for k in range(len(alignment.listOfSeqs)):
                    res[k] += alignment[k, orig]
                orig += 1
        res[len(alignment.listOfSeqs)] = align2.listseqs[1]
        return MyAlign(res, alignment.listOfSeqs)

    def align_consensus(self):
        self.alignpars.needleman_Wunsh_Algorithsm(self.seqs[0], self.seqs[1])
        res = self.alignpars.recover_align()
        for i in range(2, len(self.seqs)):
            res = self.add_sequence_alignment(res, self .seqs[i])
        return res

def test():
    s1 = MySeq("ATAGC")
    s2 = MySeq("AACC")
    s3 = MySeq("ATGAC")
    sm = SubstMatrix()
    sm.create_submat(1, -1, "ACGT")
    aseq = PairwiseAlignment(sm, -1)
    ma = MultipleAlignment([s1, s2, s3], aseq)
    al = ma.align_consensus()
    print(al)