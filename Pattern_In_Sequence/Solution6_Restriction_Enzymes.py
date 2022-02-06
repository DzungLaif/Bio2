class Regex_Restr_Enzymes:
    def __init__(self, enzyme, sequence):
        self.seq = sequence
        self.enz = enzyme

    def iub_to_RE(self, iub):
        dic = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[GA]", "Y": "[CT]","M": "[AC]", "K": "[GT]", "S": "[GC]", "W": "[AT]", "B": "[CGT]", "D":"[AGT]", "H":"[ACT]", "V":"[ACG]", "N":"[ACGT]"}
        site = iub.replace("^", "")
        regexp = ""
        for c in site:
            regexp += dic[c]
        return regexp

    # detect the cut position of an enzyme and
    # print the sub-sequence after the cut (restriction map)
    def cut_positions(self):
        from re import finditer
        cutpos = self.enz.find("^")
        regexp = self.iub_to_RE(self.enz)
        matches = finditer(regexp, self.seq)
        locs = []
        for m in matches:
            locs.append(m.start() + cutpos)
        return locs

    def cut_subsequences(self, locs):
        res = []
        positions = locs
        positions.insert(0, 0)
        positions.append(len(self.seq))
        for i in range(len(positions)-1):
            res.append(self.seq[positions[i]:positions[i + 1]])
            return res

def test():
    seq = input("Input sequence:")
    enz = input("Input enzyme:")
    new = Regex_Restr_Enzymes(seq, enz)

    res1 = new.cut_positions()
    res2 = new.cut_subsequences(res1)
    print("Cut position is: ", res1)
    print("The sub-sequence is: ", res2)

test()