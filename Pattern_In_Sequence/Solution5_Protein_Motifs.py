import re

class Regex_Protein_Motifs():
    def __init__(self, seq, prof):
        self.sequence = seq
        self.profile = prof

    # Find Zinc finger RING-type signature” (PS00518) motif
    # which is represented by “C-x-H-x-[LIVMFY]-C-x(2)-C-[LIVMYA]
    def find_zync_finger(self):
        from re import search
        regexp = "C.H.[LIVMFY]C.{2}C[LIVMYA]"
        mo = search(regexp, self.sequence)
        if (mo != None):
            return mo.span()[0]
        else:
            return -1
    # Transform Prosite pattern into an Regex form
    # and find it in the sequence
    def find_prosite(self):
        from re import search
        regexp = self.profile.replace("−", "")
        regexp = self.profile.replace("x", ".")
        regexp = self.profile.replace("(", "{")
        regexp = self.profile.replace(")", "}")
        mo = search(regexp, self.sequence)
        if(mo != None):
            return mo.span()[0]
        else:
            return -1

def test():
    seq = input("Input sequence:")
    pat = input("Input prosite form:")
    new = Regex_Protein_Motifs(seq,pat)

    res1 = new.find_zync_finger(seq, pat)
    res2 = new.find_prosite(seq,pat)
    print("The zync finger pattern found in the sequence: ", res1)
    print("Regex form of the given prosite: ", res2)

test()
