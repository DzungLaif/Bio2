import re

class Regex_solution_ADN():
    @classmethod
    def __init__(self,seq,pattern):
        self.sequence = seq
        self.pattern = pattern

    def find_pattern_regex(seq, pat):
        from re import search
        mo = search(pat, seq)
        if (mo != None):
            return mo.span()[0]
        else :
            return -1

    def find_all_occurrences_re(seq, pat):
        from re import finditer
        mos = finditer(pat, seq)
        res = []
        for x in mos:
            res.append(x.span()[0])
        return res

def test():
    seq = input("Input sequence:")
    pat = input("Input pattern (as a regular expression):")
    new = Regex_solution_ADN(seq,pat)

    res = new.find_all_occurrences_re(seq, pat)
    if res >= 0:
        print("Pattern found in position: ", res)
    else :
        print("Pattern not found")

    all_res = new.find_all_occurrences_re(seq, pat)
    if len(all_res) > 0:
        print("Pattern found in positions: ", all_res)
    else:
        print("Pattern not found")
test()
