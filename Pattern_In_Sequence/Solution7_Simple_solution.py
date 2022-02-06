class naivie_guy():
    def search_first_occ(self, seq, pattern):
        found = False
        i = 0
        while i <= len(seq) - len(pattern) and not found:
            j = 0
            while j < len(pattern) and pattern[j] == seq[i+j]:
                j = j+1
            if j == len(pattern):
                found = True
            else:
                i += 1
        if found:
            return i
        else:
            return -1

    def search_all_occs(self, seq, pattern):
        res = []
        for i in range(len(seq) - len(pattern) +1):
            j = 0
            while j < len(pattern) and pattern[i] == seq[i+j]:
                j = j+1
            if j == len(pattern):
                res.append(i)
        return res

def test():
    seq = input("Input sequence: ")
    pat = input("Input pattern: ")
    _ = naivie_guy.search_all_occs(seq, pat)

    print(pat, " can be found in given sequence at position:", )
    print(_)

test()