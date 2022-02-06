import re
import os

class BoyerMoore:
    # init function
    def __init__(self, alphabet, pattern):
        self.alphabet = alphabet
        self.pattern = pattern
        self.preprocess()

    def preprocess(self):
        self.process_bad_character()
        self.process_good_suffix()

    def process_bad_character(self):
        # define a dictionary with key is alphabet symbol
        # and the value is right most position of the symbol in pattern
        # if a symbol do not include in alphabet, value is set to -1
        self.occ={}
        for symbol in self.alphabet:
            self.occ[symbol]=-1
        for j in range(len(self.pattern)):
            c = self.pattern[j]
            self.occ[c] = j

    def process_good_suffix(self):
        self.f = [0]*(len(self.pattern)+1)
        self.s = [0]*(len(self.pattern)+1)
        i = len(self.pattern)
        j = len(self.pattern)+1
        self.f[i] = j
        while i > 0:
            while j <= len(self.pattern) and self.pattern[i-1] != self.pattern[j-1]:
                if self.s[j] == 0: self.s[i] = j-1
                j = self.f[i]
            j -= 1
            j -= 1
            self.f[i] = j
        j = self.f[0]
        for i in range(len(self.pattern)):
            if self.s[i] == 0: self.s[i] = j
            if i == j: j = self.f[j]

    def search_pattern(self, text):
        res = []
        i = 0
        while i <= len(text) - len( self.pattern):
            j = len(self.pattern) - 1
            while j >= 0 and self.pattern[j] == text[j+i]: j -= 1
            if(j < 0):
                res.append(i)
                i += self.s[0]
            else:
                c = text[j + i]
                i += max(self.s[j + 1], j - self.occ[c])
        return res

def test():
    bm = BoyerMoore("ACTG", "ACCA")
    print (bm.search_pattern("ATAGAACCAATGAACCATGATGAACCATGGATACCCAACCACC"))

test()
