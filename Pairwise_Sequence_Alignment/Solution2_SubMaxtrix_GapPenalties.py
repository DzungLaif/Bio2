#function

# Create customs substitution matrices
def create_submat(match, mismatch, alphabet):
    sm = {}
    for c1 in alphabet:
        for c2 in alphabet:
            if c1 == c2:
                sm[c1+c2] = match
            else:
                sm[c1+c2] = mismatch
    return sm

# Load these matrices from files
def read_submat_file(filename):
    sm = {}
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
            k = alphabet[i]+alphabet[j]
            sm[k] = int(tokens[j])
    return sm

# Calculate the alignment score with g parameter is the score
# given to a gap in the sequence
def score_pos (c1, c2, sm, g):
    if c1 == "−" or c2 == "−":
        return g
    else:
        return sm[c1+c2]

def cal_score_align(seq1, seq2, sm, g):
    res = 0;
    for i in range(len(seq1)):
        res += score_pos (seq1[i], seq2[i], sm, g)
    return res

