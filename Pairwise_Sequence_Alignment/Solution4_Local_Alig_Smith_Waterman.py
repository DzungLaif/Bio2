# read the complex substitution matrices from file
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
        line = f.readline()
        tokens = line.split("\t")
        for j in range(0, len (tokens)):
            k = alphabet[i]+alphabet[j]
            sm[k] = int (tokens[j])
    return sm

# calculate the score of each cells in matrix
def score_pos (sq_pos1, sq_pos2, sm, g):
    if sq_pos1 == "−" or sq_pos2=="−":
        return g
    else :
        return sm[sq_pos1+sq_pos2]

# Traceback filling rule
def max3traceback(v1, v2, v3):
    if v1 > v2:
        if v1 > v3: return 1
        else: return 3
    else:
        if v2 > v3: return 2
        else: return 3

# S matrix has n+1 rows and m+1 columns, the cell S(i,j) is filled by the following rule
# S(i,j) = max{ S(i-1,j-1) + sm[a(i),b(j)],S(i-1,j) + g, S(i,j-1)+g, with 0<i<= n, 0<j<=m }
# where sm(c1, c2) gives the value of the substitution matrix for symbols c1 and c2
# g is penalty for a gap
# T matrix is traceback matrix, which has n+1 rows and m+1 columns
def smith_Water_Algorithsm (sequence1, sequence2, submax, gap):
    S = [[0]]
    T = [[0]]
    maxscore = 0
    # initialize gap's row
    for j in range(1, len(sequence2)+1):
        (S[0][j]).append(0)
        (T[0][j]).append(0)
    # initialize gap's column
    for i in range(1, len(sequence1)+1):
        (S[i][0]).append(0)
        (T[i][0]).append(0)
    # use a recurrence relation to fill the remaining of the matrix
    for i in range(0, len(sequence1)):
        for j in range(len(sequence2)):
            s1 = S[i][j] + score_pos(sequence1[i], sequence2[j], submax, gap)
            s2 = S[i][j+1] + gap
            s3 = S[i+1][j] + gap
            b = max(s1, s2, s3)
            if b <= 0:
                S[i + 1].append(0)
                T[i + 1].append(0)
            else:
                S[i+1].append(b)
                T[i+1].append(max3traceback(s1, s2, s3))
                if b > maxscore:
                    maxscore = b
    return (S, T, maxscore)

def max_mat(mat):
    maxval = mat[0][0]
    maxrow = 0
    maxcol = 0
    for i in range (0, len (mat)):
        for j in range (0, len (mat[i])):
            if mat[i][j] > maxval:
                maxval = mat[i][j]
                maxrow = i
                maxcol = j
    return (maxrow,maxcol)

# restore the most similar sequence after traceback
def recover_align_local (S, T, seq1, seq2):
    res = ["", ""]
    i, j = max_mat(S)
    while T[i][j]>0:
        if T[i][j]==1:
            res[0] = seq1[i-1] + res [0]
            res[1] = seq2[j-1] + res [1]
            i -= 1
            j -= 1
        elif T[i][j] == 3:
            res[0] = "−" + res[0];
            res[1] = seq2[j-1] + res [1]
            j -= 1
        elif T[i][j] == 2:
            res[0] = seq1[i-1] + res [0]
            res[1] = "−" + res[1]
            i -= 1
    return res

# test function
def test_Local_Alignment():
    sm = read_submat_file("blosum62.mat")
    seq1 = "HGWAGWAGG"
    seq2 = "PHSWGGAGH"
    res = smith_Water_Algorithsm(seq1, seq2, sm, -8)
    S = res[0]
    T = res[1]
    print ("Score of optimal alignment:", res[2])
    print(S)
    print(T)
    alignment = recover_align_local(S, T, seq1, seq2)
    print(alignment[0])
    print(alignment[1])

test_Local_Alignment()