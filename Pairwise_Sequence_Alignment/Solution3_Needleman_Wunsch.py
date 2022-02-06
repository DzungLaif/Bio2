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
        for j in range(0, len(tokens)):
            k = alphabet[i]+alphabet[j]
            sm[k] = int(tokens[j])
    return sm

# calculate the score of each cells in matrix
def score_pos (sq_pos1, sq_pos2, sm, g):
    if sq_pos1 == "−" or sq_pos2 == "−":
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
def needleman_Wunsh_Algorithsm (sequence1, sequence2, submax, gap):
    S = [[0]]
    T = [[0]]
    # initialize gap's row
    for j in range(1, len(sequence2)+1):
        (S[0][j]).append(gap*j)
        (T[0][j]).append(3)
    # initialize gap's column
    for i in range(1, len(sequence1)+1):
        (S[i][0]).append([gap*i])
        (T[i][0]).append(2)
    # use a recurrence relation to fill the remaining of the matrix
    for i in range(0, len(sequence1)):
        for j in range(len(sequence2)):
            s1 = S[i][j] + score_pos(sequence1[i], sequence2[j], submax, gap)
            s2 = S[i][j+1] + gap
            s3 = S[i+1][j] + gap
            S[i+1][j].append(max(s1, s2, s3))
            T[i+1][j].append(max3traceback(s1, s2, s3))
    return (S, T)

# restore the most similar sequence after traceback
def align_recover (trb_matrix, seq1, seq2):
    res = ["", ""]
    i = len(seq1)
    j = len(seq2)
    while i > 0 or j > 0:
        if trb_matrix[i][j] == 1:
            res[0] = seq1[i-1] + res[0]
            res[1] = seq2[j-1] + res[1]
            i -= 1
            j -= 1
        elif trb_matrix[i][j] == 3:
            res[0] = "−" + res[0]
            res[1] = seq2[j-1] + res[1]
            j -= 1
        else:
            res[0] = seq1[i-1] + res[0]
            res[1] = "−" + res[1]
            i -= 1
    return res

# print the result
def print_matrix(mtrx):
    for i in range(0,len(mtrx)):
        print(mtrx[i])

# test function
def test_global_alig():
    sm = read_submat_file("blosum62.mat")
    seq1 = "PHSWG"
    seq2 = "HGWAG"
    res = needleman_Wunsh_Algorithsm(seq1, seq2, sm, -8)
    S = res[0]
    T = res[1]

    print ("Score of optimal alignment:", S[ len (seq1)][ len (seq2)])
    print_matrix(S)
    print_matrix(T)
    alig = align_recover(T, seq1, seq2)
    print (alig[0])
    print (alig[1])

test_global_alig()